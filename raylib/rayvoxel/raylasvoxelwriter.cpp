// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file implements I/O functions for writing a VoxelGrid to various file formats,
// including AMAPVox, text, and NetCDF.

#include "raylib/rayvoxel/raylasvoxelwriter.h"
#include "raylib/rayvoxel/raylasvegmetrics.h"
#include "raylib/rayvoxel/raylasbailey.h"
#include "raylib/rayvoxel/rayvox.h" // For writing the AMAPVox .vox format
#include "raylib/rayutils.h"
#include "raylib/rayunused.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <algorithm>
#include <cassert>

#if RAYLIB_WITH_NETCDF
#include <netcdf>
#endif // RAYLIB_WITH_NETCDF

#ifdef _MSC_VER
#include <intrin.h>
#define popcount __popcnt64
#else
#define popcount __builtin_popcountll
#endif

namespace ray
{

// ==================================================================================
// Local Helper Functions
// ==================================================================================
namespace { // Use an anonymous namespace for local helpers

// Helper to get the dominant point classification. Returns -1 if no hits.
int getDominantPointClassification(const std::array<float, 256>& classification_hits) {
    int dominant = -1;
    float max_hits = -1.0f;
    for (int c = 0; c < 256; ++c) {
        if (classification_hits[c] > max_hits) {
            max_hits = classification_hits[c];
            dominant = c;
        }
    }
    return (max_hits > 0.0f) ? dominant : -1;
}

// Helper to get the classification if and only if all hits in the voxel are of that class. Returns -1 otherwise.
int getAbsolutePointClassification(const std::array<float, 256>& classification_hits) {
    int found = -1;
    int nonzero_count = 0;
    for (int c = 0; c < 256; ++c) {
        if (classification_hits[c] > 0.0f) {
            ++nonzero_count;
            found = c;
            if (nonzero_count > 1) return -1;
        }
    }
    return (nonzero_count == 1) ? found : -1;
}

// Helper to format the classification hits map to a string
std::string formatClassificationHits(const std::array<float, 256>& hits) {
    bool any = false;
    for (int c = 0; c < 256; ++c) if (hits[c] > 0.0f) { any = true; break; }
    if (!any) return "{}";
    std::stringstream ss;
    ss << "{";
    bool first = true;
    for (int c = 0; c < 256; ++c) {
        if (hits[c] > 0.0f) {
            if (!first) ss << ";";
            ss << c << ":" << std::fixed << std::setprecision(2) << hits[c];
            first = false;
        }
    }
    ss << "}";
    return ss.str();
}

// Helper to parse comma-separated class strings
std::vector<int> parseClasses(const std::string& classes_str) {
    std::vector<int> classes;
    if (classes_str.empty()) {
        return classes;
    }
    std::stringstream ss(classes_str);
    std::string item;
    while (std::getline(ss, item, ',')) {
        try {
            classes.push_back(std::stoi(item));
        } catch (const std::invalid_argument& e) {
            std::cerr << "Warning: Invalid class code '" << item << "' in class list. Ignoring." << std::endl;
        }
    }
    return classes;
}

// Helper to parse comma-separated LAD parameter string
void parseLadParams(const std::string& lad_params_str, double& param1, double& param2) {
    param1 = 0.0;
    param2 = 0.0;
    if (lad_params_str.empty()) {
        return;
    }
    std::stringstream ss(lad_params_str);
    std::string item;
    if (std::getline(ss, item, ',')) {
        try {
            param1 = std::stod(item);
        } catch (const std::invalid_argument&) {}
    }
    if (std::getline(ss, item, ',')) {
        try {
            param2 = std::stod(item);
        } catch (const std::invalid_argument&) {}
    }
}

} // end anonymous namespace


// ==================================================================================
// Centralized Metric Calculator Implementation
// ==================================================================================

// Per-voxel attenuation-coefficient (lambda) estimators reproducing the AMAPVox
// PAD estimator family.
//   Stage 1 (echo weight only): free_path_length / num_beams_weighted
//   Stage 2 (beam metrics):     bs_free_path / bs_entering
//   PPL inputs:                 sum_hit_delta, sum_miss_delta (mean chord for hit/miss beams)
//
// Pimont, F., Allard, D., Soma, M. & Dupuy, J.-L. (2018). Estimators and confidence
// intervals for plant area density at voxel scale with T-LiDAR.
// Remote Sensing of Environment, 215, 343-370.
// DOI: 10.1016/j.rse.2018.06.024
// Source reference: https://github.com/umr-amap/AMAPVox
//
// "fpl":           lambda = intercepted / total traversed path (mean free path length).
// "transmittance": lambda = -ln(T), T = (entering - intercepted) / entering (Beer-Lambert).
// "ppl":           bias-corrected MLE separating mean hit/miss path lengths; solved by
//                  bisection; unbiased for > 2 beams per voxel (Pimont eq. for PPL).
double computeLambda(const VoxelGrid::Voxel& v, const std::string& method)
{
    const double eps = 1e-10;

    if (method == "fpl") {
        if (v.bs_free_path > eps)
            return static_cast<double>(v.bs_intercepted) / static_cast<double>(v.bs_free_path);
        if (v.free_path_length > eps)
            return static_cast<double>(v.num_hits) / static_cast<double>(v.free_path_length);
        return 0.0;
    }
    if (method == "transmittance") {
        if (v.bs_entering > eps) {
            double T = std::max(eps, static_cast<double>(v.bs_entering - v.bs_intercepted)
                                    / static_cast<double>(v.bs_entering));
            return -std::log(T);
        }
        if (v.num_beams_weighted > eps) {
            double hf = std::min(1.0 - eps,
                                 static_cast<double>(v.num_hits) / static_cast<double>(v.num_beams_weighted));
            return -std::log(1.0 - hf);
        }
        return 0.0;
    }
    if (method == "ppl") {
        double n      = static_cast<double>(v.num_hits);
        double m      = static_cast<double>(std::max(0, v.num_beams - static_cast<int32_t>(v.num_hits)));
        double dbar_n = (n > eps) ? static_cast<double>(v.sum_hit_delta)  / n : 0.0;
        double dbar_m = (m > eps) ? static_cast<double>(v.sum_miss_delta) / m : 0.0;
        if (n < eps || dbar_n < eps) {
            if (v.bs_entering > eps) {
                double T = std::max(eps, static_cast<double>(v.bs_entering - v.bs_intercepted)
                                        / static_cast<double>(v.bs_entering));
                return -std::log(T);
            }
            return 0.0;
        }
        if (m < eps || dbar_m < eps) {
            return 50.0 / dbar_n;
        }
        // Bisection: solve  m·δ̄_m = n·δ̄_n·exp(-λδ̄_n)/(1−exp(-λδ̄_n))
        double lo = eps;
        double hi = 50.0 / std::min(dbar_n, dbar_m);
        for (int iter = 0; iter < 60; ++iter) {
            double mid  = 0.5 * (lo + hi);
            double expm = std::exp(-mid * dbar_n);
            double den  = 1.0 - expm;
            if (den < eps) { hi = mid; continue; }
            if (m * dbar_m < n * dbar_n * expm / den) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }
    // unknown method — FPL simplified fallback
    if (v.path_length > eps)
        return static_cast<double>(v.num_hits) / static_cast<double>(v.path_length);
    return 0.0;
}

MetricResultsMap calculateOutputMetrics(const VoxelGrid& grid, const VoxelizationParameters& params,
                                         const HeightField* dtm, const ClassTable& class_table,
                                         const PerTreeIadMap& per_tree_iad,
                                         const PredominantTreeTable& predominant_tree)
{
    MetricResultsMap results;
    const double voxel_volume = grid.getVoxelWidth() * grid.getVoxelWidth() * grid.getVoxelWidth();
    const double vox_w = grid.getVoxelWidth();
    const Eigen::Vector3d& bmin = grid.getBounds().min_bound_;
    const auto& dims = grid.getDimensions();

    // Strip optional "field:" prefix before parsing numeric class codes
    auto strip_field_prefix = [](const std::string& s) -> std::string {
        auto colon = s.find(':');
        return (colon != std::string::npos) ? s.substr(colon + 1) : s;
    };
    const std::vector<int> leaf_classes = parseClasses(strip_field_prefix(params.leaf_classes_str));
    const std::vector<int> wood_classes = parseClasses(strip_field_prefix(params.wood_classes_str));
    double lad_param1, lad_param2;
    parseLadParams(params.lad_params_str, lad_param1, lad_param2);

    auto computeLambdaV = [&](const VoxelGrid::Voxel& v) -> double {
        return computeLambda(v, params.attenuation_methods[0]);
    };

    // Lambda to fill a VoxelOutputData from a voxel + its flat index
    auto populateData = [&](int64_t flat_idx, int64_t ci, int64_t cj, int64_t ck,
                            const VoxelGrid::Voxel& v) -> VoxelOutputData {
        VoxelOutputData data;
        data.i = ci; data.j = cj; data.k = ck;
        data.x = bmin.x() + vox_w * (static_cast<double>(ci) + 0.5);
        data.y = bmin.y() + vox_w * (static_cast<double>(cj) + 0.5);
        data.z = bmin.z() + vox_w * (static_cast<double>(ck) + 0.5);
        data.state = grid.getVoxelState(ci, cj, ck);
        data.num_hits = v.num_hits;
        data.num_beams = v.num_beams;
        data.num_beams_weighted = v.num_beams_weighted;
        data.path_length_raw = v.path_length_raw;
        data.path_length = v.path_length;
        data.free_path_length = v.free_path_length;
        data.effective_free_path_length = v.effective_free_path_length;
        data.num_rays_occluded = v.num_rays_occluded;
        data.path_length_occluded = v.path_length_occluded;
        data.num_unbound_rays = v.num_unbound_rays;
        data.path_length_unbound = v.path_length_unbound;
        data.num_miss_rays = v.num_miss_rays;

        // Classification from the post-traversal table
        auto cit = class_table.find(flat_idx);
        if (cit != class_table.end()) {
            data.classification_hits = cit->second;
            data.dominant_class = getDominantPointClassification(data.classification_hits);
            data.absolute_class  = getAbsolutePointClassification(data.classification_hits);
        }

        data.pad_g0_5      = v.pad_g0_5();
        data.surface_area  = data.pad_g0_5 * voxel_volume;
        data.mean_zenith_angle_rad  = (v.num_beams_weighted > 0) ? (v.sum_of_angles / v.num_beams_weighted) : 0.0;
        if (v.num_beams_weighted > 0) {
          const double mean_sin = v.sum_sin_azimuth / v.num_beams_weighted;
          const double mean_cos = v.sum_cos_azimuth / v.num_beams_weighted;
          double az = std::atan2(mean_sin, mean_cos);
          if (az < 0.0) az += 2.0 * kPi;
          data.mean_azimuth_rad = az;
          data.azimuth_concentration = std::sqrt(mean_sin * mean_sin + mean_cos * mean_cos);
        }
        data.mean_laser_dist = (v.num_beams_weighted > 0) ? (v.sum_of_laser_distances / v.num_beams_weighted) : 0.0;

        if (dtm && dtm->isValid()) {
            double ground_height;
            if (dtm->getHeight(data.x, data.y, ground_height))
                data.distance_from_ground = data.z - ground_height;
        }

        // Leaf/wood hit counts shared by the G=0.5, parametric, and empirical metric blocks.
        float leaf_hits = 0.0f, wood_hits = 0.0f;
        if (params.calc_veg_metrics || params.calc_inclination_dist || params.has_leaf || params.has_wood) {
            for (int code : leaf_classes)
                if (code >= 0 && code <= 255) leaf_hits += data.classification_hits[code];
            for (int code : wood_classes)
                if (code >= 0 && code <= 255) wood_hits += data.classification_hits[code];
        }
        data.num_hit_leaf = static_cast<int32_t>(leaf_hits);
        data.num_hit_wood = static_cast<int32_t>(wood_hits);

        if (v.num_hits > 0) {
            const double hit_total = static_cast<double>(v.num_hits);
            if (params.has_leaf) data.lad_g0_5 = data.pad_g0_5 * (leaf_hits / hit_total);
            if (params.has_wood) data.wad_g0_5 = data.pad_g0_5 * (wood_hits / hit_total);
        }

        if (params.calc_veg_metrics && v.path_length > 0) {
            // G evaluated at mean beam zenith angle — first-order approximation.
            // Use --veg_metrics (IAD active by default) for angle-integrated G_eff via pad/lad/wad.
            double g_theta = computeG(data.mean_zenith_angle_rad, params.lad, lad_param1, lad_param2);
            if (g_theta > 0) {
                double lambda    = computeLambdaV(v);
                double hit_total = std::max(1e-10, static_cast<double>(v.num_hits));
                data.pad_g_corrected = lambda / g_theta;
                data.pad_leaf        = lambda * (leaf_hits / hit_total) / g_theta;
                data.pad_wood        = lambda * (wood_hits / hit_total) / g_theta;
            }
        }

        if (params.calc_inclination_dist) {
          // Source the voxel's G correction from its predominant tree's per-tree IAD. -1 (no
          // tree_id data) falls back to the computeG-at-mean-zenith model below.
          int32_t tid = -1;
          auto pit = predominant_tree.find(flat_idx);
          if (pit != predominant_tree.end()) tid = pit->second;
          data.predominant_tree = tid;

          const IadData* iad = nullptr;
          if (tid >= 0) {
            auto iit = per_tree_iad.find(tid);
            if (iit != per_tree_iad.end()) iad = &iit->second;
          }

          if (iad) {
            data.liad_dewit = iad->liad_dewit;
            data.wiad_dewit = iad->wiad_dewit;
            data.piad_dewit = iad->piad_dewit;
            data.g_leaf  = iad->leaf_g;
            data.g_wood  = iad->wood_g;
            data.g_plant = iad->plant_g;
          }

          if (iad && v.path_length > 0) {
            const double hit_total = std::max(1e-10, static_cast<double>(v.num_hits));
            const double leaf_hits = static_cast<double>(iad->leaf_hits);
            const double wood_hits = static_cast<double>(iad->wood_hits);
            for (const auto& method : params.attenuation_methods) {
              if (method == "bailey") {
                // Bailey (2017) eq.10: per-class G from triangle facets (Eq.4).
                double pad_v = 0.0, lad_v = 0.0, wad_v = 0.0;
                if (iad->bailey_g_leaf > 0 && leaf_hits > 0)
                  lad_v = solveBaileyPadEq10(v.path_length, v.num_beams_weighted,
                                             v.num_hits * (leaf_hits / hit_total), iad->bailey_g_leaf);
                if (iad->bailey_g_wood > 0 && wood_hits > 0)
                  wad_v = solveBaileyPadEq10(v.path_length, v.num_beams_weighted,
                                             v.num_hits * (wood_hits / hit_total), iad->bailey_g_wood);
                if      (leaf_hits == 0 && wood_hits  > 0) pad_v = wad_v;
                else if (leaf_hits  > 0 && wood_hits == 0) pad_v = lad_v;
                else if (leaf_hits  > 0 && wood_hits  > 0 && iad->plant_g > 0)
                  pad_v = solveBaileyPadEq10(v.path_length, v.num_beams_weighted,
                                             v.num_hits, iad->plant_g);
                data.pad_per_method[method] = pad_v;
                data.lad_per_method[method] = lad_v;
                data.wad_per_method[method] = wad_v;
              } else {
                // Vicari et al. (2019) path: angle-integrated G from the tree's empirical PIAD.
                double lambda = computeLambda(v, method);
                if (iad->plant_g > 0) data.pad_per_method[method] = lambda / iad->plant_g;
                if (iad->leaf_g  > 0) data.lad_per_method[method] = lambda * (leaf_hits / hit_total) / iad->leaf_g;
                if (iad->wood_g  > 0) data.wad_per_method[method] = lambda * (wood_hits / hit_total) / iad->wood_g;
              }
            }
          } else if (v.path_length > 0) {
            // Backward-compat fallback (no tree_id): G evaluated at the voxel's mean beam zenith
            // via the existing LAD model, so pad/lad/wad columns are still populated.
            const double g_theta = computeG(data.mean_zenith_angle_rad, params.lad, lad_param1, lad_param2);
            data.g_plant = g_theta;
            if (g_theta > 0) {
              const double hit_total = std::max(1e-10, static_cast<double>(v.num_hits));
              for (const auto& method : params.attenuation_methods) {
                const double lambda = computeLambda(v, method);
                data.pad_per_method[method] = lambda / g_theta;
                if (params.has_leaf) data.lad_per_method[method] = lambda * (leaf_hits / hit_total) / g_theta;
                if (params.has_wood) data.wad_per_method[method] = lambda * (wood_hits / hit_total) / g_theta;
              }
            }
          }
        }

        // Always populate AMAPVox core columns (flag-free).
        // NOTE: effective_free_path_length / bs_effective_free_path are NOT persisted to
        // shards (see raylasbinaryio.h kShardVersion); like free_path_length / bs_free_path
        // they are lost in out-of-core (OOC) mode and read back as 0.
        data.bs_entering    = v.bs_entering;
        data.bs_intercepted = v.bs_intercepted;
        data.bs_potential   = v.bs_potential;
        data.transmittance  = v.transmittance();
        data.bs_free_path    = v.bs_free_path;
        data.lMeanTotal     = (v.num_beams > 0)
                              ? static_cast<double>(v.path_length_raw) / v.num_beams
                              : 0.0;
        data.lMeanFreeTotal = (v.num_beams > 0)
                              ? static_cast<double>(v.free_path_length) / v.num_beams
                              : 0.0;
        data.lMeanEffectiveFreeTotal = (v.num_beams > 0)
                              ? static_cast<double>(v.effective_free_path_length) / v.num_beams
                              : 0.0;
        data.sd_length = (v.num_beams > 1)
                         ? std::sqrt(std::max(0.0, static_cast<double>(v.path_length_sq_raw) / v.num_beams
                                              - data.lMeanTotal * data.lMeanTotal))
                         : 0.0;
        {
            const double fpl = computeLambda(v, "fpl");
            const double fpl_bias = (v.num_beams > 1 && data.lMeanTotal > 0.0)
                                    ? fpl * fpl * data.lMeanTotal / v.num_beams
                                    : 0.0;
            data.attenuation_fpl_biased     = fpl;
            data.attenuation_fpl_correction = fpl_bias;
            data.attenuation_fpl_unbiased   = fpl - fpl_bias;
            data.weighted_fpl               = v.bs_free_path;
            data.weighted_effective_fpl     = v.bs_effective_free_path;
            data.attenuation_ppl            = computeLambda(v, "ppl");
        }

        if (params.subvoxel_split > 0) {
            int total_subvoxels = params.subvoxel_split * params.subvoxel_split * params.subvoxel_split;
            uint64_t bitmap = 0;
            int set_bits = 0;
            for (int i = 0; i < total_subvoxels; i++) {
                if (v.subvoxel_counts[i] >= static_cast<uint8_t>(params.subvoxel_min_beams)) {
                    bitmap |= (1ULL << i);
                    set_bits++;
                }
            }
            data.exploration_rate = (total_subvoxels > 0) ? static_cast<double>(set_bits) / total_subvoxels : 0.0;
            data.subvoxel_bitmap = bitmap;
        }
        return data;
    };

    if (grid.isFlat()) {
        const int64_t total = dims[0] * dims[1] * dims[2];
        const int64_t dimX = dims[0], dimY = dims[1];
        for (int64_t flat_idx = 0; flat_idx < total; ++flat_idx) {
            const VoxelGrid::Voxel& v = grid.voxelAt(flat_idx);
            if (v.num_hits == 0 && v.num_beams_weighted == 0.0f && v.num_rays_occluded == 0.0f) continue;
            const int64_t ci = flat_idx % dimX;
            const int64_t cj = (flat_idx / dimX) % dimY;
            const int64_t ck = flat_idx / (dimX * dimY);
            results[{ci, cj, ck}] = populateData(flat_idx, ci, cj, ck, v);
        }
    } else {
        // Sparse fallback: iterate only occupied voxels
        for (const auto& pair : grid.getSparseVoxels()) {
            const VoxelCoord& coord = pair.first;
            const VoxelGrid::Voxel& v = pair.second;
            const int64_t flat_idx = grid.flatIndex(coord.x, coord.y, coord.z);
            results[coord] = populateData(flat_idx, coord.x, coord.y, coord.z, v);
        }
    }
    return results;
}

// ==================================================================================
// Refactored Writer Implementations
// ==================================================================================

bool writeAmapVoxFile(const std::string& out_name_stub, const VoxelGrid& grid, const MetricResultsMap& metrics,
                        int padding, const Cuboid& user_bounds, const VoxelizationParameters& params, bool filled_only)
{
  ray::VoxelSpace space;

  Eigen::Vector3d user_extent = user_bounds.max_bound_ - user_bounds.min_bound_;
  Eigen::Matrix<int64_t, 3, 1> user_dims = (user_extent / grid.getVoxelWidth()).array().ceil().cast<int64_t>();

  space.header["min_corner"] = format_vec_string(user_bounds.min_bound_);
  space.header["max_corner"] = format_vec_string(user_bounds.max_bound_);
  space.header["split"] = format_vec_string(Eigen::Vector3i(static_cast<int>(user_dims.x()),
                                                            static_cast<int>(user_dims.y()),
                                                            static_cast<int>(user_dims.z())));
  Eigen::Vector3d res(user_extent.x() / static_cast<double>(user_dims.x()),
                      user_extent.y() / static_cast<double>(user_dims.y()),
                      user_extent.z() / static_cast<double>(user_dims.z()));
  space.header["res"] = format_vec_string(res);
  if (params.subvoxel_split > 0) {
    space.header["subvoxel_split"] = std::to_string(params.subvoxel_split);
    space.header["subvoxel_min_beams"] = std::to_string(params.subvoxel_min_beams);
  }

  std::string colnames = "i j k num_hits";
  if (params.has_leaf) colnames += " num_hit_leaf";
  if (params.has_wood) colnames += " num_hit_wood";
  colnames += " num_beams num_miss_rays num_unbound_rays num_beams_occluded";
  colnames += " path_length path_length_weight free_path_length effective_free_path_length"
              " path_length_occluded path_length_unbound";
  colnames += " voxel_size surface_area";
  if (params.subvoxel_split > 0) colnames += " subvoxel_split";
  colnames += " mean_zenith_angle mean_azimuth_angle_deg azimuth_concentration mean_laser_dist";
  if (!params.dtm_file.empty() || params.dtm_from_class >= 0) colnames += " distance_from_ground";
  if (params.calc_beam_metrics) colnames += " transmittance bs_entering bs_intercepted";
  if (params.subvoxel_split > 0) colnames += " exploration_rate subvoxel_bitmap";
  if (params.calc_inclination_dist) {
    colnames += " predominant_tree piad_dewit";
    if (params.has_leaf) colnames += " liad_dewit";
    if (params.has_wood) colnames += " wiad_dewit";
    colnames += " g_plant";
    if (params.has_leaf) colnames += " g_leaf";
    if (params.has_wood) colnames += " g_wood";
    for (const auto& method : params.attenuation_methods) {
      colnames += " pad_" + method;
      if (params.has_leaf) colnames += " lad_" + method;
      if (params.has_wood) colnames += " wad_" + method;
    }
  }
  space.header["colnames"] = colnames;

  auto process_voxel = [&](int64_t i, int64_t j, int64_t k, const VoxelOutputData* data) {
    VoxelGrid::VoxelState state = data ? data->state : VoxelGrid::VoxelState::UNOBSERVED;
    if (filled_only && state != VoxelGrid::VoxelState::FILLED) return;
    if (!params.write_empty_voxels && !filled_only && state == VoxelGrid::VoxelState::UNOBSERVED) return;

    ray::VoxelData v_data;
    v_data.i = i - padding;
    v_data.j = j - padding;
    v_data.k = k - padding;
    v_data.variables.push_back(std::to_string(data ? data->num_hits : 0));
    if (params.has_leaf) v_data.variables.push_back(std::to_string(data ? data->num_hit_leaf : 0));
    if (params.has_wood) v_data.variables.push_back(std::to_string(data ? data->num_hit_wood : 0));
    v_data.variables.push_back(std::to_string(data ? data->num_beams : 0));
    v_data.variables.push_back(std::to_string(data ? data->num_miss_rays : 0));
    v_data.variables.push_back(std::to_string(data ? data->num_unbound_rays : 0));
    v_data.variables.push_back(std::to_string(data ? data->num_rays_occluded : 0));
    v_data.variables.push_back(std::to_string(data ? data->path_length_raw : 0.0f));
    v_data.variables.push_back(std::to_string(data ? data->path_length : 0.0f));
    v_data.variables.push_back(std::to_string(data ? data->free_path_length : 0.0f));
    v_data.variables.push_back(std::to_string(data ? data->effective_free_path_length : 0.0f));
    v_data.variables.push_back(std::to_string(data ? data->path_length_occluded : 0.0f));
    v_data.variables.push_back(std::to_string(data ? data->path_length_unbound : 0.0f));
    v_data.variables.push_back(std::to_string(grid.getVoxelWidth()));
    v_data.variables.push_back(std::to_string(data ? data->surface_area : 0.0));
    if (params.subvoxel_split > 0) v_data.variables.push_back(std::to_string(params.subvoxel_split));
    v_data.variables.push_back(std::to_string(data ? data->mean_zenith_angle_rad * 180.0 / kPi : 0.0));
    v_data.variables.push_back(std::to_string(data ? data->mean_azimuth_rad * 180.0 / kPi : 0.0));
    v_data.variables.push_back(std::to_string(data ? data->azimuth_concentration : 0.0));
    v_data.variables.push_back(std::to_string(data ? data->mean_laser_dist : 0.0));
    if (!params.dtm_file.empty() || params.dtm_from_class >= 0) {
      v_data.variables.push_back(
        (data && data->distance_from_ground != std::numeric_limits<double>::lowest())
          ? std::to_string(data->distance_from_ground) : "nan");
    }
    if (params.calc_beam_metrics) {
      v_data.variables.push_back(std::to_string(data ? data->transmittance : 1.0));
      v_data.variables.push_back(std::to_string(data ? data->bs_entering : 0.0f));
      v_data.variables.push_back(std::to_string(data ? data->bs_intercepted : 0.0f));
    }
    if (params.subvoxel_split > 0) {
      v_data.variables.push_back(std::to_string(data ? data->exploration_rate : 0.0));
      v_data.variables.push_back(std::to_string(data ? data->subvoxel_bitmap : uint64_t(0)));
    }
    if (params.calc_inclination_dist) {
      v_data.variables.push_back(std::to_string(data ? data->predominant_tree : -1));
      v_data.variables.push_back(data ? data->piad_dewit : "");
      if (params.has_leaf) v_data.variables.push_back(data ? data->liad_dewit : "");
      if (params.has_wood) v_data.variables.push_back(data ? data->wiad_dewit : "");
      v_data.variables.push_back(std::to_string(data ? data->g_plant : 0.0));
      if (params.has_leaf) v_data.variables.push_back(std::to_string(data ? data->g_leaf : 0.0));
      if (params.has_wood) v_data.variables.push_back(std::to_string(data ? data->g_wood : 0.0));
      for (const auto& method : params.attenuation_methods) {
        auto lookup = [&](const std::unordered_map<std::string, double>& m) -> double {
          auto it = m.find(method); return it != m.end() ? it->second : 0.0;
        };
        v_data.variables.push_back(std::to_string(data ? lookup(data->pad_per_method) : 0.0));
        if (params.has_leaf) v_data.variables.push_back(std::to_string(data ? lookup(data->lad_per_method) : 0.0));
        if (params.has_wood) v_data.variables.push_back(std::to_string(data ? lookup(data->wad_per_method) : 0.0));
      }
    }
    space.voxels.push_back(v_data);
  };

  if (params.write_empty_voxels && !filled_only) {
    auto dims = grid.getDimensions();
    for (int64_t k = padding; k < dims.z() - padding; ++k) {
    for (int64_t j = padding; j < dims.y() - padding; ++j) {
    for (int64_t i = padding; i < dims.x() - padding; ++i) {
        auto it = metrics.find({i, j, k});
        process_voxel(i, j, k, it != metrics.end() ? &it->second : nullptr);
    }}}
  } else {
    for (const auto& pair : metrics) {
      if ( (pair.first.x < padding || pair.first.x >= grid.getDimensions().x() - padding) ||
           (pair.first.y < padding || pair.first.y >= grid.getDimensions().y() - padding) ||
           (pair.first.z < padding || pair.first.z >= grid.getDimensions().z() - padding) ) continue;
      process_voxel(pair.first.x, pair.first.y, pair.first.z, &pair.second);
    }
  }

  std::string out_name = out_name_stub + ".vox";
  std::cout << "Writing AMAPVox output to " << out_name << "..." << std::endl;
  return ray::writeVox(out_name, space);
}


bool writeTextFile(const std::string& out_name_stub, const VoxelGrid& grid, const MetricResultsMap& metrics,
                     int padding, const Cuboid& user_bounds, const VoxelizationParameters& params, bool filled_only)
{
  std::string filename = out_name_stub + ".txt";
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
    return false;
  }
  outfile << std::fixed << std::setprecision(6);
  std::string header = "i j k x y z num_hits";
  if (params.has_leaf) header += " num_hit_leaf";
  if (params.has_wood) header += " num_hit_wood";
  header += " num_beams num_miss_rays num_unbound_rays num_beams_occluded";
  header += " path_length path_length_weight free_path_length effective_free_path_length"
            " path_length_occluded path_length_unbound";
  header += " voxel_size surface_area";
  if (params.subvoxel_split > 0) header += " subvoxel_split";
  header += " mean_zenith_angle mean_azimuth_angle_deg azimuth_concentration mean_laser_dist";
  if (!params.dtm_file.empty() || params.dtm_from_class >= 0) header += " distance_from_ground";
  if (params.calc_beam_metrics) header += " transmittance bs_entering bs_intercepted";
  if (params.subvoxel_split > 0) header += " exploration_rate subvoxel_bitmap";
  if (params.calc_inclination_dist) {
    header += " predominant_tree piad_dewit";
    if (params.has_leaf) header += " liad_dewit";
    if (params.has_wood) header += " wiad_dewit";
    header += " g_plant";
    if (params.has_leaf) header += " g_leaf";
    if (params.has_wood) header += " g_wood";
    for (const auto& method : params.attenuation_methods) {
      header += " pad_" + method;
      if (params.has_leaf) header += " lad_" + method;
      if (params.has_wood) header += " wad_" + method;
    }
  }
  header += "\n";
  outfile << header;

  long long point_count = 0;

  auto write_line = [&](const VoxelOutputData& data) {
    outfile << (data.i - padding) << " " << (data.j - padding) << " " << (data.k - padding) << " "
            << data.x << " " << data.y << " " << data.z << " "
            << data.num_hits;
    if (params.has_leaf) outfile << " " << data.num_hit_leaf;
    if (params.has_wood) outfile << " " << data.num_hit_wood;
    outfile << " " << data.num_beams << " " << data.num_miss_rays
            << " " << data.num_unbound_rays << " " << data.num_rays_occluded;
    outfile << " " << data.path_length_raw << " " << data.path_length
            << " " << data.free_path_length << " " << data.effective_free_path_length
            << " " << data.path_length_occluded << " " << data.path_length_unbound;
    outfile << " " << grid.getVoxelWidth() << " " << data.surface_area;
    if (params.subvoxel_split > 0) outfile << " " << params.subvoxel_split;
    outfile << " " << (data.mean_zenith_angle_rad * 180.0 / kPi)
            << " " << (data.mean_azimuth_rad * 180.0 / kPi)
            << " " << data.azimuth_concentration << " " << data.mean_laser_dist;
    if (!params.dtm_file.empty() || params.dtm_from_class >= 0) {
      if (data.distance_from_ground != std::numeric_limits<double>::lowest())
        outfile << " " << data.distance_from_ground;
      else
        outfile << " nan";
    }
    if (params.calc_beam_metrics) outfile << " " << data.transmittance << " " << data.bs_entering << " " << data.bs_intercepted;
    if (params.subvoxel_split > 0) outfile << " " << data.exploration_rate << " " << data.subvoxel_bitmap;
    if (params.calc_inclination_dist) {
      outfile << " " << data.predominant_tree << " " << data.piad_dewit;
      if (params.has_leaf) outfile << " " << data.liad_dewit;
      if (params.has_wood) outfile << " " << data.wiad_dewit;
      outfile << " " << data.g_plant;
      if (params.has_leaf) outfile << " " << data.g_leaf;
      if (params.has_wood) outfile << " " << data.g_wood;
      for (const auto& method : params.attenuation_methods) {
        auto lookup = [&](const std::unordered_map<std::string, double>& m) -> double {
          auto it = m.find(method); return it != m.end() ? it->second : 0.0;
        };
        outfile << " " << lookup(data.pad_per_method);
        if (params.has_leaf) outfile << " " << lookup(data.lad_per_method);
        if (params.has_wood) outfile << " " << lookup(data.wad_per_method);
      }
    }
    outfile << "\n";
    point_count++;
  };

  if (params.write_empty_voxels && !filled_only) {
    auto dims = grid.getDimensions();
    for (int64_t k = padding; k < dims.z() - padding; ++k) {
    for (int64_t j = padding; j < dims.y() - padding; ++j) {
    for (int64_t i = padding; i < dims.x() - padding; ++i) {
      auto it = metrics.find({i, j, k});
      if (it != metrics.end()) {
        write_line(it->second);
      } else {
        // UNOBSERVED voxel
        VoxelOutputData empty_data;
        empty_data.i = i; empty_data.j = j; empty_data.k = k;
        empty_data.x = user_bounds.min_bound_.x() + grid.getVoxelWidth() * (static_cast<double>(i - padding) + 0.5);
        empty_data.y = user_bounds.min_bound_.y() + grid.getVoxelWidth() * (static_cast<double>(j - padding) + 0.5);
        empty_data.z = user_bounds.min_bound_.z() + grid.getVoxelWidth() * (static_cast<double>(k - padding) + 0.5);
        write_line(empty_data);
      }
    }}}
  } else {
    for (const auto& pair : metrics) {
      if ( (pair.first.x < padding || pair.first.x >= grid.getDimensions().x() - padding) ||
           (pair.first.y < padding || pair.first.y >= grid.getDimensions().y() - padding) ||
           (pair.first.z < padding || pair.first.z >= grid.getDimensions().z() - padding) ) continue;
      if (filled_only && pair.second.state != VoxelGrid::VoxelState::FILLED) continue;
      write_line(pair.second);
    }
  }

  outfile.close();
  std::cout << "Wrote " << point_count << " voxels to " << filename << std::endl;
  return true;
}

bool writeNetcdfFile(const std::string& out_name_stub, const VoxelGrid& grid, const MetricResultsMap& metrics,
                       int padding, const Cuboid& user_bounds, const VoxelizationParameters& params,
                       bool filled_only)
{
#if RAYLIB_WITH_NETCDF
  try {
    std::string filename = out_name_stub + ".nc";
    netCDF::NcFile dataFile(filename, netCDF::NcFile::replace, netCDF::NcFile::nc4);

    std::vector<VoxelOutputData> data_to_write;
    if (params.write_empty_voxels && !filled_only) {
        auto dims = grid.getDimensions();
        size_t non_padded_count = (dims.x() - 2*padding) * (dims.y() - 2*padding) * (dims.z() - 2*padding);
        data_to_write.reserve(non_padded_count);
        for (int64_t k = padding; k < dims.z() - padding; ++k) {
        for (int64_t j = padding; j < dims.y() - padding; ++j) {
        for (int64_t i = padding; i < dims.x() - padding; ++i) {
          auto it = metrics.find({i, j, k});
          if (it != metrics.end()) {
            data_to_write.push_back(it->second);
          } else {
            VoxelOutputData empty_data;
            empty_data.i = i; empty_data.j = j; empty_data.k = k;
            data_to_write.push_back(empty_data);
          }
        }}}
    } else {
        for (const auto& pair : metrics) {
            if ( (pair.first.x < padding || pair.first.x >= grid.getDimensions().x() - padding) ||
                 (pair.first.y < padding || pair.first.y >= grid.getDimensions().y() - padding) ||
                 (pair.first.z < padding || pair.first.z >= grid.getDimensions().z() - padding) ) continue;
            if (filled_only && pair.second.state != VoxelGrid::VoxelState::FILLED) continue;
            data_to_write.push_back(pair.second);
        }
    }

    long long point_count = data_to_write.size();
    if (point_count == 0) {
        std::cout << "Wrote 0 voxels to " << filename << std::endl;
        return true;
    }
    long long total_classification_hits = 0;
    for(const auto& data : data_to_write) {
        for (int c = 0; c < 256; ++c) if (data.classification_hits[c] > 0.0f) ++total_classification_hits;
    }

    auto nPoints = dataFile.addDim("nPoints", point_count);
    auto nClassificationHits = dataFile.addDim("nClassificationHits", total_classification_hits > 0 ? total_classification_hits : 1);

    std::map<std::string, netCDF::NcVar> vars;
    vars["i"] = dataFile.addVar("i", netCDF::ncInt, {nPoints});
    vars["j"] = dataFile.addVar("j", netCDF::ncInt, {nPoints});
    vars["k"] = dataFile.addVar("k", netCDF::ncInt, {nPoints});
    vars["voxel_state"] = dataFile.addVar("voxel_state", netCDF::ncInt, {nPoints});
    vars["pointclass"] = dataFile.addVar("pointclass", netCDF::ncInt, {nPoints});
    vars["absolute_pointclass"] = dataFile.addVar("absolute_pointclass", netCDF::ncInt, {nPoints});
    vars["num_hits"] = dataFile.addVar("num_hits", netCDF::ncFloat, {nPoints});
    vars["num_beams_weighted"] = dataFile.addVar("num_beams_weighted", netCDF::ncFloat, {nPoints});
    vars["pad_g0_5"] = dataFile.addVar("pad_g0_5", netCDF::ncDouble, {nPoints});
    if (params.has_leaf) vars["lad_g0_5"] = dataFile.addVar("lad_g0_5", netCDF::ncDouble, {nPoints});
    if (params.has_wood) vars["wad_g0_5"] = dataFile.addVar("wad_g0_5", netCDF::ncDouble, {nPoints});
    vars["surface_area"] = dataFile.addVar("surface_area", netCDF::ncDouble, {nPoints});
    vars["mean_zenith_angle_rad"] = dataFile.addVar("mean_zenith_angle_rad", netCDF::ncDouble, {nPoints});
    vars["mean_azimuth_rad"] = dataFile.addVar("mean_azimuth_rad", netCDF::ncDouble, {nPoints});
    vars["azimuth_concentration"] = dataFile.addVar("azimuth_concentration", netCDF::ncDouble, {nPoints});
    vars["mean_laser_dist"] = dataFile.addVar("mean_laser_dist", netCDF::ncDouble, {nPoints});
    if (!params.dtm_file.empty() || params.dtm_from_class >= 0) {
      vars["distance_from_ground"] = dataFile.addVar("distance_from_ground", netCDF::ncDouble, {nPoints});
      vars["distance_from_ground"].putAtt("_FillValue", netCDF::ncDouble, std::numeric_limits<double>::lowest());
    }
    if (params.calc_veg_metrics) {
      vars["pad_g_corrected"] = dataFile.addVar("pad_g_corrected", netCDF::ncDouble, {nPoints});
      vars["pad_leaf"] = dataFile.addVar("pad_leaf", netCDF::ncDouble, {nPoints});
      vars["pad_wood"] = dataFile.addVar("pad_wood", netCDF::ncDouble, {nPoints});
    }
    if (params.calc_beam_metrics) {
      vars["transmittance"] = dataFile.addVar("transmittance", netCDF::ncDouble, {nPoints});
    }
    if (params.subvoxel_split > 0) {
      vars["exploration_rate"] = dataFile.addVar("exploration_rate", netCDF::ncDouble, {nPoints});
      vars["subvoxel_bitmap"] = dataFile.addVar("subvoxel_bitmap", netCDF::ncUint64, {nPoints});
    }
    if (params.calc_inclination_dist) {
      vars["predominant_tree"] = dataFile.addVar("predominant_tree", netCDF::ncInt, {nPoints});
      if (params.has_leaf) vars["liad_dewit"] = dataFile.addVar("liad_dewit", netCDF::ncString, {nPoints});
      if (params.has_wood) vars["wiad_dewit"] = dataFile.addVar("wiad_dewit", netCDF::ncString, {nPoints});
      vars["piad_dewit"] = dataFile.addVar("piad_dewit", netCDF::ncString, {nPoints});
    }

    vars["voxel_id"] = dataFile.addVar("voxel_id", netCDF::ncInt, {nClassificationHits});
    vars["hit_classification_code"] = dataFile.addVar("hit_classification_code", netCDF::ncUbyte, {nClassificationHits});
    vars["hit_classification_count"] = dataFile.addVar("hit_classification_count", netCDF::ncFloat, {nClassificationHits});

    std::vector<int> i_data, j_data, k_data, state_data, pclass_data, abs_pclass_data;
    std::vector<float> hits_data, rays_data;
    std::vector<double> pad_data, lad_g05_data, wad_g05_data, sa_data, angle_data, azimuth_data, concentration_data, dist_data, dfg_data, pad_g_data, pad_leaf_data, pad_wood_data, transm_data, explore_data;
    std::vector<uint64_t> bitmap_data;
    std::vector<int> voxel_id_data;
    std::vector<unsigned char> hit_class_code_data;
    std::vector<float> hit_class_count_data;

    i_data.reserve(point_count); j_data.reserve(point_count); k_data.reserve(point_count);
    state_data.reserve(point_count); pclass_data.reserve(point_count); abs_pclass_data.reserve(point_count);
    hits_data.reserve(point_count); rays_data.reserve(point_count);
    pad_data.reserve(point_count);
    if (params.has_leaf) lad_g05_data.reserve(point_count);
    if (params.has_wood) wad_g05_data.reserve(point_count);
    sa_data.reserve(point_count);
    angle_data.reserve(point_count); azimuth_data.reserve(point_count); concentration_data.reserve(point_count); dist_data.reserve(point_count);
    if (!params.dtm_file.empty() || params.dtm_from_class >= 0) { dfg_data.reserve(point_count); }
    if (params.calc_veg_metrics) {
        pad_g_data.reserve(point_count);
        pad_leaf_data.reserve(point_count);
        pad_wood_data.reserve(point_count);
    }
    if (params.calc_beam_metrics) { transm_data.reserve(point_count); }
    if (params.subvoxel_split > 0) { explore_data.reserve(point_count); bitmap_data.reserve(point_count); }
    if (total_classification_hits > 0) {
        voxel_id_data.reserve(total_classification_hits);
        hit_class_code_data.reserve(total_classification_hits);
        hit_class_count_data.reserve(total_classification_hits);
    }

    std::vector<int> predominant_tree_data;
    std::vector<std::string> liad_dewit_data, wiad_dewit_data, piad_dewit_data;
    std::unordered_map<std::string, std::vector<double>> pad_iad_data, lad_data, wad_data;
    if (params.calc_inclination_dist) {
        predominant_tree_data.reserve(point_count);
        if (params.has_leaf) liad_dewit_data.reserve(point_count);
        if (params.has_wood) wiad_dewit_data.reserve(point_count);
        piad_dewit_data.reserve(point_count);
        for (const auto& method : params.attenuation_methods) {
            pad_iad_data[method].reserve(point_count);
            if (params.has_leaf) lad_data[method].reserve(point_count);
            if (params.has_wood) wad_data[method].reserve(point_count);
        }
    }
    int voxel_idx_counter = 0;
    for (const auto& data : data_to_write) {
        i_data.push_back(data.i - padding);
        j_data.push_back(data.j - padding);
        k_data.push_back(data.k - padding);
        state_data.push_back(static_cast<int>(data.state));
        pclass_data.push_back(data.dominant_class);
        abs_pclass_data.push_back(data.absolute_class);
        hits_data.push_back(data.num_hits);
        rays_data.push_back(data.num_beams_weighted);
        pad_data.push_back(data.pad_g0_5);
        if (params.has_leaf) lad_g05_data.push_back(data.lad_g0_5);
        if (params.has_wood) wad_g05_data.push_back(data.wad_g0_5);
        sa_data.push_back(data.surface_area);
        angle_data.push_back(data.mean_zenith_angle_rad);
        azimuth_data.push_back(data.mean_azimuth_rad);
        concentration_data.push_back(data.azimuth_concentration);
        dist_data.push_back(data.mean_laser_dist);
        if (!params.dtm_file.empty() || params.dtm_from_class >= 0) { dfg_data.push_back(data.distance_from_ground); }
        if (params.calc_veg_metrics) {
            pad_g_data.push_back(data.pad_g_corrected);
            pad_leaf_data.push_back(data.pad_leaf);
            pad_wood_data.push_back(data.pad_wood);
        }
        if (params.calc_beam_metrics) { transm_data.push_back(data.transmittance); }
        if (params.subvoxel_split > 0) { explore_data.push_back(data.exploration_rate); bitmap_data.push_back(data.subvoxel_bitmap); }
        if (params.calc_inclination_dist) {
            predominant_tree_data.push_back(data.predominant_tree);
            if (params.has_leaf) liad_dewit_data.push_back(data.liad_dewit);
            if (params.has_wood) wiad_dewit_data.push_back(data.wiad_dewit);
            piad_dewit_data.push_back(data.piad_dewit);
            for (const auto& method : params.attenuation_methods) {
                auto lookup = [&](const std::unordered_map<std::string, double>& m) -> double {
                    auto it = m.find(method); return it != m.end() ? it->second : 0.0;
                };
                pad_iad_data[method].push_back(lookup(data.pad_per_method));
                if (params.has_leaf) lad_data[method].push_back(lookup(data.lad_per_method));
                if (params.has_wood) wad_data[method].push_back(lookup(data.wad_per_method));
            }
        }

        for (int c = 0; c < 256; ++c) {
            if (data.classification_hits[c] > 0.0f) {
                voxel_id_data.push_back(voxel_idx_counter);
                hit_class_code_data.push_back(static_cast<unsigned char>(c));
                hit_class_count_data.push_back(data.classification_hits[c]);
            }
        }
        voxel_idx_counter++;
    }

    vars["i"].putVar(i_data.data());
    vars["j"].putVar(j_data.data());
    vars["k"].putVar(k_data.data());
    vars["voxel_state"].putVar(state_data.data());
    vars["pointclass"].putVar(pclass_data.data());
    vars["absolute_pointclass"].putVar(abs_pclass_data.data());
    vars["num_hits"].putVar(hits_data.data());
    vars["num_beams_weighted"].putVar(rays_data.data());
    vars["pad_g0_5"].putVar(pad_data.data());
    if (params.has_leaf) vars["lad_g0_5"].putVar(lad_g05_data.data());
    if (params.has_wood) vars["wad_g0_5"].putVar(wad_g05_data.data());
    vars["surface_area"].putVar(sa_data.data());
    vars["mean_zenith_angle_rad"].putVar(angle_data.data());
    vars["mean_azimuth_rad"].putVar(azimuth_data.data());
    vars["azimuth_concentration"].putVar(concentration_data.data());
    vars["mean_laser_dist"].putVar(dist_data.data());

    if (!params.dtm_file.empty() || params.dtm_from_class >= 0) { vars["distance_from_ground"].putVar(dfg_data.data()); }
    if (params.calc_veg_metrics) {
        vars["pad_g_corrected"].putVar(pad_g_data.data());
        vars["pad_leaf"].putVar(pad_leaf_data.data());
        vars["pad_wood"].putVar(pad_wood_data.data());
    }
    if (params.calc_beam_metrics) { vars["transmittance"].putVar(transm_data.data()); }
    if (params.subvoxel_split > 0) { vars["exploration_rate"].putVar(explore_data.data()); vars["subvoxel_bitmap"].putVar(bitmap_data.data()); }

    if (total_classification_hits > 0) {
        vars["voxel_id"].putVar(voxel_id_data.data());
        vars["hit_classification_code"].putVar(hit_class_code_data.data());
        vars["hit_classification_count"].putVar(hit_class_count_data.data());
    }

    if (params.calc_inclination_dist) {
        vars["predominant_tree"].putVar(predominant_tree_data.data());
        auto putStrVar = [&](const std::string& vname, const std::vector<std::string>& sv) {
            std::vector<const char*> ptrs; ptrs.reserve(sv.size());
            for (const auto& s : sv) ptrs.push_back(s.c_str());
            vars[vname].putVar(ptrs.data());
        };
        if (params.has_leaf) putStrVar("liad_dewit", liad_dewit_data);
        if (params.has_wood) putStrVar("wiad_dewit", wiad_dewit_data);
        putStrVar("piad_dewit", piad_dewit_data);
        for (const auto& method : params.attenuation_methods) {
            dataFile.addVar("pad_" + method, netCDF::ncDouble, {nPoints}).putVar(pad_iad_data.at(method).data());
            if (params.has_leaf) dataFile.addVar("lad_" + method, netCDF::ncDouble, {nPoints}).putVar(lad_data.at(method).data());
            if (params.has_wood) dataFile.addVar("wad_" + method, netCDF::ncDouble, {nPoints}).putVar(wad_data.at(method).data());
        }
    }

    std::cout << "Wrote " << point_count << " voxels to " << filename << std::endl;
    return true;

  } catch (const netCDF::exceptions::NcException& e) {
    std::cerr << "NetCDF exception: " << e.what() << std::endl;
    return false;
  }
#else
  RAYLIB_UNUSED(out_name_stub); RAYLIB_UNUSED(grid); RAYLIB_UNUSED(metrics); RAYLIB_UNUSED(padding);
  RAYLIB_UNUSED(user_bounds); RAYLIB_UNUSED(params); RAYLIB_UNUSED(filled_only);
  std::cerr << "Error: NetCDF support is not enabled in this build." << std::endl;
  return false;
#endif
}

// ==================================================================================
// Per-Tree IAD Sidecar CSV
// ==================================================================================

bool writePerTreeIadCsv(const std::string& out_name_stub, const PerTreeIadMap& per_tree_iad,
                        const VoxelizationParameters& params, bool has_stem_id)
{
  const std::string filename = out_name_stub + "_iad.csv";
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
    return false;
  }

  const bool any_bailey = std::any_of(params.attenuation_methods.begin(), params.attenuation_methods.end(),
                                      [](const std::string& m) { return m == "bailey"; });
  const bool any_vicari = std::any_of(params.attenuation_methods.begin(), params.attenuation_methods.end(),
                                      [](const std::string& m) { return m != "bailey"; });
  const int n = params.n_iad_bins;

  // Header. Column-active gating mirrors the voxel writers: liad_* only when has_leaf,
  // wiad_* only when has_wood, piad_* always. bailey-suffixed sets are added when a bailey
  // method is active.
  outfile << "tree_id";
  if (has_stem_id) outfile << ",stem_id";  // representative column; not a per-tree join key
  if (any_vicari) {
    if (params.has_leaf) outfile << ",liad_dewit";
    if (params.has_wood) outfile << ",wiad_dewit";
    outfile << ",piad_dewit";
    if (params.has_leaf) for (int b = 0; b < n; ++b) outfile << ",liad_" << b;
    if (params.has_wood) for (int b = 0; b < n; ++b) outfile << ",wiad_" << b;
    for (int b = 0; b < n; ++b) outfile << ",piad_" << b;
  }
  if (any_bailey) {
    if (params.has_leaf) for (int b = 0; b < n; ++b) outfile << ",liad_bailey_" << b;
    if (params.has_wood) for (int b = 0; b < n; ++b) outfile << ",wiad_bailey_" << b;
    for (int b = 0; b < n; ++b) outfile << ",piad_bailey_" << b;
  }
  outfile << "\n";

  outfile << std::fixed << std::setprecision(6);

  // One row per tree_id, ascending for deterministic output.
  std::vector<int32_t> tids;
  tids.reserve(per_tree_iad.size());
  for (const auto& kv : per_tree_iad) tids.push_back(kv.first);
  std::sort(tids.begin(), tids.end());

  auto write_bins = [&](const std::vector<double>& h) {
    for (int b = 0; b < n; ++b) outfile << "," << (b < static_cast<int>(h.size()) ? h[b] : 0.0);
  };

  for (int32_t tid : tids) {
    const IadData& iad = per_tree_iad.at(tid);
    outfile << tid;
    if (has_stem_id) outfile << ",-1";  // stem_id is not a per-tree aggregation key
    if (any_vicari) {
      if (params.has_leaf) outfile << "," << iad.liad_dewit;
      if (params.has_wood) outfile << "," << iad.wiad_dewit;
      outfile << "," << iad.piad_dewit;
      if (params.has_leaf) write_bins(iad.liad);
      if (params.has_wood) write_bins(iad.wiad);
      write_bins(iad.piad);
    }
    if (any_bailey) {
      if (params.has_leaf) write_bins(iad.liad_bailey);
      if (params.has_wood) write_bins(iad.wiad_bailey);
      write_bins(iad.piad_bailey);
    }
    outfile << "\n";
  }

  outfile.close();
  std::cout << "Wrote " << tids.size() << " trees to " << filename << std::endl;
  return true;
}

} // namespace ray
