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

MetricResultsMap calculateOutputMetrics(const VoxelGrid& grid, const VoxelizationParameters& params,
                                         const HeightField* dtm, const ClassTable& class_table,
                                         const IadTable& iad_table)
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
        data.num_rays_observed = v.num_rays_observed;
        data.path_length_observed = v.path_length_observed;
        data.num_rays_occluded = v.num_rays_occluded;
        data.path_length_occluded = v.path_length_occluded;

        // Classification from the post-traversal table
        auto cit = class_table.find(flat_idx);
        if (cit != class_table.end()) {
            data.classification_hits = cit->second;
            data.dominant_class = getDominantPointClassification(data.classification_hits);
            data.absolute_class  = getAbsolutePointClassification(data.classification_hits);
        }

        data.pad_bv_total  = v.pad_bv_total();
        data.surface_area  = data.pad_bv_total * voxel_volume;
        data.mean_zenith_angle_rad  = (v.num_rays_observed > 0) ? (v.sum_of_angles / v.num_rays_observed) : 0.0;
        if (v.num_rays_observed > 0) {
          const double mean_sin = v.sum_sin_azimuth / v.num_rays_observed;
          const double mean_cos = v.sum_cos_azimuth / v.num_rays_observed;
          double az = std::atan2(mean_sin, mean_cos);
          if (az < 0.0) az += 2.0 * kPi;
          data.mean_azimuth_rad = az;
          data.azimuth_concentration = std::sqrt(mean_sin * mean_sin + mean_cos * mean_cos);
        }
        data.mean_laser_dist = (v.num_rays_observed > 0) ? (v.sum_of_laser_distances / v.num_rays_observed) : 0.0;

        if (dtm && dtm->isValid()) {
            double ground_height;
            if (dtm->getHeight(data.x, data.y, ground_height))
                data.distance_from_ground = data.z - ground_height;
        }

        // Leaf/wood hit counts are shared by both the parametric and empirical metric blocks.
        float leaf_hits = 0.0f, wood_hits = 0.0f;
        if (params.calc_veg_metrics || params.calc_inclination_dist) {
            for (int code : leaf_classes)
                if (code >= 0 && code <= 255) leaf_hits += data.classification_hits[code];
            for (int code : wood_classes)
                if (code >= 0 && code <= 255) wood_hits += data.classification_hits[code];
        }

        if (params.calc_veg_metrics && v.path_length_observed > 0) {
            double g_theta = computeG(data.mean_zenith_angle_rad, params.lad, lad_param1, lad_param2);
            if (g_theta > 0) {
                data.pad_g_corrected = v.num_hits / (g_theta * v.path_length_observed);
                data.pad_leaf = leaf_hits / (g_theta * v.path_length_observed);
                data.pad_wood = wood_hits / (g_theta * v.path_length_observed);
            }
        }

        if (params.calc_inclination_dist) {
          auto iit = iad_table.find(flat_idx);
          if (iit != iad_table.end()) {
            const IadData& iad = iit->second;
            data.leaf_g  = iad.leaf_g;
            data.wood_g  = iad.wood_g;
            data.plant_g = iad.plant_g;
            data.liad = iad.liad;
            data.wiad = iad.wiad;
            data.piad = iad.piad;
            if (v.path_length_observed > 0) {
              if (iad.plant_g > 0) data.pad = v.num_hits / (iad.plant_g * v.path_length_observed);
              if (iad.leaf_g  > 0) data.lad = leaf_hits  / (iad.leaf_g  * v.path_length_observed);
              if (iad.wood_g  > 0) data.wad = wood_hits  / (iad.wood_g  * v.path_length_observed);
            }
          }
        }

        if (params.calc_beam_metrics)
            data.transmittance = v.transmittance();

        if (params.subvoxel_split > 0) {
            int set_bits = popcount(v.subvoxel_bitmap);
            int total_subvoxels = params.subvoxel_split * params.subvoxel_split * params.subvoxel_split;
            data.exploration_rate = (total_subvoxels > 0) ? static_cast<double>(set_bits) / total_subvoxels : 0.0;
        }
        return data;
    };

    if (grid.isFlat()) {
        const int64_t total = dims[0] * dims[1] * dims[2];
        const int64_t dimX = dims[0], dimY = dims[1];
        for (int64_t flat_idx = 0; flat_idx < total; ++flat_idx) {
            const VoxelGrid::Voxel& v = grid.voxelAt(flat_idx);
            if (v.num_hits == 0.0f && v.num_rays_observed == 0.0f && v.num_rays_occluded == 0.0f) continue;
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

  std::string colnames = "i j k classification nbEchos nbSampling PadBVTotal lgTotal zenithAngleMean azimuthAngleMean azimuthConcentration distLaser";
  if (params.calc_beam_metrics) {
    colnames += " bsEntering bsIntercepted";
  }
  if (params.calc_inclination_dist) {
    colnames += " leaf_g wood_g plant_g";
    for (int b = 0; b < params.n_iad_bins; ++b) colnames += " liad_" + std::to_string(b);
    for (int b = 0; b < params.n_iad_bins; ++b) colnames += " wiad_" + std::to_string(b);
    for (int b = 0; b < params.n_iad_bins; ++b) colnames += " piad_" + std::to_string(b);
    colnames += " pad lad wad";
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
    v_data.variables.push_back(std::to_string(static_cast<int>(state)));
    v_data.variables.push_back(std::to_string(static_cast<int>(data ? data->num_hits : 0.0f)));
    v_data.variables.push_back(std::to_string(static_cast<int>(data ? data->num_rays_observed : 0.0f)));
    v_data.variables.push_back(std::to_string(data ? data->pad_bv_total : 0.0));
    v_data.variables.push_back(std::to_string(data ? data->path_length_observed : 0.0f));
    v_data.variables.push_back(std::to_string(data ? data->mean_zenith_angle_rad * 180.0 / kPi : 0.0));
    v_data.variables.push_back(std::to_string(data ? data->mean_azimuth_rad * 180.0 / kPi : 0.0));
    v_data.variables.push_back(std::to_string(data ? data->azimuth_concentration : 0.0));
    v_data.variables.push_back(std::to_string(data ? data->mean_laser_dist : 0.0));

    if (params.calc_beam_metrics) {
      const VoxelGrid::Voxel& original_voxel = grid.getVoxel(i, j, k);
      v_data.variables.push_back(std::to_string(original_voxel.bs_entering));
      v_data.variables.push_back(std::to_string(original_voxel.bs_intercepted));
    }
    if (params.calc_inclination_dist) {
      v_data.variables.push_back(std::to_string(data ? data->leaf_g  : 0.0));
      v_data.variables.push_back(std::to_string(data ? data->wood_g  : 0.0));
      v_data.variables.push_back(std::to_string(data ? data->plant_g : 0.0));
      for (int b = 0; b < params.n_iad_bins; ++b)
        v_data.variables.push_back(std::to_string(data && b < static_cast<int>(data->liad.size()) ? data->liad[b] : 0.0));
      for (int b = 0; b < params.n_iad_bins; ++b)
        v_data.variables.push_back(std::to_string(data && b < static_cast<int>(data->wiad.size()) ? data->wiad[b] : 0.0));
      for (int b = 0; b < params.n_iad_bins; ++b)
        v_data.variables.push_back(std::to_string(data && b < static_cast<int>(data->piad.size()) ? data->piad[b] : 0.0));
      v_data.variables.push_back(std::to_string(data ? data->pad : 0.0));
      v_data.variables.push_back(std::to_string(data ? data->lad : 0.0));
      v_data.variables.push_back(std::to_string(data ? data->wad : 0.0));
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
  std::string header = "i j k x y z voxel_state pointclass absolute_pointclass num_hits num_rays_observed path_length_observed "
                       "num_rays_occluded path_length_occluded pad_bv_total surface_area voxel_size "
                       "mean_zenith_angle_rad mean_azimuth_rad azimuth_concentration mean_laser_dist";
  if (!params.dtm_file.empty() || params.dtm_from_class >= 0) { header += " distance_from_ground"; }
  if (params.calc_veg_metrics) header += " pad_g_corrected pad_leaf pad_wood";
  if (params.calc_beam_metrics) header += " transmittance";
  if (params.subvoxel_split > 0) header += " exploration_rate";
  if (params.calc_inclination_dist) {
    header += " leaf_g wood_g plant_g";
    for (int b = 0; b < params.n_iad_bins; ++b) header += " liad_" + std::to_string(b);
    for (int b = 0; b < params.n_iad_bins; ++b) header += " wiad_" + std::to_string(b);
    for (int b = 0; b < params.n_iad_bins; ++b) header += " piad_" + std::to_string(b);
    header += " pad lad wad";
  }
  header += " classification_hits\n";
  outfile << header;

  long long point_count = 0;

  auto write_line = [&](const VoxelOutputData& data) {
    outfile << (data.i - padding) << " " << (data.j - padding) << " " << (data.k - padding) << " "
            << data.x << " " << data.y << " " << data.z << " "
            << static_cast<int>(data.state) << " " << data.dominant_class << " " << data.absolute_class << " "
            << data.num_hits << " " << data.num_rays_observed << " " << data.path_length_observed << " "
            << data.num_rays_occluded << " " << data.path_length_occluded << " "
            << data.pad_bv_total << " " << data.surface_area << " " << grid.getVoxelWidth() << " "
            << data.mean_zenith_angle_rad << " " << data.mean_azimuth_rad << " " << data.azimuth_concentration << " " << data.mean_laser_dist;
    if (!params.dtm_file.empty() || params.dtm_from_class >= 0) {
        if (data.distance_from_ground != std::numeric_limits<double>::lowest()) {
            outfile << " " << data.distance_from_ground;
        } else {
            outfile << " " << "nan"; // Use 'nan' for no-data values
        }
    }
    if (params.calc_veg_metrics) outfile << " " << data.pad_g_corrected << " " << data.pad_leaf << " " << data.pad_wood;
    if (params.calc_beam_metrics) outfile << " " << data.transmittance;
    if (params.subvoxel_split > 0) outfile << " " << data.exploration_rate;
    if (params.calc_inclination_dist) {
      outfile << " " << data.leaf_g << " " << data.wood_g << " " << data.plant_g;
      for (int b = 0; b < params.n_iad_bins; ++b)
        outfile << " " << (b < static_cast<int>(data.liad.size()) ? data.liad[b] : 0.0);
      for (int b = 0; b < params.n_iad_bins; ++b)
        outfile << " " << (b < static_cast<int>(data.wiad.size()) ? data.wiad[b] : 0.0);
      for (int b = 0; b < params.n_iad_bins; ++b)
        outfile << " " << (b < static_cast<int>(data.piad.size()) ? data.piad[b] : 0.0);
      outfile << " " << data.pad << " " << data.lad << " " << data.wad;
    }
    outfile << " " << formatClassificationHits(data.classification_hits) << "\n";
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
                       const IadTable& iad_table, bool filled_only)
{
#if RAYLIB_WITH_NETCDF
  try {
    std::string filename = out_name_stub + ".nc";
    netCDF::NcFile dataFile(filename, netCDF::NcFile::replace);

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
    netCDF::NcDim nIADBins;
    if (params.calc_inclination_dist)
      nIADBins = dataFile.addDim("nIADBins", params.n_iad_bins);

    std::map<std::string, netCDF::NcVar> vars;
    vars["i"] = dataFile.addVar("i", netCDF::ncInt, {nPoints});
    vars["j"] = dataFile.addVar("j", netCDF::ncInt, {nPoints});
    vars["k"] = dataFile.addVar("k", netCDF::ncInt, {nPoints});
    vars["voxel_state"] = dataFile.addVar("voxel_state", netCDF::ncInt, {nPoints});
    vars["pointclass"] = dataFile.addVar("pointclass", netCDF::ncInt, {nPoints});
    vars["absolute_pointclass"] = dataFile.addVar("absolute_pointclass", netCDF::ncInt, {nPoints});
    vars["num_hits"] = dataFile.addVar("num_hits", netCDF::ncFloat, {nPoints});
    vars["num_rays_observed"] = dataFile.addVar("num_rays_observed", netCDF::ncFloat, {nPoints});
    vars["pad_bv_total"] = dataFile.addVar("pad_bv_total", netCDF::ncDouble, {nPoints});
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
    }

    vars["voxel_id"] = dataFile.addVar("voxel_id", netCDF::ncInt, {nClassificationHits});
    vars["hit_classification_code"] = dataFile.addVar("hit_classification_code", netCDF::ncUbyte, {nClassificationHits});
    vars["hit_classification_count"] = dataFile.addVar("hit_classification_count", netCDF::ncFloat, {nClassificationHits});

    std::vector<int> i_data, j_data, k_data, state_data, pclass_data, abs_pclass_data;
    std::vector<float> hits_data, rays_data;
    std::vector<double> pad_data, sa_data, angle_data, azimuth_data, concentration_data, dist_data, dfg_data, pad_g_data, pad_leaf_data, pad_wood_data, transm_data, explore_data;
    std::vector<int> voxel_id_data;
    std::vector<unsigned char> hit_class_code_data;
    std::vector<float> hit_class_count_data;

    i_data.reserve(point_count); j_data.reserve(point_count); k_data.reserve(point_count);
    state_data.reserve(point_count); pclass_data.reserve(point_count); abs_pclass_data.reserve(point_count);
    hits_data.reserve(point_count); rays_data.reserve(point_count);
    pad_data.reserve(point_count); sa_data.reserve(point_count);
    angle_data.reserve(point_count); azimuth_data.reserve(point_count); concentration_data.reserve(point_count); dist_data.reserve(point_count);
    if (!params.dtm_file.empty() || params.dtm_from_class >= 0) { dfg_data.reserve(point_count); }
    if (params.calc_veg_metrics) {
        pad_g_data.reserve(point_count);
        pad_leaf_data.reserve(point_count);
        pad_wood_data.reserve(point_count);
    }
    if (params.calc_beam_metrics) { transm_data.reserve(point_count); }
    if (params.subvoxel_split > 0) { explore_data.reserve(point_count); }
    if (total_classification_hits > 0) {
        voxel_id_data.reserve(total_classification_hits);
        hit_class_code_data.reserve(total_classification_hits);
        hit_class_count_data.reserve(total_classification_hits);
    }

    if (params.calc_inclination_dist) {
        auto binCentresVar = dataFile.addVar("iad_bin_centres_deg", netCDF::ncFloat, {nIADBins});
        std::vector<float> bcd(params.n_iad_bins);
        for (int b = 0; b < params.n_iad_bins; ++b)
            bcd[b] = static_cast<float>((b + 0.5) * 90.0 / params.n_iad_bins);
        binCentresVar.putVar(bcd.data());
    }

    std::vector<float> liad_flat, wiad_flat, piad_flat;
    std::vector<double> leaf_g_data, wood_g_data, plant_g_data, pad_iad_data, lad_data, wad_data;
    if (params.calc_inclination_dist) {
        liad_flat.reserve(metrics.size() * params.n_iad_bins);
        wiad_flat.reserve(metrics.size() * params.n_iad_bins);
        piad_flat.reserve(metrics.size() * params.n_iad_bins);
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
        rays_data.push_back(data.num_rays_observed);
        pad_data.push_back(data.pad_bv_total);
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
        if (params.subvoxel_split > 0) { explore_data.push_back(data.exploration_rate); }
        if (params.calc_inclination_dist) {
            const int64_t flat_idx = grid.flatIndex(data.i, data.j, data.k);
            auto iit = iad_table.find(flat_idx);
            if (iit != iad_table.end()) {
                for (int b = 0; b < params.n_iad_bins; ++b) liad_flat.push_back(static_cast<float>(iit->second.liad[b]));
                for (int b = 0; b < params.n_iad_bins; ++b) wiad_flat.push_back(static_cast<float>(iit->second.wiad[b]));
                for (int b = 0; b < params.n_iad_bins; ++b) piad_flat.push_back(static_cast<float>(iit->second.piad[b]));
            } else {
                liad_flat.insert(liad_flat.end(), params.n_iad_bins, 0.0f);
                wiad_flat.insert(wiad_flat.end(), params.n_iad_bins, 0.0f);
                piad_flat.insert(piad_flat.end(), params.n_iad_bins, 0.0f);
            }
            leaf_g_data.push_back(data.leaf_g);
            wood_g_data.push_back(data.wood_g);
            plant_g_data.push_back(data.plant_g);
            pad_iad_data.push_back(data.pad);
            lad_data.push_back(data.lad);
            wad_data.push_back(data.wad);
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
    vars["num_rays_observed"].putVar(rays_data.data());
    vars["pad_bv_total"].putVar(pad_data.data());
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
    if (params.subvoxel_split > 0) { vars["exploration_rate"].putVar(explore_data.data()); }

    if (total_classification_hits > 0) {
        vars["voxel_id"].putVar(voxel_id_data.data());
        vars["hit_classification_code"].putVar(hit_class_code_data.data());
        vars["hit_classification_count"].putVar(hit_class_count_data.data());
    }

    if (params.calc_inclination_dist) {
        dataFile.addVar("liad",  netCDF::ncFloat,  {nPoints, nIADBins}).putVar(liad_flat.data());
        dataFile.addVar("wiad",  netCDF::ncFloat,  {nPoints, nIADBins}).putVar(wiad_flat.data());
        dataFile.addVar("piad",  netCDF::ncFloat,  {nPoints, nIADBins}).putVar(piad_flat.data());
        dataFile.addVar("leaf_g",  netCDF::ncDouble, {nPoints}).putVar(leaf_g_data.data());
        dataFile.addVar("wood_g",  netCDF::ncDouble, {nPoints}).putVar(wood_g_data.data());
        dataFile.addVar("plant_g", netCDF::ncDouble, {nPoints}).putVar(plant_g_data.data());
        dataFile.addVar("pad",     netCDF::ncDouble, {nPoints}).putVar(pad_iad_data.data());
        dataFile.addVar("lad",     netCDF::ncDouble, {nPoints}).putVar(lad_data.data());
        dataFile.addVar("wad",     netCDF::ncDouble, {nPoints}).putVar(wad_data.data());
    }

    std::cout << "Wrote " << point_count << " voxels to " << filename << std::endl;
    return true;

  } catch (const netCDF::exceptions::NcException& e) {
    std::cerr << "NetCDF exception: " << e.what() << std::endl;
    return false;
  }
#else
  RAYLIB_UNUSED(out_name_stub); RAYLIB_UNUSED(grid); RAYLIB_UNUSED(metrics); RAYLIB_UNUSED(padding);
  RAYLIB_UNUSED(user_bounds); RAYLIB_UNUSED(params); RAYLIB_UNUSED(iad_table); RAYLIB_UNUSED(filled_only);
  std::cerr << "Error: NetCDF support is not enabled in this build." << std::endl;
  return false;
#endif
}

} // namespace ray
