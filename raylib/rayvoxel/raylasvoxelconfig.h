// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file defines the VoxelizationParameters struct, which consolidates all
// configuration options for the advanced voxelization process. This is part
// of the "Introduce Parameter Object" refactoring pattern to improve code
// clarity, reduce long parameter lists, and make the system more extensible.

#ifndef RAYLIB_RAYVOXEL_RAYLASVOXELCONFIG_H
#define RAYLIB_RAYVOXEL_RAYLASVOXELCONFIG_H

#include <string>
#include <vector>
#include <Eigen/Dense>

namespace ray
{
  struct VoxelizationParameters
  {
    // --- Input/Grid Parameters ---
    std::string cloud_name;
    double voxel_size = 0.1;
    Eigen::Vector3d grid_bounds_min = Eigen::Vector3d::Zero();
    Eigen::Vector3d grid_bounds_max = Eigen::Vector3d::Zero();
    size_t reserve_size = 0; // Pre-allocation hint for the sparse map

    // --- DTM / Ground Parameters ---
    std::string dtm_file;
    int dtm_from_class = -1;
    double dtm_cell_size = 1.0;
    double dtm_filter_distance = 0.2;

    // --- Core Processing Parameters ---
    std::string weighting_method = "equal";
    bool use_occlusion = false;
    bool apply_flat_top = false;
    int neighbour_prior_min_rays = 0;
    size_t num_threads = 0; // 0 for auto-detection

    // --- Out-of-Core Strategy Parameters ---
    bool use_ooc = false;
    size_t ram_budget_mb = 1024;

    // --- Output Control Parameters ---
    std::string output_format = "amapvox";
    // If true, the primary output will be a dense grid including UNOBSERVED voxels.
    bool write_empty_voxels = false;
    // If true, an additional, separate output file containing ONLY FILLED voxels will be created.
    bool write_filled = false;
    // If true and the primary format is not amapvox, an additional .vox file will also be created.
    bool write_amapvox_also = false;

    // --- Vegetation & Beam Metrics Parameters (Tier 1 & 2) ---
    bool calc_veg_metrics = false;
    std::string leaf_classes_str;
    std::string wood_classes_str;
    std::string lad = "spherical";
    std::string lad_params_str;

    bool calc_inclination_dist = false;
    bool output_iad = false;          // write per-bin liad/wiad/piad columns to output
    bool has_leaf = false;            // --leaf_classes was provided
    bool has_wood = false;            // --wood_classes was provided
    int  n_iad_bins = 18;
    std::vector<std::string> attenuation_methods = {"fpl"};  // each: "fpl" | "ppl" | "transmittance" | "bailey"
    int  knn_normal = 10;

    double triangle_lmax     = 0.05;   // max triangle edge length (metres)
    double iad_tile_size = 3.0;   // XY tile edge (m) for tiled parallel KNN in the IAD pass

    bool calc_beam_metrics = false;
    std::string laser_spec_name;
    Eigen::Vector2d beam_params = {0.0, 0.0};
    int subvoxel_split = 0; // N for an N x N x N subvoxel grid
  };

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASVOXELCONFIG_H
