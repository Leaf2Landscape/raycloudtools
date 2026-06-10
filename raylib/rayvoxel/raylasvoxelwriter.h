// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file declares I/O functions for writing a VoxelGrid to various file formats.
#ifndef RAYLIB_RAYVOXEL_RAYLASVOXELWRITER_H
#define RAYLIB_RAYVOXEL_RAYLASVOXELWRITER_H

#include "raylib/rayvoxel/raylasvoxelise.h" // For VoxelGrid definition
#include "raylib/raycuboid.h"
#include "raylib/rayvoxel/raylasvoxelconfig.h" // For VoxelizationParameters
#include "raylib/rayvoxel/raylasheightfield.h" // For HeightField
#include <string>
#include <array>
#include <limits>
#include <unordered_map>

namespace ray
{
  // MODIFIED: This struct now holds all final, computed metrics for a single voxel.
  // It serves as the data transfer object between the metric calculation stage
  // and the various writer functions, ensuring calculations are done only once.
  struct VoxelOutputData {
    // Identity
    int64_t i, j, k;
    double x, y, z;
    VoxelGrid::VoxelState state = VoxelGrid::VoxelState::UNOBSERVED;

    // Base Metrics from VoxelGrid::Voxel
    int dominant_class = -1;
    int absolute_class = -1;
    int32_t num_hits = 0;
    int32_t num_beams = 0;
    float num_beams_weighted = 0.0f;
    float path_length_raw = 0.0f;
    float path_length = 0.0f;
    float num_rays_occluded = 0.0f;
    float path_length_occluded = 0.0f;
    float num_unbound_rays = 0.0f;
    float path_length_unbound = 0.0f;
    int32_t num_miss_rays = 0;
    float num_hit_leaf = 0.0f;
    float num_hit_wood = 0.0f;
    std::array<float, 256> classification_hits{};

    // Calculated Metrics
    double pad_g0_5 = 0.0;
    double lad_g0_5 = 0.0;
    double wad_g0_5 = 0.0;
    double surface_area = 0.0;
    double mean_zenith_angle_rad = 0.0;
    double mean_azimuth_rad = 0.0;
    double azimuth_concentration = 0.0;
    double mean_laser_dist = 0.0;

    // DTM-derived Metric
    double distance_from_ground = std::numeric_limits<double>::lowest();

    // Advanced/Optional Metrics (defaulted to indicate "not calculated" or a neutral value)
    double pad_g_corrected = 0.0;
    double pad_leaf = 0.0;
    double pad_wood = 0.0;
    double transmittance = 1.0;
    float bs_entering = 0.0f;       // raw beam-sample accumulator (AMAPVox bsEntering)
    float bs_intercepted = 0.0f;    // raw beam-sample accumulator (AMAPVox bsIntercepted)
    float bs_free_path = 0.0f;       // beam-area-weighted clipped path (AMAPVox weightedFreepathLength)
    double lMeanTotal = 0.0;        // lgTotal / nbSampling (AMAPVox lMeanTotal)
    double lMeanFreeTotal = 0.0;    // mean free-path length: sum_free_path_weighted / num_beams_weighted
    double lMeanEffectiveFreeTotal = 0.0;  // Stage 3 mean effective free-path: effective_free_path_length / num_beams
    double sd_length = 0.0;         // SD of per-beam path lengths (not yet tracked; always 0)
    double bs_potential = 0.0;      // potential beam cross-section (not yet tracked; always 0)
    double attenuation_fpl_biased = 0.0;
    double attenuation_fpl_correction = 0.0;
    double attenuation_fpl_unbiased = 0.0;
    double weighted_fpl = 0.0;
    double weighted_effective_fpl = 0.0;
    double attenuation_ppl = 0.0;
    double exploration_rate = 0.0;
    uint64_t subvoxel_bitmap = 0;
    double leaf_g  = 0.0;
    double wood_g  = 0.0;
    double plant_g = 0.0;
    std::vector<double> liad;
    std::vector<double> wiad;
    std::vector<double> piad;
    std::unordered_map<std::string, double> pad_per_method;
    std::unordered_map<std::string, double> lad_per_method;
    std::unordered_map<std::string, double> wad_per_method;
  };

  // MODIFIED: A type alias for a map that will store the pre-calculated output data,
  // keyed by voxel coordinates for efficient lookup.
  using MetricResultsMap = std::unordered_map<VoxelCoord, VoxelOutputData, VoxelCoordHash>;

  /// Computes the extinction coefficient λ (m⁻¹) for a voxel using the chosen estimator.
  /// method is one of "fpl", "ppl", "transmittance"; unknown values fall back to FPL simplified.
  double computeLambda(const VoxelGrid::Voxel& v, const std::string& method);

  // MODIFIED: Centralized metric calculator function declaration.
  // This function will iterate over the sparse grid once and compute all
  // required output metrics, populating a MetricResultsMap.
  MetricResultsMap calculateOutputMetrics(const VoxelGrid& grid, const VoxelizationParameters& params,
                                           const HeightField* dtm, const ClassTable& class_table,
                                           const IadTable& iad_table);

  // MODIFIED: All writer function signatures are now refactored to be cleaner.
  // They take the pre-calculated MetricResultsMap and the params object,
  // clearly separating the calculation stage from the I/O stage.
  // The 'filled_only' flag is used to differentiate between the primary output and the --write_filled output.

  /// @brief Writes the grid to the AMAPVox .vox format.
  bool writeAmapVoxFile(const std::string& out_name_stub, const VoxelGrid& grid, const MetricResultsMap& metrics,
                        int padding, const Cuboid& user_bounds, const VoxelizationParameters& params, bool filled_only = false);

  /// @brief Writes the grid to a human-readable, space-delimited text file.
  bool writeTextFile(const std::string& out_name_stub, const VoxelGrid& grid, const MetricResultsMap& metrics,
                     int padding, const Cuboid& user_bounds, const VoxelizationParameters& params, bool filled_only = false);

  /// @brief Writes the grid to a NetCDF file. Requires RAYLIB_WITH_NETCDF.
  bool writeNetcdfFile(const std::string& out_name_stub, const VoxelGrid& grid, const MetricResultsMap& metrics,
                       int padding, const Cuboid& user_bounds, const VoxelizationParameters& params,
                       const IadTable& iad_table, bool filled_only = false);

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASVOXELWRITER_H
