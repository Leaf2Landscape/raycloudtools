// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file defines the VoxelProcessor class, which encapsulates the core
// logic for voxelizing a single LiDAR point and its associated rays.
// By isolating this logic, we can use this class as a "work unit" in
// single-threaded, multi-threaded, and out-of-core processing strategies.

#ifndef RAYLIB_RAYVOXEL_RAYLASVOXELPROCESSOR_H
#define RAYLIB_RAYVOXEL_RAYLASVOXELPROCESSOR_H

#include "raylib/rayvoxel/raylasvoxelise.h" // For Voxel, VoxelCoord, etc.
#include "raylib/rayvoxel/raylasatomic.h"   // For atomic_fadd / atomic_or_u64
#include "raylib/raycuboid.h"
#include "raylib/rayvoxel/raylasheightfield.h" // For HeightField
#include <array>

namespace ray
{
  // A simple struct to pass all necessary data for a single point
  // from the producer thread to the consumer threads.
  struct PointData
  {
    double x, y, z;
    uint8_t return_number;
    uint8_t number_of_returns;
    uint8_t classification;
    Eigen::Vector3d beam_origin;
    double distance_to_sensor;
    double gps_time = 0.0;
    int32_t beam_id = -1;
    uint8_t bound = 1;  ///< 1 = bound (real return), 0 = unbound (miss / floating far end)
  };

  constexpr uint8_t kMaxReturnsPerBeam = 16;

  // A bundle of all returns sharing the same pulse (grouped by gps_time),
  // passed from the producer thread to the consumer threads.
  struct BeamData
  {
    Eigen::Vector3d beam_origin;
    double gps_time;
    std::array<PointData, kMaxReturnsPerBeam> returns;
    uint8_t num_returns = 0;
  };


  class VoxelProcessor
  {
  public:
    // The map type used to accumulate voxel results within a single processor.
    using Map = std::unordered_map<VoxelCoord, VoxelGrid::Voxel, VoxelCoordHash>;

    VoxelProcessor(const Cuboid& grid_bounds, double voxel_width, const std::string& weighting_method,
                   bool use_occlusion_rays, bool use_flat_top, const std::vector<double>* peaks,
                   bool calc_beam_metrics, double beam_diameter, double tan_half_divergence, int subvoxel_split,
                   const HeightField* dtm);

    /// @brief Processes a single point, tracing its rays and accumulating results
    ///        into the processor's internal map.
    void processPoint(const PointData& p);

    /// @brief Processes a whole beam (all returns of a single pulse), tracing the
    ///        consecutive sensor->R0->R1->... segments with corrected per-segment
    ///        weights and accumulating results into the processor's internal map.
    void processBeam(const BeamData& beam);

    /// @brief Writes the contents of the internal map to a sorted binary file (shard).
    ///        This is the core operation for the out-of-core strategy.
    /// @param shard_path The full path to the temporary file to be created.
    /// @return True on success, false on failure.
    bool flushToShard(const std::string& shard_path);

    /// @brief Returns a const reference to the internal map of this processor.
    const Map& getMap() const;

    /// @brief Moves the internal map out of this processor (single-threaded optimization).
    Map&& takeMap() { return std::move(sparse_voxels_); }

    /// @brief Returns the number of voxels currently in the internal map.
    size_t size() const { return sparse_voxels_.size(); }

    /// @brief Clears the internal map, ready for the next chunk of work.
    void clear() { sparse_voxels_.clear(); }

    /// @brief Enables direct flat-array write mode. When set, walkGrid and
    ///        hit-recording write directly into the shared flat array via
    ///        atomic adds instead of accumulating into the per-thread map.
    ///        Pass nullptr to revert to map mode (e.g. for OOC path).
    void setFlatTarget(VoxelGrid::Voxel* arr, int64_t dimX, int64_t dimXY) noexcept
    {
      flat_array_ = arr;
      flat_dim_x_ = dimX;
      flat_dim_xy_ = dimXY;
    }

  private:
    /// @enum RayType
    /// @brief Differentiates between the part of the ray before and after the last hit.
    enum class RayType { OBSERVED, OCCLUDED };

    /// @brief The core Amanatides & Woo voxel traversal algorithm.
    void walkGrid(const Eigen::Vector3d &vox_start, const Eigen::Vector3d &vox_end, RayType type, double weight);

    /// @brief Amanatides & Woo traversal over a local N*N*N subvoxel grid.
    void walkSubGrid(const Eigen::Vector3d& local_start, const Eigen::Vector3d& local_end, int split, uint64_t& bitmap);

    // --- Configuration ---
    const Cuboid& bounds_;
    double voxel_width_;
    const Eigen::Matrix<int64_t, 3, 1> voxel_dims_;
    const std::string& weighting_method_;
    bool use_occlusion_rays_;
    bool use_flat_top_;
    const std::vector<double>* peaks_; // Pointer to the main grid's peaks vector
    bool calc_beam_metrics_;
    double beam_diameter_;
    double tan_half_divergence_;
    int subvoxel_split_;
    int64_t row_stride_;
    const HeightField* dtm_; // Pointer to the DTM for ground clipping

    // --- Per-ray state ---
    Eigen::Vector3d current_ray_vox_start_;
    Eigen::Vector3d current_ray_vox_dir_;
    Eigen::Vector3d current_ray_world_start_;

    // --- Flat-array write mode (optional) ---
    VoxelGrid::Voxel* flat_array_ = nullptr;  // non-null → write direct, skip map
    int64_t flat_dim_x_  = 0;
    int64_t flat_dim_xy_ = 0;

    // --- Local Storage ---
    // Each processor accumulates results into its own private map (map mode only).
    Map sparse_voxels_;
  };

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASVOXELPROCESSOR_H
