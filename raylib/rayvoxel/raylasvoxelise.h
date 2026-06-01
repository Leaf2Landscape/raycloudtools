// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
#ifndef RAYLIB_RAYVOXEL_RAYLASVOXELISE_H
#define RAYLIB_RAYVOXEL_RAYLASVOXELISE_H

#include "raylib/raycuboid.h"
#include "raylib/rayutils.h"
#include "raylib/rayvoxel/raylasvoxelconfig.h" // For VoxelizationParameters
#include "raylib/rayvoxel/raylasheightfield.h" // For HeightField
#include <cstdint>
#include <functional>
#include <map>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace ray
{
  // Local alias for the 8-bit unsigned classification code carried per point.
  using U8 = uint8_t;

  // Forward declarations
  class VoxelProcessor;
  class VoxelGrid;
  class ProcessingStrategy;
  class ShardMerger;       // Forward-declare for friendship
  class OutOfCoreStrategy; // Forward-declare for friendship
  void calculatePeaks(VoxelGrid& grid, const std::string& file_name);
  void applyNeighbourPriors(VoxelGrid& grid, int min_rays_for_density);

  /// @brief Main orchestrator function to generate an advanced voxel grid from a point cloud.
  bool RAYLIB_EXPORT generateVoxelGrid(const VoxelizationParameters& params);

  // Represents the 3D integer coordinates of a voxel.
  // This serves as the key for the sparse voxel map.
  struct VoxelCoord {
    int64_t x, y, z;

    // Equality operator required for the hash map
    bool operator==(const VoxelCoord& other) const {
      return x == other.x && y == other.y && z == other.z;
    }
  };

  // Provides the hashing mechanism for VoxelCoord for std::unordered_map.
  struct VoxelCoordHash {
    std::size_t operator()(const VoxelCoord& c) const noexcept {
      // FNV-style mix of three int64_t values
      std::size_t h = 0;
      auto mix = [&](int64_t v) {
        h ^= std::hash<int64_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
      };
      mix(c.x); mix(c.y); mix(c.z);
      return h;
    }
  };

  /// @class VoxelGrid
  /// @brief A data container for a sparse voxel grid.
  class RAYLIB_EXPORT VoxelGrid
  {
    // Grant access to private members for strategies and post-processing.
    friend class InProcessStrategy;
    friend class ShardMerger;       // Grant friendship to the Shard Merger helper class.
    friend class OutOfCoreStrategy; // Grant friendship to the OOC strategy.
    friend void calculatePeaks(VoxelGrid& grid, const std::string& file_name);
    friend void applyNeighbourPriors(VoxelGrid& grid, int min_rays_for_density);

  public:
    /// @enum VoxelState
    /// @brief Classification of a voxel's state based on ray traversal.
    enum class VoxelState
    {
      UNOBSERVED, // Voxel was never traversed by any ray part.
      OCCLUDED,   // Traversed only by "occluded" ray parts (beyond the last return).
      EMPTY,      // Traversed by an "observed" ray part but contains no points.
      FILLED      // Contains one or more LiDAR points.
    };

    /// @struct Voxel
    /// @brief Stores the accumulated metrics for a single voxel cell.
    struct Voxel
    {
      float num_hits = 0.0f;              // Sum of hit events (unweighted).
      float num_rays_observed = 0.0f;     // Sum of observed rays passing through (weighted).
      float path_length_observed = 0.0f;  // Sum of path lengths of observed rays (weighted).
      float num_rays_occluded = 0.0f;     // Sum of occluded rays passing through (unweighted).
      float path_length_occluded = 0.0f;  // Sum of path lengths of occluded rays (unweighted).
      bool is_filled = false;             // True if a point return is located in this voxel.
      std::map<U8, float> classification_hits; // Unweighted sum of hits per classification code.
      float sum_of_angles = 0.0f;         // Weighted sum of zenith angles of rays passing through.
      float sum_of_laser_distances = 0.0f;// Weighted sum of distances from sensor to voxel center for rays.
      float bs_entering = 0.0f;           // Weighted sum of entering beam cross-sectional area.
      float bs_intercepted = 0.0f;        // Weighted sum of intercepted beam cross-sectional area.
      uint64_t subvoxel_bitmap = 0;       // Bitmap for tracking subvoxel coverage (up to 4x4x4).

      /// @brief Calculates Plant Area Density (PAD), similar to AMAPVox's PadBVTotal.
      double pad_bv_total() const;

      /// @brief Calculates transmittance based on beam surface metrics.
      double transmittance() const;

      /// @brief The voxels can be summed element-wise for neighbour priors and merging.
      inline void operator+=(const Voxel &other);

      /// @brief The voxels can be multiplied by a scalar, element-wise.
      inline Voxel operator*(double scale) const;
    };

  public:
    VoxelGrid(const Cuboid &grid_bounds, double vox_width, size_t reservation_size = 0);

    /// @brief Merges the results from a VoxelProcessor into this grid's main map.
    ///        This operation is thread-safe.
    void merge(const VoxelProcessor& processor);

    /// @brief Moves the results from a VoxelProcessor into this grid's main map. Not thread-safe.
    void take(VoxelProcessor& processor);

    /// Bulk-import a moved processor map; not thread-safe — call only after join().
    /// The argument type is identical to VoxelProcessor::Map, spelled out here
    /// because VoxelProcessor is only forward-declared in this header.
    void absorbMap(std::unordered_map<VoxelCoord, Voxel, VoxelCoordHash>&& m);

    // --- Accessors ---
    VoxelState getVoxelState(int64_t i, int64_t j, int64_t k) const;
    const Voxel& getVoxel(int64_t i, int64_t j, int64_t k) const;
    const Cuboid& getBounds() const { return bounds_; }
    const Eigen::Matrix<int64_t, 3, 1>& getDimensions() const { return voxel_dims_; }
    double getVoxelWidth() const { return voxel_width_; }
    const std::vector<double>& getPeaks() const { return peaks_; }
    int64_t getIndex(int64_t i, int64_t j, int64_t k) const;

    // Public getter for read-only access to the sparse voxel map.
    const std::unordered_map<VoxelCoord, Voxel, VoxelCoordHash>& getSparseVoxels() const { return sparse_voxels_; }

  private:
    // Private getter for write access, intended only for friend classes.
    std::unordered_map<VoxelCoord, Voxel, VoxelCoordHash>& getSparseVoxels_internal() { return sparse_voxels_; }

    Cuboid bounds_;
    std::unordered_map<VoxelCoord, Voxel, VoxelCoordHash> sparse_voxels_;
    double voxel_width_;
    Eigen::Matrix<int64_t, 3, 1> voxel_dims_;
    std::vector<double> peaks_;

    // Mutex to protect the main map during concurrent merge operations.
    std::mutex merge_mutex_;
  };

  // --- Inline Voxel Operator Implementations ---

  inline void VoxelGrid::Voxel::operator+=(const VoxelGrid::Voxel &other)
  {
    num_hits += other.num_hits;
    num_rays_observed += other.num_rays_observed;
    path_length_observed += other.path_length_observed;
    num_rays_occluded += other.num_rays_occluded;
    path_length_occluded += other.path_length_occluded;
    sum_of_angles += other.sum_of_angles;
    sum_of_laser_distances += other.sum_of_laser_distances;
    bs_entering += other.bs_entering;
    bs_intercepted += other.bs_intercepted;
    is_filled = is_filled || other.is_filled;
    subvoxel_bitmap |= other.subvoxel_bitmap;

    for (const auto& pair : other.classification_hits) {
        classification_hits[pair.first] += pair.second;
    }
  }

  inline VoxelGrid::Voxel VoxelGrid::Voxel::operator*(double scale) const
  {
    Voxel v;
    v.num_hits = static_cast<float>(num_hits * scale);
    v.num_rays_observed = static_cast<float>(num_rays_observed * scale);
    v.path_length_observed = static_cast<float>(path_length_observed * scale);
    v.num_rays_occluded = static_cast<float>(num_rays_occluded * scale);
    v.path_length_occluded = static_cast<float>(path_length_occluded * scale);
    v.sum_of_angles = static_cast<float>(sum_of_angles * scale);
    v.sum_of_laser_distances = static_cast<float>(sum_of_laser_distances * scale);
    v.bs_entering = static_cast<float>(bs_entering * scale);
    v.bs_intercepted = static_cast<float>(bs_intercepted * scale);
    v.is_filled = is_filled;
    v.subvoxel_bitmap = subvoxel_bitmap;

    for (const auto& pair : classification_hits) {
        v.classification_hits[pair.first] = static_cast<float>(pair.second * scale);
    }
    return v;
  }

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASVOXELISE_H
