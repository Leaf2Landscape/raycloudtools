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
#include <array>
#include <cstdint>
#include <functional>
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
      // Murmur3-inspired: pack coords into two 64-bit words, then finalize.
      uint64_t h = static_cast<uint64_t>(c.x) * 2654435761ULL
                 ^ static_cast<uint64_t>(c.y) * 805459861ULL
                 ^ static_cast<uint64_t>(c.z) * 3674653429ULL;
      h ^= h >> 33;
      h *= 0xff51afd7ed558ccdULL;
      h ^= h >> 33;
      h *= 0xc4ceb9fe1a85ec53ULL;
      h ^= h >> 33;
      return static_cast<std::size_t>(h);
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
    /// @brief Stores the accumulated traversal metrics for a single voxel cell.
    ///        classification_hits is NOT stored here; it is computed in a
    ///        separate post-traversal pass and held in a ClassTable.
    struct Voxel
    {
      int32_t num_hits = 0;               // Count of echo returns landing in this voxel (nbEchos).
      float num_hits_weighted = 0.0f;     // Weighted hit count: Σ beam_weight per hit endpoint; used for PPL m term.
      int32_t num_beams_observed = 0;     // Count of beams traversing this voxel (nbSampling).
      float num_beams_weighted = 0.0f;    // Weighted beam traversal sum: Σ beam_weight per traversal; used for attenuation.
      float path_length_observed = 0.0f;  // Unweighted path length sum: Σ segment_length across all observed traversals.
      float path_length_weighted = 0.0f;  // Weighted path length sum: Σ (segment_length × beam_weight); used for PAD.
      float num_rays_occluded = 0.0f;     // Sum of occluded rays passing through (unweighted).
      float path_length_occluded = 0.0f;  // Sum of path lengths of occluded rays (unweighted).
      float sum_of_angles = 0.0f;         // Weighted sum of zenith angles of rays passing through.
      float sum_sin_azimuth = 0.0f;       // Weighted sum of sin(azimuth) of rays passing through.
      float sum_cos_azimuth = 0.0f;       // Weighted sum of cos(azimuth) of rays passing through.
      float sum_of_laser_distances = 0.0f;// Weighted sum of distances from sensor to voxel center.
      float bs_entering = 0.0f;           // Weighted sum of entering beam cross-sectional area.
      float bs_intercepted = 0.0f;        // Weighted sum of intercepted beam cross-sectional area.
      float sum_bs_path    = 0.0f;        // FPL: beam-section-weighted clipped path length (all observed).
      float sum_hit_delta  = 0.0f;        // PPL: weight*full_δ for terminal (hit) rays.
      float sum_miss_delta = 0.0f;        // PPL: weight*full_δ for traversing rays.
      float num_unbound_rays = 0.0f;    // weighted count of unbound (miss) rays traversing this voxel
      float path_length_unbound = 0.0f; // weighted sum of clipped path lengths for unbound rays
      uint64_t subvoxel_bitmap = 0;       // Bitmap for tracking subvoxel coverage (up to 4x4x4).

      /// @brief Calculates Plant Area Density (PAD) assuming spherical LAD (G=0.5).
      double pad_g0_5() const;

      /// @brief Calculates transmittance based on beam surface metrics.
      double transmittance() const;

      /// @brief The voxels can be summed element-wise for neighbour priors and merging.
      inline void operator+=(const Voxel &other);

      /// @brief The voxels can be multiplied by a scalar, element-wise.
      inline Voxel operator*(double scale) const;
    };

  public:
    /// ram_budget_bytes: if the flat array would exceed this, fall back to sparse map.
    VoxelGrid(const Cuboid &grid_bounds, double vox_width,
              size_t ram_budget_bytes = 2ULL * 1024 * 1024 * 1024,
              size_t sparse_reservation = 0);

    /// Merge a processor map into the grid (thread-safe via internal mutex).
    void merge(const VoxelProcessor& processor);

    /// Move-assign a processor map into the grid (single-threaded, no lock).
    void take(VoxelProcessor& processor);

    /// Bulk-import a moved processor map; not thread-safe — call only after join().
    void absorbMap(std::unordered_map<VoxelCoord, Voxel, VoxelCoordHash>&& m);

    // --- Storage mode ---
    bool isFlat() const { return !use_sparse_fallback_; }

    // Direct flat-index access (write). Call only when isFlat().
    Voxel& voxelAt(int64_t flat_idx) { return flat_voxels_[flat_idx]; }
    const Voxel& voxelAt(int64_t flat_idx) const { return flat_voxels_[flat_idx]; }

    // Flat index arithmetic: i + j*dimX + k*dimX*dimY
    int64_t flatIndex(int64_t i, int64_t j, int64_t k) const
    {
      return i + j * voxel_dims_[0] + k * voxel_dims_[0] * voxel_dims_[1];
    }

    // --- Accessors ---
    VoxelState getVoxelState(int64_t i, int64_t j, int64_t k) const;
    const Voxel& getVoxel(int64_t i, int64_t j, int64_t k) const;
    const Cuboid& getBounds() const { return bounds_; }
    const Eigen::Matrix<int64_t, 3, 1>& getDimensions() const { return voxel_dims_; }
    double getVoxelWidth() const { return voxel_width_; }
    const std::vector<double>& getPeaks() const { return peaks_; }
    int64_t getIndex(int64_t i, int64_t j, int64_t k) const;

    // Sparse map accessor — valid only when !isFlat().
    const std::unordered_map<VoxelCoord, Voxel, VoxelCoordHash>& getSparseVoxels() const { return sparse_voxels_; }

  private:
    std::unordered_map<VoxelCoord, Voxel, VoxelCoordHash>& getSparseVoxels_internal() { return sparse_voxels_; }

    Cuboid bounds_;
    std::vector<Voxel> flat_voxels_;
    bool use_sparse_fallback_ = true;
    std::unordered_map<VoxelCoord, Voxel, VoxelCoordHash> sparse_voxels_;
    double voxel_width_;
    Eigen::Matrix<int64_t, 3, 1> voxel_dims_;
    std::vector<double> peaks_;
    std::mutex merge_mutex_;
  };

  // ClassTable: flat voxel index → per-classification hit counts.
  // Populated in a separate post-traversal pass; only hit voxels have entries.
  using ClassTable = std::unordered_map<int64_t, std::array<float, 256>>;

  struct IadData {
    std::vector<double> bin_centres;       // radians, size n_iad_bins
    std::vector<double> liad, wiad, piad;  // normalized histograms, size n_iad_bins
    double leaf_g = 0.0, wood_g = 0.0, plant_g = 0.0;
    double bailey_g_leaf = 0.0;   // Bailey eq.(4) area*sin(theta)-weighted mean G for leaf facets
    double bailey_g_wood = 0.0;   // Bailey eq.(4) area*sin(theta)-weighted mean G for wood facets
    float leaf_hits = 0.0f;
    float wood_hits = 0.0f;
  };
  using IadTable = std::unordered_map<int64_t, IadData>;

  // --- Inline Voxel Operator Implementations ---

  inline void VoxelGrid::Voxel::operator+=(const VoxelGrid::Voxel &other)
  {
    num_hits += other.num_hits;
    num_hits_weighted += other.num_hits_weighted;
    num_beams_observed += other.num_beams_observed;
    num_beams_weighted += other.num_beams_weighted;
    path_length_observed += other.path_length_observed;
    path_length_weighted += other.path_length_weighted;
    num_rays_occluded += other.num_rays_occluded;
    path_length_occluded += other.path_length_occluded;
    sum_of_angles += other.sum_of_angles;
    sum_sin_azimuth += other.sum_sin_azimuth;
    sum_cos_azimuth += other.sum_cos_azimuth;
    sum_of_laser_distances += other.sum_of_laser_distances;
    bs_entering += other.bs_entering;
    bs_intercepted += other.bs_intercepted;
    sum_bs_path += other.sum_bs_path;
    sum_hit_delta += other.sum_hit_delta;
    sum_miss_delta += other.sum_miss_delta;
    num_unbound_rays += other.num_unbound_rays;
    path_length_unbound += other.path_length_unbound;
    subvoxel_bitmap |= other.subvoxel_bitmap;
  }

  inline VoxelGrid::Voxel VoxelGrid::Voxel::operator*(double scale) const
  {
    Voxel v;
    // NOTE: scaling int32_t fields by a fractional scale then truncating is only used
    // for neighbour-prior interpolation (non-critical path). Low counts may round to 0.
    v.num_hits = static_cast<int32_t>(num_hits * scale);
    v.num_hits_weighted = static_cast<float>(num_hits_weighted * scale);
    v.num_beams_observed = static_cast<int32_t>(num_beams_observed * scale);
    v.num_beams_weighted = static_cast<float>(num_beams_weighted * scale);
    v.path_length_observed = static_cast<float>(path_length_observed * scale);
    v.path_length_weighted = static_cast<float>(path_length_weighted * scale);
    v.num_rays_occluded = static_cast<float>(num_rays_occluded * scale);
    v.path_length_occluded = static_cast<float>(path_length_occluded * scale);
    v.sum_of_angles = static_cast<float>(sum_of_angles * scale);
    v.sum_sin_azimuth = static_cast<float>(sum_sin_azimuth * scale);
    v.sum_cos_azimuth = static_cast<float>(sum_cos_azimuth * scale);
    v.sum_of_laser_distances = static_cast<float>(sum_of_laser_distances * scale);
    v.bs_entering = static_cast<float>(bs_entering * scale);
    v.bs_intercepted = static_cast<float>(bs_intercepted * scale);
    v.sum_bs_path = static_cast<float>(sum_bs_path * scale);
    v.sum_hit_delta = static_cast<float>(sum_hit_delta * scale);
    v.sum_miss_delta = static_cast<float>(sum_miss_delta * scale);
    v.num_unbound_rays = static_cast<float>(num_unbound_rays * scale);
    v.path_length_unbound = static_cast<float>(path_length_unbound * scale);
    v.subvoxel_bitmap = subvoxel_bitmap;
    return v;
  }

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASVOXELISE_H
