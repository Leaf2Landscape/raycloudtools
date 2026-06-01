// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file implements post-processing refinement algorithms for VoxelGrid objects,
// such as canopy peak calculation and spatial smoothing via neighbour priors.

#include "raylib/rayvoxel/raylasvoxelrefine.h"
#include "raylib/rayparse.h"
#include "raylib/raylaz.h"

#include <iostream>
#include <iomanip>
#include <limits>
#include <memory> // For std::unique_ptr

// This is required to access the private 'sparse_voxels_' member for modification.
// It's a common pattern when helper functions need direct access to a class's internals.
#include "raylib/rayvoxel/raylasvoxelise.h"

namespace ray
{

// ==================================================================================
// Post-Processing Refinements Implementation
// ==================================================================================

// --- REFACTORED calculatePeaks ---
void calculatePeaks(VoxelGrid& grid, const std::string& file_name)
{
  // NOTE: This function now exclusively supports LAS/LAZ ray clouds, matching the
  // capabilities of the main voxelisation pipeline.
  std::string ext = getFileNameExtension(file_name);
  if (ext != "las" && ext != "laz") {
    // It's possible for this function to be called with an unsupported file type
    // if a user enables flat-top for a non-LAS file. We should inform them and return.
    std::cerr << "Warning: Flat top compensation (calculatePeaks) only supports .las and .laz files. Skipping peak calculation." << std::endl;
    return;
  }

  const auto& dims = grid.getDimensions();

  // Read the ray cloud in chunks, accumulating the highest point in each (x,y) column.
  size_t num_bounded = 0;
  ray::readLas(file_name,
    [&](std::vector<Eigen::Vector3d>& /*starts*/, std::vector<Eigen::Vector3d>& ends,
        std::vector<double>& /*times*/, std::vector<ray::RGBA>& /*colours*/) {
      for (const auto& end_pos : ends) {
        const Eigen::Vector3d vox_end = (end_pos - grid.getBounds().min_bound_) / grid.getVoxelWidth();
        int64_t ix = static_cast<int64_t>(vox_end.x()), iy = static_cast<int64_t>(vox_end.y());
        if (ix >= 0 && ix < dims[0] && iy >= 0 && iy < dims[1]) {
          // This internal access pattern is safe because the friend declaration grants it.
          int64_t peak_id = ix + iy * dims[0];
          grid.peaks_[peak_id] = std::max(grid.peaks_[peak_id], vox_end.z());
        }
      }
    }, num_bounded, 255.0, nullptr);
}

// Unchanged function
void applyNeighbourPriors(VoxelGrid& grid, int min_rays_for_density)
{
  if (min_rays_for_density <= 0) return;

  double num_hit_voxels = 0.0;
  double num_hit_voxels_unsatisfied = 0.0;

  // MODIFIED for sparse grid: Instead of copying a huge dense vector, we copy
  // only the map of existing voxels. This is memory-efficient.
  auto read_voxels = grid.sparse_voxels_;

  // MODIFIED to iterate over the existing voxels in the read-only copy.
  // This is more efficient than iterating over the entire conceptual grid.
  for (const auto& pair : read_voxels) {
    const VoxelCoord& coord = pair.first;
    const VoxelGrid::Voxel& read_voxel = pair.second;

    // Skip voxels near the padded boundary.
    const auto& dims = grid.getDimensions();
    if (coord.x < 1 || coord.x >= dims[0] - 1 ||
        coord.y < 1 || coord.y >= dims[1] - 1 ||
        coord.z < 1 || coord.z >= dims[2] - 1) {
      continue;
    }

    // This is the voxel in the live grid that we will modify.
    VoxelGrid::Voxel& center_voxel = grid.sparse_voxels_.at(coord);

    // Explicitly cast min_rays_for_density to float to avoid a conversion warning.
    bool was_undersampled = (read_voxel.num_rays_observed < static_cast<float>(min_rays_for_density));
    if (read_voxel.is_filled) {
      num_hit_voxels++;
      if (was_undersampled) {
        num_hit_voxels_unsatisfied++;
      }
    }

    float needed = static_cast<float>(min_rays_for_density) - read_voxel.num_rays_observed;
    if (needed <= 0.0f) continue;

    // --- Safe Neighbor Access ---
    // The previous dense array logic is replaced with safe calls to getVoxel().
    // This correctly handles cases where a neighbor might not exist in the map
    // (getVoxel returns a default empty voxel in that case).
    VoxelGrid::Voxel shell1_sum;
    shell1_sum += grid.getVoxel(coord.x - 1, coord.y, coord.z); shell1_sum += grid.getVoxel(coord.x + 1, coord.y, coord.z);
    shell1_sum += grid.getVoxel(coord.x, coord.y - 1, coord.z); shell1_sum += grid.getVoxel(coord.x, coord.y + 1, coord.z);
    shell1_sum += grid.getVoxel(coord.x, coord.y, coord.z - 1); shell1_sum += grid.getVoxel(coord.x, coord.y, coord.z + 1);
    if (shell1_sum.num_rays_observed > 0) {
      double borrow_ratio = std::min(1.0, static_cast<double>(needed) / shell1_sum.num_rays_observed);
      center_voxel += shell1_sum * borrow_ratio;
      needed -= static_cast<float>(shell1_sum.num_rays_observed * borrow_ratio);
    }
    if (needed <= 0.0f) continue;

    VoxelGrid::Voxel shell2_sum;
    shell2_sum += grid.getVoxel(coord.x-1, coord.y-1, coord.z); shell2_sum += grid.getVoxel(coord.x-1, coord.y+1, coord.z); shell2_sum += grid.getVoxel(coord.x+1, coord.y-1, coord.z); shell2_sum += grid.getVoxel(coord.x+1, coord.y+1, coord.z);
    shell2_sum += grid.getVoxel(coord.x-1, coord.y, coord.z-1); shell2_sum += grid.getVoxel(coord.x-1, coord.y, coord.z+1); shell2_sum += grid.getVoxel(coord.x+1, coord.y, coord.z-1); shell2_sum += grid.getVoxel(coord.x+1, coord.y, coord.z+1);
    shell2_sum += grid.getVoxel(coord.x, coord.y-1, coord.z-1); shell2_sum += grid.getVoxel(coord.x, coord.y-1, coord.z+1); shell2_sum += grid.getVoxel(coord.x, coord.y+1, coord.z-1); shell2_sum += grid.getVoxel(coord.x, coord.y+1, coord.z+1);
    if (shell2_sum.num_rays_observed > 0) {
      double borrow_ratio = std::min(1.0, static_cast<double>(needed) / shell2_sum.num_rays_observed);
      center_voxel += shell2_sum * borrow_ratio;
      needed -= static_cast<float>(shell2_sum.num_rays_observed * borrow_ratio);
    }
    if (needed <= 0.0f) continue;

    VoxelGrid::Voxel shell3_sum;
    shell3_sum += grid.getVoxel(coord.x-1, coord.y-1, coord.z-1); shell3_sum += grid.getVoxel(coord.x-1, coord.y-1, coord.z+1); shell3_sum += grid.getVoxel(coord.x-1, coord.y+1, coord.z-1); shell3_sum += grid.getVoxel(coord.x-1, coord.y+1, coord.z+1);
    shell3_sum += grid.getVoxel(coord.x+1, coord.y-1, coord.z-1); shell3_sum += grid.getVoxel(coord.x+1, coord.y-1, coord.z+1); shell3_sum += grid.getVoxel(coord.x+1, coord.y+1, coord.z-1); shell3_sum += grid.getVoxel(coord.x+1, coord.y+1, coord.z+1);
    if (shell3_sum.num_rays_observed > 0) {
      double borrow_ratio = std::min(1.0, static_cast<double>(needed) / shell3_sum.num_rays_observed);
      center_voxel += shell3_sum * borrow_ratio;
    }
  }

  if (num_hit_voxels > 0) {
    const double percentage = 100.0 * num_hit_voxels_unsatisfied / num_hit_voxels;
    std::cout << "Density calculation: " << std::fixed << std::setprecision(1) << percentage
              << "% of filled voxels had insufficient (<" << min_rays_for_density
              << ") rays." << std::endl;
    if (percentage > 50.0) {
        std::cout << "This is high. Consider using a larger voxel size or a denser cloud for more stable results." << std::endl;
    } else if (percentage < 1.0) {
        std::cout << "This is low. You could potentially use a smaller voxel size for more detail." << std::endl;
    }
  }
}

} // namespace ray
