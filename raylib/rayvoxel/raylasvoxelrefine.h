// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file declares post-processing refinement algorithms for VoxelGrid objects.
#ifndef RAYLIB_RAYVOXEL_RAYLASVOXELREFINE_H
#define RAYLIB_RAYVOXEL_RAYLASVOXELREFINE_H

#include "raylib/rayvoxel/raylasvoxelise.h" // For VoxelGrid definition
#include <string>

namespace ray
{
  /// @brief Pre-calculates the highest point in each (x,y) column, required for flat top compensation.
  /// This function modifies the internal state of the provided grid.
  /// @param grid The VoxelGrid object to process.
  /// @param file_name The point cloud file to read point data from.
  void calculatePeaks(VoxelGrid& grid, const std::string& file_name);

  /// @brief Fuses voxels with few rays with their neighbours to create more stable density estimates.
  /// This function modifies the voxel metrics within the provided grid.
  /// @param grid The VoxelGrid object to process.
  /// @param min_rays_for_density The minimum number of observed rays required before borrowing from neighbours.
  void applyNeighbourPriors(VoxelGrid& grid, int min_rays_for_density);

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASVOXELREFINE_H
