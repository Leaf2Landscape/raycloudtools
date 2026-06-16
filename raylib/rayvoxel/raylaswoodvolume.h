// Woody volume rasterisation: project tree branch cylinders (from a trees.txt
// ForestStructure produced by rayextract trees/reconstruct) into the voxel grid
// to obtain per-voxel woody material volume (m^3).
#ifndef RAYLIB_RAYVOXEL_RAYLASWOODVOLUME_H
#define RAYLIB_RAYVOXEL_RAYLASWOODVOLUME_H

#include <unordered_map>
#include <Eigen/Dense>
#include "raylib/raylibconfig.h"
#include "raylib/raycuboid.h"
#include "raylib/rayforeststructure.h"
#include "raylib/rayvoxel/raylasvoxelise.h"

namespace ray
{
/// Rasterise every tree's branch cylinders (and the spherical root segment) into the voxel grid
/// defined by @p bounds / @p voxel_size / @p dims, returning a map from voxel coordinate to the
/// woody volume (m^3) contained in that voxel.
///
/// Each branch (segment i>=1) is the cylinder from segments[parent_id].tip to segments[i].tip with
/// radius segments[i].radius — the same convention used by TreeStructure::volume() (V = pi r^2 L).
/// The volume is distributed by equal-sub-volume 3D sampling (axial x radial x angular), so a thick
/// trunk is spread across the voxels its cross-section actually occupies. Samples landing outside the
/// grid are dropped. The voxel mapping is identical to the traversal mapping used during
/// voxelisation: coord = floor((p - bounds.min_bound_) / voxel_size).
std::unordered_map<VoxelCoord, double, VoxelCoordHash>
computeWoodVolumePerVoxel(const ForestStructure &forest, const Cuboid &bounds, double voxel_size,
                          const Eigen::Matrix<int64_t, 3, 1> &dims);

}  // namespace ray

#endif  // RAYLIB_RAYVOXEL_RAYLASWOODVOLUME_H
