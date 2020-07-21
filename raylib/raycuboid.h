// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYCUBOID_H
#define RAYLIB_RAYCUBOID_H

#include "raylib/raylibconfig.h"
#include "rayutils.h"

namespace ray
{
/// Class for intersection tests on axis-aligned cuboids. These are closed intervals, their boundary is
/// included on intersection.
class RAYLIB_EXPORT Cuboid
{
public:
  Cuboid(){}
  Cuboid(const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound);

  /// if the ray defined by @c start, @c dir and @c depth intersects the cuboid, then 
  /// @c depth is updated with the intersection depth, and the function returns true
  /// @c positive_box true has outward normals and intersects with the near box faces
  /// @c positive_box false has inward normals and intersects with the far box faces
  bool rayIntersects(const Eigen::Vector3d &start, const Eigen::Vector3d &dir, double &depth, bool positive_box) const;
  /// ia negative box has inwards normals, so the ray intersects the back wall, filling in the @c depth argument
  bool intersects(const Eigen::Vector3d &pos) const;
  /// overlap of this cuboid with another, 'kissing' cuboids count as intersection. 
  bool overlaps(const Cuboid &other) const;

  Eigen::Vector3d min_bound_, max_bound_;
};
}  // namespace ray

#endif  // RAYLIB_RAYCUBOID_H