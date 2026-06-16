// Woody volume rasterisation — see raylaswoodvolume.h.
#include "raylib/rayvoxel/raylaswoodvolume.h"

#include <algorithm>
#include <cmath>
#include "raylib/rayutils.h"

namespace ray
{
namespace
{
// Axial sub-steps per voxel edge. The axial sample spacing is voxel_size / kSamplesPerVoxel, so a
// branch is split finely relative to a voxel; finer => more accurate spread at more compute.
constexpr int kSamplesPerVoxel = 3;
// Cap on radial rings / angular samples so huge trunks don't explode the sample count.
constexpr int kMaxRings = 8;
constexpr int kMaxAngular = 16;

// Build an arbitrary orthonormal pair (u, v) perpendicular to unit axis a.
void perpendicularFrame(const Eigen::Vector3d &a, Eigen::Vector3d &u, Eigen::Vector3d &v)
{
  // Pick the world axis least aligned with a to avoid a near-degenerate cross product.
  Eigen::Vector3d helper = (std::abs(a.x()) < 0.9) ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0);
  u = a.cross(helper).normalized();
  v = a.cross(u).normalized();
}

// Accumulate @p volume worth of equal-sub-volume samples for a (possibly degenerate) cylinder of
// axis p0->p1 and radius r into @p out. When length==0 (the spherical root) p0==p1 and the samples
// fill a sphere of radius r instead.
void rasteriseCylinder(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1, double r, double volume,
                       bool is_sphere, const Cuboid &bounds, double voxel_size,
                       const Eigen::Matrix<int64_t, 3, 1> &dims,
                       std::unordered_map<VoxelCoord, double, VoxelCoordHash> &out)
{
  if (volume <= 0.0 || r <= 0.0)
    return;

  Eigen::Vector3d axis = p1 - p0;
  const double length = axis.norm();

  // Sample counts: axial from length, radial/angular from radius — all relative to voxel_size so the
  // sampling resolves voxel-scale structure. A thin branch (r << voxel) collapses to one ring/angle.
  const double step = voxel_size / static_cast<double>(kSamplesPerVoxel);
  int n_axial = is_sphere ? std::max(1, static_cast<int>(std::ceil((2.0 * r) / step)))
                          : std::max(1, static_cast<int>(std::ceil(length / step)));
  int n_rings = std::max(1, std::min(kMaxRings, static_cast<int>(std::ceil(r / step))));
  int n_ang = (n_rings == 1 && r < step) ? 1
                                         : std::max(1, std::min(kMaxAngular, static_cast<int>(std::ceil((2.0 * kPi * r) / step))));

  const size_t n_samples = static_cast<size_t>(n_axial) * static_cast<size_t>(n_rings) * static_cast<size_t>(n_ang);
  const double sub_volume = volume / static_cast<double>(n_samples);

  Eigen::Vector3d a_unit, u, v;
  if (is_sphere || length == 0.0)
  {
    a_unit = Eigen::Vector3d(0, 0, 1);
  }
  else
  {
    a_unit = axis / length;
  }
  perpendicularFrame(a_unit, u, v);

  const Eigen::Vector3d &min_b = bounds.min_bound_;
  const double inv_vs = 1.0 / voxel_size;

  for (int iz = 0; iz < n_axial; ++iz)
  {
    // Centre of each axial slab (cell-centred so the set of samples is symmetric in the cylinder).
    const double t = (static_cast<double>(iz) + 0.5) / static_cast<double>(n_axial);  // in [0,1]
    Eigen::Vector3d centre;
    if (is_sphere)
    {
      // March a point along the axis spanning [-r, r]; radial extent shrinks toward the poles.
      const double z = (t * 2.0 - 1.0) * r;  // [-r, r]
      centre = p0 + a_unit * z;
    }
    else
    {
      centre = p0 + axis * t;
    }

    // For a sphere the in-plane radius at height z is sqrt(r^2 - z^2); for a cylinder it is r.
    double ring_r = r;
    if (is_sphere)
    {
      const double z = (t * 2.0 - 1.0) * r;
      const double rr = r * r - z * z;
      ring_r = rr > 0.0 ? std::sqrt(rr) : 0.0;
    }

    for (int ir = 0; ir < n_rings; ++ir)
    {
      // Equal-area radial positions: r_m = R * sqrt((m+0.5)/n_rings) so each (ring,angle) cell carries
      // the same cross-sectional area, hence every sample carries the same sub-volume.
      const double rad = ring_r * std::sqrt((static_cast<double>(ir) + 0.5) / static_cast<double>(n_rings));
      for (int ia = 0; ia < n_ang; ++ia)
      {
        const double ang = (2.0 * kPi) * (static_cast<double>(ia) + 0.5) / static_cast<double>(n_ang);
        const Eigen::Vector3d p = centre + (u * std::cos(ang) + v * std::sin(ang)) * rad;

        const int64_t cx = static_cast<int64_t>(std::floor((p.x() - min_b.x()) * inv_vs));
        const int64_t cy = static_cast<int64_t>(std::floor((p.y() - min_b.y()) * inv_vs));
        const int64_t cz = static_cast<int64_t>(std::floor((p.z() - min_b.z()) * inv_vs));
        if (cx < 0 || cy < 0 || cz < 0 || cx >= dims[0] || cy >= dims[1] || cz >= dims[2])
          continue;
        out[VoxelCoord{ cx, cy, cz }] += sub_volume;
      }
    }
  }
}
}  // namespace

std::unordered_map<VoxelCoord, double, VoxelCoordHash>
computeWoodVolumePerVoxel(const ForestStructure &forest, const Cuboid &bounds, double voxel_size,
                          const Eigen::Matrix<int64_t, 3, 1> &dims)
{
  std::unordered_map<VoxelCoord, double, VoxelCoordHash> out;

  for (const auto &tree : forest.trees)
  {
    const auto &segments = tree.segments();
    for (size_t i = 0; i < segments.size(); ++i)
    {
      const auto &seg = segments[i];
      if (seg.parent_id == -1)
      {
        // Root segment: treat as a sphere (matches TreeStructure::closestPointOnSegment). Excluded
        // from TreeStructure::volume(), included here for completeness; small relative to the branches.
        const double r = seg.radius;
        const double vol = (4.0 / 3.0) * kPi * r * r * r;
        rasteriseCylinder(seg.tip, seg.tip, r, vol, /*is_sphere=*/true, bounds, voxel_size, dims, out);
        continue;
      }
      const Eigen::Vector3d &p0 = segments[seg.parent_id].tip;
      const Eigen::Vector3d &p1 = seg.tip;
      const double r = seg.radius;
      const double length = (p1 - p0).norm();
      const double vol = kPi * r * r * length;  // same convention as TreeStructure::volume()
      rasteriseCylinder(p0, p1, r, vol, /*is_sphere=*/false, bounds, voxel_size, dims, out);
    }
  }

  return out;
}

}  // namespace ray
