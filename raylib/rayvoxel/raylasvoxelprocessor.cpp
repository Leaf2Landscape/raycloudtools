// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file implements the VoxelProcessor class. The logic here was
// carefully extracted from the original VoxelGrid class to create a
// modular and reusable processing engine.

#include "raylib/rayvoxel/raylasvoxelprocessor.h"
#include "raylib/rayvoxel/raylasbinaryio.h" // Needed for out-of-core flushing
#include "raylib/rayutils.h" // For kPi, clamped
#include "raylib/rayunused.h"
#include <cassert> // For assert
#include <algorithm> // For std::sort

namespace ray
{

VoxelProcessor::VoxelProcessor(const Cuboid& grid_bounds, double voxel_width, const std::string& weighting_method,
                               bool use_occlusion_rays, bool use_flat_top, const std::vector<double>* peaks,
                               bool calc_beam_metrics, double beam_diameter, double tan_half_divergence, int subvoxel_split,
                               const HeightField* dtm)
  : bounds_(grid_bounds),
    voxel_width_(voxel_width),
    voxel_dims_(((grid_bounds.max_bound_ - grid_bounds.min_bound_) / voxel_width).array().ceil().cast<int64_t>()),
    weighting_method_(weighting_method),
    use_occlusion_rays_(use_occlusion_rays),
    use_flat_top_(use_flat_top),
    peaks_(peaks),
    calc_beam_metrics_(calc_beam_metrics),
    beam_diameter_(beam_diameter),
    tan_half_divergence_(tan_half_divergence),
    subvoxel_split_(subvoxel_split),
    row_stride_(voxel_dims_[0]),
    dtm_(dtm)
{
}

const VoxelProcessor::Map& VoxelProcessor::getMap() const
{
  return sparse_voxels_;
}

bool VoxelProcessor::flushToShard(const std::string& shard_path)
{
  if (sparse_voxels_.empty()) {
    return true; // Nothing to flush
  }

  // 1. Copy map contents to a vector for sorting.
  // This is necessary because hash maps are unordered.
  std::vector<std::pair<VoxelCoord, VoxelGrid::Voxel>> sorted_voxels;
  sorted_voxels.reserve(sparse_voxels_.size());
  for (const auto& pair : sparse_voxels_) {
    sorted_voxels.push_back(pair);
  }

  // 2. Sort the vector by VoxelCoord. This is critical for the external merge-sort algorithm.
  // The sort order (Z, then Y, then X) is chosen for efficient merging.
  std::sort(sorted_voxels.begin(), sorted_voxels.end(),
    [](const auto& a, const auto& b) {
      if (a.first.z != b.first.z) return a.first.z < b.first.z;
      if (a.first.y != b.first.y) return a.first.y < b.first.y;
      return a.first.x < b.first.x;
    }
  );

  // 3. Write the sorted vector to the binary shard file.
  std::ofstream out(shard_path, std::ios::binary);
  if (!out.is_open()) {
    std::cerr << "Error: Could not open temporary shard file for writing: " << shard_path << std::endl;
    return false;
  }
  if (!writeShardHeader(out)) {
    std::cerr << "Error: Failed to write shard header to: " << shard_path << std::endl;
    return false;
  }

  for (const auto& pair : sorted_voxels) {
    if (!writeVoxelData(out, pair.first, pair.second)) {
      std::cerr << "Error: Failed to write to shard file: " << shard_path << std::endl;
      out.close(); // Attempt to close before returning
      return false;
    }
  }

  return true;
}


void VoxelProcessor::processBeam(const BeamData& beam)
{
  if (beam.num_returns == 0) return;

  std::vector<const PointData*> sorted;
  for (uint8_t r = 0; r < beam.num_returns; ++r) sorted.push_back(&beam.returns[r]);
  std::sort(sorted.begin(), sorted.end(),
    [](const PointData* a, const PointData* b){ return a->return_number < b->return_number; });

  const int N = static_cast<int>(beam.num_returns);
  const PointData& farthest = *sorted[N - 1];
  Eigen::Vector3d farthest_pos(farthest.x, farthest.y, farthest.z);

  // Walk the full ray from beam origin to the farthest return exactly once.
  // Weight = 1/N (equal) models each return consuming an equal energy fraction;
  // weight = 1.0 (full) treats the full pulse as a single unattenuated traversal.
  // This matches DensityGrid's one-walk-per-pulse model and eliminates the
  // per-segment ray-count inflation of the old segmented approach.
  // If farthest.bound == 0 (unbound/miss ray), traversal still sweeps through those
  // voxels — correctly marking them as observed/free — but the hit-recording loop
  // below gates on p.bound and will not count the endpoint as a hit.
  double beam_weight = 1.0;
  if (N > 1) {
    beam_weight = (weighting_method_ == "equal") ? 1.0 / N : 1.0;
  }

  Eigen::Vector3d cs = beam.beam_origin, ce = farthest_pos;
  if (bounds_.clipRay(cs, ce, 1e-10)) {
    Eigen::Vector3d vs = (cs - bounds_.min_bound_) / voxel_width_;
    Eigen::Vector3d ve = (ce - bounds_.min_bound_) / voxel_width_;
    current_ray_vox_start_   = vs;
    current_ray_vox_dir_     = (ve - vs).normalized();
    current_ray_world_start_ = beam.beam_origin;
    current_ray_unbound_ = (farthest.bound == 0 and N <= 1);
    walkGrid(vs, ve, RayType::OBSERVED, beam_weight);
  }

  // Pre-compute unit ray direction for sum_bs_path free-path corrections below.
  // Only needed when beam-section metrics are active.
  Eigen::Vector3d ray_dir_unit = Eigen::Vector3d::Zero();
  bool ray_dir_valid = false;
  if (calc_beam_metrics_) {
    const Eigen::Vector3d raw_dir = farthest_pos - beam.beam_origin;
    const double raw_len = raw_dir.norm();
    if (raw_len > 1e-12) { ray_dir_unit = raw_dir / raw_len; ray_dir_valid = true; }
  }

  // Record hits for all returns (num_hits only; traversal already counted above).
  for (int i = 0; i < N; ++i) {
    const PointData& p = *sorted[i];
    if (p.bound == 0) continue;  // unbound (miss) ray: traversed as observed, never a hit
    Eigen::Vector3d curr(p.x, p.y, p.z);
    Eigen::Vector3d vox_coord_filled = (curr - bounds_.min_bound_) / voxel_width_;
    int64_t ix = static_cast<int64_t>(std::floor(vox_coord_filled.x()));
    int64_t iy = static_cast<int64_t>(std::floor(vox_coord_filled.y()));
    int64_t iz = static_cast<int64_t>(std::floor(vox_coord_filled.z()));
    if (ix >= 0 && ix < voxel_dims_[0] && iy >= 0 && iy < voxel_dims_[1] && iz >= 0 && iz < voxel_dims_[2]) {
      if (flat_array_) {
        VoxelGrid::Voxel& v = flat_array_[ix + iy * flat_dim_x_ + iz * flat_dim_xy_];
        atomic_iadd(v.num_hits, 1);
        atomic_fadd(v.num_hits_weighted, static_cast<float>(beam_weight));
        if (calc_beam_metrics_) {
          double dist = p.distance_to_sensor;
          double r = tan_half_divergence_ * dist + 0.5 * beam_diameter_;
          atomic_fadd(v.bs_intercepted, static_cast<float>(kPi * r * r * beam_weight));
          if (ray_dir_valid) {
            // sum_bs_path correction: walkGrid added beam_area × weight × full_transit;
            // FPL needs beam_area × weight × free_path (entry→hit, not entry→exit).
            // Subtract the excess: beam_area × weight × (t_exit − t_hit).
            const Eigen::Vector3d vox_min_w = bounds_.min_bound_ + Eigen::Vector3d(static_cast<double>(ix), static_cast<double>(iy), static_cast<double>(iz)) * voxel_width_;
            double t_entry = 0.0, t_exit = 1e30;
            for (int d = 0; d < 3; ++d) {
              if (std::abs(ray_dir_unit[d]) > 1e-15) {
                const double t1 = (vox_min_w[d]                - beam.beam_origin[d]) / ray_dir_unit[d];
                const double t2 = (vox_min_w[d] + voxel_width_ - beam.beam_origin[d]) / ray_dir_unit[d];
                t_entry = std::max(t_entry, std::min(t1, t2));
                t_exit  = std::min(t_exit,  std::max(t1, t2));
              }
            }
            t_entry = std::max(0.0, t_entry);
            const double full_transit = std::max(0.0, t_exit - t_entry);
            const double excess = std::max(0.0, std::min(full_transit, t_exit - p.distance_to_sensor));
            if (excess > 1e-12) {
              const Eigen::Vector3d vox_center = vox_min_w + Eigen::Vector3d::Constant(0.5 * voxel_width_);
              const double dc = (vox_center - beam.beam_origin).norm();
              const double rc = tan_half_divergence_ * dc + 0.5 * beam_diameter_;
              atomic_fadd(v.sum_bs_path, -static_cast<float>(kPi * rc * rc * beam_weight * excess));
            }
          }
        }
      } else {
        VoxelCoord coord = {ix, iy, iz};
        VoxelGrid::Voxel& v = sparse_voxels_[coord];
        v.num_hits += 1;
        v.num_hits_weighted += static_cast<float>(beam_weight);
        if (calc_beam_metrics_) {
          double dist = p.distance_to_sensor;
          double r = tan_half_divergence_ * dist + 0.5 * beam_diameter_;
          v.bs_intercepted += static_cast<float>(kPi * r * r * beam_weight);
          if (ray_dir_valid) {
            const Eigen::Vector3d vox_min_w = bounds_.min_bound_ + Eigen::Vector3d(static_cast<double>(ix), static_cast<double>(iy), static_cast<double>(iz)) * voxel_width_;
            double t_entry = 0.0, t_exit = 1e30;
            for (int d = 0; d < 3; ++d) {
              if (std::abs(ray_dir_unit[d]) > 1e-15) {
                const double t1 = (vox_min_w[d]                - beam.beam_origin[d]) / ray_dir_unit[d];
                const double t2 = (vox_min_w[d] + voxel_width_ - beam.beam_origin[d]) / ray_dir_unit[d];
                t_entry = std::max(t_entry, std::min(t1, t2));
                t_exit  = std::min(t_exit,  std::max(t1, t2));
              }
            }
            t_entry = std::max(0.0, t_entry);
            const double full_transit = std::max(0.0, t_exit - t_entry);
            const double excess = std::max(0.0, std::min(full_transit, t_exit - p.distance_to_sensor));
            if (excess > 1e-12) {
              const Eigen::Vector3d vox_center = vox_min_w + Eigen::Vector3d::Constant(0.5 * voxel_width_);
              const double dc = (vox_center - beam.beam_origin).norm();
              const double rc = tan_half_divergence_ * dc + 0.5 * beam_diameter_;
              v.sum_bs_path -= static_cast<float>(kPi * rc * rc * beam_weight * excess);
            }
          }
        }
      }
    }
  }

  if (use_occlusion_rays_ && farthest.bound == 1) {
    Eigen::Vector3d direction = (farthest_pos - beam.beam_origin).normalized();
    double large_distance = (bounds_.max_bound_ - bounds_.min_bound_).norm() * 2.0;
    Eigen::Vector3d start_occ = farthest_pos;
    Eigen::Vector3d end_occ = start_occ + direction * large_distance;
    if (bounds_.clipRay(start_occ, end_occ, 1e-10)) {
      Eigen::Vector3d vs = (start_occ - bounds_.min_bound_) / voxel_width_;
      Eigen::Vector3d ve = (end_occ   - bounds_.min_bound_) / voxel_width_;
      current_ray_vox_start_   = vs;
      current_ray_vox_dir_     = (ve - vs).normalized();
      current_ray_world_start_ = farthest_pos;
      walkGrid(vs, ve, RayType::OCCLUDED, 1.0);
    }
  }
}


void VoxelProcessor::walkSubGrid(const Eigen::Vector3d& local_start, const Eigen::Vector3d& local_end, int split, uint64_t& bitmap)
{
    auto sub_walk_lambda =
        [&](const Eigen::Vector3i &p, const Eigen::Vector3i &/*target*/, double, double, double) -> bool {
            if (p.x() >= 0 && p.x() < split && p.y() >= 0 && p.y() < split && p.z() >= 0 && p.z() < split) {
                int sub_idx = p.x() + p.y() * split + p.z() * split * split;
                bitmap |= (1ULL << sub_idx);
            }
            return false; // Continue walking
        };

    ray::walkGrid(local_start, local_end, sub_walk_lambda);
}

void VoxelProcessor::walkGrid(const Eigen::Vector3d &vox_start, const Eigen::Vector3d &vox_end, RayType type, double weight)
{
    // Pre-calculate angle for this ray, as it's constant throughout traversal.
    // Zenith angle is the angle between the ray direction and the vertical Z-axis (0,0,1).
    // The cosine of this angle is simply the z-component of the normalized direction vector.
    double zenith_angle = acos(clamped(current_ray_vox_dir_.z(), -1.0, 1.0));
    double azimuth_angle = std::atan2(current_ray_vox_dir_.x(), current_ray_vox_dir_.y());
    if (azimuth_angle < 0.0) azimuth_angle += 2.0 * kPi;
    const double sin_az = std::sin(azimuth_angle);
    const double cos_az = std::cos(azimuth_angle);

    auto walk_lambda =
        [&](const Eigen::Vector3i &p, const Eigen::Vector3i &/*target*/, double in_length, double out_length, double max_length) -> bool {

            if (p.x() < 0 || p.x() >= voxel_dims_[0] ||
                p.y() < 0 || p.y() >= voxel_dims_[1] ||
                p.z() < 0 || p.z() >= voxel_dims_[2]) {
                return false;
            }

            double end_length = std::min(out_length, max_length);

            // Added an assertion to make the contract explicit: if flat-top
            // compensation is enabled, the peaks_ pointer must be valid.
            if (use_flat_top_ && type == RayType::OBSERVED) {
                assert(peaks_ != nullptr && "If use_flat_top_ is true, peaks_ must be valid.");
                int64_t peak_id = p.x() + p.y() * row_stride_;
                if (peak_id >= 0 && peak_id < static_cast<int64_t>(peaks_->size())) {
                    double peak = (*peaks_)[peak_id];
                    double in_height = current_ray_vox_start_.z() + current_ray_vox_dir_.z() * in_length;
                    double end_height = current_ray_vox_start_.z() + current_ray_vox_dir_.z() * end_length;
                    if (current_ray_vox_dir_.z() < 0.0 && in_height > peak && end_height <= peak) {
                        double t = (in_height - peak) / (in_height - end_height);
                        in_length += (end_length - in_length) * std::max(0.0, std::min(t, 0.99));
                    }
                }
            }

            double length_in_voxel = (end_length - in_length) * voxel_width_;

            if (type == RayType::OCCLUDED && dtm_ && dtm_->isValid()) {
                Eigen::Vector3d voxel_center_world = bounds_.min_bound_ + (p.cast<double>() + Eigen::Vector3d(0.5, 0.5, 0.5)) * voxel_width_;
                double ground_height;
                if (dtm_->getHeightNearest(voxel_center_world.x(), voxel_center_world.y(), ground_height)) {
                    if (voxel_center_world.z() < ground_height) {
                        return false;
                    }
                }
            }

            if (flat_array_) {
                // Direct atomic writes into the shared flat array — no per-thread map.
                VoxelGrid::Voxel& v = flat_array_[p.x() + p.y() * flat_dim_x_ + p.z() * flat_dim_xy_];
                if (type == RayType::OBSERVED) {
                    atomic_iadd(v.num_beams_observed, 1);
                    atomic_fadd(v.num_beams_weighted, static_cast<float>(weight));
                    atomic_fadd(v.path_length_observed, static_cast<float>(length_in_voxel));
                    atomic_fadd(v.path_length_weighted, static_cast<float>(length_in_voxel * weight));
                    atomic_fadd(v.sum_of_angles, static_cast<float>(zenith_angle * weight));
                    atomic_fadd(v.sum_sin_azimuth, static_cast<float>(sin_az * weight));
                    atomic_fadd(v.sum_cos_azimuth, static_cast<float>(cos_az * weight));

                    Eigen::Vector3d voxel_center_world = bounds_.min_bound_ + (p.cast<double>() + Eigen::Vector3d(0.5, 0.5, 0.5)) * voxel_width_;
                    double dist_to_center = (voxel_center_world - current_ray_world_start_).norm();
                    atomic_fadd(v.sum_of_laser_distances, static_cast<float>(dist_to_center * weight));

                    if (calc_beam_metrics_) {
                        double beam_radius = tan_half_divergence_ * dist_to_center + 0.5 * beam_diameter_;
                        atomic_fadd(v.bs_entering, static_cast<float>(kPi * beam_radius * beam_radius * weight));
                        atomic_fadd(v.sum_bs_path, static_cast<float>(kPi * beam_radius * beam_radius * weight * length_in_voxel));
                    }
                    if (subvoxel_split_ > 0) {
                        Eigen::Vector3d ls = (current_ray_vox_start_ + current_ray_vox_dir_ * in_length  - p.cast<double>()) * subvoxel_split_;
                        Eigen::Vector3d le = (current_ray_vox_start_ + current_ray_vox_dir_ * end_length - p.cast<double>()) * subvoxel_split_;
                        uint64_t bits = 0;
                        walkSubGrid(ls, le, subvoxel_split_, bits);
                        if (bits) atomic_or_u64(v.subvoxel_bitmap, bits);
                    }
                    {
                        double full_delta = (out_length - in_length) * voxel_width_;
                        if (end_length < out_length) {
                            atomic_fadd(v.sum_hit_delta,  static_cast<float>(weight * full_delta));
                        } else {
                            atomic_fadd(v.sum_miss_delta, static_cast<float>(weight * full_delta));
                        }
                    }
                    if (current_ray_unbound_) {
                        atomic_fadd(v.num_unbound_rays, 1.0f);
                        atomic_fadd(v.path_length_unbound, static_cast<float>(length_in_voxel * weight));
                    }
                } else {
                    atomic_fadd(v.num_rays_occluded, static_cast<float>(weight));
                    atomic_fadd(v.path_length_occluded, static_cast<float>(length_in_voxel * weight));
                }
            } else {
                // Per-thread sparse map — OOC path and sparse-fallback mode.
                VoxelCoord coord = {p.x(), p.y(), p.z()};
                VoxelGrid::Voxel& v = sparse_voxels_[coord];
                if (type == RayType::OBSERVED) {
                    v.num_beams_observed += 1;
                    v.num_beams_weighted += static_cast<float>(weight);
                    v.path_length_observed += static_cast<float>(length_in_voxel);
                    v.path_length_weighted += static_cast<float>(length_in_voxel * weight);

                    Eigen::Vector3d voxel_center_world = bounds_.min_bound_ + (p.cast<double>() + Eigen::Vector3d(0.5, 0.5, 0.5)) * voxel_width_;
                    double dist_to_center = (voxel_center_world - current_ray_world_start_).norm();

                    v.sum_of_angles += static_cast<float>(zenith_angle * weight);
                    v.sum_sin_azimuth += static_cast<float>(sin_az * weight);
                    v.sum_cos_azimuth += static_cast<float>(cos_az * weight);
                    v.sum_of_laser_distances += static_cast<float>(dist_to_center * weight);

                    if (calc_beam_metrics_) {
                        double beam_radius = tan_half_divergence_ * dist_to_center + 0.5 * beam_diameter_;
                        v.bs_entering += static_cast<float>(kPi * beam_radius * beam_radius * weight);
                        v.sum_bs_path += static_cast<float>(kPi * beam_radius * beam_radius * weight * length_in_voxel);
                    }
                    if (subvoxel_split_ > 0) {
                        Eigen::Vector3d ls = (current_ray_vox_start_ + current_ray_vox_dir_ * in_length  - p.cast<double>()) * subvoxel_split_;
                        Eigen::Vector3d le = (current_ray_vox_start_ + current_ray_vox_dir_ * end_length - p.cast<double>()) * subvoxel_split_;
                        walkSubGrid(ls, le, subvoxel_split_, v.subvoxel_bitmap);
                    }
                    {
                        double full_delta = (out_length - in_length) * voxel_width_;
                        if (end_length < out_length) {
                            v.sum_hit_delta  += static_cast<float>(weight * full_delta);
                        } else {
                            v.sum_miss_delta += static_cast<float>(weight * full_delta);
                        }
                    }
                    if (current_ray_unbound_) {
                        v.num_unbound_rays += static_cast<float>(weight);
                        v.path_length_unbound += static_cast<float>(length_in_voxel * weight);
                    }
                } else {
                    v.num_rays_occluded += static_cast<float>(weight);
                    v.path_length_occluded += static_cast<float>(length_in_voxel * weight);
                }
            }

            return false;
        };

    ray::walkGrid(vox_start, vox_end, walk_lambda);
}

} // namespace ray
