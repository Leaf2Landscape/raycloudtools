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
#include <cmath>     // For std::log1p

namespace ray
{

namespace
{
  // Stage 3 effective free path: eff(z) = −ln(1 − λ₁·z) / λ₁, degrading to z when λ₁ = 0.
  // λ₁·z is clamped below 1 to avoid the log singularity / NaN. std::log1p improves
  // numerical stability for small arguments.
  static inline double effFreePath(double z, double lambda1) {
      if (lambda1 <= 0.0) return z;
      const double x = lambda1 * z;
      const double x_clamped = (x < 1.0 - 1e-9) ? x : (1.0 - 1e-9);
      return -std::log1p(-x_clamped) / lambda1;
  }
} // namespace

VoxelProcessor::VoxelProcessor(const Cuboid& grid_bounds, double voxel_width, const std::string& weighting_method,
                               bool use_occlusion_rays, bool use_flat_top, const std::vector<double>* peaks,
                               bool calc_beam_metrics, double beam_diameter, double tan_half_divergence, int subvoxel_split,
                               const HeightField* dtm, double lambda1)
  : bounds_(grid_bounds),
    voxel_width_(voxel_width),
    voxel_dims_(((grid_bounds.max_bound_ - grid_bounds.min_bound_) / voxel_width).array().ceil().cast<int64_t>()),
    weight_method_(parseWeightMethod(weighting_method)),
    use_occlusion_rays_(use_occlusion_rays),
    use_flat_top_(use_flat_top),
    peaks_(peaks),
    calc_beam_metrics_(calc_beam_metrics),
    beam_diameter_(beam_diameter),
    tan_half_divergence_(tan_half_divergence),
    subvoxel_split_(subvoxel_split),
    row_stride_(voxel_dims_[0]),
    dtm_(dtm),
    lambda1_(lambda1)
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

  // Per-echo weights echo_w[k] and per-segment traversal weights seg_w[k].
  // seg_w[k] is a suffix sum of echo_w: the remaining beam fraction entering segment k.
  float echo_w[kMaxReturnsPerBeam];
  float seg_w[kMaxReturnsPerBeam];
  computeEchoWeights(weight_method_, sorted.data(), N, echo_w);
  computeSegmentWeights(echo_w, N, seg_w);

  // last_nonzero_seg: last k where seg_w[k] > 0 (always true: k <= this index).
  // seg_w is a non-increasing suffix sum, so scanning backward finds the cutoff.
  int last_nonzero_seg = -1;
  for (int k = N - 1; k >= 0; --k) {
    if (seg_w[k] > 0.0f) { last_nonzero_seg = k; break; }
  }

  // --- Mechanism 1: unweighted full-ray walk (weight = 1.0) ---
  // One traversal from beam origin to the last weighted return (sorted[last_nonzero_seg]).
  // This accumulates the unweighted counts (num_beams, path_length,
  // num_miss_rays, subvoxel bitmap). If that return's bound == 0 (unbound/miss ray),
  // traversal still sweeps those voxels — marking them observed/free — but the
  // hit-recording loop below gates on p.bound and will not count the endpoint as a hit.
  const double ray_length = (farthest_pos - beam.beam_origin).norm();
  if (last_nonzero_seg >= 0 && ray_length >= 1e-6) {
    const PointData& walk_end = *sorted[last_nonzero_seg];
    Eigen::Vector3d walk_end_pos(walk_end.x, walk_end.y, walk_end.z);
    Eigen::Vector3d cs = beam.beam_origin, ce = walk_end_pos;
    if (bounds_.clipRay(cs, ce, 1e-10)) {
      Eigen::Vector3d vs = (cs - bounds_.min_bound_) / voxel_width_;
      Eigen::Vector3d ve = (ce - bounds_.min_bound_) / voxel_width_;
      current_ray_vox_start_   = vs;
      current_ray_vox_dir_     = (ve - vs).normalized();
      current_ray_world_start_ = beam.beam_origin;
      current_ray_unbound_ = (walk_end.bound == 0);
      walkGrid(vs, ve, RayType::OBSERVED, 1.0, /*weighted_only=*/false);
    }
  }

  // --- Mechanism 2: segmented weighted walks ---
  // Segment k runs from sorted[k-1] (or beam_origin for k=0) to sorted[k] and is
  // traversed with seg_w[k]. This accumulates the weighted metrics.
  {
    Eigen::Vector3d seg_start = beam.beam_origin;
    for (int k = 0; k <= last_nonzero_seg; ++k) {
      Eigen::Vector3d seg_end(sorted[k]->x, sorted[k]->y, sorted[k]->z);
      const double sw = static_cast<double>(seg_w[k]);
      if (sw > 0.0) {  // always true: k <= last_nonzero_seg
        const double seg_length = (seg_end - seg_start).norm();
        if (seg_length >= 1e-6) {
          Eigen::Vector3d cs = seg_start, ce = seg_end;
          if (bounds_.clipRay(cs, ce, 1e-10)) {
            Eigen::Vector3d vs = (cs - bounds_.min_bound_) / voxel_width_;
            Eigen::Vector3d ve = (ce - bounds_.min_bound_) / voxel_width_;
            current_ray_vox_start_   = vs;
            current_ray_vox_dir_     = (ve - vs).normalized();
            current_ray_world_start_ = seg_start;
            current_ray_unbound_ = (sorted[k]->bound == 0);
            current_seg_foliage_class_ = sorted[k]->foliage_class;
            walkGrid(vs, ve, RayType::OBSERVED, sw, /*weighted_only=*/true);
          }
        }
      }
      seg_start = seg_end;
    }
  }

  // Pre-compute unit ray direction for free-path excess corrections below.
  Eigen::Vector3d ray_dir_unit = Eigen::Vector3d::Zero();
  bool ray_dir_valid = false;
  {
    const Eigen::Vector3d raw_dir = farthest_pos - beam.beam_origin;
    const double raw_len = raw_dir.norm();
    if (raw_len > 1e-12) { ray_dir_unit = raw_dir / raw_len; ray_dir_valid = true; }
  }

  // Record hits for all returns and apply per-stage free-path corrections.
  //
  // Stage 1 (always): subtract post-hit excess from free_path_length so its
  //   net value is Σ(seg_weight × free_path), where free_path = entry→hit for hits.
  //
  // Stage 2 (calc_beam_metrics_ only): accumulate bs_intercepted and apply the same
  //   geometric correction to bs_free_path with the additional π·r² beam-area factor.
  //
  // Stage 3 (future): effective-free-path correction −log(1−λ₁·z)/λ₁ applied here.
  for (int i = 0; i < N; ++i) {
    const PointData& p = *sorted[i];
    if (p.bound == 0) continue;  // unbound ray: traversed as observed, never a hit
    if (seg_w[i] <= 0.0f) continue;  // beam already exhausted before this echo; not a counted hit
    // Weight the intercepted beam section by THIS echo's intercepted fraction (AMAPVox
    // bfIntercepted = echo_w), not the cumulative seg_w. FPL's interceptedBeamSection and bias
    // accumulator are per-echo; using seg_w over-weights early returns of multi-return pulses.
    const double hit_weight = static_cast<double>(echo_w[i]);
    Eigen::Vector3d curr(p.x, p.y, p.z);
    Eigen::Vector3d vox_coord_filled = (curr - bounds_.min_bound_) / voxel_width_;
    int64_t ix = static_cast<int64_t>(std::floor(vox_coord_filled.x()));
    int64_t iy = static_cast<int64_t>(std::floor(vox_coord_filled.y()));
    int64_t iz = static_cast<int64_t>(std::floor(vox_coord_filled.z()));
    if (ix >= 0 && ix < voxel_dims_[0] && iy >= 0 && iy < voxel_dims_[1] && iz >= 0 && iz < voxel_dims_[2]) {
      // Full voxel chord (entry→exit) of the hit voxel. Used only as the PPL hit-chord
      // (sum_hit_delta). The free-path family is NOT corrected here: the weighted segment
      // walk (Mechanism 2) terminates AT the echo, so its seed is already the entry→hit free
      // path (and, for multi-return beams, the post-hit continuation is carried by the next
      // segment at the reduced beam fraction) — matching AMAPVox's weightedFreepathLength.
      double full_transit = 0.0;
      double free_to_hit = 0.0;  // entry→hit free path within this voxel; feeds the FPL bias-correction accumulator
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
        full_transit = std::max(0.0, t_exit - t_entry);
        free_to_hit = std::max(0.0, std::min(full_transit, p.distance_to_sensor - t_entry));
      }

      if (flat_array_) {
        VoxelGrid::Voxel& v = flat_array_[ix + iy * flat_dim_x_ + iz * flat_dim_xy_];
        atomic_iadd(v.num_hits, 1);
        atomic_fadd(v.sum_hit_delta, static_cast<float>(full_transit));
        {
          const int64_t flat_idx = ix + iy * voxel_dims_[0] + iz * voxel_dims_[0] * voxel_dims_[1];
          thread_class_table_[flat_idx][p.classification] += 1.0f;
          if      (p.foliage_class == 2) thread_voxel_lw_[flat_idx].first++;
          else if (p.foliage_class == 3) thread_voxel_lw_[flat_idx].second++;
        }
        if (calc_beam_metrics_) {
          const Eigen::Vector3d vox_center = bounds_.min_bound_ + (Eigen::Vector3d(static_cast<double>(ix), static_cast<double>(iy), static_cast<double>(iz)) + Eigen::Vector3d::Constant(0.5)) * voxel_width_;
          const double dc = (vox_center - beam.beam_origin).norm();
          const double r = tan_half_divergence_ * dc + 0.5 * beam_diameter_;
          atomic_fadd(v.bs_intercepted, static_cast<float>(kPi * r * r * hit_weight));
          atomic_fadd(v.bs_eff_free_path_hits, static_cast<float>(kPi * r * r * hit_weight * effFreePath(free_to_hit, lambda1_)));
        }
      } else {
        VoxelCoord coord = {ix, iy, iz};
        VoxelGrid::Voxel& v = sparse_voxels_[coord];
        v.num_hits += 1;
        v.sum_hit_delta += static_cast<float>(full_transit);
        {
          const int64_t flat_idx = ix + iy * voxel_dims_[0] + iz * voxel_dims_[0] * voxel_dims_[1];
          thread_class_table_[flat_idx][p.classification] += 1.0f;
          if      (p.foliage_class == 2) thread_voxel_lw_[flat_idx].first++;
          else if (p.foliage_class == 3) thread_voxel_lw_[flat_idx].second++;
        }
        if (calc_beam_metrics_) {
          const Eigen::Vector3d vox_center = bounds_.min_bound_ + (Eigen::Vector3d(static_cast<double>(ix), static_cast<double>(iy), static_cast<double>(iz)) + Eigen::Vector3d::Constant(0.5)) * voxel_width_;
          const double dc = (vox_center - beam.beam_origin).norm();
          const double r = tan_half_divergence_ * dc + 0.5 * beam_diameter_;
          v.bs_intercepted += static_cast<float>(kPi * r * r * hit_weight);
          v.bs_eff_free_path_hits += static_cast<float>(kPi * r * r * hit_weight * effFreePath(free_to_hit, lambda1_));
        }
      }
    }
  }

  if (use_occlusion_rays_ && farthest.bound == 1 && ray_length >= 1e-6) {
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

  // Exact-PPL pass (gated): record intercepted beams + per-voxel miss term for the post-solve.
  if (ppl_enabled_) accumulatePpl(beam, sorted, echo_w, N);
}

void VoxelProcessor::accumulatePpl(const BeamData& beam, const std::vector<const PointData*>& sorted,
                                   const float* echo_w, int N)
{
  // Flat-array path only (ppl_enabled_ is set solely on the in-process flat runner).
  if (!flat_array_ || N <= 0) return;

  // Walk origin -> farthest return (clipped to the grid), like Mechanism 1.
  const PointData& farthest = *sorted[N - 1];
  Eigen::Vector3d cs = beam.beam_origin;
  Eigen::Vector3d ce(farthest.x, farthest.y, farthest.z);
  if (!bounds_.clipRay(cs, ce, 1e-10)) return;
  const Eigen::Vector3d vs = (cs - bounds_.min_bound_) / voxel_width_;
  const Eigen::Vector3d ve = (ce - bounds_.min_bound_) / voxel_width_;

  // Precompute each return's voxel coords (only bound returns intercept).
  int rvx[kMaxReturnsPerBeam], rvy[kMaxReturnsPerBeam], rvz[kMaxReturnsPerBeam];
  for (int r = 0; r < N; ++r) {
    const Eigen::Vector3d vc = (Eigen::Vector3d(sorted[r]->x, sorted[r]->y, sorted[r]->z) - bounds_.min_bound_) / voxel_width_;
    rvx[r] = static_cast<int>(std::floor(vc.x()));
    rvy[r] = static_cast<int>(std::floor(vc.y()));
    rvz[r] = static_cast<int>(std::floor(vc.z()));
  }

  double f_in = 1.0;                  // entering beam fraction (full beam)
  const double eps = 1e-6;
  auto ppl_lambda = [&](const Eigen::Vector3i& p, const Eigen::Vector3i& /*target*/,
                        double in_length, double out_length, double /*max_length*/) -> bool {
    if (p.x() < 0 || p.x() >= voxel_dims_[0] || p.y() < 0 || p.y() >= voxel_dims_[1] ||
        p.z() < 0 || p.z() >= voxel_dims_[2])
      return false;
    const double full_chord = (out_length - in_length) * voxel_width_;   // potential path length (full voxel chord)
    if (full_chord <= 0.0) return false;
    // Beam section at the voxel centre (constant section = 1 when beam metrics are off).
    double beam_section = 1.0;
    if (calc_beam_metrics_) {
      const Eigen::Vector3d ctr = bounds_.min_bound_ + (p.cast<double>() + Eigen::Vector3d::Constant(0.5)) * voxel_width_;
      const double dist = (ctr - beam.beam_origin).norm();
      const double rad = tan_half_divergence_ * dist + 0.5 * beam_diameter_;
      beam_section = kPi * rad * rad;
    }
    const int64_t flat = static_cast<int64_t>(p.x()) + static_cast<int64_t>(p.y()) * flat_dim_x_
                       + static_cast<int64_t>(p.z()) * flat_dim_xy_;
    // Intercepted beam fraction in this voxel: one record per bound echo (PL = full chord).
    double intercepted = 0.0;
    for (int r = 0; r < N; ++r) {
      if (sorted[r]->bound == 0) continue;
      if (rvx[r] == p.x() && rvy[r] == p.y() && rvz[r] == p.z()) {
        const double w = static_cast<double>(echo_w[r]);
        intercepted += w;
        ppl_hits_.push_back(PplHit{flat, static_cast<float>(full_chord), static_cast<float>(w * beam_section)});
      }
    }
    // Exiting (non-intercepted) fraction contributes the miss term over the full chord.
    const double exiting = f_in - intercepted;
    if (exiting > 0.0)
      atomic_fadd(flat_array_[flat].ppl_miss_wL, static_cast<float>(exiting * beam_section * full_chord));
    f_in -= intercepted;
    return f_in < eps;   // stop once the beam is exhausted
  };
  ray::walkGrid(vs, ve, ppl_lambda);
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

void VoxelProcessor::walkGrid(const Eigen::Vector3d &vox_start, const Eigen::Vector3d &vox_end, RayType type, double weight,
                              bool weighted_only)
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
            double full_chord = (out_length - in_length) * voxel_width_;

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
                    Eigen::Vector3d voxel_center_world = bounds_.min_bound_ + (p.cast<double>() + Eigen::Vector3d(0.5, 0.5, 0.5)) * voxel_width_;
                    double dist_to_center = (voxel_center_world - current_ray_world_start_).norm();

                    if (!weighted_only) {
                        // Mechanism 1: unweighted counts from the single full-ray walk.
                        atomic_iadd(v.num_beams, 1);
                        atomic_fadd(v.path_length, static_cast<float>(full_chord));
                        atomic_fadd(v.path_length_sq_raw, static_cast<float>(full_chord * full_chord));
                        if (subvoxel_split_ > 0) {
                            Eigen::Vector3d ls = (current_ray_vox_start_ + current_ray_vox_dir_ * in_length  - p.cast<double>()) * subvoxel_split_;
                            Eigen::Vector3d le = (current_ray_vox_start_ + current_ray_vox_dir_ * end_length - p.cast<double>()) * subvoxel_split_;
                            uint64_t bits = 0;
                            walkSubGrid(ls, le, subvoxel_split_, bits);
                            for (uint64_t b = bits; b; b &= b - 1)
                              atomic_inc_u8_sat(v.subvoxel_counts[__builtin_ctzll(b)]);
                        }
                        if (end_length >= out_length - 1e-6) {
                            if (!current_ray_unbound_) atomic_iadd(v.num_miss_rays, 1);
                            atomic_fadd(v.sum_miss_delta, static_cast<float>(full_chord));
                        }
                    } else {
                        // Mechanism 2: weighted metrics from the segmented walks.
                        atomic_fadd(v.num_beams_weighted, static_cast<float>(weight));
                        atomic_fadd(v.free_path_length, static_cast<float>(length_in_voxel * weight));
                        {
                            const float fpl_c = static_cast<float>(length_in_voxel * weight);
                            if      (current_seg_foliage_class_ == 1) atomic_fadd(v.free_path_length_plant, fpl_c);
                            else if (current_seg_foliage_class_ == 2) atomic_fadd(v.free_path_length_leaf,  fpl_c);
                            else if (current_seg_foliage_class_ == 3) atomic_fadd(v.free_path_length_wood,  fpl_c);
                        }
                        atomic_fadd(v.effective_free_path_length, static_cast<float>(effFreePath(length_in_voxel, lambda1_) * weight));
                        atomic_fadd(v.sum_of_angles, static_cast<float>(zenith_angle * weight));
                        atomic_fadd(v.sum_sin_azimuth, static_cast<float>(sin_az * weight));
                        atomic_fadd(v.sum_cos_azimuth, static_cast<float>(cos_az * weight));
                        atomic_fadd(v.sum_of_laser_distances, static_cast<float>(dist_to_center * weight));

                        if (calc_beam_metrics_) {
                            double beam_radius = tan_half_divergence_ * dist_to_center + 0.5 * beam_diameter_;
                            atomic_fadd(v.bs_entering, static_cast<float>(kPi * beam_radius * beam_radius * weight));
                            atomic_fadd(v.bs_free_path, static_cast<float>(kPi * beam_radius * beam_radius * weight * length_in_voxel));
                            atomic_fadd(v.bs_effective_free_path, static_cast<float>(kPi * beam_radius * beam_radius * weight * effFreePath(length_in_voxel, lambda1_)));
                            if (end_length >= out_length - 1e-6) {
                                atomic_fadd(v.bs_potential, static_cast<float>(kPi * beam_radius * beam_radius * weight));
                            }
                        }
                        if (current_ray_unbound_) {
                            atomic_iadd(v.num_unbound_rays, 1);
                            atomic_fadd(v.path_length_unbound, static_cast<float>(length_in_voxel * weight));
                        }
                    }
                } else {
                    atomic_iadd(v.num_rays_occluded, 1);
                    atomic_fadd(v.path_length_occluded, static_cast<float>(length_in_voxel * weight));
                }
            } else {
                // Per-thread sparse map — OOC path and sparse-fallback mode.
                VoxelCoord coord = {p.x(), p.y(), p.z()};
                VoxelGrid::Voxel& v = sparse_voxels_[coord];
                if (type == RayType::OBSERVED) {
                    Eigen::Vector3d voxel_center_world = bounds_.min_bound_ + (p.cast<double>() + Eigen::Vector3d(0.5, 0.5, 0.5)) * voxel_width_;
                    double dist_to_center = (voxel_center_world - current_ray_world_start_).norm();

                    if (!weighted_only) {
                        // Mechanism 1: unweighted counts from the single full-ray walk.
                        v.num_beams += 1;
                        v.path_length += static_cast<float>(full_chord);
                        v.path_length_sq_raw += static_cast<float>(full_chord * full_chord);
                        if (subvoxel_split_ > 0) {
                            Eigen::Vector3d ls = (current_ray_vox_start_ + current_ray_vox_dir_ * in_length  - p.cast<double>()) * subvoxel_split_;
                            Eigen::Vector3d le = (current_ray_vox_start_ + current_ray_vox_dir_ * end_length - p.cast<double>()) * subvoxel_split_;
                            uint64_t bits = 0;
                            walkSubGrid(ls, le, subvoxel_split_, bits);
                            for (uint64_t b = bits; b; b &= b - 1) {
                              int i = __builtin_ctzll(b);
                              v.subvoxel_counts[i] = static_cast<uint8_t>(std::min(255, static_cast<int>(v.subvoxel_counts[i]) + 1));
                            }
                        }
                        if (end_length >= out_length - 1e-6) {
                            if (!current_ray_unbound_) v.num_miss_rays += 1;
                            v.sum_miss_delta += static_cast<float>(full_chord);
                        }
                    } else {
                        // Mechanism 2: weighted metrics from the segmented walks.
                        v.num_beams_weighted += static_cast<float>(weight);
                        v.free_path_length += static_cast<float>(length_in_voxel * weight);
                        {
                            const float fpl_c = static_cast<float>(length_in_voxel * weight);
                            if      (current_seg_foliage_class_ == 1) v.free_path_length_plant += fpl_c;
                            else if (current_seg_foliage_class_ == 2) v.free_path_length_leaf  += fpl_c;
                            else if (current_seg_foliage_class_ == 3) v.free_path_length_wood  += fpl_c;
                        }
                        v.effective_free_path_length += static_cast<float>(effFreePath(length_in_voxel, lambda1_) * weight);
                        v.sum_of_angles += static_cast<float>(zenith_angle * weight);
                        v.sum_sin_azimuth += static_cast<float>(sin_az * weight);
                        v.sum_cos_azimuth += static_cast<float>(cos_az * weight);
                        v.sum_of_laser_distances += static_cast<float>(dist_to_center * weight);

                        if (calc_beam_metrics_) {
                            double beam_radius = tan_half_divergence_ * dist_to_center + 0.5 * beam_diameter_;
                            v.bs_entering += static_cast<float>(kPi * beam_radius * beam_radius * weight);
                            v.bs_free_path += static_cast<float>(kPi * beam_radius * beam_radius * weight * length_in_voxel);
                            v.bs_effective_free_path += static_cast<float>(kPi * beam_radius * beam_radius * weight * effFreePath(length_in_voxel, lambda1_));
                            if (end_length >= out_length - 1e-6) {
                                v.bs_potential += static_cast<float>(kPi * beam_radius * beam_radius * weight);
                            }
                        }
                        if (current_ray_unbound_) {
                            v.num_unbound_rays += 1;
                            v.path_length_unbound += static_cast<float>(length_in_voxel * weight);
                        }
                    }
                } else {
                    v.num_rays_occluded += 1;
                    v.path_length_occluded += static_cast<float>(length_in_voxel * weight);
                }
            }

            return false;
        };

    ray::walkGrid(vox_start, vox_end, walk_lambda);
}

} // namespace ray
