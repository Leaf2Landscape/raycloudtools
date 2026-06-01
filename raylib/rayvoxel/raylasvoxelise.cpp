// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file implements the VoxelGrid data container and the processing
// strategies for voxelization. The generateVoxelGrid function acts as the
// main orchestrator, delegating the core work to the chosen strategy:
// in-memory (parallel or single-threaded) or out-of-core.

#include "raylib/rayprogress.h"
#include "raylib/rayprogressthread.h"
#include "raylib/rayparse.h"
#include "raylib/rayutils.h"
#include "raylib/rayunused.h"
#include "raylib/rayply.h"
#include "raylib/raylaz.h"
#include "raylib/rayvoxel/rayvox.h"
#include "raylib/rayvoxel/raylasvoxelise.h"
#include "raylib/rayvoxel/raylasvoxelwriter.h"
#include "raylib/rayvoxel/raylasvoxelrefine.h"
#include "raylib/rayvoxel/raylasvegmetrics.h"
#include "raylib/rayvoxel/raylasvoxelprocessor.h"
#include "raylib/rayvoxel/raylasthreadsafequeue.h"
#include "raylib/rayvoxel/raylasheightfield.h"
#include "raylib/rayvoxel/raylasbinaryio.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <fstream>
#include <limits>
#include <cmath>
#include <mutex>
#include <thread>
#include <memory>
#include <queue>
#include <filesystem>
#include <atomic>

namespace ray
{

// The passthrough buffer produced by readLas has 10 standard LAS field bytes per point
// before any original sensor extra bytes. Field byte [0] packs return number (low nibble)
// and number of returns (high nibble); byte [2] is the extended classification.
static const uint16_t kPassthroughStdBytes = 10;

// ==================================================================================
// VoxelGrid Class Implementation
// ==================================================================================

VoxelGrid::VoxelGrid(const Cuboid &grid_bounds, double vox_width, size_t reservation_size)
    : bounds_(grid_bounds), voxel_width_(vox_width)
{
  Eigen::Vector3d extent = bounds_.max_bound_ - bounds_.min_bound_;
  Eigen::Vector3d temp = (extent / voxel_width_).array().ceil();
  voxel_dims_ = temp.cast<int64_t>();

  const int64_t max_reasonable_dim = 1000000000L;
  if (voxel_dims_.maxCoeff() > max_reasonable_dim) {
    throw std::runtime_error("VoxelGrid Error: Resulting dimensions on one or more axes are unreasonably large.");
  }

  peaks_.resize(voxel_dims_[0] * voxel_dims_[1], std::numeric_limits<double>::lowest());
  std::cout << "Initialized VoxelGrid with conceptual dimensions: " << voxel_dims_.transpose()
            << " (using sparse storage)" << std::endl;

  // Pre-allocate memory in the hash map to avoid rehashing during processing.
  if (reservation_size > 0) {
    std::cout << "Reserving space for approximately " << reservation_size << " voxels..." << std::endl;
    sparse_voxels_.reserve(reservation_size);
  }
}

void VoxelGrid::merge(const VoxelProcessor& processor)
{
  // This lock ensures that multiple threads can safely merge their results
  // into the main grid's map without causing data corruption.
  std::lock_guard<std::mutex> lock(merge_mutex_);
  for (const auto& pair : processor.getMap()) {
    sparse_voxels_[pair.first] += pair.second;
  }
}

void VoxelGrid::take(VoxelProcessor& processor)
{
  // This is a non-locking, single-threaded optimization.
  // It moves the map from the processor directly into the grid's map.
  sparse_voxels_ = processor.takeMap();
}

void VoxelGrid::absorbMap(VoxelProcessor::Map&& m)
{
  for (auto& pair : m) {
    sparse_voxels_[pair.first] += pair.second;
  }
}

VoxelGrid::VoxelState VoxelGrid::getVoxelState(int64_t i, int64_t j, int64_t k) const {
    VoxelCoord coord = {i, j, k};
    auto it = sparse_voxels_.find(coord);

    if (it == sparse_voxels_.end()) {
        // If the voxel coordinate is not in the map, it was never touched by any ray.
        return VoxelState::UNOBSERVED;
    }

    // If it is in the map, determine its state from the stored voxel data.
    const Voxel& v = it->second;
    if (v.is_filled) return VoxelState::FILLED;
    if (v.num_rays_observed > 0) return VoxelState::EMPTY;
    if (v.num_rays_occluded > 0) return VoxelState::OCCLUDED;

    // This case should not be hit if a voxel exists in the map, but it's a safe fallback.
    return VoxelState::UNOBSERVED;
}

const VoxelGrid::Voxel& VoxelGrid::getVoxel(int64_t i, int64_t j, int64_t k) const {
    // To safely handle requests for voxels that don't exist in the map (since this is a
    // const method and cannot modify the map), we return a reference to a static empty voxel.
    static const Voxel empty_voxel;
    VoxelCoord coord = {i, j, k};
    auto it = sparse_voxels_.find(coord);

    if (it != sparse_voxels_.end()) {
        return it->second; // Return a const reference to the actual voxel.
    } else {
        return empty_voxel; // Return the static empty one if not found.
    }
}

int64_t VoxelGrid::getIndex(int64_t i, int64_t j, int64_t k) const {
    return i + j * voxel_dims_[0] + k * voxel_dims_[0] * voxel_dims_[1];
}

// ==================================================================================
// Voxel Method Implementations
// ==================================================================================

double VoxelGrid::Voxel::pad_bv_total() const
{
  const double eps = 1e-10; // Avoid division by zero
  if (num_rays_observed < 2.0f) return 0.0;
  const double spherical_distribution_scale = 2.0;
  return spherical_distribution_scale * (num_rays_observed - 1.0f) * num_hits / (eps + num_rays_observed * path_length_observed);
}

double VoxelGrid::Voxel::transmittance() const
{
  if (bs_entering < 1e-10) return 1.0;
  double transmitted = bs_entering - bs_intercepted;
  return std::max(0.0, transmitted / bs_entering);
}

// ==================================================================================
// Local Helpers
// ==================================================================================

namespace {

// Build a PointData record from a single ray (start -> end) and its passthrough bytes.
// The stride is the per-point passthrough byte count (10 standard fields + sensor extras).
// distance_to_sensor is computed as the ray length (end - start).norm().
PointData makePointData(const Eigen::Vector3d& start, const Eigen::Vector3d& end,
                        double gps_time, int32_t beam_id,
                        const std::vector<uint8_t>& passthrough, size_t index, uint16_t stride)
{
  PointData pd;
  pd.beam_origin       = start;                       // sx,sy,sz already decoded by readLas
  pd.x                 = end.x();
  pd.y                 = end.y();
  pd.z                 = end.z();
  pd.gps_time          = gps_time;
  pd.beam_id           = beam_id;
  const size_t base    = index * stride;
  if (passthrough.size() >= base + stride) {
    pd.classification    = passthrough[base + 2];
    pd.return_number     = passthrough[base + 0] & 0x0F;
    pd.number_of_returns = (passthrough[base + 0] >> 4) & 0x0F;
  } else {
    pd.classification    = 0;
    pd.return_number     = 1;
    pd.number_of_returns = 1;
  }
  pd.distance_to_sensor = (end - start).norm();
  return pd;
}

} // anonymous namespace

// ==================================================================================
// Processing Strategy Implementations
// ==================================================================================

// Abstract base class for all processing strategies.
class ProcessingStrategy {
public:
  virtual ~ProcessingStrategy() = default;
  virtual bool execute(const std::string& cloud_name, VoxelGrid& grid,
                       const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                       bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                       const HeightField* dtm) = 0;
};

// --- In-Memory Strategy (Options 1 & 2) ---
class InProcessStrategy : public ProcessingStrategy {
public:
  InProcessStrategy(size_t num_threads) : num_threads_(num_threads) {}
  bool execute(const std::string& cloud_name, VoxelGrid& grid,
               const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
               bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
               const HeightField* dtm) override;
private:
  size_t num_threads_;
};

// --- Out-of-Core Strategy (Option 3) ---
class OutOfCoreStrategy : public ProcessingStrategy {
public:
  OutOfCoreStrategy(size_t num_threads, size_t ram_budget_mb) : num_threads_(num_threads), ram_budget_mb_(ram_budget_mb) {}
  bool execute(const std::string& cloud_name, VoxelGrid& grid,
               const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
               bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
               const HeightField* dtm) override;
private:
  bool createShards(const std::string& cloud_name, VoxelGrid& grid,
                    const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                    bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                    const HeightField* dtm, std::vector<std::string>& out_shard_paths);

  bool mergeShards(const std::vector<std::string>& shard_paths, VoxelGrid& grid);

  size_t num_threads_;
  size_t ram_budget_mb_;
};


bool InProcessStrategy::execute(const std::string& cloud_name, VoxelGrid& grid,
                                const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                                bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                                const HeightField* dtm)
{
  // Determine the final number of threads to use
  size_t resolved_threads = num_threads_;
  if (resolved_threads == 0) { // Auto-detect
    resolved_threads = std::thread::hardware_concurrency();
    if (resolved_threads == 0) resolved_threads = 1; // Fallback
  }

  // Common parameters for all processors
  double tan_half_divergence = calc_beam_metrics ? tan(0.5 * beam_divergence) : 0.0;
  const std::vector<double>* peaks_ptr = apply_flat_top ? &grid.getPeaks() : nullptr;

  // Determine the per-point passthrough stride from the ray cloud's extra-byte header.
  uint16_t orig_extra_size = 0;
  std::vector<uint8_t> extra_bytes_vlr;
  readLasExtraBytesVlr(cloud_name, orig_extra_size, extra_bytes_vlr);
  const uint16_t stride = static_cast<uint16_t>(kPassthroughStdBytes + orig_extra_size);

  if (resolved_threads > 1) {
    // --- OPTION 2: Parallel Producer-Consumer Implementation ---
    std::cout << "Processing point cloud using " << resolved_threads << " parallel threads (in-memory)..." << std::endl;
    ThreadSafeQueue<BeamData> beam_queue(resolved_threads * 100);
    std::vector<std::thread> threads;

    // The worker task lambda
    std::vector<VoxelProcessor::Map> worker_maps(resolved_threads);

    auto worker_task = [&](size_t thread_idx) {
      VoxelProcessor processor(grid.getBounds(), grid.getVoxelWidth(), weighting_method, use_occlusion,
                               apply_flat_top, peaks_ptr, calc_beam_metrics, beam_diameter,
                               tan_half_divergence, subvoxel_split, dtm);
      BeamData beam;
      while(beam_queue.pop(beam)) {
        processor.processBeam(beam);
      }
      worker_maps[thread_idx] = processor.takeMap();
    };

    // Launch worker threads
    for (size_t i = 0; i < resolved_threads; ++i) {
      threads.emplace_back(worker_task, i);
    }

    size_t num_bounded = 0;
    std::vector<uint8_t> passthrough;
    std::vector<int32_t> beam_ids_chunk;
    uint16_t pt_extra = 0;
    static bool not_raycloud_warned = false;
    // Beam accumulator state, persisting across readLas chunk calls.
    double pending_gps_time = std::numeric_limits<double>::quiet_NaN();
    int32_t pending_beam_id = -1;
    std::vector<PointData> pending_returns;
    Eigen::Vector3d pending_beam_origin;
    auto flush_beam = [&]() {
      if (!pending_returns.empty()) {
        BeamData beam;
        beam.beam_origin = pending_beam_origin;
        beam.gps_time    = pending_gps_time;
        beam.returns     = std::move(pending_returns);
        beam_queue.push(std::move(beam));
        pending_returns.clear();
      }
    };
    ray::readLas(cloud_name,
      [&](std::vector<Eigen::Vector3d>& starts, std::vector<Eigen::Vector3d>& ends,
          std::vector<double>& times, std::vector<ray::RGBA>& /*colours*/) {
        for (size_t i = 0; i < ends.size(); ++i) {
          // readLas only supplies ray starts for ray-cloud files (files with sx,sy,sz extra bytes).
          if (starts.empty() || starts[i] == ends[i]) {
            if (!not_raycloud_warned) {
              std::cerr << "Warning: input is not a ray cloud (no sx,sy,sz ray starts); skipping points." << std::endl;
              not_raycloud_warned = true;
            }
            continue;
          }
          const int32_t bid = (i < beam_ids_chunk.size()) ? beam_ids_chunk[i] : -1;
          PointData pd = makePointData(starts[i], ends[i], times[i], bid, passthrough, i, stride);
          const bool new_beam = beam_ids_chunk.empty()
            ? (pd.gps_time != pending_gps_time)
            : (pd.beam_id != pending_beam_id);
          if (new_beam) {
            flush_beam();
            pending_gps_time    = pd.gps_time;
            pending_beam_id     = pd.beam_id;
            pending_beam_origin = starts[i];
          }
          pending_returns.push_back(pd);
        }
        passthrough.clear();
        beam_ids_chunk.clear();
      }, num_bounded, 255.0, nullptr, 1000000, nullptr, &passthrough, &pt_extra, nullptr, nullptr, &beam_ids_chunk);

    flush_beam(); // Flush the final beam.
    beam_queue.notify_done(); // Signal that no more beams are coming
    // Join all threads
    for (auto& t : threads) { t.join(); }
    // Serial reduction — no locks needed after join
    for (size_t i = 0; i < resolved_threads; ++i) {
      grid.absorbMap(std::move(worker_maps[i]));
    }
    std::cout << "Parallel processing finished." << std::endl;

  } else {
    // --- OPTION 1: Single-Threaded Implementation ---
    std::cout << "Processing point cloud using 1 thread (in-memory)..." << std::endl;
    VoxelProcessor processor(grid.getBounds(), grid.getVoxelWidth(), weighting_method, use_occlusion,
                             apply_flat_top, peaks_ptr, calc_beam_metrics, beam_diameter,
                             tan_half_divergence, subvoxel_split, dtm);

    size_t num_bounded = 0;
    std::vector<uint8_t> passthrough;
    std::vector<int32_t> beam_ids_chunk;
    uint16_t pt_extra = 0;
    static bool not_raycloud_warned = false;
    // Beam accumulator state, persisting across readLas chunk calls.
    double pending_gps_time = std::numeric_limits<double>::quiet_NaN();
    int32_t pending_beam_id = -1;
    std::vector<PointData> pending_returns;
    Eigen::Vector3d pending_beam_origin;
    auto flush_beam = [&]() {
      if (!pending_returns.empty()) {
        BeamData beam;
        beam.beam_origin = pending_beam_origin;
        beam.gps_time    = pending_gps_time;
        beam.returns     = std::move(pending_returns);
        processor.processBeam(beam);
        pending_returns.clear();
      }
    };
    ray::readLas(cloud_name,
      [&](std::vector<Eigen::Vector3d>& starts, std::vector<Eigen::Vector3d>& ends,
          std::vector<double>& times, std::vector<ray::RGBA>& /*colours*/) {
        for (size_t i = 0; i < ends.size(); ++i) {
          if (starts.empty() || starts[i] == ends[i]) {
            if (!not_raycloud_warned) {
              std::cerr << "Warning: input is not a ray cloud (no sx,sy,sz ray starts); skipping points." << std::endl;
              not_raycloud_warned = true;
            }
            continue;
          }
          const int32_t bid = (i < beam_ids_chunk.size()) ? beam_ids_chunk[i] : -1;
          PointData pd = makePointData(starts[i], ends[i], times[i], bid, passthrough, i, stride);
          const bool new_beam = beam_ids_chunk.empty()
            ? (pd.gps_time != pending_gps_time)
            : (pd.beam_id != pending_beam_id);
          if (new_beam) {
            flush_beam();
            pending_gps_time    = pd.gps_time;
            pending_beam_id     = pd.beam_id;
            pending_beam_origin = starts[i];
          }
          pending_returns.push_back(pd);
        }
        passthrough.clear();
        beam_ids_chunk.clear();
      }, num_bounded, 255.0, nullptr, 1000000, nullptr, &passthrough, &pt_extra, nullptr, nullptr, &beam_ids_chunk);

    flush_beam(); // Flush the final beam.
    // Move results from the single processor to the main grid
    grid.take(processor);
  }

  return true;
}

// ==================================================================================
// Out-of-Core Strategy Implementation
// ==================================================================================

// Helper class for the k-way merge in Phase 2 of the out-of-core strategy.
class ShardMerger {
private:
  // An entry in the priority queue for merging.
  struct MergeEntry {
    VoxelCoord coord;
    VoxelGrid::Voxel voxel;
    size_t shard_index;

    // Custom comparator to make the priority queue a min-heap.
    bool operator>(const MergeEntry& other) const {
      if (coord.z != other.coord.z) return coord.z > other.coord.z;
      if (coord.y != other.coord.y) return coord.y > other.coord.y;
      return coord.x > other.coord.x;
    }
  };

  const std::vector<std::string>& shard_paths_;
  VoxelGrid& target_grid_;

public:
  ShardMerger(const std::vector<std::string>& paths, VoxelGrid& grid)
    : shard_paths_(paths), target_grid_(grid) {}

  bool merge() {
    std::cout << "Phase 2: Merging " << shard_paths_.size() << " temporary shards..." << std::endl;
    std::vector<std::ifstream> shard_streams;
    shard_streams.reserve(shard_paths_.size());
    for(const auto& path : shard_paths_) {
      shard_streams.emplace_back(path, std::ios::binary);
      if (!shard_streams.back().is_open()) {
        std::cerr << "Error: Could not open shard file for reading: " << path << std::endl;
        return false;
      }
    }

    // Min-priority queue to efficiently find the next voxel to merge.
    std::priority_queue<MergeEntry, std::vector<MergeEntry>, std::greater<MergeEntry>> pq;

    // Prime the queue with the first entry from each shard.
    for (size_t i = 0; i < shard_streams.size(); ++i) {
      VoxelCoord coord;
      VoxelGrid::Voxel voxel;
      if (readVoxelData(shard_streams[i], coord, voxel)) {
        pq.push({coord, voxel, i});
      }
    }

    if (pq.empty()) {
        std::cout << "No data found in shards to merge." << std::endl;
        return true;
    }

    VoxelCoord current_coord = pq.top().coord;
    VoxelGrid::Voxel accumulator = {};

    while (!pq.empty()) {
      MergeEntry entry = pq.top();
      pq.pop();

      if (entry.coord == current_coord) {
        // Accumulate data for the same voxel coordinate.
        accumulator += entry.voxel;
      } else {
        // New coordinate found; write the accumulated data for the previous one.
        target_grid_.getSparseVoxels_internal()[current_coord] = accumulator;
        current_coord = entry.coord;
        accumulator = entry.voxel;
      }

      // Read the next entry from the same shard to replace the one we just processed.
      VoxelCoord next_coord;
      VoxelGrid::Voxel next_voxel;
      if (readVoxelData(shard_streams[entry.shard_index], next_coord, next_voxel)) {
        pq.push({next_coord, next_voxel, entry.shard_index});
      }
    }

    // Write the last accumulated voxel.
    target_grid_.getSparseVoxels_internal()[current_coord] = accumulator;

    std::cout << "Merge complete. Final grid has " << target_grid_.getSparseVoxels().size() << " voxels." << std::endl;
    return true;
  }
};


bool OutOfCoreStrategy::execute(const std::string& cloud_name, VoxelGrid& grid,
                                const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                                bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                                const HeightField* dtm) {
    std::vector<std::string> shard_paths;
    std::cout << "Starting out-of-core processing..." << std::endl;

    if (!createShards(cloud_name, grid, weighting_method, use_occlusion, apply_flat_top,
                      calc_beam_metrics, beam_diameter, beam_divergence, subvoxel_split, dtm, shard_paths)) {
        std::cerr << "Error: Failed during sharding phase." << std::endl;
        return false;
    }

    if (!mergeShards(shard_paths, grid)) {
        std::cerr << "Error: Failed during merge phase." << std::endl;
        return false;
    }

    std::cout << "Cleaning up temporary shard files..." << std::endl;
    std::filesystem::path temp_dir;
    if (!shard_paths.empty()) {
        temp_dir = std::filesystem::path(shard_paths[0]).parent_path();
    }
    for (const auto& path : shard_paths) {
        std::error_code ec;
        if (std::filesystem::exists(path)) {
            std::filesystem::remove(path, ec);
            if (ec) {
                std::cerr << "Warning: Could not remove temporary file " << path << ": " << ec.message() << std::endl;
            }
        }
    }
    // Clean up temp directory if it exists and is empty
    if (!temp_dir.empty() && std::filesystem::exists(temp_dir) && std::filesystem::is_empty(temp_dir)) {
        std::filesystem::remove(temp_dir);
    }
    return true;
}

bool OutOfCoreStrategy::createShards(const std::string& cloud_name, VoxelGrid& grid,
                                     const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                                     bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                                     const HeightField* dtm, std::vector<std::string>& out_shard_paths) {
    std::cout << "Phase 1: Processing points and writing to temporary shards..." << std::endl;

    size_t resolved_threads = num_threads_ == 0 ? std::thread::hardware_concurrency() : num_threads_;
    if (resolved_threads == 0) resolved_threads = 1;

    size_t ram_per_thread_bytes = ram_budget_mb_ * 1024 * 1024;
    // Estimate size of a map entry: Key + Value + overhead (approx 2 pointers) + map node.
    size_t voxel_pair_size_approx = sizeof(VoxelCoord) + sizeof(VoxelGrid::Voxel) + sizeof(void*)*2 + sizeof(std::pair<U8,float>)*2;
    size_t max_voxels_per_thread = (ram_per_thread_bytes / voxel_pair_size_approx);

    std::cout << "Out-of-core settings: " << resolved_threads << " threads, " << ram_budget_mb_ << "MB RAM budget per thread." << std::endl;
    std::cout << "Each thread will flush to disk approx. every " << max_voxels_per_thread << " unique voxels." << std::endl;

    std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "raylasvoxel_shards";
    std::filesystem::create_directories(temp_dir);

    ThreadSafeQueue<BeamData> beam_queue(resolved_threads * 100);
    std::vector<std::thread> threads;
    std::vector<std::string> shard_paths_list;
    std::mutex shard_paths_mutex;
    std::atomic<int> chunk_counter = 0;

    double tan_half_divergence = calc_beam_metrics ? tan(0.5 * beam_divergence) : 0.0;
    const std::vector<double>* peaks_ptr = apply_flat_top ? &grid.getPeaks() : nullptr;

    auto worker_task = [&](int thread_id) {
      VoxelProcessor processor(grid.getBounds(), grid.getVoxelWidth(), weighting_method, use_occlusion,
                               apply_flat_top, peaks_ptr, calc_beam_metrics, beam_diameter,
                               tan_half_divergence, subvoxel_split, dtm);
      BeamData beam;

      while(beam_queue.pop(beam)) {
        processor.processBeam(beam);
        if (processor.size() >= max_voxels_per_thread) {
          int chunk_id = chunk_counter++;
          std::filesystem::path shard_path = temp_dir /
              ("voxel_shard_t" + std::to_string(thread_id) + "_c" + std::to_string(chunk_id) + ".bin");

          if (processor.flushToShard(shard_path.string())) {
            std::lock_guard<std::mutex> lock(shard_paths_mutex);
            shard_paths_list.push_back(shard_path.string());
          }
          processor.clear();
        }
      }

      if (processor.size() > 0) {
        int chunk_id = chunk_counter++;
        std::filesystem::path shard_path = temp_dir /
            ("voxel_shard_t" + std::to_string(thread_id) + "_c" + std::to_string(chunk_id) + ".bin");

        if (processor.flushToShard(shard_path.string())) {
          std::lock_guard<std::mutex> lock(shard_paths_mutex);
          shard_paths_list.push_back(shard_path.string());
        }
      }
    };

    for (size_t i = 0; i < resolved_threads; ++i) {
      threads.emplace_back(worker_task, i);
    }

    // --- Producer loop ---
    // Determine the per-point passthrough stride from the ray cloud's extra-byte header.
    uint16_t orig_extra_size = 0;
    std::vector<uint8_t> extra_bytes_vlr;
    readLasExtraBytesVlr(cloud_name, orig_extra_size, extra_bytes_vlr);
    const uint16_t stride = static_cast<uint16_t>(kPassthroughStdBytes + orig_extra_size);

    size_t num_bounded = 0;
    std::vector<uint8_t> passthrough;
    std::vector<int32_t> beam_ids_chunk;
    uint16_t pt_extra = 0;
    static bool not_raycloud_warned = false;
    // Beam accumulator state, persisting across readLas chunk calls.
    double pending_gps_time = std::numeric_limits<double>::quiet_NaN();
    int32_t pending_beam_id = -1;
    std::vector<PointData> pending_returns;
    Eigen::Vector3d pending_beam_origin;
    auto flush_beam = [&]() {
      if (!pending_returns.empty()) {
        BeamData beam;
        beam.beam_origin = pending_beam_origin;
        beam.gps_time    = pending_gps_time;
        beam.returns     = std::move(pending_returns);
        beam_queue.push(std::move(beam));
        pending_returns.clear();
      }
    };
    ray::readLas(cloud_name,
      [&](std::vector<Eigen::Vector3d>& starts, std::vector<Eigen::Vector3d>& ends,
          std::vector<double>& times, std::vector<ray::RGBA>& /*colours*/) {
        for (size_t i = 0; i < ends.size(); ++i) {
          if (starts.empty() || starts[i] == ends[i]) {
            if (!not_raycloud_warned) {
              std::cerr << "Warning: input is not a ray cloud (no sx,sy,sz ray starts); skipping points." << std::endl;
              not_raycloud_warned = true;
            }
            continue;
          }
          const int32_t bid = (i < beam_ids_chunk.size()) ? beam_ids_chunk[i] : -1;
          PointData pd = makePointData(starts[i], ends[i], times[i], bid, passthrough, i, stride);
          const bool new_beam = beam_ids_chunk.empty()
            ? (pd.gps_time != pending_gps_time)
            : (pd.beam_id != pending_beam_id);
          if (new_beam) {
            flush_beam();
            pending_gps_time    = pd.gps_time;
            pending_beam_id     = pd.beam_id;
            pending_beam_origin = starts[i];
          }
          pending_returns.push_back(pd);
        }
        passthrough.clear();
        beam_ids_chunk.clear();
      }, num_bounded, 255.0, nullptr, 1000000, nullptr, &passthrough, &pt_extra, nullptr, nullptr, &beam_ids_chunk);

    flush_beam(); // Flush the final beam.
    beam_queue.notify_done();
    for (auto& t : threads) { t.join(); }

    out_shard_paths = std::move(shard_paths_list);
    return true;
}

bool OutOfCoreStrategy::mergeShards(const std::vector<std::string>& shard_paths, VoxelGrid& grid) {
    ShardMerger merger(shard_paths, grid);
    return merger.merge();
}

// ==================================================================================
// Main Orchestrator Function
// ==================================================================================

bool generateVoxelGrid(const VoxelizationParameters& params)
{
    // This function is now the main entry point and orchestrator.
    // It sets up the shared components and delegates to the chosen strategy.

    auto isZeroVector = [](const Eigen::Vector3d &vec) -> bool { return vec.isApprox(Eigen::Vector3d::Zero()); };
    Cuboid user_bounds;

    double beam_diameter = 0.0;
    double beam_divergence = 0.0;

    if (params.calc_beam_metrics) {
        if (!params.laser_spec_name.empty()) {
            LaserSpecManager spec_manager;
            LaserSpecification spec;
            if (spec_manager.getSpec(params.laser_spec_name, spec)) {
                beam_diameter = spec.beam_diameter_at_exit;
                beam_divergence = spec.beam_divergence;
                std::cout << "Using predefined laser spec: '" << spec.name << "'" << std::endl;
            } else {
                std::cerr << "Error: Predefined laser spec '" << params.laser_spec_name << "' not found." << std::endl;
                return false;
            }
        } else {
            beam_diameter = params.beam_params.x();
            beam_divergence = params.beam_params.y();
            std::cout << "Using manual beam parameters: diameter=" << beam_diameter << "m, divergence=" << beam_divergence << "rad" << std::endl;
        }
    }

    if (isZeroVector(params.grid_bounds_min) && isZeroVector(params.grid_bounds_max)) {
        std::cout << "Auto-detecting grid bounds from file..." << std::endl;
        Eigen::Vector3d bounds_min = Eigen::Vector3d::Constant(std::numeric_limits<double>::max());
        Eigen::Vector3d bounds_max = Eigen::Vector3d::Constant(std::numeric_limits<double>::lowest());
        size_t num_bounded = 0;
        if (!ray::readLas(params.cloud_name,
              [&](std::vector<Eigen::Vector3d>& /*starts*/, std::vector<Eigen::Vector3d>& ends,
                  std::vector<double>& /*times*/, std::vector<ray::RGBA>& /*colours*/) {
                for (auto& e : ends) {
                  bounds_min = bounds_min.cwiseMin(e);
                  bounds_max = bounds_max.cwiseMax(e);
                }
              }, num_bounded, 255.0, nullptr)) {
            std::cerr << "Error: Could not read LAS/LAZ file to determine bounds." << std::endl;
            return false;
        }
        if (bounds_min.x() > bounds_max.x()) {
            std::cerr << "Error: No points found to determine bounds." << std::endl;
            return false;
        }
        user_bounds.min_bound_ = bounds_min;
        user_bounds.max_bound_ = bounds_max;
        std::cout << "Detected bounds min: " << format_vec_string(user_bounds.min_bound_)
                  << ", max: " << format_vec_string(user_bounds.max_bound_) << std::endl;
    } else {
        user_bounds = Cuboid(params.grid_bounds_min, params.grid_bounds_max);
    }

    // --- DTM Setup ---
    std::unique_ptr<HeightField> dtm_ptr;
    if (!params.dtm_file.empty())
    {
        std::cout << "Loading DTM from PLY mesh: " << params.dtm_file << std::endl;
        Mesh dtm_mesh;
        if (!readPlyMesh(params.dtm_file, dtm_mesh))
        {
            std::cerr << "Error: Failed to load DTM mesh file." << std::endl;
            return false;
        }
        dtm_ptr = std::make_unique<HeightField>();
        // Pass the correct DTM cell size for rasterization.
        dtm_ptr->fromMesh(dtm_mesh, user_bounds, params.dtm_cell_size);
    }
    else if (params.dtm_from_class >= 0)
    {
        std::cout << "Generating DTM from class " << params.dtm_from_class << " with cell size " << params.dtm_cell_size << "..." << std::endl;
        std::vector<Eigen::Vector3d> ground_points;

        // Read classification from the passthrough buffer (byte [2] = extended classification).
        uint16_t orig_extra_size = 0;
        std::vector<uint8_t> extra_bytes_vlr;
        readLasExtraBytesVlr(params.cloud_name, orig_extra_size, extra_bytes_vlr);
        const uint16_t stride = static_cast<uint16_t>(kPassthroughStdBytes + orig_extra_size);

        // Chunked pass collecting ground points whose passthrough classification matches.
        size_t num_bounded = 0;
        std::vector<uint8_t> passthrough;
        uint16_t pt_extra = 0;
        if (!ray::readLas(params.cloud_name,
            [&](std::vector<Eigen::Vector3d>& /*starts*/, std::vector<Eigen::Vector3d>& ends,
                std::vector<double>& /*times*/, std::vector<ray::RGBA>& /*colours*/) {
              for (size_t i = 0; i < ends.size(); ++i) {
                const size_t base = i * stride;
                if (passthrough.size() >= base + stride &&
                    passthrough[base + 2] == static_cast<uint8_t>(params.dtm_from_class)) {
                  ground_points.emplace_back(ends[i]);
                }
              }
              passthrough.clear();
            }, num_bounded, 255.0, nullptr, 1000000, nullptr, &passthrough, &pt_extra)) {
            std::cerr << "Error: Could not re-open LAS file to extract ground points for DTM." << std::endl;
            return false;
        }

        if (ground_points.empty()) {
            std::cerr << "Warning: No points found with class " << params.dtm_from_class << ". Cannot generate DTM." << std::endl;
        } else {
            dtm_ptr = std::make_unique<HeightField>();
            dtm_ptr->fromLowestPoint(ground_points, user_bounds, params.dtm_cell_size);
        }
    }

    int padding = (params.neighbour_prior_min_rays > 0) ? 1 : 0;
    Cuboid processing_bounds = user_bounds;
    if (padding > 0) {
        std::cout << "Padding grid by " << padding << " voxel(s) for neighbour prior calculations." << std::endl;
        Eigen::Vector3d padding_vec(padding * params.voxel_size, padding * params.voxel_size, padding * params.voxel_size);
        processing_bounds.min_bound_ -= padding_vec;
        processing_bounds.max_bound_ += padding_vec;
    }

    // Handle auto-detection of reservation size here.
    size_t final_reserve_size = params.reserve_size;
    if (final_reserve_size == 0 && !params.use_ooc) {
        size_t num_bounded = 0;
        ray::readLas(params.cloud_name,
            [&](std::vector<Eigen::Vector3d>& /*starts*/, std::vector<Eigen::Vector3d>& ends,
                std::vector<double>& /*times*/, std::vector<ray::RGBA>& /*colours*/) {
              final_reserve_size += ends.size();
            }, num_bounded, 255.0, nullptr);
    }
    // For OOC, we never reserve in the final grid, as it's populated at the end.
    if (params.use_ooc) { final_reserve_size = 0; }

    std::unique_ptr<VoxelGrid> grid_ptr;
    try {
        grid_ptr = std::make_unique<VoxelGrid>(processing_bounds, params.voxel_size, final_reserve_size);
    } catch (const std::exception& e) {
        std::cerr << "Error during VoxelGrid initialization: " << e.what()
                  << " Not enough memory. Consider using --out_of_core or a larger --voxel_size." << std::endl;
        return false;
    }
    VoxelGrid& grid = *grid_ptr;

    // Pre-calculate peaks if needed. This is done before the main processing.
    if (params.apply_flat_top) {
        std::cout << "Calculating peaks for flat top compensation..." << std::endl;
        calculatePeaks(grid, params.cloud_name);
    }

    // === STRATEGY SELECTION ===
    // This is where we choose which processing strategy to use.
    std::unique_ptr<ProcessingStrategy> strategy;
    if (params.use_ooc) {
        strategy = std::make_unique<OutOfCoreStrategy>(params.num_threads, params.ram_budget_mb);
    } else {
        strategy = std::make_unique<InProcessStrategy>(params.num_threads);
    }

    bool processing_success = strategy->execute(params.cloud_name, grid, params.weighting_method, params.use_occlusion, params.apply_flat_top,
                                                 params.calc_beam_metrics, beam_diameter, beam_divergence, params.subvoxel_split, dtm_ptr.get());

    if (!processing_success) {
      return false; // Strategy failed, exit early.
    }

    // === POST-PROCESSING AND OUTPUT (Same for all strategies) ===
    if (params.neighbour_prior_min_rays > 0) {
        std::cout << "Applying neighbour priors (min rays = " << params.neighbour_prior_min_rays << ")..." << std::endl;
        applyNeighbourPriors(grid, params.neighbour_prior_min_rays);
    }

    // This is the new, decoupled output pipeline.
    // Calculate all output metrics once and store them in a map.
    std::cout << "Calculating output metrics..." << std::endl;
    MetricResultsMap metrics = calculateOutputMetrics(grid, params, dtm_ptr.get());

    // Pass the pre-calculated metrics to the writer functions.
    std::string base_name_stub = getFileNameStub(params.cloud_name);

    // Create a copy of params for the primary output, as it might change write_filled_only_mode
    VoxelizationParameters primary_params = params;
    std::string primary_name_stub = base_name_stub;
    if (primary_params.write_empty_voxels) {
        primary_name_stub += "_include_empty";
    }

    bool primary_success = false;
    std::cout << "Writing primary output file(s)..." << std::endl;
    if (primary_params.output_format == "text") {
        primary_success = writeTextFile(primary_name_stub, grid, metrics, padding, user_bounds, primary_params);
    } else if (primary_params.output_format == "netcdf") {
        primary_success = writeNetcdfFile(primary_name_stub, grid, metrics, padding, user_bounds, primary_params);
    } else if (primary_params.output_format == "amapvox") {
        primary_success = writeAmapVoxFile(primary_name_stub, grid, metrics, padding, user_bounds, primary_params);
    } else {
        std::cerr << "Error: Primary output format '" << primary_params.output_format << "' not supported." << std::endl;
        primary_success = false;
    }

    bool filled_success = true;
    if (params.write_filled) {
        std::cout << "Writing additional 'filled only' output file(s)..." << std::endl;
        std::string filled_name_stub = base_name_stub + "_filled";
        // Create a dedicated params copy for the 'filled' output to control its behavior
        VoxelizationParameters filled_params = params;
        filled_params.write_empty_voxels = false; // 'filled' output is always sparse

        if (filled_params.output_format == "text") {
            filled_success = writeTextFile(filled_name_stub, grid, metrics, padding, user_bounds, filled_params, true);
        } else if (filled_params.output_format == "netcdf") {
            filled_success = writeNetcdfFile(filled_name_stub, grid, metrics, padding, user_bounds, filled_params, true);
        } else if (filled_params.output_format == "amapvox") {
            filled_success = writeAmapVoxFile(filled_name_stub, grid, metrics, padding, user_bounds, filled_params, true);
        }
    }

    bool secondary_success = true;
    if (params.write_amapvox_also && params.output_format != "amapvox") {
        std::cout << "Writing additional AMAPVox output file..." << std::endl;
        secondary_success = writeAmapVoxFile(primary_name_stub, grid, metrics, padding, user_bounds, primary_params);
    }

    return primary_success && filled_success && secondary_success;
}

} // namespace ray
