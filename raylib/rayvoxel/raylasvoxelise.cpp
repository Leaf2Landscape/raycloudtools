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
#include "raylib/raysysinfo.h"
#include "raylib/rayvoxel/rayvox.h"
#include "raylib/rayvoxel/raylasvoxelise.h"
#include "raylib/rayvoxel/raylasvoxelwriter.h"
#include "raylib/rayvoxel/raylasvoxelrefine.h"
#include "raylib/rayvoxel/raylasvegmetrics.h"
#include "raylib/rayvoxel/raylasvoxelprocessor.h"
#include "raylib/rayvoxel/raylasthreadsafequeue.h"
#include "raylib/rayvoxel/raylasheightfield.h"
#include "raylib/rayvoxel/raylasbinaryio.h"
#include "raylib/rayvoxel/raylasbailey.h"
#include "raylib/rayvoxel/raylasvegmetrics.h"

#include <nabo/nabo.h>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cctype>
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
#include <set>

namespace ray
{

// The passthrough buffer produced by readLas has 10 standard LAS field bytes per point
// before any original sensor extra bytes. Field byte [0] packs return number (low nibble)
// and number of returns (high nibble); byte [2] is the extended classification.
static const uint16_t kPassthroughStdBytes = 10;

// ==================================================================================
// VoxelGrid Class Implementation
// ==================================================================================

VoxelGrid::VoxelGrid(const Cuboid &grid_bounds, double vox_width,
                     size_t ram_budget_bytes, size_t sparse_reservation)
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

  const int64_t total = voxel_dims_[0] * voxel_dims_[1] * voxel_dims_[2];
  const size_t required = static_cast<size_t>(total) * sizeof(Voxel);

  if (total > 0 && required / sizeof(Voxel) == static_cast<size_t>(total) && required <= ram_budget_bytes) {
    flat_voxels_.assign(static_cast<size_t>(total), Voxel{});
    use_sparse_fallback_ = false;
    std::cout << "Initialized VoxelGrid: " << voxel_dims_.transpose()
              << " flat array (" << required / (1024 * 1024) << " MB)" << std::endl;
  } else {
    use_sparse_fallback_ = true;
    if (sparse_reservation > 0)
      sparse_voxels_.reserve(sparse_reservation);
    std::cout << "Initialized VoxelGrid: " << voxel_dims_.transpose()
              << " sparse map (flat would need " << required / (1024 * 1024) << " MB > budget "
              << ram_budget_bytes / (1024 * 1024) << " MB)" << std::endl;
  }
}

void VoxelGrid::merge(const VoxelProcessor& processor)
{
  std::lock_guard<std::mutex> lock(merge_mutex_);
  if (!use_sparse_fallback_) {
    for (const auto& pair : processor.getMap())
      flat_voxels_[flatIndex(pair.first.x, pair.first.y, pair.first.z)] += pair.second;
  } else {
    for (const auto& pair : processor.getMap())
      sparse_voxels_[pair.first] += pair.second;
  }
}

void VoxelGrid::take(VoxelProcessor& processor)
{
  absorbMap(processor.takeMap());
}

void VoxelGrid::absorbMap(VoxelProcessor::Map&& m)
{
  if (!use_sparse_fallback_) {
    for (auto& pair : m)
      flat_voxels_[flatIndex(pair.first.x, pair.first.y, pair.first.z)] += pair.second;
  } else {
    for (auto& pair : m)
      sparse_voxels_[pair.first] += pair.second;
  }
}

VoxelGrid::VoxelState VoxelGrid::getVoxelState(int64_t i, int64_t j, int64_t k) const {
    const Voxel& v = getVoxel(i, j, k);
    if (v.num_hits > 0) return VoxelState::FILLED;
    if (v.num_beams_weighted > 0.0f) return VoxelState::EMPTY;
    if (v.num_rays_occluded > 0.0f) return VoxelState::OCCLUDED;
    return VoxelState::UNOBSERVED;
}

const VoxelGrid::Voxel& VoxelGrid::getVoxel(int64_t i, int64_t j, int64_t k) const {
    static const Voxel empty_voxel;
    if (!use_sparse_fallback_) {
        return flat_voxels_[flatIndex(i, j, k)];
    }
    auto it = sparse_voxels_.find({i, j, k});
    return it != sparse_voxels_.end() ? it->second : empty_voxel;
}

int64_t VoxelGrid::getIndex(int64_t i, int64_t j, int64_t k) const {
    return i + j * voxel_dims_[0] + k * voxel_dims_[0] * voxel_dims_[1];
}

// ==================================================================================
// Voxel Method Implementations
// ==================================================================================

double VoxelGrid::Voxel::pad_g0_5() const
{
  // Bias-corrected MLE assuming spherical LAD (G=0.5): PAD = 2*(N-1)/N * H / L_obs.
  const double eps = 1e-10;
  if (num_beams_weighted < 2.0f) return 0.0;
  return 2.0 * (num_beams_weighted - 1.0f) * num_hits / (eps + num_beams_weighted * path_length);
}

double VoxelGrid::Voxel::transmittance() const
{
  const double eps = 1e-10;
  if (bs_free_path > eps)
    return std::exp(-static_cast<double>(bs_intercepted) / static_cast<double>(bs_free_path));
  if (bs_entering > eps) {
    double transmitted = bs_entering - bs_intercepted;
    return std::max(0.0, transmitted / static_cast<double>(bs_entering));
  }
  return 1.0;
}

// ==================================================================================
// Local Helpers
// ==================================================================================

namespace {

// Beams are queued in fixed-size batches to amortise mutex overhead.
// 32 beams/batch → 32× fewer lock acquisitions vs one-beam-per-slot.
static constexpr size_t kBeamBatchSize = 32;
struct BeamBatch {
  std::array<BeamData, kBeamBatchSize> beams;
  size_t count = 0;
};

// Build a PointData record from a single ray (start -> end) and its passthrough bytes.
// The stride is the per-point passthrough byte count (10 standard fields + sensor extras).
// distance_to_sensor is computed as the ray length (end - start).norm().
// @c alpha is the decoded per-point intensity/alpha: alpha==0 marks an unbound (miss) ray.
PointData makePointData(const Eigen::Vector3d& start, const Eigen::Vector3d& end,
                        double gps_time, int32_t beam_id, uint8_t alpha,
                        const std::vector<uint8_t>& passthrough, size_t index, uint16_t stride)
{
  PointData pd;
  pd.beam_origin       = start;                       // sx,sy,sz already decoded by readLas
  pd.x                 = end.x();
  pd.y                 = end.y();
  pd.z                 = end.z();
  pd.gps_time          = gps_time;
  pd.beam_id           = beam_id;
  pd.intensity         = alpha;                       // per-point intensity (alpha channel)
  pd.bound             = (alpha > 0) ? 1 : 0;         // alpha==0 is an unbound (miss) ray
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

// Decide whether the incoming point pd starts a new beam relative to the pending return group.
// has_beam_ids selects the authoritative beam-id grouping; otherwise points are grouped by GPS
// time, number-of-returns consistency, and monotonically increasing distance-to-sensor.
static bool isNewBeam(const PointData& pd,
                      double pending_gps_time,
                      int32_t pending_beam_id,
                      const std::vector<PointData>& pending_returns,
                      bool has_beam_ids)
{
    if (pending_returns.empty()) return true;
    if (has_beam_ids)
        return pd.beam_id != pending_beam_id
            || pd.number_of_returns != pending_returns.front().number_of_returns
            || pd.distance_to_sensor <= pending_returns.back().distance_to_sensor;

    // nor <= 1: each point is its own beam by definition
    if (pd.number_of_returns <= 1) return true;
    if (pending_returns.front().number_of_returns <= 1) return true;

    // nor > 1: group by GPS time + nor consistency + monotonic distance
    return pd.gps_time != pending_gps_time
        || pd.number_of_returns != pending_returns.front().number_of_returns
        || pd.distance_to_sensor <= pending_returns.back().distance_to_sensor;
}

// Returns true if a point is a ground hit and should be excluded from PAD/LAD/WAD.
// Two mutually exclusive modes: dtm_from_class >= 0 selects a class-match path; otherwise a
// valid DTM mesh selects the vertical-distance path (point above DTM within dtm_filter_distance).
static bool isGroundHit(double x, double y, double z, uint8_t classification,
                         int dtm_from_class, const HeightField* dtm, double dtm_filter_distance)
{
  if (dtm_from_class >= 0)
    return classification == static_cast<uint8_t>(dtm_from_class);
  if (dtm && dtm->isValid() && dtm_filter_distance > 0.0) {
    double ground_h;
    if (dtm->getHeightNearest(x, y, ground_h)) {
      const double above = z - ground_h;
      return above >= 0.0 && above <= dtm_filter_distance;
    }
  }
  return false;
}

// Byte size of each LAS extra-bytes data_type, indexed by the ASPRS data_type code.
// Identical to kExtraTypeSize in raylaz.cpp (duplicated here since it is file-local there).
static const uint8_t kExtraByteSizes[11] = { 0, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };

// Describes where a classification value lives inside a per-point passthrough record:
// a byte offset into the record and the LAS data_type used to decode it.
// Defaults select the standard Classification byte (offset 2, dtype 1 = u8).
struct ClassFieldSource { uint16_t byte_offset = 2; uint8_t las_dtype = 1; };

// Resolve a named LAS extra-byte field to its passthrough byte offset and data_type.
// An empty name (or "classification", case-insensitive) selects the standard
// Classification byte. The extra_bytes_vlr blob is a sequence of 192-byte ASPRS
// extra-bytes records (byte [2] = data_type, bytes [4..35] = null-terminated name).
// The passthrough buffer stores the sensor extra bytes after kPassthroughStdBytes
// standard bytes, in VLR order; each preceding field advances the cumulative offset
// by its data_type size. On no match, prints one warning and falls back to defaults.
ClassFieldSource resolveClassField(const std::string& field_name,
                                   const std::vector<uint8_t>& extra_bytes_vlr)
{
  if (field_name.empty()) return ClassFieldSource{};
  std::string lowered = field_name;
  for (char& c : lowered) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  if (lowered == "classification") return ClassFieldSource{};

  const size_t num_records = extra_bytes_vlr.size() / 192;
  uint16_t cumulative_offset = 0;
  for (size_t i = 0; i < num_records; ++i)
  {
    const uint8_t dtype = extra_bytes_vlr[i * 192 + 2];
    char name[33] = {};
    std::memcpy(name, &extra_bytes_vlr[i * 192 + 4], 32);
    if (field_name == name)
      return ClassFieldSource{ static_cast<uint16_t>(kPassthroughStdBytes + cumulative_offset), dtype };
    if (dtype > 0 && dtype <= 10) cumulative_offset += kExtraByteSizes[dtype];
  }

  std::cerr << "Warning: classification field '" << field_name
            << "' not found in extra bytes; using standard Classification byte." << std::endl;
  return ClassFieldSource{};
}

// Decode the integer class value at base + src.byte_offset according to src.las_dtype.
int readClassValue(const uint8_t* base, const ClassFieldSource& src)
{
  const uint8_t* p = base + src.byte_offset;
  switch (src.las_dtype)
  {
    case 1: return static_cast<int>(*p);                                            // u8
    case 2: return static_cast<int>(static_cast<int8_t>(*p));                       // i8
    case 3: return static_cast<int>(*reinterpret_cast<const uint16_t*>(p));         // u16
    case 4: return static_cast<int>(*reinterpret_cast<const int16_t*>(p));          // i16
    case 5: {                                                                       // u32
      const uint32_t v = *reinterpret_cast<const uint32_t*>(p);
      return v > static_cast<uint32_t>(std::numeric_limits<int>::max())
               ? std::numeric_limits<int>::max() : static_cast<int>(v);
    }
    case 6: return *reinterpret_cast<const int32_t*>(p);                            // i32
    default: return static_cast<int>(*p);                                          // u8 fallback
  }
}

} // anonymous namespace

// ==================================================================================
// Classification Post-Traversal Pass
// ==================================================================================

// Reads the cloud endpoints only (no ray walking) to build a flat-index→per-class
// hit count table. O(N_points) I/O pass, much cheaper than the traversal pass.
static ClassTable buildClassTable(const std::string& cloud_name, const VoxelGrid& grid,
                                  int dtm_from_class, const HeightField* dtm, double dtm_filter_distance)
{
  ClassTable class_table;

  uint16_t orig_extra_size = 0;
  std::vector<uint8_t> extra_bytes_vlr;
  readLasExtraBytesVlr(cloud_name, orig_extra_size, extra_bytes_vlr);
  const uint16_t stride = static_cast<uint16_t>(kPassthroughStdBytes + orig_extra_size);

  const Cuboid& bounds    = grid.getBounds();
  const double vox_width  = grid.getVoxelWidth();
  const auto& dims        = grid.getDimensions();

  size_t num_bounded = 0;
  std::vector<uint8_t> passthrough;
  uint16_t pt_extra = 0;
  // global_chunk_start tracks the running point offset across readLas chunk callbacks.
  // The mmap fast path pre-allocates the passthrough for all points at global indices, while
  // the sequential path appends per-chunk — both are correct when indexing with the global offset.
  size_t global_chunk_start = 0;

  ray::readLas(cloud_name,
    [&](std::vector<Eigen::Vector3d>& /*starts*/, std::vector<Eigen::Vector3d>& ends,
        std::vector<double>& /*times*/, std::vector<ray::RGBA>& colours) {
      for (size_t i = 0; i < ends.size(); ++i) {
        const size_t base = (global_chunk_start + i) * stride;
        if (passthrough.size() < base + stride) continue;
        // "bound" is encoded as bound==(alpha>0), so alpha==0 identifies an unbound ray for both
        // new files (authoritative bound field) and old files (alpha-only fallback). Unbound rays
        // are never counted as classified hits.
        const uint8_t alpha = (i < colours.size()) ? colours[i].alpha : 1;
        if (alpha == 0) continue;

        const uint8_t classification = passthrough[base + 2];
        if (isGroundHit(ends[i].x(), ends[i].y(), ends[i].z(), classification,
                        dtm_from_class, dtm, dtm_filter_distance)) continue;
        const Eigen::Vector3d vox = (ends[i] - bounds.min_bound_) / vox_width;
        const int64_t ix = static_cast<int64_t>(std::floor(vox.x()));
        const int64_t iy = static_cast<int64_t>(std::floor(vox.y()));
        const int64_t iz = static_cast<int64_t>(std::floor(vox.z()));
        if (ix < 0 || ix >= dims[0] || iy < 0 || iy >= dims[1] || iz < 0 || iz >= dims[2]) continue;
        class_table[grid.flatIndex(ix, iy, iz)][classification] += 1.0f;
      }
      global_chunk_start += ends.size();
    }, num_bounded, 255.0, nullptr, 1000000, nullptr, &passthrough, &pt_extra);

  return class_table;
}

// Per-tile accumulator for the tiled parallel KNN/IAD pass. Each worker thread owns one
// TileResult; histograms are merged serially after the thread pool joins.
struct TileResult {
  std::unordered_map<int64_t, std::vector<double>> all_hist, leaf_hist, wood_hist, beam_hist;
  std::unordered_map<int64_t, float> leaf_hit_count, wood_hit_count;
  std::unordered_map<int64_t, TriangleHistograms> triangle_histograms;
};

// Per-point surface normals estimated via KNN PCA (smallest-eigenvalue eigenvector).
// Inclination angle theta = acos(|n_z|) is binned over [0, pi/2] into LIAD/WIAD/PIAD.
// Angle-integrated G_eff is the projection kernel A(theta_beam, theta_L) weighted over
// the empirical beam-direction distribution — more accurate than single-angle evaluation.
//
// Vicari, M.B., Pisek, J. & Disney, M. (2019). New estimates of leaf angle distribution
// from terrestrial LiDAR: Comparison with measured and modelled estimates from nine
// broadleaf tree species. Agricultural and Forest Meteorology, 264, 322-333.
// DOI: 10.1016/j.agrformet.2018.10.021
//
// NOTE: this is the leaf-ANGLE paper (AgForMet). Do not confuse with Vicari et al. (2019)
// Methods Ecol. Evol. 10(5):680-694 (DOI 10.1111/2041-210X.13144), which covers
// leaf/wood point SEPARATION and is a distinct method.
static void buildClassAndIadTable(const std::string& cloud_name, const VoxelGrid& grid,
                                  const VoxelizationParameters& params, const HeightField* dtm,
                                  ClassTable& class_table_out, IadTable& iad_table_out)
{
  ClassTable class_table;
  IadTable iad_table;

  uint16_t orig_extra_size = 0;
  std::vector<uint8_t> extra_bytes_vlr;
  readLasExtraBytesVlr(cloud_name, orig_extra_size, extra_bytes_vlr);
  const uint16_t stride = static_cast<uint16_t>(kPassthroughStdBytes + orig_extra_size);

  const Cuboid& bounds    = grid.getBounds();
  const double vox_width  = grid.getVoxelWidth();
  const auto& dims        = grid.getDimensions();

  // Resolve the leaf/wood classification source fields. Syntax: "[<field>:]c1,c2,...".
  // An optional "<field>:" prefix selects a named LAS extra-byte field; without it the
  // standard Classification byte is used. Leaf and wood resolve independently.
  auto split_field_codes = [](const std::string& s, std::string& field, std::string& codes) {
    auto colon = s.find(':');
    if (colon != std::string::npos) { field = s.substr(0, colon); codes = s.substr(colon + 1); }
    else { field.clear(); codes = s; }
  };
  std::string leaf_field, leaf_codes_str, wood_field, wood_codes_str;
  split_field_codes(params.leaf_classes_str, leaf_field, leaf_codes_str);
  split_field_codes(params.wood_classes_str, wood_field, wood_codes_str);

  const ClassFieldSource leaf_src = resolveClassField(leaf_field, extra_bytes_vlr);
  const ClassFieldSource wood_src = resolveClassField(wood_field, extra_bytes_vlr);

  // Collect endpoints, classifications and flat voxel indices (same filtering as buildClassTable).
  // leaf_vals/wood_vals hold the class value read from each point's resolved field (these may
  // come from different fields, hence two separate vectors rather than one shared `classes`).
  std::vector<Eigen::Vector3d> positions;
  std::vector<int> leaf_vals;
  std::vector<int> wood_vals;
  std::vector<int64_t> flat_indices;
  std::vector<double> beam_angles;  // zenith angle [0, pi/2] of each point's inbound ray

  size_t num_bounded = 0;
  std::vector<uint8_t> passthrough;
  uint16_t pt_extra = 0;
  size_t global_chunk_start = 0;

  ray::readLas(cloud_name,
    [&](std::vector<Eigen::Vector3d>& starts, std::vector<Eigen::Vector3d>& ends,
        std::vector<double>& /*times*/, std::vector<ray::RGBA>& colours) {
      for (size_t i = 0; i < ends.size(); ++i) {
        const size_t base = (global_chunk_start + i) * stride;
        if (passthrough.size() < base + stride) continue;
        // "bound" is encoded as bound==(alpha>0), so alpha==0 identifies an unbound ray for both
        // new files (authoritative bound field) and old files (alpha-only fallback). Unbound rays
        // are never counted as classified hits.
        const uint8_t alpha = (i < colours.size()) ? colours[i].alpha : 1;
        if (alpha == 0) continue;

        const Eigen::Vector3d vox = (ends[i] - bounds.min_bound_) / vox_width;
        const int64_t ix = static_cast<int64_t>(std::floor(vox.x()));
        const int64_t iy = static_cast<int64_t>(std::floor(vox.y()));
        const int64_t iz = static_cast<int64_t>(std::floor(vox.z()));
        if (ix < 0 || ix >= dims[0] || iy < 0 || iy >= dims[1] || iz < 0 || iz >= dims[2]) continue;

        // ClassTable work (same as buildClassTable): per-voxel standard-classification counts.
        const uint8_t classification = passthrough[base + 2];
        if (isGroundHit(ends[i].x(), ends[i].y(), ends[i].z(), classification,
                        params.dtm_from_class, dtm, params.dtm_filter_distance)) continue;
        const int64_t flat_idx = grid.flatIndex(ix, iy, iz);
        class_table[flat_idx][classification] += 1.0f;

        // IadTable collection (same as buildIadTable): leaf/wood field values + flat index.
        const int lv = readClassValue(&passthrough[base], leaf_src);
        const int wv = readClassValue(&passthrough[base], wood_src);
        const Eigen::Vector3d dir = ends[i] - starts[i];
        const double len2 = dir.squaredNorm();
        const double bz = (len2 > 1e-12) ? std::acos(std::min(1.0, std::abs(dir.z() / std::sqrt(len2)))) : 0.0;
        positions.push_back(ends[i]);
        leaf_vals.push_back(lv);  // preserve sign: -1 means "neither", must not be clamped to 0
        wood_vals.push_back(wv);
        flat_indices.push_back(flat_idx);
        beam_angles.push_back(bz);
      }
      global_chunk_start += ends.size();
    }, num_bounded, 255.0, nullptr, 1000000, nullptr, &passthrough, &pt_extra);

  if (positions.size() < 2) {
    class_table_out = std::move(class_table);
    iad_table_out = std::move(iad_table);
    return;
  }

  // Parse leaf/wood class sets (read-only during tile workers).
  std::set<int> leaf_set, wood_set;
  {
    std::stringstream ss(leaf_codes_str);
    std::string item;
    while (std::getline(ss, item, ',')) { try { leaf_set.insert(std::stoi(item)); } catch (...) {} }
  }
  {
    std::stringstream ss(wood_codes_str);
    std::string item;
    while (std::getline(ss, item, ',')) { try { wood_set.insert(std::stoi(item)); } catch (...) {} }
  }

  const bool any_bailey = std::any_of(params.attenuation_methods.begin(), params.attenuation_methods.end(),
                                       [](const std::string& m){ return m == "bailey"; });

  size_t resolved_threads = std::thread::hardware_concurrency();
  if (resolved_threads == 0) resolved_threads = 1;

  // buf_m must be >= tile so 3x3 neighbour scan covers full buffer.
  const double buf_m = 1.0;
  const double tile_sz = std::max(buf_m, params.iad_tile_size);
  const double minx = bounds.min_bound_.x(), miny = bounds.min_bound_.y();
  const double maxx = bounds.max_bound_.x(), maxy = bounds.max_bound_.y();
  const int n_tx = std::max(1, (int)std::ceil((maxx - minx) / tile_sz));
  const int n_ty = std::max(1, (int)std::ceil((maxy - miny) / tile_sz));
  const int n_tiles = n_tx * n_ty;
  std::unordered_map<int64_t, std::vector<double>> all_hist, leaf_hist, wood_hist, beam_hist;
  std::unordered_map<int64_t, float> leaf_hit_count, wood_hit_count;
  std::unordered_map<int64_t, TriangleHistograms> triangle_histograms;

  if (n_tiles == 1) {
    // Small cloud: global single KD-tree, no tiling overhead.
    const int K = std::min(params.knn_normal, (int)positions.size() - 1);
    Eigen::MatrixXd points_p(3, positions.size());
    for (size_t i = 0; i < positions.size(); ++i) points_p.col(i) = positions[i];
    std::unique_ptr<Nabo::NNSearchD> nns(Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3));
    Eigen::MatrixXi indices(K, (int)positions.size());
    Eigen::MatrixXd dists2(K, (int)positions.size());
    nns->knn(points_p, indices, dists2, K, kNearestNeighbourEpsilon, 0);
    nns.reset(nullptr);

    if (any_bailey) {
      std::vector<int> class_labels_int(positions.size(), 0);
      for (size_t i = 0; i < positions.size(); ++i) {
        if (leaf_set.count(leaf_vals[i])) class_labels_int[i] = 1;
        else if (wood_set.count(wood_vals[i])) class_labels_int[i] = -1;
      }
      triangle_histograms = buildTriangleInclinationHistograms(
          positions, indices, flat_indices, class_labels_int, params.n_iad_bins, params.triangle_lmax);
    }

    for (size_t i = 0; i < positions.size(); ++i) {
      Eigen::Vector3d centroid(0, 0, 0);
      int num_neighbours = 0;
      for (int j = 0; j < K && indices(j, i) != Nabo::NNSearchD::InvalidIndex; ++j) {
        centroid += positions[indices(j, i)];
        ++num_neighbours;
      }
      if (num_neighbours < 3) continue;
      centroid /= static_cast<double>(num_neighbours);
      Eigen::Matrix3d scatter = Eigen::Matrix3d::Zero();
      for (int j = 0; j < K && indices(j, i) != Nabo::NNSearchD::InvalidIndex; ++j) {
        Eigen::Vector3d offset = positions[indices(j, i)] - centroid;
        scatter += offset * offset.transpose();
      }
      scatter /= static_cast<double>(num_neighbours);

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter);
      const Eigen::Vector3d normal = eigen_solver.eigenvectors().col(0);
      const double theta = std::acos(std::min(1.0, std::abs(normal.z())));

      int bin = static_cast<int>(theta / (kPi / 2.0) * params.n_iad_bins);
      bin = std::clamp(bin, 0, params.n_iad_bins - 1);

      const int64_t flat_idx = flat_indices[i];
      auto& ah = all_hist[flat_idx];
      if (ah.empty()) ah.assign(params.n_iad_bins, 0.0);
      ah[bin] += 1.0;

      {
        int bbin = static_cast<int>(beam_angles[i] / (kPi / 2.0) * params.n_iad_bins);
        bbin = std::clamp(bbin, 0, params.n_iad_bins - 1);
        auto& bh = beam_hist[flat_idx];
        if (bh.empty()) bh.assign(params.n_iad_bins, 0.0);
        bh[bbin] += 1.0;
      }

      if (leaf_set.count(leaf_vals[i])) {
        leaf_hit_count[flat_idx] += 1.0f;
        auto& lh = leaf_hist[flat_idx];
        if (lh.empty()) lh.assign(params.n_iad_bins, 0.0);
        lh[bin] += 1.0;
      }
      if (wood_set.count(wood_vals[i])) {
        wood_hit_count[flat_idx] += 1.0f;
        auto& wh = wood_hist[flat_idx];
        if (wh.empty()) wh.assign(params.n_iad_bins, 0.0);
        wh[bin] += 1.0;
      }
    }
  } else {
  // Large cloud: tiled parallel KNN.
  // Bucket each point by its core tile.
  std::vector<std::vector<size_t>> core_points(n_tiles);
  std::vector<int> pt_tile(positions.size());
  for (size_t i = 0; i < positions.size(); ++i) {
    int tx = std::clamp((int)((positions[i].x() - minx) / tile_sz), 0, n_tx - 1);
    int ty = std::clamp((int)((positions[i].y() - miny) / tile_sz), 0, n_ty - 1);
    int t  = ty * n_tx + tx;
    pt_tile[i] = t;
    core_points[t].push_back(i);
  }

  // For bailey: each voxel's triangle histogram is owned by a single tile to avoid
  // double-counting (TriangleHistograms holds area-weighted means, not summable counts).
  // Owner tile = tile of the lowest global point index that maps to that flat_idx.
  std::unordered_map<int64_t, int> flat_owner;
  if (any_bailey) {
    for (size_t i = 0; i < positions.size(); ++i) {
      auto it = flat_owner.find(flat_indices[i]);
      if (it == flat_owner.end()) flat_owner[flat_indices[i]] = pt_tile[i];
    }
  }

  std::vector<TileResult> results(resolved_threads);
  std::atomic<int> next_tile(0);

  auto worker = [&](size_t w) {
    TileResult& tr = results[w];
    int t;
    while ((t = next_tile.fetch_add(1)) < n_tiles) {
      if (core_points[t].empty()) continue;

      const int tx = t % n_tx;
      const int ty = t / n_tx;
      const double cx0 = minx + tx * tile_sz;
      const double cx1 = std::min(maxx, cx0 + tile_sz);
      const double cy0 = miny + ty * tile_sz;
      const double cy1 = std::min(maxy, cy0 + tile_sz);

      // Collect buffered point indices from 3x3 tile neighbourhood.
      std::vector<size_t> buf;
      for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx) {
          int nx2 = tx + dx, ny2 = ty + dy;
          if (nx2 < 0 || nx2 >= n_tx || ny2 < 0 || ny2 >= n_ty) continue;
          for (size_t gi : core_points[ny2 * n_tx + nx2]) {
            const Eigen::Vector3d& p = positions[gi];
            if (p.x() >= cx0 - buf_m && p.x() <= cx1 + buf_m &&
                p.y() >= cy0 - buf_m && p.y() <= cy1 + buf_m)
              buf.push_back(gi);
          }
        }

      const size_t Nb = buf.size();
      if (Nb < 3) continue;

      const int K = std::min(params.knn_normal, (int)Nb - 1);

      // Build tile-local KD-tree over buffered points.
      Eigen::MatrixXd pts(3, Nb);
      for (size_t c = 0; c < Nb; ++c) pts.col(c) = positions[buf[c]];
      std::unique_ptr<Nabo::NNSearchD> nns(Nabo::NNSearchD::createKDTreeLinearHeap(pts, 3));
      Eigen::MatrixXi idx(K, Nb);
      Eigen::MatrixXd d2(K, Nb);
      nns->knn(pts, idx, d2, K, kNearestNeighbourEpsilon, 0);
      nns.reset(nullptr);

      // PCA normal estimation — accumulate only for core points.
      for (size_t c = 0; c < Nb; ++c) {
        if (pt_tile[buf[c]] != t) continue;  // skip buffer-only points

        Eigen::Vector3d centroid(0, 0, 0);
        int num_neighbours = 0;
        for (int j = 0; j < K && idx(j, c) != Nabo::NNSearchD::InvalidIndex; ++j) {
          centroid += positions[buf[idx(j, c)]];
          ++num_neighbours;
        }
        if (num_neighbours < 3) continue;
        centroid /= static_cast<double>(num_neighbours);
        Eigen::Matrix3d scatter = Eigen::Matrix3d::Zero();
        for (int j = 0; j < K && idx(j, c) != Nabo::NNSearchD::InvalidIndex; ++j) {
          Eigen::Vector3d offset = positions[buf[idx(j, c)]] - centroid;
          scatter += offset * offset.transpose();
        }
        scatter /= static_cast<double>(num_neighbours);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter);
        const Eigen::Vector3d normal = eigen_solver.eigenvectors().col(0);
        const double theta = std::acos(std::min(1.0, std::abs(normal.z())));

        int bin = static_cast<int>(theta / (kPi / 2.0) * params.n_iad_bins);
        bin = std::clamp(bin, 0, params.n_iad_bins - 1);

        const int64_t flat_idx = flat_indices[buf[c]];
        auto& ah = tr.all_hist[flat_idx];
        if (ah.empty()) ah.assign(params.n_iad_bins, 0.0);
        ah[bin] += 1.0;

        {
          int bbin = static_cast<int>(beam_angles[buf[c]] / (kPi / 2.0) * params.n_iad_bins);
          bbin = std::clamp(bbin, 0, params.n_iad_bins - 1);
          auto& bh = tr.beam_hist[flat_idx];
          if (bh.empty()) bh.assign(params.n_iad_bins, 0.0);
          bh[bbin] += 1.0;
        }

        if (leaf_set.count(leaf_vals[buf[c]])) {
          tr.leaf_hit_count[flat_idx] += 1.0f;
          auto& lh = tr.leaf_hist[flat_idx];
          if (lh.empty()) lh.assign(params.n_iad_bins, 0.0);
          lh[bin] += 1.0;
        }
        if (wood_set.count(wood_vals[buf[c]])) {
          tr.wood_hit_count[flat_idx] += 1.0f;
          auto& wh = tr.wood_hist[flat_idx];
          if (wh.empty()) wh.assign(params.n_iad_bins, 0.0);
          wh[bin] += 1.0;
        }
      }

      // Bailey triangle-facet histograms (single-owner per voxel to avoid
      // double-counting area-weighted means across tile boundaries).
      if (any_bailey) {
        std::vector<Eigen::Vector3d> local_pos(Nb);
        std::vector<int64_t> local_flat(Nb);
        std::vector<int> labels(Nb, 0);
        for (size_t c = 0; c < Nb; ++c) {
          local_pos[c]  = positions[buf[c]];
          local_flat[c] = flat_indices[buf[c]];
          if (leaf_set.count(leaf_vals[buf[c]])) labels[c] = 1;
          else if (wood_set.count(wood_vals[buf[c]])) labels[c] = -1;
        }
        auto th = buildTriangleInclinationHistograms(local_pos, idx, local_flat, labels,
                                                     params.n_iad_bins, params.triangle_lmax);
        for (auto& kv : th) {
          auto oit = flat_owner.find(kv.first);
          if (oit != flat_owner.end() && oit->second == t)
            tr.triangle_histograms.emplace(kv.first, std::move(kv.second));
        }
      }
    }
  };

  const size_t n_workers = std::min(resolved_threads, (size_t)std::max(1, n_tiles));
  std::vector<std::thread> pool;
  pool.reserve(n_workers);
  for (size_t w = 0; w < n_workers; ++w) pool.emplace_back(worker, w);
  for (auto& th : pool) th.join();

  // Merge per-tile results. Histograms are additive (each point is accumulated by exactly
  // one core tile). Triangle histograms use single-owner assignment (non-additive means).
  auto merge_hist = [&](std::unordered_map<int64_t, std::vector<double>>& dst,
                        std::unordered_map<int64_t, std::vector<double>>& src) {
    for (auto& kv : src) {
      auto& d = dst[kv.first];
      if (d.empty()) d = std::move(kv.second);
      else for (int b = 0; b < params.n_iad_bins; ++b) d[b] += kv.second[b];
    }
  };

  for (auto& tr : results) {
    merge_hist(all_hist, tr.all_hist);
    merge_hist(leaf_hist, tr.leaf_hist);
    merge_hist(wood_hist, tr.wood_hist);
    merge_hist(beam_hist, tr.beam_hist);
    for (auto& kv : tr.leaf_hit_count) leaf_hit_count[kv.first] += kv.second;
    for (auto& kv : tr.wood_hit_count) wood_hit_count[kv.first] += kv.second;
    for (auto& kv : tr.triangle_histograms) triangle_histograms.emplace(kv.first, std::move(kv.second));
  }
  } // end tiled path

  // L1-normalize each histogram (leave all-zero if its sum is zero).
  auto normalize = [](std::vector<double>& h) {
    double sum = 0.0;
    for (double v : h) sum += v;
    if (sum > 0.0) for (double& v : h) v /= sum;
  };

  std::vector<double> bin_centres(params.n_iad_bins);
  for (int b = 0; b < params.n_iad_bins; ++b)
    bin_centres[b] = (b + 0.5) * (kPi / 2.0) / params.n_iad_bins;

  for (auto& pair : all_hist) {
    const int64_t flat_idx = pair.first;
    IadData iad;
    iad.bin_centres = bin_centres;
    iad.piad = pair.second;  // all points -> plant
    auto lit = leaf_hist.find(flat_idx);
    iad.liad = (lit != leaf_hist.end()) ? lit->second : std::vector<double>(params.n_iad_bins, 0.0);
    auto wit = wood_hist.find(flat_idx);
    iad.wiad = (wit != wood_hist.end()) ? wit->second : std::vector<double>(params.n_iad_bins, 0.0);

    // For the bailey method, extract triangle-facet G values without overwriting the empirical
    // liad/wiad/piad. This preserves empirical G (leaf_g/wood_g/plant_g) for non-bailey methods
    // that may be running alongside bailey in the same invocation.
    if (any_bailey) {
      auto it = triangle_histograms.find(flat_idx);
      if (it != triangle_histograms.end()) {
        iad.bailey_g_leaf = it->second.bailey_g_leaf;
        iad.bailey_g_wood = it->second.bailey_g_wood;
      }
    }

    normalize(iad.liad);
    normalize(iad.wiad);
    normalize(iad.piad);

    // Angle-integrated G: weight G(theta_beam, leaf_angles) over the empirical beam-direction
    // distribution rather than evaluating at a single mean angle.
    auto bhit = beam_hist.find(flat_idx);
    if (bhit != beam_hist.end()) {
      std::vector<double> norm_beam = bhit->second;
      normalize(norm_beam);
      double g_plant = 0.0, g_leaf = 0.0, g_wood = 0.0;
      for (int b = 0; b < params.n_iad_bins; ++b) {
        if (norm_beam[b] <= 0.0) continue;
        g_plant += norm_beam[b] * computeGFromHistogram(bin_centres[b], iad.bin_centres, iad.piad);
        g_leaf  += norm_beam[b] * computeGFromHistogram(bin_centres[b], iad.bin_centres, iad.liad);
        g_wood  += norm_beam[b] * computeGFromHistogram(bin_centres[b], iad.bin_centres, iad.wiad);
      }
      iad.plant_g = g_plant;
      iad.leaf_g  = g_leaf;
      iad.wood_g  = g_wood;
    } else {
      const int64_t ci = flat_idx % dims[0];
      const int64_t cj = (flat_idx / dims[0]) % dims[1];
      const int64_t ck = flat_idx / (dims[0] * dims[1]);
      const VoxelGrid::Voxel& vv = grid.getVoxel(ci, cj, ck);
      const double mean_zenith = (vv.num_beams_weighted > 0) ? (vv.sum_of_angles / vv.num_beams_weighted) : 0.0;
      iad.plant_g = computeGFromHistogram(mean_zenith, iad.bin_centres, iad.piad);
      iad.leaf_g  = computeGFromHistogram(mean_zenith, iad.bin_centres, iad.liad);
      iad.wood_g  = computeGFromHistogram(mean_zenith, iad.bin_centres, iad.wiad);
    }

    {
      auto lhit = leaf_hit_count.find(flat_idx);
      if (lhit != leaf_hit_count.end()) iad.leaf_hits = lhit->second;
      auto whit = wood_hit_count.find(flat_idx);
      if (whit != wood_hit_count.end()) iad.wood_hits = whit->second;
    }

    iad_table[flat_idx] = std::move(iad);
  }

  class_table_out = std::move(class_table);
  iad_table_out = std::move(iad_table);
}

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
                       const HeightField* dtm, int dtm_from_class, double dtm_filter_distance, double lambda1 = 0.0) = 0;
};

// --- In-Memory Strategy (Options 1 & 2) ---
class InProcessStrategy : public ProcessingStrategy {
public:
  InProcessStrategy(size_t num_threads) : num_threads_(num_threads) {}
  bool execute(const std::string& cloud_name, VoxelGrid& grid,
               const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
               bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
               const HeightField* dtm, int dtm_from_class, double dtm_filter_distance, double lambda1 = 0.0) override;
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
               const HeightField* dtm, int dtm_from_class, double dtm_filter_distance, double lambda1 = 0.0) override;
private:
  bool createShards(const std::string& cloud_name, VoxelGrid& grid,
                    const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                    bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                    const HeightField* dtm, int dtm_from_class, double dtm_filter_distance,
                    std::vector<std::string>& out_shard_paths, double lambda1 = 0.0);

  bool mergeShards(const std::vector<std::string>& shard_paths, VoxelGrid& grid);

  size_t num_threads_;
  size_t ram_budget_mb_;
};


bool InProcessStrategy::execute(const std::string& cloud_name, VoxelGrid& grid,
                                const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                                bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                                const HeightField* dtm, int dtm_from_class, double dtm_filter_distance, double lambda1)
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

    // Scale queue depth to available RAM: use up to 1 GB, minimum 8 batches/thread.
    const size_t avail_ram    = ray::queryAvailableMemoryBytes();
    const size_t queue_budget = std::min(avail_ram / 20, size_t(1) * 1024 * 1024 * 1024);
    const size_t queue_depth  = std::max(queue_budget / sizeof(BeamBatch),
                                         resolved_threads * 8);

    // Scale LAS read chunk: larger chunks reduce callback overhead when the queue is deep.
    // Aim for each chunk to fill ~4 queue batches per thread.
    const size_t las_chunk = std::min(
        std::max(size_t(4000000), resolved_threads * kBeamBatchSize * 4),
        size_t(16000000));

    std::cout << "Processing point cloud using " << resolved_threads
              << " parallel threads (in-memory).\n"
              << "  Queue: " << queue_depth << " batches × " << kBeamBatchSize
              << " beams = " << queue_depth * kBeamBatchSize << " beams  ("
              << (queue_depth * sizeof(BeamBatch)) / (1024 * 1024) << " MB)\n"
              << "  LAS chunk: " << las_chunk / 1000000 << "M points" << std::endl;

    ThreadSafeQueue<BeamBatch> beam_queue(queue_depth);
    std::vector<std::thread> threads;

    const bool use_flat = grid.isFlat();
    VoxelGrid::Voxel* flat_ptr = use_flat ? grid.flat_voxels_.data() : nullptr;
    const int64_t flat_dimX  = use_flat ? grid.voxel_dims_[0] : 0;
    const int64_t flat_dimXY = use_flat ? grid.voxel_dims_[0] * grid.voxel_dims_[1] : 0;

    std::vector<VoxelProcessor::Map> worker_maps(use_flat ? 0 : resolved_threads);

    auto worker_task = [&](size_t thread_idx) {
      VoxelProcessor processor(grid.getBounds(), grid.getVoxelWidth(), weighting_method, use_occlusion,
                               apply_flat_top, peaks_ptr, calc_beam_metrics, beam_diameter,
                               tan_half_divergence, subvoxel_split, dtm, lambda1);
      if (use_flat)
        processor.setFlatTarget(flat_ptr, flat_dimX, flat_dimXY);
      BeamBatch batch;
      while (beam_queue.pop(batch)) {
        for (size_t b = 0; b < batch.count; ++b)
          processor.processBeam(batch.beams[b]);
      }
      if (!use_flat)
        worker_maps[thread_idx] = processor.takeMap();
    };

    for (size_t i = 0; i < resolved_threads; ++i)
      threads.emplace_back(worker_task, i);

    size_t num_bounded = 0;
    std::vector<uint8_t> passthrough;
    std::vector<int32_t> beam_ids_chunk;
    uint16_t pt_extra = 0;
    static bool not_raycloud_warned = false;
    double pending_gps_time = std::numeric_limits<double>::quiet_NaN();
    int32_t pending_beam_id = -1;
    std::vector<PointData> pending_returns;
    Eigen::Vector3d pending_beam_origin;

    // Current batch being filled by the producer.
    BeamBatch current_batch;

    // Commit a completed beam into the current batch; push when the batch is full.
    auto flush_beam = [&]() {
      if (pending_returns.empty()) return;
      BeamData& bd = current_batch.beams[current_batch.count];
      bd.beam_origin = pending_beam_origin;
      bd.gps_time    = pending_gps_time;
      bd.num_returns = static_cast<uint8_t>(std::min(pending_returns.size(),
                                            static_cast<size_t>(kMaxReturnsPerBeam)));
      for (uint8_t r = 0; r < bd.num_returns; ++r)
        bd.returns[r] = pending_returns[r];
      ++current_batch.count;
      pending_returns.clear();
      if (current_batch.count == kBeamBatchSize) {
        BeamBatch tmp = current_batch;  // copy before reset so workers get valid data
        current_batch.count = 0;
        beam_queue.push(std::move(tmp));
      }
    };

    ray::readLas(cloud_name,
      [&](std::vector<Eigen::Vector3d>& starts, std::vector<Eigen::Vector3d>& ends,
          std::vector<double>& times, std::vector<ray::RGBA>& colours) {
        if (starts.empty() && !not_raycloud_warned) {
          std::cerr << "Warning: input is not a ray cloud (no sx,sy,sz ray starts); skipping points." << std::endl;
          not_raycloud_warned = true;
        }
        for (size_t i = 0; i < ends.size(); ++i) {
          if (starts.empty()) {
            // Non-raycloud file (start == end for every point); nothing to traverse.
            continue;
          }
          const int32_t bid = (i < beam_ids_chunk.size()) ? beam_ids_chunk[i] : -1;
          const uint8_t alpha = (i < colours.size()) ? colours[i].alpha : 1;
          PointData pd = makePointData(starts[i], ends[i], times[i], bid, alpha, passthrough, i, stride);
          if (isGroundHit(pd.x, pd.y, pd.z, pd.classification, dtm_from_class, dtm, dtm_filter_distance))
            pd.bound = 0;
          const bool new_beam = isNewBeam(pd, pending_gps_time, pending_beam_id,
                                          pending_returns, !beam_ids_chunk.empty());
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
      }, num_bounded, 255.0, nullptr, las_chunk,
         nullptr, &passthrough, &pt_extra, nullptr, nullptr, &beam_ids_chunk);

    flush_beam();  // commit the last beam
    if (current_batch.count > 0) {
      beam_queue.push(std::move(current_batch));
    }
    beam_queue.notify_done();
    for (auto& t : threads) { t.join(); }

    if (!use_flat) {
      for (size_t i = 0; i < resolved_threads; ++i)
        grid.absorbMap(std::move(worker_maps[i]));
    }
    std::cout << "Parallel processing finished." << std::endl;

  } else {
    // --- OPTION 1: Single-Threaded Implementation ---
    std::cout << "Processing point cloud using 1 thread (in-memory)..." << std::endl;
    VoxelProcessor processor(grid.getBounds(), grid.getVoxelWidth(), weighting_method, use_occlusion,
                             apply_flat_top, peaks_ptr, calc_beam_metrics, beam_diameter,
                             tan_half_divergence, subvoxel_split, dtm, lambda1);
    if (grid.isFlat())
      processor.setFlatTarget(grid.flat_voxels_.data(),
                              grid.voxel_dims_[0],
                              grid.voxel_dims_[0] * grid.voxel_dims_[1]);

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
        beam.num_returns = static_cast<uint8_t>(std::min(pending_returns.size(),
                             static_cast<size_t>(kMaxReturnsPerBeam)));
        for (uint8_t r = 0; r < beam.num_returns; ++r) {
          beam.returns[r] = pending_returns[r];
        }
        processor.processBeam(beam);
        pending_returns.clear();
      }
    };
    ray::readLas(cloud_name,
      [&](std::vector<Eigen::Vector3d>& starts, std::vector<Eigen::Vector3d>& ends,
          std::vector<double>& times, std::vector<ray::RGBA>& colours) {
        if (starts.empty() && !not_raycloud_warned) {
          std::cerr << "Warning: input is not a ray cloud (no sx,sy,sz ray starts); skipping points." << std::endl;
          not_raycloud_warned = true;
        }
        for (size_t i = 0; i < ends.size(); ++i) {
          if (starts.empty()) {
            // Non-raycloud file (start == end for every point); nothing to traverse.
            continue;
          }
          const int32_t bid = (i < beam_ids_chunk.size()) ? beam_ids_chunk[i] : -1;
          const uint8_t alpha = (i < colours.size()) ? colours[i].alpha : 1;
          PointData pd = makePointData(starts[i], ends[i], times[i], bid, alpha, passthrough, i, stride);
          if (isGroundHit(pd.x, pd.y, pd.z, pd.classification, dtm_from_class, dtm, dtm_filter_distance))
            pd.bound = 0;
          const bool new_beam = isNewBeam(pd, pending_gps_time, pending_beam_id,
                                          pending_returns, !beam_ids_chunk.empty());
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

    flush_beam();
    // Flat path: writes already landed in flat_voxels_ — nothing to move.
    // Sparse path: move the processor's map into the grid.
    if (!grid.isFlat())
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

  void writeAccumulated(const VoxelCoord& coord, const VoxelGrid::Voxel& v) {
    if (!target_grid_.use_sparse_fallback_) {
      target_grid_.voxelAt(target_grid_.flatIndex(coord.x, coord.y, coord.z)) = v;
    } else {
      target_grid_.getSparseVoxels_internal()[coord] = v;
    }
  }

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

    // Validate shard headers and prime the queue.
    for (size_t i = 0; i < shard_streams.size(); ++i) {
      if (!readShardHeader(shard_streams[i])) {
        std::cerr << "Error: Shard file has invalid or incompatible header: " << shard_paths_[i] << std::endl;
        return false;
      }
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
        accumulator += entry.voxel;
      } else {
        writeAccumulated(current_coord, accumulator);
        current_coord = entry.coord;
        accumulator = entry.voxel;
      }

      VoxelCoord next_coord;
      VoxelGrid::Voxel next_voxel;
      if (readVoxelData(shard_streams[entry.shard_index], next_coord, next_voxel)) {
        pq.push({next_coord, next_voxel, entry.shard_index});
      }
    }

    writeAccumulated(current_coord, accumulator);

    if (target_grid_.isFlat()) {
      std::cout << "Merge complete (flat grid)." << std::endl;
    } else {
      std::cout << "Merge complete. Final grid has " << target_grid_.getSparseVoxels().size() << " voxels." << std::endl;
    }
    return true;
  }
};


bool OutOfCoreStrategy::execute(const std::string& cloud_name, VoxelGrid& grid,
                                const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                                bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                                const HeightField* dtm, int dtm_from_class, double dtm_filter_distance, double lambda1) {
    std::vector<std::string> shard_paths;
    std::cout << "Starting out-of-core processing..." << std::endl;

    if (!createShards(cloud_name, grid, weighting_method, use_occlusion, apply_flat_top,
                      calc_beam_metrics, beam_diameter, beam_divergence, subvoxel_split, dtm,
                      dtm_from_class, dtm_filter_distance, shard_paths, lambda1)) {
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
                                     const HeightField* dtm, int dtm_from_class, double dtm_filter_distance,
                                     std::vector<std::string>& out_shard_paths, double lambda1) {
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
                               tan_half_divergence, subvoxel_split, dtm, lambda1);
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
        beam.num_returns = static_cast<uint8_t>(std::min(pending_returns.size(),
                             static_cast<size_t>(kMaxReturnsPerBeam)));
        for (uint8_t r = 0; r < beam.num_returns; ++r) {
          beam.returns[r] = pending_returns[r];
        }
        beam_queue.push(std::move(beam));
        pending_returns.clear();
      }
    };
    ray::readLas(cloud_name,
      [&](std::vector<Eigen::Vector3d>& starts, std::vector<Eigen::Vector3d>& ends,
          std::vector<double>& times, std::vector<ray::RGBA>& colours) {
        if (starts.empty() && !not_raycloud_warned) {
          std::cerr << "Warning: input is not a ray cloud (no sx,sy,sz ray starts); skipping points." << std::endl;
          not_raycloud_warned = true;
        }
        for (size_t i = 0; i < ends.size(); ++i) {
          if (starts.empty()) {
            // Non-raycloud file (start == end for every point); nothing to traverse.
            continue;
          }
          const int32_t bid = (i < beam_ids_chunk.size()) ? beam_ids_chunk[i] : -1;
          const uint8_t alpha = (i < colours.size()) ? colours[i].alpha : 1;
          PointData pd = makePointData(starts[i], ends[i], times[i], bid, alpha, passthrough, i, stride);
          if (isGroundHit(pd.x, pd.y, pd.z, pd.classification, dtm_from_class, dtm, dtm_filter_distance))
            pd.bound = 0;
          const bool new_beam = isNewBeam(pd, pending_gps_time, pending_beam_id,
                                          pending_returns, !beam_ids_chunk.empty());
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

    // Combined pre-scan result: a single endpoint-only readLas pass can serve both the
    // auto-bounds detection and the auto-reserve-size count, since both are lightweight
    // endpoint reads over the same file.
    struct PreScanResult {
      Eigen::Vector3d bounds_min = Eigen::Vector3d::Constant(std::numeric_limits<double>::max());
      Eigen::Vector3d bounds_max = Eigen::Vector3d::Constant(std::numeric_limits<double>::lowest());
      size_t point_count = 0;
    };

    // Whether each half of the pre-scan is needed. Auto-bounds is needed only when the user
    // supplied no explicit grid bounds; the reserve count is needed only for in-memory modes
    // without an explicit --reserve_size.
    const bool need_auto_bounds = isZeroVector(params.grid_bounds_min) && isZeroVector(params.grid_bounds_max);
    const bool need_reserve_count = (params.reserve_size == 0) && !params.use_ooc;

    auto preScanCloud = [&](const std::string& cloud_name, PreScanResult& result) -> bool {
      // Determine whether this file declares a "bound" extra attribute. When present, unbound
      // (miss) endpoints are excluded from the bounding box; old files without it include all
      // endpoints exactly as before.
      uint16_t pre_orig_extra = 0;
      std::vector<uint8_t> pre_extra_vlr;
      bool file_has_bound = false;
      readLasExtraBytesVlr(cloud_name, pre_orig_extra, pre_extra_vlr, &file_has_bound);
      size_t num_bounded = 0;
      return ray::readLas(cloud_name,
          [&](std::vector<Eigen::Vector3d>& /*starts*/, std::vector<Eigen::Vector3d>& ends,
              std::vector<double>& /*times*/, std::vector<ray::RGBA>& colours) {
            for (size_t i = 0; i < ends.size(); ++i) {
              // bound == 0 (alpha == 0) marks an unbound ray with a floating far end; exclude it
              // from the bounds. point_count still counts every point (over-reserve is fine).
              const bool is_unbound = file_has_bound && i < colours.size() && colours[i].alpha == 0;
              if (!is_unbound) {
                result.bounds_min = result.bounds_min.cwiseMin(ends[i]);
                result.bounds_max = result.bounds_max.cwiseMax(ends[i]);
              }
            }
            result.point_count += ends.size();
          }, num_bounded, 255.0, nullptr);
    };

    // Run the pre-scan once if either half is needed; reuse its results below.
    PreScanResult pre_scan;
    if (need_auto_bounds || need_reserve_count) {
      if (need_auto_bounds) {
        std::cout << "Auto-detecting grid bounds from file..." << std::endl;
      }
      if (!preScanCloud(params.cloud_name, pre_scan)) {
        std::cerr << "Error: Could not read LAS/LAZ file to pre-scan the cloud." << std::endl;
        return false;
      }
    }

    double beam_diameter = 0.0;
    double beam_divergence = 0.0;

    // Stage 3 effective free path coefficient λ₁ = 0.25·avg_leaf_area / voxel_size³.
    // Computed once; threaded to every VoxelProcessor. Zero disables the correction (eff(z)=z).
    const double lambda1 = (params.average_leaf_area > 0.0 && params.voxel_size > 0.0)
        ? 0.25 * params.average_leaf_area / (params.voxel_size * params.voxel_size * params.voxel_size)
        : 0.0;

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

    if (need_auto_bounds) {
        if (pre_scan.bounds_min.x() > pre_scan.bounds_max.x()) {
            std::cerr << "Error: No points found to determine bounds." << std::endl;
            return false;
        }
        user_bounds.min_bound_ = pre_scan.bounds_min;
        user_bounds.max_bound_ = pre_scan.bounds_max;
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

    // Handle auto-detection of reservation size here. The point count was already gathered
    // by the combined pre-scan pass above (need_reserve_count mirrors that pass's trigger).
    size_t final_reserve_size = params.reserve_size;
    if (need_reserve_count) {
        final_reserve_size += pre_scan.point_count;
    }
    // For OOC, we never reserve in the final grid, as it's populated at the end.
    if (params.use_ooc) { final_reserve_size = 0; }

    std::unique_ptr<VoxelGrid> grid_ptr;
    try {
        const size_t ram_budget_bytes = static_cast<size_t>(params.ram_budget_mb) * 1024ULL * 1024;
        grid_ptr = std::make_unique<VoxelGrid>(processing_bounds, params.voxel_size,
                                               ram_budget_bytes, final_reserve_size);
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
                                                 params.calc_beam_metrics, beam_diameter, beam_divergence, params.subvoxel_split, dtm_ptr.get(),
                                                 params.dtm_from_class, params.dtm_filter_distance, lambda1);

    if (!processing_success) {
      return false; // Strategy failed, exit early.
    }

    // === POST-PROCESSING AND OUTPUT (Same for all strategies) ===
    if (params.neighbour_prior_min_rays > 0) {
        std::cout << "Applying neighbour priors (min rays = " << params.neighbour_prior_min_rays << ")..." << std::endl;
        applyNeighbourPriors(grid, params.neighbour_prior_min_rays);
    }

    // Build the classification table in a separate O(N_points) pass (no ray walking).
    // When inclination distributions are active, the class and IAD tables share an
    // identical I/O scan, so build them together in a single readLas pass.
    ClassTable class_table;
    IadTable iad_table;
    if (params.calc_inclination_dist) {
      {
        static bool empty_classes_warned = false;
        if (!empty_classes_warned && params.leaf_classes_str.empty() && params.wood_classes_str.empty()) {
          std::cerr << "Warning: --inclination_dist active but neither --leaf_classes nor --wood_classes "
                       "is set; LAD and WAD will be zero (PAD is unaffected)." << std::endl;
          empty_classes_warned = true;
        }
      }
      std::cout << "Building classification table and inclination angle distributions..." << std::endl;
      buildClassAndIadTable(params.cloud_name, grid, params, dtm_ptr.get(), class_table, iad_table);
    } else {
      {
        static bool field_no_iad_warned = false;
        if (!field_no_iad_warned &&
            (params.leaf_classes_str.find(':') != std::string::npos ||
             params.wood_classes_str.find(':') != std::string::npos)) {
          std::cerr << "Warning: a class field prefix ('field:codes') was given but --inclination_dist "
                       "is off; the field-aware leaf/wood path is inactive, so leaf/wood counts use the "
                       "standard Classification byte and may be wrong." << std::endl;
          field_no_iad_warned = true;
        }
      }
      std::cout << "Building classification table..." << std::endl;
      class_table = buildClassTable(params.cloud_name, grid, params.dtm_from_class, dtm_ptr.get(), params.dtm_filter_distance);
    }

    std::cout << "Calculating output metrics..." << std::endl;
    MetricResultsMap metrics = calculateOutputMetrics(grid, params, dtm_ptr.get(), class_table, iad_table);

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
        primary_success = writeNetcdfFile(primary_name_stub, grid, metrics, padding, user_bounds, primary_params, iad_table);
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
            filled_success = writeNetcdfFile(filled_name_stub, grid, metrics, padding, user_bounds, filled_params, iad_table, true);
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
