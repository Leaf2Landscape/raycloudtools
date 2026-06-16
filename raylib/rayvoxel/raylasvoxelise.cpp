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
#include "raylib/raycloudreader.h"
#include "raylib/raysysinfo.h"
#include "raylib/rayvoxel/rayvox.h"
#include "raylib/rayvoxel/raylasvoxelise.h"
#include "raylib/rayvoxel/raylasvoxelwriter.h"
#include "raylib/rayvoxel/raylaswoodvolume.h"
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
#include <cassert>
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
                     size_t ram_budget_bytes, size_t sparse_reservation, bool allocate_peaks)
    : bounds_(grid_bounds), voxel_width_(vox_width)
{
  Eigen::Vector3d extent = bounds_.max_bound_ - bounds_.min_bound_;
  Eigen::Vector3d temp = (extent / voxel_width_).array().ceil();
  voxel_dims_ = temp.cast<int64_t>();

  const int64_t max_reasonable_dim = 1000000000L;
  if (voxel_dims_.maxCoeff() > max_reasonable_dim) {
    throw std::runtime_error("VoxelGrid Error: Resulting dimensions on one or more axes are unreasonably large.");
  }

  // Only allocate the (x,y) peaks store when flat-top compensation is requested. When the flat
  // array would exceed an eighth of the RAM budget, hold peaks sparsely instead of as a flat vector.
  if (allocate_peaks) {
    const size_t peaks_bytes = static_cast<size_t>(voxel_dims_[0] * voxel_dims_[1]) * sizeof(double);
    if (peaks_bytes <= ram_budget_bytes / 8) {
      peaks_.resize(static_cast<size_t>(voxel_dims_[0] * voxel_dims_[1]),
                    std::numeric_limits<double>::lowest());
      use_sparse_peaks_ = false;
    } else {
      use_sparse_peaks_ = true;
    }
  }

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

void VoxelGrid::setPeak(int64_t xy_idx, double value)
{
  if (!use_sparse_peaks_) {
    peaks_[xy_idx] = value;
  } else {
    sparse_peaks_[xy_idx] = value;
  }
}

double VoxelGrid::getPeak(int64_t xy_idx) const
{
  if (!use_sparse_peaks_) {
    return peaks_[xy_idx];
  }
  auto it = sparse_peaks_.find(xy_idx);
  return it != sparse_peaks_.end() ? it->second : std::numeric_limits<double>::lowest();
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

// Groups a cloud's returns into pulses (beams) independent of point ordering, by hashing on the
// pulse key (beam_id when present, else gps_time). A pulse is emitted as soon as all its returns
// have arrived (collected count reaches number_of_returns); single-return pulses emit immediately;
// stragglers (incomplete pulses, e.g. filtered/missing returns) are emitted at flush().
//
// Replaces the previous file-adjacency grouping, which split multi-return pulses into per-return
// beams on raycloud .laz files whose points are NOT pulse-ordered (returns sharing a gps_time are
// scattered across the file). That split each return into its own beam, re-walking the shared
// origin→hit path — inflating per-voxel beam-crossing counts (nbSampling) and the PPL miss term,
// and giving split returns full weight instead of the correct decreasing per-segment weight.
class PulseGrouper {
public:
    // Add one return. If its pulse is now complete, invokes emit(BeamData&&).
    template <class Emit>
    void add(const PointData& pd, Emit&& emit) {
        const uint8_t nor = (pd.number_of_returns < 1) ? 1 : pd.number_of_returns;
        if (nor <= 1) {  // single-return pulse → emit immediately, never buffered
            BeamData b;
            b.beam_origin = pd.beam_origin;
            b.gps_time    = pd.gps_time;
            b.num_returns = 1;
            b.returns[0]  = pd;
            emit(std::move(b));
            return;
        }
        // beam_id (when present, >= 0) is the authoritative key; otherwise gps_time (bit-exact per
        // pulse). A cloud uses one regime consistently, so the two never collide.
        const double key = (pd.beam_id >= 0) ? static_cast<double>(pd.beam_id) : pd.gps_time;
        Partial& p = pending_[key];
        if (p.returns.empty()) { p.origin = pd.beam_origin; p.gps_time = pd.gps_time; p.target = nor; }
        p.returns.push_back(pd);
        if (p.returns.size() >= p.target) { emitPulse(p, emit); pending_.erase(key); }
    }

    // Emit any pulses still incomplete at end-of-stream (missing/filtered returns).
    template <class Emit>
    void flush(Emit&& emit) {
        for (auto& kv : pending_)
            if (!kv.second.returns.empty()) emitPulse(kv.second, emit);
        pending_.clear();
    }

private:
    struct Partial {
        Eigen::Vector3d origin;
        double gps_time = 0.0;
        uint8_t target = 0;
        std::vector<PointData> returns;
    };
    template <class Emit>
    static void emitPulse(const Partial& p, Emit& emit) {
        BeamData b;
        b.beam_origin = p.origin;
        b.gps_time    = p.gps_time;
        b.num_returns = static_cast<uint8_t>(std::min(p.returns.size(), static_cast<size_t>(kMaxReturnsPerBeam)));
        for (uint8_t r = 0; r < b.num_returns; ++r) b.returns[r] = p.returns[r];
        emit(std::move(b));
    }
    std::unordered_map<double, Partial> pending_;
};

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
ClassFieldSource resolveClassField(const std::string& field_name, const ray::LasHeader& hdr)
{
  if (field_name.empty()) return ClassFieldSource{};
  std::string lowered = field_name;
  for (char& c : lowered) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  if (lowered == "classification") return ClassFieldSource{};

  const ray::LasExtraField *f = hdr.field(field_name);
  if (f && !f->is_own)
    return ClassFieldSource{ static_cast<uint16_t>(kPassthroughStdBytes + f->sensor_offset), f->dtype };

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

// Per-tile accumulator for the tiled parallel KNN/IAD pass. Each worker thread owns one
// TileResult; histograms are merged serially after the thread pool joins.
struct TileResult {
  // Histograms and hit counts are keyed by tree_id (joined across stems), not voxel index.
  std::unordered_map<int32_t, std::vector<double>> all_hist, leaf_hist, wood_hist, beam_hist;
  std::unordered_map<int32_t, float> leaf_hit_count, wood_hit_count;
  std::unordered_map<int64_t, TriangleHistograms> triangle_histograms;
  // Per-voxel tally of tree_id → point-hit count, reduced to predominant_tree during the merge.
  std::unordered_map<int64_t, std::unordered_map<int32_t, int32_t>> tree_hits;
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
                                  PerTreeIadMap& per_tree_iad_out,
                                  PredominantTreeTable& predominant_tree_out)
{
  // Compact per-tile point record. Replaces the parallel global positions/leaf_vals/wood_vals/
  // flat_indices/beam_angles arrays: routing points into per-tile buckets during the readLas pass
  // bounds peak memory at O(max single-tile point count) rather than O(N_total).
  struct TilePoint {
    Eigen::Vector3d pos;
    int leaf_val;
    int wood_val;
    int64_t flat_idx;
    double beam_angle;  // zenith angle [0, pi/2] of the point's inbound ray
    int32_t tree_id;    // raycloud tree_id of the point; -1 when absent (excluded from per-tree IAD)
  };

  PerTreeIadMap per_tree_iad;
  PredominantTreeTable predominant_tree;

  ray::CloudReader reader;
  reader.begin(cloud_name);
  const ray::LasHeader &hdr = reader.header();
  const bool has_tree_id = hdr.has("tree_id");
  const std::vector<uint8_t> extra_bytes_vlr = hdr.sensorExtraVlr();
  const uint16_t stride = static_cast<uint16_t>(kPassthroughStdBytes + hdr.sensorExtraSize());

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

  const ClassFieldSource leaf_src = resolveClassField(leaf_field, hdr);
  const ClassFieldSource wood_src = resolveClassField(wood_field, hdr);

  // Parse leaf/wood class sets before readLas so they can be used in the readLas callback
  // for per-voxel leaf/wood counting. The sets are also used later in the tile KNN loop.
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

  // Tile geometry is fixed by the grid bounds, so it can be computed before the readLas pass —
  // this lets each point be routed straight into its core-tile bucket as it is read.
  // buf_m must be >= tile so 3x3 neighbour scan covers full buffer.
  const double buf_m = 1.0;
  const double tile_sz = std::max(buf_m, params.iad_tile_size);
  const double minx = bounds.min_bound_.x(), miny = bounds.min_bound_.y();
  const double maxx = bounds.max_bound_.x(), maxy = bounds.max_bound_.y();
  const int n_tx = std::max(1, (int)std::ceil((maxx - minx) / tile_sz));
  const int n_ty = std::max(1, (int)std::ceil((maxy - miny) / tile_sz));
  const int n_tiles = n_tx * n_ty;

  // Route each bounded hit point into its core-tile bucket. For n_tiles == 1 this is a single
  // bucket processed by the fast path; otherwise tiles are processed sequentially/in parallel,
  // each discarding its points before the next is built.
  std::vector<std::vector<TilePoint>> tile_points(n_tiles);
  // Owner tile per voxel for the bailey path: the tile of the first point routed into each
  // flat_idx. Populated inline during routing to avoid a second pass over the points.
  std::unordered_map<int64_t, int> flat_owner;

  // tree_id is a raycloud attribute decoded by readLas into tree_ids_all at the global point
  // index (parallel to the passthrough buffer). When the input has no tree_id field, per-tree
  // IAD is disabled: predominant_tree stays -1 everywhere and no {stub}_iad.csv is written, but
  // the class table is still built below (it is independent of tree identity).
  if (!has_tree_id) {
    std::cerr << "Info: no tree_id field in input; per-tree IAD disabled, predominant_tree=-1." << std::endl;
  }

  size_t num_bounded = 0;
  std::vector<uint8_t> passthrough;
  std::vector<int32_t> tree_ids_all;
  size_t global_chunk_start = 0;
  size_t total_points = 0;

  reader.read(
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

        const uint8_t classification = passthrough[base + 2];
        if (isGroundHit(ends[i].x(), ends[i].y(), ends[i].z(), classification,
                        params.dtm_from_class, dtm, params.dtm_filter_distance)) continue;
        const int64_t flat_idx = grid.flatIndex(ix, iy, iz);

        // IadTable collection (same as buildIadTable): leaf/wood field values + flat index.
        const int lv = readClassValue(&passthrough[base], leaf_src);
        const int wv = readClassValue(&passthrough[base], wood_src);
        const Eigen::Vector3d dir = ends[i] - starts[i];
        const double len2 = dir.squaredNorm();
        const double bz = (len2 > 1e-12) ? std::acos(std::min(1.0, std::abs(dir.z() / std::sqrt(len2)))) : 0.0;

        // tree_id grouping key (joins across stems). -1 when the file has no tree_id field, or for
        // any unlabeled point; such points are excluded from per-tree histograms.
        const int32_t tid = has_tree_id ? tree_ids_all[global_chunk_start + i] : -1;

        // Route into the point's core tile.
        const int tx = std::clamp((int)((ends[i].x() - minx) / tile_sz), 0, n_tx - 1);
        const int ty = std::clamp((int)((ends[i].y() - miny) / tile_sz), 0, n_ty - 1);
        const int t  = ty * n_tx + tx;
        tile_points[t].push_back(TilePoint{ ends[i], lv, wv, flat_idx, bz, tid });
        ++total_points;

        // Owner tile = tile of the first point that maps to this flat_idx (lowest global index).
        if (any_bailey && flat_owner.find(flat_idx) == flat_owner.end()) flat_owner[flat_idx] = t;
      }
      global_chunk_start += ends.size();
    }, num_bounded, 255.0, nullptr, 1000000, &tree_ids_all, &passthrough);

  if (total_points < 2) {
    per_tree_iad_out = std::move(per_tree_iad);
    predominant_tree_out = std::move(predominant_tree);
    return;
  }

  size_t resolved_threads = (params.num_threads == 0)
      ? std::thread::hardware_concurrency()
      : static_cast<size_t>(params.num_threads);
  if (resolved_threads == 0) resolved_threads = 1;

  // Histograms / hit counts keyed by tree_id; tree_hits keyed by voxel for predominant_tree.
  std::unordered_map<int32_t, std::vector<double>> all_hist, leaf_hist, wood_hist, beam_hist;
  std::unordered_map<int32_t, float> leaf_hit_count, wood_hit_count;
  std::unordered_map<int64_t, TriangleHistograms> triangle_histograms;
  std::unordered_map<int64_t, std::unordered_map<int32_t, int32_t>> tree_hits;

  if (n_tiles == 1) {
    // Small cloud: single bucket, global KD-tree, no tiling overhead.
    const std::vector<TilePoint>& tp = tile_points[0];
    const int K = std::min(params.knn_normal, (int)tp.size() - 1);
    Eigen::MatrixXd points_p(3, tp.size());
    for (size_t i = 0; i < tp.size(); ++i) points_p.col(i) = tp[i].pos;
    std::unique_ptr<Nabo::NNSearchD> nns(Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3));
    Eigen::MatrixXi indices(K, (int)tp.size());
    Eigen::MatrixXd dists2(K, (int)tp.size());
    nns->knn(points_p, indices, dists2, K, kNearestNeighbourEpsilon, 0);
    nns.reset(nullptr);

    if (any_bailey) {
      std::vector<Eigen::Vector3d> positions(tp.size());
      std::vector<int64_t> flat_indices(tp.size());
      std::vector<int> class_labels_int(tp.size(), 0);
      for (size_t i = 0; i < tp.size(); ++i) {
        positions[i]    = tp[i].pos;
        flat_indices[i] = tp[i].flat_idx;
        if (leaf_set.count(tp[i].leaf_val)) class_labels_int[i] = 1;
        else if (wood_set.count(tp[i].wood_val)) class_labels_int[i] = -1;
      }
      triangle_histograms = buildTriangleInclinationHistograms(
          positions, indices, flat_indices, class_labels_int, params.n_iad_bins, params.triangle_lmax);
    }

    for (size_t i = 0; i < tp.size(); ++i) {
      Eigen::Vector3d centroid(0, 0, 0);
      int num_neighbours = 0;
      for (int j = 0; j < K && indices(j, i) != Nabo::NNSearchD::InvalidIndex; ++j) {
        centroid += tp[indices(j, i)].pos;
        ++num_neighbours;
      }
      if (num_neighbours < 3) continue;
      centroid /= static_cast<double>(num_neighbours);
      Eigen::Matrix3d scatter = Eigen::Matrix3d::Zero();
      for (int j = 0; j < K && indices(j, i) != Nabo::NNSearchD::InvalidIndex; ++j) {
        Eigen::Vector3d offset = tp[indices(j, i)].pos - centroid;
        scatter += offset * offset.transpose();
      }
      scatter /= static_cast<double>(num_neighbours);

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter);
      const Eigen::Vector3d normal = eigen_solver.eigenvectors().col(0);
      const double theta = std::acos(std::min(1.0, std::abs(normal.z())));

      int bin = static_cast<int>(theta / (kPi / 2.0) * params.n_iad_bins);
      bin = std::clamp(bin, 0, params.n_iad_bins - 1);

      // predominant_tree tally: count this point's hit toward its voxel's tree histogram.
      const int64_t flat_idx = tp[i].flat_idx;
      const int32_t tid = tp[i].tree_id;
      if (tid >= 0) tree_hits[flat_idx][tid] += 1;
      // Points without a tree_id contribute to no tree's IAD; skip the per-tree histograms.
      if (tid < 0) continue;

      auto& ah = all_hist[tid];
      if (ah.empty()) ah.assign(params.n_iad_bins, 0.0);
      ah[bin] += 1.0;

      {
        int bbin = static_cast<int>(tp[i].beam_angle / (kPi / 2.0) * params.n_iad_bins);
        bbin = std::clamp(bbin, 0, params.n_iad_bins - 1);
        auto& bh = beam_hist[tid];
        if (bh.empty()) bh.assign(params.n_iad_bins, 0.0);
        bh[bbin] += 1.0;
      }

      if (leaf_set.count(tp[i].leaf_val)) {
        leaf_hit_count[tid] += 1.0f;
        auto& lh = leaf_hist[tid];
        if (lh.empty()) lh.assign(params.n_iad_bins, 0.0);
        lh[bin] += 1.0;
      }
      if (wood_set.count(tp[i].wood_val)) {
        wood_hit_count[tid] += 1.0f;
        auto& wh = wood_hist[tid];
        if (wh.empty()) wh.assign(params.n_iad_bins, 0.0);
        wh[bin] += 1.0;
      }
    }
  } else {
  // Large cloud: tiled parallel KNN. Points are already bucketed by core tile in tile_points
  // (filled during the readLas pass), and flat_owner was populated inline during routing.

  std::vector<TileResult> results(resolved_threads);
  std::atomic<int> next_tile(0);

  auto worker = [&](size_t w) {
    TileResult& tr = results[w];
    int t;
    while ((t = next_tile.fetch_add(1)) < n_tiles) {
      if (tile_points[t].empty()) continue;

      const int tx = t % n_tx;
      const int ty = t / n_tx;
      const double cx0 = minx + tx * tile_sz;
      const double cx1 = std::min(maxx, cx0 + tile_sz);
      const double cy0 = miny + ty * tile_sz;
      const double cy1 = std::min(maxy, cy0 + tile_sz);

      // Collect buffered points from the 3x3 tile neighbourhood of the already-bucketed data.
      // buf_core marks the entries belonging to this tile's own bucket (core points); buffer-only
      // points from adjacent tiles correct the KNN at tile edges but are not accumulated.
      std::vector<const TilePoint*> buf;
      std::vector<char> buf_core;
      for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx) {
          int nx2 = tx + dx, ny2 = ty + dy;
          if (nx2 < 0 || nx2 >= n_tx || ny2 < 0 || ny2 >= n_ty) continue;
          const int nt = ny2 * n_tx + nx2;
          for (const TilePoint& q : tile_points[nt]) {
            const Eigen::Vector3d& p = q.pos;
            if (p.x() >= cx0 - buf_m && p.x() <= cx1 + buf_m &&
                p.y() >= cy0 - buf_m && p.y() <= cy1 + buf_m) {
              buf.push_back(&q);
              buf_core.push_back(nt == t ? 1 : 0);
            }
          }
        }

      const size_t Nb = buf.size();
      if (Nb < 3) continue;

      const int K = std::min(params.knn_normal, (int)Nb - 1);

      // Build tile-local KD-tree over buffered points.
      Eigen::MatrixXd pts(3, Nb);
      for (size_t c = 0; c < Nb; ++c) pts.col(c) = buf[c]->pos;
      std::unique_ptr<Nabo::NNSearchD> nns(Nabo::NNSearchD::createKDTreeLinearHeap(pts, 3));
      Eigen::MatrixXi idx(K, Nb);
      Eigen::MatrixXd d2(K, Nb);
      nns->knn(pts, idx, d2, K, kNearestNeighbourEpsilon, 0);
      nns.reset(nullptr);

      // PCA normal estimation — accumulate only for core points.
      for (size_t c = 0; c < Nb; ++c) {
        if (!buf_core[c]) continue;  // skip buffer-only points

        Eigen::Vector3d centroid(0, 0, 0);
        int num_neighbours = 0;
        for (int j = 0; j < K && idx(j, c) != Nabo::NNSearchD::InvalidIndex; ++j) {
          centroid += buf[idx(j, c)]->pos;
          ++num_neighbours;
        }
        if (num_neighbours < 3) continue;
        centroid /= static_cast<double>(num_neighbours);
        Eigen::Matrix3d scatter = Eigen::Matrix3d::Zero();
        for (int j = 0; j < K && idx(j, c) != Nabo::NNSearchD::InvalidIndex; ++j) {
          Eigen::Vector3d offset = buf[idx(j, c)]->pos - centroid;
          scatter += offset * offset.transpose();
        }
        scatter /= static_cast<double>(num_neighbours);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter);
        const Eigen::Vector3d normal = eigen_solver.eigenvectors().col(0);
        const double theta = std::acos(std::min(1.0, std::abs(normal.z())));

        int bin = static_cast<int>(theta / (kPi / 2.0) * params.n_iad_bins);
        bin = std::clamp(bin, 0, params.n_iad_bins - 1);

        // predominant_tree tally: count this core point's hit toward its voxel's tree histogram.
        const int64_t flat_idx = buf[c]->flat_idx;
        const int32_t tid = buf[c]->tree_id;
        if (tid >= 0) tr.tree_hits[flat_idx][tid] += 1;
        // Points without a tree_id contribute to no tree's IAD; skip the per-tree histograms.
        if (tid < 0) continue;

        auto& ah = tr.all_hist[tid];
        if (ah.empty()) ah.assign(params.n_iad_bins, 0.0);
        ah[bin] += 1.0;

        {
          int bbin = static_cast<int>(buf[c]->beam_angle / (kPi / 2.0) * params.n_iad_bins);
          bbin = std::clamp(bbin, 0, params.n_iad_bins - 1);
          auto& bh = tr.beam_hist[tid];
          if (bh.empty()) bh.assign(params.n_iad_bins, 0.0);
          bh[bbin] += 1.0;
        }

        if (leaf_set.count(buf[c]->leaf_val)) {
          tr.leaf_hit_count[tid] += 1.0f;
          auto& lh = tr.leaf_hist[tid];
          if (lh.empty()) lh.assign(params.n_iad_bins, 0.0);
          lh[bin] += 1.0;
        }
        if (wood_set.count(buf[c]->wood_val)) {
          tr.wood_hit_count[tid] += 1.0f;
          auto& wh = tr.wood_hist[tid];
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
          local_pos[c]  = buf[c]->pos;
          local_flat[c] = buf[c]->flat_idx;
          if (leaf_set.count(buf[c]->leaf_val)) labels[c] = 1;
          else if (wood_set.count(buf[c]->wood_val)) labels[c] = -1;
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
  auto merge_hist = [&](std::unordered_map<int32_t, std::vector<double>>& dst,
                        std::unordered_map<int32_t, std::vector<double>>& src) {
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
    for (auto& kv : tr.tree_hits) {
      auto& dst = tree_hits[kv.first];
      for (auto& tc : kv.second) dst[tc.first] += tc.second;
    }
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

  // Reduce each voxel's tree-hit tally to its predominant tree (argmax point hits). Iterate the
  // inner counts in ascending tree_id and keep the first strict maximum, so ties resolve to the
  // lowest tree_id deterministically regardless of thread/merge order.
  for (const auto& pair : tree_hits) {
    std::vector<int32_t> tids;
    tids.reserve(pair.second.size());
    for (const auto& tc : pair.second) tids.push_back(tc.first);
    std::sort(tids.begin(), tids.end());
    int32_t best_tid = -1, best_count = 0;
    for (int32_t k : tids) {
      const int32_t c = pair.second.at(k);
      if (c > best_count) { best_count = c; best_tid = k; }
    }
    predominant_tree[pair.first] = best_tid;
  }

  // Build one IadData per tree_id (joined across stems) from the tree-keyed histograms.
  for (auto& pair : all_hist) {
    const int32_t tid = pair.first;
    IadData iad;
    iad.bin_centres = bin_centres;
    iad.piad = pair.second;  // all points -> plant
    auto lit = leaf_hist.find(tid);
    iad.liad = (lit != leaf_hist.end()) ? lit->second : std::vector<double>(params.n_iad_bins, 0.0);
    auto wit = wood_hist.find(tid);
    iad.wiad = (wit != wood_hist.end()) ? wit->second : std::vector<double>(params.n_iad_bins, 0.0);

    normalize(iad.liad);
    normalize(iad.wiad);
    normalize(iad.piad);
    iad.liad_dewit = classifyDeWit(bin_centres, iad.liad);
    iad.wiad_dewit = classifyDeWit(bin_centres, iad.wiad);
    iad.piad_dewit = classifyDeWit(bin_centres, iad.piad);

    // Angle-integrated G: weight G(theta_beam, leaf_angles) over the tree's empirical beam-direction
    // distribution rather than evaluating at a single mean angle.
    auto bhit = beam_hist.find(tid);
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
      // No beam-angle samples for this tree (rare): fall back to spherical G = 0.5.
      iad.plant_g = 0.5;
      iad.leaf_g  = 0.5;
      iad.wood_g  = 0.5;
    }

    {
      auto lhit = leaf_hit_count.find(tid);
      if (lhit != leaf_hit_count.end()) iad.leaf_hits = lhit->second;
      auto whit = wood_hit_count.find(tid);
      if (whit != wood_hit_count.end()) iad.wood_hits = whit->second;
    }

    per_tree_iad[tid] = std::move(iad);
  }

  // Bailey path (Option A): triangle-facet histograms are voxel-keyed. Fold each voxel's facets
  // into its predominant tree, area-weighted, to produce per-tree bailey histograms and mean G.
  if (any_bailey) {
    std::unordered_map<int32_t, std::vector<double>> tree_tiad_leaf, tree_tiad_wood;
    std::unordered_map<int32_t, double> g_num_leaf, g_den_leaf, g_num_wood, g_den_wood;
    for (const auto& kv : triangle_histograms) {
      auto pit = predominant_tree.find(kv.first);
      if (pit == predominant_tree.end() || pit->second < 0) continue;
      const int32_t tid = pit->second;
      const TriangleHistograms& th = kv.second;
      if (!th.tiad_leaf.empty()) {
        auto& acc = tree_tiad_leaf[tid];
        if (acc.empty()) acc.assign(params.n_iad_bins, 0.0);
        double w = 0.0;
        for (int b = 0; b < params.n_iad_bins; ++b) { acc[b] += th.tiad_leaf[b]; w += th.tiad_leaf[b]; }
        g_num_leaf[tid] += th.bailey_g_leaf * w;
        g_den_leaf[tid] += w;
      }
      if (!th.tiad_wood.empty()) {
        auto& acc = tree_tiad_wood[tid];
        if (acc.empty()) acc.assign(params.n_iad_bins, 0.0);
        double w = 0.0;
        for (int b = 0; b < params.n_iad_bins; ++b) { acc[b] += th.tiad_wood[b]; w += th.tiad_wood[b]; }
        g_num_wood[tid] += th.bailey_g_wood * w;
        g_den_wood[tid] += w;
      }
    }
    for (auto& kv : per_tree_iad) {
      const int32_t tid = kv.first;
      IadData& iad = kv.second;
      std::vector<double> piad_b(params.n_iad_bins, 0.0);
      auto lit = tree_tiad_leaf.find(tid);
      if (lit != tree_tiad_leaf.end()) {
        iad.liad_bailey = lit->second;
        for (int b = 0; b < params.n_iad_bins; ++b) piad_b[b] += lit->second[b];
        normalize(iad.liad_bailey);
      }
      auto wit = tree_tiad_wood.find(tid);
      if (wit != tree_tiad_wood.end()) {
        iad.wiad_bailey = wit->second;
        for (int b = 0; b < params.n_iad_bins; ++b) piad_b[b] += wit->second[b];
        normalize(iad.wiad_bailey);
      }
      normalize(piad_b);
      iad.piad_bailey = std::move(piad_b);
      auto gdl = g_den_leaf.find(tid);
      if (gdl != g_den_leaf.end() && gdl->second > 0.0) iad.bailey_g_leaf = g_num_leaf[tid] / gdl->second;
      auto gdw = g_den_wood.find(tid);
      if (gdw != g_den_wood.end() && gdw->second > 0.0) iad.bailey_g_wood = g_num_wood[tid] / gdw->second;
    }
  }

  per_tree_iad_out = std::move(per_tree_iad);
  predominant_tree_out = std::move(predominant_tree);
}

// ==================================================================================
// Processing Implementations (file-local free functions)
// ==================================================================================

namespace {

// In-memory traversal (single-threaded or parallel producer-consumer).
bool runInProcess(const std::string& cloud_name, VoxelGrid& grid, size_t num_threads,
                  const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                  bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                  const HeightField* dtm, int dtm_from_class, double dtm_filter_distance,
                  ClassTable& class_table_out, VoxelLeafWoodTable& voxel_lw_out,
                  const std::string& leaf_classes_str, const std::string& wood_classes_str,
                  double lambda1 = 0.0,
                  bool ppl_enabled = false, std::vector<PplHit>* ppl_out = nullptr)
{
  // Determine the final number of threads to use
  size_t resolved_threads = num_threads;
  if (resolved_threads == 0) { // Auto-detect
    resolved_threads = std::thread::hardware_concurrency();
    if (resolved_threads == 0) resolved_threads = 1; // Fallback
  }

  // Common parameters for all processors
  double tan_half_divergence = calc_beam_metrics ? tan(0.5 * beam_divergence) : 0.0;
  const std::vector<double>* peaks_ptr = apply_flat_top ? &grid.getPeaks() : nullptr;

  // Determine the per-point passthrough stride from the ray cloud's extra-byte header.
  ray::CloudReader reader;
  reader.begin(cloud_name);
  const ray::LasHeader &hdr = reader.header();
  const std::vector<uint8_t> extra_bytes_vlr = hdr.sensorExtraVlr();
  const uint16_t stride = static_cast<uint16_t>(kPassthroughStdBytes + hdr.sensorExtraSize());

  // Resolve leaf/wood class field sources for foliage_class tagging
  auto split_field_codes_ = [](const std::string& s, std::string& field, std::string& codes) {
      auto colon = s.find(':');
      if (colon != std::string::npos) { field = s.substr(0, colon); codes = s.substr(colon + 1); }
      else { field.clear(); codes = s; }
  };
  std::string leaf_field_, leaf_codes_, wood_field_, wood_codes_;
  split_field_codes_(leaf_classes_str, leaf_field_, leaf_codes_);
  split_field_codes_(wood_classes_str, wood_field_, wood_codes_);
  const ClassFieldSource leaf_src = resolveClassField(leaf_field_, hdr);
  const ClassFieldSource wood_src = resolveClassField(wood_field_, hdr);
  std::set<int> leaf_set, wood_set;
  {
      std::stringstream ss(leaf_codes_);
      std::string item;
      while (std::getline(ss, item, ',')) { try { leaf_set.insert(std::stoi(item)); } catch (...) {} }
  }
  {
      std::stringstream ss(wood_codes_);
      std::string item;
      while (std::getline(ss, item, ',')) { try { wood_set.insert(std::stoi(item)); } catch (...) {} }
  }
  const bool no_foliage_classes = leaf_set.empty() && wood_set.empty();
  auto resolveFoliageClass = [&](const uint8_t* base, uint8_t cls, bool is_ground) -> uint8_t {
      if (no_foliage_classes) {
          // No --leaf_classes/--wood_classes given: treat every non-ground return as plant,
          // regardless of its class code. Ground is identified via the DTM (or --dtm_from_class)
          // and excluded; PAD is then derived from the analytic --lad G-function downstream.
          return is_ground ? 0 : 1;
      }
      if (cls < 3) return 0;
      if (leaf_set.count(readClassValue(base, leaf_src))) return 2;
      if (wood_set.count(readClassValue(base, wood_src))) return 3;
      return 1;
  };

  if (resolved_threads > 1) {
    // --- OPTION 2: Parallel Producer-Consumer Implementation ---

    // Scale queue depth to available RAM: use up to 1 GB, minimum 8 batches/thread.
    // Apply a 75% margin to the queried figure so the queue and grid share one budget.
    const size_t avail_ram    = ray::queryAvailableMemoryBytes() * 3 / 4;
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

    // RAM-fit selection (see generateVoxelGrid) guarantees a flat grid on the in-process path.
    assert(grid.isFlat());
    VoxelGrid::Voxel* flat_ptr = &grid.voxelAt(0);
    const int64_t flat_dimX  = grid.getDimensions()[0];
    const int64_t flat_dimXY = grid.getDimensions()[0] * grid.getDimensions()[1];

    std::vector<ClassTable> worker_class_tables(resolved_threads);
    std::vector<VoxelLeafWoodTable> worker_voxel_lw(resolved_threads);
    std::mutex ppl_sink_mutex;

    auto worker_task = [&](size_t thread_idx) {
      VoxelProcessor processor(grid.getBounds(), grid.getVoxelWidth(), weighting_method, use_occlusion,
                               apply_flat_top, peaks_ptr, calc_beam_metrics, beam_diameter,
                               tan_half_divergence, subvoxel_split, dtm, lambda1);
      processor.setFlatTarget(flat_ptr, flat_dimX, flat_dimXY);
      if (ppl_enabled) processor.enablePpl();
      BeamBatch batch;
      while (beam_queue.pop(batch)) {
        for (size_t b = 0; b < batch.count; ++b) {
          processor.processBeam(batch.beams[b]);
        }
      }
      worker_class_tables[thread_idx] = processor.extractClassTable();
      worker_voxel_lw[thread_idx]     = processor.extractVoxelLW();
      if (ppl_enabled && ppl_out) {
        std::vector<PplHit> h = processor.extractPplHits();
        std::lock_guard<std::mutex> lk(ppl_sink_mutex);
        ppl_out->insert(ppl_out->end(), std::make_move_iterator(h.begin()), std::make_move_iterator(h.end()));
      }
    };

    for (size_t i = 0; i < resolved_threads; ++i)
      threads.emplace_back(worker_task, i);

    size_t num_bounded = 0;
    std::vector<uint8_t> passthrough;
    std::vector<int32_t> beam_ids_chunk;
    static bool not_raycloud_warned = false;
    PulseGrouper grouper;

    // Current batch being filled by the producer.
    BeamBatch current_batch;

    // Commit a completed beam into the current batch; push when the batch is full.
    auto emit_beam = [&](BeamData&& bd) {
      current_batch.beams[current_batch.count++] = std::move(bd);
      if (current_batch.count == kBeamBatchSize) {
        BeamBatch tmp = current_batch;  // copy before reset so workers get valid data
        current_batch.count = 0;
        beam_queue.push(std::move(tmp));
      }
    };

    reader.read(
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
          const bool is_ground = isGroundHit(pd.x, pd.y, pd.z, pd.classification, dtm_from_class, dtm, dtm_filter_distance);
          if (is_ground)
            pd.bound = 0;
          const size_t base_i = i * stride;
          if (passthrough.size() >= base_i + stride)
            pd.foliage_class = resolveFoliageClass(&passthrough[base_i], pd.classification, is_ground);
          grouper.add(pd, emit_beam);
        }
        passthrough.clear();
        beam_ids_chunk.clear();
      }, num_bounded, 255.0, nullptr, las_chunk,
         nullptr, &passthrough, nullptr, &beam_ids_chunk);

    grouper.flush(emit_beam);  // emit any pulses still incomplete at end-of-stream
    if (current_batch.count > 0) {
      beam_queue.push(std::move(current_batch));
    }
    beam_queue.notify_done();
    for (auto& t : threads) { t.join(); }

    for (size_t i = 0; i < resolved_threads; ++i) {
        for (auto& kv : worker_class_tables[i])
            for (int c = 0; c < 256; ++c)
                class_table_out[kv.first][c] += kv.second[c];
        for (auto& kv : worker_voxel_lw[i]) {
            voxel_lw_out[kv.first].first  += kv.second.first;
            voxel_lw_out[kv.first].second += kv.second.second;
        }
    }
    std::cout << "Parallel processing finished." << std::endl;

  } else {
    // --- OPTION 1: Single-Threaded Implementation ---
    std::cout << "Processing point cloud using 1 thread (in-memory)..." << std::endl;
    VoxelProcessor processor(grid.getBounds(), grid.getVoxelWidth(), weighting_method, use_occlusion,
                             apply_flat_top, peaks_ptr, calc_beam_metrics, beam_diameter,
                             tan_half_divergence, subvoxel_split, dtm, lambda1);
    // RAM-fit selection (see generateVoxelGrid) guarantees a flat grid on the in-process path.
    assert(grid.isFlat());
    processor.setFlatTarget(&grid.voxelAt(0),
                            grid.getDimensions()[0],
                            grid.getDimensions()[0] * grid.getDimensions()[1]);
    if (ppl_enabled) processor.enablePpl();

    size_t num_bounded = 0;
    std::vector<uint8_t> passthrough;
    std::vector<int32_t> beam_ids_chunk;
    static bool not_raycloud_warned = false;
    PulseGrouper grouper;
    auto emit_beam = [&](BeamData&& beam) { processor.processBeam(beam); };
    reader.read(
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
          const bool is_ground = isGroundHit(pd.x, pd.y, pd.z, pd.classification, dtm_from_class, dtm, dtm_filter_distance);
          if (is_ground)
            pd.bound = 0;
          const size_t base_i = i * stride;
          if (passthrough.size() >= base_i + stride)
            pd.foliage_class = resolveFoliageClass(&passthrough[base_i], pd.classification, is_ground);
          grouper.add(pd, emit_beam);
        }
        passthrough.clear();
        beam_ids_chunk.clear();
      }, num_bounded, 255.0, nullptr, 1000000, nullptr, &passthrough, nullptr, &beam_ids_chunk);

    grouper.flush(emit_beam);
    // Flat path: writes already landed in the grid's flat array — nothing to move.
    class_table_out = processor.extractClassTable();
    voxel_lw_out    = processor.extractVoxelLW();
    if (ppl_enabled && ppl_out) {
      std::vector<PplHit> h = processor.extractPplHits();
      ppl_out->insert(ppl_out->end(), std::make_move_iterator(h.begin()), std::make_move_iterator(h.end()));
    }
  }

  return true;
}

} // anonymous namespace

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


namespace {

// File-local helpers for the out-of-core path (formerly OutOfCoreStrategy members).
bool createShards(const std::string& cloud_name, VoxelGrid& grid, size_t num_threads, size_t ram_budget_mb,
                  const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                  bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                  const HeightField* dtm, int dtm_from_class, double dtm_filter_distance,
                  ClassTable& class_table_out, VoxelLeafWoodTable& voxel_lw_out,
                  const std::string& leaf_classes_str, const std::string& wood_classes_str,
                  std::vector<std::string>& out_shard_paths, double lambda1 = 0.0);

bool mergeShards(const std::vector<std::string>& shard_paths, VoxelGrid& grid);

// Out-of-core traversal: shard to disk, then k-way merge into the grid.
bool runOutOfCore(const std::string& cloud_name, VoxelGrid& grid, size_t num_threads, size_t ram_budget_mb,
                  const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                  bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                  const HeightField* dtm, int dtm_from_class, double dtm_filter_distance,
                  ClassTable& class_table_out, VoxelLeafWoodTable& voxel_lw_out,
                  const std::string& leaf_classes_str, const std::string& wood_classes_str,
                  double lambda1 = 0.0) {
    std::vector<std::string> shard_paths;
    std::cout << "Starting out-of-core processing..." << std::endl;

    if (!createShards(cloud_name, grid, num_threads, ram_budget_mb, weighting_method, use_occlusion, apply_flat_top,
                      calc_beam_metrics, beam_diameter, beam_divergence, subvoxel_split, dtm,
                      dtm_from_class, dtm_filter_distance, class_table_out, voxel_lw_out,
                      leaf_classes_str, wood_classes_str, shard_paths, lambda1)) {
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

bool createShards(const std::string& cloud_name, VoxelGrid& grid, size_t num_threads, size_t ram_budget_mb,
                  const std::string& weighting_method, bool use_occlusion, bool apply_flat_top,
                  bool calc_beam_metrics, double beam_diameter, double beam_divergence, int subvoxel_split,
                  const HeightField* dtm, int dtm_from_class, double dtm_filter_distance,
                  ClassTable& class_table_out, VoxelLeafWoodTable& voxel_lw_out,
                  const std::string& leaf_classes_str, const std::string& wood_classes_str,
                  std::vector<std::string>& out_shard_paths, double lambda1) {
    std::cout << "Phase 1: Processing points and writing to temporary shards..." << std::endl;

    size_t resolved_threads = num_threads == 0 ? std::thread::hardware_concurrency() : num_threads;
    if (resolved_threads == 0) resolved_threads = 1;

    size_t ram_per_thread_bytes = ram_budget_mb * 1024 * 1024;
    // Estimate size of a map entry: Key + Value + overhead (approx 2 pointers) + map node.
    size_t voxel_pair_size_approx = sizeof(VoxelCoord) + sizeof(VoxelGrid::Voxel) + sizeof(void*)*2 + sizeof(std::pair<U8,float>)*2;
    size_t max_voxels_per_thread = (ram_per_thread_bytes / voxel_pair_size_approx);

    std::cout << "Out-of-core settings: " << resolved_threads << " threads, " << ram_budget_mb << "MB RAM budget per thread." << std::endl;
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

    std::vector<ClassTable> worker_class_tables(resolved_threads);
    std::vector<VoxelLeafWoodTable> worker_voxel_lw(resolved_threads);

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

      // NB: class tables are NOT cleared on flush — they persist the full run
      worker_class_tables[thread_id] = processor.extractClassTable();
      worker_voxel_lw[thread_id]     = processor.extractVoxelLW();
    };

    for (size_t i = 0; i < resolved_threads; ++i) {
      threads.emplace_back(worker_task, i);
    }

    // --- Producer loop ---
    // Determine the per-point passthrough stride from the ray cloud's extra-byte header.
    ray::CloudReader reader;
    reader.begin(cloud_name);
    const ray::LasHeader &hdr = reader.header();
    const std::vector<uint8_t> extra_bytes_vlr = hdr.sensorExtraVlr();
    const uint16_t stride = static_cast<uint16_t>(kPassthroughStdBytes + hdr.sensorExtraSize());

    // Resolve leaf/wood class field sources for foliage_class tagging
    auto split_field_codes_ = [](const std::string& s, std::string& field, std::string& codes) {
        auto colon = s.find(':');
        if (colon != std::string::npos) { field = s.substr(0, colon); codes = s.substr(colon + 1); }
        else { field.clear(); codes = s; }
    };
    std::string leaf_field_, leaf_codes_, wood_field_, wood_codes_;
    split_field_codes_(leaf_classes_str, leaf_field_, leaf_codes_);
    split_field_codes_(wood_classes_str, wood_field_, wood_codes_);
    const ClassFieldSource leaf_src = resolveClassField(leaf_field_, hdr);
    const ClassFieldSource wood_src = resolveClassField(wood_field_, hdr);
    std::set<int> leaf_set, wood_set;
    {
        std::stringstream ss(leaf_codes_);
        std::string item;
        while (std::getline(ss, item, ',')) { try { leaf_set.insert(std::stoi(item)); } catch (...) {} }
    }
    {
        std::stringstream ss(wood_codes_);
        std::string item;
        while (std::getline(ss, item, ',')) { try { wood_set.insert(std::stoi(item)); } catch (...) {} }
    }
    const bool no_foliage_classes = leaf_set.empty() && wood_set.empty();
    auto resolveFoliageClass = [&](const uint8_t* base, uint8_t cls, bool is_ground) -> uint8_t {
        if (no_foliage_classes) {
            // No --leaf_classes/--wood_classes given: treat every non-ground return as plant,
            // regardless of its class code. Ground is identified via the DTM (or --dtm_from_class)
            // and excluded; PAD is then derived from the analytic --lad G-function downstream.
            return is_ground ? 0 : 1;
        }
        if (cls < 3) return 0;
        if (leaf_set.count(readClassValue(base, leaf_src))) return 2;
        if (wood_set.count(readClassValue(base, wood_src))) return 3;
        return 1;
    };

    size_t num_bounded = 0;
    std::vector<uint8_t> passthrough;
    std::vector<int32_t> beam_ids_chunk;
    static bool not_raycloud_warned = false;
    PulseGrouper grouper;
    auto emit_beam = [&](BeamData&& beam) { beam_queue.push(std::move(beam)); };
    reader.read(
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
          const bool is_ground = isGroundHit(pd.x, pd.y, pd.z, pd.classification, dtm_from_class, dtm, dtm_filter_distance);
          if (is_ground)
            pd.bound = 0;
          const size_t base_i = i * stride;
          if (passthrough.size() >= base_i + stride)
            pd.foliage_class = resolveFoliageClass(&passthrough[base_i], pd.classification, is_ground);
          grouper.add(pd, emit_beam);
        }
        passthrough.clear();
        beam_ids_chunk.clear();
      }, num_bounded, 255.0, nullptr, 1000000, nullptr, &passthrough, nullptr, &beam_ids_chunk);

    grouper.flush(emit_beam); // emit any pulses still incomplete at end-of-stream
    beam_queue.notify_done();
    for (auto& t : threads) { t.join(); }

    for (size_t i = 0; i < resolved_threads; ++i) {
        for (auto& kv : worker_class_tables[i])
            for (int c = 0; c < 256; ++c)
                class_table_out[kv.first][c] += kv.second[c];
        for (auto& kv : worker_voxel_lw[i]) {
            voxel_lw_out[kv.first].first  += kv.second.first;
            voxel_lw_out[kv.first].second += kv.second.second;
        }
    }

    out_shard_paths = std::move(shard_paths_list);
    return true;
}

bool mergeShards(const std::vector<std::string>& shard_paths, VoxelGrid& grid) {
    ShardMerger merger(shard_paths, grid);
    return merger.merge();
}

// Exact potential-path-length (PPL) attenuation MLE, mirroring AMAPVox VoxelizationTask:
// per voxel solve for k such that  Σ_hits bsIn·L·e^{-kL}/(1-e^{-kL}) = ppl_miss_wL.
// The LHS is strictly decreasing in k (→∞ as k→0⁺, →0 as k→∞), so the root is unique; bisect on
// [0, kMaxAtt] and cap at kMaxAtt. Records are grouped by flat voxel index; result → Voxel::ppl_lambda.
void solvePplExact(VoxelGrid& grid, std::vector<PplHit>& hits)
{
  if (!grid.isFlat() || hits.empty()) return;
  const double kMaxAtt = 20.0;   // AMAPVox maximal-attenuation cap
  const double eps = 1e-12;
  std::unordered_map<int64_t, std::vector<std::pair<float, float>>> by_voxel;  // flat → [(L, bsIn)]
  by_voxel.reserve(hits.size() / 2 + 1);
  for (const PplHit& h : hits) by_voxel[h.voxel_index].push_back({ h.L, h.bsIn });

  for (auto& kv : by_voxel) {
    const std::vector<std::pair<float, float>>& recs = kv.second;
    if (recs.empty()) continue;
    VoxelGrid::Voxel& v = grid.voxelAt(kv.first);
    const double miss = static_cast<double>(v.ppl_miss_wL);
    if (miss <= eps) { v.ppl_lambda = static_cast<float>(kMaxAtt); continue; }  // no exiting beam → max attenuation
    auto f = [&](double k) -> double {
      double s = 0.0;
      for (const auto& r : recs) {
        const double L = r.first, bs = r.second;
        const double e = std::exp(k * L) - 1.0;
        if (e > eps) s += bs * L / e;
      }
      return s;
    };
    if (f(kMaxAtt) >= miss) { v.ppl_lambda = static_cast<float>(kMaxAtt); continue; }  // root beyond cap
    double lo = eps, hi = kMaxAtt;
    for (int it = 0; it < 100; ++it) {
      const double mid = 0.5 * (lo + hi);
      if (f(mid) > miss) lo = mid; else hi = mid;   // f decreasing: too-large f ⇒ move right
    }
    v.ppl_lambda = static_cast<float>(0.5 * (lo + hi));
  }
}

} // anonymous namespace

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
      ray::CloudReader reader;
      reader.begin(cloud_name);
      const ray::LasHeader &pre_hdr = reader.header();
      const bool file_has_bound = pre_hdr.has("bound");
      size_t num_bounded = 0;
      return reader.read(
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
        ray::CloudReader reader;
        reader.begin(params.cloud_name);
        const ray::LasHeader &hdr = reader.header();
        const std::vector<uint8_t> extra_bytes_vlr = hdr.sensorExtraVlr();
        const uint16_t stride = static_cast<uint16_t>(kPassthroughStdBytes + hdr.sensorExtraSize());

        // Chunked pass collecting ground points whose passthrough classification matches.
        size_t num_bounded = 0;
        std::vector<uint8_t> passthrough;
        if (!reader.read(
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
            }, num_bounded, 255.0, nullptr, 1000000, nullptr, &passthrough)) {
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

    // === RAM-FIT PATH SELECTION ===
    // Single decision: does the flat (dense) voxel array fit in a safe fraction of available RAM?
    // If not (or if the user forced --out_of_core), use the out-of-core path. The in-process path
    // always operates on a flat grid; the in-RAM sparse fallback has been removed.
    // Dimensions match the VoxelGrid ctor: ceil(extent / voxel_width) per axis.
    const Eigen::Vector3d extent = processing_bounds.max_bound_ - processing_bounds.min_bound_;
    const Eigen::Matrix<int64_t, 3, 1> dims = (extent / params.voxel_size).array().ceil().cast<int64_t>();

    // rayvoxel keeps voxels perfectly cubic (no rescaling to bounds): the grid's true max corner is
    // snapped UP to the next multiple of voxel_size so the last voxels fully enclose their points.
    // Report that voxel-aligned extent so it matches the .vox header (which uses the actual res).
    {
      const Eigen::Vector3d ub_extent = user_bounds.max_bound_ - user_bounds.min_bound_;
      const Eigen::Matrix<int64_t, 3, 1> ub_dims = (ub_extent / params.voxel_size).array().ceil().cast<int64_t>();
      const Eigen::Vector3d aligned_max = user_bounds.min_bound_ + ub_dims.cast<double>() * params.voxel_size;
      std::cout << "Voxel grid is cubic at " << params.voxel_size << " m; max corner aligned to voxel size: ("
                << aligned_max.x() << ", " << aligned_max.y() << ", " << aligned_max.z()
                << ")  (requested max (" << user_bounds.max_bound_.x() << ", "
                << user_bounds.max_bound_.y() << ", " << user_bounds.max_bound_.z() << "))" << std::endl;
    }

    const int64_t flat_cells = dims[0] * dims[1] * dims[2];
    const size_t flat_bytes = static_cast<size_t>(flat_cells) * sizeof(VoxelGrid::Voxel);
    const size_t avail_ram_bytes = ray::queryAvailableMemoryBytes();
    const bool fits_flat = (flat_bytes <= avail_ram_bytes * 3 / 4);
    const bool use_ooc_effective = params.use_ooc || !fits_flat;
    std::cout << "Path selection: " << (use_ooc_effective ? "out-of-core" : "flat (in-memory)")
              << "; flat grid = " << flat_bytes / (1024 * 1024) << " MB, available RAM = "
              << avail_ram_bytes / (1024 * 1024) << " MB"
              << (params.use_ooc ? " (--out_of_core forced)" : "") << std::endl;

    // Handle auto-detection of reservation size here. The point count was already gathered
    // by the combined pre-scan pass above (need_reserve_count mirrors that pass's trigger).
    size_t final_reserve_size = params.reserve_size;
    if (need_reserve_count) {
        final_reserve_size += pre_scan.point_count;
    }
    // For OOC, we never reserve in the final grid, as it's populated at the end.
    if (use_ooc_effective) { final_reserve_size = 0; }

    std::unique_ptr<VoxelGrid> grid_ptr;
    try {
        // In-process path: force flat allocation by giving the ctor a budget >= the flat size.
        // OOC path: keep the user's per-grid RAM budget (sparse fallback is acceptable there).
        const size_t user_budget_bytes = static_cast<size_t>(params.ram_budget_mb) * 1024 * 1024;
        const size_t ram_budget_bytes = use_ooc_effective
            ? user_budget_bytes
            : std::max(flat_bytes, user_budget_bytes);
        grid_ptr = std::make_unique<VoxelGrid>(processing_bounds, params.voxel_size,
                                               ram_budget_bytes, final_reserve_size,
                                               params.apply_flat_top);
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

    // === TRAVERSAL ===
    // ClassTable and VoxelLeafWoodTable are accumulated during the traversal pass.
    ClassTable class_table;
    VoxelLeafWoodTable voxel_lw;
    // Exact PPL is supported on the flat (in-memory) path only; OOC falls back to the
    // closed-form mean-field estimator in computeLambda("ppl").
    const bool ppl_exact = !use_ooc_effective &&
        (std::find(params.attenuation_methods.begin(), params.attenuation_methods.end(),
                   std::string("ppl")) != params.attenuation_methods.end());
    std::vector<PplHit> ppl_hits;
    bool processing_success;
    if (use_ooc_effective) {
        processing_success = runOutOfCore(params.cloud_name, grid, params.num_threads, params.ram_budget_mb,
                                          params.weighting_method, params.use_occlusion, params.apply_flat_top,
                                          params.calc_beam_metrics, beam_diameter, beam_divergence, params.subvoxel_split, dtm_ptr.get(),
                                          params.dtm_from_class, params.dtm_filter_distance,
                                          class_table, voxel_lw,
                                          params.leaf_classes_str, params.wood_classes_str, lambda1);
    } else {
        processing_success = runInProcess(params.cloud_name, grid, params.num_threads,
                                          params.weighting_method, params.use_occlusion, params.apply_flat_top,
                                          params.calc_beam_metrics, beam_diameter, beam_divergence, params.subvoxel_split, dtm_ptr.get(),
                                          params.dtm_from_class, params.dtm_filter_distance,
                                          class_table, voxel_lw,
                                          params.leaf_classes_str, params.wood_classes_str, lambda1,
                                          ppl_exact, &ppl_hits);
    }

    if (!processing_success) {
      return false; // Traversal failed, exit early.
    }

    // === OPTIONAL UNBOUND-RAY PASS ===
    // Traverse a second cloud of unbound (miss) rays into the same grid, before any post-processing.
    // Grid bounds are NOT recomputed — unbound rays clip naturally at the grid boundary via the
    // existing walkGrid bounds check. Unbound rays produce no hits, so num_hits / IAD / class
    // outputs are unaffected; only num_beams / path_length (and free-path) accumulate.
    if (!params.unbound_file.empty()) {
        // The second pass must accumulate (+=) into the already-populated grid. Only the flat
        // in-process path accumulates; the OOC merge writes by assignment and would clobber the
        // primary results, so the unbound pass is supported on the flat grid only.
        if (!grid.isFlat()) {
            std::cerr << "Warning: --unbound_file is only supported on the flat (in-memory) path; "
                         "skipping unbound traversal for: " << params.unbound_file << std::endl;
        } else {
            std::cout << "Traversing unbound-ray file: " << params.unbound_file << std::endl;
            // Throwaway accumulators — unbound rays produce no hits; classification is bound-file-only.
            ClassTable unbound_class_table;
            VoxelLeafWoodTable unbound_voxel_lw;
            bool unbound_ok = runInProcess(params.unbound_file, grid, params.num_threads,
                                           params.weighting_method, params.use_occlusion, params.apply_flat_top,
                                           params.calc_beam_metrics, beam_diameter, beam_divergence, params.subvoxel_split, dtm_ptr.get(),
                                           params.dtm_from_class, params.dtm_filter_distance,
                                           unbound_class_table, unbound_voxel_lw,
                                           params.leaf_classes_str, params.wood_classes_str, lambda1,
                                           ppl_exact, &ppl_hits);
            if (!unbound_ok) {
                std::cerr << "Warning: --unbound_file traversal failed for: " << params.unbound_file << std::endl;
                // Non-fatal — primary output is still valid.
            }
        }
    }

    // === EXACT PPL SOLVE ===
    // Solve the per-voxel potential-path-length MLE from the intercepted-beam records gathered
    // during traversal (bound + unbound passes). Results land in Voxel::ppl_lambda; the writer
    // prefers it over the closed-form fallback. Records are freed immediately after.
    if (ppl_exact) {
        std::cout << "Solving exact PPL attenuation (" << ppl_hits.size()
                  << " intercepted records)..." << std::endl;
        solvePplExact(grid, ppl_hits);
        std::vector<PplHit>().swap(ppl_hits);
    }

    // === POST-PROCESSING AND OUTPUT (Same for all strategies) ===
    if (params.neighbour_prior_min_rays > 0) {
        std::cout << "Applying neighbour priors (min rays = " << params.neighbour_prior_min_rays << ")..." << std::endl;
        applyNeighbourPriors(grid, params.neighbour_prior_min_rays);
    }

    // The classification table and per-voxel leaf/wood counts are already populated by the
    // strategy's traversal pass. Inclination angle distributions still require their own
    // O(N_points) KNN pass when requested.
    PerTreeIadMap per_tree_iad;
    PredominantTreeTable predominant_tree;
    if (params.calc_inclination_dist) {
      {
        static bool empty_classes_warned = false;
        if (!empty_classes_warned && params.leaf_classes_str.empty() && params.wood_classes_str.empty()) {
          std::cerr << "Warning: --inclination_dist active but neither --leaf_classes nor --wood_classes "
                       "is set; LAD and WAD will be zero (PAD is unaffected)." << std::endl;
          empty_classes_warned = true;
        }
      }
      std::cout << "Building inclination angle distributions..." << std::endl;
      buildClassAndIadTable(params.cloud_name, grid, params, dtm_ptr.get(),
                            per_tree_iad, predominant_tree);
    }

    std::cout << "Calculating output metrics..." << std::endl;
    MetricResultsMap metrics = calculateOutputMetrics(grid, params, dtm_ptr.get(), class_table, per_tree_iad, predominant_tree, voxel_lw);

    // Optional woody-volume rasterisation: project the trees.txt branch cylinders into the grid and
    // annotate the per-voxel output with wood_volume (m^3) and wood_volume_density (m^3/m^3). This is
    // independent of ray traversal; it only annotates voxels already present in @c metrics (observed,
    // or all with --write_empty), which covers scanned trees in practice.
    if (!params.trees_file.empty()) {
        ray::ForestStructure forest;
        if (!forest.load(params.trees_file)) {
            std::cerr << "Error: could not load --trees file: " << params.trees_file << std::endl;
            return false;
        }
        const double voxel_volume = params.voxel_size * params.voxel_size * params.voxel_size;
        auto wood_map = computeWoodVolumePerVoxel(forest, grid.getBounds(), grid.getVoxelWidth(),
                                                  grid.getDimensions());
        double total_voxel_wood = 0.0, total_tree_volume = 0.0;
        for (const auto &tree : forest.trees) total_tree_volume += tree.volume();
        for (const auto &kv : wood_map) {
            total_voxel_wood += kv.second;
            auto it = metrics.find(kv.first);
            if (it != metrics.end()) {
                it->second.wood_volume = kv.second;
                it->second.wood_volume_density = voxel_volume > 0.0 ? kv.second / voxel_volume : 0.0;
            }
        }
        std::cout << "Woody volume: rasterised " << total_voxel_wood << " m^3 across " << wood_map.size()
                  << " voxels (tree-file total branch volume " << total_tree_volume << " m^3)." << std::endl;
    }

    // Pass the pre-calculated metrics to the writer functions. All rayvoxel outputs share a
    // "_voxel" stub, e.g. {stub}_voxel.vox, {stub}_voxel_filled.vox, {stub}_voxel_iad.csv.
    std::string base_name_stub = getFileNameStub(params.cloud_name) + "_voxel";

    // Per-tree inclination angle distributions sidecar. Written whenever tree_id is present in
    // the input (per_tree_iad non-empty); --output_iad is not required. The CSV describes the
    // cloud's trees, so it is written once for the base stub regardless of voxel filtering
    // (_filled / _include_empty variants reuse the same per-tree data).
    if (params.calc_inclination_dist && !per_tree_iad.empty()) {
        ray::CloudReader iad_reader;
        iad_reader.begin(params.cloud_name);
        const ray::LasHeader &iad_hdr = iad_reader.header();
        const bool has_stem_id = iad_hdr.has("stem_id");
        writePerTreeIadCsv(base_name_stub, per_tree_iad, params, has_stem_id);
    }

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
