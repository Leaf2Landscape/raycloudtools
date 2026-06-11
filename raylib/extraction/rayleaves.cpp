// Copyright (c) 2023
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayleaves.h"
#include "../rayparse.h"
#include <nabo/nabo.h>
#include "../raycuboid.h"
#include "../rayforeststructure.h"
#include "../raymesh.h"
#include "../rayply.h"
#include "../rayrenderer.h"
#define STB_IMAGE_IMPLEMENTATION
#include <random>
#include <cmath>
#include <functional>
#include <string>
#include "raylib/imageread.h"
#include "raylib/rayvoxel/rayvox.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace ray
{
typedef std::complex<float> Cmp;


// Function to get the leaf angle distribution based on user input
std::function<double(double)> getLeafAngleDistribution(int distribution)
{
  switch (distribution)
  {
  case 1: // Uniform distribution
    return [](double) { return 1.0; };
  case 2: // Spherical distribution
    return [](double angle) { return std::sin(angle * M_PI / 180.0); };
  case 3: // Erectophile distribution
    return [](double angle) { return std::sin(2 * angle * M_PI / 180.0); };
  case 4: // Plagiophile distribution
    return [](double angle) { return std::sin(4 * angle * M_PI / 180.0); };
  case 5: // Planophile distribution
    return [](double angle) { return std::cos(angle * M_PI / 180.0); };
  case 6: // Extremophile distribution
    return [](double angle) { return std::abs(std::cos(angle * M_PI / 180.0)); };
  default:
    return [](double) { return 1.0; };  // Default to uniform distribution
  }
}

bool generateLeaves(const std::string &cloud_name, const std::string &trees_file, const std::string &leaf_file,
                    double leaf_area, double droop, int distribution, double leafAreaDensity, bool stalks,
                    const std::string &vox_file, const std::string &rayvoxel_method,
                    const std::string &leaf_classes_str)
{
  // Build leaf-point predicate from optional --leaf_classes code list (e.g. "1,2,3").
  // Strip any "field:" prefix to match rayvoxel syntax. Without --leaf_classes: alpha != 0.
  std::function<bool(uint8_t)> is_leaf;
  {
    std::string codes_str = leaf_classes_str;
    auto colon = codes_str.find(':');
    if (colon != std::string::npos)
      codes_str = codes_str.substr(colon + 1);
    if (codes_str.empty())
    {
      is_leaf = [](uint8_t a) { return a != 0; };
    }
    else
    {
      std::vector<uint8_t> codes;
      std::stringstream ss(codes_str);
      std::string tok;
      while (std::getline(ss, tok, ','))
      {
        tok.erase(0, tok.find_first_not_of(" \t"));
        tok.erase(tok.find_last_not_of(" \t") + 1);
        if (!tok.empty())
          codes.push_back(static_cast<uint8_t>(std::stoi(tok)));
      }
      is_leaf = [codes](uint8_t a) {
        return std::find(codes.begin(), codes.end(), a) != codes.end();
      };
    }
  }

  // For now we assume that woody points have been set as unbounded (alpha=0). e.g. through raycolour foliage or
  // raysplit file distance 0.2 as examples. so firstly we must calculate the foliage density across the whole map.
  const std::string cloud_stub = getFileNameStub(cloud_name);
  Cloud::Info info;
  if (!Cloud::getInfo(cloud_name, info))
  {
    return false;
  }
  const Cuboid bounds = info.ends_bound;
  const Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
  const double vox_width = 1;
  Eigen::Vector3i dims =
    (extent / vox_width).cast<int>() + Eigen::Vector3i(2, 2, 2);  // so that we have extra space to convolve
  Cuboid grid_bounds = bounds;
  grid_bounds.min_bound_ -= Eigen::Vector3d(vox_width, vox_width, vox_width);
  DensityGrid grid(grid_bounds, vox_width, dims);
  grid.calculateDensities(cloud_name);
  grid.addNeighbourPriors();

  // Optional per-voxel leaf area density (LAD) and leaf inclination angle distribution (LIAD)
  // loaded from a rayvoxel file. Supports both AMAPVox .vox format (VOXEL SPACE header, nbEchos)
  // and the extended .txt format (num_hits, subvoxel_bitmap, lad_*, liad_*).
  struct VoxKey {
    long i, j, k;
    bool operator==(const VoxKey& o) const { return i==o.i && j==o.j && k==o.k; }
  };
  struct VoxKeyHash {
    size_t operator()(const VoxKey& v) const {
      size_t h = std::hash<long>{}(v.i);
      h ^= std::hash<long>{}(v.j) + 0x9e3779b9 + (h<<6) + (h>>2);
      h ^= std::hash<long>{}(v.k) + 0x9e3779b9 + (h<<6) + (h>>2);
      return h;
    }
  };
  // nb_echos: total echoes in this voxel from the full lidar scan (used for fraction computation)
  struct VoxLeafData { double lad; std::vector<double> liad; uint64_t bitmap = 0; int nb_echos = 0; };
  using VoxMap = std::unordered_map<VoxKey, VoxLeafData, VoxKeyHash>;

  VoxMap vox_map;
  Eigen::Vector3d vox_min(0,0,0), vox_res(1,1,1);
  int n_iad_bins = 0;
  // per-DensityGrid-cell LIAD accumulator (keyed by flat 1m-grid index)
  std::unordered_map<int, std::vector<double>> cell_liad_sum;
  bool has_lad_vox = false;
  bool has_liad_vox = false;
  bool has_bitmap_vox = false;

  if (!vox_file.empty())
  {
    // Detect format: AMAPVox "VOXEL SPACE" vs extended text
    bool is_amap_format = false;
    {
      std::ifstream probe(vox_file);
      std::string first_line;
      if (std::getline(probe, first_line))
      {
        auto s = first_line.find_first_not_of(" \t\r\n");
        auto e = first_line.find_last_not_of(" \t\r\n");
        std::string trimmed = (s != std::string::npos) ? first_line.substr(s, e - s + 1) : "";
        is_amap_format = (trimmed == "VOXEL SPACE");
      }
    }

    auto notify_columns = [&]() {
      if (!has_lad_vox)
        std::cout << "Note: rayvoxel file has no lad column for method '" << rayvoxel_method
                  << "'; using --leaf_density as fallback." << std::endl;
      if (!has_liad_vox)
        std::cout << "Note: rayvoxel file has no liad_* columns; using --leaf_angle distribution as fallback." << std::endl;
      if (has_lad_vox)
        std::cout << "Note: --rayvoxel active; --leaf_density used only as fallback for uncovered voxels." << std::endl;
      if (has_liad_vox)
        std::cout << "Note: --rayvoxel active; --leaf_angle used only as fallback when liad data is absent." << std::endl;
      if (has_bitmap_vox)
        std::cout << "Note: subvoxel_bitmap found; using bitmap-guided leaf placement with three-state occupancy." << std::endl;
    };

    if (is_amap_format)
    {
      // ---- AMAPVox VOXEL SPACE format ----
      ray::VoxelSpace space;
      if (!ray::readVox(vox_file, space))
      {
        std::cerr << "Error: cannot read rayvoxel file: " << vox_file << std::endl;
        return false;
      }

      auto parse_vec3 = [](const std::string& s, Eigen::Vector3d& out) {
        std::istringstream ss(s);
        return static_cast<bool>(ss >> out.x() >> out.y() >> out.z());
      };
      auto it_min = space.header.find("min_corner");
      auto it_res = space.header.find("res");
      if (it_min == space.header.end() || it_res == space.header.end() ||
          !parse_vec3(it_min->second, vox_min) || !parse_vec3(it_res->second, vox_res))
      {
        std::cerr << "Error: rayvoxel file missing min_corner or res header." << std::endl;
        return false;
      }

      auto it_col = space.header.find("colnames");
      std::vector<std::string> colnames;
      if (it_col != space.header.end())
      {
        std::istringstream ss(it_col->second);
        std::string tok;
        while (ss >> tok) colnames.push_back(tok);
      }
      auto col_idx = [&](const std::string& name) -> int {
        for (int c = 3; c < (int)colnames.size(); ++c)
          if (colnames[c] == name) return c - 3;
        return -1;
      };

      int lad_col = col_idx("lad_" + rayvoxel_method);
      if (lad_col < 0) lad_col = col_idx("lad_fpl");
      if (lad_col < 0) lad_col = col_idx("ladG0.5");
      has_lad_vox = (lad_col >= 0);

      std::vector<int> liad_cols;
      for (int b = 0; ; ++b) {
        int c = col_idx("liad_" + std::to_string(b));
        if (c < 0) break;
        liad_cols.push_back(c);
      }
      n_iad_bins = (int)liad_cols.size();
      has_liad_vox = (n_iad_bins > 0);

      int bitmap_col = col_idx("subvoxel_bitmap");
      has_bitmap_vox = (bitmap_col >= 0);

      int nb_echos_col = col_idx("nbEchos");

      notify_columns();

      auto safe_val = [](const VoxelData& v, int col) -> double {
        if (col < 0 || col >= (int)v.variables.size()) return 0.0;
        try { return std::stod(v.variables[col]); } catch (...) { return 0.0; }
      };
      for (auto& vd : space.voxels)
      {
        VoxLeafData ld;
        ld.lad = has_lad_vox ? safe_val(vd, lad_col) : 0.0;
        ld.nb_echos = (nb_echos_col >= 0 && nb_echos_col < (int)vd.variables.size())
                        ? (int)safe_val(vd, nb_echos_col) : 0;
        if (has_liad_vox) {
          ld.liad.resize(n_iad_bins);
          for (int b = 0; b < n_iad_bins; ++b)
            ld.liad[b] = safe_val(vd, liad_cols[b]);
        }
        if (has_bitmap_vox) {
          try {
            ld.bitmap = bitmap_col < (int)vd.variables.size()
                          ? static_cast<uint64_t>(std::stoull(vd.variables[bitmap_col])) : 0;
          } catch (...) { ld.bitmap = 0; }
        }
        vox_map[{vd.i, vd.j, vd.k}] = std::move(ld);
      }
    }
    else
    {
      // ---- Extended text format (i j k x y z voxel_state ... num_hits ... subvoxel_bitmap ...) ----
      std::ifstream ifs(vox_file);
      if (!ifs.is_open()) {
        std::cerr << "Error: cannot open rayvoxel file: " << vox_file << std::endl;
        return false;
      }
      std::string header_line;
      if (!std::getline(ifs, header_line)) {
        std::cerr << "Error: empty rayvoxel file: " << vox_file << std::endl;
        return false;
      }
      std::vector<std::string> cols;
      { std::istringstream hss(header_line); std::string t; while (hss >> t) cols.push_back(t); }

      auto fc = [&](const std::string& name) -> int {
        for (int i = 0; i < (int)cols.size(); ++i)
          if (cols[i] == name) return i;
        return -1;
      };
      int i_col = fc("i"), j_col = fc("j"), k_col = fc("k");
      int x_col = fc("x"), y_col = fc("y"), z_col = fc("z");
      int vs_col = fc("voxel_size");
      if (i_col < 0 || j_col < 0 || k_col < 0) {
        std::cerr << "Error: rayvoxel txt file missing i/j/k columns." << std::endl;
        return false;
      }

      int nb_echos_col_txt = fc("num_hits");  // extended txt uses num_hits
      int bitmap_col_txt   = fc("subvoxel_bitmap");

      int lad_col_txt = fc("lad_" + rayvoxel_method);
      if (lad_col_txt < 0) lad_col_txt = fc("lad_fpl");
      if (lad_col_txt < 0) lad_col_txt = fc("lad_g0.5");
      if (lad_col_txt < 0) lad_col_txt = fc("ladG0.5");
      has_lad_vox = (lad_col_txt >= 0);

      std::vector<int> liad_cols_txt;
      for (int b = 0; ; ++b) {
        int c = fc("liad_" + std::to_string(b));
        if (c < 0) break;
        liad_cols_txt.push_back(c);
      }
      n_iad_bins = (int)liad_cols_txt.size();
      has_liad_vox = (n_iad_bins > 0);
      has_bitmap_vox = (bitmap_col_txt >= 0);

      notify_columns();

      bool vox_meta_set = false;
      std::string line;
      while (std::getline(ifs, line))
      {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> vals;
        { std::istringstream ss(line); std::string t; while (ss >> t) vals.push_back(t); }

        int needed = std::max({i_col, j_col, k_col});
        if ((int)vals.size() <= needed) continue;

        try {
          long vi = std::stol(vals[i_col]);
          long vj = std::stol(vals[j_col]);
          long vk = std::stol(vals[k_col]);

          // Infer vox_min and vox_res from x/y/z/voxel_size columns on the first valid row
          if (!vox_meta_set && x_col >= 0 && y_col >= 0 && z_col >= 0 && vs_col >= 0 &&
              std::max({x_col, y_col, z_col, vs_col}) < (int)vals.size())
          {
            double vs = std::stod(vals[vs_col]);
            vox_res = Eigen::Vector3d(vs, vs, vs);
            vox_min = Eigen::Vector3d(std::stod(vals[x_col]) - (vi + 0.5) * vs,
                                     std::stod(vals[y_col]) - (vj + 0.5) * vs,
                                     std::stod(vals[z_col]) - (vk + 0.5) * vs);
            vox_meta_set = true;
          }

          auto safe_txt = [&](int col) -> double {
            if (col < 0 || col >= (int)vals.size()) return 0.0;
            try { return std::stod(vals[col]); } catch (...) { return 0.0; }
          };

          VoxLeafData ld;
          ld.lad = has_lad_vox ? safe_txt(lad_col_txt) : 0.0;
          ld.nb_echos = (nb_echos_col_txt >= 0 && nb_echos_col_txt < (int)vals.size())
                          ? (int)std::stoi(vals[nb_echos_col_txt]) : 0;
          if (has_liad_vox) {
            ld.liad.resize(n_iad_bins);
            for (int b = 0; b < n_iad_bins; ++b) ld.liad[b] = safe_txt(liad_cols_txt[b]);
          }
          if (has_bitmap_vox) {
            try {
              ld.bitmap = (bitmap_col_txt >= 0 && bitmap_col_txt < (int)vals.size())
                            ? static_cast<uint64_t>(std::stoull(vals[bitmap_col_txt])) : 0;
            } catch (...) { ld.bitmap = 0; }
          }
          vox_map[{vi, vj, vk}] = std::move(ld);
        } catch (...) { continue; }
      }

      if (!vox_meta_set) {
        std::cerr << "Error: could not derive voxel origin from rayvoxel txt file (missing x/y/z/voxel_size columns)." << std::endl;
        return false;
      }
    }

    // Aggregate LIAD into 1m DensityGrid cells
    if (has_liad_vox)
    {
      for (auto& [key, ld] : vox_map)
      {
        Eigen::Vector3d world_centre = vox_min + vox_res.cwiseProduct(
          Eigen::Vector3d(key.i + 0.5, key.j + 0.5, key.k + 0.5));
        if (!((world_centre.array() >= grid_bounds.min_bound_.array()).all() &&
              (world_centre.array() < grid_bounds.max_bound_.array()).all()))
          continue;
        int cell = grid.getIndexFromPos(world_centre);
        auto& acc = cell_liad_sum[cell];
        if (acc.empty()) acc.assign(n_iad_bins, 0.0);
        for (int b = 0; b < n_iad_bins; ++b)
          acc[b] += ld.liad[b];
      }
      // Build per-cell CDFs (normalize then prefix-sum); stored in-place in cell_liad_sum
      for (auto& [cell, hist] : cell_liad_sum)
      {
        double sum = 0.0;
        for (double v : hist) sum += v;
        if (sum <= 0.0) { hist.clear(); continue; }
        for (double& v : hist) v /= sum;
        for (int b = 1; b < (int)hist.size(); ++b) hist[b] += hist[b-1];
      }
    }
  }

  // we want to find the few branches that are nearest to each voxel
  // possibly we want there to be no maximum distance... which is weird, but more robust I guess.
  // so the best option is to use knn, and match voxel centres to tree segment centres I guess. This has the advantage
  // that it tends not to align leaves to really thick trunks.
  std::vector<int> tree_ids;
  std::vector<int> segment_ids;
  std::vector<std::vector<int>> neighbour_segments;  // this looks up into the above two structures
  ForestStructure forest;
  std::vector<int> dense_voxel_indices(grid.voxels().size(), -1);
  {  // Tim: this block looks for the closest cylindrical branch segments to each voxel, in order to give the leaves a
     // 'direction' value
    // The reason I use knn (K-nearest neighbour search) is that there is no maximum distance to worry about, and it is
    // fast
    if (!forest.load(trees_file))
    {
      return false;
    }

    size_t num_segments = 0;
    for (auto &tree : forest.trees)
    {
      num_segments += tree.segments().size() - 1;
    }
    size_t num_dense_voxels = 0;
    int i = 0;
    for (auto &vox : grid.voxels())
    {
      if (vox.density() > 0.0)
      {
        dense_voxel_indices[i] = (int)num_dense_voxels;
        num_dense_voxels++;
      }
      i++;
    }

    const int search_size = 12;  // find the twelve nearest branch segments. For larger voxels a larger value here would be helpful
    size_t p_size = num_segments;
    size_t q_size = num_dense_voxels;
    Eigen::MatrixXd points_p(3, p_size);
    i = 0;
    // 1. get branch centre positions
    for (int tree_id = 0; tree_id < (int)forest.trees.size(); tree_id++)
    {
      auto &tree = forest.trees[tree_id];
      for (int segment_id = 0; segment_id < (int)tree.segments().size(); segment_id++)
      {
        auto &segment = tree.segments()[segment_id];
        if (segment.parent_id == -1)
        {
          continue;
        }
        points_p.col(i++) = (segment.tip + tree.segments()[segment.parent_id].tip) / 2.0;
        tree_ids.push_back(tree_id);
        segment_ids.push_back(segment_id);
      }
    }
    // 2. get
    neighbour_segments.resize(grid.voxels().size());
    Eigen::MatrixXd points_q(3, q_size);
    int c = 0;
    for (int k = 0; k < dims[2]; k++)
    {
      for (int j = 0; j < dims[1]; j++)
      {
        for (int i = 0; i < dims[0]; i++)
        {
          int index = grid.getIndex(Eigen::Vector3i(i, j, k));
          double density = grid.voxels()[index].density();
          if (density > 0.0)
          {
            points_q.col(c++) =
              grid_bounds.min_bound_ + vox_width * Eigen::Vector3d((double)i + 0.5, (double)j + 0.5, (double)k + 0.5);
          }
        }
      }
    }
    Nabo::NNSearchD *nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);
    Eigen::MatrixXi indices;
    Eigen::MatrixXd dists2;
    indices.resize(search_size, q_size);
    dists2.resize(search_size, q_size);
    const double max_distance = 2.0;
    nns->knn(points_q, indices, dists2, search_size, kNearestNeighbourEpsilon, 0, max_distance);
    delete nns;

    // Convert these set of nearest neighbours into surfels
    for (int i = 0; i < (int)grid.voxels().size(); i++)
    {
      int id = dense_voxel_indices[i];
      if (id != -1)
      {
        for (int j = 0; j < search_size && indices(j, id) != Nabo::NNSearchD::InvalidIndex; j++)
        {
          neighbour_segments[i].push_back(indices(j, id));
        }
      }
    }
  }


  // the density is now stored in grid.voxels()[grid.getIndex(Eigen::Vector3i )].density().
  struct Leaf
  {
    Eigen::Vector3d centre;
    Eigen::Vector3d direction;
    Eigen::Vector3d origin;
    double grad0;
  };

  std::vector<Leaf> leaves;
  std::vector<double> leaf_counter(grid.voxels().size());
  std::srand(1);
  for (size_t i = 0; i < grid.voxels().size(); i++)
  {
    leaf_counter[i] =
      (double)(std::rand() % 10000) / 10000.0;  // a random start stops regions of low density have 0 leaves
  }

  // Get the leaf angle distribution function based on user input
  auto leafAngleDistribution = getLeafAngleDistribution(distribution);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  std::uniform_real_distribution<> bin_dis(0.0, 1.0);  // for CDF inversion within bin

  // Helper: sample a leaf angle from per-cell LIAD CDF or fall back to analytic distribution
  auto sample_angle = [&](int grid_idx) -> double {
    if (!vox_file.empty() && has_liad_vox)
    {
      auto cit = cell_liad_sum.find(grid_idx);
      if (cit != cell_liad_sum.end() && !cit->second.empty())
      {
        double u = dis(gen);
        const auto &cdf = cit->second;
        int bin = (int)cdf.size() - 1;
        for (int b = 0; b < (int)cdf.size(); ++b)
          if (u <= cdf[b]) { bin = b; break; }
        return (bin + bin_dis(gen)) * 90.0 / n_iad_bins;
      }
    }
    double angle;
    do { angle = dis(gen) * 90.0; } while (dis(gen) > leafAngleDistribution(angle));
    return angle;
  };

  // Helper: find nearest branch to a point; returns false if point is inside a branch cylinder
  auto find_nearest_branch = [&](const Eigen::Vector3d &pos, int grid_idx,
                                 Eigen::Vector3d &closest_out) -> bool {
    double min_dist = 1e10;
    for (auto &ind : neighbour_segments[grid_idx])
    {
      auto &tree = forest.trees[tree_ids[ind]];
      Eigen::Vector3d line_closest;
      Eigen::Vector3d closest = tree.closestPointOnSegment(segment_ids[ind], pos, line_closest);
      double dist = (closest - pos).norm();
      double radius = tree.segments()[segment_ids[ind]].radius;
      if (dist <= radius)
        return false;
      if (dist < min_dist)
      {
        min_dist = dist;
        closest_out = closest;
      }
    }
    return min_dist < 1e10;
  };

  // Helper: build a Leaf given position and closest branch point
  auto build_leaf = [&](const Eigen::Vector3d &centre, const Eigen::Vector3d &closest_branch,
                        int grid_idx) -> Leaf {
    Leaf lf;
    lf.centre = centre;
    lf.origin = closest_branch;

    Eigen::Vector3d branch_dir = (centre - closest_branch).normalized();
    double angle = sample_angle(grid_idx);
    double angle_rad = angle * M_PI / 180.0;

    Eigen::Vector3d rot_axis = branch_dir.cross(Eigen::Vector3d::UnitZ()).normalized();
    if (rot_axis.norm() < 1e-6)
      rot_axis = branch_dir.cross(Eigen::Vector3d::UnitY()).normalized();

    lf.direction = Eigen::AngleAxisd(angle_rad, rot_axis) * branch_dir;

    // Apply droop
    Eigen::Vector3d flat = lf.direction;
    flat[2] = 0.0;
    double dist_h = flat.norm();
    lf.direction[2] -= droop * dist_h * dist_h;
    lf.direction.normalize();

    lf.grad0 = std::tan(angle_rad);
    return lf;
  };

  const bool use_bitmap_path = !vox_file.empty() && has_bitmap_vox;

  if (use_bitmap_path)
  {
    // Detect subvoxel split N from max set bit across all rayvox bitmaps
    int max_bit_seen = 0;
    for (auto &[key, ld] : vox_map)
      for (int b = 63; b >= 0; --b)
        if ((ld.bitmap >> b) & 1) { if (b > max_bit_seen) max_bit_seen = b; break; }
    const int split_n = (max_bit_seen >= 27) ? 4 : (max_bit_seen >= 8) ? 3 : 2;
    const int split_n3 = split_n * split_n * split_n;

    // Build input bitmaps and per-voxel point counts from the input cloud.
    // A rayvox bitmap bit=1 means a ray explored that sub-cell (explored/not-occluded).
    // An input bitmap bit=1 means a foliage point from our input cloud fell in that sub-cell.
    // Three states: hit (input=1), empty (rayvox=1 & input=0), occluded (rayvox=0 & input=0).
    std::unordered_map<VoxKey, uint64_t, VoxKeyHash> input_bitmaps;
    std::unordered_map<VoxKey, int, VoxKeyHash> input_pts_count;
    // Per hit sub-cell: sum of input point positions and count, to compute average for leaf placement.
    std::unordered_map<VoxKey, std::unordered_map<int, std::pair<Eigen::Vector3d, int>>, VoxKeyHash> hit_pos_accum;

    auto build_input = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                           std::vector<double> &, std::vector<ray::RGBA> &colours) {
      for (size_t i = 0; i < ends.size(); i++)
      {
        if (!is_leaf(colours[i].alpha)) continue;
        const Eigen::Vector3d rel = vox_res.array().inverse() * (ends[i] - vox_min).array();
        long vi = (long)std::floor(rel.x());
        long vj = (long)std::floor(rel.y());
        long vk = (long)std::floor(rel.z());
        VoxKey key{vi, vj, vk};
        input_pts_count[key]++;
        double fx = rel.x() - vi, fy = rel.y() - vj, fz = rel.z() - vk;
        int sx = std::min(std::max((int)(fx * split_n), 0), split_n - 1);
        int sy = std::min(std::max((int)(fy * split_n), 0), split_n - 1);
        int sz = std::min(std::max((int)(fz * split_n), 0), split_n - 1);
        int bit = sx + sy * split_n + sz * split_n * split_n;
        input_bitmaps[key] |= (1ULL << bit);
        auto &accum = hit_pos_accum[key][bit];
        accum.first += ends[i];
        accum.second++;
      }
    };
    if (!ray::Cloud::read(cloud_name, build_input))
      return false;

    const double vox_vol = vox_res[0] * vox_res[1] * vox_res[2];

    for (auto &[key, ld] : vox_map)
    {
      uint64_t rayvox_bmp = ld.bitmap;
      auto inp_it = input_bitmaps.find(key);
      uint64_t inp_bmp = (inp_it != input_bitmaps.end()) ? inp_it->second : 0;

      if (rayvox_bmp == 0 && inp_bmp == 0)
        continue;

      // Fraction: how much of the full-lidar echoes are represented in our input cloud.
      // If the voxel had no recorded echoes (occluded / not in vox file), assume fraction=1.
      auto cnt_it = input_pts_count.find(key);
      int n_input_pts = (cnt_it != input_pts_count.end()) ? cnt_it->second : 0;
      double fraction = (ld.nb_echos > 0) ? std::min(1.0, (double)n_input_pts / ld.nb_echos) : 1.0;

      double vox_lad = has_lad_vox ? ld.lad : leafAreaDensity;
      double expected_leaves = vox_lad * fraction * vox_vol / leaf_area;
      if (expected_leaves <= 0.0)
        continue;

      // Partition valid sub-cells into: hit (input=1), occluded (rayvox=0 & input=0).
      // Empty sub-cells (rayvox=1 & input=0) are skipped — no leaf material in explored-but-empty space.
      std::vector<int> hit_bits, occ_bits;
      for (int bit = 0; bit < split_n3; ++bit)
      {
        bool rayvox_bit = (rayvox_bmp >> bit) & 1;
        bool inp_bit    = (inp_bmp    >> bit) & 1;
        if (inp_bit)
          hit_bits.push_back(bit);   // confirmed foliage point
        else if (!rayvox_bit)
          occ_bits.push_back(bit);   // occluded: could contain foliage
        // else: empty (ray passed through, no hit) — skip
      }

      int n_active = (int)hit_bits.size() + (int)occ_bits.size();
      if (n_active == 0)
        continue;

      // Distribute expected leaves evenly across active (non-empty) sub-cells,
      // iterating hits first so that when expected_leaves < n_hit, all leaves land in hits.
      const double leaves_per_active_bit = expected_leaves / n_active;

      for (int pass = 0; pass < 2; ++pass)
      {
        const auto &bits = (pass == 0) ? hit_bits : occ_bits;
        for (int bit : bits)
        {
          int sx = bit % split_n;
          int sy = (bit / split_n) % split_n;
          int sz = bit / (split_n * split_n);

          Eigen::Vector3d sub_centre = vox_min + vox_res.cwiseProduct(
            Eigen::Vector3d(key.i + (sx + 0.5) / split_n,
                            key.j + (sy + 0.5) / split_n,
                            key.k + (sz + 0.5) / split_n));

          if (!((sub_centre.array() >= grid_bounds.min_bound_.array()).all() &&
                (sub_centre.array() < grid_bounds.max_bound_.array()).all()))
            continue;
          int grid_idx = grid.getIndexFromPos(sub_centre);

          double &count = leaf_counter[grid_idx];
          count += leaves_per_active_bit;
          if (count < 1.0)
            continue;
          count -= 1.0;

          if (neighbour_segments[grid_idx].empty())
            continue;

          // For hit sub-cells use the average input scan position (non-grid-aligned, within the sub-cell).
          Eigen::Vector3d leaf_pos = sub_centre;
          if (pass == 0)
          {
            auto hpa_it = hit_pos_accum.find(key);
            if (hpa_it != hit_pos_accum.end())
            {
              auto bp_it = hpa_it->second.find(bit);
              if (bp_it != hpa_it->second.end() && bp_it->second.second > 0)
                leaf_pos = bp_it->second.first / bp_it->second.second;
            }
          }

          Eigen::Vector3d closest;
          if (!find_nearest_branch(leaf_pos, grid_idx, closest))
            continue;

          leaves.push_back(build_leaf(leaf_pos, closest, grid_idx));
        }
      }
    }
  }
  else
  {
    auto add_leaves = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                          std::vector<double> &, std::vector<ray::RGBA> &colours) {
      for (size_t i = 0; i < ends.size(); i++)
      {
        if (!is_leaf(colours[i].alpha))
          continue;
        int index = grid.getIndexFromPos(ends[i]);
        auto &voxel = grid.voxels()[index];

        const double vox_vol = vox_width * vox_width * vox_width;
        double desired_leaf_area = leafAreaDensity * vox_vol;
        if (desired_leaf_area <= 0.0)
          continue;
        double num_leaves_d = desired_leaf_area / leaf_area;
        double num_points = (double)voxel.numHits();
        double &count = leaf_counter[index];
        count += num_leaves_d / num_points;
        if (count < 1.0)
          continue;
        count -= 1.0;

        Eigen::Vector3d closest;
        if (!find_nearest_branch(ends[i], index, closest))
          continue;

        leaves.push_back(build_leaf(ends[i], closest, index));
      }
    };

    if (!ray::Cloud::read(cloud_name, add_leaves))
      return false;
  }

  Mesh leaf_mesh;
  // could read it from file at this point
  auto &leaf_verts = leaf_mesh.vertices();
  auto &leaf_inds =
    leaf_mesh.indexList();  // one per triangle, gives the index into the vertices_ array for each corner
  auto &leaf_uvs = leaf_mesh.uvList();
  Eigen::Vector3d leaf_root(0, 0, 0);

  double leaf_width = std::sqrt(leaf_area / 2.0);
  if (leaf_file.empty())  // generate diamond leaf
  {
    // generate a 2-triangle leaf along y axis
    leaf_verts.push_back(
      Eigen::Vector3d(0, -leaf_width, -leaf_width * leaf_width * droop));  // should leaf droop just vertically?
    leaf_verts.push_back(Eigen::Vector3d(-leaf_width / 2.0, 0, 0));
    leaf_verts.push_back(Eigen::Vector3d(leaf_width / 2.0, 0, 0));
    leaf_verts.push_back(Eigen::Vector3d(0, leaf_width, -leaf_width * leaf_width * droop));
    leaf_inds.push_back(Eigen::Vector3i(0, 2, 1));
    leaf_inds.push_back(Eigen::Vector3i(2, 3, 1));
    leaf_root = leaf_verts[0];
  }
  else if (leaf_file.substr(leaf_file.length() - 4) == ".ply")  // load leaf(s) from .ply mesh
  {
    readPlyMesh(leaf_file, leaf_mesh);
    // work out its total area:
    double total_area = 0.0;
    for (auto &tri : leaf_inds)
    {
      Eigen::Vector3d side = (leaf_verts[tri[1]] - leaf_verts[tri[0]]).cross(leaf_verts[tri[2]] - leaf_verts[tri[0]]);
      total_area += side.norm() / 2.0;
    }
    double scale = std::sqrt(leaf_area / total_area);
    for (auto &vert : leaf_verts)
    {
      vert *= scale;
    }
    leaf_root = leaf_verts[0];
  }
  else if (leaf_file.substr(leaf_file.length() - 4) == ".png")  // generate leaf(s) from image
  {
    stbi_set_flip_vertically_on_load(1);
    int width, height, num_channels;
    unsigned char *image_data = stbi_load(leaf_file.c_str(), &width, &height, &num_channels, 0);
    if (!image_data)
    {
      std::cerr << "Error: cannot load file: " << leaf_file << std::endl;
      return false;
    }
    if (num_channels != 4)
    {
      std::cerr << "Error: png file has no alpha channel and no leaves are rectangles: " << leaf_file
                << ", num channels: " << num_channels << std::endl;
      return false;
    }
    double total_alpha = 0.0;
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        const int index = num_channels * (x + width * y);
        total_alpha += ((double)image_data[index + 3]) / 255.0;
      }
    }
    stbi_image_free(image_data);
    double image_area = (double)width * (double)height;
    total_alpha /= image_area;
    std::cout << "image file: " << leaf_file << " is " << total_alpha * 100.0 << "% opaque (leaf)" << std::endl;

    // generate a 4-triangle rectangle along y axis...
    // rescale the size so actual leaf coverage matches specified area
    double scale = std::sqrt(leaf_area / image_area) / total_alpha;
    double w = 0.5 * (double)width * scale;
    double h = 0.5 * (double)height * scale;
    leaf_verts.push_back(Eigen::Vector3d(-h, -w, -w * w * droop));
    leaf_verts.push_back(Eigen::Vector3d(h, -w, -w * w * droop));
    leaf_verts.push_back(Eigen::Vector3d(-h, 0, 0));
    leaf_verts.push_back(Eigen::Vector3d(h, 0, 0));
    leaf_verts.push_back(Eigen::Vector3d(-h, w, -w * w * droop));
    leaf_verts.push_back(Eigen::Vector3d(h, w, -w * w * droop));
    leaf_inds.push_back(Eigen::Vector3i(0, 1, 2));
    leaf_inds.push_back(Eigen::Vector3i(2, 1, 3));
    leaf_inds.push_back(Eigen::Vector3i(2, 3, 4));
    leaf_inds.push_back(Eigen::Vector3i(4, 3, 5));
    leaf_uvs.push_back(Eigen::Vector3cf(Cmp(0, 0), Cmp(0, 1), Cmp(0.5, 0)));
    leaf_uvs.push_back(Eigen::Vector3cf(Cmp(0.5, 0), Cmp(0, 1), Cmp(0.5, 1)));
    leaf_uvs.push_back(Eigen::Vector3cf(Cmp(0.5, 0), Cmp(0.5, 1), Cmp(1, 0)));
    leaf_uvs.push_back(Eigen::Vector3cf(Cmp(1, 0), Cmp(0.5, 1), Cmp(1, 1)));
    leaf_mesh.textureName() = leaf_file;
    leaf_root = Eigen::Vector3d(0, -w, -w * w * droop);
  }
  else
  {
    std::cerr << "Error: leaf file type unsupported: " << leaf_file << std::endl;
    return false;
  }
  Mesh mesh;
  auto &verts = mesh.vertices();
  auto &inds = mesh.indexList();  // one per triangle, gives the index into the vertices_ array for each corner
  auto &uvs = mesh.uvList();
  mesh.textureName() = leaf_mesh.textureName();

  for (auto &leaf : leaves)
  {
    // 1. convert direction into a transformation matrix...
    Eigen::Matrix3d mat;
    mat.col(1) = leaf.direction;
    mat.col(0) = leaf.direction.cross(Eigen::Vector3d(0, 0, 1)).normalized();
    mat.col(2) = mat.col(0).cross(mat.col(1));

    int num_verts = (int)verts.size();
    for (auto &tri : leaf_inds)
    {
      inds.push_back(tri + Eigen::Vector3i(num_verts, num_verts, num_verts));
    }
    for (auto &uv : leaf_uvs)
    {
      uvs.push_back(uv);  // if UVs are present in the input, they are unchanged
    }
    for (auto &vert : leaf_verts)
    {
      verts.push_back(mat * vert + leaf.centre);
      mesh.colours().push_back(RGBA::leaves());
    }
    num_verts = (int)verts.size();
    if (stalks)
    {
      if (!uvs.empty())
      {
        std::cerr << "Error: multiple textures in one mesh are unsupported, so either turn off stalks or remove "
                     "uvs/texture from leaves"
                  << std::endl;
        return false;
      }
      Eigen::Vector3d start = leaf.origin;
      Eigen::Vector3d leaf_start = mat * leaf_root + leaf.centre;
      Eigen::Vector3d flat = (leaf_start - leaf.origin);
      flat[2] = 0.0;
      double length = flat.norm();
      flat /= length;
      Eigen::Vector3d side(-flat[1], flat[0], flat[2]);
      side *= leaf_width / 16.0;
      const int num_segs = 4;
      for (int i = 0; i < num_segs; i++)
      {
        double x = (double)i / (double)(num_segs - 1);
        x *= length;
        double h = leaf.grad0 * x - droop * x * x;
        Eigen::Vector3d pos = (i == num_segs - 1) ? leaf_start : start + Eigen::Vector3d(0, 0, h) + flat * x;
        verts.push_back(pos - side);
        verts.push_back(pos + side);
        mesh.colours().push_back(RGBA::treetrunk());
        mesh.colours().push_back(RGBA::treetrunk());
        if (i != num_segs - 1)
        {
          int j = 2 * i;
          inds.push_back(Eigen::Vector3i(num_verts, num_verts, num_verts) + Eigen::Vector3i(j, j + 2, j + 1));
          inds.push_back(Eigen::Vector3i(num_verts, num_verts, num_verts) + Eigen::Vector3i(j + 3, j + 1, j + 2));
        }
      }
    }
  }
  writePlyMesh(cloud_stub + "_leaves.ply", mesh);
  return true;
}
}  // namespace ray
