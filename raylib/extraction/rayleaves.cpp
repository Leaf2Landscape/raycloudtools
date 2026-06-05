// Copyright (c) 2023
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayleaves.h"
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

bool generateLeaves(const std::string &cloud_stub, const std::string &trees_file, const std::string &leaf_file,
                    double leaf_area, double droop, int distribution, double leafAreaDensity, bool stalks,
                    const std::string &vox_file, const std::string &rayvoxel_method)
{
  // For now we assume that woody points have been set as unbounded (alpha=0). e.g. through raycolour foliage or
  // raysplit file distance 0.2 as examples. so firstly we must calculate the foliage density across the whole map.
  std::string cloud_name = cloud_stub + ".ply";
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
  // loaded from a rayvoxel .vox file. When present, these override the scalar --leaf_density
  // and analytic --leaf_angle values on a per-voxel / per-cell basis.
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
  struct VoxLeafData { double lad; std::vector<double> liad; };
  using VoxMap = std::unordered_map<VoxKey, VoxLeafData, VoxKeyHash>;

  VoxMap vox_map;
  Eigen::Vector3d vox_min(0,0,0), vox_res(1,1,1);
  int n_iad_bins = 0;
  // per-DensityGrid-cell LIAD accumulator (keyed by flat 1m-grid index)
  std::unordered_map<int, std::vector<double>> cell_liad_sum;
  bool has_lad_vox = false;
  bool has_liad_vox = false;

  if (!vox_file.empty())
  {
    ray::VoxelSpace space;
    if (!ray::readVox(vox_file, space))
    {
      std::cerr << "Error: cannot read rayvoxel file: " << vox_file << std::endl;
      return false;
    }

    // Parse min_corner and res from header (format: "x y z")
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

    // Parse colnames to find lad and liad column indices (0-based after i j k)
    auto it_col = space.header.find("colnames");
    std::vector<std::string> colnames;
    if (it_col != space.header.end())
    {
      std::istringstream ss(it_col->second);
      std::string tok;
      while (ss >> tok) colnames.push_back(tok);
    }
    // Skip the first 3 (i j k)
    auto col_idx = [&](const std::string& name) -> int {
      for (int c = 3; c < (int)colnames.size(); ++c)
        if (colnames[c] == name) return c - 3; // 0-based into variables[]
      return -1;
    };

    // LAD column: user method -> fpl -> ladG0.5 -> scalar
    int lad_col = col_idx("lad_" + rayvoxel_method);
    if (lad_col < 0) lad_col = col_idx("lad_fpl");
    if (lad_col < 0) lad_col = col_idx("ladG0.5");
    has_lad_vox = (lad_col >= 0);
    if (!has_lad_vox)
      std::cout << "Note: rayvoxel file has no lad column for method '" << rayvoxel_method
                << "'; using --leaf_density as fallback." << std::endl;

    // LIAD columns: liad_0, liad_1, ...
    std::vector<int> liad_cols;
    for (int b = 0; ; ++b)
    {
      int c = col_idx("liad_" + std::to_string(b));
      if (c < 0) break;
      liad_cols.push_back(c);
    }
    n_iad_bins = (int)liad_cols.size();
    has_liad_vox = (n_iad_bins > 0);
    if (!has_liad_vox)
      std::cout << "Note: rayvoxel file has no liad_* columns; using --leaf_angle distribution as fallback." << std::endl;

    // One-time notification if user also provided --leaf_density or --leaf_angle
    // (caller already passed leafAreaDensity and distribution; just note fallback role)
    if (has_lad_vox)
      std::cout << "Note: --rayvoxel active; --leaf_density used only as fallback for uncovered voxels." << std::endl;
    if (has_liad_vox)
      std::cout << "Note: --rayvoxel active; --leaf_angle used only as fallback when liad data is absent." << std::endl;

    // Build VoxMap
    auto safe_val = [](const VoxelData& v, int col) -> double {
      if (col < 0 || col >= (int)v.variables.size()) return 0.0;
      try { return std::stod(v.variables[col]); } catch (...) { return 0.0; }
    };
    for (auto& vd : space.voxels)
    {
      VoxLeafData ld;
      ld.lad = has_lad_vox ? safe_val(vd, lad_col) : 0.0;
      if (has_liad_vox)
      {
        ld.liad.resize(n_iad_bins);
        for (int b = 0; b < n_iad_bins; ++b)
          ld.liad[b] = safe_val(vd, liad_cols[b]);
      }
      vox_map[{vd.i, vd.j, vd.k}] = std::move(ld);
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
        // convert to CDF
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

  auto add_leaves = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends, std::vector<double> &,
                        std::vector<ray::RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); i++)
    {
      if (colours[i].alpha == 0)
        continue;
      int index = grid.getIndexFromPos(ends[i]);
      auto &voxel = grid.voxels()[index];

      double desired_leaf_area;
      if (!vox_file.empty() && has_lad_vox)
      {
        Eigen::Vector3d frac = (ends[i] - vox_min).cwiseQuotient(vox_res);
        VoxKey vk{ static_cast<long>(std::floor(frac.x())),
                   static_cast<long>(std::floor(frac.y())),
                   static_cast<long>(std::floor(frac.z())) };
        auto it = vox_map.find(vk);
        if (it != vox_map.end())
          desired_leaf_area = it->second.lad * vox_res.prod();
        else
          desired_leaf_area = leafAreaDensity * vox_width * vox_width * vox_width; // fallback
      }
      else
      {
        desired_leaf_area = leafAreaDensity * vox_width * vox_width * vox_width;
      }
      if (desired_leaf_area <= 0.0)
        continue;
      double num_leaves_d = desired_leaf_area / leaf_area;
      double num_points = (double)voxel.numHits();
      double &count = leaf_counter[index];
      count += num_leaves_d / num_points;
      bool add_leaf = false;
      if (count >= 1.0)
      {
        add_leaf = true;
        count--;
      }

      if (add_leaf)
      {
        Leaf new_leaf;
        new_leaf.centre = ends[i];

        double min_dist = 1e10;
        Eigen::Vector3d closest_point_on_branch(0, 0, 0);
        for (auto &ind : neighbour_segments[index])
        {
          auto &tree = forest.trees[tree_ids[ind]];
          Eigen::Vector3d line_closest;
          Eigen::Vector3d closest = tree.closestPointOnSegment(segment_ids[ind], ends[i], line_closest);
          double dist = (closest - ends[i]).norm();
          double radius = tree.segments()[segment_ids[ind]].radius;
          if (dist <= radius)
          {
            min_dist = 1e10;
            break;
          }
          if (dist < min_dist)
          {
            min_dist = dist;
            closest_point_on_branch = closest;
          }
        }
        if (min_dist == 1e10)
        {
          continue;
        }        // Calculate leaf direction using the user-specified leaf angle distribution
        Eigen::Vector3d branch_direction = (new_leaf.centre - closest_point_on_branch).normalized();
        
        // Generate a random angle using the distribution
        double angle;
        bool used_liad = false;
        if (!vox_file.empty() && has_liad_vox)
        {
          auto cit = cell_liad_sum.find(index);
          if (cit != cell_liad_sum.end() && !cit->second.empty())
          {
            // Inverse-transform sample from the per-cell CDF
            double u = dis(gen);
            int bin = 0;
            const auto& cdf = cit->second;
            for (int b = 0; b < (int)cdf.size(); ++b)
              if (u <= cdf[b]) { bin = b; break; }
            // Uniform within the bin
            angle = (bin + bin_dis(gen)) * 90.0 / n_iad_bins;
            used_liad = true;
          }
        }
        if (!used_liad)
        {
          do {
            angle = dis(gen) * 90.0; // Random angle between 0 and 90 degrees
          } while (dis(gen) > leafAngleDistribution(angle));
        }

        // Convert angle to radians
        double angle_rad = angle * M_PI / 180.0;

        // Create a rotation axis perpendicular to the branch direction
        Eigen::Vector3d rotation_axis = branch_direction.cross(Eigen::Vector3d::UnitZ()).normalized();
        if (rotation_axis.norm() < 1e-6) {
          rotation_axis = branch_direction.cross(Eigen::Vector3d::UnitY()).normalized();
        }

        // Create rotation matrix
        Eigen::AngleAxisd rotation(angle_rad, rotation_axis);
        
        // Apply rotation to branch direction to get leaf direction
        new_leaf.direction = rotation * branch_direction;

        // Apply droop
        Eigen::Vector3d flat = new_leaf.direction;
        flat[2] = 0.0;
        double dist = flat.norm();
        new_leaf.direction[2] -= droop * dist * dist;
        new_leaf.direction.normalize();

        new_leaf.origin = closest_point_on_branch;
        new_leaf.grad0 = std::tan(angle_rad);

        leaves.push_back(new_leaf);
      }  
    }
  };

  if (!ray::Cloud::read(cloud_name, add_leaves))
    return false;

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
