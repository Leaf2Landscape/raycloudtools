// Copyright (c) 2024
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#include "raytreeslib.h"

#include "../rayforeststructure.h"
#include "../raymesh.h"
#include "../rayply.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <vector>

namespace ray {

// --------------------------------------------------------------------------
// ray::segment()
// --------------------------------------------------------------------------
SegmentResult segment(Cloud &cloud, const Eigen::Vector3d &offset,
                      const Mesh &mesh, const TreesParams &params,
                      bool verbose, const std::string &name_stub)
{
  Trees trees_obj(cloud, offset, mesh, params, verbose);
  // cloud.tree_ids is now populated by segmentCloud (phase 15)

  SegmentResult result;
  result.seeds    = seedsFromCloud(cloud, params);  // derived from labeled cloud, local coords
  result.tree_ids = cloud.tree_ids;

  cloud.stem_ids.assign(cloud.ends.size(), -1);
  for (size_t i = 0; i < cloud.tree_ids.size(); ++i)
    if (cloud.tree_ids[i] != -1)
      cloud.stem_ids[i] = 0;
  result.stem_ids = cloud.stem_ids;

  return result;
}

// --------------------------------------------------------------------------
// ray::saveSeeds()
// --------------------------------------------------------------------------
bool saveSeeds(const std::string &filename, const SeedList &seeds,
               const Eigen::Vector3d &offset)
{
  std::ofstream ofs(filename);
  if (!ofs.is_open())
  {
    std::cerr << "raytreeslib: could not write seeds file: " << filename << std::endl;
    return false;
  }
  ofs << std::setprecision(6) << std::fixed;
  ofs << "# tree seeds file:\n";
  ofs << "x,y,z,radius,tree_id,stem_id,tree_height\n";
  for (const auto &s : seeds)
  {
    // Write in world coords (base + offset)
    ofs << s.base.x() + offset.x() << ", "
        << s.base.y() + offset.y() << ", "
        << s.base.z() + offset.z() << ", "
        << s.radius << ", "
        << s.tree_id << ", "
        << s.stem_id << ", "
        << s.tree_height << "\n";
  }
  return true;
}

// --------------------------------------------------------------------------
// ray::loadSeeds()
// --------------------------------------------------------------------------
bool loadSeeds(const std::string &filename, SeedList &seeds,
               const Eigen::Vector3d &offset)
{
  std::ifstream ifs(filename);
  if (!ifs.is_open())
  {
    std::cerr << "raytreeslib: could not open seeds file: " << filename << std::endl;
    return false;
  }
  seeds.clear();
  std::string line;
  bool header_passed = false;
  while (std::getline(ifs, line))
  {
    if (line.empty() || line[0] == '#') continue;
    // Skip the column header line
    if (!header_passed)
    {
      header_passed = true;
      continue;
    }
    // Replace commas with spaces for sscanf-like parsing
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream ss(line);
    double x, y, z, radius, tree_height;
    int32_t tree_id, stem_id;
    if (!(ss >> x >> y >> z >> radius >> tree_id >> stem_id >> tree_height))
      continue;
    StemSeed s;
    s.base         = Eigen::Vector3d(x - offset.x(), y - offset.y(), z - offset.z());
    s.radius       = radius;
    s.tree_id      = tree_id;
    s.stem_id      = stem_id;
    s.tree_height  = tree_height;
    s.source_filename = "";
    seeds.push_back(s);
  }
  return !seeds.empty();
}

// --------------------------------------------------------------------------
// ray::seedsFromCloud()
// --------------------------------------------------------------------------
SeedList seedsFromCloud(const Cloud &cloud, const TreesParams &params)
{
  using Key = std::pair<int32_t, int32_t>;

  struct GroupData {
    double          min_z     = std::numeric_limits<double>::max();
    double          max_z     = std::numeric_limits<double>::lowest();
    Eigen::Vector3d min_z_pos = Eigen::Vector3d::Zero();
    Eigen::Vector3d base_sum  = Eigen::Vector3d::Zero();
    int             base_count = 0;
  };

  std::map<Key, GroupData> groups;
  const bool has_stem = !cloud.stem_ids.empty();

  // Pass 1: height range and lowest-point position per group.
  for (size_t i = 0; i < cloud.ends.size(); ++i)
  {
    if (!cloud.rayBounded(i)) continue;
    if (i >= cloud.tree_ids.size()) continue;
    int32_t tid = cloud.tree_ids[i];
    if (tid < 0) continue;
    int32_t sid = (has_stem && i < cloud.stem_ids.size()) ? cloud.stem_ids[i] : 0;
    if (sid < 0) sid = 0;

    auto &g = groups[{tid, sid}];
    double z = cloud.ends[i].z();
    if (z < g.min_z) { g.min_z = z; g.min_z_pos = cloud.ends[i]; }
    g.max_z = std::max(g.max_z, z);
  }

  // Pass 2: accumulate base candidates (bottom 5 % of height range, min 0.1 m window).
  for (size_t i = 0; i < cloud.ends.size(); ++i)
  {
    if (!cloud.rayBounded(i)) continue;
    if (i >= cloud.tree_ids.size()) continue;
    int32_t tid = cloud.tree_ids[i];
    if (tid < 0) continue;
    int32_t sid = (has_stem && i < cloud.stem_ids.size()) ? cloud.stem_ids[i] : 0;
    if (sid < 0) sid = 0;

    auto &g = groups[{tid, sid}];
    double threshold = g.min_z + std::max(0.05 * (g.max_z - g.min_z), 0.1);
    if (cloud.ends[i].z() <= threshold)
    {
      g.base_sum   += cloud.ends[i];
      g.base_count += 1;
    }
  }

  SeedList seeds;
  seeds.reserve(groups.size());
  for (auto &kv : groups)
  {
    auto &g = kv.second;
    StemSeed s;
    s.tree_id         = kv.first.first;
    s.stem_id         = kv.first.second;
    s.tree_height     = std::max(g.max_z - g.min_z, 1.0);
    s.radius          = 0.5 * params.max_diameter;
    s.source_filename = "";
    s.base = (g.base_count > 0) ? (g.base_sum / static_cast<double>(g.base_count))
                                 : g.min_z_pos;
    seeds.push_back(s);
  }

  std::sort(seeds.begin(), seeds.end(), [](const StemSeed &a, const StemSeed &b) {
    return a.tree_id < b.tree_id || (a.tree_id == b.tree_id && a.stem_id < b.stem_id);
  });
  return seeds;
}

// --------------------------------------------------------------------------
// ray::reconstruct()
// --------------------------------------------------------------------------
bool reconstruct(Cloud &cloud, const Eigen::Vector3d &offset,
                 const Mesh &mesh, const TreesParams &params,
                 bool verbose, const std::string &name_stub)
{
  // For PLY input: tree identity is color-encoded; decode into cloud.tree_ids.
  // Black (0,0,0) decodes to -1 (unassigned) via convertColourToInt.
  if (cloud.tree_ids.empty() && !cloud.colours.empty())
  {
    cloud.tree_ids.resize(cloud.ends.size(), -1);
    for (size_t i = 0; i < cloud.ends.size(); ++i)
    {
      const int id = convertColourToInt(cloud.colours[i]);
      cloud.tree_ids[i] = static_cast<int32_t>(id); // -1 for black = unassigned
    }
  }

  // Use the pre-labeled constructor (per-tree Dijkstra; no seeds file required).
  Trees trees_obj(cloud, offset, mesh, params, verbose, PreLabeledTag{});

  // Build id_map from the cloud's tree_id/stem_id labels.
  auto id_map = trees_obj.buildLabelIdMap();

  trees_obj.save(name_stub + "_trees.txt", offset, verbose, id_map);

  ray::ForestStructure forest;
  if (!forest.load(name_stub + "_trees.txt")) return false;
  ray::Mesh tree_mesh;
  forest.generateSmoothMesh(tree_mesh, -1, 1, 1, 1);
  ray::writePlyMesh(name_stub + "_trees_mesh.ply", tree_mesh, true);
  return true;
}

// --------------------------------------------------------------------------
// ray::trees() — backward-compatible composition
// --------------------------------------------------------------------------
void trees(Cloud &cloud, const Eigen::Vector3d &offset, const Mesh &mesh,
           const TreesParams &params, bool verbose,
           const std::string &name_stub, const std::string &input_ext,
           bool save_paths)
{
  Trees trees_obj(cloud, offset, mesh, params, verbose);
  trees_obj.save(name_stub + "_trees.txt", offset, verbose);
  if (save_paths)
    trees_obj.saveShortestPaths(name_stub + "_shortest_paths.ply", offset);
  cloud.translate(offset);
  cloud.save(name_stub + "_segmented." + input_ext);
  ray::ForestStructure forest;
  if (!forest.load(name_stub + "_trees.txt")) return;
  ray::Mesh tree_mesh;
  forest.generateSmoothMesh(tree_mesh, -1, 1, 1, 1);
  ray::writePlyMesh(name_stub + "_trees_mesh.ply", tree_mesh, true);
}

} // namespace ray
