// Copyright (c) 2024
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#pragma once
#include "raysegmentresult.h"
#include "../raycloud.h"
#include "../raymesh.h"
#include "raytrees.h"
#include <string>

namespace ray {

// Run segmentation (full Trees pipeline); extract per-point labels and seed list.
// cloud.tree_ids and cloud.stem_ids are populated on return.
SegmentResult segment(Cloud &cloud, const Eigen::Vector3d &offset,
                      const Mesh &mesh, const TreesParams &params,
                      bool verbose, const std::string &name_stub);

// Write cloud_seeds.txt from a SeedList.
// Seeds are written in world coords (offset re-added).
bool saveSeeds(const std::string &filename, const SeedList &seeds,
               const Eigen::Vector3d &offset);

// Load cloud_seeds.txt into a SeedList.
// Seeds are returned in local coords (offset subtracted).
bool loadSeeds(const std::string &filename, SeedList &seeds,
               const Eigen::Vector3d &offset);

// Synthesize a SeedList from a cloud's tree_id/stem_id labels (local coords, offset already removed).
// For each unique (tree_id, stem_id) pair with tree_id >= 0, produces one StemSeed whose base is the
// centroid of points in the bottom 5% of that group's height range (min 0.1 m window), with a
// conservative radius prior of 0.5 * params.max_diameter.  Returns seeds sorted by (tree_id, stem_id).
SeedList seedsFromCloud(const Cloud &cloud, const TreesParams &params);

// Run reconstruction from a pre-labeled cloud.
// Derives tree topology directly from cloud.tree_ids / cloud.stem_ids using per-tree Dijkstra.
// Writes name_stub + "_trees.txt" and name_stub + "_trees_mesh.ply".
bool reconstruct(Cloud &cloud, const Eigen::Vector3d &offset,
                 const Mesh &mesh, const TreesParams &params,
                 bool verbose, const std::string &name_stub);

// Backward-compatible composition: segment + reconstruct in one call.
// Writes all outputs: _segmented.{ext}, _trees.txt, _trees_mesh.ply.
void trees(Cloud &cloud, const Eigen::Vector3d &offset, const Mesh &mesh,
           const TreesParams &params, bool verbose,
           const std::string &name_stub, const std::string &input_ext,
           bool save_paths);

} // namespace ray
