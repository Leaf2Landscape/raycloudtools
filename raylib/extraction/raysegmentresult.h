// Copyright (c) 2024
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#pragma once
#include <vector>
#include <cstdint>
#include <string>
#include <Eigen/Core>

namespace ray {

struct StemSeed {
    int32_t         tree_id;          // inventory tree ID; 0-based; -1 is not valid
    int32_t         stem_id;          // stem within tree; 0 = single-stem default
    Eigen::Vector3d base;             // trunk base, local coords (after removeStartPos)
    double          radius;           // trunk radius at girth window; > 0; 0 if unknown (PLY mask form)
    double          tree_height;      // max(z) - min(z) over all stem points; > 0
    std::string     source_filename;  // PLY form: "{tree_id}_{stem_id}.ply"; LAS form: ""
};

using SeedList = std::vector<StemSeed>;

// Returned by ray::segment(); passed to ray::reconstruct()
struct SegmentResult {
    SeedList              seeds;     // one entry per unique (tree_id, stem_id)
    std::vector<int32_t>  tree_ids;  // cloud.ends-parallel; -1 = unassigned
    std::vector<int32_t>  stem_ids;  // cloud.ends-parallel; -1 = unassigned
};

} // namespace ray
