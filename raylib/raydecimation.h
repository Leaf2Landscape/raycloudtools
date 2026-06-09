// Copyright (c) 2024
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYDECIMATION_H
#define RAYLIB_RAYDECIMATION_H

#include <iostream>
#include <limits>
#include "raycloud.h"
#include "rayutils.h"

namespace ray
{
/// @brief subsample to 1 point per @c vox_width wide voxel in metres
/// This is a spatially even subsampling, but also emphasises outlier as a side-effect
bool RAYLIB_EXPORT decimateSpatial(const std::string &file_name, double vox_width);

/// @brief subsample to every @c num_rays rays
/// This is an unbiased subsampling, but will be over-sampled in stationary areas as a side-effect
/// Note that while this is called temporal decimation, it decimates evenly in file order, which isn't
/// necessarily temporal order. Though it typically is stored that way on single scans.
bool RAYLIB_EXPORT decimateTemporal(const std::string &file_name, int num_rays);

/// @brief subsample to @c num_rays rays (temporally decimated) for each @c vox_width wide voxel
/// This allows a more even distribution of points while maintaining details better than pure spatial decimation
bool RAYLIB_EXPORT decimateSpatioTemporal(const std::string &file_name, double vox_width, int num_rays);

/// @brief Maintains a maximum number of rays intersecting each voxel. This has some ambiguity, but is a useful routine
/// as it maintains the integrity of the full ray cloud including free space, so is better for combine operations
/// By contrast, standard spatial decimation removes free space whenever the end points coincide
bool RAYLIB_EXPORT decimateRaysSpatial(const std::string &file_name, double vox_width);


/// @brief decimate to no more than 1 point per voxel of width @c radius_per_length x ray length.
/// This is used when error is proportional to ray length, prioritising closer measurements and leaving distant areas sparse
bool RAYLIB_EXPORT decimateAngular(const std::string &file_name, double radius_per_length);


/// Field a tiebreak comparison reads when picking the best point per voxel for @c deduplicateVoxel.
enum class TiebreakKind { Reflectance, Range, Time, ExtraByte };

/// One resolved tiebreak criterion: which field to read, the sort direction, and (for ExtraByte)
/// how to decode the value out of the per-point passthrough slice.
struct RAYLIB_EXPORT ResolvedTiebreaker
{
  TiebreakKind kind;
  bool ascending;        ///< true = prefer lower value, false = prefer higher value
  uint16_t byte_offset = 0;  ///< within the sensor-extras slice (after the 10-byte fixed prefix)
  uint8_t byte_size = 0;
  uint8_t dtype = 0;     ///< LAS extra-byte type code (same table as raycombine.cpp's parseSensorAttrs)
};

/// @brief best-wins voxel dedup post-step over an already-combined ray cloud file.
/// Two passes over @c file_name: pass 1 finds the global winner point per @c vox_width voxel cell
/// by comparing @c spec fields in priority order; pass 2 emits only the winners. The result is
/// written to a temp file and renamed over @c file_name in place. Handles PLY input (no passthrough)
/// gracefully. Returns false on read/write failure.
bool RAYLIB_EXPORT deduplicateVoxel(const std::string &file_name, double vox_width,
                                    const std::vector<ResolvedTiebreaker> &spec);


struct Subsampler
{
  inline bool operator()(const Eigen::Vector3i &p, const Eigen::Vector3i &/*target*/, double /*in_length*/, double /*out_length*/, double /*max_length*/)
  {
    if (voxel_set.insert(p).second)
    {
      subsample.push_back(index);
      return true;
    }
    return false;
  }
  std::vector<int64_t> subsample;
  std::set<Eigen::Vector3i, ray::Vector3iLess> voxel_set;
  int index;
};
}  // namespace ray

#endif  // RAYLIB_RAYDECIMATION_H
