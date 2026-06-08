// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file provides utility functions for serializing VoxelCoord and Voxel
// data structures to and from binary streams. This is a critical component
// for the out-of-core processing strategy (sharding), as it allows for
// efficient, low-overhead storage of intermediate results on disk.

#ifndef RAYLIB_RAYVOXEL_RAYLASBINARYIO_H
#define RAYLIB_RAYVOXEL_RAYLASBINARYIO_H

#include "raylib/rayvoxel/raylasvoxelise.h"
#include <fstream>

namespace ray
{
  // Helper to write Plain Old Data (POD) types to a binary stream.
  template <typename T>
  inline void writeBinary(std::ofstream& out, const T& value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
  }

  // Helper to read Plain Old Data (POD) types from a binary stream.
  template <typename T>
  inline void readBinary(std::ifstream& in, T& value) {
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
  }

  static constexpr uint32_t kShardMagic   = 0x564F584C; // "VOXL"
  static constexpr uint16_t kShardVersion = 6;          // v6: added num_hits_weighted + path_length_weighted fields

  /// @brief Writes shard-file header. Call once at the start of each shard.
  inline bool writeShardHeader(std::ofstream& out) {
    writeBinary(out, kShardMagic);
    writeBinary(out, kShardVersion);
    return out.good();
  }

  /// @brief Reads and validates the shard-file header.
  inline bool readShardHeader(std::ifstream& in) {
    uint32_t magic; uint16_t version;
    readBinary(in, magic);
    readBinary(in, version);
    return in.good() && magic == kShardMagic && version == kShardVersion;
  }

  /// @brief Serializes a Voxel and its VoxelCoord to a binary output stream.
  inline bool writeVoxelData(std::ofstream& out, const VoxelCoord& coord, const VoxelGrid::Voxel& voxel)
  {
    writeBinary(out, coord.x);
    writeBinary(out, coord.y);
    writeBinary(out, coord.z);
    writeBinary(out, voxel.num_hits);
    writeBinary(out, voxel.num_hits_weighted);
    writeBinary(out, voxel.num_beams_observed);
    writeBinary(out, voxel.num_beams_weighted);
    writeBinary(out, voxel.path_length_observed);
    writeBinary(out, voxel.path_length_weighted);
    writeBinary(out, voxel.num_rays_occluded);
    writeBinary(out, voxel.path_length_occluded);
    writeBinary(out, voxel.sum_of_angles);
    writeBinary(out, voxel.sum_sin_azimuth);
    writeBinary(out, voxel.sum_cos_azimuth);
    writeBinary(out, voxel.sum_of_laser_distances);
    writeBinary(out, voxel.bs_entering);
    writeBinary(out, voxel.bs_intercepted);
    writeBinary(out, voxel.num_unbound_rays);
    writeBinary(out, voxel.path_length_unbound);
    writeBinary(out, voxel.subvoxel_bitmap);
    return out.good();
  }

  /// @brief Deserializes a Voxel and its VoxelCoord from a binary input stream.
  inline bool readVoxelData(std::ifstream& in, VoxelCoord& coord, VoxelGrid::Voxel& voxel)
  {
    readBinary(in, coord.x);
    readBinary(in, coord.y);
    readBinary(in, coord.z);
    if (!in.good()) return false;

    readBinary(in, voxel.num_hits);
    readBinary(in, voxel.num_hits_weighted);
    readBinary(in, voxel.num_beams_observed);
    readBinary(in, voxel.num_beams_weighted);
    readBinary(in, voxel.path_length_observed);
    readBinary(in, voxel.path_length_weighted);
    readBinary(in, voxel.num_rays_occluded);
    readBinary(in, voxel.path_length_occluded);
    readBinary(in, voxel.sum_of_angles);
    readBinary(in, voxel.sum_sin_azimuth);
    readBinary(in, voxel.sum_cos_azimuth);
    readBinary(in, voxel.sum_of_laser_distances);
    readBinary(in, voxel.bs_entering);
    readBinary(in, voxel.bs_intercepted);
    readBinary(in, voxel.num_unbound_rays);
    readBinary(in, voxel.path_length_unbound);
    readBinary(in, voxel.subvoxel_bitmap);
    return in.good();
  }

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASBINARYIO_H
