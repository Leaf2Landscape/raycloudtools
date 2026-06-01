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

  /// @brief Serializes a Voxel and its VoxelCoord to a binary output stream.
  inline bool writeVoxelData(std::ofstream& out, const VoxelCoord& coord, const VoxelGrid::Voxel& voxel)
  {
    writeBinary(out, coord.x);
    writeBinary(out, coord.y);
    writeBinary(out, coord.z);
    writeBinary(out, voxel.num_hits);
    writeBinary(out, voxel.num_rays_observed);
    writeBinary(out, voxel.path_length_observed);
    writeBinary(out, voxel.num_rays_occluded);
    writeBinary(out, voxel.path_length_occluded);
    writeBinary(out, voxel.is_filled);
    writeBinary(out, voxel.sum_of_angles);
    writeBinary(out, voxel.sum_of_laser_distances);
    writeBinary(out, voxel.bs_entering);
    writeBinary(out, voxel.bs_intercepted);
    writeBinary(out, voxel.subvoxel_bitmap);

    // Sparse-on-disk: count nonzero entries, then (code, value) pairs
    uint16_t nonzero = 0;
    for (int c = 0; c < 256; ++c) if (voxel.classification_hits[c] != 0.0f) ++nonzero;
    writeBinary(out, nonzero);
    for (int c = 0; c < 256; ++c) {
      if (voxel.classification_hits[c] != 0.0f) {
        U8 code = static_cast<U8>(c);
        writeBinary(out, code);
        writeBinary(out, voxel.classification_hits[c]);
      }
    }
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
    readBinary(in, voxel.num_rays_observed);
    readBinary(in, voxel.path_length_observed);
    readBinary(in, voxel.num_rays_occluded);
    readBinary(in, voxel.path_length_occluded);
    readBinary(in, voxel.is_filled);
    readBinary(in, voxel.sum_of_angles);
    readBinary(in, voxel.sum_of_laser_distances);
    readBinary(in, voxel.bs_entering);
    readBinary(in, voxel.bs_intercepted);
    readBinary(in, voxel.subvoxel_bitmap);

    voxel.classification_hits.fill(0.0f);
    uint16_t nonzero;
    readBinary(in, nonzero);
    for (uint16_t i = 0; i < nonzero; ++i) {
      U8 code; float value;
      readBinary(in, code);
      readBinary(in, value);
      if (!in.good()) return false;
      voxel.classification_hits[code] = value;
    }
    return in.good();
  }

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASBINARYIO_H
