// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// C++ library for reading and writing AMAPVox .vox files.
// Style is consistent with rayply.h and raylaz.h.

#ifndef RAYLIB_RAYVOXEL_RAYVOX_H
#define RAYLIB_RAYVOXEL_RAYVOX_H

#include "raylib/raylibconfig.h"

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <functional>

namespace ray
{

  /// @brief Formats a 3-vector as a space-delimited string "x y z".
  inline std::string format_vec_string(const Eigen::Vector3d& v) {
    return std::to_string(v.x()) + " " + std::to_string(v.y()) + " " + std::to_string(v.z());
  }
  inline std::string format_vec_string(const Eigen::Vector3i& v) {
    return std::to_string(v.x()) + " " + std::to_string(v.y()) + " " + std::to_string(v.z());
  }

  /// @struct VoxelData
  /// @brief Represents a single voxel entry with its grid indices and variable data.
  struct VoxelData
  {
    long i, j, k;
    std::vector<std::string> variables;
  };

  /// @struct VoxelSpace
  /// @brief Represents the entire content of a .vox file, including header and data.
  struct VoxelSpace
  {
    std::map<std::string, std::string> header;
    std::vector<VoxelData> voxels;

    void clear() {
      header.clear();
      voxels.clear();
    }
  };

  /// @brief Writes a VoxelSpace object to an AMAPVox .vox file.
  /// @param file_name The path to the output .vox file.
  /// @param space The VoxelSpace object containing data to write.
  /// @return True on success, false on failure.
  bool RAYLIB_EXPORT writeVox(const std::string& file_name, const VoxelSpace& space);

  /// @brief Reads an entire AMAPVox .vox file into a VoxelSpace object.
  /// @param file_name The path to the input .vox file.
  /// @param space_out The VoxelSpace object to populate with data from the file.
  /// @return True on success, false on failure.
  bool RAYLIB_EXPORT readVox(const std::string& file_name, VoxelSpace& space_out);

  /// @brief Reads a .vox file in chunks, applying a function to each chunk.
  /// @param file_name The path to the input .vox file.
  /// @param header_out A map to be populated with the file's header information.
  /// @param apply A function called for each chunk of voxels read from the file.
  /// @param chunk_size The number of voxels to read into memory at a time.
  /// @return True on success, false on failure.
  bool RAYLIB_EXPORT readVox(const std::string& file_name,
                           std::map<std::string, std::string>& header_out,
                           std::function<void(const std::vector<VoxelData>&)> apply,
                           size_t chunk_size = 1000000);

  // --- Chunked Writing Functions ---

  /// @brief Writes the header of a .vox file and prepares the stream for data chunks.
  /// @param out The output file stream, opened in text mode.
  /// @param header The map of header properties. Must contain a "colnames" key.
  /// @return True on success, false on failure.
  bool RAYLIB_EXPORT writeVoxChunkStart(std::ofstream& out, const std::map<std::string, std::string>& header);

  /// @brief Writes a chunk of voxel data to an already opened file stream.
  /// @param out The output file stream.
  /// @param voxels A vector of VoxelData to write.
  /// @return True on success, false on failure.
  bool RAYLIB_EXPORT writeVoxChunk(std::ofstream& out, const std::vector<VoxelData>& voxels);

  /// @brief Finalizes a chunked write operation. (No-op for .vox, included for API consistency).
  /// @param out The output file stream.
  void RAYLIB_EXPORT writeVoxChunkEnd(std::ofstream& out);

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYVOX_H
