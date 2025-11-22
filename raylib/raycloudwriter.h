// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYCLOUDWRITER_H
#define RAYLIB_RAYCLOUDWRITER_H

#include "raylib/raylibconfig.h"
#include "rayply.h"
#include "raylaz.h" // Forward declare LasWriter
#include <vector>

namespace ray
{
class LasWriter; // Forward declaration
class Cloud;     // Forward declaration

/// This helper class is for writing a ray cloud to a file, one chunk at a time
/// These chunks can be any size, even 0
class RAYLIB_EXPORT CloudWriter
{
public:
  // --- START OF MODIFICATION ---
  CloudWriter();
  ~CloudWriter();
  /// Open the file to write to, determining format from extension
  bool begin(const std::string &file_name, const std::string &output_ext = ".ply");
  // --- END OF MODIFICATION ---

  /// write a set of rays to the file
  bool writeChunk(const Cloud &chunk);

  /// write a set of rays to the file, direct arguments
  // --- START OF MODIFICATION ---
  bool writeChunk(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                  const std::vector<RGBA> &colours, const std::vector<uint8_t> &classifications, const std::vector<uint16_t> &branch_ids);
  // --- END OF MODIFICATION ---
  
  /// finish writing, and adjust the vertex count at the start.
  void end();

  /// return the stored file name
  const std::string &fileName() const { return file_name_; }

private:
  /// store the output file stream
  std::ofstream ofs_;
  /// store the file name, in order to provide a clear 'saved' message on end()
  std::string file_name_;
  
  // --- START OF MODIFICATION ---
  /// ray buffer to avoid repeated reallocations, now a dynamic char buffer
  std::vector<char> buffer_;
  /// whether a warning has been issued or not. This prevents multiple warnings.
  bool has_warned_;
  
  // Members for flexible writing
  enum WriterType { PLY, LAZ, NONE };
  WriterType writer_type_;
  LasWriter* las_writer_;
  size_t ply_vertex_byte_size_;
  bool ply_header_written_;
  // --- END OF MODIFICATION ---
};

}  // namespace ray

#endif  // RAYLIB_RAYCLOUDWRITER_H