// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYCLOUDWRITER_H
#define RAYLIB_RAYCLOUDWRITER_H

#include "raylib/raylibconfig.h"
#include "raylaz.h"
#include "rayply.h"

namespace ray
{
/// This helper class is for writing a ray cloud to a file, one chunk at a time.
/// The format (PLY or LAS/LAZ) is determined automatically from the file name extension.
/// These chunks can be any size, even 0.
class RAYLIB_EXPORT CloudWriter
{
public:
  ~CloudWriter() { end(); }

  /// Open the file to write to. Format is inferred from the extension (.las/.laz → LAS, else PLY).
  bool begin(const std::string &file_name);

  /// write a set of rays to the file
  bool writeChunk(const class Cloud &chunk);

  /// write a set of rays to the file, direct arguments
  bool writeChunk(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times,
                  std::vector<RGBA> &colours);

  /// finish writing, flush and close the file
  void end();

  /// return the stored file name
  const std::string &fileName() { return file_name_; }

private:
  std::string file_name_;
  bool use_las_ = false;

  // PLY path
  std::ofstream ofs_;
  RayPlyBuffer buffer_;
  bool has_warned_ = false;

  // LAS path
  LasRayCloudWriter *las_writer_ = nullptr;
};

}  // namespace ray

#endif  // RAYLIB_RAYCLOUDWRITER_H
