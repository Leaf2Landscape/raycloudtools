// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYCLOUDREADER_H
#define RAYLIB_RAYCLOUDREADER_H

#include "raylib/raylibconfig.h"
#include "raylaz.h"
#include "rayply.h"

namespace ray
{
/// Unified chunked reader for ray clouds in PLY, LAS, or LAZ format.
/// For LAS/LAZ files the LAS header is read once on begin(); the same header
/// is then available via header() without a second file open. PLY files have
/// no LAS header (header() returns a default-constructed LasHeader).
class RAYLIB_EXPORT CloudReader
{
public:
  ~CloudReader() { end(); }

  /// Open the file. For LAS/LAZ this reads the header immediately so that
  /// header() is valid before read() is called. Returns false on error.
  bool begin(const std::string &file_name);

  /// LAS/LAZ header, populated after begin(). Always default-constructed for PLY.
  const LasHeader &header() const { return header_; }

  /// True iff the file is LAS or LAZ.
  bool isLas() const { return is_las_; }

  /// Read all points, calling apply() once per chunk.
  /// Non-null output pointers receive the corresponding per-point arrays
  /// (tree_ids, stem_ids, beam_ids, passthrough). For PLY these are always
  /// empty. chunk_size == 0 uses computeReadChunkSize().
  bool read(std::function<void(std::vector<Eigen::Vector3d> &starts,
                               std::vector<Eigen::Vector3d> &ends,
                               std::vector<double> &times,
                               std::vector<RGBA> &colours)>
              apply,
            size_t &num_bounded,
            double max_intensity = 1.0,
            Eigen::Vector3d *offset_to_remove = nullptr,
            size_t chunk_size = 0,
            std::vector<int32_t> *tree_ids_out = nullptr,
            std::vector<uint8_t> *passthrough_out = nullptr,
            std::vector<int32_t> *stem_ids_out = nullptr,
            std::vector<int32_t> *beam_ids_out = nullptr);

  void end() { file_name_.clear(); is_las_ = false; header_ = LasHeader{}; }

  const std::string &fileName() const { return file_name_; }

private:
  std::string file_name_;
  bool is_las_ = false;
  LasHeader header_;
};

}  // namespace ray

#endif  // RAYLIB_RAYCLOUDREADER_H
