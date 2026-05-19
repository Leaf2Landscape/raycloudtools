// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYLAZ_H
#define RAYLIB_RAYLAZ_H

#include "raylib/raylibconfig.h"
#include "rayutils.h"

#if RAYLIB_WITH_LAS
#include <laszip/laszip_api.h>
#endif  // RAYLIB_WITH_LAS


namespace ray
{
/// Read a laz or las file, into the fields passed by reference.
bool RAYLIB_EXPORT readLas(std::string file_name, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times,
                           std::vector<RGBA> &colours, double max_intensity,
                           Eigen::Vector3d *offset_to_remove = nullptr);

/// Chunk-based version of readLas. This calls @c apply for every @c chunk_size points loaded.
/// When @c tree_ids_out is non-null and the file contains a tree_id extra attribute, tree IDs are appended to it.
/// When @c passthrough_out is non-null, 8 standard LAS fields + original sensor extra bytes are appended per point.
/// @c orig_extra_size_out receives the per-point byte count of original sensor extra bytes (may be 0).
/// @c extra_bytes_vlr_out receives the raw EXTRA_BYTES VLR payload for the original sensor attributes.
bool RAYLIB_EXPORT readLas(const std::string &file_name,
                           std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                              std::vector<double> &times, std::vector<RGBA> &colours)>
                             apply,
                           size_t &num_bounded, double max_intensity, Eigen::Vector3d *offset_to_remove,
                           size_t chunk_size = 1000000, std::vector<int32_t> *tree_ids_out = nullptr,
                           std::vector<uint8_t> *passthrough_out = nullptr,
                           uint16_t *orig_extra_size_out = nullptr,
                           std::vector<uint8_t> *extra_bytes_vlr_out = nullptr);


/// Read only the EXTRA_BYTES VLR from a las/laz file header without reading any point data.
/// On return, @c orig_extra_size_out is the total per-point byte count of non-raycloud extra attributes,
/// and @c extra_bytes_vlr_out contains the raw 192-byte VLR records for those attributes.
/// Returns false if the file cannot be opened or LAS support is not compiled in.
bool RAYLIB_EXPORT readLasExtraBytesVlr(const std::string &file_name, uint16_t &orig_extra_size_out,
                                        std::vector<uint8_t> &extra_bytes_vlr_out);

/// Write to a laz or las file. The intensity is the only part that is extracted from the @c colours argument.
bool RAYLIB_EXPORT writeLas(std::string file_name, const std::vector<Eigen::Vector3d> &points,
                            const std::vector<double> &times, const std::vector<RGBA> &colours);

/// Write a ray cloud to a las/laz file. Ray starts are stored as float32 extra bytes.
/// RGBA colour is fully preserved: RGB in the LAS colour fields, alpha in intensity.
/// When @c tree_ids is non-empty an additional int32 "tree_id" extra byte attribute is written per point.
/// When @c extra_bytes_vlr is non-empty, the original sensor extra-byte attributes are registered and
/// copied from @c passthrough (bytes [8..] per point at stride @c extra_bytes_size).
bool RAYLIB_EXPORT writeLasRayCloud(const std::string &file_name, const std::vector<Eigen::Vector3d> &starts,
                                    const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                                    const std::vector<RGBA> &colours,
                                    const std::vector<int32_t> &tree_ids = {},
                                    const std::vector<uint8_t> &passthrough = {},
                                    const std::vector<uint8_t> &extra_bytes_vlr = {});

/// Class for chunked writing of las/laz files.
class RAYLIB_EXPORT LasWriter
{
public:
  /// construct the class with a file name, which is stored
  LasWriter(const std::string &file_name);
  /// the destructor
  ~LasWriter();
  /// write a chunk of points to the file, described by the vector arguments
  bool writeChunk(const std::vector<Eigen::Vector3d> &points, const std::vector<double> &times,
                  const std::vector<RGBA> &colours);

private:
  const std::string &file_name_;
#if RAYLIB_WITH_LAS
  laszip_POINTER writer_handle_;
  laszip_header_struct *header_;
  laszip_point_struct *point_;
  uint64_t points_written_;
#endif  // RAYLIB_WITH_LAS
};

/// Class for chunked writing of las/laz ray cloud files.
/// Ray starts are stored as three float32 extra bytes (start - end offset).
/// RGBA is fully preserved: RGB in LAS colour fields, alpha in intensity.
/// When @c with_tree_id is true, a fourth int32 "tree_id" extra attribute is added.
/// When @c extra_bytes_vlr is non-empty, the original sensor extra-byte attributes are
/// registered from its 192-byte EXTRA_BYTES VLR records and written per point from passthrough[8+].
class RAYLIB_EXPORT LasRayCloudWriter
{
public:
  explicit LasRayCloudWriter(const std::string &file_name, bool with_tree_id = false,
                             const std::vector<uint8_t> &extra_bytes_vlr = {});
  ~LasRayCloudWriter();
  bool writeChunk(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends,
                  const std::vector<double> &times, const std::vector<RGBA> &colours,
                  const std::vector<int32_t> &tree_ids = {},
                  const std::vector<uint8_t> &passthrough = {});
  unsigned long pointCount() const { return points_written_; }

private:
  std::string file_name_;
  uint64_t points_written_ = 0;
  bool with_tree_id_ = false;
  uint16_t orig_extra_size_ = 0;   ///< per-point original sensor extra bytes
  uint16_t passthrough_stride_ = 8; ///< 8 + orig_extra_size_
#if RAYLIB_WITH_LAS
  laszip_POINTER writer_handle_;
  laszip_point_struct *point_;
#endif  // RAYLIB_WITH_LAS
};

}  // namespace ray

#endif  // RAYLIB_RAYLAZ_H
