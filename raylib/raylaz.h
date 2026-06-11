// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYLAZ_H
#define RAYLIB_RAYLAZ_H

#include "raylib/raylibconfig.h"
#include "rayutils.h"
#include <string>
#include <unordered_map>
#include <vector>

#if RAYLIB_WITH_LAS
#include <laszip/laszip_api.h>
#endif  // RAYLIB_WITH_LAS


namespace ray
{

/// A single extra-byte attribute declared in the LAS EXTRA_BYTES VLR.
struct RAYLIB_EXPORT LasExtraField
{
  char     name[33] = {};     ///< null-terminated field name (32-byte VLR name field)
  uint16_t offset   = 0;     ///< byte offset within the extra_bytes block per point
  uint8_t  dtype    = 0;     ///< LAS data_type (1-10)
  uint16_t size     = 0;     ///< per-point byte count derived from dtype
  bool     is_own   = false; ///< true iff this field belongs to the raycloudtools schema
  uint8_t  vlr_record[192] = {}; ///< raw 192-byte EXTRA_BYTES VLR record for this field
};

/// Complete summary of a LAS/LAZ file header, including a generic table of all extra-byte
/// attributes discovered from the EXTRA_BYTES VLR. Callers query fields by name via has() or
/// field() without hard-coding known attribute names in their own logic.
struct RAYLIB_EXPORT LasHeader
{
  uint8_t  point_format        = 0;
  bool     is_raycloud         = false;  ///< written by raycloudtools
  bool     has_rgb             = false;  ///< point format carries native RGB fields
  bool     is_compressed       = false;  ///< .laz file
  uint64_t point_count         = 0;
  double   scale[3]            = { 1.0, 1.0, 1.0 };
  double   offset_xyz[3]       = { 0.0, 0.0, 0.0 };
  uint16_t point_record_length = 0;
  uint32_t point_data_offset   = 0;

  std::vector<LasExtraField>                    extras;   ///< all extra fields, in VLR order
  std::unordered_map<std::string, std::size_t>  by_name;  ///< name -> index into extras

  /// Returns true iff an extra field with this name exists.
  bool has(const std::string &name) const { return by_name.count(name) > 0; }

  /// Returns a pointer to the field descriptor, or nullptr if absent.
  const LasExtraField *field(const std::string &name) const
  {
    auto it = by_name.find(name);
    return (it != by_name.end()) ? &extras[it->second] : nullptr;
  }

  /// Total per-point byte count of non-raycloudtools extra attributes.
  uint16_t sensorExtraSize() const
  {
    uint16_t sz = 0;
    for (const auto &f : extras) if (!f.is_own) sz += f.size;
    return sz;
  }

  /// Raw 192-byte VLR records for non-raycloudtools extra attributes (for passthrough on write).
  std::vector<uint8_t> sensorExtraVlr() const
  {
    std::vector<uint8_t> out;
    for (const auto &f : extras)
      if (!f.is_own)
        out.insert(out.end(), f.vlr_record, f.vlr_record + 192);
    return out;
  }
};

/// Read the header (VLRs, point count, geometry, extra-byte table) of a LAS/LAZ file without
/// reading any point data. Returns false if the file cannot be opened or LAS support is not
/// compiled in. On success, @c header_out contains a complete field table for the file.
bool RAYLIB_EXPORT readLasHeader(const std::string &file_name, LasHeader &header_out);


/// Read a laz or las file, into the fields passed by reference.
bool RAYLIB_EXPORT readLas(std::string file_name, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times,
                           std::vector<RGBA> &colours, double max_intensity,
                           Eigen::Vector3d *offset_to_remove = nullptr);

/// Chunk-based version of readLas. This calls @c apply for every @c chunk_size points loaded.
/// When @c tree_ids_out is non-null and the file contains a tree_id extra attribute, tree IDs are appended to it.
/// When @c stem_ids_out is non-null and the file contains a stem_id extra attribute, stem IDs are appended to it.
/// When @c passthrough_out is non-null, 10 standard LAS field bytes + original sensor extra bytes are appended per point.
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
                           std::vector<uint8_t> *extra_bytes_vlr_out = nullptr,
                           std::vector<int32_t> *stem_ids_out = nullptr,
                           std::vector<int32_t> *beam_ids_out = nullptr);


/// Read only the EXTRA_BYTES VLR from a las/laz file header without reading any point data.
/// On return, @c orig_extra_size_out is the total per-point byte count of non-raycloud extra attributes,
/// and @c extra_bytes_vlr_out contains the raw 192-byte VLR records for those attributes.
/// When @c has_bound_out is non-null, it is set to true iff the file declares a "bound" extra attribute
/// (absent in older files written before the bound field existed).
/// When @c has_rgb_out is non-null, it is set to true iff the file uses point format 7+ (native RGB fields).
/// When @c has_tree_id_out / @c has_stem_id_out is non-null, set to true iff the file has those attributes.
/// Returns false if the file cannot be opened or LAS support is not compiled in.
bool RAYLIB_EXPORT readLasExtraBytesVlr(const std::string &file_name, uint16_t &orig_extra_size_out,
                                        std::vector<uint8_t> &extra_bytes_vlr_out,
                                        bool *has_bound_out = nullptr, bool *has_rgb_out = nullptr,
                                        bool *has_tree_id_out = nullptr, bool *has_stem_id_out = nullptr);

/// Write to a laz or las file. The intensity is the only part that is extracted from the @c colours argument.
bool RAYLIB_EXPORT writeLas(std::string file_name, const std::vector<Eigen::Vector3d> &points,
                            const std::vector<double> &times, const std::vector<RGBA> &colours);

/// Write a ray cloud to a las/laz file. Ray starts are stored as float32 extra bytes.
/// RGBA colour is fully preserved: RGB in the LAS colour fields, alpha in intensity.
/// When @c tree_ids is non-empty an additional int32 "tree_id" extra byte attribute is written per point.
/// When @c stem_ids is non-empty an additional int32 "stem_id" extra byte attribute is written per point.
/// When @c extra_bytes_vlr is non-empty, the original sensor extra-byte attributes are registered and
/// copied from @c passthrough (bytes [8..] per point at stride @c extra_bytes_size).
bool RAYLIB_EXPORT writeLasRayCloud(const std::string &file_name, const std::vector<Eigen::Vector3d> &starts,
                                    const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                                    const std::vector<RGBA> &colours,
                                    const std::vector<int32_t> &tree_ids = {},
                                    const std::vector<int32_t> &stem_ids = {},
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
/// When @c with_rgb is true, point format 7 is used (native RGB fields); otherwise format 6 (no RGB).
/// When @c with_tree_id is true, a fourth int32 "tree_id" extra attribute is added.
/// When @c with_stem_id is true, a fifth int32 "stem_id" extra attribute is added (requires with_tree_id).
/// When @c with_beam_id is true, an int32 "beam_id" extra attribute is added (per-pulse beam ID).
/// When @c extra_bytes_vlr is non-empty, the original sensor extra-byte attributes are
/// registered from its 192-byte EXTRA_BYTES VLR records and written per point from passthrough[8+].
///
/// Authoritative extra-byte write order (see writeChunk):
///   sx, sy, sz, [tree_id], [stem_id], [beam_id], alpha, bound, <sensor extras>
/// Bracketed attributes are present only when the corresponding with_*_id flag is set.
class RAYLIB_EXPORT LasRayCloudWriter
{
public:
  explicit LasRayCloudWriter(const std::string &file_name, bool with_tree_id = false,
                             bool with_stem_id = false, bool with_beam_id = false,
                             const std::vector<uint8_t> &extra_bytes_vlr = {},
                             bool with_rgb = false);
  ~LasRayCloudWriter();
  bool writeChunk(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends,
                  const std::vector<double> &times, const std::vector<RGBA> &colours,
                  const std::vector<int32_t> &tree_ids = {},
                  const std::vector<int32_t> &stem_ids = {},
                  const std::vector<uint8_t> &passthrough = {},
                  const std::vector<int32_t> &beam_ids = {});
  unsigned long pointCount() const { return points_written_; }

private:
  std::string file_name_;
  uint64_t points_written_ = 0;
  bool with_tree_id_ = false;
  bool with_stem_id_ = false;
  bool with_beam_id_ = false;
  bool with_rgb_ = false;
  Eigen::Vector3d bbox_min_{ std::numeric_limits<double>::max(),  std::numeric_limits<double>::max(),  std::numeric_limits<double>::max()  };
  Eigen::Vector3d bbox_max_{ std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest() };
  uint16_t orig_extra_size_ = 0;   ///< per-point original sensor extra bytes
  uint16_t passthrough_stride_ = 10; ///< 10 + orig_extra_size_
#if RAYLIB_WITH_LAS
  laszip_POINTER writer_handle_;
  laszip_point_struct *point_;
#endif  // RAYLIB_WITH_LAS
};

}  // namespace ray

#endif  // RAYLIB_RAYLAZ_H
