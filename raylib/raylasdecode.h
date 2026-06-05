// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYLASDECODE_H
#define RAYLIB_RAYLASDECODE_H

#include "raylib/raylibconfig.h"
#include "rayutils.h"

#if RAYLIB_WITH_LAS
#include <cstring>
#include <laszip/laszip_api.h>

namespace ray
{
// LAS EXTRA_BYTES data_type -> per-point byte size (types 0 and >10 are skipped).
static const uint16_t kDecodeExtraTypeSize[11] = { 0, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };

/// Per-file parameters needed to decode a single LAS/LAZ point record into ray-cloud fields.
/// Populated once from the LAS header + EXTRA_BYTES VLR, then reused for every point so that the
/// per-point decode is a pure function of (point record, context). This makes the decode safe to
/// run in parallel and guarantees bit-identical output across all read paths.
struct DecodeContext
{
  uint8_t format = 0;        ///< point_data_format
  bool using_colour = false; ///< whether the format carries RGB
  bool is_raycloud = false;  ///< whether this file was written by raycloudtools
  double max_intensity = 1.0;

  uint16_t local_skip_size = 0;   ///< bytes of our own extra attributes before original sensor data
  uint16_t local_orig_extra = 0;  ///< bytes of original sensor data per point
  uint16_t tree_id_offset = 0;    uint8_t tree_id_dtype = 0;  ///< 0 = absent
  uint16_t stem_id_offset = 0;    uint8_t stem_id_dtype = 0;  ///< 0 = absent
  uint16_t beam_id_offset = 0;    uint8_t beam_id_dtype = 0;  ///< 0 = absent
  uint16_t sx_offset = 0;       ///< byte offset of "sx" within extra_bytes (default: 0)
  uint16_t sy_offset = 4;       ///< byte offset of "sy" within extra_bytes (default: 4)
  uint16_t sz_offset = 8;       ///< byte offset of "sz" within extra_bytes (default: 8)
  uint16_t alpha_offset = 12;     ///< default: sx+sy+sz only; overwritten when "alpha" VLR found
  int32_t bound_offset = -1;      ///< -1 = absent (old file); byte offset of the "bound" extra attribute

  // Fields below are only required by the raw-record decode (mmap / laz-perf paths). They mirror
  // the LAS header scaling and record geometry so a fixed-layout point record can be decoded
  // without the laszip reader.
  double scale[3] = { 1.0, 1.0, 1.0 };  ///< x/y/z scale factors from the LAS header
  double offset[3] = { 0.0, 0.0, 0.0 }; ///< x/y/z offsets from the LAS header
  uint16_t point_record_length = 0;     ///< total bytes per point record on disk
  uint16_t extra_bytes_offset = 0;      ///< byte offset of the extra-bytes block within a record
  uint16_t extra_bytes_total = 0;       ///< total extra-byte count per record (ours + original)
};

/// Per-format byte offset, within a fixed LAS point record, of the start of the extra-bytes block.
/// Returns 0 for unrecognised formats (caller should fall back to the laszip path).
inline uint16_t lasBaseRecordSize(uint8_t format)
{
  switch (format)
  {
    case 0: return 20;
    case 1: return 28;
    case 2: return 26;
    case 3: return 34;
    case 4: return 57;
    case 5: return 63;
    case 6: return 30;
    case 7: return 36;
    case 8: return 38;
    case 9: return 59;
    case 10: return 67;
    default: return 0;
  }
}

/// Populate a laszip_point_struct from a fixed-layout LAS point record @c rec, mirroring exactly
/// what laszip_read_point would yield. Only the fields consumed by decodePointRecord are filled.
/// @c extra_scratch must have capacity for ctx.extra_bytes_total; it is pointed to by pt.extra_bytes.
/// Returns false for formats whose fixed layout is not supported (caller falls back to laszip).
inline bool fillPointFromRecord(const uint8_t *rec, const DecodeContext &ctx, laszip_point_struct &pt,
                                uint8_t *extra_scratch)
{
  const uint8_t format = ctx.format;
  if (lasBaseRecordSize(format) == 0)
    return false;

  // X/Y/Z int32 little-endian at bytes 0..11.
  int32_t xi, yi, zi;
  memcpy(&xi, rec + 0, 4);
  memcpy(&yi, rec + 4, 4);
  memcpy(&zi, rec + 8, 4);
  pt.X = xi;
  pt.Y = yi;
  pt.Z = zi;

  // Intensity uint16 LE at bytes 12..13.
  memcpy(&pt.intensity, rec + 12, 2);

  if (format >= 6)
  {
    // LAS 1.4 point formats 6-10.
    const uint8_t b14 = rec[14];  // return_number (bits 0-3) | number_of_returns (bits 4-7)
    const uint8_t b15 = rec[15];  // class_flags (0-3) | scanner_channel (4-5) | scan_dir (6) | edge (7)
    pt.extended_return_number = b14 & 0x0F;
    pt.extended_number_of_returns = (b14 >> 4) & 0x0F;
    pt.extended_classification_flags = b15 & 0x0F;
    pt.extended_scanner_channel = (b15 >> 4) & 0x03;
    pt.scan_direction_flag = (b15 >> 6) & 0x01;
    pt.edge_of_flight_line = (b15 >> 7) & 0x01;
    pt.extended_classification = rec[16];
    pt.user_data = rec[17];
    int16_t scan_angle;
    memcpy(&scan_angle, rec + 18, 2);
    pt.extended_scan_angle = scan_angle;
    memcpy(&pt.point_source_ID, rec + 20, 2);
    // GPS time is mandatory for formats 6-10, at byte 22.
    memcpy(&pt.gps_time, rec + 22, 8);
    // RGB present in formats 7, 8, 10 immediately after GPS time (byte 30).
    if (format == 7 || format == 8 || format == 10)
    {
      memcpy(&pt.rgb[0], rec + 30, 2);
      memcpy(&pt.rgb[1], rec + 32, 2);
      memcpy(&pt.rgb[2], rec + 34, 2);
    }
  }
  else
  {
    // Legacy LAS 1.0-1.3 point formats 0-5.
    const uint8_t b14 = rec[14];  // return_number (0-2) | number_of_returns (3-5) | scan_dir (6) | edge (7)
    const uint8_t b15 = rec[15];  // classification (0-4) | synthetic (5) | keypoint (6) | withheld (7)
    pt.return_number = b14 & 0x07;
    pt.number_of_returns = (b14 >> 3) & 0x07;
    pt.scan_direction_flag = (b14 >> 6) & 0x01;
    pt.edge_of_flight_line = (b14 >> 7) & 0x01;
    pt.classification = b15 & 0x1F;
    pt.synthetic_flag = (b15 >> 5) & 0x01;
    pt.keypoint_flag = (b15 >> 6) & 0x01;
    pt.withheld_flag = (b15 >> 7) & 0x01;
    pt.scan_angle_rank = static_cast<laszip_I8>(rec[16]);
    pt.user_data = rec[17];
    memcpy(&pt.point_source_ID, rec + 18, 2);
    // GPS time in formats 1, 3, 4, 5 at byte 20.
    if (format == 1 || format == 3 || format == 4 || format == 5)
      memcpy(&pt.gps_time, rec + 20, 8);
    // RGB in format 2 at byte 20; in formats 3, 5 at byte 28.
    if (format == 2)
    {
      memcpy(&pt.rgb[0], rec + 20, 2);
      memcpy(&pt.rgb[1], rec + 22, 2);
      memcpy(&pt.rgb[2], rec + 24, 2);
    }
    else if (format == 3 || format == 5)
    {
      memcpy(&pt.rgb[0], rec + 28, 2);
      memcpy(&pt.rgb[1], rec + 30, 2);
      memcpy(&pt.rgb[2], rec + 32, 2);
    }
  }

  // Extra bytes: copy the per-record extra block into the caller-provided scratch buffer.
  pt.num_extra_bytes = ctx.extra_bytes_total;
  pt.extra_bytes = extra_scratch;
  if (ctx.extra_bytes_total > 0)
    memcpy(extra_scratch, rec + ctx.extra_bytes_offset, ctx.extra_bytes_total);
  return true;
}

// Read a tree/stem ID extra-byte field of the given LAS data_type and normalise its
// per-type "unassigned" sentinel to int32_t -1. uint types use 0 as unassigned;
// int types use -1. IDs exceeding INT32_MAX are narrowed (known limitation).
inline int32_t decodeLasIdField(const uint8_t *extra, uint16_t off, uint8_t dtype)
{
  switch (dtype)
  {
    case 1: { uint8_t  v; memcpy(&v, extra + off, 1); return v  == 0u ? -1 : static_cast<int32_t>(v); }
    case 3: { uint16_t v; memcpy(&v, extra + off, 2); return v  == 0u ? -1 : static_cast<int32_t>(v); }
    case 5: { uint32_t v; memcpy(&v, extra + off, 4); return v  == 0u ? -1 : static_cast<int32_t>(v); }
    case 7: { uint64_t v; memcpy(&v, extra + off, 8); return v  == 0u ? -1 : static_cast<int32_t>(v); }
    case 2: { int8_t   v; memcpy(&v, extra + off, 1); return v  == -1 ? -1 : static_cast<int32_t>(v); }
    case 4: { int16_t  v; memcpy(&v, extra + off, 2); return v  == -1 ? -1 : static_cast<int32_t>(v); }
    case 6: { int32_t  v; memcpy(&v, extra + off, 4); return v; }  // raycloudtools native; identity
    case 8: { int64_t  v; memcpy(&v, extra + off, 8); return v  == -1 ? -1 : static_cast<int32_t>(v); }
    default: return -1;  // absent (0), float (9/10), or unknown
  }
}

/// Decode a single LAS/LAZ point (as a laszip_point_struct) into the ray-cloud field vectors.
/// This is a verbatim extraction of the per-point decode that previously lived inline in the
/// readLas chunk loop; behaviour is identical. The caller supplies @c position (the decoded
/// coordinates) so this routine never touches the laszip reader handle.
///
/// @c num_bounded is incremented for every point whose decoded alpha/intensity is non-zero.
/// Optional outputs (tree/stem/beam IDs, passthrough) are appended only when their pointer is
/// non-null and the matching attribute is present, exactly as in the original loop.
inline void decodePointRecord(const laszip_point_struct *point, const DecodeContext &ctx,
                              const Eigen::Vector3d &position, std::vector<Eigen::Vector3d> &starts,
                              std::vector<Eigen::Vector3d> &ends, std::vector<double> &times,
                              std::vector<RGBA> &colours, std::vector<uint8_t> &intensities,
                              size_t &num_bounded, std::vector<int32_t> *tree_ids_out,
                              std::vector<uint8_t> *passthrough_out, std::vector<int32_t> *stem_ids_out,
                              std::vector<int32_t> *beam_ids_out)
{
  const uint8_t format = ctx.format;

  ends.push_back(position);

  if (ctx.is_raycloud && point->num_extra_bytes >= 12)
  {
    // Reconstruct ray start from the stored (start - end) float32 offset.
    float sx, sy, sz;
    memcpy(&sx, point->extra_bytes + ctx.sx_offset, 4);
    memcpy(&sy, point->extra_bytes + ctx.sy_offset, 4);
    memcpy(&sz, point->extra_bytes + ctx.sz_offset, 4);
    starts.push_back({ position[0] + sx, position[1] + sy, position[2] + sz });
    if (tree_ids_out && ctx.tree_id_dtype != 0 &&
        point->num_extra_bytes >= ctx.tree_id_offset + kDecodeExtraTypeSize[ctx.tree_id_dtype])
      tree_ids_out->push_back(decodeLasIdField(point->extra_bytes, ctx.tree_id_offset, ctx.tree_id_dtype));
    if (stem_ids_out && ctx.stem_id_dtype != 0 &&
        point->num_extra_bytes >= ctx.stem_id_offset + kDecodeExtraTypeSize[ctx.stem_id_dtype])
      stem_ids_out->push_back(decodeLasIdField(point->extra_bytes, ctx.stem_id_offset, ctx.stem_id_dtype));
    else if (stem_ids_out && ctx.tree_id_dtype != 0 && ctx.stem_id_dtype == 0)
      stem_ids_out->push_back(0);
    if (beam_ids_out && ctx.beam_id_dtype != 0 &&
        point->num_extra_bytes >= ctx.beam_id_offset + kDecodeExtraTypeSize[ctx.beam_id_dtype])
      beam_ids_out->push_back(decodeLasIdField(point->extra_bytes, ctx.beam_id_offset, ctx.beam_id_dtype));
  }
  else
  {
    starts.push_back(position);
  }

  // Pack 10 bytes of LAS fields per point into passthrough, followed by sensor extras.
  // Layout: [0] ext_return[0:3]|ext_num_returns[4:7]
  //         [1] class_flags[0:3]|scanner_chan[4:5]|scan_dir[6]|edge[7]
  //         [2] extended_classification
  //         [3] user_data
  //         [4-5] extended_scan_angle (int16 LE, 0.006 deg units)
  //         [6-7] point_source_ID (uint16 LE)
  //         [8-9] original intensity (uint16 LE)
  if (passthrough_out)
  {
    uint8_t b0, b1, b2;
    int16_t ext_angle;
    if (format >= 6)
    {
      b0 = static_cast<uint8_t>(point->extended_return_number & 0x0F) |
           static_cast<uint8_t>((point->extended_number_of_returns & 0x0F) << 4);
      b1 = static_cast<uint8_t>(point->extended_classification_flags & 0x0F) |
           static_cast<uint8_t>((point->extended_scanner_channel & 0x03) << 4) |
           static_cast<uint8_t>((point->scan_direction_flag & 0x1) << 6) |
           static_cast<uint8_t>((point->edge_of_flight_line & 0x1) << 7);
      b2 = point->extended_classification;
      ext_angle = point->extended_scan_angle;
    }
    else
    {
      // Convert LAS 1.2 legacy fields to LAS 1.4 extended layout.
      b0 = static_cast<uint8_t>(point->return_number & 0x0F) |
           static_cast<uint8_t>((point->number_of_returns & 0x0F) << 4);
      b1 = static_cast<uint8_t>(point->synthetic_flag & 0x1) |
           static_cast<uint8_t>((point->keypoint_flag & 0x1) << 1) |
           static_cast<uint8_t>((point->withheld_flag & 0x1) << 2) |
           static_cast<uint8_t>((point->scan_direction_flag & 0x1) << 6) |
           static_cast<uint8_t>((point->edge_of_flight_line & 0x1) << 7);
      b2 = static_cast<uint8_t>(point->classification & 0x1F);
      // scan_angle_rank is integer degrees; extended_scan_angle is 0.006 deg units
      ext_angle = static_cast<int16_t>(static_cast<int>(point->scan_angle_rank) * 167);
    }
    passthrough_out->push_back(b0);
    passthrough_out->push_back(b1);
    passthrough_out->push_back(b2);
    passthrough_out->push_back(point->user_data);
    passthrough_out->push_back(static_cast<uint8_t>(static_cast<uint16_t>(ext_angle) & 0xFFu));
    passthrough_out->push_back(static_cast<uint8_t>(static_cast<uint16_t>(ext_angle) >> 8));
    passthrough_out->push_back(static_cast<uint8_t>(point->point_source_ID & 0xFFu));
    passthrough_out->push_back(static_cast<uint8_t>(point->point_source_ID >> 8));
    // Original intensity (uint16 LE) — preserved so the output keeps the sensor value.
    passthrough_out->push_back(static_cast<uint8_t>(point->intensity & 0xFFu));
    passthrough_out->push_back(static_cast<uint8_t>(point->intensity >> 8));
    // Append original sensor extra bytes (after skipping our raycloud-owned attributes).
    if (ctx.local_orig_extra > 0 && point->num_extra_bytes >= ctx.local_skip_size + ctx.local_orig_extra)
      passthrough_out->insert(passthrough_out->end(),
                              point->extra_bytes + ctx.local_skip_size,
                              point->extra_bytes + ctx.local_skip_size + ctx.local_orig_extra);
    else if (ctx.local_orig_extra > 0)
      passthrough_out->insert(passthrough_out->end(), ctx.local_orig_extra, 0);
  }

  if (ctx.using_colour)
  {
    RGBA col;
    // RGB stored as uint8 * 257 in the 16-bit field; low byte recovers original value.
    col.red = static_cast<uint8_t>(point->rgb[0]);
    col.green = static_cast<uint8_t>(point->rgb[1]);
    col.blue = static_cast<uint8_t>(point->rgb[2]);
    colours.push_back(col);
  }
  times.push_back(point->gps_time);

  uint8_t intensity;
  if (ctx.is_raycloud)
  {
    // Alpha is stored in extra_bytes at position 12 (or 16 with tree_id, or 20 with tree_id+stem_id).
    // Prefer extra_bytes so the intensity field is free to carry the original sensor value.
    const uint16_t alpha_pos = ctx.alpha_offset;
    intensity = (point->num_extra_bytes > alpha_pos)
                  ? point->extra_bytes[alpha_pos]
                  : static_cast<uint8_t>(point->intensity);  // fallback for old files
    // When the explicit bound field is present (new-format files), treat it as authoritative.
    // Old files (bound_offset == -1) fall through to the alpha > 0 sentinel unchanged.
    if (ctx.bound_offset >= 0 &&
        point->num_extra_bytes > static_cast<uint16_t>(ctx.bound_offset))
    {
      const uint8_t b = point->extra_bytes[ctx.bound_offset];
      if (b == 0 && intensity > 0)
        intensity = 0;  // file says unbound; suppress stray alpha
      else if (b != 0 && intensity == 0)
        intensity = 1;  // file says bound but alpha was zero; mark as bounded
    }
  }
  else
  {
    const double normalised = (ctx.max_intensity > 0) ? (255.0 * point->intensity) / ctx.max_intensity : 255.0;
    intensity = static_cast<uint8_t>(std::min(normalised, 255.0));
    // Ensure any non-zero raw intensity maps to at least alpha=1 (bounded ray).
    if (intensity == 0 && point->intensity > 0)
      intensity = 1;
  }
  if (intensity > 0)
    num_bounded++;
  intensities.push_back(intensity);
}

/// Outputs of the indexed decode, all pre-sized so points can be written in any order / in parallel.
/// A *_active flag is false when the corresponding attribute is absent from the file, mirroring the
/// conditional appends in decodePointRecord (which depend only on per-file, not per-point, state).
struct IndexedDecodeBuffers
{
  Eigen::Vector3d *starts = nullptr;   ///< [n]
  Eigen::Vector3d *ends = nullptr;     ///< [n]
  double *times = nullptr;             ///< [n]
  RGBA *colours = nullptr;             ///< [n] (only written when ctx.using_colour)
  uint8_t *intensities = nullptr;      ///< [n]
  uint8_t *passthrough = nullptr;      ///< [n * passthrough_stride], may be null
  uint16_t passthrough_stride = 0;     ///< 10 + ctx.local_orig_extra
  int32_t *tree_ids = nullptr;         ///< [n], may be null
  int32_t *stem_ids = nullptr;         ///< [n], may be null
  int32_t *beam_ids = nullptr;         ///< [n], may be null
  bool tree_active = false;            ///< whether tree_id is decoded for this file
  bool stem_active = false;            ///< whether stem_id (or its 0 fallback) is decoded
  bool beam_active = false;            ///< whether beam_id is decoded for this file
};

/// Index-addressed twin of decodePointRecord: identical field math, but writes to fixed slots so
/// the decode can run in parallel. @c num_bounded_out is set to 1 when the point is bounded so the
/// caller can sum/reduce it. The conditional ID appends become unconditional indexed writes guarded
/// by the per-file *_active flags, preserving the exact values produced by the sequential path.
inline void decodePointRecordIndexed(const laszip_point_struct *point, const DecodeContext &ctx,
                                     const Eigen::Vector3d &position, size_t i,
                                     IndexedDecodeBuffers &buf, uint8_t &num_bounded_out)
{
  const uint8_t format = ctx.format;
  num_bounded_out = 0;

  buf.ends[i] = position;

  if (ctx.is_raycloud && point->num_extra_bytes >= 12)
  {
    float sx, sy, sz;
    memcpy(&sx, point->extra_bytes + ctx.sx_offset, 4);
    memcpy(&sy, point->extra_bytes + ctx.sy_offset, 4);
    memcpy(&sz, point->extra_bytes + ctx.sz_offset, 4);
    buf.starts[i] = Eigen::Vector3d(position[0] + sx, position[1] + sy, position[2] + sz);
    if (buf.tree_ids && buf.tree_active)
      buf.tree_ids[i] = decodeLasIdField(point->extra_bytes, ctx.tree_id_offset, ctx.tree_id_dtype);
    if (buf.stem_ids && buf.stem_active)
    {
      if (ctx.stem_id_dtype != 0 &&
          point->num_extra_bytes >= ctx.stem_id_offset + kDecodeExtraTypeSize[ctx.stem_id_dtype])
        buf.stem_ids[i] = decodeLasIdField(point->extra_bytes, ctx.stem_id_offset, ctx.stem_id_dtype);
      else
        buf.stem_ids[i] = 0;  // tree_id present but no stem_id attribute: matches sequential fallback
    }
    if (buf.beam_ids && buf.beam_active)
      buf.beam_ids[i] = decodeLasIdField(point->extra_bytes, ctx.beam_id_offset, ctx.beam_id_dtype);
  }
  else
  {
    buf.starts[i] = position;
  }

  if (buf.passthrough)
  {
    uint8_t *out = buf.passthrough + i * buf.passthrough_stride;
    uint8_t b0, b1, b2;
    int16_t ext_angle;
    if (format >= 6)
    {
      b0 = static_cast<uint8_t>(point->extended_return_number & 0x0F) |
           static_cast<uint8_t>((point->extended_number_of_returns & 0x0F) << 4);
      b1 = static_cast<uint8_t>(point->extended_classification_flags & 0x0F) |
           static_cast<uint8_t>((point->extended_scanner_channel & 0x03) << 4) |
           static_cast<uint8_t>((point->scan_direction_flag & 0x1) << 6) |
           static_cast<uint8_t>((point->edge_of_flight_line & 0x1) << 7);
      b2 = point->extended_classification;
      ext_angle = point->extended_scan_angle;
    }
    else
    {
      b0 = static_cast<uint8_t>(point->return_number & 0x0F) |
           static_cast<uint8_t>((point->number_of_returns & 0x0F) << 4);
      b1 = static_cast<uint8_t>(point->synthetic_flag & 0x1) |
           static_cast<uint8_t>((point->keypoint_flag & 0x1) << 1) |
           static_cast<uint8_t>((point->withheld_flag & 0x1) << 2) |
           static_cast<uint8_t>((point->scan_direction_flag & 0x1) << 6) |
           static_cast<uint8_t>((point->edge_of_flight_line & 0x1) << 7);
      b2 = static_cast<uint8_t>(point->classification & 0x1F);
      ext_angle = static_cast<int16_t>(static_cast<int>(point->scan_angle_rank) * 167);
    }
    out[0] = b0;
    out[1] = b1;
    out[2] = b2;
    out[3] = point->user_data;
    out[4] = static_cast<uint8_t>(static_cast<uint16_t>(ext_angle) & 0xFFu);
    out[5] = static_cast<uint8_t>(static_cast<uint16_t>(ext_angle) >> 8);
    out[6] = static_cast<uint8_t>(point->point_source_ID & 0xFFu);
    out[7] = static_cast<uint8_t>(point->point_source_ID >> 8);
    out[8] = static_cast<uint8_t>(point->intensity & 0xFFu);
    out[9] = static_cast<uint8_t>(point->intensity >> 8);
    if (ctx.local_orig_extra > 0 && point->num_extra_bytes >= ctx.local_skip_size + ctx.local_orig_extra)
      memcpy(out + 10, point->extra_bytes + ctx.local_skip_size, ctx.local_orig_extra);
    else if (ctx.local_orig_extra > 0)
      memset(out + 10, 0, ctx.local_orig_extra);
  }

  if (ctx.using_colour)
  {
    RGBA col;
    col.red = static_cast<uint8_t>(point->rgb[0]);
    col.green = static_cast<uint8_t>(point->rgb[1]);
    col.blue = static_cast<uint8_t>(point->rgb[2]);
    buf.colours[i] = col;
  }
  buf.times[i] = point->gps_time;

  uint8_t intensity;
  if (ctx.is_raycloud)
  {
    const uint16_t alpha_pos = ctx.alpha_offset;
    intensity = (point->num_extra_bytes > alpha_pos)
                  ? point->extra_bytes[alpha_pos]
                  : static_cast<uint8_t>(point->intensity);
    if (ctx.bound_offset >= 0 &&
        point->num_extra_bytes > static_cast<uint16_t>(ctx.bound_offset))
    {
      const uint8_t b = point->extra_bytes[ctx.bound_offset];
      if (b == 0 && intensity > 0)
        intensity = 0;
      else if (b != 0 && intensity == 0)
        intensity = 1;
    }
  }
  else
  {
    const double normalised = (ctx.max_intensity > 0) ? (255.0 * point->intensity) / ctx.max_intensity : 255.0;
    intensity = static_cast<uint8_t>(std::min(normalised, 255.0));
    if (intensity == 0 && point->intensity > 0)
      intensity = 1;
  }
  if (intensity > 0)
    num_bounded_out = 1;
  buf.intensities[i] = intensity;
}

}  // namespace ray
#endif  // RAYLIB_WITH_LAS

#endif  // RAYLIB_RAYLASDECODE_H
