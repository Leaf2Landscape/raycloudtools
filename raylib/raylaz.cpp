// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylaz.h"
#include <algorithm>
#include <fstream>
#include <limits>
#include <memory>
#include <thread>
#include "raylib/raylasdecode.h"
#include "raylib/rayprogress.h"
#include "raylib/rayprogressthread.h"
#include "raylib/rayvoxel/raylasthreadsafequeue.h"
#include "rayunused.h"

#if RAYLIB_WITH_LAS
#include <laszip/laszip_api.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif  // _WIN32
#if RAYLIB_WITH_LAZPERF
#include <lazperf/readers.hpp>
#endif  // RAYLIB_WITH_LAZPERF
#endif  // RAYLIB_WITH_LAS

namespace ray
{
bool readLas(const std::string &file_name,
             std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                std::vector<double> &times, std::vector<RGBA> &colours)>
               apply,
             size_t &num_bounded, double max_intensity, Eigen::Vector3d *offset_to_remove, size_t chunk_size,
             std::vector<int32_t> *tree_ids_out, std::vector<uint8_t> *passthrough_out,
             uint16_t *orig_extra_size_out, std::vector<uint8_t> *extra_bytes_vlr_out,
             std::vector<int32_t> *stem_ids_out, std::vector<int32_t> *beam_ids_out)
{
#if RAYLIB_WITH_LAS
  std::cout << "readLas: filename: " << file_name << std::endl;

  laszip_POINTER reader;
  if (laszip_create(&reader))
  {
    std::cerr << "readLas: failed to create LASzip reader" << std::endl;
    return false;
  }

  laszip_BOOL is_compressed;
  if (laszip_open_reader(reader, file_name.c_str(), &is_compressed))
  {
    laszip_CHAR *error;
    laszip_get_error(reader, &error);
    std::cerr << "readLas: failed to open stream: " << error << std::endl;
    laszip_destroy(reader);
    return false;
  }

  laszip_header_struct *header;
  laszip_get_header_pointer(reader, &header);

  Eigen::Vector3d offset(header->x_offset, header->y_offset, header->z_offset);
  if (offset_to_remove)
  {
    *offset_to_remove = offset;
    std::cout << "offset to remove: " << offset.transpose() << std::endl;
  }

  // LAS 1.4 uses a 64-bit point count; legacy uses the 32-bit field
  const size_t number_of_points =
    (header->version_minor >= 4 && header->extended_number_of_point_records > 0)
      ? static_cast<size_t>(header->extended_number_of_point_records)
      : static_cast<size_t>(header->number_of_point_records);

  const uint8_t format = header->point_data_format;
  // Formats 1,3,4,5 have GPS time in LAS 1.0-1.3; formats 6-10 always have GPS time (LAS 1.4)
  const bool using_time = (format == 1 || format == 3 || format == 4 || format == 5 || format >= 6);
  // Formats 2,3,5 have RGB in LAS 1.0-1.3; formats 7,8,10 have RGB in LAS 1.4
  const bool using_colour = (format == 2 || format == 3 || format == 5 || format == 7 || format == 8 || format == 10);

  if (!using_time)
  {
    std::cerr << "No timestamps found on laz file, these are required" << std::endl;
    laszip_close_reader(reader);
    laszip_destroy(reader);
    return false;
  }

  laszip_point_struct *point;
  laszip_get_point_pointer(reader, &point);

  // Detect rayclouds written by raycloudtools: look for the custom VLR marker.
  bool is_raycloud = false;
  for (laszip_U32 v = 0; v < header->number_of_variable_length_records; v++)
  {
    if (strncmp(header->vlrs[v].user_id, "raycloudtools", 16) == 0 && header->vlrs[v].record_id == 1)
    {
      is_raycloud = true;
      break;
    }
  }
  // Fallback for older RCT files written without the VLR marker: detect by "sx" attribute name.
  if (!is_raycloud)
  {
    for (laszip_U32 v = 0; v < header->number_of_variable_length_records; v++)
    {
      auto &vlr = header->vlrs[v];
      if (strcmp(vlr.user_id, "LASF_Spec") != 0 || vlr.record_id != 4)
        continue;
      const int num_attrs = vlr.record_length_after_header / 192;
      for (int a = 0; a < num_attrs; a++)
      {
        char attr_name[33] = {};
        memcpy(attr_name, vlr.data + a * 192 + 4, 32);
        if (strcmp(attr_name, "sx") == 0) { is_raycloud = true; break; }
      }
      break;
    }
  }

  // LAS EXTRA_BYTES data_type → per-point byte size (types 0 and >10 are skipped)
  static const uint16_t kExtraTypeSize[11] = { 0, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };
  // Names of raycloud-owned extra attributes (these are skipped when extracting original data)
  static const char *kRayCloudAttrs[] = { "sx", "sy", "sz", "alpha", "bound", "tree_id", "stem_id", "beam_id" };

  uint16_t local_skip_size = 0;   // bytes of our own extra attributes before original data
  uint16_t local_orig_extra = 0;  // bytes of original sensor data per point
  uint16_t own_offset      = 0;   // running byte cursor through our own attrs (in declared order)
  uint16_t tree_id_offset  = 0;   uint8_t tree_id_dtype = 0;  // 0 = absent
  uint16_t stem_id_offset  = 0;   uint8_t stem_id_dtype = 0;  // 0 = absent
  uint16_t beam_id_offset  = 0;   uint8_t beam_id_dtype = 0;  // 0 = absent
  uint16_t alpha_offset    = 12;  // default: sx+sy+sz only; overwritten when "alpha" VLR found
  int32_t  bound_offset    = -1;  // -1 = absent (old file); set when "bound" VLR found
  std::vector<uint8_t> local_orig_vlr;

  for (laszip_U32 v = 0; v < header->number_of_variable_length_records; v++)
  {
    auto &vlr = header->vlrs[v];
    if (strcmp(vlr.user_id, "LASF_Spec") != 0 || vlr.record_id != 4)
      continue;
    const int num_attrs = vlr.record_length_after_header / 192;
    for (int a = 0; a < num_attrs; a++)
    {
      const uint8_t *rec = vlr.data + a * 192;
      const uint8_t dtype = rec[2];
      const uint16_t attr_size = (dtype > 0 && dtype <= 10) ? kExtraTypeSize[dtype] : 0;
      if (attr_size == 0)
        continue;
      char attr_name[33] = {};
      memcpy(attr_name, rec + 4, 32);
      bool is_ours = false;
      if (is_raycloud)
      {
        for (const char *own : kRayCloudAttrs)
          if (strcmp(attr_name, own) == 0) { is_ours = true; break; }
        if (strcmp(attr_name, "tree_id") == 0) { tree_id_offset = own_offset; tree_id_dtype = dtype; }
        if (strcmp(attr_name, "stem_id") == 0) { stem_id_offset = own_offset; stem_id_dtype = dtype; }
        if (strcmp(attr_name, "beam_id") == 0) { beam_id_offset = own_offset; beam_id_dtype = dtype; }
        if (strcmp(attr_name, "alpha")   == 0) { alpha_offset   = own_offset; }
        if (strcmp(attr_name, "bound")   == 0) { bound_offset   = own_offset; }
      }
      if (is_ours)
      {
        own_offset       += attr_size;
        local_skip_size  += attr_size;
      }
      else
      {
        local_orig_extra += attr_size;
        local_orig_vlr.insert(local_orig_vlr.end(), rec, rec + 192);
      }
    }
    break; // only one EXTRA_BYTES VLR
  }

  if (orig_extra_size_out)
    *orig_extra_size_out = local_orig_extra;
  if (extra_bytes_vlr_out)
    *extra_bytes_vlr_out = local_orig_vlr;

  // Gather the per-file decode parameters into a single context so the per-point decode is a pure
  // function of (point record, context). Used by decodePointRecord below and by the parallel paths.
  DecodeContext ctx;
  ctx.format = format;
  ctx.using_colour = using_colour;
  ctx.is_raycloud = is_raycloud;
  ctx.max_intensity = max_intensity;
  ctx.local_skip_size = local_skip_size;
  ctx.local_orig_extra = local_orig_extra;
  ctx.tree_id_offset = tree_id_offset;  ctx.tree_id_dtype = tree_id_dtype;
  ctx.stem_id_offset = stem_id_offset;  ctx.stem_id_dtype = stem_id_dtype;
  ctx.beam_id_offset = beam_id_offset;  ctx.beam_id_dtype = beam_id_dtype;
  ctx.alpha_offset = alpha_offset;
  ctx.bound_offset = bound_offset;
  ctx.scale[0] = header->x_scale_factor;
  ctx.scale[1] = header->y_scale_factor;
  ctx.scale[2] = header->z_scale_factor;
  ctx.offset[0] = header->x_offset;
  ctx.offset[1] = header->y_offset;
  ctx.offset[2] = header->z_offset;
  ctx.point_record_length = header->point_data_record_length;
  ctx.extra_bytes_total = static_cast<uint16_t>(local_skip_size + local_orig_extra);
  // The extra-bytes block sits at the tail of each fixed record, immediately after the base fields.
  ctx.extra_bytes_offset = static_cast<uint16_t>(header->point_data_record_length - ctx.extra_bytes_total);

  ray::Progress progress;
  ray::ProgressThread progress_thread(progress);
  const size_t num_chunks = (number_of_points + (chunk_size - 1)) / chunk_size;
  chunk_size = std::min(number_of_points, chunk_size);
  progress.begin("read and process", num_chunks);

  std::vector<Eigen::Vector3d> starts;
  std::vector<Eigen::Vector3d> ends;
  std::vector<double> times;
  std::vector<RGBA> colours;
  std::vector<uint8_t> intensities;
  starts.reserve(chunk_size);
  ends.reserve(chunk_size);
  times.reserve(chunk_size);
  intensities.reserve(chunk_size);
  colours.reserve(chunk_size);

  num_bounded = 0;

  // Index-addressed buffers shared by the parallel/decompressed fast paths. They are pre-sized to
  // number_of_points so a decode can fill any slot, then drained through @c apply in the same chunk
  // windows as the sequential path so the callback contract is identical.
  std::vector<Eigen::Vector3d> all_starts, all_ends;
  std::vector<double> all_times;
  std::vector<RGBA> all_colours;
  std::vector<uint8_t> all_intensities;
  IndexedDecodeBuffers buf;

  // Populate @c buf and pre-size every active output for an index-addressed decode of all points.
  auto setup_indexed_buffers = [&]() {
    all_starts.resize(number_of_points);
    all_ends.resize(number_of_points);
    all_times.resize(number_of_points);
    all_colours.resize(using_colour ? number_of_points : 0);
    all_intensities.resize(number_of_points);
    buf.starts = all_starts.data();
    buf.ends = all_ends.data();
    buf.times = all_times.data();
    buf.colours = using_colour ? all_colours.data() : nullptr;
    buf.intensities = all_intensities.data();

    const uint16_t pstride = static_cast<uint16_t>(10 + local_orig_extra);
    if (passthrough_out)
    {
      passthrough_out->resize(static_cast<size_t>(pstride) * number_of_points);
      buf.passthrough = passthrough_out->data();
      buf.passthrough_stride = pstride;
    }
    // ID outputs are present iff the sequential path would have produced them. These conditions
    // depend only on per-file state (num_extra_bytes is constant), so they hold for every point.
    const bool tree_ok =
      tree_ids_out && is_raycloud && tree_id_dtype != 0 &&
      ctx.extra_bytes_total >= tree_id_offset + kExtraTypeSize[tree_id_dtype];
    const bool stem_ok = stem_ids_out && is_raycloud && tree_id_dtype != 0;  // pushes value or 0
    const bool beam_ok =
      beam_ids_out && is_raycloud && beam_id_dtype != 0 &&
      ctx.extra_bytes_total >= beam_id_offset + kExtraTypeSize[beam_id_dtype];
    if (tree_ok) { tree_ids_out->resize(number_of_points); buf.tree_ids = tree_ids_out->data(); buf.tree_active = true; }
    if (stem_ok) { stem_ids_out->resize(number_of_points); buf.stem_ids = stem_ids_out->data(); buf.stem_active = true; }
    if (beam_ok) { beam_ids_out->resize(number_of_points); buf.beam_ids = beam_ids_out->data(); buf.beam_active = true; }
  };

  // Drain the index-addressed buffers through @c apply in the same chunk windows (and with the same
  // colour/alpha handling and progress cadence) as the sequential path.
  //
  // Producer-consumer pipeline: a background reader thread assembles each chunk window (the per-chunk
  // copy out of the bulk-decoded arrays, colourByTime fallback, and alpha merge) into a ChunkBuffer
  // and hands it to the main thread via a bounded double-buffer queue. The main thread pops each
  // buffer and invokes @c apply. This overlaps chunk assembly with @c apply while guaranteeing
  // @c apply only ever runs on the main thread (callers mutate non-thread-safe external state in it).
  auto flush_indexed = [&]() {
    // The set of fields apply consumes per chunk; ownership passes from reader thread to main thread.
    struct ChunkBuffer
    {
      std::vector<Eigen::Vector3d> starts, ends;
      std::vector<double> times;
      std::vector<RGBA> colours;
      bool is_last = false;
    };
    ThreadSafeQueue<std::shared_ptr<ChunkBuffer>> queue(2);  // double-buffer

    std::thread reader([&]() {
      for (size_t base = 0; base < number_of_points; base += chunk_size)
      {
        const size_t n = std::min(chunk_size, number_of_points - base);
        auto cb = std::make_shared<ChunkBuffer>();
        cb->starts.assign(all_starts.begin() + base, all_starts.begin() + base + n);
        cb->ends.assign(all_ends.begin() + base, all_ends.begin() + base + n);
        cb->times.assign(all_times.begin() + base, all_times.begin() + base + n);
        if (using_colour)
          cb->colours.assign(all_colours.begin() + base, all_colours.begin() + base + n);
        else
          colourByTime(cb->times, cb->colours);
        for (size_t j = 0; j < cb->colours.size(); j++)
          cb->colours[j].alpha = all_intensities[base + j];
        cb->is_last = (base + n >= number_of_points);
        queue.push(std::move(cb));
      }
      queue.notify_done();
    });

    std::shared_ptr<ChunkBuffer> cb;
    while (queue.pop(cb))
    {
      apply(cb->starts, cb->ends, cb->times, cb->colours);
      progress.increment();
      cb.reset();  // free the buffer back to the allocator before popping the next
    }
    reader.join();
  };

  // Fast path: for uncompressed LAS with a fixed-layout record we know how to decode, mmap the file
  // and decode all records in parallel, then drive @c apply in the same chunk windows as the
  // sequential path. Any unsupported condition or mmap failure falls through to the laszip loop.
  bool fast_path_done = false;
  if (!is_compressed && number_of_points > 0 && lasBaseRecordSize(format) != 0 &&
      !getenv("RAYLAS_NO_MMAP") &&
      ctx.point_record_length >= lasBaseRecordSize(format) + ctx.extra_bytes_total)
  {
#ifdef _WIN32
    HANDLE fh = CreateFileA(file_name.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING,
                            FILE_ATTRIBUTE_NORMAL, nullptr);
    HANDLE map = nullptr;
    const uint8_t *mmap_ptr = nullptr;
    size_t file_size = 0;
    if (fh != INVALID_HANDLE_VALUE)
    {
      LARGE_INTEGER sz;
      if (GetFileSizeEx(fh, &sz))
      {
        file_size = static_cast<size_t>(sz.QuadPart);
        map = CreateFileMappingA(fh, nullptr, PAGE_READONLY, 0, 0, nullptr);
        if (map)
          mmap_ptr = static_cast<const uint8_t *>(MapViewOfFile(map, FILE_MAP_READ, 0, 0, 0));
      }
    }
#else
    const int fd = ::open(file_name.c_str(), O_RDONLY);
    const uint8_t *mmap_ptr = nullptr;
    size_t file_size = 0;
    if (fd >= 0)
    {
      struct stat st;
      if (::fstat(fd, &st) == 0)
      {
        file_size = static_cast<size_t>(st.st_size);
        void *p = ::mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
        if (p != MAP_FAILED)
          mmap_ptr = static_cast<const uint8_t *>(p);
      }
    }
#endif
    const size_t record_len = ctx.point_record_length;
    const size_t data_off = static_cast<size_t>(header->offset_to_point_data);
    const bool bounds_ok =
      mmap_ptr != nullptr && file_size >= data_off + record_len * number_of_points;

    if (bounds_ok)
    {
      const uint8_t *raw = mmap_ptr + data_off;
      setup_indexed_buffers();

      size_t bounded_count = 0;
#pragma omp parallel reduction(+ : bounded_count)
      {
        laszip_point_struct pt;
        std::memset(&pt, 0, sizeof(pt));
        std::vector<uint8_t> extra_scratch(ctx.extra_bytes_total ? ctx.extra_bytes_total : 1);
#pragma omp for schedule(static)
        for (laszip_I64 i = 0; i < static_cast<laszip_I64>(number_of_points); ++i)
        {
          const uint8_t *rec = raw + static_cast<size_t>(i) * record_len;
          fillPointFromRecord(rec, ctx, pt, extra_scratch.data());
          int32_t xi = pt.X, yi = pt.Y, zi = pt.Z;
          Eigen::Vector3d position(xi * ctx.scale[0] + ctx.offset[0], yi * ctx.scale[1] + ctx.offset[1],
                                   zi * ctx.scale[2] + ctx.offset[2]);
          uint8_t bounded;
          decodePointRecordIndexed(&pt, ctx, position, static_cast<size_t>(i), buf, bounded);
          bounded_count += bounded;
        }
      }
      num_bounded = bounded_count;

      flush_indexed();
      fast_path_done = true;
    }

    // Release the mapping regardless of whether decode succeeded.
#ifdef _WIN32
    if (mmap_ptr) UnmapViewOfFile(const_cast<uint8_t *>(mmap_ptr));
    if (map) CloseHandle(map);
    if (fh != INVALID_HANDLE_VALUE) CloseHandle(fh);
#else
    if (mmap_ptr) ::munmap(const_cast<uint8_t *>(mmap_ptr), file_size);
    if (fd >= 0) ::close(fd);
#endif
  }

#if RAYLIB_WITH_LAZPERF
  // Fast path for compressed (LAZ) files: decompress records via laz-perf instead of the laszip
  // per-point loop, then decode each fixed record with the shared decode core. Single-threaded
  // decompression here; the producer-consumer pipeline overlaps it with apply in a later phase.
  if (is_compressed && !fast_path_done && number_of_points > 0 && lasBaseRecordSize(format) != 0 &&
      !getenv("RAYLAS_NO_LAZPERF") &&
      ctx.point_record_length >= lasBaseRecordSize(format) + ctx.extra_bytes_total)
  {
    try
    {
      lazperf::reader::named_file lazf(file_name);
      const auto &lazhdr = lazf.header();
      if (lazhdr.point_record_length == ctx.point_record_length &&
          lazhdr.point_count == number_of_points)
      {
        setup_indexed_buffers();
        std::vector<uint8_t> record(ctx.point_record_length);
        std::vector<uint8_t> extra_scratch(ctx.extra_bytes_total ? ctx.extra_bytes_total : 1);
        laszip_point_struct pt;
        std::memset(&pt, 0, sizeof(pt));
        size_t bounded_count = 0;
        for (size_t i = 0; i < number_of_points; ++i)
        {
          lazf.readPoint(reinterpret_cast<char *>(record.data()));
          fillPointFromRecord(record.data(), ctx, pt, extra_scratch.data());
          int32_t xi = pt.X, yi = pt.Y, zi = pt.Z;
          Eigen::Vector3d position(xi * ctx.scale[0] + ctx.offset[0], yi * ctx.scale[1] + ctx.offset[1],
                                   zi * ctx.scale[2] + ctx.offset[2]);
          uint8_t bounded;
          decodePointRecordIndexed(&pt, ctx, position, i, buf, bounded);
          bounded_count += bounded;
        }
        num_bounded = bounded_count;
        flush_indexed();
        fast_path_done = true;
      }
    }
    catch (const std::exception &e)
    {
      // Any laz-perf failure falls through to the laszip per-point loop below.
      std::cerr << "readLas: laz-perf decode failed (" << e.what() << "), using laszip" << std::endl;
    }
  }
#endif  // RAYLIB_WITH_LAZPERF

  for (size_t i = 0; !fast_path_done && i < number_of_points; i++)
  {
    if (laszip_read_point(reader))
    {
      laszip_CHAR *error;
      laszip_get_error(reader, &error);
      std::cerr << "readLas: error reading point " << i << ": " << error << std::endl;
      break;
    }

    laszip_F64 coords[3];
    laszip_get_coordinates(reader, coords);
    Eigen::Vector3d position(coords[0], coords[1], coords[2]);

    decodePointRecord(point, ctx, position, starts, ends, times, colours, intensities, num_bounded,
                      tree_ids_out, passthrough_out, stem_ids_out, beam_ids_out);

    if (ends.size() == chunk_size || i == number_of_points - 1)
    {
      if (colours.empty())
        colourByTime(times, colours);
      for (size_t j = 0; j < colours.size(); j++)
        colours[j].alpha = intensities[j];
      apply(starts, ends, times, colours);
      starts.clear();
      ends.clear();
      times.clear();
      colours.clear();
      intensities.clear();
      progress.increment();
    }
  }

  progress.end();
  progress_thread.requestQuit();
  progress_thread.join();

  laszip_close_reader(reader);
  laszip_destroy(reader);

  std::cout << "loaded " << file_name << " with " << number_of_points << " points" << std::endl;
  return true;
#else   // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(max_intensity);
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(apply);
  RAYLIB_UNUSED(num_bounded);
  RAYLIB_UNUSED(chunk_size);
  RAYLIB_UNUSED(tree_ids_out);
  RAYLIB_UNUSED(passthrough_out);
  RAYLIB_UNUSED(orig_extra_size_out);
  RAYLIB_UNUSED(extra_bytes_vlr_out);
  RAYLIB_UNUSED(stem_ids_out);
  RAYLIB_UNUSED(beam_ids_out);
  std::cerr << "readLas: cannot read file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

bool readLasExtraBytesVlr(const std::string &file_name, uint16_t &orig_extra_size_out,
                           std::vector<uint8_t> &extra_bytes_vlr_out, bool *has_bound_out)
{
#if RAYLIB_WITH_LAS
  if (has_bound_out)
    *has_bound_out = false;
  laszip_POINTER reader;
  if (laszip_create(&reader))
    return false;

  laszip_BOOL is_compressed;
  if (laszip_open_reader(reader, file_name.c_str(), &is_compressed))
  {
    laszip_destroy(reader);
    return false;
  }

  laszip_header_struct *header;
  laszip_get_header_pointer(reader, &header);

  bool is_raycloud = false;
  for (laszip_U32 v = 0; v < header->number_of_variable_length_records; v++)
  {
    if (strncmp(header->vlrs[v].user_id, "raycloudtools", 16) == 0 && header->vlrs[v].record_id == 1)
    {
      is_raycloud = true;
      break;
    }
  }
  if (!is_raycloud)
  {
    for (laszip_U32 v = 0; v < header->number_of_variable_length_records; v++)
    {
      auto &vlr = header->vlrs[v];
      if (strcmp(vlr.user_id, "LASF_Spec") != 0 || vlr.record_id != 4)
        continue;
      const int num_attrs = vlr.record_length_after_header / 192;
      for (int a = 0; a < num_attrs; a++)
      {
        char attr_name[33] = {};
        memcpy(attr_name, vlr.data + a * 192 + 4, 32);
        if (strcmp(attr_name, "sx") == 0) { is_raycloud = true; break; }
      }
      break;
    }
  }

  static const uint16_t kExtraTypeSize[11] = { 0, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };
  static const char *kRayCloudAttrs[] = { "sx", "sy", "sz", "alpha", "bound", "tree_id", "stem_id", "beam_id" };

  uint16_t local_orig_extra = 0;
  std::vector<uint8_t> local_orig_vlr;

  for (laszip_U32 v = 0; v < header->number_of_variable_length_records; v++)
  {
    auto &vlr = header->vlrs[v];
    if (strcmp(vlr.user_id, "LASF_Spec") != 0 || vlr.record_id != 4)
      continue;
    const int num_attrs = vlr.record_length_after_header / 192;
    for (int a = 0; a < num_attrs; a++)
    {
      const uint8_t *rec = vlr.data + a * 192;
      const uint8_t dtype = rec[2];
      const uint16_t attr_size = (dtype > 0 && dtype <= 10) ? kExtraTypeSize[dtype] : 0;
      if (attr_size == 0)
        continue;
      char attr_name[33] = {};
      memcpy(attr_name, rec + 4, 32);
      if (has_bound_out && strcmp(attr_name, "bound") == 0)
        *has_bound_out = true;
      bool is_ours = false;
      if (is_raycloud)
        for (const char *own : kRayCloudAttrs)
          if (strcmp(attr_name, own) == 0) { is_ours = true; break; }
      if (!is_ours)
      {
        local_orig_extra += attr_size;
        local_orig_vlr.insert(local_orig_vlr.end(), rec, rec + 192);
      }
    }
    break;
  }

  laszip_close_reader(reader);
  laszip_destroy(reader);

  orig_extra_size_out = local_orig_extra;
  extra_bytes_vlr_out = std::move(local_orig_vlr);
  return true;
#else
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(orig_extra_size_out);
  RAYLIB_UNUSED(extra_bytes_vlr_out);
  RAYLIB_UNUSED(has_bound_out);
  return false;
#endif
}

bool readLas(std::string file_name, std::vector<Eigen::Vector3d> &positions, std::vector<double> &times,
             std::vector<RGBA> &colours, double max_intensity, Eigen::Vector3d *offset_to_remove)
{
  std::vector<Eigen::Vector3d> starts;  // dummy as lax just reads in point clouds, not ray clouds
  auto apply = [&](std::vector<Eigen::Vector3d> &start_points, std::vector<Eigen::Vector3d> &end_points,
                   std::vector<double> &time_points, std::vector<RGBA> &colour_values)
  {
    starts.insert(starts.end(), start_points.begin(), start_points.end());
    positions.insert(positions.end(), end_points.begin(), end_points.end());
    times.insert(times.end(), time_points.begin(), time_points.end());
    colours.insert(colours.end(), colour_values.begin(), colour_values.end());
  };
  size_t num_bounded;
  bool success =
    readLas(file_name, apply, num_bounded, max_intensity, offset_to_remove, std::numeric_limits<size_t>::max());
  if (num_bounded == 0)
  {
    std::cout << "warning: all laz file intensities are 0, which would make all rays unbounded. Setting them to 1."
              << std::endl;
    for (auto &c : colours) c.alpha = 255;
  }
  return success;
}

bool RAYLIB_EXPORT writeLas(std::string file_name, const std::vector<Eigen::Vector3d> &points,
                            const std::vector<double> &times, const std::vector<RGBA> &colours)
{
#if RAYLIB_WITH_LAS
  std::cout << "saving LAZ file" << std::endl;

  laszip_POINTER writer;
  if (laszip_create(&writer))
  {
    std::cerr << "writeLas: failed to create LASzip writer" << std::endl;
    return false;
  }

  laszip_header_struct *header;
  laszip_get_header_pointer(writer, &header);

  header->version_major = 1;
  header->version_minor = 4;
  header->header_size = 375;  // LAS 1.4 header is 375 bytes
  header->point_data_format = 6;  // LAS 1.4: GPS time only
  const double scale = 1e-4;
  header->x_scale_factor = scale;
  header->y_scale_factor = scale;
  header->z_scale_factor = scale;
  header->x_offset = 0.0;
  header->y_offset = 0.0;
  header->z_offset = 0.0;
  header->extended_number_of_point_records = static_cast<laszip_U64>(points.size());
  header->offset_to_point_data = 375;  // LAS 1.4 header, no VLRs

  const bool is_laz = file_name.find(".laz") != std::string::npos;
  std::cout << "Saving points to " << file_name << std::endl;

  if (laszip_open_writer(writer, file_name.c_str(), is_laz ? 1 : 0))
  {
    laszip_CHAR *error;
    laszip_get_error(writer, &error);
    std::cerr << "writeLas: failed to open file for writing: " << error << std::endl;
    laszip_destroy(writer);
    return false;
  }

  laszip_point_struct *point;
  laszip_get_point_pointer(writer, &point);

  for (size_t i = 0; i < points.size(); i++)
  {
    laszip_F64 coords[3] = { points[i][0], points[i][1], points[i][2] };
    laszip_set_coordinates(writer, coords);
    point->intensity = colours[i].alpha;
    if (!times.empty())
      point->gps_time = times[i];
    laszip_write_point(writer);
  }

  laszip_update_inventory(writer);
  laszip_close_writer(writer);
  laszip_destroy(writer);
  return true;
#else   // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(points);
  RAYLIB_UNUSED(times);
  RAYLIB_UNUSED(colours);
  std::cerr << "writeLas: cannot write file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true"
            << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

#if RAYLIB_WITH_LAS
LasWriter::LasWriter(const std::string &file_name)
  : file_name_(file_name)
  , writer_handle_(nullptr)
  , header_(nullptr)
  , point_(nullptr)
  , points_written_(0)
{
  if (laszip_create(&writer_handle_))
  {
    std::cerr << "LasWriter: failed to create LASzip writer" << std::endl;
    writer_handle_ = nullptr;
    return;
  }

  laszip_get_header_pointer(writer_handle_, &header_);

  header_->version_major = 1;
  header_->version_minor = 4;
  header_->header_size = 375;  // LAS 1.4 header is 375 bytes
  header_->point_data_format = 6;  // LAS 1.4: GPS time only
  const double scale = 1e-4;
  header_->x_scale_factor = scale;
  header_->y_scale_factor = scale;
  header_->z_scale_factor = scale;
  header_->x_offset = 0.0;
  header_->y_offset = 0.0;
  header_->z_offset = 0.0;
  header_->offset_to_point_data = 375;  // LAS 1.4 header, no VLRs

  const bool is_laz = file_name_.find(".laz") != std::string::npos;
  std::cout << "Saving points to " << file_name_ << std::endl;

  if (laszip_open_writer(writer_handle_, file_name_.c_str(), is_laz ? 1 : 0))
  {
    laszip_CHAR *error;
    laszip_get_error(writer_handle_, &error);
    std::cerr << "LasWriter: failed to open file for writing: " << error << std::endl;
    laszip_destroy(writer_handle_);
    writer_handle_ = nullptr;
    return;
  }

  laszip_get_point_pointer(writer_handle_, &point_);
}
#else   // RAYLIB_WITH_LAS
LasWriter::LasWriter(const std::string &file_name)
  : file_name_(file_name)
{
  RAYLIB_UNUSED(file_name);
  std::cerr << "writeLas: cannot write file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true"
            << std::endl;
}
#endif  // RAYLIB_WITH_LAS

LasWriter::~LasWriter()
{
#if RAYLIB_WITH_LAS
  if (writer_handle_)
  {
    laszip_update_inventory(writer_handle_);
    laszip_close_writer(writer_handle_);
    laszip_destroy(writer_handle_);
    // laszip_close_writer clobbers point counts for streaming writes. Patch both the legacy
    // 32-bit count (offset 107) and the LAS 1.4 64-bit extended count (offset 247) on disk.
    if (points_written_ > 0)
    {
      std::fstream f(file_name_, std::ios::in | std::ios::out | std::ios::binary);
      if (f.is_open())
      {
        const laszip_U32 legacy = static_cast<laszip_U32>(
          std::min<uint64_t>(points_written_, std::numeric_limits<laszip_U32>::max()));
        f.seekp(107);
        f.write(reinterpret_cast<const char *>(&legacy), sizeof(legacy));
        const laszip_U64 extended = static_cast<laszip_U64>(points_written_);
        f.seekp(247);
        f.write(reinterpret_cast<const char *>(&extended), sizeof(extended));
      }
    }
  }
#else
  std::cerr << "writeLas: cannot write file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true"
            << std::endl;
#endif
}

bool LasWriter::writeChunk(const std::vector<Eigen::Vector3d> &points, const std::vector<double> &times,
                           const std::vector<RGBA> &colours)
{
#if RAYLIB_WITH_LAS
  if (points.size() == 0)
  {
    return true;  // this is acceptable behaviour. It avoids calling function checking for emptiness each time
  }
  if (!writer_handle_ || !point_)
  {
    std::cerr << "Error: cannot open " << file_name_ << " for writing." << std::endl;
    return false;
  }
  for (size_t i = 0; i < points.size(); i++)
  {
    laszip_F64 coords[3] = { points[i][0], points[i][1], points[i][2] };
    laszip_set_coordinates(writer_handle_, coords);
    point_->intensity = colours[i].alpha;
    if (!times.empty())
      point_->gps_time = times[i];
    laszip_write_point(writer_handle_);
  }
  points_written_ += points.size();
  return true;
#else   // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(points);
  RAYLIB_UNUSED(times);
  RAYLIB_UNUSED(colours);
  std::cerr << "writeLas: cannot write file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true"
            << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

bool RAYLIB_EXPORT writeLasRayCloud(const std::string &file_name, const std::vector<Eigen::Vector3d> &starts,
                                    const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                                    const std::vector<RGBA> &colours, const std::vector<int32_t> &tree_ids,
                                    const std::vector<int32_t> &stem_ids,
                                    const std::vector<uint8_t> &passthrough,
                                    const std::vector<uint8_t> &extra_bytes_vlr)
{
#if RAYLIB_WITH_LAS
  LasRayCloudWriter writer(file_name, !tree_ids.empty(), !stem_ids.empty(), extra_bytes_vlr);
  return writer.writeChunk(starts, ends, times, colours, tree_ids, stem_ids, passthrough);
#else   // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(starts);
  RAYLIB_UNUSED(ends);
  RAYLIB_UNUSED(times);
  RAYLIB_UNUSED(colours);
  RAYLIB_UNUSED(stem_ids);
  std::cerr << "writeLasRayCloud: WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

#if RAYLIB_WITH_LAS
LasRayCloudWriter::LasRayCloudWriter(const std::string &file_name, bool with_tree_id, bool with_stem_id,
                                     const std::vector<uint8_t> &extra_bytes_vlr, bool with_beam_id)
  : file_name_(file_name)
  , points_written_(0)
  , with_tree_id_(with_tree_id)
  , with_stem_id_(with_stem_id)
  , with_beam_id_(with_beam_id)
  , orig_extra_size_(0)
  , passthrough_stride_(10)
  , writer_handle_(nullptr)
  , point_(nullptr)
{
  if (laszip_create(&writer_handle_))
  {
    std::cerr << "LasRayCloudWriter: failed to create LASzip writer" << std::endl;
    writer_handle_ = nullptr;
    return;
  }

  laszip_header_struct *header;
  laszip_get_header_pointer(writer_handle_, &header);

  header->version_major = 1;
  header->version_minor = 4;
  header->header_size = 375;
  header->point_data_format = 7;
  const double scale = 1e-4;
  header->x_scale_factor = scale;
  header->y_scale_factor = scale;
  header->z_scale_factor = scale;
  header->x_offset = 0.0;
  header->y_offset = 0.0;
  header->z_offset = 0.0;

  // LASzip API type 8 = F32, type 5 = INT32, type 0 = U8.
  bool attr_err =
    laszip_add_attribute(writer_handle_, 8, "sx", "ray start x offset", 1.0, 0.0) ||
    laszip_add_attribute(writer_handle_, 8, "sy", "ray start y offset", 1.0, 0.0) ||
    laszip_add_attribute(writer_handle_, 8, "sz", "ray start z offset", 1.0, 0.0);
  if (!attr_err && with_tree_id_)
    attr_err = laszip_add_attribute(writer_handle_, 5, "tree_id", "per-point tree ID", 1.0, 0.0);
  if (!attr_err && with_stem_id_)
    attr_err = laszip_add_attribute(writer_handle_, 5, "stem_id", "per-point stem ID", 1.0, 0.0);
  if (!attr_err && with_beam_id_)
    attr_err = laszip_add_attribute(writer_handle_, 5, "beam_id", "per-pulse beam ID", 1.0, 0.0);
  if (!attr_err)
    attr_err = laszip_add_attribute(writer_handle_, 0, "alpha", "intensity 1-255", 1.0, 0.0);
  if (!attr_err)
    attr_err = laszip_add_attribute(writer_handle_, 0, "bound", "1=bound, 0=unbound", 1.0, 0.0);
  if (attr_err)
  {
    laszip_CHAR *error;
    laszip_get_error(writer_handle_, &error);
    std::cerr << "LasRayCloudWriter: failed to add extra attributes: " << error << std::endl;
    laszip_destroy(writer_handle_);
    writer_handle_ = nullptr;
    return;
  }

  // Register original sensor extra-byte attributes from the EXTRA_BYTES VLR payload.
  // LAS spec data_type -> size: 1=u8(1), 2=i8(1), 3=u16(2), 4=i16(2), 5=u32(4), 6=i32(4),
  //                             7=u64(8), 8=i64(8), 9=float(4), 10=double(8)
  // LASzip API type = LAS spec data_type - 1.
  static const uint16_t kTypeSize[11] = { 0, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };
  for (size_t off = 0; off + 192 <= extra_bytes_vlr.size(); off += 192)
  {
    const uint8_t *rec = extra_bytes_vlr.data() + off;
    const uint8_t dtype = rec[2];
    if (dtype == 0 || dtype > 10)
      continue;
    char name[33] = {}, desc[33] = {};
    std::memcpy(name, rec + 4,   32);
    std::memcpy(desc, rec + 160, 32);
    const uint8_t opts = rec[3];
    double scale_v = 1.0, offset_v = 0.0;
    if (opts & 0x08) std::memcpy(&scale_v,  rec + 112, 8);
    if (opts & 0x10) std::memcpy(&offset_v, rec + 136, 8);
    laszip_add_attribute(writer_handle_, static_cast<laszip_U32>(dtype - 1), name, desc, scale_v, offset_v);
    orig_extra_size_ += kTypeSize[dtype];
  }
  passthrough_stride_ = static_cast<uint16_t>(10 + orig_extra_size_);

  // LAS 1.4 format 7 base = 36 bytes.
  uint16_t extra = 12; // sx, sy, sz
  if (with_tree_id_) extra += 4;
  if (with_stem_id_) extra += 4;
  if (with_beam_id_) extra += 4;
  extra += 1; // alpha
  extra += 1; // bound
  extra += orig_extra_size_;
  const laszip_U16 record_size = static_cast<laszip_U16>(36 + extra);
  if (laszip_set_point_type_and_size(writer_handle_, 7, record_size))
  {
    laszip_CHAR *error;
    laszip_get_error(writer_handle_, &error);
    std::cerr << "LasRayCloudWriter: failed to set point type/size: " << error << std::endl;
    laszip_destroy(writer_handle_);
    writer_handle_ = nullptr;
    return;
  }

  if (laszip_add_vlr(writer_handle_, "raycloudtools", 1, 0, "raycloud", nullptr))
  {
    laszip_CHAR *error;
    laszip_get_error(writer_handle_, &error);
    std::cerr << "LasRayCloudWriter: failed to add VLR: " << error << std::endl;
    laszip_destroy(writer_handle_);
    writer_handle_ = nullptr;
    return;
  }

  // LASzip uses LAS 1.2 default header size (227) when computing offset_to_point_data
  // during VLR additions. Correct for the LAS 1.4 header (375 bytes).
  laszip_header_struct *hdr;
  laszip_get_header_pointer(writer_handle_, &hdr);
  hdr->offset_to_point_data += (375 - 227);

  const bool is_laz = file_name_.find(".laz") != std::string::npos;
  std::cout << "Saving ray cloud to " << file_name_ << std::endl;

  if (laszip_open_writer(writer_handle_, file_name_.c_str(), is_laz ? 1 : 0))
  {
    laszip_CHAR *error;
    laszip_get_error(writer_handle_, &error);
    std::cerr << "LasRayCloudWriter: failed to open file for writing: " << error << std::endl;
    laszip_destroy(writer_handle_);
    writer_handle_ = nullptr;
    return;
  }

  laszip_get_point_pointer(writer_handle_, &point_);
}
#else   // RAYLIB_WITH_LAS
LasRayCloudWriter::LasRayCloudWriter(const std::string &file_name, bool with_tree_id, bool with_stem_id,
                                     const std::vector<uint8_t> &extra_bytes_vlr, bool with_beam_id)
  : file_name_(file_name)
  , with_tree_id_(with_tree_id)
  , with_stem_id_(with_stem_id)
  , with_beam_id_(with_beam_id)
{
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(with_tree_id);
  RAYLIB_UNUSED(with_stem_id);
  RAYLIB_UNUSED(extra_bytes_vlr);
  RAYLIB_UNUSED(with_beam_id);
  std::cerr << "LasRayCloudWriter: WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
}
#endif  // RAYLIB_WITH_LAS

LasRayCloudWriter::~LasRayCloudWriter()
{
#if RAYLIB_WITH_LAS
  if (writer_handle_)
  {
    laszip_update_inventory(writer_handle_);
    laszip_close_writer(writer_handle_);
    laszip_destroy(writer_handle_);
    // laszip_close_writer clobbers point counts for streaming writes. Patch both the legacy
    // 32-bit count (offset 107) and the LAS 1.4 64-bit extended count (offset 247) on disk.
    if (points_written_ > 0)
    {
      std::fstream f(file_name_, std::ios::in | std::ios::out | std::ios::binary);
      if (f.is_open())
      {
        const laszip_U32 legacy = static_cast<laszip_U32>(
          std::min<uint64_t>(points_written_, std::numeric_limits<laszip_U32>::max()));
        f.seekp(107);
        f.write(reinterpret_cast<const char *>(&legacy), sizeof(legacy));
        const laszip_U64 extended = static_cast<laszip_U64>(points_written_);
        f.seekp(247);
        f.write(reinterpret_cast<const char *>(&extended), sizeof(extended));
      }
    }
  }
#else
  std::cerr << "LasRayCloudWriter: WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
#endif
}

bool LasRayCloudWriter::writeChunk(const std::vector<Eigen::Vector3d> &starts,
                                   const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                                   const std::vector<RGBA> &colours, const std::vector<int32_t> &tree_ids,
                                   const std::vector<int32_t> &stem_ids,
                                   const std::vector<uint8_t> &passthrough,
                                   const std::vector<int32_t> &beam_ids)
{
#if RAYLIB_WITH_LAS
  if (ends.empty())
    return true;
  if (!writer_handle_ || !point_)
  {
    std::cerr << "Error: LasRayCloudWriter not open for writing to " << file_name_ << std::endl;
    return false;
  }
  for (size_t i = 0; i < ends.size(); i++)
  {
    laszip_F64 coords[3] = { ends[i][0], ends[i][1], ends[i][2] };
    laszip_set_coordinates(writer_handle_, coords);
    point_->gps_time = times[i];
    point_->intensity = colours[i].alpha;
    point_->rgb[0] = static_cast<laszip_U16>(colours[i].red) * 257u;
    point_->rgb[1] = static_cast<laszip_U16>(colours[i].green) * 257u;
    point_->rgb[2] = static_cast<laszip_U16>(colours[i].blue) * 257u;
    // Restore LAS fields from passthrough if available.
    // Layout: [0-7] standard LAS fields, [8-9] original intensity (uint16 LE), [10..] sensor extras.
    if (passthrough.size() >= (i + 1) * passthrough_stride_)
    {
      const uint8_t *p = passthrough.data() + i * passthrough_stride_;
      point_->extended_return_number        = p[0] & 0x0Fu;
      point_->extended_number_of_returns    = (p[0] >> 4) & 0x0Fu;
      point_->extended_classification_flags = p[1] & 0x0Fu;
      point_->extended_scanner_channel      = (p[1] >> 4) & 0x03u;
      point_->scan_direction_flag           = (p[1] >> 6) & 0x1u;
      point_->edge_of_flight_line           = (p[1] >> 7) & 0x1u;
      point_->extended_classification       = p[2];
      point_->user_data                     = p[3];
      int16_t ext_angle;
      std::memcpy(&ext_angle, p + 4, 2);
      point_->extended_scan_angle  = ext_angle;
      point_->point_source_ID      = static_cast<laszip_U16>(p[6]) | (static_cast<laszip_U16>(p[7]) << 8);
      // Restore original intensity from p[8..9], overriding the alpha default set above.
      std::memcpy(&point_->intensity, p + 8, 2);
      // Original sensor extra bytes at p[10..].
      if (orig_extra_size_ > 0)
      {
        uint16_t orig_start = 14; // sx+sy+sz+alpha+bound
        if (with_tree_id_) orig_start += 4;
        if (with_stem_id_) orig_start += 4;
        if (with_beam_id_) orig_start += 4;
        std::memcpy(point_->extra_bytes + orig_start, p + 10, orig_extra_size_);
      }
    }
    // Store start - end as three float32 extra bytes so starts can be reconstructed.
    const float sx = static_cast<float>(starts[i][0] - ends[i][0]);
    const float sy = static_cast<float>(starts[i][1] - ends[i][1]);
    const float sz = static_cast<float>(starts[i][2] - ends[i][2]);
    std::memcpy(point_->extra_bytes,     &sx, 4);
    std::memcpy(point_->extra_bytes + 4, &sy, 4);
    std::memcpy(point_->extra_bytes + 8, &sz, 4);
    uint16_t off = 12;
    if (with_tree_id_) {
      const int32_t tid = (i < tree_ids.size()) ? tree_ids[i] : -1;
      std::memcpy(point_->extra_bytes + off, &tid, 4);
      off += 4;
    }
    if (with_stem_id_) {
      const int32_t sid = (i < stem_ids.size()) ? stem_ids[i] : -1;
      std::memcpy(point_->extra_bytes + off, &sid, 4);
      off += 4;
    }
    if (with_beam_id_) {
      const int32_t bid = (i < beam_ids.size()) ? beam_ids[i] : -1;
      std::memcpy(point_->extra_bytes + off, &bid, 4);
      off += 4;
    }
    point_->extra_bytes[off] = colours[i].alpha;
    point_->extra_bytes[off + 1] = (colours[i].alpha > 0) ? 1 : 0;  // bound: 1=bound, 0=unbound
    laszip_write_point(writer_handle_);
  }
  points_written_ += ends.size();
  return true;
#else   // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(starts);
  RAYLIB_UNUSED(ends);
  RAYLIB_UNUSED(times);
  RAYLIB_UNUSED(colours);
  RAYLIB_UNUSED(tree_ids);
  RAYLIB_UNUSED(stem_ids);
  RAYLIB_UNUSED(passthrough);
  RAYLIB_UNUSED(beam_ids);
  std::cerr << "LasRayCloudWriter: WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

}  // namespace ray
