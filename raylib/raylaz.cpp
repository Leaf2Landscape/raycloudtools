// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylaz.h"
#include <algorithm>
#include <fstream>
#include <limits>
#include "raylib/rayprogress.h"
#include "raylib/rayprogressthread.h"
#include "rayunused.h"

#if RAYLIB_WITH_LAS
#include <laszip/laszip_api.h>
#endif  // RAYLIB_WITH_LAS

namespace ray
{
bool readLas(const std::string &file_name,
             std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                std::vector<double> &times, std::vector<RGBA> &colours)>
               apply,
             size_t &num_bounded, double max_intensity, Eigen::Vector3d *offset_to_remove, size_t chunk_size,
             std::vector<int32_t> *tree_ids_out, std::vector<uint8_t> *passthrough_out,
             uint16_t *orig_extra_size_out, std::vector<uint8_t> *extra_bytes_vlr_out)
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

  // LAS EXTRA_BYTES data_type → per-point byte size (types 0 and >10 are skipped)
  static const uint16_t kExtraTypeSize[11] = { 0, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };
  // Names of raycloud-owned extra attributes (these are skipped when extracting original data)
  static const char *kRayCloudAttrs[] = { "sx", "sy", "sz", "alpha", "tree_id" };

  uint16_t local_skip_size = 0;   // bytes of our own extra attributes before original data
  uint16_t local_orig_extra = 0;  // bytes of original sensor data per point
  bool has_tree_id_attr = false;  // true only if VLR explicitly declares "tree_id"
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
        if (strcmp(attr_name, "tree_id") == 0)
          has_tree_id_attr = true;
      }
      if (is_ours)
        local_skip_size += attr_size;
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
  for (size_t i = 0; i < number_of_points; i++)
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

    ends.push_back(position);

    if (is_raycloud && point->num_extra_bytes >= 12)
    {
      // Reconstruct ray start from the stored (start - end) float32 offset.
      float sx, sy, sz;
      memcpy(&sx, point->extra_bytes, 4);
      memcpy(&sy, point->extra_bytes + 4, 4);
      memcpy(&sz, point->extra_bytes + 8, 4);
      starts.push_back({ position[0] + sx, position[1] + sy, position[2] + sz });
      if (tree_ids_out && has_tree_id_attr)
      {
        int32_t tid;
        memcpy(&tid, point->extra_bytes + 12, 4);
        tree_ids_out->push_back(tid);
      }
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
      if (local_orig_extra > 0 && point->num_extra_bytes >= local_skip_size + local_orig_extra)
        passthrough_out->insert(passthrough_out->end(),
                                point->extra_bytes + local_skip_size,
                                point->extra_bytes + local_skip_size + local_orig_extra);
      else if (local_orig_extra > 0)
        passthrough_out->insert(passthrough_out->end(), local_orig_extra, 0);
    }

    if (using_colour)
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
    if (is_raycloud)
    {
      // Alpha is stored in extra_bytes at position 12 (or 16 with tree_id).
      // Prefer extra_bytes so the intensity field is free to carry the original sensor value.
      const uint16_t alpha_pos = has_tree_id_attr ? 16u : 12u;
      intensity = (point->num_extra_bytes > alpha_pos)
                    ? point->extra_bytes[alpha_pos]
                    : static_cast<uint8_t>(point->intensity);  // fallback for old files
    }
    else
    {
      const double normalised = (max_intensity > 0) ? (255.0 * point->intensity) / max_intensity : 255.0;
      intensity = static_cast<uint8_t>(std::min(normalised, 255.0));
      // Ensure any non-zero raw intensity maps to at least alpha=1 (bounded ray).
      if (intensity == 0 && point->intensity > 0)
        intensity = 1;
    }
    if (intensity > 0)
      num_bounded++;
    intensities.push_back(intensity);

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
  std::cerr << "readLas: cannot read file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

bool readLasExtraBytesVlr(const std::string &file_name, uint16_t &orig_extra_size_out,
                           std::vector<uint8_t> &extra_bytes_vlr_out)
{
#if RAYLIB_WITH_LAS
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

  static const uint16_t kExtraTypeSize[11] = { 0, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };
  static const char *kRayCloudAttrs[] = { "sx", "sy", "sz", "alpha", "tree_id" };

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
                                    const std::vector<uint8_t> &passthrough,
                                    const std::vector<uint8_t> &extra_bytes_vlr)
{
#if RAYLIB_WITH_LAS
  LasRayCloudWriter writer(file_name, !tree_ids.empty(), extra_bytes_vlr);
  return writer.writeChunk(starts, ends, times, colours, tree_ids, passthrough);
#else   // RAYLIB_WITH_LAS
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(starts);
  RAYLIB_UNUSED(ends);
  RAYLIB_UNUSED(times);
  RAYLIB_UNUSED(colours);
  std::cerr << "writeLasRayCloud: WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

#if RAYLIB_WITH_LAS
LasRayCloudWriter::LasRayCloudWriter(const std::string &file_name, bool with_tree_id,
                                     const std::vector<uint8_t> &extra_bytes_vlr)
  : file_name_(file_name)
  , points_written_(0)
  , with_tree_id_(with_tree_id)
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
    laszip_add_attribute(writer_handle_, 8, "sz", "ray start z offset", 1.0, 0.0) ||
    (with_tree_id_ && laszip_add_attribute(writer_handle_, 5, "tree_id", "per-point tree ID", 1.0, 0.0)) ||
    laszip_add_attribute(writer_handle_, 0, "alpha", "intensity 1-255", 1.0, 0.0);
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
  const laszip_U16 record_size =
    static_cast<laszip_U16>(36 + 12 + (with_tree_id_ ? 4 : 0) + 1 + orig_extra_size_);
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
LasRayCloudWriter::LasRayCloudWriter(const std::string &file_name, bool with_tree_id,
                                     const std::vector<uint8_t> &extra_bytes_vlr)
  : file_name_(file_name)
  , with_tree_id_(with_tree_id)
{
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(with_tree_id);
  RAYLIB_UNUSED(extra_bytes_vlr);
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
                                   const std::vector<uint8_t> &passthrough)
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
        const uint16_t orig_start = static_cast<uint16_t>(with_tree_id_ ? 17 : 13);
        std::memcpy(point_->extra_bytes + orig_start, p + 10, orig_extra_size_);
      }
    }
    // Store start - end as three float32 extra bytes so starts can be reconstructed.
    const float sx = static_cast<float>(starts[i][0] - ends[i][0]);
    const float sy = static_cast<float>(starts[i][1] - ends[i][1]);
    const float sz = static_cast<float>(starts[i][2] - ends[i][2]);
    std::memcpy(point_->extra_bytes, &sx, 4);
    std::memcpy(point_->extra_bytes + 4, &sy, 4);
    std::memcpy(point_->extra_bytes + 8, &sz, 4);
    if (with_tree_id_)
    {
      const int32_t tid = (i < tree_ids.size()) ? tree_ids[i] : -1;
      std::memcpy(point_->extra_bytes + 12, &tid, 4);
    }
    point_->extra_bytes[with_tree_id_ ? 16 : 12] = colours[i].alpha;
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
  RAYLIB_UNUSED(passthrough);
  std::cerr << "LasRayCloudWriter: WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
}

}  // namespace ray
