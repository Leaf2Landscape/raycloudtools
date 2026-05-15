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
             std::vector<int32_t> *tree_ids_out, std::vector<uint8_t> *passthrough_out)
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
      if (tree_ids_out && point->num_extra_bytes >= 16)
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

    // Pack 8 bytes of LAS 1.4 extended fields per point into passthrough.
    // Layout: [0] ext_return[0:3]|ext_num_returns[4:7]
    //         [1] class_flags[0:3]|scanner_chan[4:5]|scan_dir[6]|edge[7]
    //         [2] extended_classification
    //         [3] user_data
    //         [4-5] extended_scan_angle (int16 LE, 0.006 deg units)
    //         [6-7] point_source_ID (uint16 LE)
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
      // Alpha is stored directly in the intensity field (0-255).
      intensity = static_cast<uint8_t>(point->intensity);
    }
    else
    {
      const double normalised = (max_intensity > 0) ? (255.0 * point->intensity) / max_intensity : 255.0;
      intensity = static_cast<uint8_t>(std::min(normalised, 255.0));
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
  RAYLIB_UNUSED(max_intensity);
  RAYLIB_UNUSED(tree_ids_out);
  RAYLIB_UNUSED(passthrough_out);
  std::cerr << "readLas: cannot read file as WITHLAS not enabled. Enable using: cmake .. -DWITH_LAS=true" << std::endl;
  return false;
#endif  // RAYLIB_WITH_LAS
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
  header->point_data_format = 6;  // LAS 1.4: GPS time only
  const double scale = 1e-4;
  header->x_scale_factor = scale;
  header->y_scale_factor = scale;
  header->z_scale_factor = scale;
  header->x_offset = 0.0;
  header->y_offset = 0.0;
  header->z_offset = 0.0;
  header->extended_number_of_point_records = static_cast<laszip_U64>(points.size());

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
  header_->point_data_format = 6;  // LAS 1.4: GPS time only
  const double scale = 1e-4;
  header_->x_scale_factor = scale;
  header_->y_scale_factor = scale;
  header_->z_scale_factor = scale;
  header_->x_offset = 0.0;
  header_->y_offset = 0.0;
  header_->z_offset = 0.0;

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
                                    const std::vector<uint8_t> &passthrough)
{
#if RAYLIB_WITH_LAS
  LasRayCloudWriter writer(file_name, !tree_ids.empty());
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
LasRayCloudWriter::LasRayCloudWriter(const std::string &file_name, bool with_tree_id)
  : file_name_(file_name)
  , points_written_(0)
  , with_tree_id_(with_tree_id)
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
  header->point_data_format = 7;  // LAS 1.4: GPS time + RGB (equivalent of format 3)
  const double scale = 1e-4;
  header->x_scale_factor = scale;
  header->y_scale_factor = scale;
  header->z_scale_factor = scale;
  header->x_offset = 0.0;
  header->y_offset = 0.0;
  header->z_offset = 0.0;

  // Add three float32 extra attributes for ray start offset (start - end).
  // LASzip API type 8 = F32. Scale/offset unused for floating-point types.
  // Optionally add a fourth int32 attribute for the per-point tree ID.
  // LASzip API type 5 = INT32.
  bool attr_err =
    laszip_add_attribute(writer_handle_, 8, "sx", "ray start x offset", 1.0, 0.0) ||
    laszip_add_attribute(writer_handle_, 8, "sy", "ray start y offset", 1.0, 0.0) ||
    laszip_add_attribute(writer_handle_, 8, "sz", "ray start z offset", 1.0, 0.0) ||
    (with_tree_id_ && laszip_add_attribute(writer_handle_, 5, "tree_id", "per-point tree ID", 1.0, 0.0));
  if (attr_err)
  {
    laszip_CHAR *error;
    laszip_get_error(writer_handle_, &error);
    std::cerr << "LasRayCloudWriter: failed to add extra attributes: " << error << std::endl;
    laszip_destroy(writer_handle_);
    writer_handle_ = nullptr;
    return;
  }

  // Explicitly set point type 7 (LAS 1.4 GPS time + RGB) and total record size.
  // LAS 1.4 format 7 base = 36 bytes. laszip_add_attribute uses a default base, so override here.
  const laszip_U16 record_size = static_cast<laszip_U16>(36 + 12 + (with_tree_id_ ? 4 : 0));
  if (laszip_set_point_type_and_size(writer_handle_, 7, record_size))
  {
    laszip_CHAR *error;
    laszip_get_error(writer_handle_, &error);
    std::cerr << "LasRayCloudWriter: failed to set point type/size: " << error << std::endl;
    laszip_destroy(writer_handle_);
    writer_handle_ = nullptr;
    return;
  }

  // Custom VLR marks the file as a ray cloud so readers can detect it.
  if (laszip_add_vlr(writer_handle_, "raycloudtools", 1, 0, "raycloud", nullptr))
  {
    laszip_CHAR *error;
    laszip_get_error(writer_handle_, &error);
    std::cerr << "LasRayCloudWriter: failed to add VLR: " << error << std::endl;
    laszip_destroy(writer_handle_);
    writer_handle_ = nullptr;
    return;
  }

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
LasRayCloudWriter::LasRayCloudWriter(const std::string &file_name, bool with_tree_id)
  : file_name_(file_name)
  , with_tree_id_(with_tree_id)
{
  RAYLIB_UNUSED(file_name);
  RAYLIB_UNUSED(with_tree_id);
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
    // Restore LAS 1.4 extended fields from 8-byte passthrough if available.
    // Layout mirrors readLas passthrough packing (see readLas for byte definitions).
    if (passthrough.size() >= (i + 1) * 8)
    {
      const uint8_t *p = passthrough.data() + i * 8;
      point_->extended_return_number       = p[0] & 0x0Fu;
      point_->extended_number_of_returns   = (p[0] >> 4) & 0x0Fu;
      point_->extended_classification_flags= p[1] & 0x0Fu;
      point_->extended_scanner_channel     = (p[1] >> 4) & 0x03u;
      point_->scan_direction_flag          = (p[1] >> 6) & 0x1u;
      point_->edge_of_flight_line          = (p[1] >> 7) & 0x1u;
      point_->extended_classification      = p[2];
      point_->user_data                    = p[3];
      int16_t ext_angle;
      memcpy(&ext_angle, p + 4, 2);
      point_->extended_scan_angle          = ext_angle;
      point_->point_source_ID              = static_cast<laszip_U16>(p[6]) | (static_cast<laszip_U16>(p[7]) << 8);
    }
    // Store start - end as three float32 extra bytes so starts can be reconstructed.
    const float sx = static_cast<float>(starts[i][0] - ends[i][0]);
    const float sy = static_cast<float>(starts[i][1] - ends[i][1]);
    const float sz = static_cast<float>(starts[i][2] - ends[i][2]);
    memcpy(point_->extra_bytes, &sx, 4);
    memcpy(point_->extra_bytes + 4, &sy, 4);
    memcpy(point_->extra_bytes + 8, &sz, 4);
    if (with_tree_id_)
    {
      const int32_t tid = (i < tree_ids.size()) ? tree_ids[i] : -1;
      memcpy(point_->extra_bytes + 12, &tid, 4);
    }
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
