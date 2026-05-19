// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raysplitter.h"
#include <iostream>
#include <limits>
#include <map>
#include "extraction/rayforest.h"
#include "raycloudwriter.h"
#include "raycuboid.h"
#include "extraction/raytrees.h"
#include "rayparse.h"
#include "raylaz.h"

namespace ray
{
/// This is a helper function to aid in splitting the cloud while chunk-loading it. The purpose is to be able to
/// split clouds of any size, without running out of main memory.
bool split(const std::string &file_name, const std::string &in_name, const std::string &out_name,
           std::function<bool(const Cloud &cloud, int i)> is_outside)
{
  const std::string ext = getFileNameExtension(file_name);
  const bool is_las = (ext == "las" || ext == "laz");

  std::vector<uint8_t> extra_bytes_vlr;
  if (is_las)
  {
    uint16_t orig_extra = 0;
    readLasExtraBytesVlr(file_name, orig_extra, extra_bytes_vlr);
  }

  Cloud cloud_buffer;
  CloudWriter in_writer, out_writer;
  if (!in_writer.begin(in_name, extra_bytes_vlr))
    return false;
  if (!out_writer.begin(out_name, extra_bytes_vlr))
    return false;
  Cloud in_chunk, out_chunk;

  std::vector<uint8_t> passthrough_buf;

  // Move each ray into either the in_chunk or out_chunk, depending on the condition is_outside.
  // Passthrough (extra sensor bytes) is routed per-point via Cloud::addRay(const Cloud&, size_t).
  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                       std::vector<double> &times, std::vector<RGBA> &colours) {
    cloud_buffer.starts  = starts;
    cloud_buffer.ends    = ends;
    cloud_buffer.times   = times;
    cloud_buffer.colours = colours;
    if (!passthrough_buf.empty())
    {
      cloud_buffer.passthrough   = std::move(passthrough_buf);
      passthrough_buf.clear();
      cloud_buffer.extra_bytes_size = static_cast<uint16_t>(cloud_buffer.passthrough.size() / starts.size());
      cloud_buffer.extra_bytes_vlr  = extra_bytes_vlr;
    }

    for (int i = 0; i < (int)cloud_buffer.ends.size(); i++)
    {
      Cloud &chunk = is_outside(cloud_buffer, i) ? out_chunk : in_chunk;
      chunk.addRay(cloud_buffer, static_cast<size_t>(i));
    }
    in_writer.writeChunk(in_chunk);
    out_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
    cloud_buffer.passthrough.clear();
  };

  bool res;
  if (is_las)
  {
    size_t num_bounded;
    res = readLas(file_name, per_chunk, num_bounded, 1.0, nullptr, 1000000, nullptr, &passthrough_buf);
  }
  else
  {
    res = Cloud::read(file_name, per_chunk);
  }
  if (!res)
    return false;
  in_writer.end();
  out_writer.end();
  return true;
}

/// Special case for splitting around a plane.
bool splitPlane(const std::string &file_name, const std::string &in_name, const std::string &out_name,
                const Eigen::Vector3d &plane)
{
  const std::string ext = getFileNameExtension(file_name);
  const bool is_las = (ext == "las" || ext == "laz");

  std::vector<uint8_t> extra_bytes_vlr;
  if (is_las)
  {
    uint16_t orig_extra = 0;
    readLasExtraBytesVlr(file_name, orig_extra, extra_bytes_vlr);
  }

  CloudWriter inside_writer, outside_writer;
  if (!inside_writer.begin(in_name, extra_bytes_vlr))
    return false;
  if (!outside_writer.begin(out_name, extra_bytes_vlr))
    return false;
  Cloud in_chunk, out_chunk;
  Cloud cloud_buffer;
  std::vector<uint8_t> passthrough_buf;

  auto copy_pass = [&](Cloud &dst, size_t src_i) {
    const uint16_t stride = cloud_buffer.extra_bytes_size;
    if (stride > 0 && cloud_buffer.passthrough.size() >= (src_i + 1) * stride)
    {
      const uint8_t *p = cloud_buffer.passthrough.data() + src_i * stride;
      dst.passthrough.insert(dst.passthrough.end(), p, p + stride);
      dst.extra_bytes_size = stride;
      if (dst.extra_bytes_vlr.empty())
        dst.extra_bytes_vlr = extra_bytes_vlr;
    }
  };

  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                       std::vector<double> &times, std::vector<RGBA> &colours) {
    cloud_buffer.starts  = starts;
    cloud_buffer.ends    = ends;
    cloud_buffer.times   = times;
    cloud_buffer.colours = colours;
    if (!passthrough_buf.empty())
    {
      cloud_buffer.passthrough      = std::move(passthrough_buf);
      passthrough_buf.clear();
      cloud_buffer.extra_bytes_size = static_cast<uint16_t>(cloud_buffer.passthrough.size() / starts.size());
      cloud_buffer.extra_bytes_vlr  = extra_bytes_vlr;
    }
    const Eigen::Vector3d plane_vec = plane / plane.dot(plane);
    for (size_t i = 0; i < ends.size(); i++)
    {
      const double d1 = starts[i].dot(plane_vec) - 1.0;
      const double d2 = ends[i].dot(plane_vec) - 1.0;
      if (d1 * d2 > 0.0)  // start and end are on the same side of the plane, so don't split...
      {
        Cloud &chunk = d1 > 0.0 ? out_chunk : in_chunk;
        chunk.addRay(cloud_buffer, i);
      }
      else  // split the ray...
      {
        RGBA col = colours[i];
        col.red = col.green = col.blue = col.alpha = 0;
        const Eigen::Vector3d mid = starts[i] + (ends[i] - starts[i]) * d1 / (d1 - d2);
        if (d1 > 0.0)
        {
          out_chunk.addRay(starts[i], mid, times[i], col);
          copy_pass(out_chunk, i);
          in_chunk.addRay(mid, ends[i], times[i], colours[i]);
          copy_pass(in_chunk, i);
        }
        else
        {
          in_chunk.addRay(starts[i], mid, times[i], col);
          copy_pass(in_chunk, i);
          out_chunk.addRay(mid, ends[i], times[i], colours[i]);
          copy_pass(out_chunk, i);
        }
      }
    }
    inside_writer.writeChunk(in_chunk);
    outside_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
    cloud_buffer.passthrough.clear();
  };

  bool res;
  if (is_las)
  {
    size_t num_bounded;
    res = readLas(file_name, per_chunk, num_bounded, 1.0, nullptr, 1000000, nullptr, &passthrough_buf);
  }
  else
  {
    res = Cloud::read(file_name, per_chunk);
  }
  if (!res)
    return false;
  inside_writer.end();
  outside_writer.end();
  return true;
}

/// Special case for splitting a capsule.
bool splitCapsule(const std::string &file_name, const std::string &in_name, const std::string &out_name,
                  const Eigen::Vector3d &end1, const Eigen::Vector3d &end2, double radius)
{
  const std::string ext = getFileNameExtension(file_name);
  const bool is_las = (ext == "las" || ext == "laz");

  std::vector<uint8_t> extra_bytes_vlr;
  if (is_las)
  {
    uint16_t orig_extra = 0;
    readLasExtraBytesVlr(file_name, orig_extra, extra_bytes_vlr);
  }

  CloudWriter inside_writer, outside_writer;
  if (!inside_writer.begin(in_name, extra_bytes_vlr))
    return false;
  if (!outside_writer.begin(out_name, extra_bytes_vlr))
    return false;
  Cloud in_chunk, out_chunk;
  Cloud cloud_buffer;
  std::vector<uint8_t> passthrough_buf;

  Eigen::Vector3d dir = end2 - end1;
  double length = dir.norm();
  if (length > 0.0)
  {
    dir /= length;
  }

  auto copy_pass = [&](Cloud &dst, size_t src_i) {
    const uint16_t stride = cloud_buffer.extra_bytes_size;
    if (stride > 0 && cloud_buffer.passthrough.size() >= (src_i + 1) * stride)
    {
      const uint8_t *p = cloud_buffer.passthrough.data() + src_i * stride;
      dst.passthrough.insert(dst.passthrough.end(), p, p + stride);
      dst.extra_bytes_size = stride;
      if (dst.extra_bytes_vlr.empty())
        dst.extra_bytes_vlr = extra_bytes_vlr;
    }
  };

  // splitting per chunk
  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                       std::vector<double> &times, std::vector<RGBA> &colours) {
    cloud_buffer.starts  = starts;
    cloud_buffer.ends    = ends;
    cloud_buffer.times   = times;
    cloud_buffer.colours = colours;
    if (!passthrough_buf.empty())
    {
      cloud_buffer.passthrough      = std::move(passthrough_buf);
      passthrough_buf.clear();
      cloud_buffer.extra_bytes_size = static_cast<uint16_t>(cloud_buffer.passthrough.size() / starts.size());
      cloud_buffer.extra_bytes_vlr  = extra_bytes_vlr;
    }
    for (size_t i = 0; i < ends.size(); i++)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      Eigen::Vector3d ray = end - start;
      // The approach is to find the d value (ratio along ray) for the first and second intersection with 
      // the capsule. This can be found by breaking it into a cylinder and two spheres, and
      // a few min/maxs.

      // cylinder part:
      double cylinder_intersection1 = 1e10;
      double cylinder_intersection2 = -1e10;
      Eigen::Vector3d up = dir.cross(ray);
      double mag = up.norm();
      if (mag > 0.0) // two rays are not inline
      {
        up /= mag;
        double gap = std::abs((start - end1).dot(up));
        if (gap >= radius)
        {
          out_chunk.addRay(start, end, times[i], colours[i]);
          copy_pass(out_chunk, i);
          continue; // if it doesn't hit the endless cylinder it won't hit the capsule
        }
        Eigen::Vector3d lateral_dir = ray - dir * ray.dot(dir);
        double lateral_length = lateral_dir.norm();
        double d_mid = (end1 - start).dot(lateral_dir) / ray.dot(lateral_dir);
        
        double shift = std::sqrt(radius*radius - gap*gap) / lateral_length;
        double d_min = d_mid - shift;
        double d_max = d_mid + shift;
        double d1 = (start + ray*d_min - end1).dot(dir) / length;
        double d2 = (start + ray*d_max - end1).dot(dir) / length;
        if (d1 > 0.0 && d1 < 1.0)
        {
          cylinder_intersection1 = d_min;
        }
        if (d2 > 0.0 && d2 < 1.0)
        {
          cylinder_intersection2 = d_max;
        }
      }

      // the spheres part:
      double ray_length = ray.norm();
      double sphere_intersection1[2] = {1e10, 1e10};
      double sphere_intersection2[2] = {-1e10, -1e10};
      Eigen::Vector3d ends[2] = {end1, end2};
      for (int e = 0; e<2; e++)
      {
        double mid_d = (ends[e] - start).dot(ray) / (ray_length * ray_length);
        Eigen::Vector3d shortest_dir = (ends[e] - start) - ray * mid_d;
        double shortest_sqr = shortest_dir.squaredNorm();
        if (shortest_sqr < radius*radius)
        {
          double shift = std::sqrt(radius * radius - shortest_sqr) / ray_length;
          sphere_intersection1[e] = mid_d - shift;
          sphere_intersection2[e] = mid_d + shift;
        }
      }

      // combining together
      double closest_d = std::min( {cylinder_intersection1, sphere_intersection1[0], sphere_intersection1[1]} );
      double farthest_d = std::max( {cylinder_intersection2, sphere_intersection2[0], sphere_intersection2[1]} );

      RGBA black;
      black.red = black.green = black.blue = black.alpha = 0;
      // easy case
      if (closest_d >= 1.0 || farthest_d <= 0.0)
      {
        out_chunk.addRay(start, end, times[i], colours[i]);
        copy_pass(out_chunk, i);
        continue;
      }
      // first outside ray
      if (closest_d > 0.0)
      {
        out_chunk.addRay(start, start + ray*closest_d, times[i], black);
        copy_pass(out_chunk, i);
      }
      // second outside ray
      if (farthest_d < 1.0)
      {
        out_chunk.addRay(start + ray * farthest_d, end, times[i], colours[i]);
        copy_pass(out_chunk, i);
      }
      // inside ray
      if (farthest_d < 1.0)
      {
        in_chunk.addRay(start + ray * std::max(0.0, closest_d), start + ray*std::min(farthest_d, 1.0), times[i], black);
        copy_pass(in_chunk, i);
      }
      else
      {
        in_chunk.addRay(start + ray * std::max(0.0, closest_d), end, times[i], colours[i]);
        copy_pass(in_chunk, i);
      }
    }
    inside_writer.writeChunk(in_chunk);
    outside_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
    cloud_buffer.passthrough.clear();
  };

  bool res;
  if (is_las)
  {
    size_t num_bounded;
    res = readLas(file_name, per_chunk, num_bounded, 1.0, nullptr, 1000000, nullptr, &passthrough_buf);
  }
  else
  {
    res = Cloud::read(file_name, per_chunk);
  }
  if (!res)
    return false;
  inside_writer.end();
  outside_writer.end();
  return true;
}


/// Special case for splitting a box.
bool splitBox(const std::string &file_name, const std::string &in_name, const std::string &out_name,
              const Eigen::Vector3d &centre, const Eigen::Vector3d &extents)
{
  const std::string ext = getFileNameExtension(file_name);
  const bool is_las = (ext == "las" || ext == "laz");

  std::vector<uint8_t> extra_bytes_vlr;
  if (is_las)
  {
    uint16_t orig_extra = 0;
    readLasExtraBytesVlr(file_name, orig_extra, extra_bytes_vlr);
  }

  CloudWriter inside_writer, outside_writer;
  if (!inside_writer.begin(in_name, extra_bytes_vlr))
    return false;
  if (!outside_writer.begin(out_name, extra_bytes_vlr))
    return false;
  Cloud in_chunk, out_chunk;
  Cloud cloud_buffer;
  std::vector<uint8_t> passthrough_buf;

  auto copy_pass = [&](Cloud &dst, size_t src_i) {
    const uint16_t stride = cloud_buffer.extra_bytes_size;
    if (stride > 0 && cloud_buffer.passthrough.size() >= (src_i + 1) * stride)
    {
      const uint8_t *p = cloud_buffer.passthrough.data() + src_i * stride;
      dst.passthrough.insert(dst.passthrough.end(), p, p + stride);
      dst.extra_bytes_size = stride;
      if (dst.extra_bytes_vlr.empty())
        dst.extra_bytes_vlr = extra_bytes_vlr;
    }
  };

  // splitting per chunk
  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                       std::vector<double> &times, std::vector<RGBA> &colours) {
    cloud_buffer.starts  = starts;
    cloud_buffer.ends    = ends;
    cloud_buffer.times   = times;
    cloud_buffer.colours = colours;
    if (!passthrough_buf.empty())
    {
      cloud_buffer.passthrough      = std::move(passthrough_buf);
      passthrough_buf.clear();
      cloud_buffer.extra_bytes_size = static_cast<uint16_t>(cloud_buffer.passthrough.size() / starts.size());
      cloud_buffer.extra_bytes_vlr  = extra_bytes_vlr;
    }
    const Cuboid cuboid(centre - extents, centre + extents);
    for (size_t i = 0; i < ends.size(); i++)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      if (cuboid.clipRay(start, end))  // true if ray intersects the cuboid
      {
        RGBA col = colours[i];
        if (!cuboid.intersects(ends[i]))  // mark as unbounded for the in_chunk
        {
          col.red = col.green = col.blue = col.alpha = 0;
        }
        in_chunk.addRay(start, end, times[i], col);
        copy_pass(in_chunk, i);
        if (start != starts[i])  // start part is clipped
        {
          col.red = col.green = col.blue = col.alpha = 0;
          out_chunk.addRay(starts[i], start, times[i], col);
          copy_pass(out_chunk, i);
        }
        if (ends[i] != end)  // end part is clipped
        {
          out_chunk.addRay(end, ends[i], times[i], colours[i]);
          copy_pass(out_chunk, i);
        }
      }
      else  // no intersection
      {
        out_chunk.addRay(starts[i], ends[i], times[i], colours[i]);
        copy_pass(out_chunk, i);
      }
    }
    inside_writer.writeChunk(in_chunk);
    outside_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
    cloud_buffer.passthrough.clear();
  };

  bool res;
  if (is_las)
  {
    size_t num_bounded;
    res = readLas(file_name, per_chunk, num_bounded, 1.0, nullptr, 1000000, nullptr, &passthrough_buf);
  }
  else
  {
    res = Cloud::read(file_name, per_chunk);
  }
  if (!res)
    return false;
  inside_writer.end();
  outside_writer.end();
  return true;
}

/// Special case for splitting based on a grid.
bool splitGrid(const std::string &file_name, const std::string &cloud_name_stub, const Eigen::Vector3d &cell_width,
               double overlap)
{
  return splitGrid(file_name, cloud_name_stub, Eigen::Vector4d(cell_width[0], cell_width[1], cell_width[2], 0),
                   overlap);
}

/// Special case for splitting based on a grid.
bool splitGrid(const std::string &file_name, const std::string &cloud_name_stub, const Eigen::Vector4d &cell_width,
               double overlap)
{
  const std::string ext = getFileNameExtension(file_name);
  const bool is_las = (ext == "las" || ext == "laz");

  std::vector<uint8_t> extra_bytes_vlr;
  if (is_las)
  {
    uint16_t orig_extra = 0;
    readLasExtraBytesVlr(file_name, orig_extra, extra_bytes_vlr);
  }

  overlap /= 2.0;  // it now means overlap relative to grid edge
  Cloud::Info info;
  Cloud::getInfo(file_name, info);
  const Eigen::Vector3d &min_bound = info.rays_bound.min_bound_;
  const Eigen::Vector3d &max_bound = info.rays_bound.max_bound_;

  Eigen::Vector4d width = cell_width;
  for (int i = 0; i < 4; i++)
  {
    if (width[i] ==
        0.0)  // '0' is a convenience for the highest possible cell width, meaning that this axis isn't split
      width[i] = std::numeric_limits<double>::max();
  }

  const Eigen::Vector3d minID(std::floor(0.5 + min_bound[0] / width[0]), std::floor(0.5 + min_bound[1] / width[1]),
                              std::floor(0.5 + min_bound[2] / width[2]));
  const Eigen::Vector3d maxID(std::ceil(0.5 + max_bound[0] / width[0]), std::ceil(0.5 + max_bound[1] / width[1]),
                              std::ceil(0.5 + max_bound[2] / width[2]));
  const Eigen::Vector3i min_index = minID.cast<int>();
  const Eigen::Vector3i max_index = maxID.cast<int>();
  const Eigen::Vector3i dimensions = max_index - min_index;

  // I do something different for time, because the absolute values are so large that integer overflow is possible (e.g.
  // gridding every 10th of a second on a v short scan) A double can store consecutive integers up to 9e15, compared to
  // an int which is 2e9, so I cast to a long int to get the larger range
  const long int min_time = static_cast<long int>(std::floor(0.5 + info.min_time / width[3]));
  const long int max_time = static_cast<long int>(std::ceil(0.5 + info.max_time / width[3]));
  const int time_dimension = static_cast<int>(
    max_time - min_time);  // the difference won't overflow integers. We don't scan for 20 years straight.

  const int length = dimensions[0] * dimensions[1] * dimensions[2] * time_dimension;
  std::cout << "splitting into maximum of: " << length << " files" << std::endl;
  if (length > 50000)
  {
    std::cerr << "error: output of over 50,000 files is probably a mistake, exiting" << std::endl;
    return false;
  }
  const int max_open_files = 256;
  if (length > max_open_files)
  {
    std::cout << "Warning: nominally more than " << max_open_files << " file pointers will be open at once." << std::endl;
    std::cout << "Diving the operation into " << 1+length/max_open_files << " passes" << std::endl;
  }
  const std::string grid_ext = getFileNameExtension(file_name);
  for (int pass = 0; pass<length; pass+=max_open_files)
  {
    if (pass > 0)
    {
      std::cout << "Running pass " << 1 + pass/max_open_files << " / " << 1+length/max_open_files << std::endl;
    }
    std::vector<CloudWriter> cells(max_open_files);
    std::vector<Cloud> chunks(max_open_files);
    Cloud cloud_buffer;
    std::vector<uint8_t> passthrough_buf;

    // splitting performed per chunk
    auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts,
                         std::vector<Eigen::Vector3d> &ends, std::vector<double> &times,
                         std::vector<RGBA> &colours) {
      cloud_buffer.starts  = starts;
      cloud_buffer.ends    = ends;
      cloud_buffer.times   = times;
      cloud_buffer.colours = colours;
      if (!passthrough_buf.empty())
      {
        cloud_buffer.passthrough      = std::move(passthrough_buf);
        passthrough_buf.clear();
        cloud_buffer.extra_bytes_size = static_cast<uint16_t>(cloud_buffer.passthrough.size() / starts.size());
        cloud_buffer.extra_bytes_vlr  = extra_bytes_vlr;
      }
      for (size_t i = 0; i < ends.size(); i++)
      {
        // get set of cells that the ray may intersect
        const Eigen::Vector3d from(0.5 + starts[i][0] / width[0], 0.5 + starts[i][1] / width[1], 0.5 + starts[i][2] / width[2]);
        const Eigen::Vector3d to(0.5 + ends[i][0] / width[0], 0.5 + ends[i][1] / width[1], 0.5 + ends[i][2] / width[2]);
        const Eigen::Vector3d pos0 = minVector(from, to) - Eigen::Vector3d(overlap, overlap, 0.0);
        const Eigen::Vector3d pos1 = maxVector(from, to) + Eigen::Vector3d(overlap, overlap, 0.0);
        Eigen::Vector3i minI = Eigen::Vector3d(std::floor(pos0[0]), std::floor(pos0[1]), std::floor(pos0[2])).cast<int>();
        Eigen::Vector3i maxI = Eigen::Vector3d(std::ceil(pos1[0]), std::ceil(pos1[1]), std::ceil(pos1[2])).cast<int>();
        if (overlap > 0.0)
        {
          minI = maxVector(minI, min_index);
          maxI = minVector(maxI, max_index);
        }
        const long int t = static_cast<long int>(std::floor(0.5 + times[i] / width[3]));
        for (int x = minI[0]; x < maxI[0]; x++)
        {
          for (int y = minI[1]; y < maxI[1]; y++)
          {
            for (int z = minI[2]; z < maxI[2]; z++)
            {
              const int time_dif = static_cast<int>(t - min_time);
              int index = (x - min_index[0]) + dimensions[0] * (y - min_index[1]) +
                                dimensions[0] * dimensions[1] * (z - min_index[2]) +
                                dimensions[0] * dimensions[1] * dimensions[2] * time_dif;
              if (index < 0 || index >= length)
              {
                std::cout << "Error: bad index: " << index << std::endl;  // this should not happen
                return;
              }
              index -= pass;
              if (index < 0 || index >= max_open_files)
                continue;
              // do actual clipping here....
              const Eigen::Vector3d box_min(((double)x - 0.5) * width[0] - overlap,
                                            ((double)y - 0.5) * width[1] - overlap, ((double)z - 0.5) * width[2]);
              const Eigen::Vector3d box_max(((double)x + 0.5) * width[0] + overlap,
                                            ((double)y + 0.5) * width[1] + overlap, ((double)z + 0.5) * width[2]);
              const Cuboid cuboid(box_min, box_max);
              Eigen::Vector3d start = starts[i];
              Eigen::Vector3d end = ends[i];

              if (cuboid.clipRay(start, end))
              {
                RGBA col = colours[i];
                if (cells[index].fileName().empty())  // first time in this cell, so start writing to a new file
                {
                  std::stringstream name;
                  name << cloud_name_stub;
                  if (cell_width[0] > 0.0)
                    name << "_" << x;
                  if (cell_width[1] > 0.0)
                    name << "_" << y;
                  if (cell_width[2] > 0.0)
                    name << "_" << z;
                  if (cell_width[3] > 0.0)
                    name << "_" << t;
                  name << "." << grid_ext;
                  cells[index].begin(name.str(), extra_bytes_vlr);
                }
                if (!cuboid.intersects(ends[i]))  // end point is outside, so mark an unbounded ray
                {
                  col.red = col.green = col.blue = col.alpha = 0;
                }
                chunks[index].addRay(start, end, times[i], col);
                // copy passthrough bytes from original ray
                const uint16_t stride = cloud_buffer.extra_bytes_size;
                if (stride > 0 && cloud_buffer.passthrough.size() >= (i + 1) * stride)
                {
                  const uint8_t *p = cloud_buffer.passthrough.data() + i * stride;
                  chunks[index].passthrough.insert(chunks[index].passthrough.end(), p, p + stride);
                  chunks[index].extra_bytes_size = stride;
                  if (chunks[index].extra_bytes_vlr.empty())
                    chunks[index].extra_bytes_vlr = extra_bytes_vlr;
                }
              }
            }
          }
        }
      }
      for (int i = 0; i < max_open_files; i++)
      {
        if (chunks[i].ends.size() > 0)
        {
          cells[i].writeChunk(chunks[i]);
          chunks[i].clear();
        }
      }
      cloud_buffer.passthrough.clear();
    };

    bool res;
    if (is_las)
    {
      size_t num_bounded;
      res = readLas(file_name, per_chunk, num_bounded, 1.0, nullptr, 1000000, nullptr, &passthrough_buf);
    }
    else
    {
      res = Cloud::read(file_name, per_chunk);
    }
    if (!res)
      return false;

    for (int i = 0; i < max_open_files; i++)
    {
      cells[i].end();  // has no effect on writers where begin has not been called
    }
  }
  return true;
}

class RGBALess
{
public:
  bool operator()(const RGBA &a, const RGBA &b) const
  {
    if (a.red != b.red)
      return a.red < b.red;
    if (a.green != b.green)
      return a.green < b.green;
    return a.blue < b.blue;
  }
};

/// Special case for splitting based on a colour
bool splitColour(const std::string &file_name, const std::string &cloud_name_stub, bool seg_colour)
{
  const std::string file_ext = getFileNameExtension(file_name);
  const bool is_las = (file_ext == "las" || file_ext == "laz");

  std::vector<uint8_t> extra_bytes_vlr;
  if (is_las)
  {
    uint16_t orig_extra = 0;
    readLasExtraBytesVlr(file_name, orig_extra, extra_bytes_vlr);
  }

  std::map<RGBA, int, RGBALess> vox_map;
  // firstly, find out how many different colours there are
  int num_colours = 0;
  auto count_colours = [&vox_map, &num_colours](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &,
                                                std::vector<double> &, std::vector<ray::RGBA> &colours) {
    for (auto &colour : colours)
    {
      if (vox_map.find(colour) == vox_map.end())
      {
        vox_map.insert(std::pair<RGBA, int>(colour, num_colours++));
      }
    }
  };
  if (!ray::Cloud::read(file_name, count_colours))
    return false;

  const int max_total_files = 50000; // raysplit colour more likely to be a mistake in this case
  if (num_colours > max_total_files)
  {
    std::cerr << "Error: " << num_colours << " colours generates more than the maximum number of files: " << max_total_files << std::endl;
    return false;
  }
  const int max_files_at_once = 512;  // operating systems will fail with too many open file pointers.
  std::cout << "splitting into: " << num_colours << " files" << std::endl;
  if (num_colours > max_files_at_once)
  {
    std::cout << "Warning: cloud has more unique colours than allowed for simultaneous files " << max_files_at_once << " so using multiple passes." << std::endl;
  }

  int chunk_size = std::min(num_colours, max_files_at_once);
  const std::string colour_ext = getFileNameExtension(file_name);

  for (int batch = 0; batch < num_colours; batch+=max_files_at_once)
  {
    if (batch > 0)
    {
      std::cout << "Running pass " << 1+batch/max_files_at_once << " / " << 1 + num_colours / max_files_at_once << std::endl;
    }
    std::vector<CloudWriter> cells(chunk_size);
    std::vector<Cloud> chunks(chunk_size);
    Cloud cloud_buffer;
    std::vector<uint8_t> passthrough_buf;

    int batch_max = std::min(num_colours, batch + max_files_at_once);
    if (num_colours > max_files_at_once)
    {
      std::cout << "batch processing colours " << batch << ", to " << batch_max << std::endl;
    }
    // splitting performed per chunk
    auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                         std::vector<double> &times, std::vector<RGBA> &colours) {
      cloud_buffer.starts  = starts;
      cloud_buffer.ends    = ends;
      cloud_buffer.times   = times;
      cloud_buffer.colours = colours;
      if (!passthrough_buf.empty())
      {
        cloud_buffer.passthrough      = std::move(passthrough_buf);
        passthrough_buf.clear();
        cloud_buffer.extra_bytes_size = static_cast<uint16_t>(cloud_buffer.passthrough.size() / starts.size());
        cloud_buffer.extra_bytes_vlr  = extra_bytes_vlr;
      }
      for (size_t i = 0; i < ends.size(); i++)
      {
        RGBA colour = colours[i];
        const auto &vox = vox_map.find(colour);
        if (vox != vox_map.end())
        {
          int index = vox->second - batch;
          if (index < 0 || index >= max_files_at_once)
          {
            continue; // ignores as this ray will not go into this batch
          }

          if (cells[index].fileName().empty())  // first time in this cell, so start writing to a new file
          {
            std::stringstream name;
            if (seg_colour)
            {
              name << cloud_name_stub << "_" << convertColourToInt(colour) << "." << colour_ext;
            }
            else
            {
              name << cloud_name_stub << "_" << (int)colour.red << "_" << (int)colour.green << "_" << (int)colour.blue << "." << colour_ext;
            }
            cells[index].begin(name.str(), extra_bytes_vlr);
          }
          chunks[index].addRay(cloud_buffer, i);
        }
      }
      for (int i = 0; i < chunk_size; i++)
      {
        if (chunks[i].ends.size() > 0)
        {
          cells[i].writeChunk(chunks[i]);
          chunks[i].clear();
        }
      }
      cloud_buffer.passthrough.clear();
    };

    bool res;
    if (is_las)
    {
      size_t num_bounded;
      res = readLas(file_name, per_chunk, num_bounded, 1.0, nullptr, 1000000, nullptr, &passthrough_buf);
    }
    else
    {
      res = Cloud::read(file_name, per_chunk);
    }
    if (!res)
      return false;

    for (int i = 0; i < chunk_size; i++)
    {
      cells[i].end();  // has no effect on writers where begin has not been called
    }
  }
  return true;
}

}  // namespace ray
