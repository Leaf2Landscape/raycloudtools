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
#include "raylaz.h" // For LasWriter

namespace ray
{
/// This is a helper function to aid in splitting the cloud while chunk-loading it. The purpose is to be able to
/// split clouds of any size, without running out of main memory.
bool split(const std::string &file_name, const std::string &in_name, const std::string &out_name,
           std::function<bool(const Cloud &cloud, int i)> is_outside,
           const std::string& output_ext)
{
  Cloud cloud_buffer;
  CloudWriter in_writer, out_writer;
  if (!in_writer.begin(in_name, output_ext)) return false;
  if (!out_writer.begin(out_name, output_ext)) return false;
  Cloud in_chunk, out_chunk;

  /// move each ray into either the in_chunk or out_chunk, depending on the condition function is_outside
  auto per_chunk = [&](
                     std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                     std::vector<double> &times, std::vector<RGBA> &colours,
                     std::vector<uint8_t> &classifications, std::vector<uint16_t> &branch_ids) {
    // I move these into the cloud buffer, so that they can be indexed easily in is_outside (by index).
    cloud_buffer.starts = std::move(starts);
    cloud_buffer.ends = std::move(ends);
    cloud_buffer.times = std::move(times);
    cloud_buffer.colours = std::move(colours);
    cloud_buffer.classifications = std::move(classifications);
    cloud_buffer.branch_ids = std::move(branch_ids);

    for (int i = 0; i < (int)cloud_buffer.ends.size(); i++)
    {
      Cloud &chunk = is_outside(cloud_buffer, i) ? out_chunk : in_chunk;
      chunk.addRay(cloud_buffer, i);
    }
    in_writer.writeChunk(in_chunk);
    out_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
  };
  if (!Cloud::read(file_name, per_chunk)) return false;
  
  in_writer.end();
  out_writer.end();
  return true;
}

/// Special case for splitting around a plane.
bool splitPlane(const std::string &file_name, const std::string &in_name, const std::string &out_name,
                const Eigen::Vector3d &plane, const std::string& output_ext)
{
  CloudWriter inside_writer, outside_writer;
  if (!inside_writer.begin(in_name, output_ext)) return false;
  if (!outside_writer.begin(out_name, output_ext)) return false;
  Cloud in_chunk, out_chunk;

  // the split operation
  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                       std::vector<double> &times, std::vector<RGBA> &colours,
                       std::vector<uint8_t> &classifications, std::vector<uint16_t> &branch_ids) {
    const Eigen::Vector3d plane_vec = plane / plane.dot(plane);
    for (size_t i = 0; i < ends.size(); i++)
    {
      const double d1 = starts[i].dot(plane_vec) - 1.0;
      const double d2 = ends[i].dot(plane_vec) - 1.0;
      if (d1 * d2 > 0.0)
      {
        Cloud &chunk = d1 > 0.0 ? out_chunk : in_chunk;
        chunk.addRay(starts[i], ends[i], times[i], colours[i], classifications[i], branch_ids[i]);
      }
      else
      {
        RGBA col = colours[i];
        col.alpha = 0; // Mark as unbounded
        const Eigen::Vector3d mid = starts[i] + (ends[i] - starts[i]) * d1 / (d1 - d2);
        if (d1 > 0.0)
        {
          out_chunk.addRay(starts[i], mid, times[i], col, classifications[i], branch_ids[i]);
          in_chunk.addRay(mid, ends[i], times[i], colours[i], classifications[i], branch_ids[i]);
        }
        else
        {
          in_chunk.addRay(starts[i], mid, times[i], col, classifications[i], branch_ids[i]);
          out_chunk.addRay(mid, ends[i], times[i], colours[i], classifications[i], branch_ids[i]);
        }
      }
    }
    inside_writer.writeChunk(in_chunk);
    outside_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
  };
  if (!readPly(file_name, true, per_chunk, 0)) return false;

  inside_writer.end();
  outside_writer.end();
  return true;
}

/// Special case for splitting a capsule.
bool splitCapsule(const std::string &file_name, const std::string &in_name, const std::string &out_name,
                  const Eigen::Vector3d &end1, const Eigen::Vector3d &end2, double radius, const std::string& output_ext)
{
    CloudWriter inside_writer, outside_writer;
    if (!inside_writer.begin(in_name, output_ext)) return false;
    if (!outside_writer.begin(out_name, output_ext)) return false;
    Cloud in_chunk, out_chunk;

    Eigen::Vector3d dir = end2 - end1;
    double length = dir.norm();
    if (length > 0.0) dir /= length;

    auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                         std::vector<double> &times, std::vector<RGBA> &colours,
                         std::vector<uint8_t> &classifications, std::vector<uint16_t> &branch_ids) {
    for (size_t i = 0; i < ends.size(); i++)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      Eigen::Vector3d ray = end - start;
      
      double cylinder_intersection1 = 1e10;
      double cylinder_intersection2 = -1e10;
      Eigen::Vector3d up = dir.cross(ray);
      double mag = up.norm();
      if (mag > 0.0)
      {
        up /= mag;
        double gap = std::abs((start - end1).dot(up));
        if (gap >= radius)
        {
          out_chunk.addRay(start, end, times[i], colours[i], classifications[i], branch_ids[i]);
          continue;
        }
        Eigen::Vector3d lateral_dir = ray - dir * ray.dot(dir);
        double lateral_length = lateral_dir.norm();
        if (lateral_length == 0.0) continue;
        double d_mid = (end1 - start).dot(lateral_dir) / ray.dot(lateral_dir);
        
        double shift = std::sqrt(radius*radius - gap*gap) / lateral_length;
        double d_min = d_mid - shift;
        double d_max = d_mid + shift;
        double d1 = (start + ray*d_min - end1).dot(dir) / length;
        double d2 = (start + ray*d_max - end1).dot(dir) / length;
        if (d1 > 0.0 && d1 < 1.0) cylinder_intersection1 = d_min;
        if (d2 > 0.0 && d2 < 1.0) cylinder_intersection2 = d_max;
      }

      double ray_length = ray.norm();
      if(ray_length == 0.0) continue;
      double sphere_intersection1[2] = {1e10, 1e10};
      double sphere_intersection2[2] = {-1e10, -1e10};
      Eigen::Vector3d sphere_ends[2] = {end1, end2};
      for (int e = 0; e<2; e++)
      {
        double mid_d = (sphere_ends[e] - start).dot(ray) / (ray_length * ray_length);
        Eigen::Vector3d shortest_dir = (sphere_ends[e] - start) - ray * mid_d;
        double shortest_sqr = shortest_dir.squaredNorm();
        if (shortest_sqr < radius*radius)
        {
          double shift = std::sqrt(radius * radius - shortest_sqr) / ray_length;
          sphere_intersection1[e] = mid_d - shift;
          sphere_intersection2[e] = mid_d + shift;
        }
      }

      double closest_d = std::min({cylinder_intersection1, sphere_intersection1[0], sphere_intersection1[1]});
      double farthest_d = std::max({cylinder_intersection2, sphere_intersection2[0], sphere_intersection2[1]});

      RGBA black = {0,0,0,0};
      
      if (closest_d >= 1.0 || farthest_d <= 0.0)
      {
        out_chunk.addRay(start, end, times[i], colours[i], classifications[i], branch_ids[i]);
        continue;
      }
      
      if (closest_d > 0.0)
      {
        out_chunk.addRay(start, start + ray*closest_d, times[i], black, classifications[i], branch_ids[i]);
      }
      
      if (farthest_d < 1.0)
      {
        out_chunk.addRay(start + ray * farthest_d, end, times[i], colours[i], classifications[i], branch_ids[i]);
      }
      
      if (farthest_d < 1.0)
      {
        in_chunk.addRay(start + ray * std::max(0.0, closest_d), start + ray*std::min(farthest_d, 1.0), times[i], black, classifications[i], branch_ids[i]);
      }
      else
      {
        in_chunk.addRay(start + ray * std::max(0.0, closest_d), end, times[i], colours[i], classifications[i], branch_ids[i]);
      }
    }
    inside_writer.writeChunk(in_chunk);
    outside_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
    };
    if (!readPly(file_name, true, per_chunk, 0)) return false;

    inside_writer.end();
    outside_writer.end();
    return true;
}


/// Special case for splitting a box.
bool splitBox(const std::string &file_name, const std::string &in_name, const std::string &out_name,
              const Eigen::Vector3d &centre, const Eigen::Vector3d &extents, const std::string& output_ext)
{
  CloudWriter inside_writer, outside_writer;
  if (!inside_writer.begin(in_name, output_ext)) return false;
  if (!outside_writer.begin(out_name, output_ext)) return false;
  Cloud in_chunk, out_chunk;

  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                       std::vector<double> &times, std::vector<RGBA> &colours,
                       std::vector<uint8_t> &classifications, std::vector<uint16_t> &branch_ids) {
    const Cuboid cuboid(centre - extents, centre + extents);
    for (size_t i = 0; i < ends.size(); i++)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      if (cuboid.clipRay(start, end))
      {
        RGBA col = colours[i];
        if (!cuboid.intersects(ends[i])) col.alpha = 0;
        in_chunk.addRay(start, end, times[i], col, classifications[i], branch_ids[i]);
        if (start != starts[i])
        {
          RGBA out_col = {0,0,0,0};
          out_chunk.addRay(starts[i], start, times[i], out_col, classifications[i], branch_ids[i]);
        }
        if (ends[i] != end)
        {
          out_chunk.addRay(end, ends[i], times[i], colours[i], classifications[i], branch_ids[i]);
        }
      }
      else
      {
        out_chunk.addRay(starts[i], ends[i], times[i], colours[i], classifications[i], branch_ids[i]);
      }
    }
    inside_writer.writeChunk(in_chunk);
    outside_writer.writeChunk(out_chunk);
    in_chunk.clear();
    out_chunk.clear();
  };
  if (!readPly(file_name, true, per_chunk, 0)) return false;

  inside_writer.end();
  outside_writer.end();
  return true;
}

bool splitGrid(const std::string &file_name, const std::string &cloud_name_stub, const Eigen::Vector3d &cell_width,
               double overlap, const std::string& output_ext)
{
  return splitGrid(file_name, cloud_name_stub, Eigen::Vector4d(cell_width[0], cell_width[1], cell_width[2], 0), overlap, output_ext);
}

bool splitGrid(const std::string &file_name, const std::string &cloud_name_stub, const Eigen::Vector4d &cell_width,
               double overlap, const std::string& output_ext)
{
  overlap /= 2.0;
  Cloud::Info info;
  Cloud::getInfo(file_name, info);
  const Eigen::Vector3d &min_bound = info.rays_bound.min_bound_;
  const Eigen::Vector3d &max_bound = info.rays_bound.max_bound_;

  Eigen::Vector4d width = cell_width;
  for (int i = 0; i < 4; i++) {
    if (width[i] == 0.0) width[i] = std::numeric_limits<double>::max();
  }

  const Eigen::Vector3i min_index = Eigen::Vector3d(std::floor(0.5 + min_bound[0] / width[0]), std::floor(0.5 + min_bound[1] / width[1]), std::floor(0.5 + min_bound[2] / width[2])).cast<int>();
  const Eigen::Vector3i max_index = Eigen::Vector3d(std::ceil(0.5 + max_bound[0] / width[0]), std::ceil(0.5 + max_bound[1] / width[1]), std::ceil(0.5 + max_bound[2] / width[2])).cast<int>();
  const Eigen::Vector3i dimensions = max_index - min_index;
  
  const long int min_time = static_cast<long int>(std::floor(0.5 + info.min_time / width[3]));
  const long int max_time = static_cast<long int>(std::ceil(0.5 + info.max_time / width[3]));
  const int time_dimension = static_cast<int>(max_time - min_time);

  const int length = dimensions[0] * dimensions[1] * dimensions[2] * time_dimension;
  std::cout << "splitting into maximum of: " << length << " files" << std::endl;
  if (length > 50000) { std::cerr << "error: output of over 50,000 files is probably a mistake, exiting" << std::endl; return false; }
  
  const int max_open_files = 256;
  for (int pass = 0; pass<length; pass+=max_open_files)
  {
    std::vector<CloudWriter> cells(max_open_files);
    std::vector<Cloud> chunks(max_open_files);

    auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                         std::vector<double> &times, std::vector<RGBA> &colours,
                         std::vector<uint8_t> &classifications, std::vector<uint16_t> &branch_ids) {
      for (size_t i = 0; i < ends.size(); i++)
      {
        // --- START OF CORRECTION ---
        Eigen::Vector4d start_coords(starts[i][0], starts[i][1], starts[i][2], times[i]);
        Eigen::Vector4d end_coords(ends[i][0], ends[i][1], ends[i][2], times[i]);
        Eigen::Vector4d from = start_coords.array() / width.array() + 0.5;
        Eigen::Vector4d to = end_coords.array() / width.array() + 0.5;

        Eigen::Vector4d pos0 = from.cwiseMin(to);
        Eigen::Vector4d pos1 = from.cwiseMax(to);

        // Apply overlap to X and Y spatial dimensions only, not Z or Time, as per original intent
        pos0[0] -= overlap;
        pos0[1] -= overlap;
        pos1[0] += overlap;
        pos1[1] += overlap;

        Eigen::Vector4i minI = pos0.array().floor().cast<int>();
        Eigen::Vector4i maxI = pos1.array().ceil().cast<int>();
        // --- END OF CORRECTION ---
        
        for (int x = minI[0]; x < maxI[0]; x++)
        for (int y = minI[1]; y < maxI[1]; y++)
        for (int z = minI[2]; z < maxI[2]; z++)
        for (int t_idx = minI[3]; t_idx < maxI[3]; t_idx++)
        {
            const int time_dif = t_idx - min_time;
            int index = (x - min_index[0]) + dimensions[0] * (y - min_index[1]) + dimensions[0] * dimensions[1] * (z - min_index[2]) + dimensions[0] * dimensions[1] * dimensions[2] * time_dif;
            if (index < pass || index >= pass + max_open_files) continue;
            index -= pass;

            const Eigen::Vector3d box_min = (Eigen::Vector3d(x - 0.5, y - 0.5, z - 0.5).array() * width.head<3>().array()) - overlap;
            const Eigen::Vector3d box_max = (Eigen::Vector3d(x + 0.5, y + 0.5, z + 0.5).array() * width.head<3>().array()) + overlap;
            const Cuboid cuboid(box_min, box_max);
            Eigen::Vector3d start = starts[i];
            Eigen::Vector3d end = ends[i];

            if (cuboid.clipRay(start, end))
            {
                if (cells[index].fileName().empty())
                {
                    std::stringstream name;
                    name << cloud_name_stub;
                    if (cell_width[0] > 0.0) name << "_" << x;
                    if (cell_width[1] > 0.0) name << "_" << y;
                    if (cell_width[2] > 0.0) name << "_" << z;
                    if (cell_width[3] > 0.0) name << "_" << (long int)t_idx;
                    name << output_ext;
                    cells[index].begin(name.str(), output_ext);
                }
                RGBA col = colours[i];
                if (!cuboid.intersects(ends[i])) col.alpha = 0;
                chunks[index].addRay(start, end, times[i], col, classifications[i], branch_ids[i]);
            }
        }
      }
      for (int i = 0; i < max_open_files; i++)
      {
        if (chunks[i].rayCount() > 0)
        {
          cells[i].writeChunk(chunks[i]);
          chunks[i].clear();
        }
      }
    };
    if (!Cloud::read(file_name, per_chunk)) return false;
    for (int i = 0; i < max_open_files; i++) cells[i].end();
  }
  return true;
}

class RGBALess { public: bool operator()(const RGBA &a, const RGBA &b) const { if (a.red != b.red) return a.red < b.red; if (a.green != b.green) return a.green < b.green; return a.blue < b.blue; } };

bool splitColour(const std::string &file_name, const std::string &cloud_name_stub, bool seg_colour, const std::string& output_ext)
{
  std::map<RGBA, int, RGBALess> vox_map;
  int num_colours = 0;
  auto count_colours = [&vox_map, &num_colours](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &,
                                                std::vector<double> &, std::vector<ray::RGBA> &colours,
                                                std::vector<uint8_t>&, std::vector<uint16_t>&) {
    for (auto &colour : colours)
    {
      if (vox_map.find(colour) == vox_map.end())
      {
        vox_map.insert(std::pair<RGBA, int>(colour, num_colours++));
      }
    }
  };
  if (!ray::Cloud::read(file_name, count_colours)) return false;

  const int max_total_files = 50000;
  if (num_colours > max_total_files) { std::cerr << "Error: " << num_colours << " colours generates more than the maximum number of files: " << max_total_files << std::endl; return false; }
  
  const int max_files_at_once = 512;  // operating systems will fail with too many open file pointers.
  std::cout << "splitting into: " << num_colours << " files" << std::endl;
  if (num_colours > max_files_at_once) { std::cout << "Warning: cloud has more unique colours than allowed for simultaneous files " << max_files_at_once << " so using multiple passes." << std::endl; }
  
  for (int batch = 0; batch < num_colours; batch+=max_files_at_once)
  {
    std::vector<CloudWriter> cells(max_files_at_once);
    std::vector<Cloud> chunks(max_files_at_once);
    
    auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                         std::vector<double> &times, std::vector<RGBA> &colours,
                         std::vector<uint8_t> &classifications, std::vector<uint16_t> &branch_ids) {
      for (size_t i = 0; i < ends.size(); i++)
      {
        RGBA colour = colours[i];
        const auto &vox = vox_map.find(colour);
        if (vox != vox_map.end())
        {
          int index = vox->second - batch;
          if (index < 0 || index >= max_files_at_once) continue;

          if (cells[index].fileName().empty())
          {
            std::stringstream name;
            if (seg_colour) name << cloud_name_stub << "_" << convertColourToInt(colour) << output_ext;
            else name << cloud_name_stub << "_" << (int)colour.red << "_" << (int)colour.green << "_" << (int)colour.blue << output_ext;
            cells[index].begin(name.str(), output_ext);
          }
          chunks[index].addRay(starts[i], ends[i], times[i], colours[i], classifications[i], branch_ids[i]);
        }
      }
      for (size_t i = 0; i < chunks.size(); i++)
      {
        if (chunks[i].rayCount() > 0)
        {
          cells[i].writeChunk(chunks[i]);
          chunks[i].clear();
        }
      }
    };
    if (!Cloud::read(file_name, per_chunk)) return false;
    for (size_t i = 0; i < cells.size(); i++) cells[i].end();
  }
  return true;
}

}  // namespace ray