// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayply.h"
#include "raylib/rayprogress.h"
#include "raylib/rayprogressthread.h"
#include "raymesh.h"

#include <fstream>
#include <iostream>
// #define OUTPUT_MOMENTS // useful when setting up unit test expected ray clouds

namespace ray
{
namespace
{
// these are set once and are constant after that
unsigned long chunk_header_length = 0;
unsigned long point_cloud_chunk_header_length = 0;
unsigned long vertex_size_pos = 0;
unsigned long point_cloud_vertex_size_pos = 0;

enum DataType
{
  kDTfloat,
  kDTdouble,
  kDTushort,
  kDTuchar,
  kDTint,
  kDTnone
};
}  // namespace

bool writeRayCloudChunkStart(const std::string &file_name, std::ofstream &out)
{
  int num_zeros = std::numeric_limits<unsigned long>::digits10;
  out.open(file_name, std::ios::binary | std::ios::out);
  if (out.fail())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing." << std::endl;
    return false;
  }
  out << "ply" << std::endl;
  out << "format binary_little_endian 1.0" << std::endl;
  out << "comment generated by raycloudtools library" << std::endl;
  out << "element vertex ";
  for (int i = 0; i < num_zeros; i++)
    out << "0";  // fill in with zeros. I will replace rightmost characters later, to give actual number
  vertex_size_pos = out.tellp();
  out << std::endl;
#if RAYLIB_DOUBLE_RAYS
  out << "property double x" << std::endl;
  out << "property double y" << std::endl;
  out << "property double z" << std::endl;
#else
  out << "property float x" << std::endl;
  out << "property float y" << std::endl;
  out << "property float z" << std::endl;
#endif
  out << "property double time" << std::endl;
#if RAYLIB_WITH_NORMAL_FIELD
  out << "property float nx" << std::endl;
  out << "property float ny" << std::endl;
  out << "property float nz" << std::endl;
#else
  out << "property float rayx" << std::endl;
  out << "property float rayy" << std::endl;
  out << "property float rayz" << std::endl;
#endif
  out << "property uchar red" << std::endl;
  out << "property uchar green" << std::endl;
  out << "property uchar blue" << std::endl;
  out << "property uchar alpha" << std::endl;
  out << "end_header" << std::endl;
  chunk_header_length = out.tellp();
  return true;
}

bool writeRayCloudChunk(std::ofstream &out, RayPlyBuffer &vertices, const std::vector<Eigen::Vector3d> &starts,
                        const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                        const std::vector<RGBA> &colours, bool &has_warned)
{
  if (ends.size() == 0)
  {
    // this is not an error. Allowing empty chunks avoids wrapping every call to writeRayCloudChunk in a condition
    return true;
  }
  if (out.tellp() < (long)chunk_header_length)
  {
    std::cerr << "Error: file header has not been written, use writeRayCloudChunkStart" << std::endl;
    return false;
  }
  vertices.resize(ends.size());

  for (size_t i = 0; i < ends.size(); i++)
  {
    if (!has_warned)
    {
      if (!(ends[i] == ends[i]))
      {
        std::cout << "WARNING: nans in point: " << i << ": " << ends[i].transpose() << std::endl;
        has_warned = true;
      }
#if !RAYLIB_DOUBLE_RAYS
      if (std::abs(ends[i][0]) > 100000.0)
      {
        std::cout << "WARNING: very large point location at: " << i << ": " << ends[i].transpose() << ", suspicious"
                  << std::endl;
        has_warned = true;
      }
#endif
      bool b = starts[i] == starts[i];
      if (!b)
      {
        std::cout << "WARNING: nans in start: " << i << ": " << starts[i].transpose() << std::endl;
        has_warned = true;
      }
    }
    Eigen::Vector3d n = starts[i] - ends[i];
    union U  // TODO: this is nasty, better to just make vertices an unsigned char vector
    {
      float f[2];
      double d;
    };
    U u;
    u.d = times[i];

#if RAYLIB_DOUBLE_RAYS
    U end0, end1, end2;
    end0.d = ends[i][0];
    end1.d = ends[i][1];
    end2.d = ends[i][2];
    vertices[i] << end0.f[0], end0.f[1], end1.f[0], end1.f[1], end2.f[0], end2.f[1], u.f[0], u.f[1], (float)n[0],
      (float)n[1], (float)n[2], (float &)colours[i];
#else
    vertices[i] << (float)ends[i][0], (float)ends[i][1], (float)ends[i][2], u.f[0], u.f[1], (float)n[0], (float)n[1],
      (float)n[2], (float &)colours[i];
#endif
  }
  out.write((const char *)&vertices[0], sizeof(RayPlyEntry) * vertices.size());
  if (!out.good())
  {
    std::cerr << "error writing to file" << std::endl;
    return false;
  }
  return true;
}

unsigned long writeRayCloudChunkEnd(std::ofstream &out)
{
  const unsigned long size = static_cast<unsigned long>(out.tellp()) - chunk_header_length;
  const unsigned long number_of_rays = size / sizeof(RayPlyEntry);
  std::stringstream stream;
  stream << number_of_rays;
  std::string str = stream.str();
  out.seekp(vertex_size_pos - str.length());
  out << str;
  return number_of_rays;
}

// Save the polygon file to disk
bool writePlyRayCloud(const std::string &file_name, const std::vector<Eigen::Vector3d> &starts,
                      const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                      const std::vector<RGBA> &colours)
{
  std::vector<RGBA> rgb(times.size());
  if (colours.size() > 0)
    rgb = colours;
  else
    colourByTime(times, rgb);

  std::ofstream ofs;
  if (!writeRayCloudChunkStart(file_name, ofs))
    return false;
  RayPlyBuffer buffer;
  bool has_warned = false;
  // TODO: could split this into chunks aswell, it would allow saving out files roughly twice as large
  if (!writeRayCloudChunk(ofs, buffer, starts, ends, times, rgb, has_warned))
  {
    return false;
  }
  const unsigned long num_rays = ray::writeRayCloudChunkEnd(ofs);
  std::cout << num_rays << " rays saved to " << file_name << std::endl;
  return true;
}

// Point cloud chunked writing

bool writePointCloudChunkStart(const std::string &file_name, std::ofstream &out)
{
  int num_zeros = std::numeric_limits<unsigned long>::digits10;
  std::cout << "saving to " << file_name << " ..." << std::endl;
  out.open(file_name, std::ios::binary | std::ios::out);
  if (out.fail())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing." << std::endl;
    return false;
  }
  out << "ply" << std::endl;
  out << "format binary_little_endian 1.0" << std::endl;
  out << "comment generated by raycloudtools library" << std::endl;
  out << "element vertex ";
  for (int i = 0; i < num_zeros; i++)
    out << "0";  // fill in with zeros. I will replace rightmost characters later, to give actual number
  point_cloud_vertex_size_pos = out.tellp();  // same value as for ray cloud, so we can use the same varia
  out << std::endl;
#if RAYLIB_DOUBLE_RAYS
  out << "property double x" << std::endl;
  out << "property double y" << std::endl;
  out << "property double z" << std::endl;
#else
  out << "property float x" << std::endl;
  out << "property float y" << std::endl;
  out << "property float z" << std::endl;
#endif
  out << "property double time" << std::endl;
  out << "property uchar red" << std::endl;
  out << "property uchar green" << std::endl;
  out << "property uchar blue" << std::endl;
  out << "property uchar alpha" << std::endl;
  out << "end_header" << std::endl;
  point_cloud_chunk_header_length = out.tellp();
  return true;
}

bool writePointCloudChunk(std::ofstream &out, PointPlyBuffer &vertices, const std::vector<Eigen::Vector3d> &points,
                          const std::vector<double> &times, const std::vector<RGBA> &colours, bool &has_warned)
{
  if (points.size() == 0)
  {
    std::cerr << "Error: saving out ray file chunk with zero rays" << std::endl;
    return false;
  }
  if (out.tellp() < (long)point_cloud_chunk_header_length)
  {
    std::cerr << "Error: file header has not been written, use writeRayCloudChunkStart" << std::endl;
    return false;
  }
  vertices.resize(points.size());  // allocates the chunk size the first time, and nullop on subsequent chunks

  for (size_t i = 0; i < points.size(); i++)
  {
    if (!has_warned)
    {
      if (!(points[i] == points[i]))
      {
        std::cout << "WARNING: nans in point: " << i << ": " << points[i].transpose() << std::endl;
        has_warned = true;
      }
#if !RAYLIB_DOUBLE_RAYS
      if (std::abs(points[i][0]) > 100000.0)
      {
        std::cout << "WARNING: very large point location at: " << i << ": " << points[i].transpose() << ", suspicious"
                  << std::endl;
        has_warned = true;
      }
#endif
    }
    union U  // TODO: this is nasty, better to just make vertices an unsigned char vector
    {
      float f[2];
      double d;
    };
    U u;
    u.d = times[i];
#if RAYLIB_DOUBLE_RAYS
    U end0, end1, end2;
    end0.d = points[i][0];
    end1.d = points[i][1];
    end2.d = points[i][2];
    vertices[i] << end0.f[0], end0.f[1], end1.f[0], end1.f[1], end2.f[0], end2.f[1], u.f[0], u.f[1],
      (float &)colours[i];
#else
    vertices[i] << (float)points[i][0], (float)points[i][1], (float)points[i][2], (float)u.f[0], (float)u.f[1],
      (float &)colours[i];
#endif
  }
  out.write((const char *)&vertices[0], sizeof(PointPlyEntry) * vertices.size());
  if (!out.good())
  {
    std::cerr << "error writing to file" << std::endl;
    return false;
  }
  return true;
}

void writePointCloudChunkEnd(std::ofstream &out)
{
  const unsigned long size = static_cast<unsigned long>(out.tellp()) - point_cloud_chunk_header_length;
  const unsigned long number_of_points = size / sizeof(PointPlyEntry);
  std::stringstream stream;
  stream << number_of_points;
  std::string str = stream.str();
  out.seekp(point_cloud_vertex_size_pos - str.length());
  out << str;
  std::cout << "... saved out " << number_of_points << " points." << std::endl;
}

// Save the polygon point-cloud file to disk
bool writePlyPointCloud(const std::string &file_name, const std::vector<Eigen::Vector3d> &points,
                        const std::vector<double> &times, const std::vector<RGBA> &colours)
{
  std::vector<RGBA> rgb(times.size());
  if (colours.size() > 0)
  {
    rgb = colours;
  }
  else
  {
    colourByTime(times, rgb);
  }

  std::ofstream ofs;
  if (!writePointCloudChunkStart(file_name, ofs))
    return false;
  PointPlyBuffer buffer;
  bool has_warned = false;
  // TODO: could split this into chunks aswell, it would allow saving out files roughly twice as large
  if (!writePointCloudChunk(ofs, buffer, points, times, rgb, has_warned))
  {
    return false;
  }
  writePointCloudChunkEnd(ofs);
  return true;
}

bool readPly(const std::string &file_name, bool is_ray_cloud,
             std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                std::vector<double> &times, std::vector<RGBA> &colours)>
               apply, 
             double max_intensity, bool times_optional, size_t chunk_size)
{
  std::cout << "reading: " << file_name << std::endl;
  std::ifstream input(file_name.c_str(), std::ios::in | std::ios::binary);
  if (input.fail())
  {
    std::cerr << "Couldn't open file: " << file_name << std::endl;
    return false;
  }
  std::string line;
  int row_size = 0;
  int offset = -1, normal_offset = -1, time_offset = -1, colour_offset = -1;
  int intensity_offset = -1;
  bool time_is_float = false;
  bool pos_is_float = false;
  bool normal_is_float = false;
  DataType intensity_type = kDTnone;
  int rowsteps[] = { int(sizeof(float)), int(sizeof(double)), int(sizeof(unsigned short)), int(sizeof(unsigned char)), int(sizeof(int)),
                     0 };  // to match each DataType enum

  while (line != "end_header\r" && line != "end_header")
  {
    if (!getline(input, line))
    {
      break;
    }

    if (line.find("format ascii 1.0") != std::string::npos)
    {
      std::cerr << "ASCII PLY not supported " << file_name << std::endl;
      return false;
    }

    // support multiple data types
    DataType data_type = kDTnone;
    if (line.find("property float") != std::string::npos)
      data_type = kDTfloat;
    else if (line.find("property double") != std::string::npos)
      data_type = kDTdouble;
    else if (line.find("property uchar") != std::string::npos || line.find("property uint8") != std::string::npos)
      data_type = kDTuchar;
    else if (line.find("property ushort") != std::string::npos)
      data_type = kDTushort;    
    else if (line.find("property int") != std::string::npos)
      data_type = kDTint;

    if (line == "property float x" || line == "property double x")
    {
      offset = row_size;
      if (line.find("float") != std::string::npos)
        pos_is_float = true;
    }
    if (line == "property float rayx" || line == "property double rayx")
    {
#if RAYLIB_WITH_NORMAL_FIELD
      if (normal_offset == -1)
#endif
      {
        normal_offset = row_size;
        normal_is_float = line.find("float") != std::string::npos;
      }
    }
    // Support both standard PLY format and CloudCompare format with scalar_ prefixes
    if (line == "property float nx" || line == "property double nx" || 
        line == "property float scalar_nx" || line == "property double scalar_nx")
    {
#if !RAYLIB_WITH_NORMAL_FIELD
      if (normal_offset == -1)
#endif
      {
        normal_offset = row_size;
        normal_is_float = line.find("float") != std::string::npos;
      }
    }
    if (line.find("time") != std::string::npos || line.find("scalar_time") != std::string::npos)
    {
      time_offset = row_size;
      if (line.find("float") != std::string::npos)
        time_is_float = true;
    }
    if (line.find("intensity") != std::string::npos || line.find("scalar_alpha") != std::string::npos)
    {
      intensity_offset = row_size;
      intensity_type = data_type;
    }
    if (line == "property uchar red" || line == "property uint8 red")
      colour_offset = row_size;

    row_size += rowsteps[data_type];
  }
  if (offset == -1)
  {
    std::cerr << "could not find position properties of file: " << file_name << std::endl;
    return false;
  }
  if (is_ray_cloud && normal_offset == -1)
  {
    std::cerr << "could not find normal properties of file: " << file_name << std::endl;
    std::cerr << "ray clouds store the ray starts using the normal field" << std::endl;
    return false;
  }

  std::streampos start = input.tellg();
  input.seekg(0, input.end);
  size_t length = input.tellg() - start;
  input.seekg(start);
  size_t size = length / row_size;

  ray::Progress progress;
  ray::ProgressThread progress_thread(progress);
  size_t num_chunks = (size + (chunk_size - 1)) / chunk_size;
  progress.begin("read and process", num_chunks);

  std::vector<unsigned char> vertices(row_size);
  bool warning_set = false;
  if (size == 0)
  {
    std::cerr << "no entries found in ply file" << std::endl;
    return false;
  }
  if (time_offset == -1)
  {
    if (times_optional)
    {
      std::cout << "Warning: no times provided in file, applying 1 second difference per ray consecutively, starting at 0 seconds" << std::endl;
    }
    else
    {
      std::cerr << "error: no time information found in " << file_name << std::endl;
      return false;
    }
  }
  if (colour_offset == -1)
  {
    std::cout << "warning: no colour information found in " << file_name
              << ", setting colours red->green->blue based on time" << std::endl;
  }
  if (!is_ray_cloud && intensity_offset != -1)
  {
    if (colour_offset != -1)
    {
      std::cout << "warning: intensity and colour information both found in file. Replacing alpha with intensity value."
                << std::endl;
    }
    else
    {
      std::cout << "intensity information found in file, storing this in the ray cloud 8-bit alpha channel."
                << std::endl;
    }
  }

  // pre-reserving avoids memory fragmentation
  std::vector<Eigen::Vector3d> ends;
  std::vector<Eigen::Vector3d> starts;
  std::vector<double> times;
  std::vector<ray::RGBA> colours;
  std::vector<uint8_t> intensities;
  size_t reserve_size = std::min(chunk_size, size);
  ends.reserve(reserve_size);
  starts.reserve(reserve_size);
  if (time_offset != -1)
    times.reserve(reserve_size);
  if (colour_offset != -1)
    colours.reserve(reserve_size);
  if (intensity_offset != -1)
    intensities.reserve(reserve_size);
  bool any_returns = false;
  int identical_times = 0;
  double last_time = std::numeric_limits<double>::lowest();
  double last_unique_time = std::numeric_limits<double>::lowest();
  
  for (size_t i = 0; i < size; i++)
  {
    input.read((char *)&vertices[0], row_size);
    Eigen::Vector3d end;
    if (pos_is_float)
    {
      Eigen::Vector3f e = (Eigen::Vector3f &)vertices[offset];
      end = Eigen::Vector3d(e[0], e[1], e[2]);
    }
    else
    {
      end = (Eigen::Vector3d &)vertices[offset];
    }
    bool end_valid = end == end;
    if (!warning_set)
    {
      if (!end_valid)
      {
        std::cout << "warning, NANs in point " << i << ", removing all NANs." << std::endl;
        warning_set = true;
      }
      if (std::abs(end[0]) > 100000.0)
      {
        std::cout << "warning: very large data in point " << i << ", suspicious: " << end.transpose() << std::endl;
        warning_set = true;
      }
    }
    if (!end_valid)
      continue;

    Eigen::Vector3d normal(0, 0, 0);
    if (is_ray_cloud)
    {
      if (normal_is_float)
      {
        Eigen::Vector3f n = (Eigen::Vector3f &)vertices[normal_offset];
        normal = Eigen::Vector3d(n[0], n[1], n[2]);
      }
      else
      {
        normal = (Eigen::Vector3d &)vertices[normal_offset];
      }
      bool norm_valid = normal == normal;
      if (!warning_set)
      {
        if (!norm_valid)
        {
          std::cout << "warning, NANs in raystart stored in normal " << i << ", removing all such rays." << std::endl;
          warning_set = true;
        }
      }
      if (!norm_valid)
        continue;
      if (std::abs(normal[0]) > 100000.0 && !warning_set)
      {
        std::cerr << "Error: very large ray length in ray index " << i << " " << normal.transpose() << ", bad input." << std::endl;
        std::cerr << "Use rayexport then rayimport the exported point cloud with a fixed trajectory file" << std::endl;
        warning_set = true;
      }        
    }

    starts.push_back(end + normal);
    ends.push_back(end);
    if (time_offset != -1)
    {
      double time;
      if (time_is_float)
      {
        time = (double)((float &)vertices[time_offset]);
      }
      else
      {  
        time = (double &)vertices[time_offset];
      }
      if (!is_ray_cloud)
      {
        if (time==last_unique_time)
        {
          const double time_delta = 1e-6; // this is a sufficient difference for rayrestore (see time_eps in rayrestore.cpp)
          time = last_time + time_delta;
          identical_times++;
        }
        else
        {
          last_unique_time = time;
        }
        last_time = time;
      }
      times.push_back(time);
    }

    if (colour_offset != -1)
    {
      RGBA colour = (RGBA &)vertices[colour_offset];
      colours.push_back(colour);
    }
    if (!is_ray_cloud)
    {
      if (intensity_offset != -1)
      {
        double intensity;
        if (intensity_type == kDTfloat)
          intensity = (double)((float &)vertices[intensity_offset]);
        else if (intensity_type == kDTdouble)
          intensity = (double &)vertices[intensity_offset];
        else  // (intensity_type == kDTushort)
          intensity = (double)((unsigned short &)vertices[intensity_offset]);
        if (intensity >= 0.0)
        {
          // only intensity exactly 0 will be used for alpha=0 in uint_8 format.
          intensity = std::ceil(255.0 * clamped(intensity / max_intensity, 0.0, 1.0));  
        }
        // support for special codes for out of range cases, defined by intensity:
        // -1 non-return of unknown length
        // -2 the object is within minimum range, so range is not certain but small
        // -3 outside maximum range, so range is uncertain but large
        else if (intensity == -1.0) 
        {
          intensity = 0.0;
        }
        else // here a range is specified, just low certainty. We choose to this range.
        {
          intensity = 1.0;
        }
        intensities.push_back(static_cast<uint8_t>(intensity));
      }
    }
    if (ends.size() == chunk_size || i == size - 1)
    {
      if (time_offset == -1)
      {
        times.resize(ends.size());
        for (size_t j = 0; j < times.size(); j++) 
        {
          times[j] = (double)(i + j);
        }
      }
      if (colour_offset == -1)
      {
        colourByTime(times, colours);
      }
      if (!is_ray_cloud)
      {
        if (intensity_offset != -1)
        {
          for (size_t j = 0; j < intensities.size(); j++)
          {
            colours[j].alpha = intensities[j];
            // colour zero-intensity rays black. This is a helpful debug tool.
            if (intensities[j] == 0)
            {
              colours[j].red = colours[j].green = colours[j].blue = 0;
            }
            else
            {
              any_returns = true;
            }
          }
        }
        else
        {
          for (size_t j = 0; j < colours.size(); j++)
          {
            if (colours[j].alpha == 0)
            {
              // colour zero-intensity rays black. This is a helpful debug tool.
              colours[j].red = colours[j].green = colours[j].blue = 0;
            }
            else
            {
              any_returns = true;
            }
          }
        }
      }
      apply(starts, ends, times, colours);
      starts.clear();
      ends.clear();
      times.clear();
      colours.clear();
      intensities.clear();
      progress.increment();

      if (!is_ray_cloud && i==size-1 && identical_times > 0)
      {
        std::cout << std::endl;
        std::cout << "warning: " << identical_times << "/" << size << " rays have identical times," << std::endl;
        std::cout << "since rayrestore relies on unique time stamps, a 1 microsecond increment has been applied for these times." << std::endl;
      }
    }
  }
  progress.end();
  progress_thread.requestQuit();
  progress_thread.join();

  if (!is_ray_cloud && any_returns == false) // no return rays
  {
    std::cerr << "Error: ray cloud has no identified points; all rays are zero-intensity non-returns," << std::endl;
    std::cerr << "many functions will not operate on this degerenate case." << std::endl;
    std::cerr << "Use raycolour cloud.ply alpha 1 to force all rays to have end points and full intensity." << std::endl;
    std::cerr << "Or re-import using rayimport points.ply --max_intensity 0." << std::endl;
  }

  return true;
}

bool readPly(const std::string &file_name, std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
             std::vector<double> &times, std::vector<RGBA> &colours, bool is_ray_cloud, double max_intensity)
{
  // Note: this lambda function assumes that the passed in vectors are end-of-life, and can be moved
  // this is true for the readPly function, with maximum chunk size.
  auto apply = [&](std::vector<Eigen::Vector3d> &start_points, std::vector<Eigen::Vector3d> &end_points,
                   std::vector<double> &time_points, std::vector<RGBA> &colour_values) 
  {
    starts.insert(starts.end(), start_points.begin(), start_points.end());
    ends.insert(ends.end(), end_points.begin(), end_points.end());
    times.insert(times.end(), time_points.begin(), time_points.end());
    colours.insert(colours.end(), colour_values.begin(), colour_values.end());
  };
  return readPly(file_name, is_ray_cloud, apply, max_intensity, std::numeric_limits<size_t>::max());
}

bool writePlyMesh(const std::string &file_name, const Mesh &mesh, bool flip_normals)
{
  std::cout << "saving to " << file_name << ", " << mesh.vertices().size() << " vertices." << std::endl;
  if (mesh.indexList().size() == 0 || mesh.vertices().size() < 3)
  {
    std::cout << "Warning: mesh is empty or too small to save. Num vertices: " << mesh.vertices().size() << std::endl;
    return false;
  }

#if RAYLIB_DOUBLE_RAYS
  int row_size = 28;
  struct Vert
  {
    Eigen::Vector3d pos;
    ray::RGBA colour;
  };
#else 
  int row_size = 16;
  struct Vert
  {
    Eigen::Vector3f pos;
    ray::RGBA colour;
  };
#endif 
  std::vector<unsigned char> vertices(row_size * mesh.vertices().size());  // 4d to give space for colour
  for (size_t i = 0; i < mesh.vertices().size(); i++)
  {
    Vert vert;
#if RAYLIB_DOUBLE_RAYS
    vert.pos = mesh.vertices()[i];
#else
    vert.pos = mesh.vertices()[i].cast<float>();
#endif
    vert.colour = mesh.colours().empty() ? ray::RGBA(127,127,127,255) : mesh.colours()[i];
    memcpy(&vertices[row_size * i], &vert, row_size);
  }

  FILE *fid = fopen(file_name.c_str(), "w+");
  if (!fid)
  {
    std::cerr << "error opening file " << file_name << " for writing." << std::endl;
    return false;
  }
  fprintf(fid, "ply\n");
  fprintf(fid, "format binary_little_endian 1.0\n");
  fprintf(fid, "comment SDK generated\n");  // TODO: add version here
  if (!mesh.uvList().empty())
  {
    std::string tex_name = mesh.textureName() == "" ? "wood_texture.png" : mesh.textureName();
    fprintf(fid, "comment TextureFile %s\n", tex_name.c_str()); 
  }
  fprintf(fid, "element vertex %u\n", unsigned(mesh.vertices().size()));
#if RAYLIB_DOUBLE_RAYS
  fprintf(fid, "property double x\n");
  fprintf(fid, "property double y\n");
  fprintf(fid, "property double z\n");
#else
  fprintf(fid, "property float x\n");
  fprintf(fid, "property float y\n");
  fprintf(fid, "property float z\n");
#endif
  fprintf(fid, "property uchar red\n");
  fprintf(fid, "property uchar green\n");
  fprintf(fid, "property uchar blue\n");
  fprintf(fid, "property uchar alpha\n");
  fprintf(fid, "element face %u\n", (unsigned)mesh.indexList().size());
  fprintf(fid, "property list int int vertex_indices\n");
  if (!mesh.uvList().empty())
  {
    fprintf(fid, "property list int float texcoord\n");
    fprintf(fid, "property int texnumber\n");
  }
  fprintf(fid, "end_header\n");

  fwrite(&vertices[0], sizeof(unsigned char), vertices.size(), fid);

  auto &list = mesh.indexList();
  size_t written;
  if (mesh.uvList().empty())
  {
    std::vector<Eigen::Vector4i> triangles(mesh.indexList().size());
    if (flip_normals)
      for (size_t i = 0; i < list.size(); i++) triangles[i] = Eigen::Vector4i(3, list[i][2], list[i][1], list[i][0]);
    else
      for (size_t i = 0; i < list.size(); i++) triangles[i] = Eigen::Vector4i(3, list[i][0], list[i][1], list[i][2]);
    written = fwrite(&triangles[0], sizeof(Eigen::Vector4i), triangles.size(), fid);
  }
  else
  {
    struct Face
    {
      int num_corners;
      Eigen::Vector3i ids;
      int num_coords;
      float uvs[6];
      int texnumber;
    };
    auto &uvs = mesh.uvList();
    std::vector<Face> faces(uvs.size());
    for (size_t i = 0; i < list.size(); i++) 
    {
      faces[i].num_corners = 3;
      if (flip_normals)
        faces[i].ids = Eigen::Vector3i(list[i][2], list[i][1], list[i][0]);
      else
        faces[i].ids = Eigen::Vector3i(list[i][0], list[i][1], list[i][2]);
      
      faces[i].num_coords = 6;
      faces[i].uvs[0] = uvs[i][0].real();
      faces[i].uvs[1] = uvs[i][0].imag();
      faces[i].uvs[2] = uvs[i][1].real();
      faces[i].uvs[3] = uvs[i][1].imag();
      faces[i].uvs[4] = uvs[i][2].real();
      faces[i].uvs[5] = uvs[i][2].imag();
      faces[i].texnumber = 0;
    }
    written = fwrite(&faces[0], sizeof(Face), faces.size(), fid);
  }

  fclose(fid);
  if (written == 0)
  {
    std::cerr << "Error writing to file " << file_name << std::endl;
    return false;
  }

#if defined OUTPUT_MOMENTS // Only used to supply data to unit tests
  Eigen::Array<double, 6, 1> mom = mesh.getMoments();
  std::cout << "stats: " << std::endl;
  for (int i = 0; i < mom.rows(); i++) 
  { 
    std::cout << ", " << mom[i];
  }
  std::cout << std::endl;
#endif  // defined OUTPUT_MOMENTS
  return true;
}

bool readPlyMesh(const std::string &file, Mesh &mesh)
{
  std::ifstream input(file.c_str(), std::ios::in | std::ios::binary);
  if (input.fail())
  {
    std::cerr << "Couldn't open file: " << file << std::endl;
    return false;
  }
  std::string line;
  unsigned row_size = 0;
  unsigned number_of_faces = 0;
  unsigned number_of_vertices = 0;
  char char1[100], char2[100], char3[100], char4[100], char5[100];
  int pos_offset = -1;
  bool pos_is_float = false;
  int order_size = 0;
  int vertex_index_size = 0;
  int uv_order_size = 0;
  int uv_size = 0;
  int texnumber_size = 0;
  while (line != "end_header\r" && line != "end_header")
  {
    if (!getline(input, line))
    {
      break;
    }
    if (line.find("format ascii") != std::string::npos)
    {
      std::cerr << "error: can only read in binary ply mesh files" << std::endl;
      return false;
    }
    if (line.find("property float x") != std::string::npos)
    {
      pos_offset = row_size;
      pos_is_float = true;
    }
    if (line.find("property double x") != std::string::npos)
    {
      pos_offset = row_size;
    }
    if (line.find("property float") != std::string::npos)
      row_size += int(sizeof(float));
    else if (line.find("property double") != std::string::npos)
      row_size += int(sizeof(double));
    else if (line.find("property uchar") != std::string::npos || line.find("property uint8") != std::string::npos)
      row_size += int(sizeof(unsigned char));
    else if (line.find("property ushort") != std::string::npos)
      row_size += int(sizeof(unsigned short));

    if (line.find("element vertex") != std::string::npos)
    {
      sscanf(line.c_str(), "%s %s %u", char1, char2, &number_of_vertices);
    }
    if (line.find("element face") != std::string::npos)
    {
      sscanf(line.c_str(), "%s %s %u", char1, char2, &number_of_faces);
    }
    if (line.find("property list") != std::string::npos && line.find("vertex_indices") != std::string::npos)
    {
      sscanf(line.c_str(), "%s %s %s %s %s", char1, char2, char3, char4, char5);
      std::string type1(char3);
      std::string type2(char4);
      order_size        = type1 == "uchar" ? int(sizeof(unsigned char)) : (type1 == "ushort" ? int(sizeof(unsigned short)) : int(sizeof(int)));
      vertex_index_size = type2 == "uchar" ? int(sizeof(unsigned char)) : (type2 == "ushort" ? int(sizeof(unsigned short)) : int(sizeof(int)));
    }
    if (line.find("property list") != std::string::npos && line.find("texcoord") != std::string::npos)
    {
      sscanf(line.c_str(), "%s %s %s %s %s", char1, char2, char3, char4, char5);
      std::string type1(char3);
      std::string type2(char4);
      uv_order_size     = type1 == "uchar" ? int(sizeof(unsigned char)) : (type1 == "ushort" ? int(sizeof(unsigned short)) : int(sizeof(int)));
      uv_size = type2 == "float" ? int(sizeof(float)) : int(sizeof(double));
    }
    if (line.find("property int texnumber") != std::string::npos)
    {
      texnumber_size = 4;
    }
  }
  if (pos_offset == -1)
  {
    std::cerr << "error, mesh file does not contain x, y, z data" << std::endl;
    return false;
  }

  mesh.vertices().resize(number_of_vertices);
  std::vector<unsigned char> vertices(row_size);
  std::cout << "row size: " << row_size << std::endl;
  for (unsigned int i = 0; i<number_of_vertices; i++)
  {
    input.read((char *)&vertices[0], row_size);
    if (pos_is_float)
    {
      Eigen::Vector3f e = (Eigen::Vector3f &)vertices[pos_offset];
      mesh.vertices()[i] = Eigen::Vector3d(e[0], e[1], e[2]);
    }
    else
    {
      mesh.vertices()[i] = (Eigen::Vector3d &)vertices[pos_offset];
    }
  }
  row_size = order_size + 3*vertex_index_size + uv_order_size + 6*uv_size + texnumber_size;
  std::vector<unsigned char> triangles(number_of_faces * row_size);
  input.read((char *)&triangles[0], triangles.size());

  mesh.indexList().resize(number_of_faces);
  if (uv_size > 0)
  {
    mesh.uvList().resize(number_of_faces);
  }
  for (unsigned int i = 0; i < number_of_faces; i++)
  {
    for (int j = 0; j<3; j++)
    {
      if (vertex_index_size == int(sizeof(unsigned char))) 
      {
        unsigned char index = (unsigned char &)triangles[row_size*i + order_size + j*vertex_index_size];
        mesh.indexList()[i][j] = int(index);
      }
      else if (vertex_index_size == int(sizeof(unsigned short))) 
      {
        unsigned short index = (unsigned short &)triangles[row_size*i + order_size + j*vertex_index_size];
        mesh.indexList()[i][j] = int(index);
      }
      else if (vertex_index_size == int(sizeof(int))) 
      {
        mesh.indexList()[i][j] = (int &)triangles[row_size*i + order_size + j*vertex_index_size];
      }
    }
    if (uv_size > 0)
    {
      for (int j = 0; j<3; j++)
      {
        if (uv_size == 4)
        {
          float u = (float &)triangles[row_size*i + order_size + 3*vertex_index_size + uv_order_size + 2*j*uv_size];
          float v = (float &)triangles[row_size*i + order_size + 3*vertex_index_size + uv_order_size + (2*j + 1)*uv_size];
          mesh.uvList()[i][j] = std::complex<float>(u, v);
        }
        else // uv_size == 8
        {
          double u = (double &)triangles[row_size*i + order_size + 3*vertex_index_size + uv_order_size + 2*j*uv_size];
          double v = (double &)triangles[row_size*i + order_size + 3*vertex_index_size + uv_order_size + (2*j + 1)*uv_size];
          mesh.uvList()[i][j] = std::complex<float>((float)u, (float)v);
        }      
      }
    }
  }
  std::cout << "reading from " << file << ", " << mesh.indexList().size() << " triangles." << std::endl;
  return true;
}

bool convertCloud(const std::string &in_name, const std::string &out_name,
                  std::function<void(Eigen::Vector3d &start, Eigen::Vector3d &ends, double &time, RGBA &colour)> apply)
{
  std::ofstream ofs;
  if (!writeRayCloudChunkStart(out_name, ofs))
  {
    return false;
  }
  ray::RayPlyBuffer buffer;

  bool has_warned = false;
  // run the function 'apply' on each ray as it is read in, and write it out, one chunk at a time
  auto applyToChunk = [&apply, &buffer, &ofs, &has_warned](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                              std::vector<double> &times, std::vector<ray::RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); i++)
    {
      // We can adjust the applyToChunk arguments directly as they are non-const and their modification doesn't have
      // side effects
      apply(starts[i], ends[i], times[i], colours[i]);
    }
    ray::writeRayCloudChunk(ofs, buffer, starts, ends, times, colours, has_warned);
  };
  if (!ray::readPly(in_name, true, applyToChunk, 0))
  {
    return false;
  }
  ray::writeRayCloudChunkEnd(ofs);
  return true;
}

}  // namespace ray