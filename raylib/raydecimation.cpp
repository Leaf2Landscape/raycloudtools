// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raydecimation.h"
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include "raycloudwriter.h"
#include "rayparse.h"
#include "raylaz.h"
#include "raysysinfo.h"

namespace ray
{
// Open writer and pre-read extra-bytes VLR from a LAS/LAZ file so original sensor attributes
// are registered in the output before any points are written.
static bool beginWriter(CloudWriter &writer, const std::string &out_file,
                        const std::string &in_file, const std::string &ext,
                        uint16_t &pass_stride_out, std::vector<uint8_t> &extra_bytes_vlr_out,
                        bool &has_tree_id_out, bool &has_stem_id_out)
{
  has_tree_id_out = false;
  has_stem_id_out = false;
  pass_stride_out = 10;
  if (ext == "las" || ext == "laz")
  {
    LasHeader hdr;
    if (readLasHeader(in_file, hdr))
    {
      extra_bytes_vlr_out = hdr.sensorExtraVlr();
      pass_stride_out     = static_cast<uint16_t>(10 + hdr.sensorExtraSize());
      has_tree_id_out     = hdr.has("tree_id");
      has_stem_id_out     = hdr.has("stem_id");
    }
  }
  return writer.begin(out_file, extra_bytes_vlr_out, false, has_tree_id_out, has_stem_id_out);
}

// Read a LAS/LAZ file with per-point passthrough (and optionally tree/stem IDs), or fall back to
// Cloud::read for PLY. passthrough_buf is appended to before each callback, and cleared inside the
// callback after consuming the current chunk's bytes. tree_ids_out/stem_ids_out accumulate across
// chunks and must NOT be cleared per-chunk; callers track a base offset instead.
static bool readWithPassthrough(const std::string &file_name, const std::string &ext,
                                std::function<void(std::vector<Eigen::Vector3d> &,
                                                   std::vector<Eigen::Vector3d> &,
                                                   std::vector<double> &,
                                                   std::vector<RGBA> &)> apply,
                                std::vector<uint8_t> *passthrough_buf,
                                std::vector<int32_t> *tree_ids_out = nullptr,
                                std::vector<int32_t> *stem_ids_out = nullptr)
{
  if ((ext == "las" || ext == "laz") && passthrough_buf)
  {
    size_t num_bounded;
    return readLas(file_name, apply, num_bounded, 1.0, nullptr, computeReadChunkSize(),
                   tree_ids_out, passthrough_buf, nullptr, nullptr, stem_ids_out);
  }
  return Cloud::read(file_name, apply);
}

bool decimateSpatial(const std::string &file_name, double vox_width)
{
  const std::string stub = getFileNameStub(file_name);
  const std::string ext = getFileNameExtension(file_name);

  ray::CloudWriter writer;
  uint16_t pass_stride = 8;
  std::vector<uint8_t> extra_bytes_vlr;
  bool has_tree_id = false, has_stem_id = false;
  if (!beginWriter(writer, stub + "_decimated." + ext, file_name, ext, pass_stride, extra_bytes_vlr,
                   has_tree_id, has_stem_id))
    return false;

  ray::Cloud chunk;
  std::vector<int64_t> subsample;
  std::set<Eigen::Vector3i, ray::Vector3iLess> voxel_set;
  std::vector<uint8_t> passthrough_buf;
  std::vector<int32_t> tree_ids_buf, stem_ids_buf;
  size_t tree_ids_base = 0;

  auto decimate = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                      std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    double width = 0.01 * vox_width;
    subsample.clear();
    voxelSubsample(ends, width, subsample, voxel_set);
    chunk.resize(subsample.size());
    std::vector<uint8_t> chunk_pass;
    std::vector<int32_t> chunk_tree, chunk_stem;
    if (!passthrough_buf.empty())
    {
      chunk_pass.reserve(subsample.size() * pass_stride);
      for (int64_t i = 0; i < (int64_t)subsample.size(); i++)
      {
        const uint8_t *src = passthrough_buf.data() + subsample[i] * pass_stride;
        chunk_pass.insert(chunk_pass.end(), src, src + pass_stride);
      }
      passthrough_buf.clear();
    }
    for (int64_t i = 0; i < (int64_t)subsample.size(); i++)
    {
      int64_t id = subsample[i];
      chunk.starts[i] = starts[id];
      chunk.ends[i] = ends[id];
      chunk.colours[i] = colours[id];
      chunk.times[i] = times[id];
      if (!tree_ids_buf.empty())
        chunk_tree.push_back(tree_ids_buf[tree_ids_base + id]);
      if (!stem_ids_buf.empty())
        chunk_stem.push_back(stem_ids_buf[tree_ids_base + id]);
    }
    tree_ids_base += ends.size();
    writer.writeChunk(chunk.starts, chunk.ends, chunk.times, chunk.colours, chunk_pass, {}, chunk_tree, chunk_stem);
  };

  if (!readWithPassthrough(file_name, ext, decimate, &passthrough_buf,
                           has_tree_id ? &tree_ids_buf : nullptr,
                           has_stem_id ? &stem_ids_buf : nullptr))
    return false;
  writer.end();
  return true;
}

bool decimateTemporal(const std::string &file_name, int num_rays)
{
  const std::string stub = getFileNameStub(file_name);
  const std::string ext = getFileNameExtension(file_name);

  ray::CloudWriter writer;
  uint16_t pass_stride = 8;
  std::vector<uint8_t> extra_bytes_vlr;
  bool has_tree_id = false, has_stem_id = false;
  if (!beginWriter(writer, stub + "_decimated." + ext, file_name, ext, pass_stride, extra_bytes_vlr,
                   has_tree_id, has_stem_id))
    return false;

  ray::Cloud chunk;
  std::vector<uint8_t> passthrough_buf;
  std::vector<int32_t> tree_ids_buf, stem_ids_buf;
  size_t tree_ids_base = 0;

  auto decimate = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                      std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    size_t decimation = (size_t)num_rays;
    size_t count = (ends.size() + decimation - 1) / decimation;
    chunk.resize(count);
    std::vector<uint8_t> chunk_pass;
    std::vector<int32_t> chunk_tree, chunk_stem;
    if (!passthrough_buf.empty())
    {
      chunk_pass.reserve(count * pass_stride);
      for (size_t i = 0; i < ends.size(); i += decimation)
      {
        const uint8_t *src = passthrough_buf.data() + i * pass_stride;
        chunk_pass.insert(chunk_pass.end(), src, src + pass_stride);
      }
      passthrough_buf.clear();
    }
    for (size_t i = 0, c = 0; i < ends.size(); i += decimation, c++)
    {
      chunk.starts[c] = starts[i];
      chunk.ends[c] = ends[i];
      chunk.times[c] = times[i];
      chunk.colours[c] = colours[i];
      if (!tree_ids_buf.empty())
        chunk_tree.push_back(tree_ids_buf[tree_ids_base + i]);
      if (!stem_ids_buf.empty())
        chunk_stem.push_back(stem_ids_buf[tree_ids_base + i]);
    }
    tree_ids_base += ends.size();
    writer.writeChunk(chunk.starts, chunk.ends, chunk.times, chunk.colours, chunk_pass, {}, chunk_tree, chunk_stem);
  };

  if (!readWithPassthrough(file_name, ext, decimate, &passthrough_buf,
                           has_tree_id ? &tree_ids_buf : nullptr,
                           has_stem_id ? &stem_ids_buf : nullptr))
    return false;
  writer.end();
  return true;
}

bool decimateSpatioTemporal(const std::string &file_name, double vox_width, int num_rays)
{
  const std::string stub = getFileNameStub(file_name);
  const std::string ext = getFileNameExtension(file_name);

  ray::CloudWriter writer;
  uint16_t pass_stride = 8;
  std::vector<uint8_t> extra_bytes_vlr;
  bool has_tree_id = false, has_stem_id = false;
  if (!beginWriter(writer, stub + "_decimated." + ext, file_name, ext, pass_stride, extra_bytes_vlr,
                   has_tree_id, has_stem_id))
    return false;

  std::map<Eigen::Vector3i, Eigen::Vector2i, ray::Vector3iLess> voxel_map;
  std::vector<Eigen::Vector3i> samples;
  double voxel_width = 0.01 * vox_width;

  // Pass 1: count occupancy per voxel (no writing, no passthrough needed).
  auto count_voxels = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                          std::vector<double> &, std::vector<ray::RGBA> &)
  {
    for (size_t i = 0; i < ends.size(); i++)
    {
      Eigen::Vector3d coords = ends[i] / voxel_width;
      Eigen::Vector3i coordsi = Eigen::Vector3d(std::floor(coords[0]), std::floor(coords[1]),
                                                std::floor(coords[2])).cast<int>();
      auto found = voxel_map.find(coordsi);
      if (found == voxel_map.end())
      {
        voxel_map.insert({ coordsi, Eigen::Vector2i(1, 0) });
        samples.push_back(coordsi);
      }
      else
      {
        found->second[0]++;
      }
    }
  };
  if (!ray::Cloud::read(file_name, count_voxels))
    return false;

  for (auto &pos : samples)
  {
    int max_num = 0;
    for (int x = pos[0] - 1; x <= pos[0] + 1; x++)
      for (int y = pos[1] - 1; y <= pos[1] + 1; y++)
        for (int z = pos[2] - 1; z <= pos[2] + 1; z++)
        {
          auto found = voxel_map.find(Eigen::Vector3i(x, y, z));
          if (found != voxel_map.end())
            max_num = std::max(max_num, found->second[0]);
        }
    voxel_map.find(pos)->second[1] = max_num;
  }

  // Pass 2: emit selected points with passthrough preserved.
  std::vector<uint8_t> passthrough_buf;
  std::vector<int32_t> tree_ids_buf, stem_ids_buf;
  size_t tree_ids_base = 0;
  auto finalise = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                      std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    std::vector<Eigen::Vector3d> out_starts, out_ends;
    std::vector<double> out_times;
    std::vector<ray::RGBA> out_colours;
    std::vector<uint8_t> chunk_pass;
    std::vector<int32_t> chunk_tree, chunk_stem;
    for (size_t i = 0; i < ends.size(); i++)
    {
      Eigen::Vector3d coords = ends[i] / voxel_width;
      Eigen::Vector3i coordsi = Eigen::Vector3d(std::floor(coords[0]), std::floor(coords[1]),
                                                std::floor(coords[2])).cast<int>();
      auto found = voxel_map.find(coordsi);
      if (found == voxel_map.end())
        continue;
      int num = found->second[1];
      double segmentation = std::max(1.0, (double)num / (double)num_rays);
      int &ends_left = found->second[0];
      if (std::fmod((double)ends_left + 1.0, segmentation) <= std::fmod((double)ends_left, segmentation))
      {
        out_starts.push_back(starts[i]);
        out_ends.push_back(ends[i]);
        out_colours.push_back(colours[i]);
        out_times.push_back(times[i]);
        if (!passthrough_buf.empty())
        {
          const uint8_t *src = passthrough_buf.data() + i * pass_stride;
          chunk_pass.insert(chunk_pass.end(), src, src + pass_stride);
        }
        if (!tree_ids_buf.empty())
          chunk_tree.push_back(tree_ids_buf[tree_ids_base + i]);
        if (!stem_ids_buf.empty())
          chunk_stem.push_back(stem_ids_buf[tree_ids_base + i]);
      }
      ends_left--;
    }
    tree_ids_base += ends.size();
    passthrough_buf.clear();
    writer.writeChunk(out_starts, out_ends, out_times, out_colours, chunk_pass, {}, chunk_tree, chunk_stem);
  };
  if (!readWithPassthrough(file_name, ext, finalise, &passthrough_buf,
                           has_tree_id ? &tree_ids_buf : nullptr,
                           has_stem_id ? &stem_ids_buf : nullptr))
    return false;
  writer.end();
  return true;
}


bool decimateRaysSpatial(const std::string &file_name, double vox_width)
{
  const std::string stub = getFileNameStub(file_name);
  const std::string ext = getFileNameExtension(file_name);

  ray::CloudWriter writer;
  uint16_t pass_stride = 8;
  std::vector<uint8_t> extra_bytes_vlr;
  bool has_tree_id = false, has_stem_id = false;
  if (!beginWriter(writer, stub + "_decimated." + ext, file_name, ext, pass_stride, extra_bytes_vlr,
                   has_tree_id, has_stem_id))
    return false;

  ray::Cloud chunk;
  Subsampler subsampler;
  std::vector<uint8_t> passthrough_buf;
  std::vector<int32_t> tree_ids_buf, stem_ids_buf;
  size_t tree_ids_base = 0;

  auto decimate = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                      std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    double width = 0.01 * vox_width;
    subsampler.subsample.clear();
    for (int i = 0; i < (int)ends.size(); i++)
    {
      subsampler.index = i;
      #define END_FIRST // Testing on building.ply it finds more rays, so deemed to be more successful at filling space
      #if defined END_FIRST
      walkGrid(ends[i] / width, starts[i] / width, subsampler);
      #else
      walkGrid(starts[i] / width, ends[i] / width, subsampler);
      #endif
    }
    chunk.resize(subsampler.subsample.size());
    std::vector<uint8_t> chunk_pass;
    std::vector<int32_t> chunk_tree, chunk_stem;
    if (!passthrough_buf.empty())
    {
      chunk_pass.reserve(subsampler.subsample.size() * pass_stride);
      for (int64_t i = 0; i < (int64_t)subsampler.subsample.size(); i++)
      {
        const uint8_t *src = passthrough_buf.data() + subsampler.subsample[i] * pass_stride;
        chunk_pass.insert(chunk_pass.end(), src, src + pass_stride);
      }
      passthrough_buf.clear();
    }
    for (int64_t i = 0; i < (int64_t)subsampler.subsample.size(); i++)
    {
      int64_t id = subsampler.subsample[i];
      chunk.starts[i] = starts[id];
      chunk.ends[i] = ends[id];
      chunk.colours[i] = colours[id];
      chunk.times[i] = times[id];
      if (!tree_ids_buf.empty())
        chunk_tree.push_back(tree_ids_buf[tree_ids_base + id]);
      if (!stem_ids_buf.empty())
        chunk_stem.push_back(stem_ids_buf[tree_ids_base + id]);
    }
    tree_ids_base += ends.size();
    writer.writeChunk(chunk.starts, chunk.ends, chunk.times, chunk.colours, chunk_pass, {}, chunk_tree, chunk_stem);
  };

  if (!readWithPassthrough(file_name, ext, decimate, &passthrough_buf,
                           has_tree_id ? &tree_ids_buf : nullptr,
                           has_stem_id ? &stem_ids_buf : nullptr))
    return false;
  writer.end();
  return true;
}

bool decimateAngular(const std::string &file_name, double radius_per_length)
{
  const std::string stub = getFileNameStub(file_name);
  const std::string ext = getFileNameExtension(file_name);

  ray::CloudWriter writer;
  uint16_t pass_stride = 8;
  std::vector<uint8_t> extra_bytes_vlr;
  bool has_tree_id = false, has_stem_id = false;
  if (!beginWriter(writer, stub + "_decimated." + ext, file_name, ext, pass_stride, extra_bytes_vlr,
                   has_tree_id, has_stem_id))
    return false;

  int min_index = -20;
  int max_index = 50;
  std::vector<std::set<Eigen::Vector3i, ray::Vector3iLess>> voxel_sets(max_index + 1 - min_index);
  std::vector<std::set<Eigen::Vector3i, ray::Vector3iLess>> visiteds(max_index + 1 - min_index);
  std::vector<int> candidate_indices;
  const double root2 = std::sqrt(2.0);
  const double logroot2 = std::log(root2);
  std::vector<double> voxel_widths(voxel_sets.size());
  for (int i = 0; i < (int)voxel_widths.size(); i++)
    voxel_widths[i] = std::pow(root2, (double)(i + min_index));
  int index = -1;

  // Pass 1: identify candidate indices (no passthrough needed).
  auto identify = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                      std::vector<double> &, std::vector<ray::RGBA> &)
  {
    for (size_t i = 0; i < ends.size(); i++)
    {
      index++;
      double radius = (starts[i] - ends[i]).norm() * 0.01 * radius_per_length;
      int map_index = std::max(min_index, std::min((int)std::round(std::log(2.0 * radius) / logroot2), max_index));
      Eigen::Vector3d coords = ends[i] / voxel_widths[map_index - min_index];
      Eigen::Vector3i coordsi = Eigen::Vector3d(std::floor(coords[0]), std::floor(coords[1]),
                                                std::floor(coords[2])).cast<int>();
      int ind = map_index - min_index;
      if (visiteds[ind].find(coordsi) != visiteds[ind].end())
        continue;
      if (voxel_sets[ind].insert(coordsi).second)
      {
        candidate_indices.push_back(index);
        Eigen::Vector3i pos = coordsi;
        double scale = root2;
        pos = Eigen::Vector3d(std::floor((double)coordsi[0] / scale), std::floor((double)coordsi[1] / scale),
                              std::floor((double)coordsi[2] / scale)).cast<int>();
        ind++;
        while (ind < (int)visiteds.size() && visiteds[ind].insert(pos).second)
        {
          ind++;
          scale *= root2;
          pos = Eigen::Vector3d(std::floor((double)coordsi[0] / scale), std::floor((double)coordsi[1] / scale),
                                std::floor((double)coordsi[2] / scale)).cast<int>();
        }
      }
    }
  };
  if (!ray::Cloud::read(file_name, identify))
    return false;

  std::cout << "finalising" << std::endl;
  for (auto &map : voxel_sets)
    map.clear();
  index = -1;
  int head = 0;

  // Pass 2: emit kept points with passthrough preserved.
  std::vector<uint8_t> passthrough_buf;
  std::vector<int32_t> tree_ids_buf, stem_ids_buf;
  size_t tree_ids_base = 0;
  auto finalise = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                      std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    std::vector<Eigen::Vector3d> out_starts, out_ends;
    std::vector<double> out_times;
    std::vector<ray::RGBA> out_colours;
    std::vector<uint8_t> chunk_pass;
    std::vector<int32_t> chunk_tree, chunk_stem;
    for (size_t i = 0; i < ends.size(); i++)
    {
      index++;
      if (head >= (int)candidate_indices.size() || index != candidate_indices[head])
        continue;
      head++;
      double radius = (starts[i] - ends[i]).norm() * 0.01 * radius_per_length;
      int map_index = std::max(min_index, std::min((int)std::round(std::log(2.0 * radius) / logroot2), max_index));
      Eigen::Vector3d coords = ends[i] / voxel_widths[map_index - min_index];
      Eigen::Vector3i coordsi = Eigen::Vector3d(std::floor(coords[0]), std::floor(coords[1]),
                                                std::floor(coords[2])).cast<int>();
      int ind = map_index - min_index;
      if (visiteds[ind].find(coordsi) == visiteds[ind].end())
      {
        out_starts.push_back(starts[i]);
        out_ends.push_back(ends[i]);
        out_colours.push_back(colours[i]);
        out_times.push_back(times[i]);
        if (!passthrough_buf.empty())
        {
          const uint8_t *src = passthrough_buf.data() + i * pass_stride;
          chunk_pass.insert(chunk_pass.end(), src, src + pass_stride);
        }
        if (!tree_ids_buf.empty())
          chunk_tree.push_back(tree_ids_buf[tree_ids_base + i]);
        if (!stem_ids_buf.empty())
          chunk_stem.push_back(stem_ids_buf[tree_ids_base + i]);
      }
    }
    tree_ids_base += ends.size();
    passthrough_buf.clear();
    writer.writeChunk(out_starts, out_ends, out_times, out_colours, chunk_pass, {}, chunk_tree, chunk_stem);
  };
  if (!readWithPassthrough(file_name, ext, finalise, &passthrough_buf,
                           has_tree_id ? &tree_ids_buf : nullptr,
                           has_stem_id ? &stem_ids_buf : nullptr))
    return false;
  writer.end();
  return true;
}

namespace
{
// Decode one extra-byte sensor value from a passthrough slice into a double, using the LAS
// extra-byte type code (same table as raycombine.cpp's parseSensorAttrs). @c ptr points at the
// start of this point's passthrough record; the value lives at the 10-byte fixed prefix + offset.
// Returns NaN when all bytes at the field position are the 0xFF missing-field sentinel written by
// the combiner for inputs that lack this attribute. NaN propagates into beats() as a universal loser.
double decodeExtraByte(const uint8_t *ptr, const ResolvedTiebreaker &tb)
{
  const uint8_t *p = ptr + 10 + tb.byte_offset;
  bool all_ff = true;
  for (uint8_t b = 0; b < tb.byte_size; ++b)
    if (p[b] != 0xFF) { all_ff = false; break; }
  if (all_ff)
    return std::numeric_limits<double>::quiet_NaN();
  switch (tb.dtype)
  {
  case 1: { uint8_t v;  std::memcpy(&v, p, sizeof(v)); return static_cast<double>(v); }
  case 2: { int8_t v;   std::memcpy(&v, p, sizeof(v)); return static_cast<double>(v); }
  case 3: { uint16_t v; std::memcpy(&v, p, sizeof(v)); return static_cast<double>(v); }
  case 4: { int16_t v;  std::memcpy(&v, p, sizeof(v)); return static_cast<double>(v); }
  case 5: { uint32_t v; std::memcpy(&v, p, sizeof(v)); return static_cast<double>(v); }
  case 6: { int32_t v;  std::memcpy(&v, p, sizeof(v)); return static_cast<double>(v); }
  case 7: { uint64_t v; std::memcpy(&v, p, sizeof(v)); return static_cast<double>(v); }
  case 8: { int64_t v;  std::memcpy(&v, p, sizeof(v)); return static_cast<double>(v); }
  case 9: { float v;    std::memcpy(&v, p, sizeof(v)); return static_cast<double>(v); }
  case 10:{ double v;   std::memcpy(&v, p, sizeof(v)); return v; }
  default: return 0.0;
  }
}

// Resolve the field values for one point in spec order. When @c pass_ptr is null (no passthrough:
// PLY input) all criteria return NaN. NaN is also returned by decodeExtraByte when the extra-byte
// bytes are the 0xFF missing-field sentinel. beats() treats any NaN value as a universal loser.
void resolveValues(const std::vector<ResolvedTiebreaker> &spec, const Eigen::Vector3d &start,
                   const Eigen::Vector3d &end, double time, const RGBA &colour,
                   const uint8_t *pass_ptr, std::vector<double> &out)
{
  static const double kMissing = std::numeric_limits<double>::quiet_NaN();
  out.resize(spec.size());
  for (size_t s = 0; s < spec.size(); ++s)
  {
    switch (spec[s].kind)
    {
    case TiebreakKind::Reflectance: out[s] = static_cast<double>(colour.alpha); break;
    case TiebreakKind::Range:       out[s] = (end - start).norm(); break;
    case TiebreakKind::Time:        out[s] = time; break;
    case TiebreakKind::ExtraByte:
      out[s] = pass_ptr ? decodeExtraByte(pass_ptr, spec[s]) : kMissing; break;
    // Fixed 10-byte LAS passthrough fields (bytes 0-9, always present for LAS/LAZ inputs).
    case TiebreakKind::ReturnNumber:
      out[s] = pass_ptr ? static_cast<double>(pass_ptr[0] & 0x0Fu) : kMissing; break;
    case TiebreakKind::NumberOfReturns:
      out[s] = pass_ptr ? static_cast<double>((pass_ptr[0] >> 4) & 0x0Fu) : kMissing; break;
    case TiebreakKind::Classification:
      out[s] = pass_ptr ? static_cast<double>(pass_ptr[2]) : kMissing; break;
    case TiebreakKind::UserData:
      out[s] = pass_ptr ? static_cast<double>(pass_ptr[3]) : kMissing; break;
    case TiebreakKind::ScanAngle: {
      if (!pass_ptr) { out[s] = kMissing; break; }
      int16_t v; std::memcpy(&v, pass_ptr + 4, 2);
      out[s] = static_cast<double>(v); break;
    }
    case TiebreakKind::PointSourceId: {
      if (!pass_ptr) { out[s] = kMissing; break; }
      uint16_t v; std::memcpy(&v, pass_ptr + 6, 2);
      out[s] = static_cast<double>(v); break;
    }
    }
  }
}

// Lexicographic comparison of candidate values against the current winner. Returns true iff
// @c cand should replace @c best under @c spec (first differing field with the wrong order loses).
// NaN (missing-field sentinel) always loses: a NaN candidate never beats a real best, and a real
// candidate always beats a NaN best. Two NaN values tie and fall through to the next criterion.
bool beats(const std::vector<ResolvedTiebreaker> &spec, const std::vector<double> &cand,
           const std::vector<double> &best)
{
  for (size_t s = 0; s < spec.size(); ++s)
  {
    if (cand[s] == best[s])
      continue;
    const bool cand_nan = std::isnan(cand[s]);
    const bool best_nan = std::isnan(best[s]);
    if (cand_nan && best_nan) continue;
    if (cand_nan) return false;
    if (best_nan) return true;
    const bool cand_lower = cand[s] < best[s];
    return spec[s].ascending ? cand_lower : !cand_lower;
  }
  return false;  // equal across all criteria: keep the existing (earlier) winner
}

Eigen::Vector3i voxelKey(const Eigen::Vector3d &end, double vox_width)
{
  return Eigen::Vector3i(static_cast<int>(std::floor(end[0] / vox_width)),
                         static_cast<int>(std::floor(end[1] / vox_width)),
                         static_cast<int>(std::floor(end[2] / vox_width)));
}
}  // namespace

bool deduplicateVoxel(const std::string &file_name, double vox_width,
                      const std::vector<ResolvedTiebreaker> &spec)
{
  const std::string ext = getFileNameExtension(file_name);

  const bool is_las = (ext == "las" || ext == "laz");

  // Inspect the combined file's schema so the rewritten file preserves all of its columns:
  // sensor extra-bytes (passthrough), tree_id/stem_id labels, and native RGB.
  uint16_t pass_stride = 10;
  std::vector<uint8_t> extra_bytes_vlr;
  bool has_rgb = false;
  if (is_las)
  {
    LasHeader hdr;
    if (readLasHeader(file_name, hdr))
    {
      extra_bytes_vlr = hdr.sensorExtraVlr();
      pass_stride     = static_cast<uint16_t>(10 + hdr.sensorExtraSize());
      has_rgb         = hdr.has_rgb;
    }
  }

  // Per-chunk buffers shared by both passes. readLas appends tree_id/stem_id and passthrough across
  // chunks; the callbacks clear them each chunk, so each callback sees only its own chunk's slice.
  std::vector<uint8_t> passthrough_buf;
  std::vector<int32_t> tree_ids_buf, stem_ids_buf;
  // tree_id/stem_id are only populated by readLas when the file declares them; detect presence by
  // whether the first pass yielded any labels.
  bool has_tree_ids = false, has_stem_ids = false;

  // Drive @c apply over the file with labels + passthrough for LAS/LAZ, or plain Cloud::read for PLY.
  auto readAll = [&](std::function<void(std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &,
                                        std::vector<double> &, std::vector<RGBA> &)> apply) -> bool
  {
    if (is_las)
    {
      size_t num_bounded;
      return readLas(file_name, apply, num_bounded, 1.0, nullptr, computeReadChunkSize(),
                     &tree_ids_buf, &passthrough_buf, nullptr, nullptr, &stem_ids_buf);
    }
    return Cloud::read(file_name, apply);
  };

  // Pass 1: find the global winner point index per voxel cell.
  std::map<Eigen::Vector3i, int64_t, ray::Vector3iLess> winners;
  std::map<Eigen::Vector3i, std::vector<double>, ray::Vector3iLess> winner_vals;
  int64_t index = -1;
  std::vector<double> vals;
  auto find_winners = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                          std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    const size_t n_pts = ends.size();
    const bool have_pass = passthrough_buf.size() >= n_pts * pass_stride;
    if (!tree_ids_buf.empty()) has_tree_ids = true;
    if (!stem_ids_buf.empty()) has_stem_ids = true;
    for (size_t i = 0; i < n_pts; ++i)
    {
      ++index;
      const uint8_t *pass_ptr = have_pass ? passthrough_buf.data() + i * pass_stride : nullptr;
      resolveValues(spec, starts[i], ends[i], times[i], colours[i], pass_ptr, vals);
      const Eigen::Vector3i key = voxelKey(ends[i], vox_width);
      auto found = winners.find(key);
      if (found == winners.end())
      {
        winners.insert({ key, index });
        winner_vals.insert({ key, vals });
      }
      else if (beats(spec, vals, winner_vals[key]))
      {
        found->second   = index;
        winner_vals[key] = vals;
      }
    }
    passthrough_buf.clear();
    tree_ids_buf.clear();
    stem_ids_buf.clear();
  };
  if (!readAll(find_winners))
    return false;
  winner_vals.clear();

  // Pass 2: stream again, emitting only the winning point per voxel, preserving every column.
  // Keep the original extension on the temp file so CloudWriter infers the same format (LAS vs PLY).
  const std::string tmp_file = getFileNameStub(file_name) + "_dedup_tmp." + ext;
  ray::CloudWriter writer;
  if (!writer.begin(tmp_file, extra_bytes_vlr, /*with_beam_id=*/false,
                    /*with_tree_id=*/has_tree_ids, /*with_stem_id=*/has_stem_ids, /*with_rgb=*/has_rgb))
    return false;

  index = -1;
  auto emit_winners = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                          std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    const size_t n_pts = ends.size();
    const bool have_pass   = passthrough_buf.size() >= n_pts * pass_stride;
    const bool have_tree   = tree_ids_buf.size() >= n_pts;
    const bool have_stem   = stem_ids_buf.size() >= n_pts;
    std::vector<Eigen::Vector3d> out_starts, out_ends;
    std::vector<double> out_times;
    std::vector<ray::RGBA> out_colours;
    std::vector<uint8_t> chunk_pass;
    std::vector<int32_t> chunk_tree, chunk_stem;
    for (size_t i = 0; i < n_pts; ++i)
    {
      ++index;
      const Eigen::Vector3i key = voxelKey(ends[i], vox_width);
      auto found = winners.find(key);
      if (found == winners.end() || found->second != index)
        continue;
      out_starts.push_back(starts[i]);
      out_ends.push_back(ends[i]);
      out_times.push_back(times[i]);
      out_colours.push_back(colours[i]);
      if (have_pass)
      {
        const uint8_t *src = passthrough_buf.data() + i * pass_stride;
        chunk_pass.insert(chunk_pass.end(), src, src + pass_stride);
      }
      if (have_tree) chunk_tree.push_back(tree_ids_buf[i]);
      if (have_stem) chunk_stem.push_back(stem_ids_buf[i]);
    }
    passthrough_buf.clear();
    tree_ids_buf.clear();
    stem_ids_buf.clear();
    writer.writeChunk(out_starts, out_ends, out_times, out_colours, chunk_pass, {}, chunk_tree, chunk_stem);
  };
  if (!readAll(emit_winners))
    return false;
  writer.end();

  if (std::rename(tmp_file.c_str(), file_name.c_str()) != 0)
  {
    std::cerr << "Error: could not replace " << file_name << " with deduplicated output" << std::endl;
    return false;
  }
  return true;
}
}  // namespace ray
