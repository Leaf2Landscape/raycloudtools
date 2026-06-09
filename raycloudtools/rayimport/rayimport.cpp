// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>

#include "raylib/raycloud.h"
#include "raylib/raycloudwriter.h"
#include "raylib/raylaz.h"
#include "raylib/raysysinfo.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/rayutils.h"

#include "raylib/rayriegl.h"
#include "raylib/raytrajectory.h"

namespace
{
// Collect the comma-separated tokens of one --filters value, reading from a file if the value
// names an openable file, otherwise treating the value itself as the inline list.
void appendFilterTokens(const std::string &value, std::vector<std::string> &tokens)
{
  std::ifstream file(value);
  if (file.is_open())
  {
    std::string line;
    while (std::getline(file, line))
    {
      std::stringstream ss(line);
      std::string tok;
      while (std::getline(ss, tok, ','))
        if (!tok.empty())
          tokens.push_back(tok);
    }
  }
  else
  {
    std::stringstream ss(value);
    std::string tok;
    while (std::getline(ss, tok, ','))
      if (!tok.empty())
        tokens.push_back(tok);
  }
}

std::vector<ray::FieldFilter> parseFilterArgs(int argc, char *argv[])
{
  std::vector<std::string> tokens;
  for (int i = 1; i < argc; i++)
  {
    if ((std::strcmp(argv[i], "--filters") == 0 || std::strcmp(argv[i], "-f") == 0) && i + 1 < argc)
    {
      appendFilterTokens(argv[i + 1], tokens);
      ++i;
    }
  }

  std::vector<ray::FieldFilter> filters;
  if (tokens.size() % 3 != 0)
  {
    std::cout << "warning: --filters expects groups of 3 (field,min,max); dropping "
              << (tokens.size() % 3) << " leftover token(s)" << std::endl;
  }
  const size_t groups = tokens.size() / 3;
  for (size_t g = 0; g < groups; g++)
  {
    const std::string &name = tokens[g * 3 + 0];
    char *end_min = nullptr;
    char *end_max = nullptr;
    const double min_val = std::strtod(tokens[g * 3 + 1].c_str(), &end_min);
    const double max_val = std::strtod(tokens[g * 3 + 2].c_str(), &end_max);
    if (end_min == tokens[g * 3 + 1].c_str() || end_max == tokens[g * 3 + 2].c_str())
    {
      std::cout << "warning: --filters could not parse min/max for field '" << name << "'; skipping"
                << std::endl;
      continue;
    }
    ray::FieldFilter f;
    f.name = name;
    f.min_val = min_val;
    f.max_val = max_val;
    filters.push_back(f);
  }
  return filters;
}

std::vector<char *> stripFilterArgs(int argc, char *argv[])
{
  std::vector<char *> out;
  out.reserve(argc);
  for (int i = 0; i < argc; i++)
  {
    if (i > 0 && (std::strcmp(argv[i], "--filters") == 0 || std::strcmp(argv[i], "-f") == 0))
    {
      ++i;  // also skip the following value
      continue;
    }
    out.push_back(argv[i]);
  }
  return out;
}

void resolveFiltersLas(std::vector<ray::FieldFilter> &filters, const std::vector<uint8_t> &extra_bytes_vlr,
                       uint16_t /*orig_extra_size*/)
{
  struct StdField
  {
    const char *name;
    int offset;
    int size;
    bool is_signed;
    bool is_float;
    double scale;
  };
  static const StdField kStdFields[] = {
    { "classification", 2, 1, false, false, 1.0 },
    { "user_data", 3, 1, false, false, 1.0 },
    { "scan_angle", 4, 2, true, false, 0.006 },
    { "point_source_id", 6, 2, false, false, 1.0 },
    { "intensity", 8, 2, false, false, 1.0 },
  };
  // Mirror of kDecodeExtraTypeSize in raylasdecode.h: LAS EXTRA_BYTES data_type -> byte size.
  static const uint16_t kExtraTypeSize[11] = { 0, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };

  for (ray::FieldFilter &f : filters)
  {
    bool matched = false;
    for (const StdField &sf : kStdFields)
    {
      if (f.name == sf.name)
      {
        f.pass_offset = sf.offset;
        f.pass_size = sf.size;
        f.is_signed = sf.is_signed;
        f.is_float = sf.is_float;
        f.scale = sf.scale;
        f.resolved = true;
        matched = true;
        break;
      }
    }
    if (matched)
      continue;

    int byte_offset = 10;  // sensor extras follow the 10-byte standard passthrough block
    const size_t num_attrs = extra_bytes_vlr.size() / 192;
    for (size_t a = 0; a < num_attrs; a++)
    {
      const uint8_t *rec = extra_bytes_vlr.data() + a * 192;
      const uint8_t dtype = rec[2];
      const uint16_t attr_size = (dtype > 0 && dtype <= 10) ? kExtraTypeSize[dtype] : 0;
      if (attr_size == 0)
        continue;
      char attr_name[33] = {};
      std::memcpy(attr_name, rec + 4, 32);
      if (f.name == attr_name)
      {
        f.pass_offset = byte_offset;
        f.pass_size = attr_size;
        f.is_signed = (dtype == 2 || dtype == 4 || dtype == 6 || dtype == 8 || dtype == 9 || dtype == 10);
        f.is_float = (dtype == 9 || dtype == 10);
        // LAS EXTRA_BYTES VLR: options bit 3 = scale relevant, bit 4 = offset relevant.
        // scale at bytes 112-119 (first double), offset at bytes 136-143 (first double).
        const uint8_t opts = rec[3];
        f.scale  = (opts & 0x08u) ? *reinterpret_cast<const double *>(rec + 112) : 1.0;
        f.offset = (opts & 0x10u) ? *reinterpret_cast<const double *>(rec + 136) : 0.0;
        f.resolved = true;
        matched = true;
        break;
      }
      byte_offset += attr_size;
    }
    if (!matched)
      std::cout << "Filter field '" << f.name << "' not found in input — skipping" << std::endl;
  }
}

double readPassthroughField(const uint8_t *p, int byte_offset, int byte_size, bool is_signed, bool is_float)
{
  const uint8_t *src = p + byte_offset;
  if (is_float)
  {
    if (byte_size == 4)
    {
      float v;
      std::memcpy(&v, src, 4);
      return static_cast<double>(v);
    }
    double v;
    std::memcpy(&v, src, 8);
    return v;
  }
  if (is_signed)
  {
    if (byte_size == 1)
    {
      int8_t v;
      std::memcpy(&v, src, 1);
      return static_cast<double>(v);
    }
    if (byte_size == 2)
    {
      int16_t v;
      std::memcpy(&v, src, 2);
      return static_cast<double>(v);
    }
    if (byte_size == 4)
    {
      int32_t v;
      std::memcpy(&v, src, 4);
      return static_cast<double>(v);
    }
    int64_t v;
    std::memcpy(&v, src, 8);
    return static_cast<double>(v);
  }
  if (byte_size == 1)
  {
    uint8_t v;
    std::memcpy(&v, src, 1);
    return static_cast<double>(v);
  }
  if (byte_size == 2)
  {
    uint16_t v;
    std::memcpy(&v, src, 2);
    return static_cast<double>(v);
  }
  if (byte_size == 4)
  {
    uint32_t v;
    std::memcpy(&v, src, 4);
    return static_cast<double>(v);
  }
  uint64_t v;
  std::memcpy(&v, src, 8);
  return static_cast<double>(v);
}
}  // namespace

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Import a point cloud and trajectory file into a ray cloud" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayimport pointcloudfile trajectoryfile  - pointcloudfile can be a .laz, .las, .ply or .rxp file" << std::endl;
  std::cout << "                                           trajectoryfile is a text file using 'time x y z' format per line" << std::endl;
  std::cout << "                                           trajectoryfile may instead be a 4x4 transform matrix" << std::endl;
  std::cout << "                                           (16 values over 4 lines); auto-detected, or force with --transform/-t" << std::endl;
  std::cout << "rayimport pointcloudfile 0,0,0           - use 0,0,0 as the sensor location" << std::endl;
  std::cout << "rayimport pointcloudfile ray 0,0,-10     - use 0,0,-10 as the constant ray vector from start to point" << std::endl;
  std::cout << "                                          --max_intensity 100 - specify maximum intensity value (default: 65535 for .las/.laz, 100 otherwise)." << std::endl;
  std::cout << "                                                              0 sets all to full intensity (bounded rays)." << std::endl;
  std::cout << "                                        --remove_start_pos  - translate so first point is at 0,0,0" << std::endl;
  std::cout << "                                        --beam_id           - assign a per-pulse beam_id extra attribute" << std::endl;
  std::cout << "                                        --filters/-f \"field,min,max[,field2,min2,max2,...]\"" << std::endl;
  std::cout << "                                                            - keep only LAS/LAZ points whose field is in [min,max]." << std::endl;
  std::cout << "                                                              Groups of 3 comma-separated values, or a path to a text" << std::endl;
  std::cout << "                                                              file (one field,min,max per line). Standard fields:" << std::endl;
  std::cout << "                                                              intensity, classification, scan_angle (degrees)," << std::endl;
  std::cout << "                                                              user_data, point_source_id. Sensor extra fields by VLR name." << std::endl;
  std::cout << "rayimport pointcloudfile unbound transformfile - load unbound data (pulses that missed) from RIEGL .rxp file" << std::endl;
  std::cout << "                                               transformfile is a text file containing a 4x4 transformation matrix" << std::endl;
  std::cout << "The output is a _raycloud.las/.laz file (preserving .laz if the input is .laz)." << std::endl;
  // clang-format on
  exit(exit_code);
}


int rayImport(int argc, char *argv[])
{
  ray::DoubleArgument max_intensity(0.0, 1e8, 100.0);
  ray::Vector3dArgument position, ray_vec;
  ray::TextArgument ray_text("ray");
  ray::TextArgument unbound_text("unbound");
  ray::OptionalKeyValueArgument max_intensity_option("max_intensity", 'm', &max_intensity);
  ray::OptionalFlagArgument remove("remove_start_pos", 'r');
  ray::OptionalFlagArgument beam_id_opt("beam_id", 'b');
  ray::OptionalFlagArgument transform_flag("transform", 't');
  ray::FileArgument cloud_file, trajectory_file, transform_file;
  std::vector<ray::FieldFilter> filters = parseFilterArgs(argc, argv);
  std::vector<char *> filt_argv = stripFilterArgs(argc, argv);
  int filt_argc = static_cast<int>(filt_argv.size());
  bool standard_format =
    ray::parseCommandLine(filt_argc, filt_argv.data(), { &cloud_file, &trajectory_file }, { &max_intensity_option, &remove, &beam_id_opt, &transform_flag });
  bool position_format =
    ray::parseCommandLine(filt_argc, filt_argv.data(), { &cloud_file, &position }, { &max_intensity_option, &remove, &beam_id_opt });
  bool ray_format =
    ray::parseCommandLine(filt_argc, filt_argv.data(), { &cloud_file, &ray_text, &ray_vec }, { &max_intensity_option, &remove, &beam_id_opt });
  bool unbound_format = ray::parseCommandLine(filt_argc, filt_argv.data(), { &cloud_file, &unbound_text, &transform_file }, { &remove, &beam_id_opt });
  if (!standard_format && !position_format && !ray_format && !unbound_format)
    usage();

  ray::Cloud cloud;
  const std::string &traj_file = trajectory_file.name();
  const std::string &trans_file = transform_file.name();
  double maximum_intensity = max_intensity.value();
  if (!max_intensity_option.isSet())
  {
    const std::string ext = cloud_file.nameExt();
    maximum_intensity = (ext == "las" || ext == "laz") ? 65535.0 : 100.0;
  }

  // init transformation
  std::vector<double> transformation;
  Eigen::Matrix4d transform_matrix = Eigen::Matrix4d::Identity();
  bool transform_format = false;
  // load the trajectory first, it should fit into main memory

  ray::Trajectory trajectory;
  if (standard_format)
  {
    const std::string traj_end = traj_file.substr(traj_file.size() - 4);
    const bool is_cloud_traj = (traj_end == ".ply" || traj_end == ".las" || traj_end == ".laz");

    bool ambiguous = false;
    if (!is_cloud_traj)
    {
      if (transform_flag.isSet())
        transform_format = true;
      else if (ray::looksLikeTransformMatrix(traj_file, &ambiguous))
        transform_format = true;
    }

    if (transform_format)
    {
      std::cout << "parsing " << traj_file << " as transform matrix" << std::endl;
      const auto vals = ray::readNumericFile(traj_file);
      if (vals.size() != 16)
        usage();
      for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
          transform_matrix(r, c) = vals[r * 4 + c];
    }
    else
    {
      if (ambiguous)
        std::cout << "warning: " << traj_file
                  << " is 4x4 but its first column is monotonically increasing; "
                     "parsing as trajectory (pass --transform to force matrix)" << std::endl;
      std::cout << "parsing " << traj_file << " as trajectory" << std::endl;

      if (is_cloud_traj)
      {
        std::vector<Eigen::Vector3d> starts;
        std::vector<Eigen::Vector3d> ends;
        std::vector<double> times;
        std::vector<ray::RGBA> colours;
        if (traj_end == ".ply")
        {
          if (!ray::readPly(traj_file, starts, ends, times, colours, false))
            return false;
        }
        else
        {
          if (!ray::readLas(traj_file, ends, times, colours, maximum_intensity))
            return false;
        }
        trajectory.points() = std::move(ends);
        trajectory.times() = std::move(times);
      }
      else if (!trajectory.load(traj_file))
        usage();
    }
  }

  std::string save_file = cloud_file.nameStub() + "_raycloud";
  const std::string in_ext = cloud_file.nameExt();
  const std::string save_ext = (in_ext == "laz") ? "laz" : "las";
  size_t num_bounded = 0;
  uint8_t max_alpha_seen = 0;

  // Pre-read original sensor extra-byte attributes from the input LAS/LAZ header so the writer
  // can register and preserve them before opening the output file.
  std::vector<uint8_t> input_extra_bytes_vlr;
  bool input_has_rgb = false;
  uint16_t orig_extra = 0;
  if (in_ext == "laz" || in_ext == "las")
  {
    ray::readLasExtraBytesVlr(cloud_file.name(), orig_extra, input_extra_bytes_vlr, nullptr, &input_has_rgb);
    resolveFiltersLas(filters, input_extra_bytes_vlr, orig_extra);
  }

  // Pre-scan: build a global GPS-time -> beam_id map so that all returns of one pulse
  // (same timestamp) get the same beam_id regardless of chunk boundaries or spatial sorting.
  std::unordered_map<double, int32_t> beam_id_map;
  if (beam_id_opt.isSet() && (in_ext == "las" || in_ext == "laz"))
  {
    int32_t next_id = 0;
    size_t dummy_bounded = 0;
    ray::readLas(cloud_file.name(), [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &,
                                        std::vector<double> &scan_times, std::vector<ray::RGBA> &) {
      for (const double t : scan_times)
        if (beam_id_map.emplace(t, next_id).second)
          ++next_id;
    }, dummy_bounded, maximum_intensity, nullptr, ray::computeReadChunkSize(), nullptr, nullptr);
  }

  ray::CloudWriter writer;
  if (!writer.begin(save_file + "." + save_ext, input_extra_bytes_vlr, beam_id_opt.isSet(), false, false, input_has_rgb))
    usage();
  Eigen::Vector3d start_pos(0, 0, 0);
  bool first_chunk_done = false;
  double min_time = std::numeric_limits<double>::max();
  double max_time = std::numeric_limits<double>::lowest();
  std::vector<uint8_t> all_passthrough;
  size_t prev_pass_size = 0;
  int32_t current_beam_id = -1;
  double last_beam_time = std::numeric_limits<double>::quiet_NaN();
  Eigen::Vector3d last_beam_start(0, 0, 0);
  bool warned_beam_fallback = false;
  std::vector<int32_t> chunk_beam_ids;
  auto add_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                       std::vector<double> &times, std::vector<ray::RGBA> &colours) {
    // Capture the first point position once, used only by --remove_start_pos.
    if (!first_chunk_done)
    {
      start_pos = ends[0];
      first_chunk_done = true;
    }
    // user provides a single sensor location (e.g. for static scanners)
    if (position_format)
    {
      starts = ends;
      Eigen::Vector3d pos = position.value();
      for (auto &start : starts)
      {
        start = pos;
      }
    }
    // user provides a constant ray vector
    // e.g. for an overhead aerial scan, if no trajectory is available
    else if (ray_format)
    {
      starts = ends;
      Eigen::Vector3d offset = -ray_vec.value();
      for (auto &start : starts)
      {
        start += offset;
      }
    }
    // unbound data: starts are already set correctly by readRXP with the transformation applied
    else if (unbound_format)
    {
    }
    // a 4x4 rigid transform: move points from scanner frame to world frame
    else if (transform_format)
    {
      // sensor origin = translation column; points transformed from scanner to world frame
      const Eigen::Matrix3d R = transform_matrix.block<3, 3>(0, 0);
      const Eigen::Vector3d t = transform_matrix.block<3, 1>(0, 3);
      starts.resize(ends.size());
      for (size_t i = 0; i < ends.size(); i++)
      {
        ends[i]   = R * ends[i] + t;
        starts[i] = t;
      }
    }
    // otherwise, a trajectory has been passed in
    else
    {
      // find the corresponding sensor locations for each point in the cloud
      trajectory.calculateStartPoints(times, starts);
      for (size_t i = 0; i < colours.size(); i++)
      {
        min_time = std::min(min_time, times[i]);
        max_time = std::max(max_time, times[i]);
        if (colours[i].alpha == 0 && ends[i][2] < starts[i][2])  // a nonreturn, we need to remove downward ones
        {
          Eigen::Vector3d dir = (ends[i] - starts[i]).normalized();
          const double minimal_distance_for_nonreturns = 0.1;
          ends[i] = starts[i] + dir * minimal_distance_for_nonreturns;
        }
      }
    }
    // option to remove the start position, for data that is in a global frame
    // this is particularly useful if we are storing the ray cloud positions using floats
    if (remove.isSet())
    {
      for (auto &end : ends)
      {
        end -= start_pos;
      }
      for (auto &start : starts)
      {
        start -= start_pos;
      }
    }
    for (const auto &c : colours)
      max_alpha_seen = std::max(max_alpha_seen, c.alpha);
    if (maximum_intensity == 0.0)
    {
      for (auto &c : colours)
      {
        c.alpha = 255;
      }
    }
    if (unbound_format)
    {
      for (auto &c : colours)
      {
        c.alpha = uint8_t(0);
      }
    }
    std::vector<uint8_t> chunk_pass(all_passthrough.begin() + prev_pass_size, all_passthrough.end());
    prev_pass_size = all_passthrough.size();
    chunk_beam_ids.clear();
    if (beam_id_opt.isSet())
    {
      chunk_beam_ids.resize(times.size());
      if (!beam_id_map.empty())
      {
        // LAS/LAZ: look up beam_id from the globally pre-assigned map keyed by GPS time.
        // All returns of a pulse share one timestamp and thus one beam_id, even when the
        // file is spatially sorted and returns of the same pulse span different chunks.
        for (size_t i = 0; i < times.size(); i++)
        {
          const auto it = beam_id_map.find(times[i]);
          chunk_beam_ids[i] = (it != beam_id_map.end()) ? it->second : -1;
        }
      }
      else
      {
        // PLY/RXP: detect new beam by ray-origin change or GPS-time change.
        for (size_t i = 0; i < times.size(); i++)
        {
          bool new_beam = false;
          if (!starts.empty())
          {
            const Eigen::Vector3d &prev = (i == 0) ? last_beam_start : starts[i - 1];
            new_beam = !starts[i].isApprox(prev, 0.005);
          }
          else
          {
            if (times[i] != 0.0 || i > 0)
            {
              new_beam = (times[i] != last_beam_time);
            }
            else
            {
              if (!warned_beam_fallback)
              {
                std::cout << "warning: no sensor positions or GPS timestamps detected; "
                             "beam_id will be one per point" << std::endl;
                warned_beam_fallback = true;
              }
              new_beam = true;
            }
          }
          if (new_beam)
          {
            ++current_beam_id;
            last_beam_time = times[i];
            if (!starts.empty())
              last_beam_start = starts[i];
          }
          chunk_beam_ids[i] = current_beam_id;
        }
      }
    }
    const size_t pass_stride = 10u + orig_extra;
    if (!filters.empty() && !chunk_pass.empty())
    {
      std::vector<size_t> keep;
      keep.reserve(ends.size());
      for (size_t i = 0; i < ends.size(); i++)
      {
        bool accept = true;
        for (const auto &f : filters)
        {
          if (!f.resolved)
            continue;
          double val = readPassthroughField(chunk_pass.data() + i * pass_stride, f.pass_offset, f.pass_size,
                                            f.is_signed, f.is_float);
          val = val * f.scale + f.offset;
          if (val < f.min_val || val > f.max_val)
          {
            accept = false;
            break;
          }
        }
        if (accept)
          keep.push_back(i);
      }
      if (keep.size() < ends.size())
      {
        auto compact_vec = [&keep](auto &vec) {
          std::remove_reference_t<decltype(vec)> tmp;
          tmp.reserve(keep.size());
          for (size_t idx : keep)
            tmp.push_back(std::move(vec[idx]));
          vec = std::move(tmp);
        };
        compact_vec(starts);
        compact_vec(ends);
        compact_vec(times);
        compact_vec(colours);
        if (!chunk_beam_ids.empty())
          compact_vec(chunk_beam_ids);
        std::vector<uint8_t> tmp_pass;
        tmp_pass.reserve(keep.size() * pass_stride);
        for (size_t idx : keep)
          tmp_pass.insert(tmp_pass.end(), chunk_pass.begin() + idx * pass_stride,
                          chunk_pass.begin() + idx * pass_stride + pass_stride);
        chunk_pass = std::move(tmp_pass);
      }
    }
    if (!writer.writeChunk(starts, ends, times, colours, chunk_pass, chunk_beam_ids))
      usage();
  };

  if (unbound_format)
  {
    transformation = ray::readNumericFile(trans_file);
    if (transformation.empty())
      usage();
  }

  if (!unbound_format)
    std::cout << "max_intensity: " << maximum_intensity << std::endl;
  if (!filters.empty() && in_ext != "las" && in_ext != "laz")
    std::cout << "warning: --filters is only supported for LAS/LAZ input; skipping filters" << std::endl;
  Eigen::Vector3d *offset = remove.isSet() ? &start_pos : nullptr;
  if (cloud_file.nameExt() == "ply")
  {
    bool can_times_be_missing = position_format || ray_format;
    if (!ray::readPly(cloud_file.name(), false, add_chunk, maximum_intensity,
                      can_times_be_missing))  // special case of reading a non-ray-cloud ply
    {
      usage();
    }
  }
  else if (cloud_file.nameExt() == "laz" || cloud_file.nameExt() == "las")
  {
    if (!ray::readLas(cloud_file.name(), add_chunk, num_bounded, maximum_intensity, offset, ray::computeReadChunkSize(), nullptr,
                      &all_passthrough))
    {
      usage();
    }
  }
  else if (cloud_file.nameExt() == "rxp")
  {
    if (!ray::readRXP(cloud_file.name(), add_chunk, num_bounded, maximum_intensity, transformation))
    {
      usage();
    }
  }
  else
  {
    std::cout << "Error converting unknown type: " << cloud_file.name() << std::endl;
    usage();
  }
  if (standard_format && !transform_format)
  {
    const float grace_period = 30.0;
    if (trajectory.times()[0] < min_time - grace_period)
    {
      std::cout << "trajectory begins " << min_time - trajectory.times()[0] << " s before first point cloud time"
                << std::endl;
    }
    if (trajectory.times().back() > max_time + grace_period)
    {
      std::cout << "trajectory ends " << trajectory.times().back() - max_time << " s after last point cloud time"
                << std::endl;
    }
    if (min_time < trajectory.times()[0] - grace_period || max_time > trajectory.times().back() + grace_period ||
        min_time > trajectory.times().back() || max_time < trajectory.times()[0])
    {
      std::cerr.precision(10);
      std::cerr << "Error: trajectory times " << trajectory.times()[0] << "-" << trajectory.times().back()
                << " do not span the point cloud times " << min_time << "-" << max_time << std::endl;
      usage();
    }
  }
  if (num_bounded == 0 && maximum_intensity > 0 && !unbound_format)
  {
    std::cout << "warning: all point cloud intensities are 0." << std::endl;
    std::cout << "If your sensor lacks intensity information, set them to full using:" << std::endl;
    std::cout << "rayimport <point cloud> <trajectory file> --max_intensity 0" << std::endl;
  }
  if (!unbound_format && maximum_intensity > 0.0 && max_alpha_seen > 0)
  {
    if (max_alpha_seen == 255)
      std::cout << "warning: intensity values hit or exceeded max_intensity (" << maximum_intensity
                << "); data may be clipped. Consider increasing --max_intensity." << std::endl;
    else if (max_alpha_seen < 26)
      std::cout << "warning: peak intensity is <10% of max_intensity (" << maximum_intensity
                << "); consider reducing --max_intensity for better precision." << std::endl;
  }
  writer.end();
  // if we remove the start position, then it is useful to print this value that is removed
  // so that the user hasn't lost information
  if (remove.isSet())
  {
 //   std::cout.precision(10);
    std::cout << "start position: " << std::setprecision(5) << std::fixed << start_pos.transpose() << " removed from all points" << std::endl;
  }
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayImport, argc, argv);
}