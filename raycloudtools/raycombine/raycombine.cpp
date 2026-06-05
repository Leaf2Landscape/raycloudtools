// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymerger.h"
#include "raylib/raymesh.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/rayprogressthread.h"
#include "raylib/raythreads.h"
#include "raylib/raycloudwriter.h"
#include "raylib/raylaz.h"
#include "raylib/raysysinfo.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Combines multiple ray clouds. Clouds are not moved but rays are omitted in the combined cloud according to the merge type specified." << std::endl;
  std::cout << "Outputs the combined cloud and the residual cloud of differences." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raycombine all raycloud1 raycloud2 ... raycloudN   - concatenate all the rays in the _combined.ply cloud ('all' is optional)" << std::endl;
  std::cout << "           min raycloud1 ... raycloudN 20 rays - combines into one cloud with minimal objects at differences" << std::endl;
  std::cout << "                                                 20 is the number of pass through rays to define " << std::endl;
  std::cout << "           max    - maximal objects included. This is a form of volume intersection (rather than min: union)." << std::endl;
  std::cout << "           oldest - keeps the oldest geometry when there is a difference in later ray clouds." << std::endl;
  std::cout << "           newest - uses the newest geometry when there is a difference in newer ray clouds." << std::endl;
  std::cout << "           order  - conflicts are resolved in argument order, with the first taking priority." << std::endl;
  std::cout << "raycombine basecloud min raycloud1 raycloud2 20 rays - 3-way merge, choses the changed geometry (from basecloud) at any differences. " << std::endl;
  std::cout << "                                                       For merge conflicts it uses the specified merge type." << std::endl;
  std::cout << "        --output raycloud_combined.ply               - optionally specify the output file name." << std::endl;
  // clang-format on
  exit(exit_code);
}

// Combines multiple clouds together
int rayCombine(int argc, char *argv[])
{
  ray::KeyChoice merge_type({ "min", "max", "oldest", "newest", "order" });
  ray::FileArgumentList cloud_files(2);
  ray::DoubleArgument num_rays(0.0, 100.0);
  ray::TextArgument rays_text("rays"), all_text("all");

  // Below: false = allow unusual file extensions, for auto-merging, which occurs on non-standard temporary file names
  ray::FileArgument base_cloud(false), cloud_1(false), cloud_2(false), output_file(false);
  ray::OptionalKeyValueArgument output("output", 'o', &output_file);

  // three-way merge option
  bool standard_format = ray::parseCommandLine(argc, argv, { &merge_type, &cloud_files, &num_rays, &rays_text }, { &output });
  bool concatenate_all = ray::parseCommandLine(argc, argv, { &all_text, &cloud_files }, { &output });
  bool threeway = ray::parseCommandLine(
    argc, argv, { &base_cloud, &merge_type, &cloud_1, &cloud_2, &num_rays, &rays_text }, { &output });
  bool threeway_concatenate =
    ray::parseCommandLine(argc, argv, { &base_cloud, &all_text, &cloud_1, &cloud_2 }, { &output });
  if (!standard_format && !concatenate_all && !threeway && !threeway_concatenate)
  {
    concatenate_all = ray::parseCommandLine(argc, argv, { &cloud_files }, { &output }); // a bit more ambiguous, so only try if the other formats failed
    if (!concatenate_all)
    {
      usage();
    }
  }

  // we know there is at least one file, as we specified a minimum number in FileArgumentList
  std::string file_stub =
    (threeway || threeway_concatenate) ? base_cloud.nameStub() : cloud_files.files()[0].nameStub();

  std::vector<ray::Cloud> clouds;
  if (threeway || threeway_concatenate)
  {
    clouds.resize(2);
    if (!clouds[0].load(cloud_1.name(), false))
      usage();
    if (!clouds[1].load(cloud_2.name(), false))
      usage();
  }
  else if (!concatenate_all)
  {
    clouds.resize(cloud_files.files().size());
    for (int i = 0; i < (int)cloud_files.files().size(); i++)
      if (!clouds[i].load(cloud_files.files()[i].name()))
        usage();
  }

  ray::Threads::init();
  ray::MergerConfig config;
  config.voxel_size = 0.0;  // Infer voxel size
  config.num_rays_filter_threshold = num_rays.value();
  config.merge_type = ray::MergeType::Mininum;

  if (merge_type.selectedKey() == "order")
  {
    config.merge_type = ray::MergeType::Order;
  }
  if (merge_type.selectedKey() == "oldest")
  {
    config.merge_type = ray::MergeType::Oldest;
  }
  if (merge_type.selectedKey() == "newest")
  {
    config.merge_type = ray::MergeType::Newest;
  }
  if (merge_type.selectedKey() == "min")
  {
    config.merge_type = ray::MergeType::Mininum;
  }
  if (merge_type.selectedKey() == "max")
  {
    config.merge_type = ray::MergeType::Maximum;
  }
  if (threeway_concatenate || concatenate_all)
  {
    config.merge_type = ray::MergeType::All;
  }
  const std::string combine_ext = ray::getFileNameExtension(
    (threeway || threeway_concatenate) ? base_cloud.name() : cloud_files.files()[0].name());
  std::string combined_file = output.isSet() ? output_file.name() : file_stub + "_combined." + combine_ext;
  if (concatenate_all)
  {
    // Each sensor-extra attribute declared in a file's EXTRA_BYTES VLR.
    struct SensorAttr
    {
      std::string name;
      uint16_t size   = 0;
      uint16_t offset = 0;  ///< cumulative byte offset within this file's sensor-extras slice
      std::array<uint8_t, 192> record{};  ///< raw 192-byte EXTRA_BYTES VLR record
    };

    // Parse a stripped EXTRA_BYTES VLR payload (as returned by readLasExtraBytesVlr) into
    // named attributes with cumulative byte offsets.
    auto parseSensorAttrs = [](const std::vector<uint8_t> &vlr)
    {
      constexpr uint16_t kTypeSize[11] = { 0, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };
      std::vector<SensorAttr> result;
      uint16_t off = 0;
      const int n = static_cast<int>(vlr.size()) / 192;
      for (int a = 0; a < n; ++a)
      {
        const uint8_t *rec = vlr.data() + a * 192;
        const uint8_t dtype = rec[2];
        const uint16_t sz = (dtype > 0 && dtype <= 10) ? kTypeSize[dtype] : 0;
        if (sz == 0)
          continue;
        SensorAttr attr;
        char name_buf[33] = {};
        std::memcpy(name_buf, rec + 4, 32);
        attr.name   = name_buf;
        attr.size   = sz;
        attr.offset = off;
        std::memcpy(attr.record.data(), rec, 192);
        result.push_back(std::move(attr));
        off += sz;
      }
      return result;
    };

    // Pre-pass: read the EXTRA_BYTES VLR from every LAS/LAZ input file.
    struct FileSchema
    {
      uint16_t orig_extra = 0;
      std::vector<SensorAttr> attrs;
    };
    const int nfiles = static_cast<int>(cloud_files.files().size());
    std::vector<FileSchema> schemas(nfiles);
    for (int f = 0; f < nfiles; ++f)
    {
      const std::string &fn = cloud_files.files()[f].name();
      const std::string fe  = ray::getFileNameExtension(fn);
      if (fe == "las" || fe == "laz")
      {
        std::vector<uint8_t> file_vlr;
        ray::readLasExtraBytesVlr(fn, schemas[f].orig_extra, file_vlr);
        schemas[f].attrs = parseSensorAttrs(file_vlr);
      }
    }

    // Build the union sensor-attr schema: name-deduplicated, ordered by first appearance.
    // Any attr absent in a particular file will have its union-layout slot filled with 0xFF
    // (the -1 sentinel for signed interpretations) during passthrough reformatting below.
    std::vector<SensorAttr> union_attrs;
    std::vector<uint8_t> union_vlr;
    uint16_t union_sensor_size = 0;
    for (const auto &schema : schemas)
    {
      for (const auto &attr : schema.attrs)
      {
        bool already = false;
        for (const auto &ua : union_attrs)
          if (ua.name == attr.name) { already = true; break; }
        if (!already)
        {
          SensorAttr ua = attr;
          ua.offset     = union_sensor_size;
          union_vlr.insert(union_vlr.end(), ua.record.begin(), ua.record.end());
          union_sensor_size += ua.size;
          union_attrs.push_back(std::move(ua));
        }
      }
    }

    // Per-file remap: for each union attr present in this file, record src_offset → dst_offset.
    // Attrs absent in a file are left at the 0xFF sentinel initialised during reformatting.
    struct AttrRemap { uint16_t src_off, dst_off, size; };
    const uint16_t writer_pass_stride = static_cast<uint16_t>(10 + union_sensor_size);
    std::vector<std::vector<AttrRemap>> remaps(nfiles);
    std::vector<bool> needs_remap(nfiles, false);
    for (int f = 0; f < nfiles; ++f)
    {
      for (const auto &ua : union_attrs)
        for (const auto &fa : schemas[f].attrs)
          if (fa.name == ua.name)
          {
            remaps[f].push_back({ fa.offset, ua.offset, std::min(fa.size, ua.size) });
            break;
          }
      const uint16_t fp = static_cast<uint16_t>(10 + schemas[f].orig_extra);
      needs_remap[f] = (fp != writer_pass_stride) || (remaps[f].size() < union_attrs.size());
      if (!needs_remap[f])
        for (const auto &r : remaps[f])
          if (r.src_off != r.dst_off) { needs_remap[f] = true; break; }
    }

    ray::CloudWriter writer;
    if (!writer.begin(combined_file, union_vlr))
      usage();

    for (int i = 0; i < nfiles; ++i)
    {
      const std::string &fname = cloud_files.files()[i].name();
      const std::string fext   = ray::getFileNameExtension(fname);
      const uint16_t file_pass_stride = static_cast<uint16_t>(10 + schemas[i].orig_extra);
      std::vector<uint8_t> passthrough_buf;
      size_t passthrough_cursor = 0;  ///< byte offset into passthrough_buf for the next chunk

      auto concatenate = [&, i, file_pass_stride](
          std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
          std::vector<double> &times, std::vector<ray::RGBA> &colours)
      {
        // Slice this chunk's passthrough at the cursor position.
        // Sequential path: passthrough_buf grows per chunk; cursor stays aligned.
        // laz-perf path: passthrough_buf is pre-allocated for all N points; cursor advances
        //                through it so each chunk gets the correct absolute slice.
        const size_t n_pts       = starts.size();
        const size_t chunk_bytes = n_pts * file_pass_stride;
        std::vector<uint8_t> chunk_pass;
        if (passthrough_buf.size() >= passthrough_cursor + chunk_bytes)
          chunk_pass.assign(
            passthrough_buf.begin() + static_cast<ptrdiff_t>(passthrough_cursor),
            passthrough_buf.begin() + static_cast<ptrdiff_t>(passthrough_cursor + chunk_bytes));
        passthrough_cursor += chunk_bytes;

        // Reformat sensor-extra bytes from this file's VLR layout to the union layout.
        // Missing attrs keep their 0xFF sentinel; matching attrs are mapped by name.
        if (needs_remap[i] && !chunk_pass.empty())
        {
          const size_t n = chunk_pass.size() / file_pass_stride;
          std::vector<uint8_t> out(n * writer_pass_stride, 0xFFu);
          for (size_t j = 0; j < n; ++j)
          {
            // Preserve the 10 fixed LAS bytes (return number, scan angle, intensity, …).
            std::memcpy(out.data()       + j * writer_pass_stride,
                        chunk_pass.data() + j * file_pass_stride, 10);
            // Map each named sensor attr from this file's layout to the union layout.
            for (const auto &r : remaps[i])
              std::memcpy(out.data()       + j * writer_pass_stride + 10 + r.dst_off,
                          chunk_pass.data() + j * file_pass_stride  + 10 + r.src_off, r.size);
          }
          chunk_pass = std::move(out);
        }
        writer.writeChunk(starts, ends, times, colours, chunk_pass);
      };

      if (fext == "las" || fext == "laz")
      {
        size_t num_bounded;
        if (!ray::readLas(fname, concatenate, num_bounded, 1.0, nullptr,
                          ray::computeReadChunkSize(), nullptr, &passthrough_buf))
          usage();
      }
      else
      {
        if (!ray::Cloud::read(fname, concatenate))
          usage();
      }
    }
    writer.end();
    return 0;
  }

  ray::Merger merger(config);
  ray::Progress progress;
  ray::ProgressThread progress_thread(progress);
  ray::Cloud concatenated_cloud;
  const ray::Cloud *fixed_cloud = &merger.fixedCloud();

  if (threeway || threeway_concatenate)
  {
    ray::Cloud base_cloud;
    if (!base_cloud.load(argv[1], false))
      usage();
    merger.mergeThreeWay(base_cloud, clouds[0], clouds[1], &progress);
  }
  else
  {
    merger.mergeMultiple(clouds, &progress);
    std::cout << merger.differenceCloud().rayCount() << " transients, " << merger.fixedCloud().rayCount()
              << " fixed rays." << std::endl;
    merger.differenceCloud().save(file_stub + "_differences." + combine_ext);
  }

  progress_thread.join();
  fixed_cloud->save(combined_file);
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayCombine, argc, argv);
}
