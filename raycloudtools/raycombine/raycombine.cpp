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
#include "raylib/raycloudreader.h"
#include "raylib/raydecimation.h"
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
  std::cout << "        --dedup 0.02                                 - 'all' mode only: deduplicate to one point per 0.02 m voxel cell (metres; first-seen wins)." << std::endl;
  std::cout << "        --dedup 0.02 --tiebreak \">reflectance,<range\" - 'all' mode only: keep the best point per voxel by field priority ('>' prefers higher, '<' prefers lower)." << std::endl;
  std::cout << "                    built-in fields: reflectance, range, time" << std::endl;
  std::cout << "                    fixed LAS fields: return_number, number_of_returns, classification, user_data, scan_angle, point_source_id" << std::endl;
  std::cout << "                    any named sensor extra-byte field from any input file (files lacking it use 0 as sentinel)." << std::endl;
  // clang-format on
  exit(exit_code);
}

// One parsed --tiebreak token: a field name and its preferred sort direction.
struct TiebreakToken
{
  std::string name;
  bool ascending;  ///< '<' prefers lower (ascending), '>' prefers higher
};

// Parse a --tiebreak spec like ">reflectance,<range" into ordered tokens. Each comma-separated
// token must start with '>' or '<'; anything malformed (empty token, missing/invalid prefix,
// empty name) is a hard error that prints usage and exits.
std::vector<TiebreakToken> parseTiebreakTokens(const std::string &spec)
{
  std::vector<TiebreakToken> tokens;
  size_t start = 0;
  while (start <= spec.size())
  {
    const size_t comma = spec.find(',', start);
    const std::string tok = spec.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
    if (tok.size() < 2 || (tok[0] != '>' && tok[0] != '<'))
    {
      std::cerr << "Malformed --tiebreak token: '" << tok << "' (expected '>name' or '<name')" << std::endl;
      usage();
    }
    tokens.push_back({ tok.substr(1), tok[0] == '<' });
    if (comma == std::string::npos)
      break;
    start = comma + 1;
  }
  return tokens;
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

  // 'all' mode voxel deduplication. --dedup is the voxel width in metres; --tiebreak (requires
  // --dedup) selects the best point per voxel by a comma-separated field-priority spec.
  ray::DoubleArgument dedup_width(0.0, 1e10);
  ray::StringArgument tiebreak_spec;
  ray::OptionalKeyValueArgument dedup_option("dedup", '\0', &dedup_width);
  ray::OptionalKeyValueArgument tiebreak_option("tiebreak", '\0', &tiebreak_spec);

  // three-way merge option
  bool standard_format = ray::parseCommandLine(argc, argv, { &merge_type, &cloud_files, &num_rays, &rays_text }, { &output });
  bool concatenate_all = ray::parseCommandLine(argc, argv, { &all_text, &cloud_files }, { &output, &dedup_option, &tiebreak_option });
  bool threeway = ray::parseCommandLine(
    argc, argv, { &base_cloud, &merge_type, &cloud_1, &cloud_2, &num_rays, &rays_text }, { &output });
  bool threeway_concatenate =
    ray::parseCommandLine(argc, argv, { &base_cloud, &all_text, &cloud_1, &cloud_2 }, { &output });
  if (!standard_format && !concatenate_all && !threeway && !threeway_concatenate)
  {
    concatenate_all = ray::parseCommandLine(argc, argv, { &cloud_files }, { &output, &dedup_option, &tiebreak_option }); // a bit more ambiguous, so only try if the other formats failed
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

    // Pre-pass: read the full header from every input file via CloudReader. A single begin() call
    // per file populates the field table, label flags, and RGB flag (PLY yields an empty header).
    struct FileSchema
    {
      ray::LasHeader hdr;
      std::vector<SensorAttr> attrs;  ///< sensor (non-own) extras in VLR order
    };
    const int nfiles = static_cast<int>(cloud_files.files().size());
    std::vector<FileSchema> schemas(nfiles);
    std::vector<bool> has_labels(nfiles, false);
    bool union_has_labels = false;
    bool union_has_rgb = false;
    for (int f = 0; f < nfiles; ++f)
    {
      const std::string &fn = cloud_files.files()[f].name();
      ray::CloudReader reader;
      reader.begin(fn);
      schemas[f].hdr = reader.header();
      for (const auto &ef : schemas[f].hdr.extras)
      {
        if (ef.is_own) continue;
        SensorAttr attr;
        attr.name = ef.name; attr.size = ef.size; attr.offset = ef.sensor_offset;
        std::memcpy(attr.record.data(), ef.vlr_record, 192);
        schemas[f].attrs.push_back(attr);
      }
      if (schemas[f].hdr.has("tree_id") || schemas[f].hdr.has("stem_id"))
        { has_labels[f] = true; union_has_labels = true; }
      if (schemas[f].hdr.has_rgb)
        union_has_rgb = true;
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
      const uint16_t fp = static_cast<uint16_t>(10 + schemas[f].hdr.sensorExtraSize());
      needs_remap[f] = (fp != writer_pass_stride) || (remaps[f].size() < union_attrs.size());
      if (!needs_remap[f])
        for (const auto &r : remaps[f])
          if (r.src_off != r.dst_off) { needs_remap[f] = true; break; }
    }

    // Voxel deduplication ('all' mode only). When --dedup is set, points sharing a voxel cell are
    // collapsed to the first-seen point, streaming across all input files. The set persists across
    // every chunk and every file below, so it is the cross-file streaming dedup state.
    if (tiebreak_option.isSet() && !dedup_option.isSet())
    {
      std::cerr << "--tiebreak requires --dedup" << std::endl;
      usage();
    }
    std::set<Eigen::Vector3i, ray::Vector3iLess> dedup_vox_set;

    // Resolve --tiebreak tokens against the built-in fields and the union sensor-attr schema.
    // The resolved spec drives the best-wins post-step (deduplicateVoxel) after the file is written.
    std::vector<ray::ResolvedTiebreaker> resolved;
    if (tiebreak_option.isSet())
    {
      const std::vector<TiebreakToken> tokens = parseTiebreakTokens(tiebreak_spec.text());
      for (const auto &tok : tokens)
      {
        // Built-in computed fields.
        if (tok.name == "reflectance")      { resolved.push_back({ ray::TiebreakKind::Reflectance,      tok.ascending, 0, 0, 0 }); continue; }
        if (tok.name == "range")            { resolved.push_back({ ray::TiebreakKind::Range,            tok.ascending, 0, 0, 0 }); continue; }
        if (tok.name == "time")             { resolved.push_back({ ray::TiebreakKind::Time,             tok.ascending, 0, 0, 0 }); continue; }
        // Fixed LAS passthrough fields (bytes 0-9, present in every LAS/LAZ input).
        if (tok.name == "return_number")    { resolved.push_back({ ray::TiebreakKind::ReturnNumber,     tok.ascending, 0, 0, 0 }); continue; }
        if (tok.name == "number_of_returns"){ resolved.push_back({ ray::TiebreakKind::NumberOfReturns,  tok.ascending, 0, 0, 0 }); continue; }
        if (tok.name == "classification")   { resolved.push_back({ ray::TiebreakKind::Classification,   tok.ascending, 0, 0, 0 }); continue; }
        if (tok.name == "user_data")        { resolved.push_back({ ray::TiebreakKind::UserData,         tok.ascending, 0, 0, 0 }); continue; }
        if (tok.name == "scan_angle")       { resolved.push_back({ ray::TiebreakKind::ScanAngle,        tok.ascending, 0, 0, 0 }); continue; }
        if (tok.name == "point_source_id")  { resolved.push_back({ ray::TiebreakKind::PointSourceId,    tok.ascending, 0, 0, 0 }); continue; }
        // Named sensor extra-byte fields. Any input that has the field contributes its value;
        // inputs that lack it have those bytes zeroed (sentinel = 0).
        bool found = false;
        for (const auto &ua : union_attrs)
          if (ua.name == tok.name)
          {
            resolved.push_back({ ray::TiebreakKind::ExtraByte, tok.ascending,
                                 ua.offset, static_cast<uint8_t>(ua.size), ua.record[2] });
            found = true;
            break;
          }
        if (!found)
        {
          std::cerr << "Unknown tiebreak field: '" << tok.name << "' (field names are case-sensitive)" << std::endl;
          std::cerr << "  built-in: reflectance, range, time" << std::endl;
          std::cerr << "  fixed LAS: return_number, number_of_returns, classification, user_data, scan_angle, point_source_id" << std::endl;
          if (!union_attrs.empty())
          {
            std::cerr << "  sensor extra-byte fields found in input files:";
            for (const auto &ua : union_attrs)
              std::cerr << " " << ua.name;
            std::cerr << std::endl;
          }
          usage();
        }
      }
    }

    ray::CloudWriter writer;
    if (!writer.begin(combined_file, union_vlr, /*with_beam_id=*/false,
                      /*with_tree_id=*/union_has_labels, /*with_stem_id=*/union_has_labels,
                      /*with_rgb=*/union_has_rgb))
      usage();

    for (int i = 0; i < nfiles; ++i)
    {
      const std::string &fname = cloud_files.files()[i].name();
      const uint16_t file_pass_stride = static_cast<uint16_t>(10 + schemas[i].hdr.sensorExtraSize());
      std::vector<uint8_t> passthrough_buf;
      size_t passthrough_cursor = 0;  ///< byte offset into passthrough_buf for the next chunk
      // Label buffers: readLas appends tree_id/stem_id for files that declare them. Sliced
      // per chunk through label_cursor, mirroring the passthrough cursor above.
      std::vector<int32_t> tree_ids_buf, stem_ids_buf;
      size_t label_cursor = 0;  ///< point offset into the label buffers for the next chunk

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
        // Missing attrs keep their 0xFF sentinel; decodeExtraByte detects all-0xFF bytes and
        // returns NaN, which beats() treats as a universal loser in tiebreak comparisons.
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

        // Carry tree_id/stem_id as first-class label columns when the union output has labels.
        std::vector<int32_t> chunk_tree_ids, chunk_stem_ids;
        if (union_has_labels)
        {
          if (has_labels[i])
          {
            // Forward this file's labels, sliced at the label cursor.
            if (tree_ids_buf.size() >= label_cursor + n_pts)
              chunk_tree_ids.assign(tree_ids_buf.begin() + static_cast<ptrdiff_t>(label_cursor),
                                    tree_ids_buf.begin() + static_cast<ptrdiff_t>(label_cursor + n_pts));
            if (stem_ids_buf.size() >= label_cursor + n_pts)
              chunk_stem_ids.assign(stem_ids_buf.begin() + static_cast<ptrdiff_t>(label_cursor),
                                    stem_ids_buf.begin() + static_cast<ptrdiff_t>(label_cursor + n_pts));
          }
          else
          {
            // File has no labels: fill sentinels (-1 = unassigned, per raycloud.h).
            chunk_tree_ids.assign(n_pts, -1);
            chunk_stem_ids.assign(n_pts, -1);
          }
        }
        label_cursor += n_pts;

        // Streaming first-wins voxel dedup: keep only points whose voxel cell has not yet been
        // seen across all chunks/files. Compact every parallel array (geometry, passthrough,
        // labels) by the kept indices so they stay aligned. Skipped when --tiebreak is active:
        // the best-wins post-step needs every candidate point, so it does the voxel collapse itself.
        if (dedup_option.isSet() && !tiebreak_option.isSet())
        {
          const double w = dedup_width.value();
          const bool have_pass   = !chunk_pass.empty();
          const bool have_labels = !chunk_tree_ids.empty();
          size_t keep = 0;
          for (size_t j = 0; j < n_pts; ++j)
          {
            const Eigen::Vector3i key(static_cast<int>(std::floor(ends[j][0] / w)),
                                      static_cast<int>(std::floor(ends[j][1] / w)),
                                      static_cast<int>(std::floor(ends[j][2] / w)));
            if (!dedup_vox_set.insert(key).second)
              continue;
            if (keep != j)
            {
              starts[keep]  = starts[j];
              ends[keep]    = ends[j];
              times[keep]   = times[j];
              colours[keep] = colours[j];
              if (have_pass)
                std::memcpy(chunk_pass.data() + keep * writer_pass_stride,
                            chunk_pass.data() + j * writer_pass_stride, writer_pass_stride);
              if (have_labels)
              {
                chunk_tree_ids[keep] = chunk_tree_ids[j];
                chunk_stem_ids[keep] = chunk_stem_ids[j];
              }
            }
            ++keep;
          }
          starts.resize(keep);
          ends.resize(keep);
          times.resize(keep);
          colours.resize(keep);
          if (have_pass)
            chunk_pass.resize(keep * writer_pass_stride);
          if (have_labels)
          {
            chunk_tree_ids.resize(keep);
            chunk_stem_ids.resize(keep);
          }
        }

        writer.writeChunk(starts, ends, times, colours, chunk_pass, {}, chunk_tree_ids, chunk_stem_ids);
      };

      ray::CloudReader reader;
      if (!reader.begin(fname))
        usage();
      size_t num_bounded;
      // Only request label output for files that declare labels; others stay empty and are
      // sentinel-filled in concatenate when the union output carries labels.
      std::vector<int32_t> *tree_ids_out = has_labels[i] ? &tree_ids_buf : nullptr;
      std::vector<int32_t> *stem_ids_out = has_labels[i] ? &stem_ids_buf : nullptr;
      if (!reader.read(concatenate, num_bounded, 1.0, nullptr,
                       ray::computeReadChunkSize(), tree_ids_out, &passthrough_buf, stem_ids_out))
        usage();
    }
    writer.end();

    // Best-wins post-step: re-resolve each voxel to its highest-priority point per the --tiebreak
    // spec, rewriting the combined file in place. The streaming first-wins pass above already
    // collapsed each voxel to one point; this replaces that representative with the spec winner.
    if (dedup_option.isSet() && !resolved.empty())
    {
      if (!ray::deduplicateVoxel(combined_file, dedup_width.value(), resolved))
        usage();
    }
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
