// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This program provides the 'rayvoxel' command-line interface.
// It generates an AMAPVox-style Beer-Lambert voxel grid from a
// raycloudtools ray cloud.

#include "raylib/rayparse.h"
#include "raylib/rayutils.h"
#include "raylib/rayvoxel/raylasvoxelise.h"
#include "raylib/rayvoxel/raylasvoxelconfig.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace ray;

void usage()
{
  std::cout << "Usage: rayvoxel <cloud_file> [options]" << std::endl;
  std::cout << "Generates an advanced voxel grid from a ray cloud, calculating density and occlusion metrics." << std::endl;
  std::cout << "Input must be a raycloudtools .las/.laz ray cloud (with sx,sy,sz ray-start extra bytes)." << std::endl;
  std::cout << "The distance_to_sensor used for beam metrics is computed as the ray length (end - start).norm()." << std::endl << std::endl;
  std::cout << "Required Arguments:" << std::endl;
  std::cout << "  <cloud_file>                    Input ray cloud (.las, .laz)." << std::endl << std::endl;
  std::cout << "Processing Strategy:" << std::endl;
  std::cout << "  --parallel, -p                  Enable parallel in-memory processing. (Default: on)" << std::endl;
  std::cout << "  --no_parallel                   Disable parallel processing and run single-threaded." << std::endl;
  std::cout << "  --threads, -t <N>               Specify number of worker threads. (Default: all available cores)" << std::endl;
  std::cout << "  --out_of_core                   Enable out-of-core processing for datasets larger than RAM." << std::endl;
  std::cout << "  --ram_budget_mb <MB>            Approx. RAM limit per thread for out-of-core mode. Default: 1024." << std::endl;
  std::cout << "  --reserve_size <num>            Manually specify number of voxels to reserve for in-memory modes." << std::endl << std::endl;
  std::cout << "General Options:" << std::endl;
  std::cout << "  --grid_bounds_min <x,y,z>       Minimum bounds (corner) of the voxel grid. (Default: auto-detect)" << std::endl;
  std::cout << "  --grid_bounds_max <x,y,z>       Maximum bounds (corner) of the voxel grid. (Default: auto-detect)" << std::endl;
  std::cout << "  --voxel_size <val>              Voxel size (in metres). Default: 0.1." << std::endl;
  std::cout << "  --output_format <format>        Primary output format (text, netcdf, amapvox). Default: amapvox." << std::endl;
  std::cout << "  --weighting_method <method>     Method for weighting ray contributions ('equal', 'full'). Default: equal." << std::endl;
  std::cout << "  --occlusion, -o                 Enable occlusion mapping (traces rays beyond last hits)." << std::endl << std::endl;
  std::cout << "Output Content Control:" << std::endl;
  std::cout << "  --write_empty, -e               Primary output includes all observed voxels (empty and filled)." << std::endl;
  std::cout << "  --write_filled, -w              Creates an *additional* output file containing only FILLED voxels." << std::endl;
  std::cout << "  --write_amapvox_also, -a        In addition to primary output, also write an AMAPVox .vox file." << std::endl << std::endl;
  std::cout << "DTM / Ground Control Options:" << std::endl;
  std::cout << "  --dtm <file.ply>                Use an external PLY mesh as a DTM for ground clipping and metrics." << std::endl;
  std::cout << "  --dtm_from_class <class>        Generate a DTM internally from points with the given classification code." << std::endl;
  std::cout << "  --dtm_cell_size <val>           Cell size (in metres) for the DTM rasterization. Default: 1.0." << std::endl << std::endl;
  std::cout << "Post-Processing Options:" << std::endl;
  std::cout << "  --flat_top_compensation         Apply density correction for flat-topped canopies (e.g., crops)." << std::endl;
  std::cout << "  --neighbour_priors <min_rays>   Apply spatial smoothing to voxels with fewer than <min_rays>. Default: 0 (off)." << std::endl << std::endl;
  std::cout << "Vegetation & Beam Metrics (Tier 1 & 2):" << std::endl;
  std::cout << "  --veg_metrics                   Enable calculation of Tier 1 (angle, dist) and advanced vegetation metrics." << std::endl;
  std::cout << "  --leaf_classes [<field>:]c1,c2  Classification codes for leaves (PAD_leaf). Optional <field>: selects a named LAS extra-byte field." << std::endl;
  std::cout << "  --wood_classes [<field>:]c1,c2  Classification codes for wood/stems (PAD_wood). Optional <field>: selects a named LAS extra-byte field." << std::endl;
  std::cout << "  --lad <type>                    Leaf Angle Distribution for G-function (e.g., spherical, ellipsoidal). Default: spherical." << std::endl;
  std::cout << "  --lad_params <p1,p2>            Comma-separated parameters for the LAD (e.g., chi for ellipsoidal)." << std::endl;
  std::cout << "  --beam_metrics                  Enable Tier 2 beam-based metrics (transmittance, etc.)." << std::endl;
  std::cout << "  --laser_spec <name>             Select a predefined laser specification (e.g., VZ-400)." << std::endl;
  std::cout << "  --beam_params <diam,div>        Manually specify beam diameter (m) and divergence (rad)." << std::endl;
  std::cout << "  --subvoxel_split <N>            Enable exploration rate calculation with an N x N x N grid (N=2,3,4). Default: 0 (off)." << std::endl;
  std::cout << "  --no_inclination_dist           Disable the inclination-distribution pass (skips KNN normal estimation) while keeping --veg_metrics." << std::endl;
  std::cout << "  --output_iad                    Write per-bin LIAD/WIAD/PIAD histogram columns to output (default: off; G scalars are always written)." << std::endl;
  std::cout << "  --n_iad_bins <N>                Number of inclination-angle histogram bins over [0, pi/2]. Default: 18." << std::endl;
  std::cout << "  --attenuation_method <methods>  Comma-separated PAD/LAD/WAD estimators: fpl (default), ppl, transmittance, bailey." << std::endl;
  std::cout << "  --knn_normal <N>                Number of nearest neighbours used for per-point normal estimation. Default: 10." << std::endl;
  std::cout << "  --triangle_lmax <m>             Max triangle edge length for Bailey facets (only used with --attenuation_method bailey). Default: 0.05." << std::endl;
  exit(1);
}

int main_function(int argc, char *argv[])
{
  // --- Argument Definitions ---
  FileArgument cloud_file;

  // Processing Strategy
  OptionalFlagArgument parallel_flag("parallel", 'p');
  OptionalFlagArgument no_parallel_flag("no_parallel", '\0');
  IntArgument num_threads_val(1, 256, 0);
  OptionalKeyValueArgument num_threads("threads", 't', &num_threads_val);
  OptionalFlagArgument out_of_core_flag("out_of_core", '\0');
  IntArgument ram_budget_mb_val(64, 65536, 1024);
  OptionalKeyValueArgument ram_budget_mb("ram_budget_mb", '\0', &ram_budget_mb_val);
  IntArgument reserve_size_val(0, 2000000000, 0);
  OptionalKeyValueArgument reserve_size("reserve_size", '\0', &reserve_size_val);

  // General Options
  DoubleArgument voxel_size_val(0.001, 1000.0, 0.1);
  OptionalKeyValueArgument voxel_size("voxel_size", 's', &voxel_size_val);
  Vector3dArgument grid_bounds_min_val;
  OptionalKeyValueArgument grid_bounds_min("grid_bounds_min", '\0', &grid_bounds_min_val);
  Vector3dArgument grid_bounds_max_val;
  OptionalKeyValueArgument grid_bounds_max("grid_bounds_max", '\0', &grid_bounds_max_val);
  StringArgument output_format_val("amapvox");
  OptionalKeyValueArgument output_format("output_format", 'f', &output_format_val);
  StringArgument weighting_method_val("equal");
  OptionalKeyValueArgument weighting_method("weighting_method", '\0', &weighting_method_val);
  OptionalFlagArgument occlusion("occlusion", 'o');

  // Output Content
  OptionalFlagArgument write_amapvox_also("write_amapvox_also", 'a');
  OptionalFlagArgument write_empty("write_empty", 'e');
  OptionalFlagArgument write_filled("write_filled", 'w');

  // DTM Options
  FileArgument dtm_file_val;
  OptionalKeyValueArgument dtm_file("dtm", '\0', &dtm_file_val);
  IntArgument dtm_from_class_val(0, 255, -1);
  OptionalKeyValueArgument dtm_from_class("dtm_from_class", '\0', &dtm_from_class_val);
  DoubleArgument dtm_cell_size_val(0.1, 100.0, 1.0);
  OptionalKeyValueArgument dtm_cell_size("dtm_cell_size", '\0', &dtm_cell_size_val);

  // Post-Processing
  OptionalFlagArgument flat_top_compensation("flat_top_compensation", '\0');
  IntArgument neighbour_priors_val(0, 1000, 0);
  OptionalKeyValueArgument neighbour_priors("neighbour_priors", 'n', &neighbour_priors_val);

  // Metrics
  OptionalFlagArgument veg_metrics("veg_metrics", '\0');
  StringArgument leaf_classes_val("");
  OptionalKeyValueArgument leaf_classes("leaf_classes", '\0', &leaf_classes_val);
  StringArgument wood_classes_val("");
  OptionalKeyValueArgument wood_classes("wood_classes", '\0', &wood_classes_val);
  StringArgument lad_val("spherical");
  OptionalKeyValueArgument lad("lad", '\0', &lad_val);
  StringArgument lad_params_val("");
  OptionalKeyValueArgument lad_params("lad_params", '\0', &lad_params_val);
  OptionalFlagArgument beam_metrics("beam_metrics", '\0');
  StringArgument laser_spec_val("");
  OptionalKeyValueArgument laser_spec("laser_spec", '\0', &laser_spec_val);
  Vector2dArgument beam_params_val;
  OptionalKeyValueArgument beam_params("beam_params", '\0', &beam_params_val);
  IntArgument subvoxel_split_val(0, 4, 0);
  OptionalKeyValueArgument subvoxel_split("subvoxel_split", '\0', &subvoxel_split_val);
  OptionalFlagArgument inclination_dist("inclination_dist", '\0');
  OptionalFlagArgument no_inclination_dist("no_inclination_dist", '\0');
  OptionalFlagArgument output_iad("output_iad", '\0');
  IntArgument n_iad_bins_val(1, 180, 18);
  OptionalKeyValueArgument n_iad_bins("n_iad_bins", '\0', &n_iad_bins_val);
  StringArgument attenuation_method_val("fpl");
  OptionalKeyValueArgument attenuation_method("attenuation_method", '\0', &attenuation_method_val);
  IntArgument knn_normal_val(2, 1000, 10);
  OptionalKeyValueArgument knn_normal("knn_normal", '\0', &knn_normal_val);
  DoubleArgument triangle_lmax_val(0.001, 10.0, 0.05);
  OptionalKeyValueArgument triangle_lmax("triangle_lmax", '\0', &triangle_lmax_val);

  // --- Parse Command Line ---
  std::vector<FixedArgument *> fixed_args = { &cloud_file };
  std::vector<OptionalArgument *> optional_args = {
      &parallel_flag, &no_parallel_flag, &num_threads, &out_of_core_flag, &ram_budget_mb, &reserve_size,
      &voxel_size, &grid_bounds_min, &grid_bounds_max, &output_format, &weighting_method,
      &occlusion, &write_amapvox_also, &write_empty, &write_filled,
      &dtm_file, &dtm_from_class, &dtm_cell_size,
      &flat_top_compensation, &neighbour_priors,
      &veg_metrics, &leaf_classes, &wood_classes, &lad, &lad_params,
      &beam_metrics, &laser_spec, &beam_params, &subvoxel_split,
      &inclination_dist, &no_inclination_dist, &output_iad, &n_iad_bins, &attenuation_method, &knn_normal,
      &triangle_lmax };

  if (!parseCommandLine(argc, argv, fixed_args, optional_args)) {
    usage();
  }

  // --- Validate Arguments ---
  if (cloud_file.name().empty()) {
      std::cerr << "Error: An input cloud file must be specified." << std::endl; usage();
  }
  if (grid_bounds_min.isSet() != grid_bounds_max.isSet()) {
      std::cerr << "Error: You must specify both --grid_bounds_min and --grid_bounds_max, or neither." << std::endl; usage();
  }
  if (laser_spec.isSet() && beam_params.isSet()) {
      std::cerr << "Error: --laser_spec and --beam_params are mutually exclusive." << std::endl; usage();
  }
  if (beam_metrics.isSet() && !laser_spec.isSet() && !beam_params.isSet()) {
      std::cerr << "Error: --beam_metrics requires either --laser_spec or --beam_params." << std::endl; usage();
  }
  if ((parallel_flag.isSet() || num_threads.isSet()) && no_parallel_flag.isSet()) {
      std::cerr << "Error: --no_parallel cannot be used with --parallel or --threads." << std::endl; usage();
  }
  if (out_of_core_flag.isSet() && reserve_size.isSet()) {
      std::cerr << "Warning: --reserve_size is ignored when using --out_of_core mode." << std::endl;
  }
  if (dtm_file.isSet() && dtm_from_class.isSet()) {
      std::cerr << "Error: --dtm and --dtm_from_class are mutually exclusive. Please specify only one." << std::endl; usage();
  }
  std::vector<std::string> parsed_methods;
  if (attenuation_method.isSet()) {
      std::stringstream ss(attenuation_method_val.text());
      std::string token;
      while (std::getline(ss, token, ',')) {
          token.erase(0, token.find_first_not_of(" \t"));
          token.erase(token.find_last_not_of(" \t") + 1);
          for (char& c : token) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
          if (!token.empty()) parsed_methods.push_back(token);
      }
      for (const auto& m : parsed_methods) {
          if (m != "fpl" && m != "ppl" && m != "transmittance" && m != "bailey") {
              std::cerr << "Error: --attenuation_method: unknown method '" << m << "'. Must be fpl, ppl, transmittance, or bailey." << std::endl;
              return 1;
          }
      }
  }
  const bool any_bailey = std::any_of(parsed_methods.begin(), parsed_methods.end(),
                                       [](const std::string& m){ return m == "bailey"; });
  if (any_bailey && (leaf_classes_val.text().empty() || wood_classes_val.text().empty())) {
      std::cerr << "Error: --attenuation_method bailey requires both --leaf_classes and --wood_classes." << std::endl;
      return 1;
  }

  // --- Populate Parameters Struct ---
  // All command-line arguments are now consolidated into a single
  // VoxelizationParameters struct. This simplifies the call to the core
  // library function and makes the code much cleaner.
  VoxelizationParameters params;
  params.cloud_name = cloud_file.name();
  params.voxel_size = voxel_size_val.value();
  if (grid_bounds_min.isSet()) {
    params.grid_bounds_min = grid_bounds_min_val.value();
    params.grid_bounds_max = grid_bounds_max_val.value();
  }
  params.use_occlusion = occlusion.isSet();
  params.write_empty_voxels = write_empty.isSet();
  params.output_format = output_format_val.text();
  params.weighting_method = weighting_method_val.text();
  params.write_amapvox_also = write_amapvox_also.isSet();
  params.apply_flat_top = flat_top_compensation.isSet();
  params.neighbour_prior_min_rays = neighbour_priors_val.value();
  params.write_filled = write_filled.isSet();
  params.calc_veg_metrics = veg_metrics.isSet();
  // --inclination_dist implies --veg_metrics: if passed standalone, silently enable veg_metrics.
  if (inclination_dist.isSet()) params.calc_veg_metrics = true;
  params.leaf_classes_str = leaf_classes_val.text();
  params.wood_classes_str = wood_classes_val.text();
  // Derive has_leaf/has_wood by stripping any "field:" prefix before checking for content.
  {
    auto strip_prefix = [](const std::string& s) {
      auto c = s.find(':'); return (c != std::string::npos) ? s.substr(c + 1) : s;
    };
    params.has_leaf = !strip_prefix(leaf_classes_val.text()).empty();
    params.has_wood = !strip_prefix(wood_classes_val.text()).empty();
  }
  if (!params.has_leaf && !params.has_wood && (veg_metrics.isSet() || inclination_dist.isSet()))
    std::cerr << "Info: no --leaf_classes or --wood_classes specified; only pad_* (plant) columns will be written.\n";
  else if (!params.has_leaf && params.has_wood)
    std::cerr << "Info: no --leaf_classes specified; lad_* columns will be omitted.\n";
  else if (params.has_leaf && !params.has_wood)
    std::cerr << "Info: no --wood_classes specified; wad_* columns will be omitted.\n";
  params.lad = lad_val.text();
  params.lad_params_str = lad_params_val.text();
  params.calc_beam_metrics = beam_metrics.isSet();
  params.laser_spec_name = laser_spec_val.text();
  if (beam_params.isSet()) {
    params.beam_params = beam_params_val.value();
  }
  params.subvoxel_split = subvoxel_split_val.value();
  // IAD is on by default whenever vegetation metrics are active; --no_inclination_dist opts out.
  params.calc_inclination_dist = (veg_metrics.isSet() || inclination_dist.isSet()) && !no_inclination_dist.isSet();
  params.output_iad = output_iad.isSet();
  params.n_iad_bins = n_iad_bins_val.value();
  if (!parsed_methods.empty()) params.attenuation_methods = parsed_methods;
  params.knn_normal = knn_normal_val.value();
  params.reserve_size = static_cast<size_t>(reserve_size_val.value());
  params.triangle_lmax     = triangle_lmax_val.value();
  if (any_bailey) params.calc_inclination_dist = true;  // ensures KNN matrix exists

  // DTM parameters
  // The dtm_cell_size now applies to both DTM creation methods.
  if (dtm_file.isSet()) {
    params.dtm_file = dtm_file_val.name();
  }
  if (dtm_from_class.isSet()) {
    params.dtm_from_class = dtm_from_class_val.value();
  }
  params.dtm_cell_size = dtm_cell_size_val.value();

  // Determine number of threads to use. 0 means auto-detect.
  // Strategy parameters
  params.use_ooc = out_of_core_flag.isSet();
  params.ram_budget_mb = static_cast<size_t>(ram_budget_mb_val.value());
  if (no_parallel_flag.isSet()) {
    params.num_threads = 1;
  } else if (num_threads.isSet()) {
    params.num_threads = num_threads_val.value();
  } else {
    params.num_threads = 0; // 0 means auto-detect
  }

  // --- Execute Voxelization ---
  std::cout << "Starting advanced voxelization process..." << std::endl;
  bool success = generateVoxelGrid(params);

  if (success) {
    std::cout << "Voxelization process completed successfully." << std::endl;
    return 0;
  } else {
    std::cerr << "Voxelization process failed." << std::endl;
    return 1;
  }
}

// The main function is a wrapper for runWithMemoryCheck to catch out-of-memory exceptions
int main(int argc, char* argv[])
{
  return runWithMemoryCheck(main_function, argc, argv);
}
