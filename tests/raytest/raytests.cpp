// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include "raycloud.h"
#include "raymesh.h"
#include "rayply.h"
#include "rayforeststructure.h"
#include "raylasdecode.h"
#include <vector>
#include <fstream>
#include <string>
#include <gtest/gtest.h>
#include <cstdlib>
#include "raylaz.h"
#include "raysysinfo.h"
#include "rayvoxel/raylasvoxelprocessor.h"
#include "raycuboid.h"
#include "rayvoxel/raylasbailey.h"

/// Raycloud testing framework. In each test, the statistics of the resulting clouds are compared to the statistics
/// of the cloud when it was confirmed to be operating correctly. 
namespace raytest
{
  /// Issues the specified system command, including the required prefix on non-windows systems.
  int command(const std::string &system_command)
  {
    #ifdef _WIN32
    return system(system_command);
    #else
    return system(("./" + system_command).c_str());
    #endif // _WIN32
  }

  /// Issues the command to copy a file, which is a platform dependent system command.
  int copy(const std::string &copy_command)
  {
    #ifdef _WIN32
    return system("copy " + copy_command);
    #else
    return system(("cp " + copy_command).c_str());
    #endif // _WIN32
  }

  /// Compare the statistical (1st and 2nd order) moments of the two ray clouds. This almost surely
  /// detects differing clouds, and always equal clouds, given a tolerance @c eps.
  void compareMoments(const Eigen::ArrayXd &m1, const std::vector<double> &m2, double eps = 0.1)
  {
    for (size_t i = 0; i<m2.size(); i++)
    {
      EXPECT_GT(m1[i], m2[i]-eps);
      EXPECT_LT(m1[i], m2[i]+eps);
    }
  }

  /// Compare the statistical (1st and 2nd order) moments of the two ray clouds. This almost surely
  /// detects differing clouds, and always equal clouds, given a tolerance @c eps.
  void compareMomentsPercentageError(const Eigen::ArrayXd &m1, const std::vector<double> &m2, double percentage = 5.0)
  {
    double eps = 0.01 * percentage;
    for (size_t i = 0; i<m2.size(); i++)
    {
      EXPECT_GE(m1[i], m2[i] * (1.0 - eps));
      EXPECT_LE(m1[i], m2[i] * (1.0 + eps));
    }
  }

  /// Creates two copies of the same room with a rotational difference, then aligns the first onto the second 
  TEST(Basic, RayAlign)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(copy("room.las room2.las"), 0);
    EXPECT_EQ(command("rayrotate room2.las 0,0,35"), 0);
    EXPECT_EQ(command("rayalign room.las room2.las"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_aligned.las"));
    compareMoments(cloud.getMoments(), {-0.0618268, -0.077552, 0.0531072, 7.58334e-08, 7.97642e-08, 1.93877e-08, -0.180532, -0.219257, 0.0654452, 2.47241, 2.08183, 1.28226, 17.539, 10.1994, 0.304682, 0.761892, 0.429502, 0.987362, 0.318932, 0.225742, 0.389901, 0.111705});  }

  /// Colours a room according to the normal direction of the surfaces, comparing to the expected results
  TEST(Basic, RayColour)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(command("raycolour room.las normal"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_coloured.las"));
    compareMoments(cloud.getMoments(), {-0.108066, -0.0410134, 0.052168, 7.05134e-08, 8.45038e-08, 1.93877e-08, -0.276144, -0.0760758, 0.065631, 2.42455, 2.13738, 1.28226, 17.539, 10.1994, 0.304682, 0.761892, 0.429502, 0.987362, 0.248361, 0.203648, 0.385192, 0.111705});
  }
  
  /// Creates two rooms, with different transformations, then combines them, and compares to the expected result.
  TEST(Basic, RayCombine)
  {
    EXPECT_EQ(command("./raycreate room 1"), 0);
    EXPECT_EQ(copy("room.las room2.las"), 0);
    EXPECT_EQ(command("./raytranslate room2.las 0,0,1"), 0);
    EXPECT_EQ(command("./rayrotate room2.las 0,0,35"), 0);
    EXPECT_EQ(command("./raycombine min room.las room2.las 1 rays"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_combined.las"));
    compareMoments(cloud.getMoments(), {-0.0867714, -0.0679941, 0.546619, 0.0215326, 0.0272819, 0.499969, -0.305657, -0.186353, 0.582642, 2.95777, 2.47531, 1.63323, 17.4967, 10.1789, 0.305355, 0.763356, 0.427376, 0.979005, 0.318409, 0.225661, 0.389366, 0.143369});
  }
  
  /// Creates a building with random seed 1, and compares to the expected results
  TEST(Basic, RayCreate)
  {
    EXPECT_EQ(command("raycreate building 1"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("building.las"));
    compareMoments(cloud.getMoments(), {-3.168, 16.5472, 7.04175, 3.28539, 18.0871, 2.96654, -3.19551, 16.564, 7.30826, 4.14359, 18.2463, 3.34431, 935.715, 540.236, 0.499997, 0.500408, 0.429416, 0.998894, 0.372134, 0.372389, 0.390499, 0.0332398});
  }
  
  /// Creates a forest and decimates it, comparing to the expected result 
  TEST(Basic, RayDecimate)
  {
    EXPECT_EQ(command("raycreate forest 1"), 0);
    EXPECT_EQ(command("raydecimate forest.las 10 cm"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("forest_decimated.las"));
    // Below does not compare the time values (or the colour values, which are based on time here)
    // because spatial decimation does not constraint which time it picks points from.
    compareMoments(cloud.getMoments(), {-0.222571, 1.08156, 1.67264, 6.00755, 5.78731, 0.508713, -0.202668, 1.09517, 2.6238, 6.0285, 5.85715, 3.22093, 69.0574, 35.2775, 0.48969, 0.498403, 0.443549, 1, 0.379062, 0.366963, 0.389535, 0});
  }

  /// Creates a room, and calls denoise using a fixed distance threshols, and compares to expected result
  TEST(Basic, RayDenoise)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(command("raydenoise room.las 3 cm"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_denoised.las"));
    compareMoments(cloud.getMoments(), {-0.108066, -0.0410134, 0.052168, 8.67026e-08, 8.81787e-08, 2.24394e-08, -0.464107, -0.113806, 0.161496, 2.82122, 2.34281, 1.35279, 17.81, 10.2005, 0.297047, 0.758802, 0.440232, 0.975166, 0.317215, 0.226682, 0.390971, 0.155618});
  }

  /// Creates two rooms, the second is decimated and transformed, then rayrestore is called to apply this transformation to
  /// the first (high resolution) room
  TEST(Basic, RayRestore)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(copy("room.las room2.las"), 0);
    EXPECT_EQ(command("raydecimate room2.las 10 cm"), 0);
    EXPECT_EQ(command("raytranslate room2_decimated.las 1,2,3"), 0);
    EXPECT_EQ(command("rayrotate room2_decimated.las 0,0,-50"), 0);
    EXPECT_EQ(command("rayrestore room2_decimated.las 10 cm room.las"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_restored.las"));
    compareMoments(cloud.getMoments(), {2.07399, 0.575952, 3.05217, 7.85442e-08, 7.70963e-08, 1.93877e-08, 1.9391, 0.682169, 3.06563, 2.10068, 2.45642, 1.28226, 17.539, 10.1994, 0.304682, 0.761892, 0.429502, 0.987362, 0.318932, 0.225742, 0.389901, 0.111705});
  }  

  /// Creates a forest and rotates it in all three axes, comparing to the expected result
  TEST(Basic, RayRotate)
  {
    EXPECT_EQ(command("raycreate forest 1"), 0);
    EXPECT_EQ(command("rayrotate forest.las 10,20,30"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("forest.las"));
    compareMoments(cloud.getMoments(), {-0.254879, 0.846076, 2.02322, 5.62648, 5.68622, 2.56306, 0.266873, 0.772016, 3.2888, 5.54595, 5.69021, 4.28394, 62.683, 36.1903, 0.514327, 0.504407, 0.413534, 1, 0.372377, 0.365965, 0.391709, 0});
  }  

  /// Creates a room and smooths this ray cloud, comparing to the expected result
  TEST(Basic, RaySmooth)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(command("raysmooth room.las"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_smooth.las"));
    compareMoments(cloud.getMoments(), {-0.108066, -0.0410134, 0.052168, 7.05134e-08, 8.45038e-08, 1.93877e-08, -0.27615, -0.0761079, 0.0656267, 2.42413, 2.13691, 1.28163, 17.539, 10.1994, 0.304682, 0.761892, 0.429502, 0.987362, 0.318932, 0.225742, 0.389901, 0.111705});
  }  

  /// Creates a room, then splits it around a plane, comparing agaisnt the expected result
  TEST(Basic, RaySplit)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(command("raysplit room.las plane 0,0.1,1.5"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_outside.las"));
    compareMoments(cloud.getMoments(), {-0.467731, 1.05075, 1.43662, 2.20441, 1.60162, 0.106775, -0.77974, 1.03139, 1.57353, 3.67521, 2.64766, 0.485084, 17.3995, 10.279, 0.311066, 0.759795, 0.425206, 0.951355, 0.321609, 0.226785, 0.39073, 0.215125});
  }  

  /// Creates a room and runs raytransients, comparing the identified transients ray cloud to the expected results
  TEST(Basic, RayTransients)
  {
    EXPECT_EQ(command("raycreate room 2"), 0);
    EXPECT_EQ(command("raytransients min room.las 1 rays"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room_transient.las"));
    compareMoments(cloud.getMoments(), {-1.05406, -0.240721, -0.0629182, 5.05649e-08, 3.32941e-08, 2.54759e-08, 0.268724, -0.136746, -0.596782, 1.04798, 0.921776, 0.527205, 32.1452, 6.7491, 0.205871, 0.395641, 0.884296, 1, 0.225501, 0.296487, 0.153923, 0});
  }  

  /// Creates a forest and translates it in all three axes, comparing to the expected result
  TEST(Basic, RayTranslate)
  {
    EXPECT_EQ(command("raycreate forest 1"), 0);
    EXPECT_EQ(command("raytranslate forest.las 10,20,30"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("forest.las"));
    compareMoments(cloud.getMoments(), {9.66298, 21.3454, 31.7177, 6.0926, 5.75511, 0.56438, 9.69155, 21.3605, 33.0883, 6.10555, 5.82564, 3.20507, 62.683, 36.1903, 0.514327, 0.504407, 0.413534, 1, 0.372377, 0.365965, 0.391709, 0});
  }

#if RAYLIB_WITH_QHULL
  /// Creates a terrain ray cloud, then wraps it from below, comparing the mesh to the expected results
  TEST(Basic, RayWrap)
  {
    EXPECT_EQ(command("raycreate terrain 1"), 0);
    EXPECT_EQ(command("raywrap terrain.las upwards 1.0"), 0);
    ray::Mesh mesh;
    EXPECT_TRUE(ray::readPlyMesh("terrain_mesh.ply", mesh));
    compareMoments(mesh.getMoments(), {0.0386662, -1.52168, -0.139079, 3.30621, 3.35391, 0.705937});
  }  

  /// Tests extraction of terrain and extraction of trees
  TEST(Basic, RayExtract)
  {
    EXPECT_EQ(command("raycreate forest 2"), 0);
    EXPECT_EQ(command("rayextract terrain forest.las"), 0);
    ray::Mesh mesh;
    EXPECT_TRUE(ray::readPlyMesh("forest_mesh.ply", mesh));
    compareMoments(mesh.getMoments(), {-0.00147491, -0.00191917, -0.0617946, 5.77995, 5.80266, 0.0426993});
    EXPECT_EQ(command("rayextract trees forest.las forest_mesh.ply"), 0);

    ray::ForestStructure forest;
    EXPECT_TRUE(forest.load("forest_trees.txt"));
    // Values re-baselined for .las input (raycreate seed 2): LAS 1mm quantization
    // slightly shifts the KNN graph used in tree reconstruction; algorithm unchanged.
    compareMomentsPercentageError(forest.getMoments(), {21, 20.0515, 1062.27, 1.4679, 0.11764, 2.28854, 21918, 0, 98.7857});

    EXPECT_EQ(command("rayextract forest forest.las --ground forest_mesh.ply"), 0);

    ray::ForestStructure forest2;
    EXPECT_TRUE(forest2.load("forest_forest.txt"));
    compareMoments(forest2.getMoments(), {11, 8.43829, 586.427, 1.40054, 0.200644, 0, 16697, 8.49419, 3.09917});

    EXPECT_EQ(command("rayextract trunks forest.las"), 0);

    ray::ForestStructure forest3;
    EXPECT_TRUE(forest3.load("forest_trunks.txt"));
    compareMoments(forest3.getMoments(), {21, 20.0797, 1124.61, 1.60427, 0.135159, 0, 0, 0, 0});
  }

  /// Tests rayextract segment + reconstruct split pipeline (new in refactor).
  /// Uses forest.las since raycreate now writes .las by default.
  TEST(Basic, RayExtractRefactor)
  {
    // Create forest data (las format from current raycreate)
    EXPECT_EQ(command("raycreate forest 2"), 0);
    EXPECT_EQ(command("rayextract terrain forest.las"), 0);
    { ray::Mesh m; EXPECT_TRUE(ray::readPlyMesh("forest_mesh.ply", m)); }

    // --- Test 1: rayextract segment produces segmented cloud + seeds file ---
    EXPECT_EQ(command("rayextract segment forest.las --ground forest_mesh.ply"), 0);

    ray::Cloud seg_cloud;
    EXPECT_TRUE(seg_cloud.load("forest_segmented.las"));
    EXPECT_FALSE(seg_cloud.tree_ids.empty());
    // Confirm at least some points have valid tree_ids (non -1)
    int labeled = 0;
    for (auto tid : seg_cloud.tree_ids)
      if (tid != -1) labeled++;
    EXPECT_GT(labeled, 0);

    // Seeds file must exist and have sensible content
    {
      std::ifstream ifs("forest_segmented_seeds.txt");
      EXPECT_TRUE(ifs.is_open());
    }

    // --- Test 2: rayextract reconstruct from segmented cloud ---
    EXPECT_EQ(command("rayextract reconstruct forest_segmented.las forest_mesh.ply"), 0);

    ray::ForestStructure recon_forest;
    EXPECT_TRUE(recon_forest.load("forest_segmented_trees.txt"));
    // Should reconstruct at least some trees
    EXPECT_GT(recon_forest.getMoments()[0], 0);

    // The new format must contain tree_id,stem_id header
    {
      std::ifstream ifs("forest_segmented_trees.txt");
      std::string line;
      bool found_id_header = false;
      while (std::getline(ifs, line))
      {
        if (line.find("tree_id") != std::string::npos &&
            line.find("stem_id") != std::string::npos)
        {
          found_id_header = true;
          break;
        }
      }
      EXPECT_TRUE(found_id_header);
    }

    // --- Test 3: legacy rayextract trees (no mask) produces OLD format ---
    EXPECT_EQ(command("rayextract trees forest.las forest_mesh.ply"), 0);
    {
      std::ifstream ifs("forest_trees.txt");
      EXPECT_TRUE(ifs.is_open());
      std::string line;
      bool has_old_header = false;
      while (std::getline(ifs, line))
      {
        if (line.find('#') == std::string::npos && !line.empty())
        {
          // Old format header has no tree_id/stem_id prefix
          if (line.find("tree_id") == std::string::npos)
            has_old_header = true;
          break;
        }
      }
      EXPECT_TRUE(has_old_header);
    }
  }
  /// Tests that rayextract reconstruct works without a seeds file by synthesizing seeds from
  /// the tree_id/stem_id labels already present in the segmented cloud.
  TEST(Basic, RayReconstructNoSeeds)
  {
    EXPECT_EQ(command("raycreate forest 2"), 0);
    EXPECT_EQ(command("rayextract terrain forest.las"), 0);
    EXPECT_EQ(command("rayextract segment forest.las --ground forest_mesh.ply"), 0);

    // Remove the auto-generated seeds file to force the synthesis path.
    EXPECT_EQ(remove("forest_segmented_seeds.txt"), 0);

    // Reconstruct must succeed by synthesizing seeds from cloud labels alone.
    EXPECT_EQ(command("rayextract reconstruct forest_segmented.las forest_mesh.ply"), 0);

    ray::ForestStructure recon;
    EXPECT_TRUE(recon.load("forest_segmented_trees.txt"));
    EXPECT_GT(recon.getMoments()[0], 0);

    std::ifstream mesh_ifs("forest_segmented_trees_mesh.ply");
    EXPECT_TRUE(mesh_ifs.is_open());
  }
#endif  // RAYLIB_WITH_QHULL

#if RAYLIB_WITH_LAS
  // Verify that the explicit `bound` extra-byte field takes precedence over alpha when both are
  // present. This exercises the fix for old files (bound_offset == -1) vs new files.
  TEST(BoundDecode, BoundAuthoritativeWhenPresent)
  {
    using namespace ray;

    // Extra-bytes layout: [sx(4) sy(4) sz(4) alpha(1) bound(1)] = 14 bytes total.
    // alpha_offset = 12, bound_offset = 13.
    DecodeContext ctx;
    ctx.is_raycloud = true;
    ctx.alpha_offset = 12;
    ctx.bound_offset = 13;
    ctx.sx_offset = 0;
    ctx.sy_offset = 4;
    ctx.sz_offset = 8;

    auto makePoint = [](uint8_t alpha, uint8_t bound, std::array<uint8_t, 14> &extra) -> laszip_point_struct {
      std::fill(extra.begin(), extra.end(), 0);
      extra[12] = alpha;
      extra[13] = bound;
      laszip_point_struct pt = {};
      pt.extra_bytes = extra.data();
      pt.num_extra_bytes = static_cast<laszip_I32>(extra.size());
      pt.gps_time = 0.0;
      return pt;
    };

    // Case 1: alpha > 0, bound == 0  →  bound wins: ray must be unbound (intensity == 0).
    {
      std::array<uint8_t, 14> extra;
      laszip_point_struct pt = makePoint(5, 0, extra);
      std::vector<Eigen::Vector3d> starts, ends;
      std::vector<double> times;
      std::vector<RGBA> colours;
      std::vector<uint8_t> intensities;
      size_t num_bounded = 0;
      decodePointRecord(&pt, ctx, Eigen::Vector3d::Zero(), starts, ends, times, colours, intensities,
                        num_bounded, nullptr, nullptr, nullptr, nullptr);
      EXPECT_EQ(intensities.back(), 0u) << "bound=0 must suppress stray alpha";
      EXPECT_EQ(num_bounded, 0u);
    }

    // Case 2: alpha == 0, bound == 1  →  bound wins: ray must be bounded (intensity == 1).
    {
      std::array<uint8_t, 14> extra;
      laszip_point_struct pt = makePoint(0, 1, extra);
      std::vector<Eigen::Vector3d> starts, ends;
      std::vector<double> times;
      std::vector<RGBA> colours;
      std::vector<uint8_t> intensities;
      size_t num_bounded = 0;
      decodePointRecord(&pt, ctx, Eigen::Vector3d::Zero(), starts, ends, times, colours, intensities,
                        num_bounded, nullptr, nullptr, nullptr, nullptr);
      EXPECT_GT(intensities.back(), 0u) << "bound=1 must mark ray as bounded when alpha is zero";
      EXPECT_EQ(num_bounded, 1u);
    }

    // Case 3: old file (bound_offset == -1) — alpha alone drives boundedness.
    {
      DecodeContext old_ctx = ctx;
      old_ctx.bound_offset = -1;

      std::array<uint8_t, 14> extra;
      laszip_point_struct pt = makePoint(7, 0, extra);  // bound byte says unbound, but field absent
      std::vector<Eigen::Vector3d> starts, ends;
      std::vector<double> times;
      std::vector<RGBA> colours;
      std::vector<uint8_t> intensities;
      size_t num_bounded = 0;
      decodePointRecord(&pt, old_ctx, Eigen::Vector3d::Zero(), starts, ends, times, colours,
                        intensities, num_bounded, nullptr, nullptr, nullptr, nullptr);
      EXPECT_EQ(intensities.back(), 7u) << "old file: alpha is sole source of truth";
      EXPECT_EQ(num_bounded, 1u);
    }
  }
#endif  // RAYLIB_WITH_LAS

#if RAYLIB_WITH_LAS
  TEST(RayLas, RoundTrip)
  {
    // raycreate writes .las directly on this branch — load it to exercise the write→read path.
    EXPECT_EQ(command("raycreate room 1"), 0);
    ray::Cloud cloud;
    EXPECT_TRUE(cloud.load("room.las"));
    EXPECT_FALSE(cloud.ends.empty());
    EXPECT_FALSE(cloud.starts.empty());
    EXPECT_EQ(cloud.ends.size(), cloud.starts.size());
    size_t bounded = 0;
    for (const auto &c : cloud.colours)
      if (c.alpha > 0) ++bounded;
    EXPECT_GT(bounded, 0u);
  }
#endif

#if RAYLIB_WITH_LAS
  TEST(RayLas, BeamId)
  {
    EXPECT_EQ(command("raycreate room 1"), 0);
    EXPECT_EQ(command("rayimport room.las 0,0,0 --beam_id"), 0);

    std::vector<Eigen::Vector3d> starts, ends;
    std::vector<double> times;
    std::vector<ray::RGBA> colours;
    std::vector<int32_t> beam_ids;
    size_t num_bounded = 0;
    bool ok = ray::readLas("room_raycloud.las",
      [&](std::vector<Eigen::Vector3d> &s, std::vector<Eigen::Vector3d> &e,
          std::vector<double> &t, std::vector<ray::RGBA> &c) {
        starts.insert(starts.end(), s.begin(), s.end());
        ends.insert(ends.end(), e.begin(), e.end());
        times.insert(times.end(), t.begin(), t.end());
        colours.insert(colours.end(), c.begin(), c.end());
      },
      num_bounded, 100.0, nullptr, ray::computeReadChunkSize(),
      nullptr, nullptr, nullptr, nullptr, nullptr, &beam_ids);

    EXPECT_TRUE(ok);
    EXPECT_FALSE(beam_ids.empty());
    EXPECT_EQ(beam_ids.size(), ends.size());
    // At least some beam IDs must be non-negative (valid assignment).
    bool has_valid = false;
    for (auto id : beam_ids)
      if (id >= 0) { has_valid = true; break; }
    EXPECT_TRUE(has_valid);
    // Beam IDs must be non-decreasing (same pulse gets same id, new pulse increments).
    for (size_t i = 1; i < beam_ids.size(); ++i)
      EXPECT_GE(beam_ids[i], beam_ids[i - 1]);
  }
#endif

#if RAYLIB_WITH_QHULL && RAYLIB_WITH_LAS
  TEST(RayLas, LabelledRoundTrip)
  {
    EXPECT_EQ(command("raycreate forest 2"), 0);
    EXPECT_EQ(command("rayextract terrain forest.las"), 0);
    EXPECT_EQ(command("rayextract segment forest.las --ground forest_mesh.ply"), 0);

    ray::Cloud seg;
    EXPECT_TRUE(seg.load("forest_segmented.las"));
    EXPECT_FALSE(seg.tree_ids.empty());
    EXPECT_EQ(seg.tree_ids.size(), seg.ends.size());
    bool has_labelled = false;
    for (auto tid : seg.tree_ids)
      if (tid != -1) { has_labelled = true; break; }
    EXPECT_TRUE(has_labelled);
    EXPECT_FALSE(seg.stem_ids.empty());
    EXPECT_EQ(seg.stem_ids.size(), seg.ends.size());
  }
#endif

#if RAYLIB_WITH_LAS
  TEST(RayVoxel, Smoke)
  {
    // raycreate writes .las; feed it directly to rayvoxel (it's already a raycloud).
    EXPECT_EQ(command("raycreate forest 2"), 0);
    EXPECT_EQ(command("rayvoxel forest.las --voxel_size 0.5"), 0);
    std::ifstream vox("forest.vox");
    EXPECT_TRUE(vox.is_open());
    vox.seekg(0, std::ios::end);
    EXPECT_GT(static_cast<long>(vox.tellg()), 0);
  }
#endif

#if RAYLIB_WITH_QHULL && RAYLIB_WITH_LAS
  TEST(RayCombine, LabelUnion)
  {
    EXPECT_EQ(command("raycreate forest 2"), 0);
    EXPECT_EQ(command("rayextract terrain forest.las"), 0);
    EXPECT_EQ(command("rayextract segment forest.las --ground forest_mesh.ply"), 0);
    EXPECT_EQ(copy("forest_segmented.las forest_segmented2.las"), 0);
    EXPECT_EQ(command("raycombine all forest_segmented.las forest_segmented2.las"), 0);

    ray::Cloud combined;
    EXPECT_TRUE(combined.load("forest_segmented_combined.las"));
    EXPECT_FALSE(combined.ends.empty());
    EXPECT_FALSE(combined.tree_ids.empty());
    EXPECT_EQ(combined.tree_ids.size(), combined.ends.size());
    // Combined should have roughly twice the points of the original.
    ray::Cloud original;
    EXPECT_TRUE(original.load("forest_segmented.las"));
    EXPECT_GT(combined.ends.size(), original.ends.size());
  }
#endif

  // Verify that an unbound (miss) beam never registers a hit, but does mark
  // voxels as observed (free-space traversal).
  TEST(RayVoxel, UnboundBeamNoHits)
  {
    const double voxel_size = 1.0;
    // weighting_method_ is stored as const std::string& — must outlive the processor.
    const std::string weighting = "full";
    ray::Cuboid bounds(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(10, 10, 10));
    ray::VoxelProcessor vp(bounds, voxel_size, weighting,
                           /*use_occlusion_rays=*/false,
                           /*use_flat_top=*/false, /*peaks=*/nullptr,
                           /*calc_beam_metrics=*/false, 0.0, 0.0,
                           /*subvoxel_split=*/0, /*dtm=*/nullptr);

    ray::BeamData beam;
    beam.beam_origin = Eigen::Vector3d(0.5, 0.5, 0.5);
    beam.gps_time    = 0.0;
    beam.num_returns = 1;
    beam.returns[0].x = 8.5;  beam.returns[0].y = 0.5;  beam.returns[0].z = 0.5;
    beam.returns[0].beam_origin     = beam.beam_origin;
    beam.returns[0].return_number   = 1;
    beam.returns[0].number_of_returns = 1;
    beam.returns[0].distance_to_sensor = 8.0;
    beam.returns[0].bound = 0;  // unbound / floating endpoint

    vp.processBeam(beam);

    const ray::VoxelProcessor::Map& m = vp.getMap();
    EXPECT_FALSE(m.empty()) << "traversal should have populated at least one voxel";
    float total_hits     = 0.0f;
    float total_observed = 0.0f;
    float total_unbound_rays = 0.0f;
    float total_path_length_unbound = 0.0f;
    for (const auto& kv : m)
    {
      total_hits     += kv.second.num_hits;
      total_observed += kv.second.num_beams_weighted;
      total_unbound_rays += kv.second.num_unbound_rays;
      total_path_length_unbound += kv.second.path_length_unbound;
    }
    EXPECT_EQ(total_hits, 0.0f)  << "unbound endpoint must not register as a hit";
    EXPECT_GT(total_observed, 0.0f) << "ray path must be marked as observed/free";
    EXPECT_GT(total_unbound_rays, 0.0f)      << "unbound ray must increment num_unbound_rays";
    EXPECT_GT(total_path_length_unbound, 0.0f) << "unbound ray must accumulate path_length_unbound";
  }

  // Verify that a bound (hit) beam registers exactly one endpoint hit and never
  // populates the unbound accumulators.
  TEST(RayVoxel, BoundBeamPathLength)
  {
    const double voxel_size = 1.0;
    const std::string weighting = "full";
    ray::Cuboid bounds(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(3, 1, 1));
    ray::VoxelProcessor vp(bounds, voxel_size, weighting,
                           /*use_occlusion_rays=*/false,
                           /*use_flat_top=*/false, /*peaks=*/nullptr,
                           /*calc_beam_metrics=*/false, 0.0, 0.0,
                           /*subvoxel_split=*/0, /*dtm=*/nullptr);

    ray::BeamData beam;
    beam.beam_origin = Eigen::Vector3d(0.5, 0.5, 0.5);
    beam.gps_time    = 0.0;
    beam.num_returns = 1;
    beam.returns[0].x = 2.5;  beam.returns[0].y = 0.5;  beam.returns[0].z = 0.5;
    beam.returns[0].beam_origin     = beam.beam_origin;
    beam.returns[0].return_number   = 1;
    beam.returns[0].number_of_returns = 1;
    beam.returns[0].distance_to_sensor = 2.0;
    beam.returns[0].bound = 1;  // bound / real return

    vp.processBeam(beam);

    const ray::VoxelProcessor::Map& m = vp.getMap();
    EXPECT_FALSE(m.empty()) << "traversal should have populated at least one voxel";
    float total_hits     = 0.0f;
    float total_path_length_observed = 0.0f;
    float total_unbound_rays = 0.0f;
    float total_path_length_unbound = 0.0f;
    int   hit_voxels = 0;
    for (const auto& kv : m)
    {
      total_hits += kv.second.num_hits;
      total_path_length_observed += kv.second.path_length_observed;
      total_unbound_rays += kv.second.num_unbound_rays;
      total_path_length_unbound += kv.second.path_length_unbound;
      if (kv.second.num_hits == 1.0f) ++hit_voxels;
    }
    EXPECT_EQ(hit_voxels, 1) << "exactly one (endpoint) voxel should register a hit";
    EXPECT_EQ(total_hits, 1.0f) << "a single bound return is exactly one hit";
    EXPECT_GT(total_path_length_observed, 0.0f) << "traversal must accumulate observed path length";
    EXPECT_EQ(total_unbound_rays, 0.0f) << "bound ray must not populate num_unbound_rays";
    EXPECT_EQ(total_path_length_unbound, 0.0f) << "bound ray must not populate path_length_unbound";
  }

  // Verify multi-return weighting: with weighting="equal" each return contributes
  // weight 1/N per voxel traversal; with weighting="full" each contributes 1.0.
  TEST(RayVoxel, MultiReturnWeighting)
  {
    const double voxel_size = 1.0;
    ray::Cuboid bounds(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(10, 1, 1));

    ray::BeamData beam;
    beam.beam_origin = Eigen::Vector3d(0.5, 0.5, 0.5);
    beam.gps_time    = 0.0;
    beam.num_returns = 2;
    beam.returns[0].x = 3.5;  beam.returns[0].y = 0.5;  beam.returns[0].z = 0.5;
    beam.returns[0].beam_origin       = beam.beam_origin;
    beam.returns[0].return_number     = 1;
    beam.returns[0].number_of_returns = 2;
    beam.returns[0].distance_to_sensor = 3.0;
    beam.returns[0].bound = 0;
    beam.returns[1].x = 7.5;  beam.returns[1].y = 0.5;  beam.returns[1].z = 0.5;
    beam.returns[1].beam_origin       = beam.beam_origin;
    beam.returns[1].return_number     = 2;
    beam.returns[1].number_of_returns = 2;
    beam.returns[1].distance_to_sensor = 7.0;
    beam.returns[1].bound = 0;

    // weighting="equal": per-voxel num_beams_weighted is a multiple of 1/N = 0.5.
    {
      const std::string weighting = "equal";
      ray::VoxelProcessor vp(bounds, voxel_size, weighting,
                             /*use_occlusion_rays=*/false,
                             /*use_flat_top=*/false, /*peaks=*/nullptr,
                             /*calc_beam_metrics=*/false, 0.0, 0.0,
                             /*subvoxel_split=*/0, /*dtm=*/nullptr);
      vp.processBeam(beam);
      const ray::VoxelProcessor::Map& m = vp.getMap();
      EXPECT_FALSE(m.empty()) << "traversal should have populated at least one voxel";
      for (const auto& kv : m)
      {
        EXPECT_LT(std::fmod(kv.second.num_beams_weighted, 0.5f), 1e-4f)
          << "equal weighting must contribute multiples of 1/N (0.5) per voxel";
      }
    }

    // weighting="full": per-voxel num_beams_weighted is a multiple of 1.0.
    {
      const std::string weighting = "full";
      ray::VoxelProcessor vp(bounds, voxel_size, weighting,
                             /*use_occlusion_rays=*/false,
                             /*use_flat_top=*/false, /*peaks=*/nullptr,
                             /*calc_beam_metrics=*/false, 0.0, 0.0,
                             /*subvoxel_split=*/0, /*dtm=*/nullptr);
      vp.processBeam(beam);
      const ray::VoxelProcessor::Map& m = vp.getMap();
      EXPECT_FALSE(m.empty()) << "traversal should have populated at least one voxel";
      for (const auto& kv : m)
      {
        EXPECT_LT(std::fmod(kv.second.num_beams_weighted, 1.0f), 1e-4f)
          << "full weighting must contribute multiples of 1.0 per voxel";
      }
    }
  }

  // --- RayVoxelAttenuation: Voxel::pad_g0_5 ---

  // Guard: fewer than 2 observed rays returns 0 regardless of hits/path.
  TEST(RayVoxelAttenuation, PadG05Guard)
  {
    ray::VoxelGrid::Voxel v{};
    v.num_beams_weighted = 1.0f;
    v.num_hits = 1.0f;
    v.path_length_observed = 1.0f;
    EXPECT_EQ(v.pad_g0_5(), 0.0);
  }

  // Zero hits produces zero PAD regardless of the observed-ray count.
  TEST(RayVoxelAttenuation, PadG05ZeroHits)
  {
    ray::VoxelGrid::Voxel v{};
    v.num_beams_weighted = 10.0f;
    v.num_hits = 0.0f;
    v.path_length_observed = 5.0f;
    EXPECT_EQ(v.pad_g0_5(), 0.0);
  }

  // N=4, H=2, L=3.0: 2*(3/4)*(2/3.0) = 1.0 (eps in denominator is negligible).
  TEST(RayVoxelAttenuation, PadG05KnownValue)
  {
    ray::VoxelGrid::Voxel v{};
    v.num_beams_weighted = 4.0f;
    v.num_hits = 2.0f;
    v.path_length_observed = 3.0f;
    EXPECT_NEAR(v.pad_g0_5(), 1.0, 1e-6);
  }

  // N=10, H=1, L=2.0: 2*(9/10)*(1/2.0) = 0.9.
  TEST(RayVoxelAttenuation, PadG05SingleHit)
  {
    ray::VoxelGrid::Voxel v{};
    v.num_beams_weighted = 10.0f;
    v.num_hits = 1.0f;
    v.path_length_observed = 2.0f;
    EXPECT_NEAR(v.pad_g0_5(), 0.9, 1e-6);
  }

  // --- RayVoxelAttenuation: Voxel::transmittance ---

  // Guard: bs_entering below threshold returns 1.0 (no interception data).
  TEST(RayVoxelAttenuation, TransmittanceGuard)
  {
    ray::VoxelGrid::Voxel v{};
    v.bs_entering = 0.0f;
    v.bs_intercepted = 0.0f;
    EXPECT_EQ(v.transmittance(), 1.0);
  }

  // Full transmittance: nothing intercepted.
  TEST(RayVoxelAttenuation, TransmittanceFull)
  {
    ray::VoxelGrid::Voxel v{};
    v.bs_entering = 1.0f;
    v.bs_intercepted = 0.0f;
    EXPECT_NEAR(v.transmittance(), 1.0, 1e-10);
  }

  // Full interception: entering equals intercepted → 0.
  TEST(RayVoxelAttenuation, TransmittanceZero)
  {
    ray::VoxelGrid::Voxel v{};
    v.bs_entering = 1.0f;
    v.bs_intercepted = 1.0f;
    EXPECT_NEAR(v.transmittance(), 0.0, 1e-10);
  }

  // Partial: (4-1)/4 = 0.75.
  TEST(RayVoxelAttenuation, TransmittancePartial)
  {
    ray::VoxelGrid::Voxel v{};
    v.bs_entering = 4.0f;
    v.bs_intercepted = 1.0f;
    EXPECT_NEAR(v.transmittance(), 0.75, 1e-10);
  }

  // Over-interception clamps to 0 via max(0,...).
  TEST(RayVoxelAttenuation, TransmittanceOverInterception)
  {
    ray::VoxelGrid::Voxel v{};
    v.bs_entering = 1.0f;
    v.bs_intercepted = 2.0f;
    EXPECT_NEAR(v.transmittance(), 0.0, 1e-10);
  }

  // --- RayVoxelAttenuation: solveBaileyPadEq10 ---

  // G=0: function must return 0.
  TEST(RayVoxelAttenuation, BaileyGuardZeroG)
  {
    EXPECT_EQ(ray::solveBaileyPadEq10(2.0, 5.0, 2.0, 0.0), 0.0);
  }

  // num_rays < 1: function must return 0.
  TEST(RayVoxelAttenuation, BaileyGuardFewRays)
  {
    EXPECT_EQ(ray::solveBaileyPadEq10(2.0, 0.5, 0.0, 0.5), 0.0);
  }

  // r_bar = path_length/num_rays = 0: function must return 0.
  TEST(RayVoxelAttenuation, BaileyGuardZeroRBar)
  {
    EXPECT_EQ(ray::solveBaileyPadEq10(0.0, 5.0, 2.0, 0.5), 0.0);
  }

  // Low attenuation: a_L=0.5, G=0.5, r_bar=1.0.
  // P_bar = exp(-0.25); num_rays=100, path=100, num_hits=(1-P_bar)*100.
  TEST(RayVoxelAttenuation, BaileyLowAttenuation)
  {
    const double a_L = 0.5, G = 0.5, r_bar = 1.0;
    const double P_bar = std::exp(-a_L * G * r_bar);
    const double num_rays = 100.0;
    const double path_length = num_rays * r_bar;
    const double num_hits = (1.0 - P_bar) * num_rays;
    EXPECT_NEAR(ray::solveBaileyPadEq10(path_length, num_rays, num_hits, G), a_L, 1e-6);
  }

  // High attenuation: a_L=5.0, G=0.5, r_bar=1.0.
  // P_bar = exp(-2.5); num_rays=100, path=100, num_hits=(1-P_bar)*100.
  TEST(RayVoxelAttenuation, BaileyHighAttenuation)
  {
    const double a_L = 5.0, G = 0.5, r_bar = 1.0;
    const double P_bar = std::exp(-a_L * G * r_bar);
    const double num_rays = 100.0;
    const double path_length = num_rays * r_bar;
    const double num_hits = (1.0 - P_bar) * num_rays;
    EXPECT_NEAR(ray::solveBaileyPadEq10(path_length, num_rays, num_hits, G), a_L, 1e-4);
  }

  // Unity G: a_L=2.0, G=1.0, r_bar=0.5.
  // P_bar = exp(-1.0); num_rays=50, path=25, num_hits=(1-P_bar)*50.
  TEST(RayVoxelAttenuation, BaileyUnityG)
  {
    const double a_L = 2.0, G = 1.0, r_bar = 0.5;
    const double P_bar = std::exp(-a_L * G * r_bar);
    const double num_rays = 50.0;
    const double path_length = num_rays * r_bar;
    const double num_hits = (1.0 - P_bar) * num_rays;
    EXPECT_NEAR(ray::solveBaileyPadEq10(path_length, num_rays, num_hits, G), a_L, 1e-4);
  }

  // --- RayVoxelAttenuation: buildTriangleInclinationHistograms ---

  // A horizontal triangle (z=0 plane): area=0.5 is recorded correctly.
  // bailey_g_leaf stays 0 because the eq.(4) weight is area*sin(theta)=0 at theta=0,
  // so the denominator guard fires and G_bar is never written.
  TEST(RayVoxelAttenuation, TriHistHorizontalLeaf)
  {
    std::vector<Eigen::Vector3d> positions = {
      Eigen::Vector3d(0.0, 0.0, 0.0),
      Eigen::Vector3d(1.0, 0.0, 0.0),
      Eigen::Vector3d(0.0, 1.0, 0.0),
    };
    // K=2 neighbours per point; knn_indices is (2,3): column i lists neighbours of point i.
    Eigen::MatrixXi knn(2, 3);
    knn(0, 0) = 1;  knn(1, 0) = 2;
    knn(0, 1) = 0;  knn(1, 1) = 2;
    knn(0, 2) = 0;  knn(1, 2) = 1;
    std::vector<int64_t> flat_indices = { 0, 0, 0 };
    std::vector<int> class_labels = { 1, 1, 1 };
    auto result = ray::buildTriangleInclinationHistograms(positions, knn, flat_indices, class_labels, 4, 2.0);
    EXPECT_FALSE(result.empty()) << "horizontal leaf facet must produce a histogram entry";
    auto it = result.find(0);
    ASSERT_NE(it, result.end());
    EXPECT_FALSE(it->second.tiad_leaf.empty());
    EXPECT_NEAR(it->second.total_leaf_area, 0.5, 1e-6)
      << "area of unit right triangle = 0.5";
    // sin(theta)=0 for horizontal facet → weight=0 → denominator guard → G_bar=0.
    EXPECT_NEAR(it->second.bailey_g_leaf, 0.0, 1e-6)
      << "horizontal facet has zero eq.(4) weight; G_bar denominator guard must leave it at 0";
  }

  // A 45-degree tilted triangle has normal (1/√2, 0, 1/√2), theta=pi/4,
  // G_i=|r_hat·n_hat|=1/√2, weight=area*sin(pi/4)>0, so G_bar = 1/√2 ≈ 0.707.
  // Vertices: p0=(0,0,0), p1=(0,1,0), p2=(-1/√2, 0, 1/√2) — all within l_max=2.
  TEST(RayVoxelAttenuation, TriHist45DegreeLeaf)
  {
    const double s = 1.0 / std::sqrt(2.0);
    std::vector<Eigen::Vector3d> positions = {
      Eigen::Vector3d(0.0, 0.0, 0.0),
      Eigen::Vector3d(0.0, 1.0, 0.0),
      Eigen::Vector3d(-s,  0.0,  s),
    };
    Eigen::MatrixXi knn(2, 3);
    knn(0, 0) = 1;  knn(1, 0) = 2;
    knn(0, 1) = 0;  knn(1, 1) = 2;
    knn(0, 2) = 0;  knn(1, 2) = 1;
    std::vector<int64_t> flat_indices = { 0, 0, 0 };
    std::vector<int> class_labels = { 1, 1, 1 };
    auto result = ray::buildTriangleInclinationHistograms(positions, knn, flat_indices, class_labels, 4, 2.0);
    auto it = result.find(0);
    ASSERT_NE(it, result.end());
    EXPECT_NEAR(it->second.bailey_g_leaf, s, 1e-6)
      << "45-degree facet: G_bar = G_i = |r_hat·n_hat| = 1/sqrt(2)";
  }

  // A single vertical triangle (yz-plane) has normal (1,0,0), G_i=0.0,
  // so bailey_g_leaf should be 0.0.
  TEST(RayVoxelAttenuation, TriHistVerticalLeaf)
  {
    std::vector<Eigen::Vector3d> positions = {
      Eigen::Vector3d(0.0, 0.0, 0.0),
      Eigen::Vector3d(0.0, 1.0, 0.0),
      Eigen::Vector3d(0.0, 0.0, 1.0),
    };
    Eigen::MatrixXi knn(2, 3);
    knn(0, 0) = 1;  knn(1, 0) = 2;
    knn(0, 1) = 0;  knn(1, 1) = 2;
    knn(0, 2) = 0;  knn(1, 2) = 1;
    std::vector<int64_t> flat_indices = { 0, 0, 0 };
    std::vector<int> class_labels = { 1, 1, 1 };
    auto result = ray::buildTriangleInclinationHistograms(positions, knn, flat_indices, class_labels, 4, 2.0);
    EXPECT_FALSE(result.empty()) << "vertical leaf facet must produce a histogram entry";
    auto it = result.find(0);
    ASSERT_NE(it, result.end());
    EXPECT_NEAR(it->second.bailey_g_leaf, 0.0, 1e-6)
      << "vertical facet normal is horizontal so G_i = |r_hat.(1,0,0)| = 0";
  }

  // Unknown class (label=0) points are skipped — result must be empty.
  TEST(RayVoxelAttenuation, TriHistUnknownClassSkipped)
  {
    std::vector<Eigen::Vector3d> positions = {
      Eigen::Vector3d(0.0, 0.0, 0.0),
      Eigen::Vector3d(1.0, 0.0, 0.0),
      Eigen::Vector3d(0.0, 1.0, 0.0),
    };
    Eigen::MatrixXi knn(2, 3);
    knn(0, 0) = 1;  knn(1, 0) = 2;
    knn(0, 1) = 0;  knn(1, 1) = 2;
    knn(0, 2) = 0;  knn(1, 2) = 1;
    std::vector<int64_t> flat_indices = { 0, 0, 0 };
    std::vector<int> class_labels = { 0, 0, 0 };  // all unknown
    auto result = ray::buildTriangleInclinationHistograms(positions, knn, flat_indices, class_labels, 4, 2.0);
    EXPECT_TRUE(result.empty()) << "unknown-class points must not produce any facets";
  }

} // raytest
