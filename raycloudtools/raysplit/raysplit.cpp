// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymesh.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/raysplitter.h"
#include "raylib/rayforeststructure.h"
// --- START OF MODIFICATION ---
#include "raylib/raylaz.h"
// --- END OF MODIFICATION ---

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Split a ray cloud relative to the supplied geometry, generating two cropped ray clouds" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raysplit raycloud plane 10,0,0           - splits around plane at 10 m along x axis" << std::endl;
  std::cout << "                  colour                 - splits by colour, one cloud per colour" << std::endl;
  std::cout << "                  colour 0.5,0,0         - splits by colour, around half red component" << std::endl;
  std::cout << "                  single_colour 255,0,0  - splits out a single colour, in 0-255 units" << std::endl;
  std::cout << "                  seg_colour             - splits to one cloud per colour, converting _segmented.ply colours to their index suffix" << std::endl;
  std::cout << "                  alpha 0.0              - splits out unbounded rays, which have zero intensity" << std::endl;
  std::cout << "                  file distance 0.2      - splits raycloud at 0.2m from the (ply mesh or trees) file surface" << std::endl;
  std::cout << "                  raydir 0,0,0.8         - splits based on ray direction, here around nearly vertical rays" << std::endl;
  std::cout << "                  range 10               - splits out rays more than 10 m long" << std::endl;
  std::cout << "                  time 1000 (or time 3 %)- splits at given time stamp (or percentage along)" << std::endl;
  std::cout << "                  box x,y,z rx,ry,rz     - splits around a given XYZ centred axis-aligned box of the given radii" << std::endl;  
  std::cout << "                  gap 0.1                - splits into largest cloud connected within this gap, and the remainder." << std::endl;
  std::cout << "                  grid wx,wy,wz          - splits into a 0,0,0 centred grid of files, cell width wx,wy,wz. 0 for unused axes." << std::endl;
  std::cout << "                  grid wx,wy,wz 1        - same as above, but with a 1 metre overlap between cells." << std::endl;
  std::cout << "                  grid wx,wy,wz,wt       - splits into a grid of files, cell width wx,wy,wz and period wt. 0 for unused axes." << std::endl;
  std::cout << "                  capsule 1,2,3 10,11,12 5  - splits within a capsule using start, end and radius" << std::endl;
  // --- START OF MODIFICATION ---
  std::cout << "                  --output_format las     - save output files as .las or .laz instead of .ply" << std::endl;
  // --- END OF MODIFICATION ---
  // clang-format on
  exit(exit_code);
}

// Decimates the ray cloud, spatially or in time
int raySplit(int argc, char *argv[])
{
  ray::FileArgument cloud_file;
  double max_val = std::numeric_limits<double>::max();
  ray::Vector3dArgument plane, colour(0.0, 1.0), single_colour(0.0, 255.0), raydir(-1.0, 1.0),
    box_centre, box_radius(0.0, max_val), cell_width(0.0, max_val), capsule_start, capsule_end;
  ray::Vector4dArgument cell_width2(0.0, max_val);
  ray::DoubleArgument overlap(0.0, 10000.0);
  ray::DoubleArgument time, alpha(0.0, 1.0), range(0.0, 1000.0), capsule_radius(0.001, 1000.0), gap(0.000001, 10000.0);
  ray::KeyValueChoice choice({ "plane", "time", "colour", "single_colour", "alpha", "raydir", "range", "gap" },
                             { &plane, &time, &colour, &single_colour, &alpha, &raydir, &range, &gap });
  ray::FileArgument mesh_file, tree_file;
  ray::TextArgument distance_text("distance"), time_text("time"), percent_text("%");
  ray::TextArgument box_text("box"), grid_text("grid"), colour_text("colour"), seg_colour_text("seg_colour"), capsule_text("capsule");
  ray::DoubleArgument mesh_offset;
  // --- START OF FIX: Replaced TextArgument with a flexible FileArgument ---
  ray::FileArgument output_format(false); // Use FileArgument, disable extension checking
  output_format.name() = "ply";           // Manually set the default value
  ray::OptionalKeyValueArgument output_format_option("output_format", 'f', &output_format);
  // --- END OF FIX ---

  bool standard_format = ray::parseCommandLine(argc, argv, { &cloud_file, &choice }, {&output_format_option});
  bool colour_format = ray::parseCommandLine(argc, argv, { &cloud_file, &colour_text }, {&output_format_option});
  bool seg_colour_format = ray::parseCommandLine(argc, argv, { &cloud_file, &seg_colour_text }, {&output_format_option});
  bool time_percent = ray::parseCommandLine(argc, argv, { &cloud_file, &time_text, &time, &percent_text }, {&output_format_option});
  bool box_format = ray::parseCommandLine(argc, argv, { &cloud_file, &box_text, &box_centre, &box_radius }, {&output_format_option});
  bool grid_format = ray::parseCommandLine(argc, argv, { &cloud_file, &grid_text, &cell_width }, {&output_format_option});
  bool grid_format2 = ray::parseCommandLine(argc, argv, { &cloud_file, &grid_text, &cell_width2 }, {&output_format_option});
  bool grid_format3 = ray::parseCommandLine(argc, argv, { &cloud_file, &grid_text, &cell_width, &overlap }, {&output_format_option});
  bool mesh_split = ray::parseCommandLine(argc, argv, { &cloud_file, &mesh_file, &distance_text, &mesh_offset }, {&output_format_option});
  bool capsule_split =
    ray::parseCommandLine(argc, argv, { &cloud_file, &capsule_text, &capsule_start, &capsule_end, &capsule_radius }, {&output_format_option});

  if (!standard_format && !colour_format && !seg_colour_format && !box_format && !grid_format && !grid_format2 && !grid_format3 &&
      !mesh_split && !time_percent && !capsule_split)
  {
    usage();
  }

  // --- START OF MODIFICATION: Dynamic output filenames ---
  const std::string out_ext = std::string(".") + output_format.name();
  const std::string in_name = cloud_file.nameStub() + "_inside" + out_ext;
  const std::string out_name = cloud_file.nameStub() + "_outside" + out_ext;
  // --- END OF MODIFICATION ---
  const std::string rc_name = cloud_file.name();
  bool res = true;

  if (capsule_split)
  {
    res = ray::splitCapsule(rc_name, in_name, out_name, capsule_start.value(), capsule_end.value(), capsule_radius.value(), out_ext);
  }
  else if (colour_format)
  {
    res = ray::splitColour(cloud_file.name(), cloud_file.nameStub(), false, out_ext);
  } 
  else if (seg_colour_format)
  {
    res = ray::splitColour(cloud_file.name(), cloud_file.nameStub(), true, out_ext);
  }
  else if (mesh_split) 
  {
    if (mesh_file.nameExt() == "ply")
    {
      ray::Mesh mesh;
      ray::readPlyMesh(mesh_file.name(), mesh);
      if (!mesh.splitCloud(rc_name, mesh_offset.value(), in_name, out_name, out_ext))
      {
        usage();
      }
    }
    else if (mesh_file.nameExt() == "txt")
    {
      ray::ForestStructure forest;
      forest.load(mesh_file.name());
      ray::Cloud cloud;
      if (!cloud.load(rc_name)) usage();
      
      ray::Cloud inside, outside;
      forest.splitCloud(cloud, mesh_offset.value(), inside, outside);

      if (output_format.name() == "las" || output_format.name() == "laz") {
          ray::LasWriter inside_writer(in_name);
          inside_writer.writeChunk(inside.ends, inside.times, inside.colours, inside.classifications, inside.branch_ids);
          ray::LasWriter outside_writer(out_name);
          outside_writer.writeChunk(outside.ends, outside.times, outside.colours, outside.classifications, outside.branch_ids);
      } else {
          inside.save(in_name, true);
          outside.save(out_name, true);
      }
    }
  }
  else if (time_percent)
  {
    double min_time = std::numeric_limits<double>::max();
    double max_time = std::numeric_limits<double>::lowest();
    auto time_bounds = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &, std::vector<double> &times,
                           std::vector<ray::RGBA> &, std::vector<uint8_t>&, std::vector<uint16_t>&) {
      for (auto &time : times)
      {
        min_time = std::min(min_time, time);
        max_time = std::max(max_time, time);
      }
    };
    if (!ray::Cloud::read(cloud_file.name(), time_bounds))
      usage();
    
    const double time_thresh = min_time + (max_time - min_time) * time.value() / 100.0;
    res = ray::split(rc_name, in_name, out_name,
                     [&](const ray::Cloud &cloud, int i) -> bool { return cloud.times[i] > time_thresh; }, out_ext);
  }
  else if (box_format)
  {
    Eigen::Vector3d extents = box_radius.value();
    for (int i = 0; i<3; i++)
    {
      if (extents[i] == 0.0) extents[i] = 1e7;
    }
    res = ray::splitBox(rc_name, in_name, out_name, box_centre.value(), extents, out_ext);
  }
  else if (grid_format)
  {
    res = ray::splitGrid(rc_name, cloud_file.nameStub(), cell_width.value(), 0.0, out_ext);
  }
  else if (grid_format2)
  {
    res = ray::splitGrid(rc_name, cloud_file.nameStub(), cell_width2.value(), 0.0, out_ext);
  }
  else if (grid_format3)
  {
    res = ray::splitGrid(rc_name, cloud_file.nameStub(), cell_width.value(), overlap.value(), out_ext);
  }
  else
  {
    const std::string &parameter = choice.selectedKey();
    if (parameter == "time")
    {
      res = ray::split(rc_name, in_name, out_name,
                       [&](const ray::Cloud &cloud, int i) -> bool { return cloud.times[i] > time.value(); }, out_ext);
    }
    else if (parameter == "alpha")
    {
      uint8_t c = uint8_t(255.0 * alpha.value());
      res = ray::split(rc_name, in_name, out_name,
                       [&](const ray::Cloud &cloud, int i) -> bool { return cloud.colours[i].alpha > c; }, out_ext);
    }
    else if (parameter == "plane")
    {
      ray::splitPlane(rc_name, in_name, out_name, plane.value(), out_ext);
    }
    else if (parameter == "raydir")
    {
      Eigen::Vector3d vec = raydir.value() / raydir.value().squaredNorm();
      res = ray::split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool {
        Eigen::Vector3d ray_dir = (cloud.ends[i] - cloud.starts[i]).normalized();
        return ray_dir.dot(vec) > 1.0;
      }, out_ext);
    }
    else if (parameter == "colour")
    {
      Eigen::Vector3d vec = colour.value() / colour.value().squaredNorm();
      res = ray::split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool {
        Eigen::Vector3d col((double)cloud.colours[i].red / 255.0, (double)cloud.colours[i].green / 255.0,
                            (double)cloud.colours[i].blue / 255.0);
        return col.dot(vec) > 1.0;
      }, out_ext);
    }
    else if (parameter == "single_colour")
    {
      ray::RGBA col;
      col.red = (uint8_t)single_colour.value()[0];
      col.green = (uint8_t)single_colour.value()[1];
      col.blue = (uint8_t)single_colour.value()[2];
      res = ray::split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool {
        return !(cloud.colours[i].red == col.red && cloud.colours[i].green == col.green &&
                 cloud.colours[i].blue == col.blue);
      }, out_ext);
    }
    else if (parameter == "range")
    {
      res = ray::split(rc_name, in_name, out_name, [&](const ray::Cloud &cloud, int i) -> bool {
        return (cloud.starts[i] - cloud.ends[i]).norm() > range.value();
      }, out_ext);
    }
    else if (parameter == "gap")
    {
      ray::Cloud cloud;
      if (!cloud.load(rc_name)) usage();
      Eigen::Vector3d offset = cloud.removeStartPos();

      Eigen::MatrixXi neighbour_indices;
      int search_size = 10;
      cloud.getSurfels(search_size, nullptr, nullptr, nullptr, nullptr, &neighbour_indices, gap.value(), false);

      std::vector<bool> visited(neighbour_indices.cols(), false);
      std::vector<std::vector<int>> clusters;
      int largest_cluster = -1;
      int largest_cluster_size = 0;
      for (int i = 0; i<(int)neighbour_indices.cols(); i++)
      {
        if (visited[i]) continue;
        clusters.push_back(std::vector<int>());
        auto &cluster = clusters.back();
        cluster.push_back(i);
        visited[i] = true;
        for (size_t j = 0; j<cluster.size(); j++)
        {
          int id = cluster[j];
          for (int k = 0; k<search_size; k++)
          {
            int ind = neighbour_indices(k, id);
            if (ind == -1) break;
            if (!visited[ind])
            {
              cluster.push_back(ind);
              visited[ind] = true;
            }
          }
        }
        if ((int)cluster.size() > largest_cluster_size)
        {
          largest_cluster_size = (int)cluster.size();
          largest_cluster = (int)clusters.size()-1;
        }
      }
      cloud.translate(offset);

      for (auto &ind: clusters[largest_cluster])
        visited[ind] = false;

      ray::Cloud inside, outside;
      for (size_t i = 0; i<visited.size(); i++)
      {
        if (visited[i]) {
            outside.addRay(cloud, i);
        } else {
            inside.addRay(cloud, i);
        }
      }

      if (output_format.name() == "las" || output_format.name() == "laz") {
          ray::LasWriter inside_writer(in_name);
          inside_writer.writeChunk(inside.ends, inside.times, inside.colours, inside.classifications, inside.branch_ids);
          ray::LasWriter outside_writer(out_name);
          outside_writer.writeChunk(outside.ends, outside.times, outside.colours, outside.classifications, outside.branch_ids);
      } else {
          inside.save(in_name, true);
          outside.save(out_name, true);
      }
    }
  }
  if (!res)
    usage();
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(raySplit, argc, argv);
}