// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <iomanip> // For std::setprecision
#include "raylib/raycloud.h"
#include "raylib/raylaz.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"
#include "raylib/raytrajectory.h"
#include "raylib/extraction/raytrees.h" // For convertIntToColour

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Export a ray cloud into a point cloud amd trajectory file" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayexport raycloudfile.ply pointcloud.ply/.laz/.las/.txt/.xyz trajectoryfile.ply/.txt - output in the chosen point cloud and trajectory formats" << std::endl;
  std::cout << "                           --traj_delta 0.1 - trajectory temporal decimation period in s. Default is 0.1" << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayExport(int argc, char *argv[])
{
  ray::FileArgument raycloud_file, pointcloud_file, trajectory_file;
  ray::DoubleArgument traj_delta(0.0, 10000, 0.1);
  ray::OptionalKeyValueArgument delta_option("traj_delta", 't', &traj_delta);
  if (!ray::parseCommandLine(argc, argv, { &raycloud_file, &pointcloud_file, &trajectory_file }, { &delta_option }))
    usage();

  // --- START OF MODIFICATION: Handle flexible output ---
  const std::string ext = pointcloud_file.nameExt();

  if (ext == "laz" || ext == "las")
  {
    ray::LasWriter las_writer(pointcloud_file.name());
    auto add_chunk = [&las_writer](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                                   std::vector<double> &times, std::vector<ray::RGBA> &colours,
                                   std::vector<uint8_t> &classifications, std::vector<uint16_t> &branch_ids) {
        
        // If branch_ids were not in the input file, fall back to converting from color
        if (branch_ids.empty() || (branch_ids.size() > 0 && branch_ids[0] == 0)) {
            branch_ids.resize(colours.size());
            for(size_t i = 0; i < colours.size(); ++i) {
                int id = ray::convertColourToInt(colours[i]);
                branch_ids[i] = (id >= 0) ? static_cast<uint16_t>(id) : 0;
            }
        }
        las_writer.writeChunk(ends, times, colours, classifications, branch_ids); 
    };
    if (!ray::readPly(raycloud_file.name(), true, add_chunk, 0))
      usage();
  }
  else if (ext == "ply")
  {
    std::ofstream ofs;
    bool header_written = false;
    std::vector<char> buffer;
    size_t vertex_byte_size = 0;
    
    auto add_chunk = [&](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                         std::vector<double> &times, std::vector<ray::RGBA> &colours,
                         std::vector<uint8_t> &classifications, std::vector<uint16_t> &branch_ids) {
        if (!header_written) {
            bool has_classification = !classifications.empty();
            bool has_point_source_id = !branch_ids.empty();
            ray::writePointCloudChunkStart(pointcloud_file.name(), ofs, has_classification, has_point_source_id);

            #if RAYLIB_DOUBLE_RAYS
            vertex_byte_size += sizeof(double) * 3;
            #else
            vertex_byte_size += sizeof(float) * 3;
            #endif
            vertex_byte_size += sizeof(double); // time
            vertex_byte_size += sizeof(ray::RGBA);
            if (has_classification) vertex_byte_size += sizeof(uint8_t);
            if (has_point_source_id) vertex_byte_size += sizeof(uint16_t);

            header_written = true;
        }

        bool has_warned = false;
        ray::writePointCloudChunk(ofs, buffer, ends, times, colours, classifications, branch_ids, has_warned);
    };
    if (!ray::readPly(raycloud_file.name(), true, add_chunk, 0))
      usage();
      
    if (header_written) {
        ray::writePointCloudChunkEnd(ofs, vertex_byte_size);
    }
  }
  else if (ext == "xyz" || ext == "txt")
  {
    std::ofstream ofs;
    ofs.open(pointcloud_file.name(), std::ios::out);
    if (ofs.fail()) usage();
    
    ofs << std::setprecision(4) << std::fixed;
    bool txt = (ext == "txt");
    if (txt)
    {
      ofs << "# x,y,z,time,red,green,blue,alpha,classification,branch_id" << std::endl;
    }

    auto add_chunk = [&ofs, txt](std::vector<Eigen::Vector3d> &, std::vector<Eigen::Vector3d> &ends,
                                 std::vector<double> &times, std::vector<ray::RGBA> &colours,
                                 std::vector<uint8_t> &classifications, std::vector<uint16_t> &branch_ids) {
      for (size_t i = 0; i < ends.size(); i++)
      {
        if (txt)
        {
          ofs << ends[i][0] << "," << ends[i][1] << "," << ends[i][2] << "," << times[i] << "," 
              << (int)colours[i].red << "," << (int)colours[i].green << "," << (int)colours[i].blue << "," << (int)colours[i].alpha << ","
              << (int)classifications[i] << "," << branch_ids[i] << std::endl;
        }
        else // xyz format, just coordinates
        {
          ofs << ends[i][0] << " " << ends[i][1] << " " << ends[i][2] << std::endl;
        }
      }    
    };
    if (!ray::readPly(raycloud_file.name(), true, add_chunk, 0))
      usage();
    ofs.close();
  }
  else
  {
    usage();
  }
  // --- END OF MODIFICATION ---

  const double time_step = traj_delta.value();
  std::set<int64_t> time_slots;
  int64_t last_time_slot = std::numeric_limits<int64_t>::min();

  if (trajectory_file.nameExt() == "ply")
  {
    std::ofstream ofs;
    // Trajectory PLY files are simple point clouds, so we write without extra fields for compatibility.
    ray::writePointCloudChunkStart(trajectory_file.name(), ofs, false, false);
    std::vector<char> buffer;
    ray::Cloud chunk;
    bool has_warned = false;
    size_t vertex_byte_size = sizeof(float) * 3 + sizeof(double) + sizeof(ray::RGBA);
    #if RAYLIB_DOUBLE_RAYS
    vertex_byte_size = sizeof(double) * 3 + sizeof(double) + sizeof(ray::RGBA);
    #endif

    auto decimate_time = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                             std::vector<double> &times, std::vector<ray::RGBA> &colours,
                             std::vector<uint8_t>&, std::vector<uint16_t>&) {
      chunk.clear();
      for (size_t i = 0; i < ends.size(); i++)
      {
        const int64_t time_slot = static_cast<int64_t>(std::floor(times[i] / time_step));
        if (time_slot == last_time_slot) continue;
        if (time_slots.insert(time_slot).second)
        {
          chunk.starts.push_back(starts[i]); // Using starts as the points for trajectory
          chunk.times.push_back(times[i]);
          chunk.colours.push_back(colours[i]);
        }
        last_time_slot = time_slot;
      }
      ray::writePointCloudChunk(ofs, buffer, chunk.starts, chunk.times, chunk.colours, {}, {}, has_warned);
    };
    if (!ray::readPly(raycloud_file.name(), true, decimate_time, 0))
      usage();
    ray::writePointCloudChunkEnd(ofs, vertex_byte_size);
  }
  else if (trajectory_file.nameExt() == "txt")
  {
    std::vector<ray::TrajectoryNode> traj_nodes;
    bool sorted = true;

    auto decimate_time = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                             std::vector<double> &times, std::vector<ray::RGBA> &,
                             std::vector<uint8_t>&, std::vector<uint16_t>&) {
      for (size_t i = 0; i < ends.size(); i++)
      {
        const int64_t time_slot = static_cast<int64_t>(std::floor(times[i] / time_step));
        if (time_slot == last_time_slot) continue;
        
        if (time_slots.insert(time_slot).second)
        {
          if (!traj_nodes.empty() && times[i] < traj_nodes.back().time) sorted = false;
          
          ray::TrajectoryNode traj_node;
          traj_node.time = times[i];
          traj_node.point = starts[i];
          traj_nodes.push_back(traj_node);
        }
        last_time_slot = time_slot;
      }
    };
    if (!ray::readPly(raycloud_file.name(), true, decimate_time, 0))
    {
      usage();
    }

    if (!sorted)
    {
      std::sort(traj_nodes.begin(), traj_nodes.end(),
                [](const ray::TrajectoryNode &a, const ray::TrajectoryNode &b) { return a.time < b.time; });
    }

    ray::saveTrajectory(traj_nodes, trajectory_file.name());
  }
  else
    usage();
    
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayExport, argc, argv);
}