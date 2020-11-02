// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "raylib/raycloud.h"
#include "raylib/rayply.h"
#include "raylib/raylaz.h"
#include "raylib/rayparse.h"
#include "raylib/raytrajectory.h"

void usage(int exit_code = 1)
{
  std::cout << "Import a point cloud and trajectory file into a ray cloud" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "rayimport pointcloudfile trajectoryfile  - pointcloudfile can be a .laz, .las or .ply file" << std::endl;
  std::cout << "                                           trajectoryfile is a text file in time,x,y,z format"
       << std::endl;
  std::cout << "The output is a .ply file of the same name (or with suffix _raycloud if the input was a .ply file)." << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  ray::FileArgument cloud_file, trajectory_file;
  if (!ray::parseCommandLine(argc, argv, {&cloud_file, &trajectory_file}))
    usage();

  ray::Cloud cloud;
  std::string point_cloud = cloud_file.name();
  std::string traj_file = trajectory_file.name();

  // load the trajectory first, it should fit into main memory
  ray::Trajectory trajectory;
  std::string traj_end = traj_file.substr(traj_file.size() - 4);
  if (traj_end == ".ply")
  {
    std::vector<Eigen::Vector3d> starts;
    std::vector<Eigen::Vector3d> ends;
    std::vector<double> times;
    std::vector<ray::RGBA> colours;
    if (!ray::readPly(traj_file, starts, ends, times, colours, false))
      return false;
    trajectory.points_ = std::move(ends);
    trajectory.times_ = std::move(times);
  }
  else if (!trajectory.load(traj_file))
    usage();

  std::string save_file = cloud_file.nameStub();
  if (cloud_file.nameExt() == "ply")
    save_file += "_raycloud";


  std::string name_end = point_cloud.substr(point_cloud.size() - 4);
  ray::RayPlyBuffer buffer;
  std::ofstream ofs;
  if (!ray::writePlyChunkStart(save_file + ".ply", ofs))
    usage();
  auto add_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<ray::RGBA> &colours)
  {
    trajectory.calculateStartPoints(times, starts);
    ray::writePlyChunk(ofs, buffer, starts, ends, times, colours);
  };
  if (name_end == ".ply")
  {
    if (!ray::readPly(point_cloud, false, add_chunk)) // special case of reading a non-ray-cloud ply
      usage();
  }
  else if (name_end == ".laz" || name_end == ".las")
  {
    if (!ray::readLas(point_cloud, add_chunk))
      usage();
  }
  else
  {
    std::cout << "Error converting unknown type: " << point_cloud << std::endl;
    usage();
  }
  ray::writePlyChunkEnd(ofs);
  return 0;
}
