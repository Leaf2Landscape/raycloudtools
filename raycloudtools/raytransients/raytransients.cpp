// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raycloud.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/raydebugdraw.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace ray;

void usage(int exit_code = 0)
{
  cout << "Splits a raycloud into the transient rays and the fixed part" << endl;
  cout << "usage:" << endl;
  cout << "raytransients min raycloud 20 rays - splits out positive transients (objects that have since moved)."
       << endl;
  cout << "                                     20 is number of pass through rays to classify as transient." << endl;
  cout << "              max    - finds negative transients, such as a hallway exposed when a door opens." << endl;
  cout << "              oldest - keeps the oldest geometry when there is a difference over time." << endl;
  cout << "              newest - uses the newest geometry when there is a difference over time." << endl;
  cout << " --colour     - also colours the clouds, to help tweak numRays. red: opacity, green: pass throughs, blue: "
          "planarity."
       << endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  DebugDraw::init(argc, argv, "raytransients");
  if (argc != 5 && argc != 6)
    usage();

  bool colour = false;
  if (argc == 6)
  {
    if (string(argv[5]) != "--colour" && string(argv[5]) != "-c")
      usage();
    colour = true;
  }
  double num_rays = stod(argv[3]);
  string merge_type = argv[1];
  if (merge_type != "min" && merge_type != "max" && merge_type != "oldest" && merge_type != "newest")
    usage();
  string file = argv[2];
  Cloud cloud;
  cloud.load(file);

  Cloud transient;
  Cloud fixed;
  cloud.findTransients(transient, fixed, merge_type, num_rays, colour);

  string file_stub = file;
  if (file.substr(file.length() - 4) == ".ply")
    file_stub = file.substr(0, file.length() - 4);

  transient.save(file_stub + "_transient.ply");
  fixed.save(file_stub + "_fixed.ply");
  return true;
}
