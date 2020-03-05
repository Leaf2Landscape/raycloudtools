#pragma once
#include "rayutils.h"
#include "raytreegen.h"

namespace RAY
{
struct ForestGen
{
  void make(double randomFactor = 0.0);
  void generateRays(double rayDensity);
  std::vector<Eigen::Vector3d> getCanopy();
  std::vector<Eigen::Vector3d> getPointCloud();

  std::vector<TreeGen> trees;
};
}