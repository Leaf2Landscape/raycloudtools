// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file defines the HeightField class, which provides an efficient
// way to query terrain height from a rasterized Digital Terrain Model (DTM).
// It's used for clipping occlusion rays and calculating height-from-ground metrics.

#ifndef RAYLIB_RAYVOXEL_RAYLASHEIGHTFIELD_H
#define RAYLIB_RAYVOXEL_RAYLASHEIGHTFIELD_H

#include "raylib/raylibconfig.h"
#include "raylib/raymesh.h"
#include <Eigen/Dense>

namespace ray
{
  class RAYLIB_EXPORT HeightField
  {
  public:
    HeightField();

    /// @brief Constructs a HeightField from a mesh.
    /// @param mesh The DTM mesh to rasterize.
    /// @param bounds The world-coordinate bounds for the height field grid.
    /// @param cell_size The desired resolution (width and length) of each grid cell.
    void fromMesh(const Mesh& mesh, const Cuboid& bounds, double cell_size);

    /// @brief Constructs a HeightField by finding the lowest point within each grid cell.
    /// @param points Vector of ground points.
    /// @param bounds The world-coordinate bounds for the height field grid.
    /// @param cell_size The desired resolution of each grid cell.
    void fromLowestPoint(const std::vector<Eigen::Vector3d>& points, const Cuboid& bounds, double cell_size);

    /// @brief Retrieves the terrain height at a given world coordinate.
    /// @param world_x The x-coordinate in world space.
    /// @param world_y The y-coordinate in world space.
    /// @param out_height Reference to store the resulting height.
    /// @return True if a valid height was found at the direct location, false otherwise.
    bool getHeight(double world_x, double world_y, double& out_height) const;

    /// @brief Retrieves the terrain height, using nearest-neighbor search if the direct location has no data.
    /// @param world_x The x-coordinate in world space.
    /// @param world_y The y-coordinate in world space.
    /// @param out_height Reference to store the resulting height.
    /// @return True if a valid height was found (either direct or nearby), false if the grid is empty.
    bool getHeightNearest(double world_x, double world_y, double& out_height) const;

    /// @brief Checks if the HeightField has been successfully initialized.
    /// @return True if the grid contains data, false otherwise.
    bool isValid() const { return !grid_.isZero(); }

  private:
    Eigen::ArrayXXd grid_;
    Eigen::Vector2d min_corner_;
    double cell_size_;
    int grid_dim_x_;
    int grid_dim_y_;
    const double no_data_value_;
  };

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASHEIGHTFIELD_H
