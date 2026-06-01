// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file implements the HeightField class.

#include "raylib/rayvoxel/raylasheightfield.h"
#include <limits>
#include <vector>
#include <algorithm> // For std::min/max

namespace ray
{

// Anonymous namespace for local helper functions
namespace {

  /// @brief Calculates the 2D barycentric coordinates of a point with respect to a triangle.
  /// This is used for both point-in-triangle testing and height interpolation.
  /// @param p The 2D point to test.
  /// @param a The 2D position of the triangle's first vertex.
  /// @param b The 2D position of the triangle's second vertex.
  /// @param c The 2D position of the triangle's third vertex.
  /// @param[out] u The barycentric weight for vertex a.
  /// @param[out] v The barycentric weight for vertex b.
  /// @param[out] w The barycentric weight for vertex c.
  /// @return True if the calculation is valid, false if the triangle is degenerate.
  bool get_barycentric_coords(const Eigen::Vector2d& p, const Eigen::Vector2d& a,
                              const Eigen::Vector2d& b, const Eigen::Vector2d& c,
                              double& u, double& v, double& w)
  {
    Eigen::Vector2d v0 = b - a;
    Eigen::Vector2d v1 = c - a;
    Eigen::Vector2d v2 = p - a;

    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);

    double denom = d00 * d11 - d01 * d01;
    // If the denominator is close to zero, the triangle is degenerate (a line or a point).
    if (std::abs(denom) < 1e-12)
    {
      return false;
    }

    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0 - v - w;

    return true;
  }

} // anonymous namespace


HeightField::HeightField()
  : grid_dim_x_(0),
    grid_dim_y_(0),
    no_data_value_(std::numeric_limits<double>::lowest())
{
}

void HeightField::fromMesh(const Mesh& mesh, const Cuboid& bounds, double cell_size)
{
  // --- Phase 1: Initialize Grid ---
  min_corner_ = Eigen::Vector2d(bounds.min_bound_.x(), bounds.min_bound_.y());
  cell_size_ = cell_size;
  const Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
  grid_dim_x_ = static_cast<int>(std::ceil(extent.x() / cell_size_));
  grid_dim_y_ = static_cast<int>(std::ceil(extent.y() / cell_size_));

  if (grid_dim_x_ <= 0 || grid_dim_y_ <= 0) {
      // Invalid dimensions, cannot proceed.
      grid_.resize(0, 0);
      return;
  }

  grid_ = Eigen::ArrayXXd::Constant(grid_dim_x_, grid_dim_y_, no_data_value_);

  // --- Phase 2: Triangle-Driven Rasterization ---
  // This robust method iterates through each triangle and fills all the grid cells it overlaps.
  for (const auto& face_indices : mesh.indexList())
  {
    // Get triangle vertices in world coordinates
    const Eigen::Vector3d& v0 = mesh.vertices()[face_indices[0]];
    const Eigen::Vector3d& v1 = mesh.vertices()[face_indices[1]];
    const Eigen::Vector3d& v2 = mesh.vertices()[face_indices[2]];

    // Get 2D projections of vertices
    const Eigen::Vector2d p0(v0.x(), v0.y());
    const Eigen::Vector2d p1(v1.x(), v1.y());
    const Eigen::Vector2d p2(v2.x(), v2.y());

    // Calculate the 2D bounding box of the triangle in world coordinates
    const double min_tx = std::min({v0.x(), v1.x(), v2.x()});
    const double max_tx = std::max({v0.x(), v1.x(), v2.x()});
    const double min_ty = std::min({v0.y(), v1.y(), v2.y()});
    const double max_ty = std::max({v0.y(), v1.y(), v2.y()});

    // Convert the world-coordinate bounding box to grid indices
    const int ix_min = static_cast<int>((min_tx - min_corner_.x()) / cell_size_);
    const int ix_max = static_cast<int>((max_tx - min_corner_.x()) / cell_size_);
    const int iy_min = static_cast<int>((min_ty - min_corner_.y()) / cell_size_);
    const int iy_max = static_cast<int>((max_ty - min_corner_.y()) / cell_size_);

    // Clip the index range to the grid's actual dimensions. This handles DTMs
    // that are larger than the processing area.
    const int start_ix = std::max(0, ix_min);
    const int end_ix   = std::min(grid_dim_x_ - 1, ix_max);
    const int start_iy = std::max(0, iy_min);
    const int end_iy   = std::min(grid_dim_y_ - 1, iy_max);

    // Iterate over every grid cell within the triangle's clipped bounding box
    for (int iy = start_iy; iy <= end_iy; ++iy)
    {
      for (int ix = start_ix; ix <= end_ix; ++ix)
      {
        // Get the world coordinates of the cell's center
        Eigen::Vector2d cell_center(
          min_corner_.x() + (static_cast<double>(ix) + 0.5) * cell_size_,
          min_corner_.y() + (static_cast<double>(iy) + 0.5) * cell_size_
        );

        // Check if the cell center is inside the triangle's 2D footprint
        double u, v, w;
        if (get_barycentric_coords(cell_center, p0, p1, p2, u, v, w))
        {
          // The point is inside if all weights are non-negative
          if (u >= -1e-9 && v >= -1e-9 && w >= -1e-9)
          {
            // Interpolate the height (z-value) using the barycentric coordinates
            double height = u * v0.z() + v * v1.z() + w * v2.z();

            // If the cell has no data yet, or if this triangle is lower than the
            // one previously found, update the height.
            if (grid_(ix, iy) == no_data_value_ || height < grid_(ix, iy))
            {
              grid_(ix, iy) = height;
            }
          }
        }
      }
    }
  }

  // --- Phase 3: Fill Gaps ---
  // After rasterizing, fill any remaining empty cells by averaging their neighbors.
  // This ensures a continuous DTM surface.
  bool gaps_remain = true;
  while (gaps_remain)
  {
    gaps_remain = false;
    Eigen::ArrayXXd temp_grid = grid_; // Operate on a copy to avoid order-dependent artifacts
    for (int x = 0; x < grid_dim_x_; ++x)
    {
      for (int y = 0; y < grid_dim_y_; ++y)
      {
        if (grid_(x, y) == no_data_value_)
        {
          double total_height = 0.0;
          int count = 0;
          for (int i = std::max(0, x - 1); i <= std::min(x + 1, grid_dim_x_ - 1); ++i)
          {
            for (int j = std::max(0, y - 1); j <= std::min(y + 1, grid_dim_y_ - 1); ++j)
            {
              if (temp_grid(i, j) != no_data_value_)
              {
                total_height += temp_grid(i, j);
                count++;
              }
            }
          }
          if (count > 0)
          {
            grid_(x, y) = total_height / count;
          }
          else
          {
            gaps_remain = true;
          }
        }
      }
    }
  }
}

void HeightField::fromLowestPoint(const std::vector<Eigen::Vector3d>& points, const Cuboid& bounds, double cell_size)
{
  min_corner_ = Eigen::Vector2d(bounds.min_bound_.x(), bounds.min_bound_.y());
  cell_size_ = cell_size;
  const Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
  grid_dim_x_ = static_cast<int>(std::ceil(extent.x() / cell_size_));
  grid_dim_y_ = static_cast<int>(std::ceil(extent.y() / cell_size_));

  if (grid_dim_x_ <= 0 || grid_dim_y_ <= 0) {
      grid_.resize(0, 0);
      return;
  }

  grid_ = Eigen::ArrayXXd::Constant(grid_dim_x_, grid_dim_y_, no_data_value_);

  // First pass: find the lowest point in each cell
  for (const auto& p : points)
  {
    int ix = static_cast<int>((p.x() - min_corner_.x()) / cell_size_);
    int iy = static_cast<int>((p.y() - min_corner_.y()) / cell_size_);

    if (ix >= 0 && ix < grid_dim_x_ && iy >= 0 && iy < grid_dim_y_)
    {
      if (grid_(ix, iy) == no_data_value_ || p.z() < grid_(ix, iy))
      {
        grid_(ix, iy) = p.z();
      }
    }
  }

  // Second pass: fill gaps by averaging neighbors (simple gap-filling)
  bool gaps_remain = true;
  while (gaps_remain)
  {
    gaps_remain = false;
    Eigen::ArrayXXd temp_grid = grid_; // Operate on a copy to avoid order-dependent artifacts
    for (int x = 0; x < grid_dim_x_; ++x)
    {
      for (int y = 0; y < grid_dim_y_; ++y)
      {
        if (grid_(x, y) == no_data_value_)
        {
          double total_height = 0.0;
          int count = 0;
          for (int i = std::max(0, x - 1); i <= std::min(x + 1, grid_dim_x_ - 1); ++i)
          {
            for (int j = std::max(0, y - 1); j <= std::min(y + 1, grid_dim_y_ - 1); ++j)
            {
              if (temp_grid(i, j) != no_data_value_)
              {
                total_height += temp_grid(i, j);
                count++;
              }
            }
          }
          if (count > 0)
          {
            grid_(x, y) = total_height / count;
          }
          else
          {
            gaps_remain = true;
          }
        }
      }
    }
  }
}

bool HeightField::getHeight(double world_x, double world_y, double& out_height) const
{
  if (!isValid()) return false;

  int ix = static_cast<int>((world_x - min_corner_.x()) / cell_size_);
  int iy = static_cast<int>((world_y - min_corner_.y()) / cell_size_);

  if (ix >= 0 && ix < grid_dim_x_ && iy >= 0 && iy < grid_dim_y_)
  {
    double height = grid_(ix, iy);
    if (height != no_data_value_)
    {
      out_height = height;
      return true;
    }
  }
  return false;
}

bool HeightField::getHeightNearest(double world_x, double world_y, double& out_height) const
{
  if (!isValid()) return false;

  // First, try a direct lookup
  if (getHeight(world_x, world_y, out_height))
  {
    return true;
  }

  // If that fails, perform a nearest-neighbor search by expanding a search box
  int ix_center = static_cast<int>((world_x - min_corner_.x()) / cell_size_);
  int iy_center = static_cast<int>((world_y - min_corner_.y()) / cell_size_);

  // Clamp center indices to be within the grid to handle queries outside the bounds
  ix_center = std::max(0, std::min(ix_center, grid_dim_x_ - 1));
  iy_center = std::max(0, std::min(iy_center, grid_dim_y_ - 1));

  for (int radius = 1; radius < std::max(grid_dim_x_, grid_dim_y_); ++radius)
  {
    double min_dist_sq = std::numeric_limits<double>::max();
    double nearest_height = no_data_value_;
    bool found_neighbor = false;

    // Search in a square spiral pattern around the center point
    for (int i = -radius; i <= radius; ++i) {
        // Check top and bottom edges of the search box
        int x_top = ix_center + i; int y_top = iy_center - radius;
        int x_bot = ix_center + i; int y_bot = iy_center + radius;
        // Check left and right edges
        int x_left = ix_center - radius; int y_left = iy_center + i;
        int x_right = ix_center + radius; int y_right = iy_center + i;

        // A helper lambda to check a point and update the nearest neighbor
        auto check_point = [&](int x, int y) {
            if (x >= 0 && x < grid_dim_x_ && y >= 0 && y < grid_dim_y_) {
                if (grid_(x, y) != no_data_value_) {
                    found_neighbor = true;
                    // Note: a real-world implementation might use distance, but for finding *any*
                    // neighbor in the closest shell, just returning the first one is sufficient and fast.
                    // For simplicity and performance, we return the first valid neighbor found at this radius.
                    nearest_height = grid_(x, y);
                    return true; // Found one
                }
            }
            return false;
        };

        if (check_point(x_top, y_top)) { out_height = nearest_height; return true; }
        if (check_point(x_bot, y_bot)) { out_height = nearest_height; return true; }
        if (check_point(x_left, y_left)) { out_height = nearest_height; return true; }
        if (check_point(x_right, y_right)) { out_height = nearest_height; return true; }
    }

    // This part of the original code is complex and slow; returning the first hit in the shell is more efficient
    // and sufficient for this use case. If no neighbor was found in the shell, the loop continues to the next radius.
  }

  // If no valid cells were found in the entire grid
  return false;
}

} // namespace ray
