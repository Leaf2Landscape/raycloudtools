// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYOCCUPANCY2D_H
#define RAYLIB_RAYOCCUPANCY2D_H

#include "raylib/raylibconfig.h"
#include "raytrunks.h"
#include "../rayutils.h"

namespace ray
{
// 2D occupancy grid used in estimating free space in a ray cloud
struct Occupancy2D
{
  /// initialise for a given bounds and pixel width
  void init(const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound, double pixel_width);

  /// save the grid
  void save(const std::string &filename);

  /// load the grid
  bool load(const std::string &filename);

  static constexpr int subpixels = 4;
  /// pixel structure that is used to estimate the density of free space
  struct Pixel
  {
    uint16_t bits;
    inline double density() const { return 1.0 - ((double)bits / (double)(subpixels * subpixels)); }
  };

  Eigen::Vector3i dims_;
  Eigen::Vector3d min_bound_;
  double pixel_width_;
  std::vector<Pixel> pixels_;
  Pixel dummy_pixel_;

  /// 3D index of the pixel
  inline Eigen::Vector3i pixelIndex(const Eigen::Vector3d &pos) const
  {
    return ((pos - min_bound_) / pixel_width_).cast<int>();
  }
  /// pixel accessors
  inline const Pixel &pixel(const Eigen::Vector3i &index) const
  {
    if (index[0] < 0 || index[1] < 0 || index[0] >= dims_[0] || index[1] >= dims_[1])
      return dummy_pixel_;
    return pixels_[dims_[1] * index[0] + index[1]];
  }
  inline Pixel &pixel(const Eigen::Vector3i &index)
  {
    if (index[0] < 0 || index[1] < 0 || index[0] >= dims_[0] || index[1] >= dims_[1])
      return dummy_pixel_;
    return pixels_[dims_[1] * index[0] + index[1]];
  }
  inline const Pixel &pixel(const Eigen::Vector3d &pos) const { return pixel(pixelIndex(pos)); }
  inline Pixel &pixel(const Eigen::Vector3d &pos) { return pixel(pixelIndex(pos)); }

  /// fill in the occupancy data (pixels_) based on the rays in the cloud @c cloudname, 
  /// within a height window @c clip_min to @c clip_max
  void fillDensities(const std::string &cloudname, const Eigen::ArrayXXd &lows, double clip_min, double clip_max);

  /// draw the occupancy grid
  void draw(const std::string &filename);
};


}  // namespace ray

#endif  // RAYLIB_RAYOCCUPANCY2D_H
