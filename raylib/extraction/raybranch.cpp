// Copyright (c) 2021
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raybranch.h"
#include "../raydebugdraw.h"
#include "../raygrid.h"
#include "../raycuboid.h"
#include <map>
#include <nabo/nabo.h>
#include <queue>

namespace ray
{
const double inf = 1e10;

Branch::Branch() : centre(0,0,0), radius(0), score(0), last_score(0), length(0), dir(0,0,0), parent(-1), tree_score(inf), distance_to_ground(inf), active(true), visited(false) 
{
}

void Branch::getOverlap(const Grid<Eigen::Vector3d> &grid, std::vector<Eigen::Vector3d> &points, double spacing)
{
  Eigen::Vector3d base = centre - 0.5*length*dir;
  Eigen::Vector3d top = centre + 0.5*length*dir;
  double outer_radius = (radius + spacing) * boundary_radius_scale;
  Eigen::Vector3d rad(outer_radius, outer_radius, outer_radius);
  Cuboid cuboid(minVector(base, top) - rad, maxVector(base, top) + rad);

  Eigen::Vector3i mins = ((cuboid.min_bound_ - grid.box_min) / grid.voxel_width).cast<int>();
  Eigen::Vector3i maxs = ((cuboid.max_bound_ - grid.box_min) / grid.voxel_width).cast<int>();

  mins = maxVector(mins, Eigen::Vector3i(0,0,0));
  Eigen::Vector3i min_dims = grid.dims - Eigen::Vector3i(1,1,1);
  maxs = minVector(maxs, min_dims);

  Eigen::Vector3i ind;
  for (ind[0] = mins[0]; ind[0]<=maxs[0]; ind[0]++)
  {
    for (ind[1] = mins[1]; ind[1]<=maxs[1]; ind[1]++)
    {
      for (ind[2] = mins[2]; ind[2]<=maxs[2]; ind[2]++)
      {
        auto &cell = grid.cell(ind);
        for (auto &pos: cell.data)
        {
          Eigen::Vector3d p = pos - centre;
          double h = p.dot(dir);
          if (std::abs(h) > length*0.5)
          {
            continue;
          }
          p -= dir*h;
          double dist2 = p.squaredNorm();
          if (dist2 <= outer_radius*outer_radius)
          {
            points.push_back(pos);
          }
        }
      }
    }
  }
}   

void Branch::estimatePose(const std::vector<Eigen::Vector3d> &points)
{
  centre = mean(points);
  Eigen::Matrix3d scatter;
  scatter.setZero();
  for (auto &point: points)
    scatter += (point - centre) * (point - centre).transpose();
  scatter /= (double)points.size();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(scatter.transpose());   
  dir = eigen_solver.eigenvectors().col(2);       
}

void Branch::updateDirection(const std::vector<Eigen::Vector3d> &points)
{
  struct Accumulator
  {
    Accumulator(): weight(0), x(0), abs_x(0), y(0,0), xy(0,0), x2(0) {}

    double weight;
    double x;
    double abs_x;
    Eigen::Vector2d y;
    Eigen::Vector2d xy;
    double x2;
  };
  Eigen::Vector3d ax1 = Eigen::Vector3d(1,2,3).cross(dir).normalized();
  Eigen::Vector3d ax2 = ax1.cross(dir);
  Accumulator sum;
  for (size_t i = 0; i<points.size(); i++)
  {
    Eigen::Vector3d to_point = points[i] - centre;
    Eigen::Vector2d offset(to_point.dot(ax1), to_point.dot(ax2));
    //const double dist = offset.norm();
    double w = 1.0;// 1.0 - dist/(radius * boundary_radius_scale) // lateral fade off?
    // remove radius. If radius_removal_factor=0 then half-sided trees will have estimated branch centred on that edge
    //                If radius_removal_factor=1 then v thin branches may accidentally get a radius and it won't shrink down
    const double radius_removal_factor = 0.5;
    offset -= offset * radius_removal_factor * radius / offset.norm(); 

    double h = to_point.dot(dir);
    // lean, shift and change radius
    sum.x += h*w;
    sum.y += offset*w;
    sum.xy += h*offset*w;
    sum.x2 += h*h*w;
    sum.abs_x += std::abs(h)*w;
    sum.weight += w;      
  }
  double n = sum.weight;
  centre += dir*(sum.x / n); 

  // based on http://mathworld.wolfram.com/LeastSquaresFitting.html
  Eigen::Vector2d sXY = sum.xy - sum.x*sum.y/n;
  double sXX = sum.x2 - sum.x*sum.x/n;
  if (std::abs(sXX) > 1e-10)
    sXY /= sXX;

  dir = (dir + ax1*sXY[0] + ax2*sXY[1]).normalized();
  actual_length = 4.0*(sum.abs_x/n);
  length = std::max(actual_length, 2.0*radius * branch_height_to_width);    
}

// shift to an estimation of the centre of the cylinder's circle
void Branch::updateCentre(const std::vector<Eigen::Vector3d> &points)
{
  Eigen::Vector3d ax1 = Eigen::Vector3d(1,2,3).cross(dir).normalized();
  Eigen::Vector3d ax2 = ax1.cross(dir);

  std::vector<Eigen::Vector3d> ps(points.size());
  Eigen::Vector3d mean_p(0,0,0);
  for (size_t i = 0; i<points.size(); i++)
  {
    Eigen::Vector3d to_point = points[i] - centre;
    Eigen::Vector2d offset(to_point.dot(ax1), to_point.dot(ax2));
    Eigen::Vector2d xy = offset/radius;
    double l2 = xy.squaredNorm();
    Eigen::Vector3d point(xy[0], xy[1], 0.5*l2); // a paraboloid that has gradient 1 at 1
    ps[i] = point;
    mean_p += point;
  }
  mean_p /= (double)points.size();      
  struct Acc
  {
    Acc(){ x2 = y2 = xy = xz = yz = 0; }
    double x2, y2, xy, xz, yz;
  };
  Acc plane;
  for (auto &p: ps)
  {
    Eigen::Vector3d q = p - mean_p;
    plane.x2 += q[0]*q[0];
    plane.y2 += q[1]*q[1];
    plane.xy += q[0]*q[1];        
    plane.xz += q[0]*q[2];        
    plane.yz += q[1]*q[2];        
  }      
  double A = (plane.xz*plane.y2 - plane.yz*plane.xy) / (plane.x2*plane.y2 - plane.xy*plane.xy);
  double B = (plane.yz - A * plane.xy) / plane.y2;
  Eigen::Vector2d shift(A,B);
  double shift2 = shift.squaredNorm();
  if (shift2 > 1.0) // don't shift more than one radius each iteration
    shift /= std::sqrt(shift2);

  centre += (ax1*shift[0] + ax2*shift[1]) * radius;       
}

void Branch::updateRadiusAndScore(const std::vector<Eigen::Vector3d> &points)
{
  double rad = 0;
  std::vector<double> scores(points.size());
  for (auto &point: points)
  {
    Eigen::Vector3d to_point = point - centre;
    rad += (to_point - dir*dir.dot(to_point)).norm();
  }
  double n = (double)points.size();
  radius = rad / n;
  double raddiff = 0.0;
  for (auto &point: points)
  {
    Eigen::Vector3d to_point = point - centre;
    double rad = (to_point - dir*dir.dot(to_point)).norm();
    raddiff += std::abs(rad - radius);
  }
  double num_points = (double)points.size() - 4.0; 
  double variation = raddiff / num_points; // end part gives sample variation
  
  last_score = score;
  score = actual_length / variation; // unitless
}

} // namespace ray
