// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYPLY_H
#define RAYLIB_RAYPLY_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"

// --- START OF FIX ---
// Add missing headers for types used in this file's declarations.
#include <vector>
#include <string>
#include <functional>
#include <fstream>
// Eigen headers are included via rayutils.h, but explicitly including
// the core header here makes the dependency clear.
#include <Eigen/Core>
// --- END OF FIX ---

namespace ray
{
#if RAYLIB_DOUBLE_RAYS
using PointPlyEntry = Eigen::Matrix<float, 9, 1>;
using RayPlyEntry = Eigen::Matrix<float, 12, 1>;
#else
using PointPlyEntry = Eigen::Matrix<float, 6, 1>;
using RayPlyEntry = Eigen::Matrix<float, 9, 1>;
#endif
using PointPlyBuffer = std::vector<PointPlyEntry>;
using RayPlyBuffer = std::vector<RayPlyEntry>;

/// read in a .ply file into the fields given by reference
/// Note that @c max_intensity is only used when reading in a point cloud. Intensities are already stored in the
/// colour alpha channel in ray clouds.
// Overloaded to support new fields
bool RAYLIB_EXPORT readPly(const std::string &file_name, std::vector<Eigen::Vector3d> &starts,
                           std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<RGBA> &colours,
                           bool is_ray_cloud, double max_intensity = 0);

bool RAYLIB_EXPORT readPly(const std::string &file_name, std::vector<Eigen::Vector3d> &starts,
                           std::vector<Eigen::Vector3d> &ends, std::vector<double> &times, std::vector<RGBA> &colours,
                           std::vector<uint8_t> &classifications, std::vector<uint16_t> &branch_ids,
                           bool is_ray_cloud, double max_intensity = 0);
/// read in a .ply file that represents a triangular mesh, into the @c Mesh structure
bool RAYLIB_EXPORT readPlyMesh(const std::string &file, class Mesh &mesh);

/// write a .ply file representing a triangular mesh
bool RAYLIB_EXPORT writePlyMesh(const std::string &file_name, const class Mesh &mesh, bool flip_normals = false);

/// ready in a ray cloud or point cloud .ply file, and call the @c apply function one chunk at a time,
/// @c chunk_size is the number of rays to read at one time. This method can be used on large clouds where
/// the full set of rays is not required to be in memory at one time.
/// @c times_optional flag allows clouds to be read with no time stamps
// Updated apply function signature to be flexible
bool RAYLIB_EXPORT readPly(const std::string &file_name, bool is_ray_cloud,
                           std::function<void(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                                              std::vector<double> &times, std::vector<RGBA> &colours,
                                              std::vector<uint8_t> &classifications, std::vector<uint16_t> &branch_ids)>
                             apply, 
                           double max_intensity, bool times_optional = false, size_t chunk_size = 1000000);

/// write a .ply file representing a point cloud
// Updated signatures to handle optional extra fields
bool RAYLIB_EXPORT writePlyPointCloud(const std::string &file_name, const std::vector<Eigen::Vector3d> &points,
                                      const std::vector<double> &times, const std::vector<RGBA> &colours,
                                      const std::vector<uint8_t> &classifications,
                                      const std::vector<uint16_t> &point_source_ids,
                                      bool write_extra_fields);

/// write a .ply file representing a ray cloud
// Updated signatures to handle optional extra fields
bool RAYLIB_EXPORT writePlyRayCloud(const std::string &file_name, const std::vector<Eigen::Vector3d> &starts,
                                    const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                                    const std::vector<RGBA> &colours,
                                    const std::vector<uint8_t> &classifications,
                                    const std::vector<uint16_t> &branch_ids,
                                    bool write_extra_fields);

/// Chunked version of writePlyRayCloud
bool RAYLIB_EXPORT writeRayCloudChunkStart(const std::string &file_name, std::ofstream &out, bool has_classification, bool has_branch_id);
bool RAYLIB_EXPORT writeRayCloudChunk(std::ofstream &out, std::vector<char>& buffer,
                                      const std::vector<Eigen::Vector3d> &starts,
                                      const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                                      const std::vector<RGBA> &colours, 
                                      const std::vector<uint8_t> &classifications,
                                      const std::vector<uint16_t> &branch_ids,
                                      bool &has_warned);
unsigned long RAYLIB_EXPORT writeRayCloudChunkEnd(std::ofstream &out, size_t vertex_byte_size);

/// Chunked version of writePlyPointCloud
bool RAYLIB_EXPORT writePointCloudChunkStart(const std::string &file_name, std::ofstream &out, bool has_classification, bool has_point_source_id);
bool RAYLIB_EXPORT writePointCloudChunk(std::ofstream &out, std::vector<char>& buffer,
                                        const std::vector<Eigen::Vector3d> &points, const std::vector<double> &times,
                                        const std::vector<RGBA> &colours,
                                        const std::vector<uint8_t> &classifications,
                                        const std::vector<uint16_t> &point_source_ids,
                                        bool &has_warned);
void RAYLIB_EXPORT writePointCloudChunkEnd(std::ofstream &out, size_t vertex_byte_size);

/// Simple function for converting a ray cloud according to the per-ray function @c apply
bool convertCloud(const std::string &in_name, const std::string &out_name,
                  std::function<void(Eigen::Vector3d &start, Eigen::Vector3d &ends, double &time, RGBA &colour, uint8_t &classification, uint16_t &branch_id)> apply);
}  // namespace ray

#endif  // RAYLIB_RAYPLY_H
