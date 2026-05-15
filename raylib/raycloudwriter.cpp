// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloudwriter.h"
#include "raycloud.h"
#include "rayparse.h"

namespace ray
{
bool CloudWriter::begin(const std::string &file_name)
{
  if (file_name.empty())
  {
    std::cerr << "Error: cloud writer begin called with empty file name" << std::endl;
    return false;
  }
  file_name_ = file_name;
  const std::string ext = getFileNameExtension(file_name);
  use_las_ = (ext == "las" || ext == "laz");

  if (use_las_)
  {
    las_writer_ = new LasRayCloudWriter(file_name_);
    return true;
  }

  has_warned_ = false;
  return writeRayCloudChunkStart(file_name_, ofs_);
}

void CloudWriter::end()
{
  if (file_name_.empty())
    return;

  if (use_las_)
  {
    if (las_writer_)
    {
      std::cout << las_writer_->pointCount() << " rays saved to " << file_name_ << std::endl;
      delete las_writer_;
      las_writer_ = nullptr;
    }
  }
  else
  {
    const unsigned long num_rays = writeRayCloudChunkEnd(ofs_);
    std::cout << num_rays << " rays saved to " << file_name_ << std::endl;
    ofs_.close();
  }
  file_name_.clear();
}

bool CloudWriter::writeChunk(const Cloud &chunk)
{
  if (use_las_ && las_writer_)
    return las_writer_->writeChunk(chunk.starts, chunk.ends, chunk.times, chunk.colours,
                                   chunk.tree_ids, chunk.passthrough);
  return writeChunk(const_cast<std::vector<Eigen::Vector3d> &>(chunk.starts),
                    const_cast<std::vector<Eigen::Vector3d> &>(chunk.ends),
                    const_cast<std::vector<double> &>(chunk.times),
                    const_cast<std::vector<RGBA> &>(chunk.colours));
}

bool CloudWriter::writeChunk(std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                             std::vector<double> &times, std::vector<RGBA> &colours)
{
  if (use_las_)
  {
    if (!las_writer_)
      return false;
    return las_writer_->writeChunk(starts, ends, times, colours);
  }
  return writeRayCloudChunk(ofs_, buffer_, starts, ends, times, colours, has_warned_);
}

}  // namespace ray
