// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloudreader.h"
#include "rayparse.h"
#include "raysysinfo.h"

namespace ray
{

bool CloudReader::begin(const std::string &file_name)
{
  file_name_ = file_name;
  const std::string ext = getFileNameExtension(file_name);
  is_las_ = (ext == "las" || ext == "laz");
  header_ = LasHeader{};
  if (is_las_)
    return readLasHeader(file_name_, header_);
  return true;
}

bool CloudReader::read(std::function<void(std::vector<Eigen::Vector3d> &,
                                          std::vector<Eigen::Vector3d> &,
                                          std::vector<double> &,
                                          std::vector<RGBA> &)>
                         apply,
                       size_t &num_bounded,
                       double max_intensity,
                       Eigen::Vector3d *offset_to_remove,
                       size_t chunk_size,
                       std::vector<int32_t> *tree_ids_out,
                       std::vector<uint8_t> *passthrough_out,
                       std::vector<int32_t> *stem_ids_out,
                       std::vector<int32_t> *beam_ids_out)
{
  const size_t csize = chunk_size > 0 ? chunk_size : computeReadChunkSize();
  if (is_las_)
    return readLas(file_name_, std::move(apply), num_bounded, max_intensity, offset_to_remove, csize,
                   tree_ids_out, passthrough_out, nullptr, nullptr, stem_ids_out, beam_ids_out);
  num_bounded = 0;
  return readPly(file_name_, true, std::move(apply), max_intensity, false, csize);
}

}  // namespace ray
