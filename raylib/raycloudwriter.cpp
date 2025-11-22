// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raycloudwriter.h"
#include "raycloud.h"
#include "raylaz.h"

namespace ray
{

// --- START OF MODIFICATION ---
CloudWriter::CloudWriter()
    : has_warned_(false), writer_type_(NONE), las_writer_(nullptr), ply_vertex_byte_size_(0), ply_header_written_(false)
{
}

CloudWriter::~CloudWriter()
{
    if (writer_type_ == LAZ) {
        delete las_writer_;
        las_writer_ = nullptr;
    }
}

bool CloudWriter::begin(const std::string &file_name, const std::string &output_ext)
{
  if (file_name.empty())
  {
    std::cerr << "Error: cloud writer begin called with empty file name" << std::endl;
    return false;
  }
  has_warned_ = false;
  file_name_ = file_name;

  if (output_ext == ".las" || output_ext == ".laz")
  {
      writer_type_ = LAZ;
      las_writer_ = new LasWriter(file_name_);
  }
  else // Default to PLY
  {
      writer_type_ = PLY;
      // The header for PLY cannot be written until the first chunk,
      // as we need to know which optional fields are present.
  }
  return true;
}

void CloudWriter::end()
{
  if (file_name_.empty()) return;

  if (writer_type_ == PLY && ply_header_written_)
  {
    const unsigned long num_rays = writeRayCloudChunkEnd(ofs_, ply_vertex_byte_size_);
    std::cout << num_rays << " rays saved to " << file_name_ << std::endl;
    ofs_.close();
  }
  else if (writer_type_ == LAZ)
  {
      delete las_writer_;
      las_writer_ = nullptr;
      std::cout << "File saved to " << file_name_ << std::endl;
  }
  file_name_.clear(); // Prevent double-closing
}

bool CloudWriter::writeChunk(const Cloud &chunk)
{
    return writeChunk(chunk.starts, chunk.ends, chunk.times, chunk.colours, chunk.classifications, chunk.branch_ids);
}

bool CloudWriter::writeChunk(const std::vector<Eigen::Vector3d> &starts, const std::vector<Eigen::Vector3d> &ends, const std::vector<double> &times,
                             const std::vector<RGBA> &colours, const std::vector<uint8_t> &classifications, const std::vector<uint16_t> &branch_ids)
{
    if (ends.empty()) return true;

    if (writer_type_ == LAZ)
    {
        if (las_writer_) {
            return las_writer_->writeChunk(ends, times, colours, classifications, branch_ids);
        }
        return false;
    }
    else if (writer_type_ == PLY)
    {
        if (!ply_header_written_) {
            bool has_classification = !classifications.empty();
            bool has_branch_id = !branch_ids.empty();
            if (!writeRayCloudChunkStart(file_name_, ofs_, has_classification, has_branch_id)) {
                return false;
            }

            ply_vertex_byte_size_ = 0;
            #if RAYLIB_DOUBLE_RAYS
            ply_vertex_byte_size_ += sizeof(double) * 3;
            #else
            ply_vertex_byte_size_ += sizeof(float) * 3;
            #endif
            ply_vertex_byte_size_ += sizeof(double); // time
            ply_vertex_byte_size_ += sizeof(float) * 3; // ray
            ply_vertex_byte_size_ += sizeof(RGBA);
            if (has_classification) ply_vertex_byte_size_ += sizeof(uint8_t);
            if (has_branch_id) ply_vertex_byte_size_ += sizeof(uint16_t);
            
            ply_header_written_ = true;
        }
        return writeRayCloudChunk(ofs_, buffer_, starts, ends, times, colours, classifications, branch_ids, has_warned_);
    }

    return false; // No writer type
}
// --- END OF MODIFICATION ---

}  // namespace ray