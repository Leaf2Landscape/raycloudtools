// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
#include "raylib/rayvoxel/rayvox.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <limits>

namespace ray
{

  // Helper function to trim whitespace from both ends of a string
  static inline void trim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
      return !std::isspace(ch);
    }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
      return !std::isspace(ch);
    }).base(), s.end());
  }

  bool writeVox(const std::string& file_name, const VoxelSpace& space) {
    std::ofstream ofs(file_name);
    if (!ofs.is_open()) {
      std::cerr << "Error: cannot open " << file_name << " for writing." << std::endl;
      return false;
    }

    if (!writeVoxChunkStart(ofs, space.header)) {
      return false;
    }

    if (!writeVoxChunk(ofs, space.voxels)) {
      return false;
    }

    // writeVoxChunkEnd does not close the stream; the destructor of 'ofs' will handle it here.
    writeVoxChunkEnd(ofs);
    std::cout << "Saved voxel file " << file_name << " [OK]" << std::endl;
    return true;
  }

  bool writeVoxChunkStart(std::ofstream& out, const std::map<std::string, std::string>& header) {
    out << "VOXEL SPACE" << "\n";

    // Write header properties. std::map automatically keeps keys sorted.
    for (const auto& pair : header) {
      if (pair.first != "colnames" && pair.first != "nline") {
        out << "#" << pair.first << ":" << pair.second << "\n";
      }
    }

    auto it = header.find("colnames");
    if (it == header.end()) {
      std::cerr << "Error: Header must contain 'colnames' key." << std::endl;
      return false;
    }
    out << it->second << "\n";

    // Set precision for floating point numbers to avoid scientific notation, like scipen=999
    out << std::fixed << std::setprecision(10);

    return out.good();
  }

  bool writeVoxChunk(std::ofstream& out, const std::vector<VoxelData>& voxels) {
    if (voxels.empty()) {
      return true; // Not an error to write an empty chunk
    }

    for (const auto& voxel : voxels) {
      out << voxel.i << " " << voxel.j << " " << voxel.k;
      for (const auto& var : voxel.variables) {
        out << " " << var;
      }
      out << "\n";
    }

    if (!out.good()) {
      std::cerr << "Error writing voxel data to file." << std::endl;
      return false;
    }
    return true;
  }

  void writeVoxChunkEnd(std::ofstream& out) {
    // FIXED: Only flush the stream. Do not close it.
    // The caller who opened the stream is responsible for closing it.
    out.flush();
  }

  bool readVox(const std::string& file_name, VoxelSpace& space_out) {
    space_out.clear();

    auto apply = [&](const std::vector<VoxelData>& chunk) {
      space_out.voxels.insert(space_out.voxels.end(), chunk.begin(), chunk.end());
    };

    // FIXED: Do not pass a max chunk_size, which causes a bad_alloc on reserve.
    // Rely on the default chunk size for efficient reading.
    return readVox(file_name, space_out.header, apply);
  }

  bool readVox(const std::string& file_name,
               std::map<std::string, std::string>& header_out,
               std::function<void(const std::vector<VoxelData>&)> apply,
               size_t chunk_size) {

    std::ifstream ifs(file_name);
    if (!ifs.is_open()) {
      std::cerr << "Error: cannot open " << file_name << " for reading." << std::endl;
      return false;
    }

    std::string line;
    if (!std::getline(ifs, line)) {
      std::cerr << "Error: File is empty or could not be read." << std::endl;
      return false;
    }
    trim(line);
    if (line != "VOXEL SPACE") {
      std::cerr << "Error: Invalid voxel file. Expected 'VOXEL SPACE' identifier." << std::endl;
      return false;
    }

    header_out.clear();
    int nLineHeader = 0;
    while (std::getline(ifs, line)) {
      nLineHeader++;
      trim(line);
      if (line.empty()) continue;

      if (line[0] != '#') {
        header_out["colnames"] = line;
        break;
      }

      // Remove leading '#' characters before processing
      size_t first_char = line.find_first_not_of('#');
      if (first_char == std::string::npos) continue;

      std::string segment = line.substr(first_char);
      size_t colon_pos = segment.find(':');
      if (colon_pos == std::string::npos) {
          std::cerr << "Warning: Malformed header line, skipping: " << segment << std::endl;
          continue;
      }

      std::string key = segment.substr(0, colon_pos);
      std::string value = segment.substr(colon_pos + 1);
      trim(key);
      trim(value);
      header_out[key] = value;
    }

    if (header_out.find("colnames") == header_out.end()) {
      std::cerr << "Error: Voxel file is missing column names header." << std::endl;
      return false;
    }
    header_out["nline"] = std::to_string(nLineHeader);

    std::vector<VoxelData> chunk;
    chunk.reserve(chunk_size);

    while (std::getline(ifs, line)) {
      trim(line);
      if(line.empty()) continue;

      std::stringstream ss(line);
      VoxelData voxel;

      if (!(ss >> voxel.i >> voxel.j >> voxel.k)) {
        std::cerr << "Warning: Could not parse voxel indices from line: " << line << std::endl;
        continue;
      }

      std::string var;
      while (ss >> var) {
        voxel.variables.push_back(var);
      }
      chunk.push_back(voxel);

      if (chunk.size() >= chunk_size) {
        apply(chunk);
        chunk.clear();
      }
    }

    if (!chunk.empty()) {
      apply(chunk);
    }

    return true;
  }

} // namespace ray
