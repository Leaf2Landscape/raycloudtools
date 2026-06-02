// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
#include "raylib/raysysinfo.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <thread>

namespace ray
{
size_t queryAvailableMemoryBytes()
{
  // cgroups v2: memory.max holds a byte count, or "max" when unlimited.
  {
    std::ifstream f("/sys/fs/cgroup/memory.max");
    std::string value;
    if (f >> value && value != "max")
    {
      size_t bytes = 0;
      try
      {
        bytes = std::stoull(value);
      }
      catch (...)
      {
        bytes = 0;
      }
      if (bytes > 0)
        return bytes;
    }
  }

  // cgroups v1: memory.limit_in_bytes uses a very large sentinel when unlimited.
  {
    std::ifstream f("/sys/fs/cgroup/memory/memory.limit_in_bytes");
    size_t bytes = 0;
    if (f >> bytes && bytes != std::numeric_limits<size_t>::max() && bytes < (size_t(1) << 60))
      return bytes;
  }

  // SLURM: per-node memory in MB.
  if (const char *mem_per_node = std::getenv("SLURM_MEM_PER_NODE"))
  {
    try
    {
      const size_t mb = std::stoull(mem_per_node);
      if (mb > 0)
        return mb * 1024ULL * 1024ULL;
    }
    catch (...)
    {
    }
  }

  // SLURM: per-CPU memory (MB) times CPUs per task.
  if (const char *mem_per_cpu = std::getenv("SLURM_MEM_PER_CPU"))
  {
    if (const char *cpus_per_task = std::getenv("SLURM_CPUS_PER_TASK"))
    {
      try
      {
        const size_t mb = std::stoull(mem_per_cpu);
        const size_t cpus = std::stoull(cpus_per_task);
        if (mb > 0 && cpus > 0)
          return mb * cpus * 1024ULL * 1024ULL;
      }
      catch (...)
      {
      }
    }
  }

  // /proc/meminfo: MemAvailable is reported in kB.
  {
    std::ifstream f("/proc/meminfo");
    std::string line;
    while (std::getline(f, line))
    {
      if (line.find("MemAvailable:") == 0)
      {
        size_t kb = 0;
        try
        {
          kb = std::stoull(line.substr(std::string("MemAvailable:").size()));
        }
        catch (...)
        {
          kb = 0;
        }
        if (kb > 0)
          return kb * 1024ULL;
      }
    }
  }

  return 512ULL * 1024 * 1024;
}

size_t computeReadChunkSize(size_t num_threads)
{
  const size_t resolved = (num_threads == 0)
      ? static_cast<size_t>(std::thread::hardware_concurrency())
      : num_threads;
  // kBeamBatchSize = 32 (matches raylasvoxelise.cpp constant)
  constexpr size_t kBatchMultiplier = 32 * 4;
  const size_t scaled = resolved * kBatchMultiplier;
  const size_t clamped = std::max(size_t(4000000), std::min(scaled, size_t(16000000)));
  return clamped;
}
}  // namespace ray
