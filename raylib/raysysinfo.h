// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
#ifndef RAYLIB_RAYSYSINFO_H
#define RAYLIB_RAYSYSINFO_H

#include "raylib/raylibconfig.h"

#include <cstddef>

namespace ray
{
// Query the available memory budget in bytes. Honours cgroup limits (v2 then v1),
// SLURM allocation env vars, and /proc/meminfo:MemAvailable, in that priority order.
// Falls back to 512 MB when no source is available.
size_t RAYLIB_EXPORT queryAvailableMemoryBytes();

// Return the number of threads available for parallel work. Priority order:
// OMP_NUM_THREADS env var, SLURM_CPUS_PER_TASK env var, hardware_concurrency().
// Use this instead of hardware_concurrency() in container/HPC environments where
// the visible CPU count exceeds the allocation (e.g. Apptainer on SLURM).
size_t RAYLIB_EXPORT computeAvailableThreads();

// Compute a readLas chunk size scaled to the thread count (0 = computeAvailableThreads).
// Clamped to a floor of 4M and a ceiling of 16M points.
size_t RAYLIB_EXPORT computeReadChunkSize(size_t num_threads = 0);
}  // namespace ray

#endif  // RAYLIB_RAYSYSINFO_H
