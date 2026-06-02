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

// Compute a readLas chunk size scaled to the thread count (0 = hardware_concurrency).
// Clamped to a floor of 4M and a ceiling of 16M points.
size_t RAYLIB_EXPORT computeReadChunkSize(size_t num_threads = 0);
}  // namespace ray

#endif  // RAYLIB_RAYSYSINFO_H
