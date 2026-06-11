// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// Relaxed atomic helpers for float and uint64 — used when worker threads write
// directly into the shared flat voxel array without a per-thread map.
// Uses GCC/Clang __atomic builtins (C++17-safe, no std::atomic_ref required).

#ifndef RAYLIB_RAYVOXEL_RAYLASATOMIC_H
#define RAYLIB_RAYVOXEL_RAYLASATOMIC_H

#include <cstdint>
#include <cstring>

namespace ray
{
  // Relaxed atomic float add via compare-exchange loop.
  // On x86-64 there is no hardware lock-fadd; this CAS loop is the standard approach.
  static inline void atomic_fadd(float& dest, float delta) noexcept
  {
    static_assert(sizeof(float) == sizeof(uint32_t), "float must be 32-bit");
    uint32_t* addr = reinterpret_cast<uint32_t*>(&dest);
    uint32_t expected = __atomic_load_n(addr, __ATOMIC_RELAXED);
    uint32_t desired;
    do {
      float fval;
      std::memcpy(&fval, &expected, 4);
      fval += delta;
      std::memcpy(&desired, &fval, 4);
    } while (!__atomic_compare_exchange_n(addr, &expected, desired,
                                          /*weak=*/true,
                                          __ATOMIC_RELAXED, __ATOMIC_RELAXED));
  }

  // Relaxed atomic integer add — maps to a single lock-add instruction on x86-64.
  static inline void atomic_iadd(int32_t& dest, int32_t delta) noexcept
  {
    __atomic_fetch_add(&dest, delta, __ATOMIC_RELAXED);
  }

  // Relaxed atomic bitwise-OR for uint64 — maps to a single lock-or instruction.
  static inline void atomic_or_u64(uint64_t& dest, uint64_t bits) noexcept
  {
    __atomic_fetch_or(&dest, bits, __ATOMIC_RELAXED);
  }

  // Saturating atomic increment for uint8 — CAS loop since there is no hardware lock-incb.
  static inline void atomic_inc_u8_sat(uint8_t& dest) noexcept
  {
    uint8_t expected = __atomic_load_n(&dest, __ATOMIC_RELAXED);
    while (expected < 255) {
      uint8_t desired = expected + 1;
      if (__atomic_compare_exchange_n(&dest, &expected, desired,
                                      /*weak=*/true,
                                      __ATOMIC_RELAXED, __ATOMIC_RELAXED))
        break;
    }
  }

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASATOMIC_H
