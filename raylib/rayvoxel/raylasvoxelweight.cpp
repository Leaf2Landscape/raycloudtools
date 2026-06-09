// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// Implements per-echo weighting strategies for rayvoxel beam processing.

#include "raylib/rayvoxel/raylasvoxelweight.h"
#include "raylib/rayvoxel/raylasvoxelprocessor.h"  // For PointData

#include <stdexcept>

namespace ray
{

WeightMethod parseWeightMethod(const std::string& s)
{
  if (s == "equal") return WeightMethod::kEqual;
  if (s == "full") return WeightMethod::kFull;
  if (s == "first") return WeightMethod::kFirst;
  if (s == "relative") return WeightMethod::kRelative;
  if (s == "strongest") return WeightMethod::kStrongest;
  throw std::invalid_argument("unknown weighting_method '" + s +
                              "'. Must be equal, full, first, relative, or strongest.");
}

static inline int effectiveNumReturns(uint8_t n) { return n <= 1 ? 1 : static_cast<int>(n); }

void computeEchoWeights(WeightMethod method,
                        const PointData* const* sorted, int N,
                        float* echo_w)
{
  if (N <= 0) return;

  switch (method)
  {
    case WeightMethod::kEqual:
    {
      for (int k = 0; k < N; ++k)
        echo_w[k] = 1.0f / static_cast<float>(effectiveNumReturns(sorted[k]->number_of_returns));
      break;
    }
    case WeightMethod::kFull:
    {
      for (int k = 0; k < N; ++k) echo_w[k] = 0.0f;
      echo_w[N - 1] = 1.0f;  // only the last echo carries the full pulse
      break;
    }
    case WeightMethod::kFirst:
    {
      for (int k = 0; k < N; ++k) echo_w[k] = 0.0f;
      echo_w[0] = (sorted[0]->bound == 1) ? 1.0f : 0.0f;
      break;
    }
    case WeightMethod::kRelative:
    {
      float sum = 0.0f;
      for (int k = 0; k < N; ++k) sum += static_cast<float>(sorted[k]->intensity);
      if (sum == 0.0f)
      {
        // Fallback to kEqual when the intensity sum is degenerate.
        for (int k = 0; k < N; ++k)
          echo_w[k] = 1.0f / static_cast<float>(effectiveNumReturns(sorted[k]->number_of_returns));
      }
      else
      {
        for (int k = 0; k < N; ++k)
          echo_w[k] = static_cast<float>(sorted[k]->intensity) / sum;
      }
      break;
    }
    case WeightMethod::kStrongest:
    {
      // First occurrence of the maximum intensity. Since sorted[] is ordered by
      // return_number, the first occurrence is the lowest return_number tie-break.
      int best = 0;
      for (int k = 1; k < N; ++k)
      {
        if (sorted[k]->intensity > sorted[best]->intensity) best = k;
      }
      for (int k = 0; k < N; ++k) echo_w[k] = 0.0f;
      echo_w[best] = 1.0f;
      break;
    }
  }
}

void computeSegmentWeights(const float* echo_w, int N, float* seg_w)
{
  if (N == 0) return;
  seg_w[N - 1] = echo_w[N - 1];
  for (int k = N - 2; k >= 0; --k)
    seg_w[k] = seg_w[k + 1] + echo_w[k];
}

}  // namespace ray
