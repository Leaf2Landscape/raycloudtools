// Per-echo weighting strategies for rayvoxel beam processing.
#ifndef RAYLIB_RAYVOXEL_RAYLASVOXELWEIGHT_H
#define RAYLIB_RAYVOXEL_RAYLASVOXELWEIGHT_H

#include <string>

namespace ray {

struct PointData;  // forward declaration

enum class WeightMethod { kEqual, kFull, kFirst, kRelative, kStrongest };

// Parses CLI string → WeightMethod. Throws std::invalid_argument on unknown values.
WeightMethod parseWeightMethod(const std::string& s);

// Fills echo_w[0..N-1]. sorted[] is return_number-sorted view of beam returns.
void computeEchoWeights(WeightMethod method,
                        const PointData* const* sorted, int N,
                        float* echo_w);

// Fills seg_w[0..N-1] as suffix sums of echo_w.
// seg_w[k] = sum(echo_w[k..N-1]) = remaining beam fraction entering segment k.
void computeSegmentWeights(const float* echo_w, int N, float* seg_w);

}  // namespace ray

#endif
