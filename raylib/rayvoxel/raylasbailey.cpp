// Bailey, B.N. & Mahaffee, W.F. (2017). "Rapid, high-resolution measurement
// of leaf area and leaf orientation using terrestrial LiDAR scanning data."
// Measurement Science and Technology 28(6):064006.
// DOI: 10.1088/1361-6501/aa5cfd
//
#include "raylib/rayvoxel/raylasbailey.h"
#include "raylib/rayutils.h"  // kPi

#include <algorithm>
#include <cmath>
#include <set>

namespace ray
{

std::unordered_map<int64_t, TriangleHistograms> buildTriangleInclinationHistograms(
    const std::vector<Eigen::Vector3d>& positions,
    const Eigen::MatrixXi& knn_indices,
    const std::vector<int64_t>& flat_indices,
    const std::vector<int>& class_labels,
    int n_bins,
    double l_max)
{
  std::unordered_map<int64_t, TriangleHistograms> result;
  const int N = static_cast<int>(positions.size());
  const int K = static_cast<int>(knn_indices.rows());
  const double l_max2 = l_max * l_max;

  // Bailey & Mahaffee (2017) eq.(3): G_i = |r_hat · n_hat|; eq.(4): G_bar = Σ(G_i·A_i·sinθ_i)/Σ(A_i·sinθ_i)
  // r_hat is the mean ray direction for the voxel. A per-voxel mean ray direction is not available
  // inside this function, so we use a vertical approximation r_hat = (0,0,1).
  const Eigen::Vector3d r_hat(0.0, 0.0, 1.0);
  // Eq.(4) numerator/denominator accumulators, keyed by flat voxel index (finalized after the loop).
  std::unordered_map<int64_t, double> g_num_leaf, g_den_leaf, g_num_wood, g_den_wood;

  // Deduplicate facets by their sorted vertex-id triple.
  std::set<std::array<int, 3>> seen;

  for (int i = 0; i < N; ++i) {
    const int ci = class_labels[i];
    if (ci == 0) continue;  // unknown class cannot form a same-class facet

    // Candidate neighbours of i: same class, within l_max of i.
    for (int ja = 0; ja < K; ++ja) {
      const int a = knn_indices(ja, i);
      if (a < 0 || a == i) continue;
      if (class_labels[a] != ci) continue;
      if ((positions[a] - positions[i]).squaredNorm() > l_max2) continue;

      for (int jb = ja + 1; jb < K; ++jb) {
        const int b = knn_indices(jb, i);
        if (b < 0 || b == i || b == a) continue;
        if (class_labels[b] != ci) continue;
        if ((positions[b] - positions[i]).squaredNorm() > l_max2) continue;
        if ((positions[b] - positions[a]).squaredNorm() > l_max2) continue;

        // Deduplicate: represent the facet by its sorted vertex-id triple.
        std::array<int, 3> tri = { i, a, b };
        std::sort(tri.begin(), tri.end());
        if (!seen.insert(tri).second) continue;

        // Facet normal and area.
        const Eigen::Vector3d v0 = positions[i];
        const Eigen::Vector3d v1 = positions[a];
        const Eigen::Vector3d v2 = positions[b];
        const Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
        const double n_norm = n.norm();
        if (n_norm < 1e-12) continue;  // degenerate facet
        const double area = 0.5 * n_norm;

        // Inclination of the facet normal relative to vertical.
        const double theta = std::acos(std::clamp(std::abs(n.z()) / n_norm, 0.0, 1.0));
        int bin = static_cast<int>(theta / (kPi / 2.0) * n_bins);
        bin = std::clamp(bin, 0, n_bins - 1);

        // eq.(3) per-facet projection coefficient against the (vertical-approximated) ray direction.
        const double G_i = std::abs(r_hat.dot(n / n_norm));
        const double w = area * std::sin(theta);  // eq.(4) facet weight

        // Assign the facet to the voxel of its first vertex i.
        const int64_t fidx = flat_indices[i];
        TriangleHistograms& th = result[fidx];
        if (ci > 0) {
          if (th.tiad_leaf.empty()) th.tiad_leaf.assign(n_bins, 0.0);
          th.tiad_leaf[bin] += w;
          th.total_leaf_area += area;
          g_num_leaf[fidx] += G_i * w;
          g_den_leaf[fidx] += w;
        } else {
          if (th.tiad_wood.empty()) th.tiad_wood.assign(n_bins, 0.0);
          th.tiad_wood[bin] += w;
          th.total_wood_area += area;
          g_num_wood[fidx] += G_i * w;
          g_den_wood[fidx] += w;
        }
      }
    }
  }

  // Finalize eq.(4) per-class mean G for each voxel (guard denominator > 0).
  for (auto& pair : result) {
    TriangleHistograms& th = pair.second;
    auto lnit = g_num_leaf.find(pair.first);
    auto ldit = g_den_leaf.find(pair.first);
    if (ldit != g_den_leaf.end() && ldit->second > 0.0)
      th.bailey_g_leaf = lnit->second / ldit->second;
    auto wnit = g_num_wood.find(pair.first);
    auto wdit = g_den_wood.find(pair.first);
    if (wdit != g_den_wood.end() && wdit->second > 0.0)
      th.bailey_g_wood = wnit->second / wdit->second;
  }

  return result;
}

double solveBaileyPadEq10(double path_length_observed,
                          double num_beams_weighted,
                          double num_hits,
                          double G)
{
  // Guard: no usable geometry / projection.
  if (G <= 0.0 || num_beams_weighted < 1.0) return 0.0;
  const double r_bar = path_length_observed / num_beams_weighted;  // mean scan path through voxel
  if (r_bar <= 0.0) return 0.0;

  // Mean gap probability. Clamp away from 0 so -ln is finite; a fully-blocked voxel
  // would otherwise drive a_L -> infinity.
  const double eps = 1e-12;
  double P_bar = 1.0 - (num_hits / num_beams_weighted);
  P_bar = std::clamp(P_bar, eps, 1.0 - eps);

  // Bailey & Mahaffee (2017) Eq. 10 residual: f(a_L) = P_bar - exp(-a_L * G * r_bar).
  // Approximation: Eq. 10 specifies an exponentially-weighted mean path r̄_exp; only an
  // arithmetic mean (path_length_observed / num_beams_weighted) is available from the voxel
  // accumulators, so r̄_exp is replaced by r_bar. This biases a_L when path length varies
  // significantly within a voxel. Higher-accuracy accumulation is a tracked follow-up.
  auto residual = [&](double a_L) -> double {
    return P_bar - std::exp(-a_L * G * r_bar);
  };

  // Seed with the closed-form thin-medium solution.
  double a0 = -std::log(P_bar) / (G * r_bar);
  double a1 = a0 * 1.0001 + 1e-6;  // second seed for the secant method
  double f0 = residual(a0);
  double f1 = residual(a1);

  for (int iter = 0; iter < 50; ++iter) {
    if (std::abs(f1) < 1e-9) break;
    const double denom = (f1 - f0);
    if (std::abs(denom) < 1e-18) break;  // flat residual; keep current estimate
    const double a2 = a1 - f1 * (a1 - a0) / denom;
    a0 = a1; f0 = f1;
    a1 = a2; f1 = residual(a1);
  }

  return std::max(0.0, a1);
}

} // namespace ray
