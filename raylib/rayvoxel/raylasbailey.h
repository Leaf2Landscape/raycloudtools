// Bailey, B.N. & Mahaffee, W.F. (2017). "Rapid, high-resolution measurement
// of leaf area and leaf orientation using terrestrial LiDAR scanning data."
// Measurement Science and Technology 28(6):064006.
// DOI: 10.1088/1361-6501/aa5cfd
//
#ifndef RAYLIB_RAYVOXEL_RAYLASBAILEY_H
#define RAYLIB_RAYVOXEL_RAYLASBAILEY_H

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <Eigen/Dense>

namespace ray
{
  // Per-class, area-weighted inclination histograms accumulated from triangle facets
  // built out of nearby LiDAR returns (Bailey & Mahaffee 2017).
  struct TriangleHistograms {
    std::vector<double> tiad_leaf;   // area-weighted inclination histogram, size n_bins
    std::vector<double> tiad_wood;
    double total_leaf_area = 0.0;
    double total_wood_area = 0.0;
    double bailey_g_leaf = 0.0;   // Eq.(4) area*sin(theta) weighted mean G for leaf facets
    double bailey_g_wood = 0.0;   // Eq.(4) area*sin(theta) weighted mean G for wood facets
  };

  // For each voxel, build triangle facets from within-class triples of nearby LiDAR
  // returns and accumulate area * sin(theta) into per-class inclination histograms.
  // Returns a map from flat voxel index to TriangleHistograms.
  std::unordered_map<int64_t, TriangleHistograms> buildTriangleInclinationHistograms(
      const std::vector<Eigen::Vector3d>& positions,
      const Eigen::MatrixXi& knn_indices,   // (K, N), reuse the existing KNN matrix
      const std::vector<int64_t>& flat_indices,
      const std::vector<int>& class_labels, // per-point: +1=leaf, -1=wood, 0=unknown
      int n_bins,
      double l_max = 0.05               // max edge length to accept a facet (metres)
  );

  // Solve Bailey & Mahaffee (2017) eq.10 for plant area density.
  // Inputs: raw voxel accumulators — no new fields needed.
  //
  // Variable mapping (paper symbol -> local variable):
  //   r_bar  = path_length / num_beams_weighted  (mean scan path through voxel)
  //   P_bar  = 1 - (num_hits / num_beams_weighted)         (mean gap probability)
  //   G      = projection coefficient from IAD histogram
  //
  // Consistency check: at thin-medium limit, -ln(P_bar)/r_bar -> num_hits/path_length
  // (recovers the existing Vicari estimator). Bailey is preferred at higher attenuation.
  //
  // Iterative form: secant method on residual f(a_L) per eq.10.
  // Seed: a_L0 = -ln(P_bar) / (G * r_bar)   (closed-form thin-limit)
  // Convergence: |f| < 1e-9 or 50 iterations.
  // Guard: return 0.0 if G <= 0 or r_bar <= 0 or num_beams_weighted < 1.
  double solveBaileyPadEq10(double path_length,
                            double num_beams_weighted,
                            double num_hits,
                            double G);

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASBAILEY_H
