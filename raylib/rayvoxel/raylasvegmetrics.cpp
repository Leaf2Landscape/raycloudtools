// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file implements advanced vegetation metrics calculations.

#include "raylib/rayvoxel/raylasvegmetrics.h"
#include "raylib/rayutils.h" // For kPi
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <iomanip>

namespace ray
{

// Anonymous namespace for helper functions local to this file
namespace {

// --- Numerical integration (trapezoidal rule) ---
double trapezoidal_integral(const std::function<double(double)>& f, double a, double b, int n) {
    if (n <= 0) {
        return 0.0;
    }
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i)
    {
        sum += f(a + i * h);
    }
    return h * sum;
}

// --- Leaf Angle Distribution Probability Density Functions (PDFs) ---
// Based on amapvox's LeafAngleDistribution values

double dplanophile(double thetaL) {
    return (2.0 / kPi) * (1.0 + cos(2.0 * thetaL));
}

double derectophile(double thetaL) {
    return (2.0 / kPi) * (1.0 - cos(2.0 * thetaL));
}

double dplagiophile(double thetaL) {
    return (2.0 / kPi) * (1.0 - cos(4.0 * thetaL));
}

double dextremophile(double thetaL) {
    return (2.0 / kPi) * (1.0 + cos(4.0 * thetaL));
}

double dspherical(double thetaL) {
    return sin(thetaL);
}

double duniform(double thetaL) {
    return 2.0 / kPi;
}

double dellipsoidal(double thetaL, double chi) {
    if (chi == 1.0) return sin(thetaL);

    double epsilon, lambda;
    if (chi < 1.0) {
        epsilon = sqrt(1.0 - chi * chi);
        lambda = chi + asin(epsilon) / epsilon;
    } else {
        epsilon = sqrt(1.0 - 1.0 / (chi * chi));
        lambda = chi + log((1.0 + epsilon) / (1.0 - epsilon)) / (2.0 * epsilon * chi);
    }
    return 2.0 * (chi * chi * chi) * sin(thetaL) / (lambda * pow(cos(thetaL) * cos(thetaL) + (chi * chi * sin(thetaL) * sin(thetaL)), 2));
}

// Beta function B(a,b) using lgamma (log-gamma) for numerical stability
double beta_func(double a, double b) {
    return std::exp(std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
}

// Probability density function for the Beta distribution
double dbeta(double x, double mu, double nu) {
    if (x < 0.0 || x > 1.0) {
        return 0.0;
    }
    if (mu <= 0 || nu <= 0) {
        throw std::invalid_argument("Beta distribution parameters mu and nu must be positive.");
    }
    double B = beta_func(mu, nu);
    if (B == 0) return 0; // Should not happen for mu, nu > 0
    return (pow(x, mu - 1.0) * pow(1.0 - x, nu - 1.0)) / B;
}

double dtwoParamBeta(double thetaL, double mu, double nu) {
    double t = 2.0 * thetaL / kPi;
    return dbeta(t, mu, nu) * (2.0 / kPi);
}

// G-function projection kernel A(theta, thetaL). Extracted so it can be reused
// by computeGFromHistogram without re-implementing the projection geometry.
static double projectionKernelA(double theta, double thetaL) {
  double cotcot = 1.0 / (std::tan(theta) * std::tan(thetaL));
  if (std::abs(cotcot) > 1.0 || std::isinf(cotcot)) {
    return std::cos(theta) * std::cos(thetaL);
  }
  double acos_cotcot = std::acos(cotcot);
  return std::cos(theta) * std::cos(thetaL) * (1.0 + (2.0 / kPi) * (std::tan(acos_cotcot) - acos_cotcot));
}

} // anonymous namespace


double computeG(double theta, const std::string& lad, double param1, double param2)
{
  // For spherical LAD, G(theta) is always 0.5, a significant optimization.
  if (lad == "spherical") {
      return 0.5;
  }

  // Normalize theta to the range [0, pi/2]
  theta = fmod(theta, kPi);
  if (theta > (kPi / 2.0)) {
    theta = kPi - theta;
  }

  // Avoid tan(pi/2) issues by slightly adjusting the angle if it's exactly at the boundary
  if (theta >= kPi / 2.0) {
      theta = kPi / 2.0 - 1e-9;
  }

  // Define the G-function kernel, A(theta, thetaL)
  auto A = [theta](double thetaL){ return projectionKernelA(theta, thetaL); };

  // Define the full function to be integrated: A(theta, thetaL) * g_L(thetaL)
  std::function<double(double)> integrand;
  if (lad == "planophile") integrand = [&](double thetaL){ return A(thetaL) * dplanophile(thetaL); };
  else if (lad == "erectophile") integrand = [&](double thetaL){ return A(thetaL) * derectophile(thetaL); };
  else if (lad == "plagiophile") integrand = [&](double thetaL){ return A(thetaL) * dplagiophile(thetaL); };
  else if (lad == "extremophile") integrand = [&](double thetaL){ return A(thetaL) * dextremophile(thetaL); };
  else if (lad == "uniform") integrand = [&](double thetaL){ return A(thetaL) * duniform(thetaL); };
  else if (lad == "ellipsoidal") integrand = [&](double thetaL){ return A(thetaL) * dellipsoidal(thetaL, param1); };
  else if (lad == "twoParamBeta") integrand = [&](double thetaL){ return A(thetaL) * dtwoParamBeta(thetaL, param1, param2); };
  else {
      static bool warned = false;
      if (!warned) {
          std::cerr << "Warning: Unsupported Leaf Angle Distribution (LAD) '" << lad << "'. Defaulting to 'spherical'." << std::endl;
          warned = true;
      }
      return 0.5; // Default to spherical
  }

  // Perform numerical integration over the leaf angle domain [0, pi/2]
  // The number of steps (180) is chosen to match amapvox for consistency.
  int integration_steps = 180;
  return trapezoidal_integral(integrand, 0.0, kPi / 2.0, integration_steps);
}

double computeGFromHistogram(double theta_beam,
                             const std::vector<double>& bin_centres,
                             const std::vector<double>& liad)
{
  if (bin_centres.empty() || liad.empty() || bin_centres.size() != liad.size())
    return 0.5;
  theta_beam = std::fmod(theta_beam, kPi);
  if (theta_beam > (kPi / 2.0)) theta_beam = kPi - theta_beam;
  if (theta_beam >= kPi / 2.0)  theta_beam = kPi / 2.0 - 1e-9;
  double G = 0.0;
  for (size_t b = 0; b < bin_centres.size(); ++b)
    G += projectionKernelA(theta_beam, bin_centres[b]) * liad[b];
  return G;
}

std::string encodeIadToJson(const std::vector<double>& bin_centres_deg,
                            const std::vector<double>& values)
{
  if (bin_centres_deg.empty() || values.empty() || bin_centres_deg.size() != values.size())
    return "{}";
  std::ostringstream oss;
  oss << "{";
  for (size_t i = 0; i < bin_centres_deg.size(); ++i) {
    if (i > 0) oss << ",";
    oss << "\"" << std::fixed << std::setprecision(1) << bin_centres_deg[i] << "\":"
        << std::fixed << std::setprecision(6) << values[i];
  }
  oss << "}";
  return oss.str();
}

LaserSpecManager::LaserSpecManager()
{
    // Populate the map with predefined laser specifications from AMAPVox.
    // Note: AMAPVox uses several names for the same spec. We will support the primary ones.
    specs_["LMS-Q560"] = {"LMS_Q560", 0.0003, 0.0005};
    specs_["LMS-Q780"] = {"LMS_Q780", 0.005, 0.00025};
    specs_["VZ-400"] = {"VZ_400/VZ_400i", 0.007, 0.00035};
    specs_["VZ-400i"] = {"VZ_400/VZ_400i", 0.007, 0.00035};
    specs_["LEICA-SCANSTATION-P30-40"] = {"LEICA_SCANSTATION_P30_40", 0.0035, 0.00023};
    specs_["LEICA-SCANSTATION-C10"] = {"LEICA_SCANSTATION_C10", 0.004, 0.0001};
    specs_["FARO-FOCUS-X330"] = {"FARO_FOCUS_X330", 0.0025, 0.00019};
    specs_["MINIVUX-1UAV"] = {"miniVUX-1UAV", 0.0145, 0.00105};
    specs_["TRIMBLE-X7"] = {"TRIMBLE_X7", 0.0026, 0.0008};
    specs_["UNITARY-BEAM-SECTION"] = {"Unitary beam section", 0.0, 0.0};
}

std::string classifyDeWit(const std::vector<double>& bin_centres,
                           const std::vector<double>& hist)
{
  if (bin_centres.empty() || hist.size() != bin_centres.size()) return "";
  double total = 0.0;
  for (double v : hist) total += v;
  if (total <= 0.0) return "";

  const int n = static_cast<int>(bin_centres.size());

  using PdfFn = double(*)(double);
  const std::pair<const char*, PdfFn> candidates[] = {
    {"planophile",   dplanophile},
    {"erectophile",  derectophile},
    {"plagiophile",  dplagiophile},
    {"extremophile", dextremophile},
    {"spherical",    dspherical},
    {"uniform",      duniform},
  };

  std::string best_name;
  double best_dist = std::numeric_limits<double>::max();

  for (const auto& [name, pdf] : candidates) {
    std::vector<double> ref(n);
    double ref_sum = 0.0;
    for (int b = 0; b < n; ++b) { ref[b] = pdf(bin_centres[b]); ref_sum += ref[b]; }
    if (ref_sum <= 0.0) continue;
    double d2 = 0.0;
    for (int b = 0; b < n; ++b) {
      double diff = hist[b] / total - ref[b] / ref_sum;
      d2 += diff * diff;
    }
    if (d2 < best_dist) { best_dist = d2; best_name = name; }
  }
  return best_name;
}

bool LaserSpecManager::getSpec(const std::string& name, LaserSpecification& spec_out) const
{
    // Create a case-insensitive version of the name for matching
    std::string upper_name = name;
    std::transform(upper_name.begin(), upper_name.end(), upper_name.begin(), ::toupper);

    auto it = specs_.find(upper_name);
    if (it != specs_.end()) {
        spec_out = it->second;
        return true;
    }
    return false;
}

} // namespace ray
