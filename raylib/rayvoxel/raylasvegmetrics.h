// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file declares advanced vegetation metrics calculations, such as the
// G-function for Plant Area Density (PAD) estimation.

#ifndef RAYLIB_RAYVOXEL_RAYLASVEGMETRICS_H
#define RAYLIB_RAYVOXEL_RAYLASVEGMETRICS_H

#include <string>
#include <map>
#include <vector>

namespace ray
{
  /// @brief Computes the foliage projection ratio G(theta).
  /// This function calculates the mean projection of unit leaf area onto the plane
  /// perpendicular to the beam direction, assuming a symmetric distribution of leaf azimuth angles.
  /// @param theta The incident beam inclination angle (zenith angle), in radians, ranging from [0, pi/2].
  /// @param lad The name of the probability density function for the Leaf Angle Distribution (LAD).
  ///            Supported: "spherical", "uniform", "planophile", "erectophile", "plagiophile",
  ///            "extremophile", "ellipsoidal", "twoParamBeta".
  /// @param param1 Primary parameter for the LAD (e.g., chi for ellipsoidal, mu for beta).
  /// @param param2 Secondary parameter for the LAD (e.g., nu for beta).
  /// @return The G(theta) value.
  double computeG(double theta, const std::string& lad = "spherical", double param1 = 0.0, double param2 = 0.0);

  /// @brief Computes G(theta) via dot-product of projection kernel A against an
  /// empirical inclination histogram. Returns 0.5 for empty or degenerate input.
  double computeGFromHistogram(double theta_beam,
                               const std::vector<double>& bin_centres,
                               const std::vector<double>& liad);

  /// @brief Serializes a histogram as {"angle_deg": fraction, ...} JSON. Returns "{}" for empty.
  std::string encodeIadToJson(const std::vector<double>& bin_centres_deg,
                              const std::vector<double>& values);

  /// @struct LaserSpecification
  /// @brief Holds the physical properties of a laser scanner's beam.
  struct LaserSpecification
  {
      std::string name;
      double beam_diameter_at_exit; // in meters
      double beam_divergence;       // in radians
  };

  /// @class LaserSpecManager
  /// @brief Manages a collection of predefined laser specifications.
  class LaserSpecManager
  {
  public:
      LaserSpecManager();

      /// @brief Retrieves a laser specification by its name.
      /// @param name The name of the predefined laser spec.
      /// @param spec_out The LaserSpecification object to populate if found.
      /// @return True if the specification was found, false otherwise.
      bool getSpec(const std::string& name, LaserSpecification& spec_out) const;

  private:
      std::map<std::string, LaserSpecification> specs_;
  };

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASVEGMETRICS_H
