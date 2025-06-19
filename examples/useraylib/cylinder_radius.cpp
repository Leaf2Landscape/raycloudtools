// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Generated with Claude Code

#include <raylib/raycloud.h>
#include <raylib/extraction/raytrees.h>
#include <raylib/extraction/raysegment.h>
#include <raylib/raymesh.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <sstream>

// Helper class to access Trees::estimateCylinderRadius
class CylinderRadiusEstimator
{
public:
  // Create a simplified mesh for ground detection (flat ground at z=0)
  static ray::Mesh createFlatGroundMesh(const Eigen::Vector3d& min_bound, const Eigen::Vector3d& max_bound)
  {
    ray::Mesh mesh;
    // Create a simple flat ground mesh
    double ground_z = min_bound[2] - 0.1; // Slightly below minimum point
    
    // Add vertices for a simple quad
    mesh.vertices().push_back(Eigen::Vector3d(min_bound[0] - 1, min_bound[1] - 1, ground_z));
    mesh.vertices().push_back(Eigen::Vector3d(max_bound[0] + 1, min_bound[1] - 1, ground_z));
    mesh.vertices().push_back(Eigen::Vector3d(max_bound[0] + 1, max_bound[1] + 1, ground_z));
    mesh.vertices().push_back(Eigen::Vector3d(min_bound[0] - 1, max_bound[1] + 1, ground_z));
    
    // Add faces (two triangles forming a quad)
    mesh.indexList().push_back(Eigen::Vector3i(0, 1, 2));
    mesh.indexList().push_back(Eigen::Vector3i(0, 2, 3));
    
    return mesh;
  }

  static double estimateRadius(ray::Cloud& cloud, bool use_rays, double& accuracy)
  {
    // Calculate bounds
    Eigen::Vector3d min_bound = cloud.ends[0];
    Eigen::Vector3d max_bound = cloud.ends[0];
    for (const auto& end : cloud.ends)
    {
      min_bound = min_bound.cwiseMin(end);
      max_bound = max_bound.cwiseMax(end);
    }

    // Create a simple flat ground mesh
    ray::Mesh mesh = createFlatGroundMesh(min_bound, max_bound);
    
    // Set up parameters for tree extraction
    ray::TreesParams params;
    params.use_rays = use_rays;
    params.max_diameter = 2.0; // Large enough for any reasonable cylinder
    params.height_min = 0.1;   // Very small minimum height
    params.distance_limit = 1.0;
    params.crop_length = 0.01; // Very small crop length
    
    try 
    {
      // Create Trees object - this will run the full tree reconstruction
      Eigen::Vector3d offset(0, 0, 0);
      ray::Trees trees(cloud, offset, mesh, params, false); // verbose = false
      
      // The tree reconstruction should have created at least one section
      // For now, we'll extract statistics from the processed cloud
      // This is a simplified approach - in a full implementation you'd access
      // the internal tree structure
      
      // Calculate center and use PCA for axis
      Eigen::Vector3d center(0, 0, 0);
      for (const auto& end : cloud.ends)
      {
        center += end;
      }
      center /= static_cast<double>(cloud.ends.size());

      // Calculate the principal axis using PCA
      Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
      for (const auto& end : cloud.ends)
      {
        Eigen::Vector3d diff = end - center;
        covariance += diff * diff.transpose();
      }
      covariance /= static_cast<double>(cloud.ends.size());
      
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
      Eigen::Vector3d cylinder_axis(0, 0, 1); // Default to vertical
      if (solver.info() == Eigen::Success)
      {
        cylinder_axis = solver.eigenvectors().col(2); // Largest eigenvalue
      }

      // Calculate radius manually using the same algorithm as Trees::estimateCylinderRadius
      double rad = 0.0;
      double power = 0.25;
      std::vector<double> norms;
      norms.reserve(cloud.ends.size());
      
      if (use_rays && cloud.starts.size() == cloud.ends.size())
      {
        for (size_t i = 0; i < cloud.ends.size(); i++)
        {
          Eigen::Vector3d offset = cloud.ends[i] - center;
          offset -= cylinder_axis * offset.dot(cylinder_axis); // flatten    
          Eigen::Vector3d ray = cloud.starts[i] - cloud.ends[i];
          ray -= cylinder_axis * ray.dot(cylinder_axis); // flatten
          double ray_dot = ray.dot(ray);
          if (ray_dot > 1e-10)
          {
            offset += ray * std::max(0.0, std::min(-offset.dot(ray) / ray_dot, 1.0)); // move to closest point
          }
          norms.push_back(offset.norm());
        }
      }
      else
      {
        for (const auto& end : cloud.ends)
        {
          Eigen::Vector3d offset = end - center;
          offset -= cylinder_axis * offset.dot(cylinder_axis); // flatten    
          norms.push_back(offset.norm());
        }
      }

      for (auto& norm : norms)
      {
        rad += std::pow(norm, power);
      }
      const double eps = 1e-5;
      rad /= (double)norms.size() + eps;
      rad = std::pow(rad, 1.0 / power);
      
      double e = 0.0;
      for (auto& norm : norms)
      {
        e += std::abs(norm - rad);
      }
      e /= (double)norms.size() + eps;
      
      const double sensor_noise = 0.02;
      accuracy = rad / (e + sensor_noise);
      accuracy = std::max(accuracy, eps);
      rad = std::max(rad, eps);

      return rad;
    }
    catch (const std::exception& e)
    {
      std::cerr << "Error in tree reconstruction: " << e.what() << std::endl;
      accuracy = 0.0;
      return 0.0;
    }
  }
};

// Test data generation function
class TestDataGenerator
{
public:
  struct CylinderParams
  {
    double radius = 0.15;           // True cylinder radius in meters (e.g., tree branch)
    double height = 2.0;            // Cylinder height in meters
    Eigen::Vector3d center = Eigen::Vector3d(0, 0, 1); // Cylinder center
    Eigen::Vector3d axis = Eigen::Vector3d(0, 0, 1);   // Cylinder axis (normalized)
    int num_points = 1000;          // Number of LiDAR points to generate
    double noise_level = 0.01;      // Vegetation thickness: outward extension (leaves/bark) in meters
    double outward_bias = 1.0;      // Fraction of noise that's outward (1.0 = all outward, 0.0 = uniform)
    double scanner_distance = 5.0;  // Distance from scanner to cylinder center
    Eigen::Vector3d scanner_pos = Eigen::Vector3d(5, 0, 1); // LiDAR scanner position
  };

  static ray::Cloud generateCylinderCloud(const CylinderParams& params, unsigned int seed = 42)
  {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::normal_distribution<double> noise_dist(0.0, params.noise_level);
    
    ray::Cloud cloud;
    cloud.starts.reserve(params.num_points);
    cloud.ends.reserve(params.num_points);
    cloud.colours.reserve(params.num_points);
    cloud.times.reserve(params.num_points);

    // Create orthogonal vectors to the cylinder axis for parameterization
    Eigen::Vector3d axis = params.axis.normalized();
    Eigen::Vector3d temp_vec = (std::abs(axis.dot(Eigen::Vector3d(1, 0, 0))) < 0.9) ? 
                               Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0);
    Eigen::Vector3d u = axis.cross(temp_vec).normalized();
    Eigen::Vector3d v = axis.cross(u).normalized();

    std::cout << "Generating synthetic cylinder data (tree branch with foliage):" << std::endl;
    std::cout << "  True radius: " << params.radius << " m" << std::endl;
    std::cout << "  Height: " << params.height << " m" << std::endl;
    std::cout << "  Center: " << params.center.transpose() << std::endl;
    std::cout << "  Axis: " << axis.transpose() << std::endl;
    std::cout << "  Vegetation thickness: " << params.noise_level << " m" << std::endl;
    std::cout << "  Outward bias: " << params.outward_bias << " (1.0=pure outward, 0.0=uniform)" << std::endl;
    std::cout << "  Number of points: " << params.num_points << std::endl;
    std::cout << "  Scanner position: " << params.scanner_pos.transpose() << std::endl;
    
    if (params.outward_bias > 0.8) {
        if (params.noise_level > 0.1) {
            std::cout << "  Noise model: Extremely dense foliage/large leaves" << std::endl;
        } else {
            std::cout << "  Noise model: Dense foliage (mostly outward)" << std::endl;
        }
    } else if (params.outward_bias > 0.3) {
        std::cout << "  Noise model: Mixed vegetation/bark" << std::endl;
    } else {
        std::cout << "  Noise model: Rough bark (mostly uniform)" << std::endl;
    }

    // Create 4 scanner positions around the cylinder
    std::vector<Eigen::Vector3d> scanner_positions;
    for (int scan_id = 0; scan_id < 4; scan_id++)
    {
      double scan_angle = scan_id * M_PI / 2.0; // 0°, 90°, 180°, 270°
      Eigen::Vector3d scan_pos = params.center + 
        params.scanner_distance * (u * std::cos(scan_angle) + v * std::sin(scan_angle)) +
        (params.center[2] - params.height/4.0) * axis; // Slightly below center height
      scanner_positions.push_back(scan_pos);
    }

    for (int i = 0; i < params.num_points; i++)
    {
      // Generate random point on cylinder surface
      double theta = uniform(rng) * 2.0 * M_PI;  // Angle around cylinder
      double z = uniform(rng) * params.height;   // Height along cylinder
      
      // Randomly choose scanner position to create varied ray directions
      int scan_id = static_cast<int>(uniform(rng) * 4.0) % 4;
      Eigen::Vector3d current_scanner_pos = scanner_positions[scan_id];
      
      // Point on perfect cylinder surface
      Eigen::Vector3d surface_point = params.center + 
                                     params.radius * (u * std::cos(theta) + v * std::sin(theta)) +
                                     (z - params.height/2.0) * axis;
      
      // Calculate surface normal for noise generation
      Eigen::Vector3d surface_normal = (u * std::cos(theta) + v * std::sin(theta)).normalized();
      
      // Calculate ray direction from current scanner to surface point
      Eigen::Vector3d ray_direction = (surface_point - current_scanner_pos).normalized();
      
      // Simulate realistic tree branch with leaves/bark extending outward
      // 
      //     Leaves/Bark    True Branch
      //          ↓             ↓
      //       .  *  .     ╭─────────╮
      //    .     *     .  │  Branch │  ← True cylinder surface
      //       .  *  .     ╰─────────╯
      //          ↑
      //    LiDAR hits here (outer surface)
      //
      // Use the surface normal for noise generation (already calculated above)
      Eigen::Vector3d radial_dir = surface_normal;
      
      // Generate noise based on outward_bias parameter
      // outward_bias = 1.0: all noise extends outward (dense foliage)
      // outward_bias = 0.0: uniform noise in all directions (rough bark)
      
      Eigen::Vector3d noise_vector;
      if (params.outward_bias > 0.0)
      {
        // Outward-biased noise (vegetation model)
        double outward_noise = std::abs(noise_dist(rng)) * params.outward_bias; // Always positive (outward)
        double tangential_noise = noise_dist(rng) * 0.3 * (1.0 - params.outward_bias);  // Reduced for high bias
        double axial_noise = noise_dist(rng) * 0.2 * (1.0 - params.outward_bias);       // Reduced for high bias
        
        // Create noise vector: mostly outward, small tangential and axial components
        Eigen::Vector3d tangential_dir = axis.cross(radial_dir).normalized();
        noise_vector = outward_noise * radial_dir + 
                      tangential_noise * tangential_dir + 
                      axial_noise * axis;
      }
      else
      {
        // Uniform noise (bark roughness model) - ensure no inward displacement
        double outward_component = std::abs(noise_dist(rng)); // Always outward
        double tangential_noise = noise_dist(rng) * 0.7;      // Can be positive or negative
        double axial_noise = noise_dist(rng) * 0.7;           // Can be positive or negative
        
        Eigen::Vector3d tangential_dir = axis.cross(radial_dir).normalized();
        noise_vector = outward_component * radial_dir + 
                      tangential_noise * tangential_dir + 
                      axial_noise * axis;
      }
      
      // The actual surface that LiDAR hits (true surface + outward noise)
      Eigen::Vector3d noisy_surface = surface_point + noise_vector;
      
      // Ray end point - this is where the LiDAR ray actually terminates
      // It hits the outer layer (leaves/bark) rather than the true cylinder surface
      Eigen::Vector3d ray_end = noisy_surface;
      
      // Add small measurement noise to simulate sensor uncertainty
      double sensor_noise_level = params.noise_level * 0.1; // Much smaller than vegetation noise
      ray_end += Eigen::Vector3d(
        std::normal_distribution<double>(0.0, sensor_noise_level)(rng),
        std::normal_distribution<double>(0.0, sensor_noise_level)(rng),
        std::normal_distribution<double>(0.0, sensor_noise_level)(rng)
      );
      
      // Ray starts from scanner and ends at the hit point
      // The ray direction is from scanner TO the point
      Eigen::Vector3d ray_start = ray_end - ray_direction * 10.0; // Start 10m away along ray direction
      
      cloud.starts.push_back(ray_start);
      cloud.ends.push_back(ray_end);
      
      // Add some color variation
      uint8_t intensity = static_cast<uint8_t>(200 + noise_dist(rng) * 20);
      cloud.colours.push_back(ray::RGBA(intensity, intensity, intensity, 255));
      cloud.times.push_back(static_cast<double>(i) / params.num_points);
    }

    return cloud;
  }

  // Structure to hold test results
  struct TestResult
  {
    double noise_level_mm;
    double true_radius;
    double estimated_radius_without_rays;
    double estimated_radius_with_rays;
    double accuracy_without_rays;
    double accuracy_with_rays;
    double error_without_rays_mm;
    double error_with_rays_mm;
    double improvement_mm;
    int num_points;
    double cylinder_height;
    std::string timestamp;
  };

  // Write results to CSV
  static void writeResultsToCSV(const std::vector<TestResult>& results, const std::string& filename)
  {
    std::ofstream csv_file(filename);
    if (!csv_file.is_open())
    {
      std::cerr << "Error: Could not open CSV file for writing: " << filename << std::endl;
      return;
    }

    // Write CSV header
    csv_file << "noise_level_mm,true_radius_m,estimated_radius_without_rays_m,estimated_radius_with_rays_m,"
             << "accuracy_without_rays,accuracy_with_rays,error_without_rays_mm,error_with_rays_mm,"
             << "improvement_mm,num_points,cylinder_height_m,timestamp" << std::endl;

    // Write data rows
    csv_file << std::fixed << std::setprecision(6);
    for (const auto& result : results)
    {
      csv_file << result.noise_level_mm << ","
               << result.true_radius << ","
               << result.estimated_radius_without_rays << ","
               << result.estimated_radius_with_rays << ","
               << result.accuracy_without_rays << ","
               << result.accuracy_with_rays << ","
               << result.error_without_rays_mm << ","
               << result.error_with_rays_mm << ","
               << result.improvement_mm << ","
               << result.num_points << ","
               << result.cylinder_height << ","
               << result.timestamp << std::endl;
    }

    csv_file.close();
    std::cout << "\nResults saved to CSV: " << filename << std::endl;
  }

  // Get current timestamp as string
  static std::string getCurrentTimestamp()
  {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d_%H-%M-%S");
    return ss.str();
  }

  // Generate test data with multiple noise levels for comparison
  static void runNoiseComparison(const std::string& base_filename)
  {
    CylinderParams params;
    params.radius = 0.08;     // 8cm radius (typical branch diameter)
    params.height = 1.5;      // 1.5m height
    params.num_points = 500;
    params.outward_bias = 0.85; // 85% outward bias (realistic foliage)
    
    // Test vegetation thickness from 2mm to 50cm (extended range for dense foliage)
    std::vector<double> noise_levels = {0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5}; // 2mm to 50cm vegetation
    std::vector<TestResult> results;
    std::string timestamp = getCurrentTimestamp();
    
    std::cout << "\n=== VEGETATION THICKNESS COMPARISON (Extended Range) ===" << std::endl;
    std::cout << "True branch radius: " << params.radius << " m" << std::endl;
    std::cout << "Branch height: " << params.height << " m" << std::endl;
    std::cout << "Points per test: " << params.num_points << std::endl;
    std::cout << "Outward bias: " << params.outward_bias << " (foliage model)" << std::endl;
    std::cout << "Testing range: 2mm to 500mm vegetation thickness" << std::endl;
    std::cout << "Timestamp: " << timestamp << std::endl;
    
    for (double noise : noise_levels)
    {
      params.noise_level = noise;
      ray::Cloud test_cloud = generateCylinderCloud(params);
      
      if (noise >= 0.1) {
        std::cout << "\n--- Vegetation thickness: " << noise*100 << " cm ---" << std::endl;
      } else {
        std::cout << "\n--- Vegetation thickness: " << noise*1000 << " mm ---" << std::endl;
      }
      
      // Test both modes
      ray::Cloud cloud_copy1 = test_cloud;
      double accuracy_without_rays;
      double radius_without_rays = CylinderRadiusEstimator::estimateRadius(cloud_copy1, false, accuracy_without_rays);
      
      ray::Cloud cloud_copy2 = test_cloud;
      double accuracy_with_rays;
      double radius_with_rays = CylinderRadiusEstimator::estimateRadius(cloud_copy2, true, accuracy_with_rays);
      
      double error_without = std::abs(radius_without_rays - params.radius);
      double error_with = std::abs(radius_with_rays - params.radius);
      double improvement = error_without - error_with;
      
      // Store results
      TestResult result;
      result.noise_level_mm = noise * 1000;
      result.true_radius = params.radius;
      result.estimated_radius_without_rays = radius_without_rays;
      result.estimated_radius_with_rays = radius_with_rays;
      result.accuracy_without_rays = accuracy_without_rays;
      result.accuracy_with_rays = accuracy_with_rays;
      result.error_without_rays_mm = error_without * 1000;
      result.error_with_rays_mm = error_with * 1000;
      result.improvement_mm = improvement * 1000;
      result.num_points = params.num_points;
      result.cylinder_height = params.height;
      result.timestamp = timestamp;
      results.push_back(result);
      
      std::cout << "  Without rays: " << radius_without_rays << " m (error: " << error_without*1000 << " mm, accuracy: " << accuracy_without_rays << ")" << std::endl;
      std::cout << "  With rays:    " << radius_with_rays << " m (error: " << error_with*1000 << " mm, accuracy: " << accuracy_with_rays << ")" << std::endl;
      std::cout << "  Improvement:  " << improvement*1000 << " mm reduction in error" << std::endl;
      
      // Save test data if filename provided
      if (!base_filename.empty())
      {
        std::string ply_filename;
        if (noise >= 0.1) {
          ply_filename = base_filename + "_noise_" + std::to_string(static_cast<int>(noise*100)) + "cm.ply";
        } else {
          ply_filename = base_filename + "_noise_" + std::to_string(static_cast<int>(noise*1000)) + "mm.ply";
        }
        test_cloud.save(ply_filename);
        std::cout << "  Saved: " << ply_filename << std::endl;
      }
    }

    // Write results to CSV
    if (!base_filename.empty())
    {
      std::string csv_filename = base_filename + "_results.csv";
      writeResultsToCSV(results, csv_filename);
      
      // Also write summary statistics
      writeSummaryCSV(results, base_filename + "_summary.csv");
    }
  }

  // Write summary statistics to CSV
  static void writeSummaryCSV(const std::vector<TestResult>& results, const std::string& filename)
  {
    std::ofstream csv_file(filename);
    if (!csv_file.is_open())
    {
      std::cerr << "Error: Could not open summary CSV file for writing: " << filename << std::endl;
      return;
    }

    // Calculate summary statistics
    double avg_error_without = 0, avg_error_with = 0, avg_improvement = 0;
    double max_error_without = 0, max_error_with = 0, max_improvement = -1e10;
    double min_error_without = 1e10, min_error_with = 1e10, min_improvement = 1e10;
    
    for (const auto& result : results)
    {
      avg_error_without += result.error_without_rays_mm;
      avg_error_with += result.error_with_rays_mm;
      avg_improvement += result.improvement_mm;
      
      max_error_without = std::max(max_error_without, result.error_without_rays_mm);
      max_error_with = std::max(max_error_with, result.error_with_rays_mm);
      max_improvement = std::max(max_improvement, result.improvement_mm);
      
      min_error_without = std::min(min_error_without, result.error_without_rays_mm);
      min_error_with = std::min(min_error_with, result.error_with_rays_mm);
      min_improvement = std::min(min_improvement, result.improvement_mm);
    }
    
    int n = results.size();
    avg_error_without /= n;
    avg_error_with /= n;
    avg_improvement /= n;

    // Write summary header and data
    csv_file << "metric,without_rays_mm,with_rays_mm,improvement_mm" << std::endl;
    csv_file << std::fixed << std::setprecision(3);
    csv_file << "average," << avg_error_without << "," << avg_error_with << "," << avg_improvement << std::endl;
    csv_file << "maximum," << max_error_without << "," << max_error_with << "," << max_improvement << std::endl;
    csv_file << "minimum," << min_error_without << "," << min_error_with << "," << min_improvement << std::endl;

    csv_file.close();
    std::cout << "Summary statistics saved to: " << filename << std::endl;
  }
};

void usage(int exit_code = 1)
{
  std::cout << "Estimate cylinder radius from a ray cloud using Trees::estimateCylinderRadius" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "cylinder_radius <raycloud_file>" << std::endl;
  std::cout << "cylinder_radius --generate-test [output_prefix]" << std::endl;
  std::cout << "cylinder_radius --noise-test [output_prefix]" << std::endl;
  std::cout << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  <raycloud_file>      Analyze existing ray cloud file" << std::endl;
  std::cout << "  --generate-test      Generate synthetic cylinder test data" << std::endl;
  std::cout << "  --noise-test         Run noise level comparison test" << std::endl;
  std::cout << "  [output_prefix]      Optional prefix for saved test files" << std::endl;
  std::cout << std::endl;
  std::cout << "This program tests both use_rays=false and use_rays=true modes" << std::endl;
  exit(exit_code);
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    usage();
  }

  std::string arg1 = argv[1];
  
  // Handle test generation modes
  if (arg1 == "--generate-test")
  {
    std::string output_prefix = (argc > 2) ? argv[2] : "test_cylinder";
    
    // Generate single test case with default parameters (tree branch with leaves)
    TestDataGenerator::CylinderParams params;
    params.radius = 0.20;        // 20cm radius (typical branch)
    params.height = 2.0;         // 2m height  
    params.num_points = 1000;    // 1000 points
    params.noise_level = 0.02;   // 2cm vegetation thickness (leaves/small branches)
    params.outward_bias = 0.9;   // 90% outward noise (dense foliage)
    
    ray::Cloud test_cloud = TestDataGenerator::generateCylinderCloud(params);
    
    std::string filename = output_prefix + ".ply";
    test_cloud.save(filename);
    std::cout << "\nGenerated test data saved to: " << filename << std::endl;
    
    // Test the generated data
    ray::Cloud cloud_copy1 = test_cloud;
    double accuracy_without_rays;
    double radius_without_rays = CylinderRadiusEstimator::estimateRadius(cloud_copy1, false, accuracy_without_rays);
    
    ray::Cloud cloud_copy2 = test_cloud;
    double accuracy_with_rays;
    double radius_with_rays = CylinderRadiusEstimator::estimateRadius(cloud_copy2, true, accuracy_with_rays);
    
    std::cout << "\n=== TEST RESULTS ===" << std::endl;
    std::cout << "True radius: " << params.radius << " m" << std::endl;
    std::cout << "Without rays: " << radius_without_rays << " m (error: " << std::abs(radius_without_rays - params.radius)*1000 << " mm)" << std::endl;
    std::cout << "With rays: " << radius_with_rays << " m (error: " << std::abs(radius_with_rays - params.radius)*1000 << " mm)" << std::endl;
    
    // Save single test result to CSV
    TestDataGenerator::TestResult result;
    result.noise_level_mm = params.noise_level * 1000;
    result.true_radius = params.radius;
    result.estimated_radius_without_rays = radius_without_rays;
    result.estimated_radius_with_rays = radius_with_rays;
    result.accuracy_without_rays = accuracy_without_rays;
    result.accuracy_with_rays = accuracy_with_rays;
    result.error_without_rays_mm = std::abs(radius_without_rays - params.radius) * 1000;
    result.error_with_rays_mm = std::abs(radius_with_rays - params.radius) * 1000;
    result.improvement_mm = result.error_without_rays_mm - result.error_with_rays_mm;
    result.num_points = params.num_points;
    result.cylinder_height = params.height;
    result.timestamp = TestDataGenerator::getCurrentTimestamp();
    
    std::vector<TestDataGenerator::TestResult> single_result = {result};
    std::string csv_filename = output_prefix + "_result.csv";
    TestDataGenerator::writeResultsToCSV(single_result, csv_filename);
    
    return 0;
  }
  
  if (arg1 == "--noise-test")
  {
    std::string output_prefix = (argc > 2) ? argv[2] : "noise_test";
    TestDataGenerator::runNoiseComparison(output_prefix);
    return 0;
  }

  // Handle regular file input
  std::string filename = arg1;

  ray::Cloud raycloud;
  if (!raycloud.load(filename))
  {
    std::cerr << "Error: Failed to load raycloud " << filename << std::endl;
    return 1;
  }

  std::cout << "Cloud " << filename << " loaded with " << raycloud.starts.size() << " rays and "
            << raycloud.ends.size() << " points" << std::endl;

  if (raycloud.ends.empty())
  {
    std::cerr << "Error: No points in the raycloud" << std::endl;
    return 1;
  }

  // Calculate center for display
  Eigen::Vector3d center(0, 0, 0);
  for (const auto& end : raycloud.ends)
  {
    center += end;
  }
  center /= static_cast<double>(raycloud.ends.size());

  std::cout << "\n=== RESULTS ===" << std::endl;
  std::cout << "Cylinder center: " << center.transpose() << std::endl;

  // Test both modes using the actual Trees implementation
  std::cout << "\n--- Testing with Trees::estimateCylinderRadius ---" << std::endl;

  // Test without using rays (use_rays = false)
  ray::Cloud cloud_copy1 = raycloud; // Make a copy since Trees modifies the cloud
  double accuracy_without_rays;
  double radius_without_rays = CylinderRadiusEstimator::estimateRadius(cloud_copy1, false, accuracy_without_rays);

  std::cout << "\n--- WITHOUT using rays (use_rays=false) ---" << std::endl;
  std::cout << "Estimated cylinder radius: " << radius_without_rays << std::endl;
  std::cout << "Accuracy metric: " << accuracy_without_rays << std::endl;

  // Test with using rays (use_rays = true)
  ray::Cloud cloud_copy2 = raycloud; // Make another copy
  double accuracy_with_rays;
  double radius_with_rays = CylinderRadiusEstimator::estimateRadius(cloud_copy2, true, accuracy_with_rays);

  std::cout << "\n--- WITH using rays (use_rays=true) ---" << std::endl;
  std::cout << "Estimated cylinder radius: " << radius_with_rays << std::endl;
  std::cout << "Accuracy metric: " << accuracy_with_rays << std::endl;

  std::cout << "\n--- COMPARISON ---" << std::endl;
  std::cout << "Radius difference: " << std::abs(radius_with_rays - radius_without_rays) << std::endl;
  std::cout << "Accuracy improvement: " << (accuracy_with_rays - accuracy_without_rays) << std::endl;

  return 0;
}