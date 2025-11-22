// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raysegment.h"
#include <nabo/nabo.h>
#include "rayterrain.h"
#include <queue>
#include <fstream> // For debug output
#include <set>
#include "raytrees.h" // For convertIntToColour

namespace ray
{
// Helper function to save debug PLY files for root consolidation
void saveRootDebugPly(const std::string& filename,
                      const std::string& description,
                      const std::vector<Vertex>& all_points,
                      const Eigen::Vector3d& offset,
                      const std::vector<int>& point_indices_to_save, // Pass indices of points we care about
                      const std::vector<int>& cluster_ids_for_points) // Map from point index to cluster id
{
    std::ofstream ofs(filename.c_str());
    if (!ofs.is_open()) {
        std::cerr << "Error: Cannot open debug root file " << filename << std::endl;
        return;
    }

    // Header
    ofs << "ply\nformat ascii 1.0\n";
    ofs << "comment " << description << "\n";
    ofs << "element vertex " << point_indices_to_save.size() << "\n";
    ofs << "property float x\nproperty float y\nproperty float z\n";
    ofs << "property uchar red\nproperty uchar green\nproperty uchar blue\n";
    ofs << "end_header\n";

    // Vertex data
    for (const int point_idx : point_indices_to_save)
    {
        const Eigen::Vector3d& pos = all_points[point_idx].pos;
        int cluster_id = cluster_ids_for_points[point_idx];

        RGBA colour = {0, 0, 0, 255}; // Black for unassigned/filtered, with full alpha
        if (cluster_id != -1) {
            convertIntToColour(cluster_id, colour);
        }

        ofs << pos.x() + offset.x() << " "
            << pos.y() + offset.y() << " "
            << pos.z() + offset.z() << " "
            << (int)colour.red << " "
            << (int)colour.green << " "
            << (int)colour.blue << "\n";
    }
    std::cout << "Saved root debug file: " << filename << std::endl;
}

// Helper function to save Eigen integer arrays as a point cloud grid for debugging
void saveGridAsPly(const std::string& filename, const std::string& description,
                   const Eigen::ArrayXXi& grid,
                   const Eigen::Vector3d& box_min, double pixel_width, const Eigen::Vector3d& offset)
{
    std::ofstream ofs(filename.c_str());
    if (!ofs.is_open()) {
        std::cerr << "Error: Cannot open debug grid file " << filename << std::endl;
        return;
    }
    
    int point_count = (grid.array() > 0).count();
    
    // Header
    ofs << "ply\nformat ascii 1.0\n";
    ofs << "comment " << description << "\n";
    ofs << "element vertex " << point_count << "\n";
    ofs << "property float x\nproperty float y\nproperty float z\n";
    ofs << "property int value\n";
    ofs << "end_header\n";

    // Vertex data
    for (int i = 0; i < grid.rows(); ++i) {
        for (int j = 0; j < grid.cols(); ++j) {
            if (grid(i, j) > 0) {
                double x = box_min.x() + i * pixel_width + offset.x();
                double y = box_min.y() + j * pixel_width + offset.y();
                double z = offset.z(); // Flat grid
                ofs << x << " " << y << " " << z << " " << grid(i, j) << "\n";
            }
        }
    }
    std::cout << "Saved grid debug file: " << filename << std::endl;
}


/// nodes of priority queue used in shortest path algorithm
struct QueueNode
{
//  QueueNode() {}
  QueueNode(double distance_to_source, double score, double radius, int root, int index)
    : distance_to_source(distance_to_source)
    , score(score)
    , radius(radius)
    , root(root)
    , id(index)
  {}

  double distance_to_source; // path distance to the source (root or tip)
  double score;              // score is the modified edge length metric being minimised
  double radius;             // radius of the tree base, this acts as a score scale coefficient
  int root;                  // index of the root of the path
  int id;                    // index into the points_ array for this node
};

class QueueNodeComparator
{
public:
  bool operator()(const QueueNode &p1, const QueueNode &p2) { return p1.score > p2.score; }
};

/// Connect the supplied set of points @c points according to the shortest path to the ground (bottom-up)
void connectPointsShortestPath_BottomUp(
  std::vector<Vertex> &points,
  std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> &closest_node,
  const Eigen::MatrixXi& indices, const Eigen::MatrixXd& dists2,
  int search_size, double gravity_factor)
{
  while (!closest_node.empty())
  {
    QueueNode node = closest_node.top();
    closest_node.pop();
    if (!points[node.id].visited)
    {
      for (int i = 0; i < search_size && indices(i, node.id) != Nabo::NNSearchD::InvalidIndex; i++)
      {
        const int child = indices(i, node.id);
        const double dist2 = dists2(i, node.id);
        const double dist = std::sqrt(dist2);
        
        const Eigen::Vector3d dif = (points[child].pos - points[node.id].pos).normalized();
        Eigen::Vector3d dir(0, 0, 1);
        const int ppar = points[node.id].parent;
        if (ppar != -1)
        {
          if (points[ppar].parent != -1)
          {
            dir = (points[node.id].pos - points[points[ppar].parent].pos).normalized();
          }
          else
          {
            dir = (points[node.id].pos - points[ppar].pos).normalized();
          }
        }
        const double d = std::max(0.001, dif.dot(dir));
        double score = dist2 / (d * d * (double)points[child].weight);

        if (gravity_factor > 0.0)
        {
          Eigen::Vector3d to_node = points[node.id].pos - points[node.root].pos;
          to_node[2] = 0.0;
          const double lateral_sqr = to_node.squaredNorm();
          const double gravity_scale = 1.0 + gravity_factor * lateral_sqr;
          score *= gravity_scale;
        }

        score /= node.radius;

        double new_score = node.score + score;
        if (new_score < points[child].score)
        {
          points[child].score = new_score;
          points[child].distance_to_ground = node.distance_to_source + dist;
          points[child].parent = node.id;
          points[child].root = node.root;
          closest_node.push(
            QueueNode(points[child].distance_to_ground, points[child].score, node.radius, node.root, child));
        }
      }
      points[node.id].visited = true;
    }
  }
}

/// Connect the supplied set of points @c points according to the shortest path to the tips (top-down)
void connectPointsShortestPath_TopDown(
  std::vector<Vertex> &points,
  std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> &closest_node,
  const Eigen::MatrixXi& indices, const Eigen::MatrixXd& dists2,
  int search_size, double gravity_factor)
{
  while (!closest_node.empty())
  {
    QueueNode node = closest_node.top();
    closest_node.pop();
    if (!points[node.id].visited)
    {
      for (int i = 0; i < search_size && indices(i, node.id) != Nabo::NNSearchD::InvalidIndex; i++)
      {
        const int child = indices(i, node.id); // 'child' is now conceptually the parent
        const double dist2 = dists2(i, node.id);
        const double dist = std::sqrt(dist2);
        
        const Eigen::Vector3d dif = (points[child].pos - points[node.id].pos).normalized();
        Eigen::Vector3d dir(0, 0, -1); // Default direction is down

        const int ppar = points[node.id].parent_to_tip;
        if (ppar != -1)
        {
          if (points[ppar].parent_to_tip != -1)
          {
            dir = (points[node.id].pos - points[points[ppar].parent_to_tip].pos).normalized();
          }
          else
          {
            dir = (points[node.id].pos - points[ppar].pos).normalized();
          }
        }
        
        const double d = std::max(0.001, dif.dot(dir));
        double score = dist2 / (d * d * (double)points[child].weight);

        if (gravity_factor > 0.0)
        {
           Eigen::Vector3d to_node = points[node.id].pos - points[node.root].pos;
           double gravity_assist = 1.0 - (to_node.z() / std::max(0.1, to_node.norm()));
           to_node[2] = 0.0;
           const double lateral_sqr = to_node.squaredNorm();
           const double gravity_scale = 1.0 + gravity_factor * lateral_sqr;
           score *= gravity_scale * std::max(0.1, gravity_assist);
        }

        double new_score = node.score + score;
        if (new_score < points[child].score_to_tip)
        {
          points[child].score_to_tip = new_score;
          points[child].distance_to_tip = node.distance_to_source + dist;
          points[child].parent_to_tip = node.id;
          points[child].tip_root = node.root;
          closest_node.push(
            QueueNode(points[child].distance_to_tip, new_score, 1.0, node.root, child));
        }
      }
      points[node.id].visited = true;
    }
  }
}

/// Applies the results from the top-down path search to the main parent links.
void applyVerificationCorrections(std::vector<Vertex>& points)
{
    int changes = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i].parent_to_tip != -1) { 
            if (points[i].parent != points[i].parent_to_tip) {
                changes++;
            }
            points[i].parent = points[i].parent_to_tip;
            points[i].root = points[i].tip_root; 
        }
    }
    std::cout << "Applied " << changes << " parent changes from top-down path search." << std::endl;
}

/// Saves the current path graph to a PLY file with edges for debugging.
void savePathsToMesh(const std::string& filename, const std::string& description, const std::vector<Vertex>& points, bool use_top_down_parents, const Eigen::Vector3d& offset)
{
    std::ofstream ofs(filename.c_str(), std::ios::binary);
    if (!ofs.is_open()) {
        std::cerr << "Error: Cannot open debug path file " << filename << std::endl;
        return;
    }
    
    // Header
    ofs << "ply\nformat binary_little_endian 1.0\n";
    ofs << "comment " << description << "\n";
    ofs << "element vertex " << points.size() << "\n";
    ofs << "property float x\nproperty float y\nproperty float z\n";
    ofs << "property uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha\n";
    
    size_t edge_count = 0;
    for(const auto& p : points) {
        int parent_id = use_top_down_parents ? p.parent_to_tip : p.parent;
        if (parent_id != -1) {
            edge_count++;
        }
    }
    ofs << "element edge " << edge_count << "\n";
    ofs << "property int vertex1\nproperty int vertex2\n";
    ofs << "end_header\n";

    // Vertex data
    #pragma pack(push, 1)
    struct PlyVert { float x, y, z; ray::RGBA c; };
    #pragma pack(pop)
    
    std::vector<std::vector<int>> children(points.size());
     for (size_t i = 0; i < points.size(); ++i) {
        int parent_id = use_top_down_parents ? points[i].parent_to_tip : points[i].parent;
        if (parent_id != -1 && parent_id < (int)points.size()) { // Safety check
            children[parent_id].push_back(i);
        }
    }

    for (size_t i = 0; i < points.size(); ++i) {
        PlyVert vert;
        vert.x = static_cast<float>(points[i].pos.x() + offset.x());
        vert.y = static_cast<float>(points[i].pos.y() + offset.y());
        vert.z = static_cast<float>(points[i].pos.z() + offset.z());
        
        int parent_id = use_top_down_parents ? points[i].parent_to_tip : points[i].parent;
        if (parent_id == -1) {
             vert.c = {0, 0, 255, 255}; // Blue for sources (roots or tips)
        } else if (children[i].empty()) {
            vert.c = {255, 0, 0, 255}; // Red for sinks
        } else {
            vert.c = {0, 255, 0, 255}; // Green for path
        }
        
        ofs.write(reinterpret_cast<const char*>(&vert), sizeof(PlyVert));
    }

    // Edge data
    for (size_t i = 0; i < points.size(); ++i) {
        int parent_id = use_top_down_parents ? points[i].parent_to_tip : points[i].parent;
        if (parent_id != -1) {
            int vertex1 = static_cast<int>(i);
            int vertex2 = parent_id;
            ofs.write(reinterpret_cast<const char*>(&vertex1), sizeof(int));
            ofs.write(reinterpret_cast<const char*>(&vertex2), sizeof(int));
        }
    }
}


/// Converts a ray cloud to a set of points @c points connected by the shortest path to the ground @c mesh
/// the returned vector of index sets provides the root points for each separated tree
std::vector<std::vector<int>> getRootsAndSegment(std::vector<Vertex> &points, const Cloud &cloud, const Mesh &mesh,
                                                 double max_diameter, double distance_limit, double height_min,
                                                 double gravity_factor, bool alpha_weighting, bool verify_paths,
                                                 bool debug_paths, const std::string& name_stub, const Eigen::Vector3d& offset,
                                                 bool root_debug)
{
  if (root_debug) {
    std::cout << "Root debug mode enabled. This will generate multiple large PLY files." << std::endl;
  }
  // first fill in the basic attributes of the points structure
  points.reserve(cloud.ends.size());
  for (unsigned int i = 0; i < cloud.ends.size(); i++)
  {
    if (cloud.rayBounded(i))
    {
      uint8_t weight = 1;
      if (alpha_weighting && cloud.colours[i].alpha > 0)
        weight = cloud.colours[i].alpha;
      points.push_back(Vertex(cloud.ends[i], cloud.starts[i], weight));
    }
  }

  const double pixel_width = max_diameter;
  Eigen::Vector3d box_min, box_max;
  cloud.calcBounds(&box_min, &box_max);
  
  // also add points for every vertex on the ground mesh.
  const int roots_start = static_cast<int>(points.size());
  for (auto &vert : mesh.vertices())
  {
    if (vert[0] >= box_min[0] && vert[1] >= box_min[1] &&
        vert[0] <= box_max[0] && vert[1] <= box_max[1])
    {
      points.push_back(Vertex(vert, vert + Eigen::Vector3d(0,0,0.01), 1)); // make small vertical ray
    }
  }

  // Perform k-NN search once and pass it to the functions
  const int search_size = std::min(20, static_cast<int>(points.size()) - 1);
  Eigen::MatrixXd points_p(3, points.size());
  for (unsigned int i = 0; i < points.size(); i++)
  {
    points_p.col(i) = points[i].pos;
  }
  Nabo::NNSearchD *nns = Nabo::NNSearchD::createKDTreeLinearHeap(points_p, 3);
  Eigen::MatrixXi indices;
  Eigen::MatrixXd dists2;
  indices.resize(search_size, points.size());
  dists2.resize(search_size, points.size());
  nns->knn(points_p, indices, dists2, search_size, kNearestNeighbourEpsilon, 0, distance_limit);
  delete nns;

  // convert the ground mesh to an easy look-up height field
  Eigen::ArrayXXd lowfield;
  mesh.toHeightField(lowfield, box_min, box_max, pixel_width);

  // set heightfield as the height of the canopy above the ground
  Eigen::ArrayXXd heightfield =
    Eigen::ArrayXXd::Constant(static_cast<int>(lowfield.rows()), static_cast<int>(lowfield.cols()), std::numeric_limits<double>::lowest());
  for (const auto &point : points)
  {
    Eigen::Vector3i index = ((point.pos - box_min) / pixel_width).cast<int>();
    heightfield(index[0], index[1]) = std::max(heightfield(index[0], index[1]), point.pos[2]);
  }
  // make heightfield relative to the ground
  for (int i = 0; i < heightfield.rows(); i++)
  {
    for (int j = 0; j < heightfield.cols(); j++)
    {
      heightfield(i, j) = std::max(1e-10, heightfield(i, j) - lowfield(i, j));
    }
  }

  // --- PASS 1: BOTTOM-UP SEARCH ---
  std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> closest_node_bottom_up;
  for(auto& p : points) { p.visited = false; } // Reset visited flags
  for (int ind = roots_start; ind < static_cast<int>(points.size()); ind++)
  {
    points[ind].distance_to_ground = 0.0;
    points[ind].score = 0.0;
    points[ind].root = ind;
    const Eigen::Vector3i index = ((points[ind].pos - box_min) / pixel_width).cast<int>();
    closest_node_bottom_up.push(QueueNode(0, 0, heightfield(index[0], index[1]), ind, ind));
  }
  connectPointsShortestPath_BottomUp(points, closest_node_bottom_up, indices, dists2, search_size, gravity_factor);

  if(debug_paths) {
      savePathsToMesh(name_stub + "_debug_paths_original.ply", "Original bottom-up paths", points, false, offset);
  }

  // --- PASS 2: TOP-DOWN SEARCH (SANDBOXED) ---
  for(auto& p : points) { p.visited = false; }

  std::vector<std::vector<int>> children(points.size());
  for (size_t i = 0; i < points.size(); ++i) {
      if(points[i].parent != -1) {
          children[points[i].parent].push_back(i);
      }
  }
  
  std::priority_queue<QueueNode, std::vector<QueueNode>, QueueNodeComparator> closest_node_top_down;
  for (size_t i = 0; i < points.size(); ++i) {
      if (points[i].parent != -1 && children[i].empty()) {
          points[i].score_to_tip = 0.0;
          points[i].tip_root = i;
          points[i].distance_to_tip = 0.0;
          closest_node_top_down.push(QueueNode(0, 0, 1.0, i, i));
      }
  }
  connectPointsShortestPath_TopDown(points, closest_node_top_down, indices, dists2, search_size, gravity_factor);

  if(debug_paths) {
      savePathsToMesh(name_stub + "_debug_paths_verified.ply", "Sandboxed top-down paths", points, true, offset);
  }

  // --- (OPTIONAL) PASS 3: APPLY CORRECTIONS ---
  if (verify_paths) {
      applyVerificationCorrections(points);
  }


  // --- ROOT CONSOLIDATION (uses final `parent` and `root` links) ---
  std::vector<int> point_to_cluster_id(points.size(), -1);
  std::vector<int> all_point_indices;
  for(size_t i = 0; i < points.size(); ++i) all_point_indices.push_back(i);
  std::vector<int> ground_point_indices;
  for (int i = roots_start; i < static_cast<int>(points.size()); i++) ground_point_indices.push_back(i);
  
  Eigen::ArrayXXi counts =
    Eigen::ArrayXXi::Constant(static_cast<int>(heightfield.rows()), static_cast<int>(heightfield.cols()), 0);
  Eigen::ArrayXXd heights =
    Eigen::ArrayXXd::Constant(static_cast<int>(heightfield.rows()), static_cast<int>(heightfield.cols()), 0);
  for (const auto &point : points)
  {
    if (point.root == -1)
    {
      continue;
    }
    const Eigen::Vector3i index = ((points[point.root].pos - box_min) / pixel_width).cast<int>();
    counts(index[0], index[1])++;
    heights(index[0], index[1]) = std::max(heights(index[0], index[1]), point.pos[2] - lowfield(index[0], index[1]));
  }

  Eigen::ArrayXXi sums = Eigen::ArrayXXi::Constant(static_cast<int>(counts.rows()), static_cast<int>(counts.cols()), 0);
  for (int i = 0; i < static_cast<int>(sums.rows()); i++)
  {
    for (int j = 0; j < static_cast<int>(sums.cols()); j++)
    {
      const int i2 = std::min(i + 1, static_cast<int>(sums.rows()) - 1);
      const int j2 = std::min(j + 1, static_cast<int>(sums.cols()) - 1);
      sums(i, j) = counts(i, j) + counts(i, j2) + counts(i2, j) + counts(i2, j2);
    }
  }
  std::vector<Eigen::Vector2i> bests(static_cast<int>(counts.rows()) * static_cast<int>(counts.cols()));
  for (int x = 0; x < static_cast<int>(sums.rows()); x++)
  {
    for (int y = 0; y < static_cast<int>(sums.cols()); y++)
    {
      Eigen::Vector2i best_index(-1, -1);
      int largest_sum = -1;
      for (int i = std::max(0, x - 1); i <= x; i++)
      {
        for (int j = std::max(0, y - 1); j <= y; j++)
        {
          if (sums(i, j) > largest_sum)
          {
            largest_sum = sums(i, j);
            best_index = Eigen::Vector2i(i, j);
          }
        }
      }
      bests[x + sums.rows() * y] = best_index;
    }
  }

  if (root_debug)
  {
      saveGridAsPly(name_stub + "_debug_grid_counts.ply", "Root density grid",
                    counts, box_min, pixel_width, offset);
      saveGridAsPly(name_stub + "_debug_grid_sums.ply", "Smoothed root density grid",
                    sums, box_min, pixel_width, offset);
      
      Eigen::ArrayXXi bests_grid(sums.rows(), sums.cols());
      for (int x = 0; x < static_cast<int>(sums.rows()); x++) {
          for (int y = 0; y < static_cast<int>(sums.cols()); y++) {
              const Eigen::Vector2i& best_idx = bests[x + sums.rows() * y];
              bests_grid(x, y) = best_idx.x() + sums.rows() * best_idx.y();
          }
      }
      saveGridAsPly(name_stub + "_debug_grid_bests.ply", "Best root neighborhood grid",
                    bests_grid, box_min, pixel_width, offset);

      // Populate the cluster ID map for all points
      for(size_t i = 0; i < points.size(); ++i)
      {
          if (points[i].root != -1)
          {
              const Eigen::Vector3i index = ((points[points[i].root].pos - box_min) / pixel_width).cast<int>();
              const Eigen::Vector2i best_index = bests[index[0] + static_cast<int>(sums.rows()) * index[1]];
              point_to_cluster_id[i] = best_index[0] + static_cast<int>(sums.rows()) * best_index[1];
          } else {
              point_to_cluster_id[i] = -1;
          }
      }

      // OUTPUT 1: All ground root points coloured by their consolidated cluster.
      saveRootDebugPly(name_stub + "_debug_roots_1_all_consolidated.ply",
                       "All ground points coloured by consolidated cluster ID",
                       points, offset, ground_point_indices, point_to_cluster_id);

      // OUTPUT 3: The entire point cloud coloured by consolidated cluster.
      saveRootDebugPly(name_stub + "_debug_roots_3_cloud_consolidated.ply",
                       "Full point cloud coloured by consolidated cluster ID",
                       points, offset, all_point_indices, point_to_cluster_id);
  }

  Eigen::ArrayXXd max_heights =
    Eigen::ArrayXXd::Constant(static_cast<int>(counts.rows()), static_cast<int>(counts.cols()), 0);
  for (int i = 0; i < static_cast<int>(sums.rows()); i++)
  {
    for (int j = 0; j < static_cast<int>(sums.cols()); j++)
    {
      Eigen::Vector2i best_index = bests[i + sums.rows() * j];
      double max_height = 0.0;
      for (int x = best_index[0]; x < std::min(best_index[0] + 2, static_cast<int>(sums.rows())); x++)
      {
        for (int y = best_index[1]; y < std::min(best_index[1] + 2, static_cast<int>(sums.cols())); y++)
        {
          if (bests[x + static_cast<int>(sums.rows()) * y] == best_index)
          {
            max_height = std::max(max_height, heights(x, y));
          }
        }
      }
      max_heights(best_index[0], best_index[1]) = max_height;
    }
  }

  std::vector<std::vector<int>> roots_lists(sums.rows() * sums.cols());
  for (int i = roots_start; i < static_cast<int>(points.size()); i++)
  {
    const Eigen::Vector3i index = ((points[i].pos - box_min) / pixel_width).cast<int>();
    const Eigen::Vector2i best_index = bests[index[0] + static_cast<int>(sums.rows()) * index[1]];
    const double max_height = max_heights(best_index[0], best_index[1]);
    if (max_height >= height_min)
    {
      const int id = best_index[0] + static_cast<int>(sums.rows()) * best_index[1];
      roots_lists[id].push_back(i);
    }
  }

  if (root_debug)
  {
      // Create a set of valid cluster IDs that passed the height filter
      std::set<int> valid_cluster_ids;
      for (int i = 0; i < static_cast<int>(sums.rows()); i++) {
          for (int j = 0; j < static_cast<int>(sums.cols()); j++) {
              int id = i + static_cast<int>(sums.rows()) * j;
              if (!roots_lists[id].empty()) {
                  valid_cluster_ids.insert(id);
              }
          }
      }

      // Create a new cluster ID map, filtering out the invalid ones
      std::vector<int> filtered_point_to_cluster_id = point_to_cluster_id;
      for (size_t i = 0; i < points.size(); ++i) {
          if (point_to_cluster_id[i] != -1 && valid_cluster_ids.find(point_to_cluster_id[i]) == valid_cluster_ids.end()) {
              filtered_point_to_cluster_id[i] = -1; // Set to black/invalid
          }
      }

      // OUTPUT 2: Ground root points after height filtering.
      saveRootDebugPly(name_stub + "_debug_roots_2_height_filtered.ply",
                       "Ground points after height filtering (removed roots are black)",
                       points, offset, ground_point_indices, filtered_point_to_cluster_id);
      
      // OUTPUT 4: The entire point cloud after height filtering.
      saveRootDebugPly(name_stub + "_debug_roots_4_cloud_filtered.ply",
                       "Full point cloud after height filtering (removed trees are black)",
                       points, offset, all_point_indices, filtered_point_to_cluster_id);
  }

  std::vector<std::vector<int>> roots_set;
  for (int i = 0; i < static_cast<int>(sums.rows()); i++)
  {
    for (int j = 0; j < static_cast<int>(sums.cols()); j++)
    {
      auto &roots = roots_lists[i + static_cast<int>(sums.rows()) * j];
      if (roots.size() > 0)
      {
        roots_set.push_back(roots);
      }
    }
  }

  return roots_set;
}

}  // namespace ray