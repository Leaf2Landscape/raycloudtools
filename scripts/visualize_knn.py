#!/usr/bin/env python3
"""
Visualize k-nearest neighbor connections from a point cloud.
Reads a PLY file and creates a new PLY file with edges connecting each point to its k nearest neighbors.
"""

import numpy as np
from scipy.spatial import cKDTree
import argparse
import sys
import struct

def read_ply(filename):
    """Read a PLY file (ASCII or binary) and return points and colors."""
    print(f"Reading {filename}...")

    # Read header in binary mode to handle both formats
    with open(filename, 'rb') as f:
        # Read header
        line = f.readline().decode('ascii').strip()
        if line != 'ply':
            raise ValueError("Not a valid PLY file")

        vertex_count = 0
        properties = []
        is_binary = False
        is_little_endian = True
        in_header = True
        header_end_pos = 0

        while in_header:
            line = f.readline().decode('ascii').strip()
            if line.startswith('format'):
                parts = line.split()
                if parts[1] == 'binary_little_endian':
                    is_binary = True
                    is_little_endian = True
                elif parts[1] == 'binary_big_endian':
                    is_binary = True
                    is_little_endian = False
                elif parts[1] == 'ascii':
                    is_binary = False
            elif line.startswith('element vertex'):
                vertex_count = int(line.split()[2])
            elif line.startswith('property'):
                parts = line.split()
                prop_type = parts[1]
                prop_name = parts[2]
                properties.append({'name': prop_name, 'type': prop_type})
            elif line == 'end_header':
                in_header = False
                header_end_pos = f.tell()

    if is_binary:
        return read_binary_ply(filename, header_end_pos, vertex_count, properties, is_little_endian)
    else:
        return read_ascii_ply(filename, vertex_count, properties)

def read_ascii_ply(filename, vertex_count, properties):
    """Read ASCII PLY format."""
    with open(filename, 'r') as f:
        # Skip header
        for line in f:
            if line.strip() == 'end_header':
                break

        # Read vertex data
        points = []
        colors = []

        prop_names = [p['name'] for p in properties]

        for i in range(vertex_count):
            line = f.readline().strip().split()

            # Extract x, y, z
            x_idx = prop_names.index('x')
            y_idx = prop_names.index('y')
            z_idx = prop_names.index('z')

            point = [float(line[x_idx]), float(line[y_idx]), float(line[z_idx])]
            points.append(point)

            # Extract colors if available
            color = [128, 128, 128]  # default gray
            try:
                r_idx = prop_names.index('red')
                g_idx = prop_names.index('green')
                b_idx = prop_names.index('blue')
                color = [int(line[r_idx]), int(line[g_idx]), int(line[b_idx])]
            except ValueError:
                pass
            colors.append(color)

    points = np.array(points)
    colors = np.array(colors)

    print(f"Loaded {len(points)} points (ASCII format)")
    return points, colors

def read_binary_ply(filename, header_end_pos, vertex_count, properties, is_little_endian):
    """Read binary PLY format."""
    endian_char = '<' if is_little_endian else '>'

    # Map PLY types to struct format characters
    type_map = {
        'float': 'f',
        'double': 'd',
        'uchar': 'B',
        'char': 'b',
        'uint': 'I',
        'int': 'i',
        'ushort': 'H',
        'short': 'h'
    }

    # Build struct format string
    fmt = endian_char
    for prop in properties:
        fmt += type_map[prop['type']]

    vertex_size = struct.calcsize(fmt)
    prop_names = [p['name'] for p in properties]

    # Read binary data
    points = []
    colors = []

    with open(filename, 'rb') as f:
        f.seek(header_end_pos)

        for i in range(vertex_count):
            data = struct.unpack(fmt, f.read(vertex_size))

            # Extract x, y, z
            x_idx = prop_names.index('x')
            y_idx = prop_names.index('y')
            z_idx = prop_names.index('z')

            point = [data[x_idx], data[y_idx], data[z_idx]]
            points.append(point)

            # Extract colors if available
            color = [128, 128, 128]  # default gray
            try:
                r_idx = prop_names.index('red')
                g_idx = prop_names.index('green')
                b_idx = prop_names.index('blue')
                color = [int(data[r_idx]), int(data[g_idx]), int(data[b_idx])]
            except ValueError:
                pass
            colors.append(color)

    points = np.array(points)
    colors = np.array(colors)

    print(f"Loaded {len(points)} points (binary format)")
    return points, colors

def write_ply_with_edges(filename, points, colors, edges):
    """Write a PLY file with vertices and edges."""
    print(f"Writing {filename}...")

    with open(filename, 'w') as f:
        # Write header
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write(f"comment K-nearest neighbor connections (k={args.k})\n")
        f.write(f"element vertex {len(points)}\n")
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property uchar red\n")
        f.write("property uchar green\n")
        f.write("property uchar blue\n")
        f.write(f"element edge {len(edges)}\n")
        f.write("property int vertex1\n")
        f.write("property int vertex2\n")
        f.write("end_header\n")

        # Write vertices
        for i in range(len(points)):
            f.write(f"{points[i, 0]:.6f} {points[i, 1]:.6f} {points[i, 2]:.6f} ")
            f.write(f"{colors[i, 0]} {colors[i, 1]} {colors[i, 2]}\n")

        # Write edges
        for edge in edges:
            f.write(f"{edge[0]} {edge[1]}\n")

    print(f"Saved {len(edges)} edges")

def compute_knn_edges(points, k=20):
    """Compute k-nearest neighbor edges for all points."""
    print(f"Building k-d tree for {len(points)} points...")
    tree = cKDTree(points)

    print(f"Finding {k} nearest neighbors for each point...")
    edges = []

    # Query k+1 neighbors (including the point itself)
    distances, indices = tree.query(points, k=k+1)

    # Build edges, avoiding duplicates
    edge_set = set()

    for i in range(len(points)):
        # Skip first neighbor (which is the point itself)
        for j in range(1, min(k+1, len(indices[i]))):
            neighbor_idx = indices[i][j]

            # Store edge in canonical form (smaller index first) to avoid duplicates
            edge = (min(i, neighbor_idx), max(i, neighbor_idx))
            edge_set.add(edge)

    edges = list(edge_set)
    print(f"Generated {len(edges)} unique edges")

    return edges

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Visualize k-nearest neighbor connections from a PLY point cloud"
    )
    parser.add_argument("input", help="Input PLY file")
    parser.add_argument("-o", "--output", help="Output PLY file (default: input_knn.ply)")
    parser.add_argument("-k", "--k", type=int, default=20,
                       help="Number of nearest neighbors (default: 20)")

    args = parser.parse_args()

    # Determine output filename
    if args.output is None:
        if args.input.endswith('.ply'):
            args.output = args.input[:-4] + '_knn.ply'
        else:
            args.output = args.input + '_knn.ply'

    try:
        # Read input PLY
        points, colors = read_ply(args.input)

        # Compute k-nearest neighbor edges
        edges = compute_knn_edges(points, k=args.k)

        # Write output PLY with edges
        write_ply_with_edges(args.output, points, colors, edges)

        print(f"\nSuccess! K-NN visualization saved to {args.output}")
        print(f"You can view this file in CloudCompare, MeshLab, or similar PLY viewers")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
