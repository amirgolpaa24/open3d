import open3d as o3d
import numpy as np


# Create a sphere mesh
sphere = o3d.geometry.TriangleMesh.create_sphere(radius=1.0)

# Compute vertex normals for better lighting effects
sphere.compute_vertex_normals()

# Optionally apply a color
sphere.paint_uniform_color([0.1, 0.8, 0.1])  # RGB values (green)

# Create a visualizer
vis = o3d.visualization.Visualizer()
vis.create_window()

# Add the sphere to the visualizer
vis.add_geometry(sphere)

# Run the visualizer
vis.run()

# Destroy the visualizer window
vis.destroy_window()
