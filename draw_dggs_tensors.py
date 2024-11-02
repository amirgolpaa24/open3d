import open3d as o3d

from icosahedron import *


# resolution
target_res = 3
tensor_res = 2
tensor_children_count = 4**(target_res - tensor_res)

icosahedron = create_subdivided_icosahedron(resolution=target_res, radius=10.0)
icosahedron_copy = create_subdivided_icosahedron(resolution=target_res, radius=10.0)
icosahedron_colored_tensors = []

for f in np.arange(0, 20 * 4**(target_res), tensor_children_count):
    tensor_color = [np.random.rand(), np.random.rand(), np.random.rand()]
    tensor_mesh = make_colored_triangles(icosahedron, selected_triangles=range(f, f + tensor_children_count), selected_color=tensor_color, remove_from_original=True)
    icosahedron_colored_tensors.append(tensor_mesh)

# edges
line_set, cylinders = make_edges(icosahedron_copy, offset=0.0, edge_thickness=0.003)

################################################

# Create a visualizer
vis = o3d.visualization.Visualizer()
vis.create_window()

# Add the meshes to the visualizer
vis.add_geometry(icosahedron)
for tensor_mesh in icosahedron_colored_tensors:
    vis.add_geometry(tensor_mesh)

# Add the edges to the visualizer
vis.add_geometry(line_set)

# Run the visualizer
vis.run()

# Destroy the visualizer window
vis.destroy_window()
