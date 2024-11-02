import open3d as o3d

import numpy as np


def subdivide_mesh(mesh):
    # Get vertices and triangles of the mesh
    vertices = np.asarray(mesh.vertices)
    triangles = np.asarray(mesh.triangles)

    # Create a dictionary to store new vertex indices
    new_vertices = {}
    new_triangles = []

    vertices_initial_len = len(vertices)

    # Loop through each triangle
    for triangle in triangles:
        # Get the vertices of the triangle
        v0, v1, v2 = vertices[triangle]

        # Find midpoints of each edge
        m01 = (v0 + v1) / 2
        m12 = (v1 + v2) / 2
        m20 = (v2 + v0) / 2

        # Add midpoints to the dictionary if they are not already there
        midpoints = [m01, m12, m20]
        midpoint_indices = []

        for m in midpoints:
            key = tuple(np.round(m, decimals=5))  # Use a tuple for dictionary keys
            if key in new_vertices:
                midpoint_indices.append(new_vertices[key])
            else:
                new_index = vertices_initial_len + len(new_vertices)
                new_vertices[key] = new_index
                midpoint_indices.append(new_index)

                # Add the midpoint to the vertex list
                vertices = np.vstack((vertices, m))

        # Create new triangles
        new_triangles.append([triangle[0], midpoint_indices[0], midpoint_indices[2]])
        new_triangles.append([triangle[1], midpoint_indices[1], midpoint_indices[0]])
        new_triangles.append([triangle[2], midpoint_indices[2], midpoint_indices[1]])
        new_triangles.append([midpoint_indices[0], midpoint_indices[1], midpoint_indices[2]])

    # Create a new mesh with the updated vertices and triangles
    new_mesh = o3d.geometry.TriangleMesh()
    new_mesh.vertices = o3d.utility.Vector3dVector(vertices)
    new_mesh.triangles = o3d.utility.Vector3iVector(new_triangles)

    return new_mesh

def subdivide_mesh_many(mesh, level=1):
    if level == 0:
        return mesh
    elif level == 1:
        return subdivide_mesh(mesh)
    else:
        return subdivide_mesh(subdivide_mesh_many(mesh, level - 1))

def normalize_mesh(mesh):
    # Get vertices of the mesh
    vertices = np.asarray(mesh.vertices)

    # Normalize vertices
    vertices /= np.linalg.norm(vertices, axis=1)[:, None]

    # Update the mesh vertices
    mesh.vertices = o3d.utility.Vector3dVector(vertices)

    return mesh

### meshes
def create_icosahedron(radius=1.0):
    # Create an icosahedron mesh
    icosahedron = o3d.geometry.TriangleMesh.create_icosahedron(radius=radius)
    # normalize the vertices
    icosahedron = normalize_mesh(icosahedron)
    # Compute vertex normals for the original mesh
    icosahedron.compute_vertex_normals()

    return icosahedron

def create_subdivided_icosahedron(resolution=0, radius=1.0):
    # initial icosahedron
    icosahedron = create_icosahedron(radius=radius)

    # Subdivide the mesh
    icosahedron_subdivided = subdivide_mesh_many(icosahedron, level=resolution)
    # normalize the vertices
    icosahedron_subdivided = normalize_mesh(icosahedron_subdivided)
    # Compute vertex normals for the new mesh
    icosahedron_subdivided.compute_vertex_normals()

    return icosahedron_subdivided

def make_colored_triangles(mesh, selected_triangles=[], selected_color=[1, 0, 0], remove_from_original=False):
    # current mesh: attributes
    vertices = np.asarray(mesh.vertices)
    triangles = np.asarray(mesh.triangles)
    vertex_colors = np.asarray(mesh.vertex_colors)
    num_verts = len(vertices)
    num_triangles = len(triangles)

    # we create a new mesh for just the new colored triangels
    mesh_new = o3d.geometry.TriangleMesh()
    # new mesh: attributes
    vertices_new = np.zeros((num_verts, 3))
    triangles_new = np.zeros((num_triangles, 3))
    vertex_colors_new = np.zeros((num_verts, 3))

    # mapping of old -> new vertex indices
    vertex_map = {}

    for triangle_idx_new, triangle_idx in enumerate(selected_triangles):
        (vert_index1, vert_index2, vert_index3) = triangles[triangle_idx]

        # mapping vertex indices
        if vert_index1 not in vertex_map:
            vertex_map[vert_index1] = len(vertex_map)
        if vert_index2 not in vertex_map:
            vertex_map[vert_index2] = len(vertex_map)
        if vert_index3 not in vertex_map:
            vertex_map[vert_index3] = len(vertex_map)
        vert_index1_new = vertex_map[vert_index1]
        vert_index2_new = vertex_map[vert_index2]
        vert_index3_new = vertex_map[vert_index3]

        # adding vertices
        vertices_new[vert_index1_new] = vertices[vert_index1]
        vertices_new[vert_index2_new] = vertices[vert_index2]
        vertices_new[vert_index3_new] = vertices[vert_index3]

        # adding triangles
        triangles_new[triangle_idx_new] = [vert_index1_new, vert_index2_new, vert_index3_new]

        # adding vertex colors
        vertex_colors_new[vert_index1_new] = selected_color
        vertex_colors_new[vert_index2_new] = selected_color
        vertex_colors_new[vert_index3_new] = selected_color

        # remove the triangle from the original mesh
        if remove_from_original:
            triangles[triangle_idx] = [-1, -1, -1]

    mesh_new.vertices = o3d.utility.Vector3dVector(vertices_new)
    mesh_new.triangles = o3d.utility.Vector3iVector(triangles_new)
    mesh_new.vertex_colors = o3d.utility.Vector3dVector(vertex_colors_new)
    mesh_new.compute_vertex_normals()
    
    return mesh_new

### edges
def create_cylinder_between_points(p1, p2, radius=0.01, resolution=20):
    cylinder = o3d.geometry.TriangleMesh.create_cylinder(radius, np.linalg.norm(p2 - p1), resolution)
    midpoint = (p1 + p2) / 2
    direction = (p2 - p1) / np.linalg.norm(p2 - p1)

    # Rotate the cylinder to align with the edge direction
    z_axis = np.array([0, 0, 1])
    axis = np.cross(z_axis, direction)
    angle = np.arccos(np.dot(z_axis, direction))

    if np.linalg.norm(axis) > 0:  # Avoid issues when direction is along the z-axis
        rotation_matrix = o3d.geometry.get_rotation_matrix_from_axis_angle(axis * angle)
        cylinder.rotate(rotation_matrix, center=(0, 0, 0))

    # Translate cylinder to the midpoint of p1 and p2
    cylinder.translate(midpoint - np.array([0, 0, np.linalg.norm(p2 - p1) / 2]))
    return cylinder

def make_edges(mesh, offset=0.01, edge_thickness=0.01):
    edges = set()  # Using a set to store unique edges

    # Each triangle has 3 edges: (v0, v1), (v1, v2), (v2, v0)
    for tri in mesh.triangles:
        edges.add(tuple(sorted([tri[0], tri[1]])))
        edges.add(tuple(sorted([tri[1], tri[2]])))
        edges.add(tuple(sorted([tri[2], tri[0]])))

    # Convert edges to a list and prepare for LineSet
    edges = list(edges)
    points = np.asarray(mesh.vertices)  # Get the vertices (points) of the mesh

    # Create a LineSet from points and edges
    line_set = o3d.geometry.LineSet()
    line_set.points = o3d.utility.Vector3dVector(points)
    line_set.lines = o3d.utility.Vector2iVector(edges)

    # Optionally, add colors to each line (same color for all lines in this example)
    line_set.colors = o3d.utility.Vector3dVector([[0, 0, 0] for _ in range(len(edges))])  # Blue edges
    
    # Cylinders
    cylinders = []
    for edge in edges:
        p1, p2 = points[edge[0]], points[edge[1]]
        p1 = (1.0 + offset) * p1  # Apply offset to slightly enlarge the mesh
        p2 = (1.0 + offset) * p2

        # Create a cylinder for each edge
        cylinder = create_cylinder_between_points(p1, p2, radius=edge_thickness)
        cylinders.append(cylinder)

    return line_set, cylinders

