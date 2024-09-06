"""
Create STL objects for building a transformations .glb file.
"""
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import os
import numpy as np
from stl import mesh


def create_flat_cuboid_arrays(
    scale: float = 1, length: float = 2, width: float = 1.3, height: float = 0.3
):
    """Get vertices and index arrays for creating a flat cuboid.

    Args:
        scale (float, optional): The scale of the cuboid. Defaults to 1.
        length (float, optional): The length of the cuboid. Defaults to 2.
        width (float, optional): The width of the cuboid. Defaults to 1.3.
        height (float, optional): The height of the cuboid. Defaults to 0.3.

    Returns:
        (np.ndarray, np.ndarray): The points and triangles array of the cuboid.
    """
    # Define the 8 vertices of the cuboid
    vertices = np.array(
        [
            [0, 0, 0],
            [length, 0, 0],
            [length, width, 0],
            [0, width, 0],
            [0, 0, height],
            [length, 0, height],
            [length, width, height],
            [0, width, height],
        ],
        dtype="float32",
    )

    # Define the 12 triangles (2 per cuboid face)
    indices = np.array(
        [
            [0, 1, 2],
            [0, 2, 3],  # Bottom face
            [4, 5, 6],
            [4, 6, 7],  # Top face
            [0, 1, 5],
            [0, 5, 4],  # Front face
            [1, 2, 6],
            [1, 6, 5],  # Right face
            [2, 3, 7],
            [2, 7, 6],  # Back face
            [3, 0, 4],
            [3, 4, 7],  # Left face
        ],
        dtype="uint32",
    )

    return indices, vertices * scale


import numpy as np


def create_inverted_frustum_cone_arrays(
    scale: float = 1,
    top_radius: float = 0.3,
    bottom_radius: float = 1.3,
    height: float = 2,
    num_sides: int = 50,
):
    """Get vertices and index arrays for creating an inverted frustum cone with a smaller opening at the bottom.

    Args:
        scale (float, optional): The scale of the cone. Defaults to 1.
        top_radius (float, optional): The radius of the top circular base. Defaults to 0.3.
        bottom_radius (float, optional): The radius of the bottom circular opening. Defaults to 1.3.
        height (float, optional): The height of the cone. Defaults to 2.
        num_sides (int, optional): The number of sides for the circular base and opening. Defaults to 50.

    Returns:
        (np.ndarray, np.ndarray): The points and triangles array of the cone.
    """
    theta = np.linspace(0, 2 * np.pi, num_sides, endpoint=False)

    # Bottom circle (smaller opening, now at the bottom)
    bottom_circle = np.stack(
        (
            top_radius * np.cos(theta),
            top_radius * np.sin(theta),
            np.zeros_like(theta),
        ),
        axis=-1,
    )
    # Top circle (larger base, now at the top)
    top_circle = np.stack(
        (
            bottom_radius * np.cos(theta),
            bottom_radius * np.sin(theta),
            np.full_like(theta, height),
        ),
        axis=-1,
    )

    # Vertices array
    vertices = np.vstack((bottom_circle, top_circle))

    # Indices array
    indices = []
    for i in range(num_sides):
        next_i = (i + 1) % num_sides
        # Side triangles
        indices.append([i, num_sides + next_i, num_sides + i])
        indices.append([i, next_i, num_sides + next_i])

    indices = np.array(indices, dtype="uint32")

    return indices, vertices * scale


def save_to_stl(indices, vertices, filename):
    """Save the mesh defined by indices and vertices to an STL file.

    Args:
        indices (np.ndarray): Array of triangle indices.
        vertices (np.ndarray): Array of vertices.
        filename (str): The filename of the output STL file.
    """
    # Create the mesh
    stl_mesh = mesh.Mesh(np.zeros(indices.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(indices):
        for j in range(3):
            stl_mesh.vectors[i][j] = vertices[f[j], :]

    # Save the mesh to file
    stl_mesh.save(filename)


if __name__ == "__main__":
    folder = os.path.dirname(__file__)
    shapes = {
        "sample": (create_flat_cuboid_arrays, (1, 2, 1.3, 0.3)),
        "beam": (create_flat_cuboid_arrays, (1, 10, 0.1, 0.2)),
        "analyser": (create_inverted_frustum_cone_arrays, (1, 0.3, 1.3, 2, 50)),
    }

    for shape, (shape_func, args) in shapes.items():
        indices, vertices = shape_func(*args)
        filename = os.path.join(folder, f"{shape}.stl")
        save_to_stl(indices, vertices, filename)
        print(f"STL file saved as '{shape}.stl'")
