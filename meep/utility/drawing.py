import numpy as np
from numpy.typing import NDArray


def _get_pos_grid(arr: NDArray, size: NDArray):
    x, y, z = [
        np.linspace(-size[ii] / 2.0, size[ii] / 2.0, arr.shape[ii]) for ii in range(3)
    ]
    xv, yv, zv = np.meshgrid(x, y, z, indexing="ij")
    return np.stack([xv, yv, zv], axis=-1)


def draw_sphere_in_array(
    arr: NDArray, size: NDArray, position: NDArray, radius: float, value: float
) -> NDArray:
    grid = _get_pos_grid(arr=arr, size=size)
    
    mask = np.sum((grid - position) ** 2, axis=-1) <= radius**2
    arr[mask] = value
    return arr


def draw_cylinder_in_array(
    arr: NDArray,
    size: NDArray,
    position: NDArray,
    normal: NDArray,
    radius: float,
    height: float,
    value: float,
) -> NDArray:
    
    normal = normal / np.linalg.norm(normal)
    grid = _get_pos_grid(arr=arr, size=size)
    normal_component = np.sum((grid - position) * normal, axis=-1)
    height_mask = np.abs(normal_component) <= height / 2.0
    grid_radial = grid - (normal_component[..., None] * normal[None, None, None, :])

    radial_mask = np.sum((grid_radial - position) ** 2, axis=-1) <= radius**2
    print("radial mask", radial_mask)
    print("height mask", height_mask)
    arr[height_mask & radial_mask] = value
    return arr


def draw_box_in_array(
    arr: NDArray,
    size: NDArray,
    position: NDArray,
    edge_1: NDArray,
    edge_2: NDArray,
    edge_3: NDArray,
    edge_lenths: NDArray,
    value: float,
) -> NDArray:
    grid = _get_pos_grid(arr=arr, size=size)
    edge_1 = edge_1 / np.linalg.norm(edge_1)
    edge_2 = edge_2 / np.linalg.norm(edge_2)
    edge_3 = edge_3 / np.linalg.norm(edge_3)
    edges = [edge_1, edge_2, edge_3]
    normals = [
        np.cross(edge_2, edge_3),
        np.cross(edge_3, edge_1),
        np.cross(edge_1, edge_2),
    ]
    lengths = np.array(
        [np.dot(edges[ii], normals[ii]) * edge_lenths[ii] for ii in range(3)]
    )
    grid = np.stack(
        [np.abs(np.sum((grid - position) * normal, axis=-1)) for normal in normals],
        axis=-1,
    )
    masks = grid <= lengths / 2.0
    arr[np.all(masks, axis=-1)] = value
    return arr
