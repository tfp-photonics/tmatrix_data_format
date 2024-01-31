import meep as mp
import numpy as np
from numpy.typing import NDArray

from utility.drawing import (
    draw_box_in_array,
    draw_cylinder_in_array,
    draw_sphere_in_array,
)


def get_epsgrid_from_geometric_objects(
    res: int,
    size: NDArray,
    meep_objects: list[mp.GeometricObject],
    eps_background=1.0,
) -> NDArray:

    eps = np.full(tuple(map(int, tuple(res * size))), eps_background)
    for object in meep_objects:
        if type(object) == mp.Sphere:
            eps = draw_sphere_in_array(
                arr=eps,
                size=size,
                position=np.array([*object.center]),
                radius=object.radius,
                value=object.material.epsilon_diag[0],
            )
        elif type(object) == mp.Cylinder:
            eps = draw_cylinder_in_array(
                arr=eps,
                size=size,
                position=np.array([*object.center]),
                normal=np.array([*object.axis]),
                radius=object.radius,
                height=object.height,
                value=object.material.epsilon_diag[0],
            )
        elif type(object) == mp.Block:
            eps = draw_box_in_array(
                arr=eps,
                size=size,
                position=np.array([*object.center]),
                edge_1=np.array([*object.e1]),
                edge_2=np.array([*object.e2]),
                edge_3=np.array([*object.e3]),
                edge_lenths=np.array([*object.size]),
                value=object.material.epsilon_diag[0],
            )
    return eps


def get_geometry_from_eps_grid(
    eps: NDArray,
    size: NDArray,
    center=np.zeros((3,)),
) -> list[mp.GeometricObject]:
    eps_matgrid: mp.MaterialGrid = get_matgrid_from_epsgrid(eps)
    geometry = [
        mp.Block(
            center=mp.Vector3(*center),
            size=mp.Vector3(*size),
            material=eps_matgrid,
        ),
    ]
    return geometry


def get_matgrid_from_epsgrid(eps: NDArray) -> mp.MaterialGrid:
    eps_min = np.min(eps)
    eps_max = np.max(eps)
    weights = (
        (eps - eps_min) / (eps_max - eps_min)
        if eps_min != eps_max
        else np.zeros_like(eps)
    )
    matgrid = mp.MaterialGrid(
        mp.Vector3(*eps.shape),
        mp.Medium(epsilon=eps_min),
        mp.Medium(epsilon=eps_max),
        do_averaging=False,
    )
    matgrid.update_weights(weights)
    return matgrid
