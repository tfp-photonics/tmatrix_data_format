import warnings
from typing import Optional

import meep as mp
import numpy as np
from numpy.typing import NDArray


def cartesian_to_spherical(xyz: NDArray):
    x = xyz[..., 0]
    y = xyz[..., 1]
    z = xyz[..., 2]
    r = np.sqrt(x * x + y * y + z * z)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return np.stack([r, theta, phi], axis=-1)


def spherical_to_catesian(rtp: NDArray):
    r = rtp[..., 0]
    theta = rtp[..., 1]
    phi = rtp[..., 2]
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.stack([x, y, z], axis=-1)


def rotation_matrix(axis: NDArray, theta: float) -> NDArray:
    """
    Return the rotation matrix associated with counterclockwise rotation around
    the given axis by theta in radians.
    """
    axis = axis / np.linalg.norm(axis)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )


def meep_field_component_to_index(meep_component) -> Optional[int]:
    if (
        meep_component == mp.Ex
        or meep_component == mp.Bx
        or meep_component == mp.Dx
        or meep_component == mp.Hx
    ):
        return 0
    elif (
        meep_component == mp.Ey
        or meep_component == mp.By
        or meep_component == mp.Dy
        or meep_component == mp.Hy
    ):
        return 1
    elif (
        meep_component == mp.Ez
        or meep_component == mp.Bz
        or meep_component == mp.Dz
        or meep_component == mp.Hz
    ):
        return 2
    else:
        warnings.warn("Meep field component not recognized.")
        return None


def index_to_meep_e_field_component(index) -> Optional[int]:
    if index == 0:
        return mp.Ex
    elif index == 1:
        return mp.Ey
    elif index == 2:
        return mp.Ez
    else:
        warnings.warn("Meep E-field index not recognized.")
        return None
