from typing import Callable, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.ndimage import rotate

from geometry.rotation_data import RotationData


def rotate_eps_grid(
    eps: NDArray, rotation_data: RotationData = RotationData(), eps_embedding=1.0
) -> NDArray:
    return _get_rotate_eps_grid(rotation_data, eps_embedding)(eps)


def _get_rotate_eps_grid(
    rotation_data: RotationData = RotationData(), eps_embedding=1.0
) -> Callable[[NDArray], NDArray]:
    def _rotate_eps_grid(eps: NDArray):
        eps = np.asarray(eps)
        angles_and_planes = (
            rotation_data.rotation_as_list_of_angles_and_rotation_planes_inverse
        )
        for anp in angles_and_planes:
            angle_in_rad=anp[0]
            axes=anp[1]
            angle_in_degrees = angle_in_rad * 180.0 / np.pi
            eps = rotate(
            eps, angle=angle_in_degrees, axes=axes, reshape=False, cval=eps_embedding, order=0
            )
        return eps

    return _rotate_eps_grid

