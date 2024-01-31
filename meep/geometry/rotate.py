from typing import Callable, Tuple

import jax.numpy as jnp
from jax import custom_vjp
from numpy.typing import NDArray
from scipy.ndimage import rotate

from geometry.rotation_data import RotationData


def rotate_eps_grid(
    eps: NDArray, rotation_data: RotationData = RotationData(), eps_embedding=1.0
) -> NDArray:
    return _get_jax_rotate_eps_grid(rotation_data, eps_embedding)(eps)


def _get_jax_rotate_eps_grid(
    rotation_data: RotationData = RotationData(), eps_embedding=1.0
) -> Callable[[NDArray], NDArray]:
    def _jax_rotate_eps_grid(eps: NDArray):
        eps = jnp.asarray(eps)
        angles_and_planes = (
            rotation_data.rotation_as_list_of_angles_and_rotation_planes_inverse
        )
        for anp in angles_and_planes:
            eps = _get_jax_ndimage_rotate(
                angle_in_rad=anp[0], axes=anp[1], cval=eps_embedding
            )(eps)
        return eps

    return _jax_rotate_eps_grid


def _get_jax_ndimage_rotate(
    angle_in_rad: float, axes: Tuple[int], cval=0.0
) -> Callable[[NDArray], NDArray]:
    angle_in_degrees = angle_in_rad * 180.0 / jnp.pi

    @custom_vjp
    def jax_rotate(im: NDArray):
        return rotate(
            im, angle=angle_in_degrees, axes=axes, reshape=False, cval=cval, order=0
        )

    def rotate_fwd(im: NDArray):
        return jax_rotate(im), None

    def rotate_bwd(res, grads):
        return (
            jnp.asarray(
                rotate(
                    grads,
                    angle=-angle_in_degrees,
                    axes=axes,
                    reshape=False,
                    cval=0.0,
                    order=0,
                )
            ),
        )

    jax_rotate.defvjp(rotate_fwd, rotate_bwd)
    return jax_rotate
