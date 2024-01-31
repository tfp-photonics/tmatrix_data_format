import jax.numpy as jnp
import numpy as np
import treams
from numpy.typing import NDArray

from config.config import Config


def get_diag_goal_T(c: Config):
    t_shape = c.t_shape
    goal_T = np.zeros(t_shape)
    mask = np.eye(t_shape[-1], dtype=bool)
    goal_T[:, mask] = 1
    return jnp.asarray(goal_T)


def get_sphere_goal_T(c: Config):
    goal_T = []
    for f in c.frequencies:
        goal_T.append(
            np.asarray(
                treams.TMatrix.sphere(
                    c.l_max,
                    f * 2 * np.pi,
                    0.4,
                    [treams.Material(1.2), treams.Material()],
                )
            )
        )
    return np.stack(goal_T)


def get_goal_T(c: Config) -> NDArray:
    return get_sphere_goal_T(c=c)
