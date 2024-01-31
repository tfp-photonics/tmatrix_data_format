import numpy as np
import treams
from numpy.typing import NDArray

from config.loss_func import compare_matrices


def delta_treams_vs_me(t_me: NDArray, t_treams: treams.PhysicsArray) -> None:
    f_index = int(t_me.shape[0] / 2.0)
    t_me = t_me[f_index]
    size = t_me.shape[0]
    t_treams = t_treams[:size, :size]
    t_treams = np.asarray(t_treams)
    delta = compare_matrices(t_me, t_treams)
    print("delta")
    print(delta)
    return delta
