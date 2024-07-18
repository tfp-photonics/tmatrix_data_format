from typing import Callable

import numpy as np
from numpy.typing import NDArray


def get_default_loss_func(goal_T: NDArray) -> Callable[[NDArray], float]:
    def default_loss_func(t: NDArray):
        shape = t.shape
        loss = 0
        for ii in range(shape[0]):
            loss += compare_matrices(t[ii, ...], goal_T[ii, : shape[1], : shape[2]])
        return loss

    return default_loss_func


def compare_matrices(A, B):
    c = (
        np.trace((A - B).conj().T @ (A - B))
        / (np.trace(A.conj().T @ A) + np.trace(B.conj().T @ B))
        / 2
    )
    return c.real
