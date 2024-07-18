import numpy as np
from numpy.typing import NDArray

from decomposition.coef_data import CoefData
from decomposition.t_data import TData


def solve_lstsq(a: NDArray, b: NDArray):
    return np.linalg.lstsq(a, b, rcond=-1)[0]


def calc_T_from_coefs(coefs: CoefData) -> TData:
    nfreq = coefs.scatter_matrix.shape[0]
    ts = np.stack(
        [
            solve_lstsq(coefs.incident_matrix[ii], coefs.scatter_matrix[ii]).T
            for ii in range(nfreq)
        ]
    )
    return TData(t=ts)
