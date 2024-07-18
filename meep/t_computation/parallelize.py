import os
import sys
from typing import Tuple
import meep as mp
import numpy as np
from mpi4py import MPI
from numpy.typing import NDArray

COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()


def start_parallel(
    n_sim: int, processes_per_group: int = 4
) -> Tuple[Tuple[int], int, range]:
    if not mp.am_master():
        sys.stdout = mp.saved_stdout
    n_groups = int(np.ceil(SIZE / processes_per_group))
    group_id = mp.divide_parallel_processes(n_groups)
    master_ids = tuple(mp.get_group_masters())
    sims_per_group = (1.0 * n_sim) / n_groups
    start_sim_id = int(group_id * sims_per_group)
    end_sim_id = (
        int((group_id + 1) * sims_per_group) if group_id < n_groups - 1 else n_sim
    )
    sim_id_range = range(start_sim_id, end_sim_id)
    return master_ids, group_id, sim_id_range


def end_parallel() -> None:
    mp.end_divide_parallel()
    if not mp.am_master():
        sys.stdout = open(os.devnull, "w")


def all_gather_matrix_concatenate(
    matrix: NDArray, master_ids: Tuple[int], group_id: int, axis=0
) -> NDArray:
    gathered_matrices = allgather_group_masters(
        data=matrix,
        master_ids=master_ids,
        group_id=group_id,
    )
    return np.concatenate([*gathered_matrices], axis=axis)

def allgather_group_masters(data: NDArray, master_ids: Tuple[int], group_id: int):
    x = COMM.allgather(data)
    return x
