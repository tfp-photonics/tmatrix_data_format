import os
import sys
from typing import Tuple

import jax.numpy as jnp
import meep as mp
import mpi4jax
import numpy as np
from jax import custom_vjp, value_and_grad
from jax.test_util import check_grads
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


def allgather_group_masters(data: NDArray, master_ids: Tuple[int], group_id: int):
    @custom_vjp
    def jax_allgather(x: NDArray):
        x, _token = mpi4jax.allgather(x, comm=COMM)
        print(f"gathered_x {RANK}", x)
        return x[master_ids, ...]

    def allgather_fwd(x: NDArray):
        return jax_allgather(x), None

    def allgather_bwd(res, grads):
        print(f"bwd {RANK}", grads[group_id, ...])
        return (grads[group_id, ...],)

    jax_allgather.defvjp(allgather_fwd, allgather_bwd)
    return jax_allgather(data)


def test_mpi4jax():
    # works only for one process per group
    master_ids, group_id, range = start_parallel(99, 2)

    def loss(x):
        x = allgather_group_masters(x, master_ids, group_id)
        return jnp.sum(x * x)

    np.random.seed(group_id)
    rand = np.random.random(3)
    x = jnp.array(rand)
    print(f"x {RANK}", x)
    check_grads(loss, (x,), order=1, modes="rev")
    f, grad = value_and_grad(loss)(x)
    print(f"val_grad {RANK}", f, grad)
    end_parallel()


test_mpi4jax()
