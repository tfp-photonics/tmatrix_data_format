import sys

import jax.numpy as jnp
import meep as mp
import numpy as np
from jax import value_and_grad
from jax.test_util import check_grads
from mpi4py import MPI

from optimization.parallelize import (
    allgather_group_masters,
    end_parallel,
    start_parallel,
)

COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()


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
