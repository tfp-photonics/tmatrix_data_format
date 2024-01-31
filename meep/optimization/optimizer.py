from typing import Callable, Tuple

import jax.numpy as jnp
import mpi4jax
import numpy as np
from jax import value_and_grad
from jax.test_util import check_grads
from mpi4py import MPI
from numpy.typing import NDArray

from config.config import Config
from config.loss_func import get_default_loss_func
from geometry.geometry_data import GeometryData
from optimization.forward_pass import calc_T
from optimization.loss_data import LossData
from optimization.parallelize import end_parallel, start_parallel
from optimization.progress_data import ProgressData

COMM = MPI.COMM_WORLD


def calculate_T(c: Config) -> NDArray:
    geo_dat = get_start_geometry(c)
    eps_grid = jnp.asarray(geo_dat.eps_grid, dtype=float)
    master_ids, group_id, sim_id_range = start_parallel(
        n_sim=c.sim_amount, processes_per_group=c.cpu_cores_per_simulation
    )
    t_dat = calc_T(c, master_ids, group_id, sim_id_range, eps_grid)
    t_dat.to_hdf(path=c.path, it=0)
    end_parallel()
    return t_dat.t


def optimize(
    c: Config, goal_T: NDArray = None, loss_func_of_T: Callable[[NDArray], float] = None
) -> NDArray:
    if loss_func_of_T is None:
        check_goal_T_shape(c, goal_T)
    loss_func_of_T = get_loss_func_of_T_considering_goal_T(goal_T, loss_func_of_T)
    start_it = get_start_it(c)
    geo_dat = get_start_geometry(c, start_it)
    eps_grid = jnp.asarray(geo_dat.eps_grid, dtype=float)
    for it in range(start_it, c.opt_iterations):
        master_ids, group_id, sim_id_range = start_parallel(
            n_sim=c.sim_amount, processes_per_group=c.cpu_cores_per_simulation
        )
        loss_func = get_loss_func_of_eps(
            c=c,
            master_ids=master_ids,
            group_id=group_id,
            sim_id_range=sim_id_range,
            it=it,
            loss_func_of_T=loss_func_of_T,
        )
        # check_grads(loss_func, (eps_grid,), order=1, modes="rev")
        loss, loss_grad = value_and_grad(loss_func)(eps_grid)
        LossData(loss).to_hdf(path=c.path, it=it)
        gathered_loss_grad, _token = mpi4jax.allgather(loss_grad, comm=COMM)
        loss_grad = jnp.sum(gathered_loss_grad[master_ids, ...], axis=0)
        eps_grid.at[geo_dat.no_pad_mask].set(
            jnp.clip(
                c.opt_eps_min,
                c.opt_eps_max,
                eps_grid[geo_dat.no_pad_mask] - 0.05 * loss_grad[geo_dat.no_pad_mask],
            )
        )
        geo_dat.eps_grid = eps_grid
        end_parallel()
        geo_dat.to_hdf(path=c.path, it=it + 1)
        ProgressData(it + 1).to_hdf(path=c.path)
    return geo_dat.eps_grid[geo_dat.no_pad_mask]


def check_goal_T_shape(c: Config, goal_T: NDArray):
    if goal_T is not None:
        t_shape = c.t_shape
        goal_shape = goal_T.shape
        if np.any((np.array([*goal_shape]) - np.array([*t_shape]) < 0)):
            raise ValueError(
                f"Provided goal_T does not suffice in size. The optimizer produces matrices of shape {t_shape}. Provided goal_T has shape {goal_shape}."
            )


def get_loss_func_of_T_considering_goal_T(
    goal_T: NDArray = None, loss_func_of_T: Callable[[NDArray], float] = None
) -> Callable[[NDArray], float]:
    if loss_func_of_T is None:
        if goal_T is None:
            raise ValueError(
                "Either loss_func_of_T or goal_T must be provided to the optimizer."
            )
        else:
            loss_func_of_T = get_default_loss_func(goal_T)
    else:
        if goal_T is not None:
            print("loss_func_of_T was provided, overwriting the given goal_T.")
    return loss_func_of_T


def get_loss_func_of_eps(
    c: Config,
    master_ids: Tuple[int],
    group_id: int,
    sim_id_range: range,
    it: int,
    loss_func_of_T: Callable[[NDArray], float],
) -> float:
    def _loss_func(eps_grid: NDArray):
        t_dat = calc_T(c, master_ids, group_id, sim_id_range, eps_grid)
        t_dat.to_hdf(path=c.path, it=it)
        # plot_treams_vs_me(path=c.path, it=it, l_max=c.l_max, f_cen=c.fcen)
        # plot_jcm_vs_me(path=c.path, it=it)
        loss = loss_func_of_T(t_dat.t)
        return loss

    return _loss_func


def get_start_it(c: Config) -> int:
    # Loading not safe when using multiple processes. Cannot read simultaneously...
    if c.continue_opt_iterations:
        try:
            return ProgressData.from_hdf(path=c.path).progress
        except Exception as e:
            return 0
    else:
        return 0


def get_start_geometry(c: Config, start_it: int = 0) -> GeometryData:
    if start_it == 0:
        geometry_resolution = np.array([*c.start_geometry.shape]) / c.object_size
        n_pad = ((c.no_pml_size - c.object_size) * geometry_resolution / 2.0).astype(
            int
        )
        no_pad_mask = np.ones_like(c.start_geometry, dtype=bool)
        no_pad_mask = np.pad(
            no_pad_mask,
            ((n_pad[0],), (n_pad[1],), (n_pad[2],)),
            constant_values=False,
        )
        eps_grid = np.pad(
            c.start_geometry,
            ((n_pad[0],), (n_pad[1],), (n_pad[2],)),
            constant_values=c.eps_embedding,
        )
        geo_dat = GeometryData(eps_grid=eps_grid, no_pad_mask=no_pad_mask)
        geo_dat.to_hdf(c.path, 0)
    else:
        # Loading not safe when using multiple processes. Cannot read simultaneously...
        geo_dat = GeometryData.from_hdf(path=c.path, it=start_it)
    return geo_dat
