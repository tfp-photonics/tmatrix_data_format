from typing import Tuple

import jax.numpy as jnp
from numpy.typing import NDArray

from config.config import Config
from decomposition.coef_data import CoefData
from decomposition.t_data import TData
from decomposition.tmat_lstsq import calc_T_from_coefs
from decomposition.wave_decomposer import (
    get_e_amp_from_e_in,
    get_inc_coefs,
    get_sca_coefs,
)
from geometry.rotate import rotate_eps_grid
from geometry.rotation_factory import get_persistent_sphere_angles
from optimization.parallelize import all_gather_matrix_concatenate
from sim.sim import do_incident_sim, do_scattered_sim
from sources.source_factory import get_plane_wave_source


def calc_T(
    c: Config,
    master_ids: Tuple[int],
    group_id: int,
    sim_id_range: range,
    eps_grid: NDArray,
) -> TData:

    source_data = get_plane_wave_source(c)
    sim_res_inc = do_incident_sim(c=c, source_data=source_data)
    e_amp_inc = get_e_amp_from_e_in(c=c, source_data=source_data, sim_res=sim_res_inc)

    incident_matrix = []
    scatter_matrix = []
    for ii in sim_id_range:
        print(f"Simulation id = {ii}")
        rotation_data = get_persistent_sphere_angles(c=c, id=ii)
        rotated_epsgrid = rotate_eps_grid(eps_grid, rotation_data, c.eps_embedding)
        sim_res_sca = do_scattered_sim(
            id=ii,
            c=c,
            source_data=source_data,
            eps_grid=rotated_epsgrid,
            sim_res_inc=sim_res_inc,
        )
        incident_coefs = get_inc_coefs(
            c=c, source_data=source_data, rotation_data=rotation_data, e_amp=e_amp_inc
        )
        scatter_coefs = get_sca_coefs(
            c=c, rotation_data=rotation_data, sim_res=sim_res_sca
        )
        incident_matrix.append(incident_coefs)
        scatter_matrix.append(scatter_coefs)
    incident_matrix = jnp.stack(incident_matrix, axis=1)
    scatter_matrix = jnp.stack(scatter_matrix, axis=1)
    gathered_inc_matrix = all_gather_matrix_concatenate(
        matrix=incident_matrix,
        master_ids=master_ids,
        group_id=group_id,
        axis=1,
    )
    gathered_sca_matrix = all_gather_matrix_concatenate(
        matrix=scatter_matrix,
        master_ids=master_ids,
        group_id=group_id,
        axis=1,
    )
    coefs = CoefData(
        incident_matrix=gathered_inc_matrix, scatter_matrix=gathered_sca_matrix
    )
    coefs.to_hdf(c.path)
    return calc_T_from_coefs(coefs=coefs)
