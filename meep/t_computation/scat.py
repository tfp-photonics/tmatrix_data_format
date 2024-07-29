from typing import Tuple
from numpy.typing import NDArray
import numpy as np
from config.config import Config
from decomposition.coef_data import CoefData
from decomposition.t_data import TData
from decomposition.tmat_lstsq import calc_T_from_coefs
from decomposition.wave_decomposer import (
    get_e_amp_from_e_in,
    get_inc_coefs,
    get_sca_coefs,
)
from t_computation.parallelize import end_parallel, start_parallel
from geometry.geometry_data import GeometryData
from geometry.rotate import rotate_eps_grid
from geometry.rotation_factory import get_persistent_sphere_angles
from t_computation.parallelize import all_gather_matrix_concatenate
from sim.sim import do_incident_sim, do_scattered_sim
from sources.source_factory import get_plane_wave_source


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


def calculate_T(c: Config) -> NDArray:
    geo_dat = get_start_geometry(c)
    eps_grid = np.asarray(geo_dat.eps_grid, dtype=float)
    master_ids, group_id, sim_id_range = start_parallel(
        n_sim=c.sim_amount, processes_per_group=c.cpu_cores_per_simulation
    )
    t_dat = core_calc_T(c, master_ids, group_id, sim_id_range, eps_grid)
    t_dat.to_hdf(path=c.path, it=0)
    end_parallel()
    return t_dat.t


def core_calc_T(
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
    rot_matrix = []
    for ii in sim_id_range:
        print(f"Simulation id = {ii}")
        rotation_data = get_persistent_sphere_angles(c=c, id=ii)
        rotated_epsgrid = rotate_eps_grid(eps_grid, rotation_data, c.eps_embedding)
        rot_array = np.array(
            [rotation_data.theta, rotation_data.phi, rotation_data.alpha]
        )
        rot_matrix.append(rot_array)
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
    incident_matrix = np.stack(incident_matrix, axis=1)
    scatter_matrix = np.stack(scatter_matrix, axis=1)
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
    tm = calc_T_from_coefs(coefs=coefs)
    coefs.to_hdf(c.path)
    return tm
