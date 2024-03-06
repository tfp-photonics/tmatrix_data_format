from typing import List

import jax.numpy as jnp
import matplotlib.pyplot as plt
import meep as mp
import meep.adjoint as mpa
import numpy as np
from numpy.typing import NDArray

from config.config import Config
from t_computation.sim_jax_wrapper import SimJaxWrapper
from sim.sim_data import SimData
from sources.source_data import SourceData
from utility.plot_lib import plot_im


def _reshape_and_complete_monitor_values(
    monitor_values: NDArray, monitors: List[mpa.FourierFields]
) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    # monitor_values.shape = (cubesides and field components like [(side 0, Ex), (side 0, Ey), ...], frequency, grid index 1, grid index 2)
    n_sides = 6
    fts_E = []
    fts_B = []
    positions = []
    weights = []
    for side in range(n_sides):
        side_values = monitor_values[side * n_sides : (side + 1) * n_sides]
        fts_E.append(jnp.moveaxis(side_values[:3, ...], 0, -1))
        fts_B.append(jnp.moveaxis(side_values[3:, ...], 0, -1))
        a_monitor_of_side = monitors[side * n_sides]
        pos, w = a_monitor_of_side.sim.get_array_metadata(
            center=a_monitor_of_side.volume.center,
            size=a_monitor_of_side.volume.size,
            return_pw=True,
        )
        positions.append(np.array(pos).reshape(w.shape + (3,)))
        weights.append(np.asarray(w))
    fts_E = jnp.stack(fts_E, axis=1)
    fts_B = jnp.stack(fts_B, axis=1)
    positions = np.stack(positions)
    weights = np.stack(weights)
    return fts_E, fts_B, positions, weights


def _get_fourier_monitors(
    sim: mp.Simulation, cube_length: float
) -> tuple[list[mpa.ObjectiveQuantity], NDArray]:
    monitors = []
    components = [mp.Ex, mp.Ey, mp.Ez, mp.Bx, mp.By, mp.Bz]
    normal_vectors = np.array(
        [
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 0.0, 1.0],
        ]
    )
    for normal in normal_vectors:
        for comp in components:
            center = mp.Vector3(*((cube_length / 2.0) * normal))
            size = mp.Vector3(*(cube_length * np.logical_not(np.abs(normal))))
            monitors.append(
                mpa.FourierFields(
                    sim=sim,
                    component=comp,
                    volume=mp.Volume(center=center, size=size),
                )
            )

    return monitors, normal_vectors


def _check_setup(c: Config, sim: mp.Simulation):
    n_plots = 20
    zs = np.linspace(-c.no_pml_size[2] / 2.0, c.no_pml_size[2] / 2.0, num=n_plots)
    for z in zs:
        sim.plot2D(
            output_plane=mp.Volume(
                center=(0, 0, z),
                size=mp.Vector3(c.no_pml_size[0], c.no_pml_size[1], 0),
            )
        )
        plt.savefig(f"setup_z={z}.png")


def _plot_eps(c: Config, sim: mp.Simulation):
    eps = sim.get_epsilon()
    for z_index in range(eps.shape[2]):
        plot_im(eps[..., z_index], "eps", centered=False)
        plt.savefig(f"eps_{z_index}.png")


def _do_sim(
    c: Config,
    source_data: SourceData,
    eps_grid: NDArray,
    check_setup=False,
) -> SimData:
    cell = mp.Vector3(*c.with_pml_size)
    pml = [mp.PML(c.dpml)]
    sim = mp.Simulation(
        cell_size=cell,
        resolution=c.resolution,
        boundary_layers=pml,
        k_point=mp.Vector3(),
        sources=source_data.meep_sources,
        default_material=mp.Medium(epsilon=c.eps_embedding),
    )

    if check_setup:
        _check_setup(c=c, sim=sim)

    cube_length = np.amin(c.no_pml_size)
    monitors, normals = _get_fourier_monitors(sim, cube_length)

    wrapped_meep = SimJaxWrapper(
        simulation=sim,
        design_volume=mp.Volume(size=mp.Vector3(*c.no_pml_size)),
        monitors=monitors,
        frequencies=c.frequencies,
    )

    monitor_values = wrapped_meep.simulate(eps_grid)

    if check_setup:
        _plot_eps(c=c, sim=sim)

    fts_E, fts_B, pos, weights = _reshape_and_complete_monitor_values(
        monitor_values, monitors
    )
    result = SimData(
        e=fts_E,
        b=fts_B,
        normals=normals,
        pos=pos,
        weights=weights,
    )
    return result


def do_incident_sim(
    c: Config,
    source_data: SourceData,
) -> SimData:
    sim_res = SimData.from_hdf_incident(
        path=c.path, identifier=c.scatter_sim_params_hash
    )
    if sim_res is None:
        print("Incident field simulation could not be loaded. Starting calculation.")
        sim_res = _do_sim(
            c=c, source_data=source_data, eps_grid=np.ones_like(c.start_geometry)
        )
        sim_res.to_hdf_incident(path=c.path, identifier=c.scatter_sim_params_hash)
    else:
        print("Incident field simulation successfully loaded.")
    return sim_res


def do_scattered_sim(
    id: int,
    c: Config,
    source_data: SourceData,
    eps_grid: NDArray,
    sim_res_inc: SimData,
) -> SimData:

    total_sim_res = None
    if c.load_simulations:
        total_sim_res = SimData.from_hdf_scattered(path=c.path, id=id)
        if total_sim_res is None:
            print(
                f"Scattered field simulation {id} could not be loaded. Starting calculation."
            )
        else:
            print(f"Scattered field simulation {id} successfully loaded.")
    if total_sim_res is None:
        total_sim_res = _do_sim(c=c, source_data=source_data, eps_grid=eps_grid)
        total_sim_res.to_hdf_scattered(path=c.path, id=id)
    sca_sim_res = SimData(
        e=total_sim_res.e - sim_res_inc.e,
        b=total_sim_res.b - sim_res_inc.b,
        normals=total_sim_res.normals,
        pos=total_sim_res.pos,
        weights=total_sim_res.weights,
    )
    return sca_sim_res
