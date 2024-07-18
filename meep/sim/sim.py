from typing import List

import matplotlib.pyplot as plt
import meep as mp
import meep.adjoint as mpa
import numpy as np
from numpy.typing import NDArray
from geometry.geometry_factory import (
    get_geometry_from_eps_grid,
    get_matgrid_from_epsgrid,
)
from config.config import Config
from sim.sim_data import SimData
from sources.source_data import SourceData
from utility.plot_lib import plot_im

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
    # print("eps_grid", eps_grid)
    cell = mp.Vector3(*c.with_pml_size)
    pml = [mp.PML(c.dpml)]
    sim = mp.Simulation(
        cell_size=cell,
        resolution=c.resolution,
        boundary_layers=pml,
        k_point=mp.Vector3(),
        sources=source_data.meep_sources,
        default_material=mp.Medium(epsilon=c.eps_embedding)
    )

    if check_setup:
        _check_setup(c=c, sim=sim)

    cube_length = np.amin(c.no_pml_size)
    # monitors, normals = _get_fourier_monitors(c, sim, cube_length)
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
            monitors.append(sim.add_dft_fields([comp], c.f_src, c.df_src, c.n_freq, center=center, size=size
            ))
    sim.geometry = get_geometry_from_eps_grid(
            eps=eps_grid,
            center=mp.Volume(size=mp.Vector3(*c.no_pml_size)).center,
            size=mp.Volume(size=mp.Vector3(*c.no_pml_size)).size,
        )
    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol= c.tol))
    results = []
    # Get the Fourier-transformed fields
    for i, normal in enumerate(normal_vectors):
        for j, comp in enumerate(components): 
            for k in range(c.n_freq):
                ind = i*len(components) +j
                center = mp.Vector3(*((cube_length / 2.0) * normal))
                size = mp.Vector3(*(cube_length * np.logical_not(np.abs(normal))))
                results.append( sim.get_dft_array(monitors[ind], comp, k))

    n_sides = 6
    results = np.array(results).reshape((n_sides * len(components), c.n_freq, np.array(results).shape[-1],  np.array(results).shape[-1] ))
    fts_E = []
    fts_B = []
    positions = []
    weights = []
    for side in range(n_sides):
        side_values = results[side * n_sides : (side + 1) * n_sides]
        fts_E.append(np.moveaxis(side_values[:3, ...], 0, -1))
        fts_B.append(np.moveaxis(side_values[3:, ...], 0, -1))
        a_monitor_of_side = monitors[side * n_sides]
        center = mp.Vector3(*((cube_length / 2.0) * normal_vectors[side]))
        size = mp.Vector3(*(cube_length * np.logical_not(np.abs(normal_vectors[side]))))
        pos, w = sim.get_array_metadata(
            center=center,
            size=size,
            return_pw=True,
        )
        positions.append(np.array(pos).reshape(np.asarray(w).shape + (3,)))
        weights.append(np.asarray(w))
    fts_E = np.stack(np.array(fts_E), axis=1)
    fts_B = np.stack(np.array(fts_B), axis=1)
    positions = np.stack(positions)
    weights = np.stack(weights)

    # Run forward simulation

    if check_setup:
        _plot_eps(c=c, sim=sim)

    result = SimData(
        e=fts_E,
        b=fts_B,
        normals=normal_vectors,
        pos=positions,
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
    print(sim_res)
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
