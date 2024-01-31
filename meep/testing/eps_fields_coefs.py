import matplotlib.pyplot as plt
import meep as mp
import numpy as np
import treams
from numpy.typing import ArrayLike, NDArray

from config.config import Config
from geometry.geometry_factory import (
    get_epsgrid_from_geometric_objects,
    get_geometry_from_eps_grid,
)
from geometry.rotation_data import RotationData
from sim.sim_data import SimData
from sources.source_data import SourceData
from testing.treams_t import get_t_treams
from utility.plot_lib import add_im_to_plot, plot_complex_versus
from typing import List

def eps_grid_meep_object_comparison():
    res = 30
    no_pml_size = np.array([2.0, 2.5, 3.0])
    meep_geometry = [
        mp.Sphere(
            0.2,
            center=mp.Vector3(),
            material=mp.Medium(epsilon=1.15),
        ),
        mp.Cylinder(
            center=mp.Vector3(0.2, 0.2, 0),
            axis=mp.Vector3(0.5, 2.0, 1.0),
            radius=0.2,
            height=1,
            material=mp.Medium(epsilon=1.25),
        ),
        mp.Block(
            center=mp.Vector3(-0.5, -0.2),
            size=mp.Vector3(0.5, 0.5, 0.5),
            e1=mp.Vector3(2.0, 1.0, 0.0),
            e2=mp.Vector3(0, 1.0, 0.0),
            e3=mp.Vector3(0, 0, 1.0),
            material=mp.Medium(epsilon=1.35),
        ),
    ]

    def get_eps(geometry):
        sim = mp.Simulation(
            cell_size=mp.Vector3(4, 4, 4),
            resolution=res,
            geometry=geometry,
            sources=[],
            boundary_layers=[mp.PML(1)],
            k_point=mp.Vector3(),
        )
        sim.init_sim()
        return sim.get_epsilon()

    eps1 = get_eps(meep_geometry)
    epsgrid = get_epsgrid_from_geometric_objects(
        res=2 * res, size=no_pml_size, meep_objects=meep_geometry
    )
    # plot_im(epsgrid[..., int(epsgrid.shape[2] / 2)])
    grid_geometry = get_geometry_from_eps_grid(eps=epsgrid, size=no_pml_size)
    eps2 = get_eps(grid_geometry)
    plot_complex_versus(
        eps1[..., int(eps1.shape[2] / 2)].T,
        eps2[..., int(eps2.shape[2] / 2)].T,
        "meep_object",
        "eps_grid",
    )
    plt.show()


def plot_field(sim_res: SimData, comps:List[int], sides:List[int], filename ='incident.png', xlabel = '', ylabel='') -> None:
    if mp.am_really_master():
        f_index = int(sim_res.e.shape[0] / 2)
        fig, axes = plt.subplots(len(comps), len(sides), figsize=(len(sides)*4, len(comps)*4))
        axes = np.asarray(axes)
        if (len(comps) == 1):
            axes = axes[None,...]
        if (len(sides) == 1):
            axes = axes[...,None]
        for ii,side in enumerate(sides):
            field = sim_res.e[f_index, side, :, :, :]
            side_str = "x-"
            if side == 1:
                side_str = "x+"
            elif side == 2:
                side_str = "y-"
            elif side == 3:
                side_str = "y+"
            elif side == 4:
                side_str = "z-"
            elif side == 5:
                side_str = "z+"
            for jj,comp in enumerate(comps):
                comp_str = "x"
                if comp == 1:
                    comp_str = "y"
                elif comp == 2:
                    comp_str = "z"
                ax = axes[jj,ii]
                field = 1.02*field/np.max(np.abs(field[...,comp]))
                add_im_to_plot(
                    fig,
                    ax,
                    np.real(field[..., comp]),
                    f"Slice of $E_{comp_str}$ on cuboid surface {side_str}",
                    cbar_ticks=[-1,0,1]
                )
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.set_xticks([-0.5,field.shape[0]/2,field.shape[0]-0.5],[-1,0,1])
                ax.set_yticks([-0.5,field.shape[1]/2,field.shape[1]-0.5],[-1,0,1])
        plt.tight_layout()
        plt.savefig(filename, dpi=400)
        plt.show(block=False)


def plot_fitted_incident_field(
    source_data: SourceData,
    sim_res: SimData,
    e_abs_amp: ArrayLike,
    e_phase: ArrayLike,
) -> None:
    if mp.am_really_master():
        f_index = int(e_abs_amp.shape[0] / 2.0)
        fig, axes = plt.subplots(3, 6, figsize=(18, 10))
        for ii in range(6):
            for jj in range(3):
                amp = e_abs_amp[f_index]
                phase = e_phase[f_index]
                field = np.real(
                    source_data.pol_src[jj]
                    * amp
                    * np.exp(
                        1j
                        * (np.dot(sim_res.pos[ii, :, :, :], source_data.k_src) - phase)
                    )
                )
                add_im_to_plot(
                    fig,
                    axes[jj, ii],
                    field,
                    f"E_{jj} cube side {ii}",
                )
        plt.show(block=True)


def plot_scatter_coefs_treams_vs_me(
    frequencies: NDArray,
    inc_coefs: NDArray,
    sca_coefs: NDArray,
    l_max=15,
) -> None:
    treams_sca = np.zeros_like(sca_coefs)
    for ii, f in enumerate(frequencies):
        t = get_t_treams(l_max=l_max, f=f)
        treams_sca[ii] = (t @ inc_coefs[ii])[: sca_coefs.shape[1]]
    treams_sca = treams_sca.T
    sca = sca_coefs.T
    plot_complex_versus(
        treams_sca, sca, title_1="treams sca_coef", title_2="my sca_coef"
    )
    plt.show()


def plot_field_treams_vs_me(
    c: Config,
    sim_res_inc: SimData,
    sim_res_sca: SimData,
    rotation_data: RotationData,
    source_data: SourceData,
    e_amp: NDArray,
    l_max=15,
) -> None:
    freq_index = int(c.frequencies.shape[0] / 2.0)
    f = c.frequencies[freq_index]
    t_treams = get_t_treams(l_max=l_max, f=f)
    k_abs = c.wave_numbers[freq_index]
    k_pw_rot = rotation_data.rotation_matrix @ (
        k_abs * source_data.k_src / np.linalg.norm(source_data.k_src)
    )
    e_amp_rot = (rotation_data.rotation_matrix @ e_amp.T).T
    pw = treams.plane_wave(
        k_pw_rot,
        [*e_amp_rot[freq_index]],
        k0=k_abs,
        material=treams.Material(),
        poltype="parity",
    )
    xyz_component = 0
    side = 0

    inc = pw.expand(t_treams.basis)
    sca = t_treams @ pw.expand(t_treams.basis)

    def _plot_field_versus(treams_object, sim_res: SimData):
        grid = sim_res.pos[side]
        field_treams = (treams_object.efield(grid) * treams_object[:, None]).sum(
            axis=-2
        )
        field_treams = field_treams[:, :, xyz_component]
        field_me = sim_res.e[freq_index, side, :, :, xyz_component]
        plot_complex_versus(
            field_treams, field_me, title_1="treams field", title_2="my field"
        )
        plt.show()

    _plot_field_versus(inc, sim_res_inc)
    _plot_field_versus(sca, sim_res_sca)
