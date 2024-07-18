import h5py
import matplotlib.pyplot as plt
import meep as mp
import numpy as np
import treams
from numpy.typing import NDArray

from decomposition.t_data import TData
from testing.treams_t import get_t_treams
from utility.plot_lib import plot_anim, plot_complex_versus


def plot_t_Anim(t, path: str) -> None:
    plot_anim(
        t,
        path=f"{path}/t.gif",
        cmap="Purples",
        title="Absolute values of the T-matrix for different frequencies",
    )
    plt.show()


def plot_me_vs_treams(
    t_me: NDArray, t_treams: treams.PhysicsArray, path: str, name: str = "tmat"
) -> None:
    if mp.am_really_master():
        f_index = int(t_me.shape[0] / 2.0)
        t_me = t_me[f_index]
        size = t_me.shape[0]
        t_treams = t_treams[:size, :size]
        plot_complex_versus(
            t_treams,
            t_me,
            title_1="$T^{\mathrm{treams}}$",
            title_2="$T^{\mathrm{meep}}$",
        )
        plt.tight_layout()
        print(path, name)
        plt.savefig(f"{path}/{name}.png", dpi=400)
        # plotTAnim(t, path)


def plot_jcm_vs_me(path: str, it: int = 0) -> None:
    if mp.am_really_master():
        t_data = TData.from_hdf(path, it).t
        f_index = int(t_data.shape[0] / 2.0)
        t_me = t_data[f_index]
        with h5py.File("../example_data/cylinder_tio2.h5", "r") as hf:
            t_jcm = np.array(hf.get("tmatrix"))[0]
        size = t_me.shape[0]
        t_jcm = t_jcm[:size, :size]
        plot_complex_versus(t_jcm, t_me, title_1="jcm T", title_2="meep T")
        plt.show()
        # plotTAnim(t, path)
