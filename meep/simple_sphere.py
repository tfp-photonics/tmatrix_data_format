import jax
import meep as mp
import numpy as np
import treams
import os
from config.config import Config
from decomposition.t_data import TData
from optimization.optimizer import calculate_T
from testing.delta import delta_treams_vs_me
from testing.t_plots import plot_me_vs_treams

jax.config.update("jax_enable_x64", True)

if __name__ == "__main__":

    # path = os.getcwd()
    path = "./"
    geometry = [
        mp.Sphere(0.4, center=mp.Vector3(), material=mp.Medium(epsilon=1.2)),
    ]
    c = Config(
        cpu_cores_per_simulation=8,
        resolution=30,
        sim_amount_mult=2,
        load_simulations=True,
        l_max=2,
        object_size=np.array([1.0, 1.0, 1.0]),
        path=path,
        start_geometry=geometry,
    )
    t = calculate_T(c)
    # t = TData.from_hdf(path, 0).t
    t_treams = treams.TMatrix.sphere(
        c.l_max, c.f_cen * 2 * np.pi, 0.4, [treams.Material(1.2), treams.Material()]
    )
    t_treams = t_treams.expand(treams.SphericalWaveBasis.default(c.l_max))
    t_treams = t_treams.changepoltype()

    t = t[:, : t_treams.shape[-2], : t_treams.shape[-1]]
    plot_me_vs_treams(t, t_treams, c.path)
    delta_treams_vs_me(t, t_treams)
 