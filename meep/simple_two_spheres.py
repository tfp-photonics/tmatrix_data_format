import jax
import meep as mp
import numpy as np
import treams

from config.config import Config
from decomposition.t_data import TData
from optimization.optimizer import calculate_T
from testing.delta import delta_treams_vs_me
from testing.t_plots import plot_me_vs_treams

jax.config.update("jax_enable_x64", True)


def get_treams_t(l_max=1, f=1) -> treams.PhysicsArray:
    objects = [
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.2, [treams.Material(1.2), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.2, [treams.Material(1.2), treams.Material()]
        ),
    ]
    positions = [
        [0.3, 0.0, 0.0],
        [-0.3, 0.0, 0.0],
    ]
    t_treams = treams.TMatrix.cluster(objects, positions).interaction.solve()
    t_treams = t_treams.expand(treams.SphericalWaveBasis.default(l_max))
    return t_treams.changepoltype()


if __name__ == "__main__":
    path = f"/scratch/local/pscherer/tmat_meep_out/horeka_results/simple_two_spheres"
    geometry = [
        mp.Sphere(0.2, center=mp.Vector3(0.3, 0, 0), material=mp.Medium(epsilon=1.2)),
        mp.Sphere(0.2, center=mp.Vector3(-0.3, 0, 0), material=mp.Medium(epsilon=1.2)),
    ]
    c = Config(
        cpu_cores_per_simulation=8,
        resolution=30,
        sim_amount_mult=2,
        load_simulations=True,
        l_max=3,
        object_size=np.array([1.2, 1.2, 1.2]),
        path=path,
        start_geometry=geometry,
    )
    # t = calculate_T(c)
    t = TData.from_hdf(path, 0).t
    t_treams = get_treams_t(c.l_max, c.f_cen)

    t = t[:, : t_treams.shape[-2], : t_treams.shape[-1]]
    plot_me_vs_treams(t, t_treams, c.path)
    delta_treams_vs_me(t, t_treams)
