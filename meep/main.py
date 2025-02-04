import numpy as np
import treams
import meep as mp
from config.config import Config
from config.config_factory import get_test_config, get_four_sphere_config
from testing.delta import delta_treams_vs_me
from testing.t_plots import plot_me_vs_treams
from t_computation.scat import calculate_T
from store_format import h5save


def cyl():
    """
    Generate reference cylinder from the T-matrix data format project.
    """
    eps_embedding = 1.0
    eps_mat = 6.25
    radius = 0.25
    height = 0.3
    geometry = [
        mp.Cylinder(radius=radius, height=height, material=mp.Medium(epsilon=eps_mat)),
    ]
    wl_min = 0.75
    wl_max = 1.25
    f_min = 1.0 / wl_max
    f_max = 1.0 / wl_min
    size = 2 * wl_max
    nf = 10
    res = 75
    dpml = 0.5 * wl_max
    c = Config(
        path="./meep_input/",
        path_output="./output/",
        object_size=np.array([size, size, size]),
        resolution=res,
        sim_amount_mult=2,
        load_simulations=False,
        cpu_cores_per_simulation=8,
        l_max=3,
        material=eps_mat,
        dpml=dpml,
        n_freq=nf,
        f_min=f_min,
        f_max=f_max,
        f_src=(f_min + f_max) / 2,
        df_src=f_max - f_min,
        eps_embedding=eps_embedding,
        start_geometry=geometry,
        shape="cylinder",
        params={"radius": radius, "height": height},
        material_names=["Air, Vacuum", "TiO2, Titanium Dioxide"],
        keywords="reciprocal, passive, lossless, mirrorxyz, czinfinity",
    )
    t = calculate_T(c)
    with open(__file__, "r", encoding="utf-8") as scriptfile:
        script = scriptfile.read()
    h5save(
        c,
        t,
        script,
        f"cylinder_res_{res}_wl_{wl_min}_{wl_max}_{nf}_domain_{size}",
        "Reference from data format repo",
    )


def tetrahedron():
    """
    Generate the reference file from T-matrix format project:
    radii: 50, 60, 70, 80 nm at distances.
    """
    c = get_four_sphere_config()
    t = calculate_T(c)
    with open(__file__, "r", encoding="utf-8") as scriptfile:
        script = scriptfile.read()
    h5save(
        c,
        t,
        script,
        f"four_spheres_res_{c.resolution}_freq_{c.f_min}_{c.f_max}_{c.n_freq}_dom_{c.object_size[0]}",
        "Reference from data format repo",
    )


def sphere() -> None:
    """
    Generate a simple sphere with permittivity 1.15 and radius 400 nm at wavelength
    """
    c = get_test_config()
    t = calculate_T(c)
    with open(__file__, "r", encoding="utf-8") as scriptfile:
        script = scriptfile.read()
    h5save(c, t, script, "test_sphere", "Test single sphere file")
    t_treams = treams.TMatrix.sphere(
        c.l_max,
        c.f_cen * 2 * np.pi,
        c.params["radius"],
        [treams.Material(c.material), treams.Material(c.eps_embedding)],
    )
    t_treams = t_treams.expand(treams.SphericalWaveBasis.default(c.l_max))
    t_treams = t_treams.changepoltype()

    t = t[len(t) // 2, : t_treams.shape[-2], : t_treams.shape[-1]]
    print(t_treams[0, 0], t[0, 0])
    plot_me_vs_treams(t, t_treams, c.path_output)
    delta_treams_vs_me(t, t_treams)


if __name__ == "__main__":
    cyl()
