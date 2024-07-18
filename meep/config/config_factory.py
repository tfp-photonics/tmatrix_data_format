import meep as mp
import numpy as np
from mpi4py import MPI
from config.config import Config

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()

def get_four_sphere_config() -> Config:
    a = 0.3
    radii = [0.05, 0.06, 0.07, 0.08]
    positions = [[-0.5, -np.sqrt(3)/6, -np.sqrt(6)/12],[ -0.5, np.sqrt(3)/6, np.sqrt(6)/12], [0.0,  np.sqrt(3)/3,  -np.sqrt(6)/12],  [0.0, 0.,  np.sqrt(6)/4]]
    geometry = [
        mp.Sphere(
            radii[0], center=a*mp.Vector3(*positions[0]), material=mp.Medium(epsilon=9.)
        ),
        mp.Sphere(
            radii[1], center=a*mp.Vector3(*positions[1]), material=mp.Medium(epsilon=9.)
        ),
        mp.Sphere(
            radii[2], center=a*mp.Vector3(*positions[2]), material=mp.Medium(epsilon=9.)
        ),
        mp.Sphere(
            radii[3], center=a*mp.Vector3(*positions[3]), material=mp.Medium(epsilon=9.)
        )
    ]
    dom = 1.2
    f_min  = 1/0.6
    f_max  = 1/0.5
    nf = 10
    c = Config(
        path="./meep_input/",
        path_output="./output/",
        resolution = 90,
        f_min = f_min,
        f_max = f_max,
        f_src = (f_min+f_max)/2,
        df_src = f_max -f_min,
        n_freq = nf,
        sim_amount_mult=2,
        load_simulations=False,
        object_size=np.array([dom, dom, dom], dtype=float),
        l_max=4,
        material = [9., 9., 9., 9.],
        tol = 1e-9,
        start_geometry=geometry,
        positions = positions, 
        shape = ["sphere", "sphere", "sphere", "sphere"], 
        params = [{"radius": radii[i]} for i in range(len(radii))],
        material_names = ["", "", "", "", ""],
        material_descriptions = ["", "", "", "", ""],
        material_keywords = ["", "", "", "", ""], 
        keywords = "reciprocal, passive, lossless",
        mu = [1, 1, 1, 1, 1]
    )

    return c


def get_test_config() -> Config:
    c = Config(
        path="./meep_input/",
        path_output="./output/",
        resolution=10,
        sim_amount_mult=2,
        load_simulations=False,
        l_max=2,
        material=1.15,
        params={"radius" : 0.4},
        shape="sphere",
        start_geometry =[
        mp.Sphere(0.4, center=mp.Vector3(), material=mp.Medium(epsilon=1.15)),
    ],
        keywords = "reciprocal, passive, lossless, mirrorxyz, czinfinity",
        eps_embedding = 1.,
        cpu_cores_per_simulation=SIZE,
        opt_iterations=5,
    )
    return c


def get_ref_config() -> Config:  
    c = get_four_sphere_config()
    c.to_hdf()
    return c
