import meep as mp
import numpy as np
from mpi4py import MPI

from config.config import Config

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()


def get_six_sphere_config() -> Config:
    geometry = [
        mp.Sphere(
            0.25, center=mp.Vector3(-0.5, 0.0, 0.0), material=mp.Medium(epsilon=1.15)
        ),
        mp.Sphere(
            0.25, center=mp.Vector3(0.5, 0.0, 0.0), material=mp.Medium(epsilon=1.15)
        ),
        mp.Sphere(
            0.25, center=mp.Vector3(0.0, -0.5, 0.0), material=mp.Medium(epsilon=1.15)
        ),
        mp.Sphere(
            0.25, center=mp.Vector3(0.0, 0.5, 0.0), material=mp.Medium(epsilon=1.15)
        ),
        mp.Sphere(
            0.25, center=mp.Vector3(0.0, 0.0, -0.5), material=mp.Medium(epsilon=1.15)
        ),
        mp.Sphere(
            0.25, center=mp.Vector3(0.0, 0.0, 0.5), material=mp.Medium(epsilon=1.15)
        ),
    ]

    c = Config(
        path="./meep_input",
        resolution=30,
        sim_amount_mult=2,
        load_simulations=False,
        object_size=np.array([1.0, 1.0, 1.0], dtype=float),
        l_max=3,
        start_geometry=geometry,
    )

    return c


def get_four_sphere_config() -> Config:
    a = 0.3
    geometry = [
        mp.Sphere(
            0.05, center=a*mp.Vector3(-0.5, -np.sqrt(3)/6, -np.sqrt(6)/12), material=mp.Medium(epsilon=1.)
        ),
        mp.Sphere(
            0.06, center=a*mp.Vector3(-0.5, np.sqrt(3)/6, np.sqrt(6)/12), material=mp.Medium(epsilon=1.)
        ),
        mp.Sphere(
            0.07, center=a*mp.Vector3(0.0,  np.sqrt(3)/3,  -np.sqrt(6)/12), material=mp.Medium(epsilon=1.)
        ),
        mp.Sphere(
            0.08, center=a*mp.Vector3(0.0, 0.,  np.sqrt(6)/4), material=mp.Medium(epsilon=1.)
        )
    ]

    c = Config(
        path="./meep_input",
        resolution=30,
        sim_amount_mult=2,
        load_simulations=False,
        object_size=np.array([1.0, 1.0, 1.0], dtype=float),
        l_max=3,
        start_geometry=geometry,
    )

    return c


def get_core_shell_config(
    path="./meep_input",
    center=mp.Vector3(),
) -> Config:

    geometry = [
        mp.Sphere(0.5, center=center, material=mp.Medium(epsilon=1.2)),
        mp.Sphere(0.25, center=center, material=mp.Medium(epsilon=1.1)),
    ]

    c = Config(
        resolution=30,
        sim_amount_mult=3,
        load_simulations=False,
        object_size=np.array([*center]) + np.array([1.0, 1.0, 1.0], dtype=float),
        path=path,
        start_geometry=geometry,
    )
    return c


def get_two_spheres_config() -> Config:
    geometry = [
        mp.Sphere(
            0.4, center=mp.Vector3(-0.4, -0.4, -0.4), material=mp.Medium(epsilon=1.1)
        ),
        mp.Sphere(
            0.3,
            center=mp.Vector3(0.5, 0.5, 0),
            material=mp.Medium(epsilon=1.2),
        ),
    ]

    c = Config(
        resolution=15,
        sim_amount_mult=2,
        load_simulations=False,
        l_max=3,
        object_size=np.array([1.7, 1.7, 1.7]),
        path="./meep_input",
        start_geometry=geometry,
    )

    return c


def get_two_spheres_best_config() -> Config:
    geometry = [
        mp.Sphere(
            0.4, center=mp.Vector3(-0.4, -0.4, -0.4), material=mp.Medium(epsilon=1.1)
        ),
        mp.Sphere(
            0.3,
            center=mp.Vector3(0.5, 0.5, 0),
            material=mp.Medium(epsilon=1.2),
        ),
    ]

    c = Config(
        resolution=30,
        sim_amount_mult=3,
        load_simulations=False,
        l_max=5,
        object_size=np.array([1.7, 1.7, 1.7]),
        path="./meep_input",
        start_geometry=geometry,
    )

    return c


def get_three_spheres_config() -> Config:
    geometry = [
        mp.Sphere(
            0.4, center=mp.Vector3(-0.4, -0.4, -0.4), material=mp.Medium(epsilon=1.1)
        ),
        mp.Sphere(
            0.3,
            center=mp.Vector3(0.5, 0.5, -0.1),
            material=mp.Medium(epsilon=1.2),
        ),
        mp.Sphere(
            0.45,
            center=mp.Vector3(-0.2, 0.3, 0.3),
            material=mp.Medium(epsilon=1.3),
        ),
    ]

    c = Config(
        resolution=20,
        sim_amount_mult=2,
        load_simulations=True,
        l_max=3,
        object_size=np.array([1.65, 1.65, 1.65]),
        path="./meep_input",
        start_geometry=geometry,
    )

    return c


def get_sphere_embedded_config() -> Config:
    c = Config(
        resolution=20,
        sim_amount_mult=2,
        load_simulations=False,
        l_max=3,
        path="./meep_input",
        opt_eps_min=1.1,
        eps_embedding=1.1,
    )
    return c


def get_embedded_cylinder_config() -> Config:
    eps_embedding = 1.44**2
    eps_mat = 11.7
    geometry = [
        mp.Cylinder(radius=0.05, height=0.06, material=mp.Medium(epsilon=eps_mat)),
    ]

    c = Config(
        path="./meep_input",
        object_size=np.array([0.1, 0.1, 0.1]),
        resolution=300,
        sim_amount_mult=3,
        load_simulations=False,
        cpu_cores_per_simulation=8,
        l_max=5,
        dpml=0.5,
        opt_eps_min=eps_embedding,
        eps_embedding=eps_embedding,
        opt_eps_max=eps_mat,
        start_geometry=geometry,
    )
    return c


def get_test_config() -> Config:
    c = Config(
        path="./meep_input",
        resolution=10,
        sim_amount_mult=2,
        load_simulations=False,
        l_max=1,
        cpu_cores_per_simulation=SIZE,
        opt_iterations=5,
    )
    return c


def get_sphere_opt_config() -> Config:
    geometry = [
        mp.Sphere(
            0.5,
            center=mp.Vector3(),
            material=mp.Medium(epsilon=1.2),
        )
    ]
    c = Config(
        path="./meep_input",
        start_geometry=geometry,
        object_size=np.array([1.2, 1.2, 1.2]),
        opt_eps_min=1,
        opt_eps_max=1.3,
        resolution=20,
        sim_amount_mult=2,
        load_simulations=False,
        continue_opt_iterations=True,
        l_max=3,
        cpu_cores_per_simulation=SIZE,
        opt_iterations=100,
    )
    return c


def get_config() -> Config:
#    c = get_test_config()
    c = get_four_sphere_config()
    c.to_hdf()
    return c
