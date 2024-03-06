import meep as mp
import numpy as np
from mpi4py import MPI

from config.config import Config

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()


def get_six_sphere_config() -> Config:

    radii = [0.25]*6

    positions = [[-0.5, 0., 0.],[ 0.5, 0., 0.], [0.0,  -0.5,  0.],  [0.0, 0.5,  0.0], [0.0, 0.0, -0.5], [0.0, 0.0, 0.5]]
    epsilons = [1.15]*6
    geometry = [
        mp.Sphere(
            radii[i], center=mp.Vector3(*positions[i]), material=mp.Medium(epsilon=epsilons[i])
        ) for i in range(len(radii))

    ]

    c = Config(
        path="./meep_input",
        resolution=30,
        sim_amount_mult=2,
        load_simulations=False,
        object_size=np.array([1.0, 1.0, 1.0], dtype=float),
        l_max=3,
        start_geometry=geometry,
        positions = positions,
        shape = ["sphere"]*6,
        params = [{"radius": radii[i]} for i in range(len(radii))],
        material = epsilons 
    )

    return c


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

    c = Config(
        path="./meep_input",
        resolution=30,
        sim_amount_mult=2,
        load_simulations=False,
        object_size=np.array([1.0, 1.0, 1.0], dtype=float),
        l_max=3,
        start_geometry=geometry,
        positions = positions, 
        shape = ["sphere", "sphere", "sphere", "sphere"], 
        params = [{"radius": radii[i]} for i in range(len(radii))]
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
    positions = [[-0.4, -0.4, -0.4], [0.5, 0.5, 0.]]
    rad1 = 0.4
    rad2 = 0.3
    epsilons = [1.1, 1.2]
    geometry = [
        mp.Sphere(
            rad1, 
            center=mp.Vector3(*positions[0]), 
            material=mp.Medium(epsilons[0])
        ),
        mp.Sphere(
            rad2,
            center=mp.Vector3(*positions[1]),
            material=mp.Medium(epsilons[1]),
        ),
    ]

    c = Config(
        resolution=15,
        sim_amount_mult=2,
        load_simulations=False,
        l_max=3,
        object_size=np.array([1.7, 1.7, 1.7]),
        path="./meep_input",
        shape=["sphere", "sphere"],
        params=[{"radius": rad1}, {"radius" :rad2} ],
        material = epsilons,
        start_geometry=geometry,
        positions = positions
    )

    return c


def get_two_spheres_best_config() -> Config:
    positions = [[-0.4, -0.4, -0.4], [0.5, 0.5, 0.]]
    rad1 = 0.4
    rad2 = 0.3
    eps1 = 1.1
    eps2 = 1.2
    geometry = [
        mp.Sphere(
            rad1, center=mp.Vector3(*positions[0]), material=mp.Medium(epsilon=eps1)
        ),
        mp.Sphere(
            rad2,
            center=mp.Vector3(*positions[1]),
            material=mp.Medium(epsilon=eps2),
        ),
    ]

    c = Config(
        resolution=30,
        sim_amount_mult=3,
        load_simulations=False,
        l_max=5,
        object_size=np.array([1.7, 1.7, 1.7]),
        path="./meep_input/",
        shape=["sphere", "sphere"],
        params=[{"radius": rad1}, {"radius" : rad2} ],
        material = [eps1, eps2],
        start_geometry=geometry,
    )

    return c


def get_three_spheres_config() -> Config:
    positions = np.array([[-0.4, -0.4, -0.4], [0.5, 0.5, -0.1], [-0.2, 0.3, 0.3]])
    radii = [0.4, 0.3, 0.45]
    epsilons = [1.1, 1.2, 1.3]
    geometry = [
        mp.Sphere(
            radii[0], center=mp.Vector3(*positions[0]), material=mp.Medium(epsilon=epsilons[0])
        ),
        mp.Sphere(
            radii[1],
            center=mp.Vector3(*positions[1]),
            material=mp.Medium(epsilon=epsilons[1]),
        ),
        mp.Sphere(
            radii[2],
            center=mp.Vector3(*positions[2]),
            material=mp.Medium(epsilon=epsilons[2]),
        ),
    ]

    c = Config(
        resolution=20,
        sim_amount_mult=2,
        load_simulations=True,
        l_max=1,
        object_size=np.array([1.65, 1.65, 1.65]),
        path="./meep_input",
        start_geometry=geometry,
        shape=["sphere"]*3,
        params=[{"radius": radii[i]} for i in range(len(radii))],
        material = epsilons,
        positions = positions
        )

    return c


def get_sphere_embedded_config() -> Config:
    c = Config(
        resolution=20,
        sim_amount_mult=2,
        load_simulations=False,
        l_max=3,
        path="./meep_input",
        opt_eps_min = 1.1,
        eps_embedding = 1.1,
        )
    return c


def get_embedded_cylinder_config() -> Config:
    eps_embedding = 1.44**2
    eps_mat = 11.7
    radius = 0.05
    height = 0.06
    geometry = [
        mp.Cylinder(radius=radius, height=height, material=mp.Medium(epsilon=eps_mat)),
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
        shape = "cylinder", 
        params={"radius": radius, "height" : height},
        material = eps_mat
    )
    return c


def get_test_config() -> Config:
    c = Config(
        path="./meep_input/",
        resolution=10,
        sim_amount_mult=2,
        load_simulations=False,
        l_max=1,
        cpu_cores_per_simulation=SIZE,
        opt_iterations=5,
    )
    return c


def get_sphere_opt_config() -> Config:
    rad = 0.5
    eps = 1.2
    geometry = [
        mp.Sphere(
            rad,
            center=mp.Vector3(),
            material=mp.Medium(epsilon=eps),
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
        shape = "sphere",
        params = {"radius": rad},
        material = eps,

        )
    return c


def get_ref_config() -> Config:
    
    c = get_four_sphere_config()
    c.to_hdf()
    return c
