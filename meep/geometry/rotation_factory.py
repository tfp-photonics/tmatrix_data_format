import numpy as np

from config.config import Config
from geometry.rotation_data import RotationData


def get_persistent_sphere_angles(c: Config, id: int) -> RotationData:
    rotation_data = None
    if c.load_simulations:
        rotation_data = RotationData.from_hdf(path=c.path, id=id)
    if rotation_data is None:
        n_fib_points = int(np.ceil(c.sim_amount / 2.0))
        np.random.seed(0)
        permutation = np.random.permutation(n_fib_points)
        fib_id = permutation[int(np.floor(id / 2.0))]
        theta, phi = _get_fib_sphere_angles(n=fib_id, n_tot=n_fib_points)
        alpha = 0.0 if id % 2 == 0 else np.pi / 2.0
        rotation_data = RotationData(theta=theta, phi=phi, alpha=alpha)
        rotation_data.to_hdf(path=c.path, id=id)
    return rotation_data


def _get_fib_sphere_angles(n: int, n_tot: int) -> tuple[float, float]:
    theta = np.arccos(1.0 - (2.0 * n) / n_tot)
    phi = np.mod(np.pi * (1.0 + 5.0**0.5) * n, 2 * np.pi)
    return theta, phi
