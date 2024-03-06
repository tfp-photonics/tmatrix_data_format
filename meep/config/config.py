from __future__ import annotations

import hashlib
import warnings
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import meep as mp
import numpy as np
from numpy.typing import ArrayLike, NDArray

from geometry.geometry_factory import get_epsgrid_from_geometric_objects
from utility.data import CONFIG_NAME, PersitableDataClass


@dataclass(frozen=True, slots=True)
class Config(PersitableDataClass):
    path: str = "./"
    load_simulations: bool = False
    cpu_cores_per_simulation: int = 4
    object_size: ArrayLike = (1.0, 1.0, 1.0)
    dpml: float = 1.0
    eps_embedding: float = 1.0
    l_max: int = 3
    sim_amount_mult: float = 2
    resolution: int = 15
    n_freq: int = 5
    f_min: float = 0.9
    f_max: float = 1.1
    f_src: float = 1.0
    df_src: float = 0.2
    shape : list[str] | str = 'sphere'
    
    params : list[dict] | dict = field(default_factory=lambda: {"radius":0.5}) 
    material: list[complex] | complex = 1.15
    positions: list[float] = field(default_factory=lambda: [0, 0, 0])
    start_geometry: list[mp.GeometricObject] | NDArray = field(
        default_factory=lambda: [
            mp.Sphere(
                0.5,
                center=mp.Vector3(),
                material=mp.Medium(epsilon=1.15),
            )
        ]
    )
    continue_opt_iterations: bool = False
    opt_iterations: int = 1
    opt_eps_min: float = None
    opt_eps_max: float = None

    def __post_init__(self):
        print('postint', self.start_geometry)
        Path(self.path).mkdir(parents=True, exist_ok=True)
        object.__setattr__(self, "object_size", np.asarray(self.object_size))
        if isinstance(self.start_geometry, list):
            object.__setattr__(
                self,
                "start_geometry",
                get_epsgrid_from_geometric_objects(
                    res=self.resolution * 2.0,
                    size=self.object_size,
                    meep_objects=self.start_geometry,
                    eps_background=self.opt_eps_min
                    if self.opt_eps_min is not None
                    else 1.0,
                ),
            )
        if self.opt_eps_min is None:
            object.__setattr__(self, "opt_eps_min", np.min(self.start_geometry))
        if self.opt_eps_max is None:
            object.__setattr__(self, "opt_eps_max", np.max(self.start_geometry))
        if self.dpml < self.max_wl / 2:
            warnings.warn("PML size is smaller than half the maximum wavelength.")

        if (min_resolution := int(np.ceil(10 / self.min_wl))) > self.resolution:
            warnings.warn(
                "The resolution of your simulation is too small for the specified "
                "frequencies and maximum permittivity. Minimum recommended resolution "
                f"is {int(np.ceil(min_resolution))}."
            )
        if self.opt_iterations <= 0:
            raise ValueError("Number of iterations must be positive.")
        if self.dpml <= 0:
            raise ValueError("PML thickness must be positive")
        if any(s <= 0 for s in self.object_size):
            raise ValueError("Object sizes must be positive.")
        if self.opt_eps_max < self.opt_eps_min:
            raise ValueError(
                "Maximal epsilon must be larger than or equal to minimal epsilon."
            )
        if self.l_max <= 0:
            raise ValueError("L_max must be positive.")
        if self.sim_amount_mult <= 0:
            raise ValueError("Simulation amount multiplier must be positive.")
        if self.resolution <= 0:
            raise ValueError("Resolution must be positive.")
        if self.n_freq <= 0:
            raise ValueError("Number of frequencies must be positive.")
        if self.f_min >= self.f_max:
            raise ValueError(
                "Minimum frequency must be smaller than maximum frequency."
            )
        if self.f_src <= 0:
            raise ValueError("Frequency of the field sources must be positive.")
        if self.df_src <= 0:
            raise ValueError("Frequency width of the field sources must be positive.")
        if (
            np.min(self.start_geometry) < self.opt_eps_min
            or np.max(self.start_geometry) > self.opt_eps_max
        ):
            raise ValueError(
                "Provided geometry does not conform to the given maximum and minimum permittivity."
            )

    @cached_property
    def scatter_sim_params_hash(self) -> str:
        relevant_information = str.encode(
            f"{self.resolution}{self.eps_embedding}{self.object_size}{self.dpml}{self.f_src}{self.df_src}"
        )
        hash_object = hashlib.sha1(relevant_information)
        return hash_object.hexdigest()

    @cached_property
    def t_shape(self) -> tuple[int, int, int]:
        return (self.n_freq, self.min_sim_amount, self.min_sim_amount)

    @cached_property
    def no_pml_size(self) -> NDArray:
        return np.full(
            (3,),
            np.sqrt(
                self.object_size[0] ** 2
                + self.object_size[1] ** 2
                + self.object_size[2] ** 2
            ),
        )

    @cached_property
    def with_pml_size(self) -> NDArray:
        return self.no_pml_size + 2 * self.dpml

    @cached_property
    def min_sim_amount(self) -> int:
        return ((self.l_max + 1) ** 2 - 1) * 2

    @cached_property
    def sim_amount(self) -> int:
        return int(np.ceil(self.sim_amount_mult * self.min_sim_amount))

    @cached_property
    def f_cen(self) -> float:
        return (self.f_min + self.f_max) / 2

    @cached_property
    def frequencies(self) -> NDArray:
        return np.linspace(self.f_min, self.f_max, self.n_freq)

    @cached_property
    def omegas(self) -> NDArray:
        return 2 * np.pi * self.frequencies

    @cached_property
    def wave_numbers(self) -> NDArray:
        eps, mu = self.eps_embedding, 1
        return self.omegas * np.sqrt(eps * mu)

    @cached_property
    def wavelengths(self) -> NDArray:
        return 1 / self.frequencies

    @cached_property
    def min_wl(self) -> float:
        return np.min(self.wavelengths) / np.sqrt(self.opt_eps_max)

    @cached_property
    def max_wl(self) -> float:
        return np.max(self.wavelengths) / np.sqrt(self.opt_eps_min)

    @staticmethod
    def _get_full_path(root_path: str) -> str:
        return f"{root_path}/{CONFIG_NAME}.h5"

    def to_hdf(self) -> None:
        self._to_hdf(self._get_full_path(self.path))

    @classmethod
    def from_hdf(cls, path: str) -> Config:
        return cls._from_hdf(cls._get_full_path(path))
