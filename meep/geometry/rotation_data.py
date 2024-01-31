from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from os.path import exists
from typing import List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from utility.data import ROTATION_DUMP_NAME, PersitableDataClass
from utility.vector_geometry import rotation_matrix


@dataclass(frozen=True, slots=True)
class RotationData(PersitableDataClass):
    theta: float = 0
    phi: float = 0
    alpha: float = 0

    @staticmethod
    def _get_full_path(path: str, id: int) -> str:
        return f"{path}/{ROTATION_DUMP_NAME}{id}.h5"

    def to_hdf(self, path: str, id: int) -> None:
        self._to_hdf(self._get_full_path(path, id))

    @classmethod
    def from_hdf(cls, path: str, id: int) -> Optional[RotationData]:
        try:
            return cls._from_hdf(cls._get_full_path(path, id))
        except Exception as e:
            return None

    @cached_property
    def rotation_matrix(self) -> NDArray:
        up = np.array([0, 0, 1])
        left = np.array([0, 1, 0])
        rotation_left = rotation_matrix(left, self.theta)
        rotation_up_phi = rotation_matrix(up, self.phi)
        rotation_up_alpha = rotation_matrix(up, self.alpha)
        return rotation_up_phi @ rotation_left @ rotation_up_alpha

    @cached_property
    def rotation_matrix_inverse(self) -> NDArray:
        down = np.array([0, 0, -1])
        right = np.array([0, -1, 0])
        rotation_left_inverse = rotation_matrix(right, self.theta)
        rotation_up_phi_inverse = rotation_matrix(down, self.phi)
        rotation_up_alpha_inverse = rotation_matrix(down, self.alpha)
        return (
            rotation_up_alpha_inverse @ rotation_left_inverse @ rotation_up_phi_inverse
        )

    @cached_property
    def rotation_as_list_of_angles_and_rotation_planes(
        self,
    ) -> List[Tuple[float, Tuple[int]]]:
        angles = [self.alpha, -self.theta, self.phi]
        planes_of_rotation = [(0, 1), (0, 2), (0, 1)]
        return list(zip(angles, planes_of_rotation))

    @cached_property
    def rotation_as_list_of_angles_and_rotation_planes_inverse(
        self,
    ) -> List[Tuple[float, Tuple[int]]]:
        angles = [-self.phi, self.theta, -self.alpha]
        planes_of_rotation = [(0, 1), (0, 2), (0, 1)]
        return list(zip(angles, planes_of_rotation))
