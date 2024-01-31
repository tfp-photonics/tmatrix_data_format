from __future__ import annotations

from dataclasses import dataclass
from os.path import exists
from typing import Optional

from numpy.typing import NDArray

from utility.data import SIM_DUMP_INC_NAME, SIM_DUMP_SCA_NAME, PersitableDataClass


@dataclass(frozen=True, slots=True)
class SimData(PersitableDataClass):
    # e.shape = b.shape = (frequency, cube side, grid index 1, grid index 2, field components xyz)
    # pos.shape = (cube side, grid index 1, grid index 2, xyz)
    # weights.shape = (cube side, grid index 1, grid index 2)
    # normals.shape = (cube side, xyz)
    # src_k.shape = (xyz)
    e: NDArray = None
    b: NDArray = None
    normals: NDArray = None
    pos: NDArray = None
    weights: NDArray = None

    @staticmethod
    def _get_full_path_incident(path: str, identifier: str) -> str:
        return f"{path}/{SIM_DUMP_INC_NAME}{identifier}.h5"

    @staticmethod
    def _get_full_path_scattered(path: str, id: int) -> str:
        return f"{path}/{SIM_DUMP_SCA_NAME}{id}.h5"

    def to_hdf_incident(self, path: str, identifier: str) -> None:
        self._to_hdf(self._get_full_path_incident(path, identifier))

    def to_hdf_scattered(self, path: str, id: int) -> None:
        self._to_hdf(self._get_full_path_scattered(path, id))

    @classmethod
    def from_hdf_incident(cls, path: str, identifier: str) -> Optional[SimData]:
        try:
            return cls._from_hdf(cls._get_full_path_incident(path, identifier))
        except Exception as e:
            return None

    @classmethod
    def from_hdf_scattered(cls, path: str, id: int) -> Optional[SimData]:
        try:
            return cls._from_hdf(cls._get_full_path_scattered(path, id))
        except Exception as e:
            return None
