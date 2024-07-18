from __future__ import annotations

from dataclasses import dataclass

from numpy.typing import NDArray

from utility.data import COEF_DUMP_NAME, PersitableDataClass


@dataclass(frozen=True, slots=True)
class CoefData(PersitableDataClass):
    # shape = (frequencies, simulation, coefficient index)
    incident_matrix: NDArray = None
    scatter_matrix: NDArray = None

    @staticmethod
    def _get_full_path(path: str) -> str:
        return f"{path}/{COEF_DUMP_NAME}.h5"

    def to_hdf(self, path: str) -> None:
        self._to_hdf(self._get_full_path(path))

    @classmethod
    def from_hdf(cls, path: str) -> CoefData:
        return cls._from_hdf(cls._get_full_path(path))
