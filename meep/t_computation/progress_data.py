from __future__ import annotations

from dataclasses import dataclass

from utility.data import N_EPSILON_CALCULATED_NAME, PersitableDataClass


@dataclass(frozen=True, slots=True)
class ProgressData(PersitableDataClass):
    progress: int = None

    @staticmethod
    def _get_full_path(path: str) -> str:
        return f"{path}/{N_EPSILON_CALCULATED_NAME}.h5"

    def to_hdf(self, path: str) -> None:
        self._to_hdf(self._get_full_path(path))

    @classmethod
    def from_hdf(cls, path: str) -> ProgressData:
        return cls._from_hdf(cls._get_full_path(path))
