from __future__ import annotations

from dataclasses import dataclass

from numpy.typing import NDArray

from utility.data import T_DUMP_NAME, PersitableDataClass


@dataclass(frozen=True, slots=True)
class TData(PersitableDataClass):
    t: NDArray = None

    @staticmethod
    def _get_full_path(path: str, it: int) -> str:
        return f"{path}/{T_DUMP_NAME}{it}.h5"

    def to_hdf(self, path: str, it: int) -> None:
        self._to_hdf(self._get_full_path(path, it))

    @classmethod
    def from_hdf(cls, path: str, it: int) -> TData:
        return cls._from_hdf(cls._get_full_path(path, it))
