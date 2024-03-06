from __future__ import annotations

from dataclasses import dataclass

from utility.data import LOSS_LOG_NAME, PersitableDataClass


@dataclass(slots=True)
class LossData(PersitableDataClass):
    loss: float = 0

    @staticmethod
    def _get_full_path(path: str, it: int) -> str:
        return f"{path}/{LOSS_LOG_NAME}{it}.h5"

    def to_hdf(self, path: str, it: int) -> None:
        self._to_hdf(self._get_full_path(path, it))

    @classmethod
    def from_hdf(cls, path: str, it: int) -> LossData:
        return cls._from_hdf(cls._get_full_path(path, it))
