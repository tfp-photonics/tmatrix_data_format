from __future__ import annotations

from abc import ABC
from dataclasses import asdict, is_dataclass

import h5py
import meep as mp
import numpy as np

# once
CONFIG_NAME = "config"
N_EPSILON_CALCULATED_NAME = "progress"
# for each iteration
EPS_DUMP_NAME = "eps_dump_"
COEF_DUMP_NAME = "coef_dump_"
T_DUMP_NAME = "t_dump_"
LOSS_LOG_NAME = "loss_"
# for each simulation
ROTATION_DUMP_NAME = "rotation_data_"
SIM_DUMP_SCA_NAME = "sim_dump_sca_"
SIM_DUMP_INC_NAME = "sim_dump_inc_"

import jax.numpy
from jax.errors import TracerArrayConversionError


def save_dataset(file: h5py.File, key: str, data) -> None:
    try:
        del file[key]
    except Exception as e:
        pass
    file.create_dataset(
        key,
        data=data,
        compression="gzip",
        compression_opts=1,
    )


class PersitableDataClass(ABC):
    def _to_hdf(self, path, mode="a") -> None:
        if not mp.am_master():
            return
        if not is_dataclass(self):
            raise TypeError(
                "The object to be persisted is not an instance of a dataclass."
            )
        try:
            with h5py.File(path, mode) as hf:
                for k, v in asdict(self).items():
                    if isinstance(v, np.ndarray):
                        save_dataset(hf, k, v)
                    elif isinstance(v, jax.numpy.ndarray):
                        try:
                            save_dataset(hf, k, np.array(v))
                        except TracerArrayConversionError:
                            print(
                                f"Saving of dataset with key {k} failed because of TracerArrayConversionError."
                            )
                    else:
                        hf.attrs[k] = v
        except Exception as e:
            pass

    @classmethod
    def _from_hdf(cls, path) -> PersitableDataClass:
        with h5py.File(path, "r") as hf:
            dset_dict = {k: np.array(v) for k, v in hf.items()}
            return cls(**dset_dict, **hf.attrs)
