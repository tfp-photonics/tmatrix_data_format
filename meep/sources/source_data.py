from dataclasses import dataclass

import meep as mp
from numpy.typing import NDArray


@dataclass
class SourceData:
    meep_sources: list[mp.Source] = None
    k_src: NDArray = None
    pol_src: NDArray = None
