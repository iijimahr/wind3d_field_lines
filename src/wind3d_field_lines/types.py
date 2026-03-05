from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class FieldLineResult:
    """磁力線トレース結果。"""

    i: NDArray[np.float64]
    j: NDArray[np.float64]
    k: NDArray[np.float64]
    lmin: NDArray[np.int32]
    lmax: NDArray[np.int32]
    lcen: int
    nx: int
    lx: int


@dataclass(frozen=True)
class OpenFieldResult:
    """開放磁場の面積充填率。"""

    f_opn: NDArray[np.float64]


@dataclass(frozen=True)
class ObservationMapResult:
    """観測高さへのマッピング結果。"""

    i_obs: NDArray[np.float64]
    j_obs: NDArray[np.float64]
    dk_obs: NDArray[np.float64]
