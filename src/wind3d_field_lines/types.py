from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class FieldLineResult:
    """Container for traced field-line coordinates and index bounds."""

    i: NDArray[np.float64]
    j: NDArray[np.float64]
    k: NDArray[np.float64]
    lmin: NDArray[np.int32]
    lmax: NDArray[np.int32]
    lcen: int
    nx: int
    lx: int
