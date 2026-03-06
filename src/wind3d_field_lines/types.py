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
    line_center: int
    num_lines: int
    line_length: int


@dataclass(frozen=True)
class CurvilinearFieldLineResult:
    """Container for traced field-line coordinates in curvilinear coordinate space.

    Coordinates are returned as physical values in the (xi, eta, zeta) space.
    Field-line tracing is performed using scaled magnetic field components as
    described in the orthogonal curvilinear coordinates theory.
    """

    xi: NDArray[np.float64]
    eta: NDArray[np.float64]
    zeta: NDArray[np.float64]
    lmin: NDArray[np.int32]
    lmax: NDArray[np.int32]
    line_center: int
    num_lines: int
    line_length: int
