from __future__ import annotations

import importlib
from typing import Any, cast

import numpy as np
from numpy.typing import NDArray

from .types import CurvilinearFieldLineResult, FieldLineResult


def trace_field_lines(
    *,
    bx: NDArray[np.floating[Any]],
    by: NDArray[np.floating[Any]],
    bz: NDArray[np.floating[Any]],
    dx: NDArray[np.floating[Any]],
    dy: NDArray[np.floating[Any]],
    dz: NDArray[np.floating[Any]],
    seed_i: NDArray[np.floating[Any]],
    seed_j: NDArray[np.floating[Any]],
    seed_k: NDArray[np.floating[Any]],
    line_center: int,
    line_length: int,
    margin: int = 0,
    n_substeps: int = 3,
) -> FieldLineResult:
    """Trace magnetic field lines using the Fortran backend via f2py.

    Parameters
    ----------
    bx, by, bz:
        Magnetic field components with shape ``(ix, jx, kx)``.
    dx, dy, dz:
        Grid spacing profiles along the vertical index with shape ``(kx,)``.
    seed_i, seed_j, seed_k:
        Per-line center coordinates with shape ``(num_lines,)``.
    line_center:
        1-based center index along the line coordinate.
    line_length:
        Number of points along each field line.
    margin:
        Number of ghost cells for periodic boundary handling.
    n_substeps:
        Number of integration substeps per line segment.

    Returns
    -------
    FieldLineResult
        Traced line coordinates and valid index ranges.

    Raises
    ------
    ValueError
        If input shapes or scalar arguments are invalid.
    ImportError
        If the compiled Fortran extension is not available.
    """

    bx64 = _as_float64_array("bx", bx, ndim=3)
    by64 = _as_float64_array("by", by, ndim=3)
    bz64 = _as_float64_array("bz", bz, ndim=3)

    if bx64.shape != by64.shape or bx64.shape != bz64.shape:
        raise ValueError("bx, by, and bz must have the same shape.")

    _, _, kx = bx64.shape

    dx64 = _as_float64_array("dx", dx, ndim=1)
    dy64 = _as_float64_array("dy", dy, ndim=1)
    dz64 = _as_float64_array("dz", dz, ndim=1)

    if dx64.shape != (kx,) or dy64.shape != (kx,) or dz64.shape != (kx,):
        raise ValueError("dx, dy, and dz must have shape (kx,).")

    seed_i64 = _as_float64_array("seed_i", seed_i, ndim=1)
    seed_j64 = _as_float64_array("seed_j", seed_j, ndim=1)
    seed_k64 = _as_float64_array("seed_k", seed_k, ndim=1)

    if not (seed_i64.shape == seed_j64.shape == seed_k64.shape):
        raise ValueError("seed_i, seed_j, and seed_k must have the same shape.")

    num_lines = seed_i64.shape[0]

    if line_length <= 0:
        raise ValueError("line_length must be greater than 0.")
    if not (1 <= line_center <= line_length):
        raise ValueError("line_center must satisfy 1 <= line_center <= line_length.")
    if margin < 0:
        raise ValueError("margin must be non-negative.")
    if n_substeps <= 0:
        raise ValueError("n_substeps must be greater than 0.")

    i_bln = np.zeros((num_lines, line_length), dtype=np.float64, order="F")
    j_bln = np.zeros((num_lines, line_length), dtype=np.float64, order="F")
    k_bln = np.zeros((num_lines, line_length), dtype=np.float64, order="F")
    lmin_bln = np.zeros(num_lines, dtype=np.int32)
    lmax_bln = np.zeros(num_lines, dtype=np.int32)

    try:
        _bbtobln = cast(Any, importlib.import_module("wind3d_field_lines._bbtobln"))
    except ImportError as exc:
        raise ImportError(
            "Failed to import the Fortran extension. Run `pip install -e '.[dev]'`."
        ) from exc

    _bbtobln.bbtobln(
        i_bln=i_bln,
        j_bln=j_bln,
        k_bln=k_bln,
        lmin_bln=lmin_bln,
        lmax_bln=lmax_bln,
        lcen_bln=int(line_center),
        icen_bln=np.asfortranarray(seed_i64),
        jcen_bln=np.asfortranarray(seed_j64),
        kcen_bln=np.asfortranarray(seed_k64),
        bx=np.asfortranarray(bx64),
        by=np.asfortranarray(by64),
        bz=np.asfortranarray(bz64),
        dx=np.asfortranarray(dx64),
        dy=np.asfortranarray(dy64),
        dz=np.asfortranarray(dz64),
        nsubstepx=int(n_substeps),
        margin=int(margin),
    )

    return FieldLineResult(
        i=i_bln,
        j=j_bln,
        k=k_bln,
        lmin=lmin_bln,
        lmax=lmax_bln,
        line_center=int(line_center),
        num_lines=int(num_lines),
        line_length=int(line_length),
    )


def trace_field_lines_curvilinear(
    *,
    bxi: NDArray[np.floating[Any]],
    bet: NDArray[np.floating[Any]],
    bzt: NDArray[np.floating[Any]],
    dxi: float,
    det: float,
    dzt: float,
    hxi: NDArray[np.floating[Any]],
    het: NDArray[np.floating[Any]],
    hzt: NDArray[np.floating[Any]],
    seed_i: NDArray[np.floating[Any]],
    seed_j: NDArray[np.floating[Any]],
    seed_k: NDArray[np.floating[Any]],
    line_center: int,
    line_length: int,
    margin: int = 0,
    n_substeps: int = 3,
) -> CurvilinearFieldLineResult:
    """Trace magnetic field lines in an orthogonal curvilinear coordinate system.

    The field-line tracing problem in orthogonal curvilinear coordinates reduces
    to the standard Cartesian streamline problem after scaling the magnetic field
    components by products of the scale factors.  See the theory documentation
    for the derivation.

    Parameters
    ----------
    bxi, bet, bzt:
        Physical components of the magnetic field (B_xi, B_eta, B_zeta) with
        shape ``(ix, jx, kx)``.
    dxi, det, dzt:
        Uniform grid spacing in the xi, eta, and zeta directions (scalars).
    hxi, het, hzt:
        Scale factors h_xi, h_eta, h_zeta with shape ``(kx,)``.  Each scale
        factor may depend on the zeta index only.
    seed_i, seed_j, seed_k:
        Per-line seed coordinates as 1-based grid indices with shape
        ``(num_lines,)``.
    line_center:
        1-based center index along the line coordinate.
    line_length:
        Number of points along each field line.
    margin:
        Number of ghost cells for periodic boundary handling.
    n_substeps:
        Number of integration substeps per line segment.

    Returns
    -------
    CurvilinearFieldLineResult
        Traced field-line coordinates in physical (xi, eta, zeta) space and
        valid index ranges.

    Raises
    ------
    ValueError
        If input shapes or scalar arguments are invalid.
    ImportError
        If the compiled Fortran extension is not available.
    """

    bxi64 = _as_float64_array("bxi", bxi, ndim=3)
    bet64 = _as_float64_array("bet", bet, ndim=3)
    bzt64 = _as_float64_array("bzt", bzt, ndim=3)

    if bxi64.shape != bet64.shape or bxi64.shape != bzt64.shape:
        raise ValueError("bxi, bet, and bzt must have the same shape.")

    _, _, kx = bxi64.shape

    hxi64 = _as_float64_array("hxi", hxi, ndim=1)
    het64 = _as_float64_array("het", het, ndim=1)
    hzt64 = _as_float64_array("hzt", hzt, ndim=1)

    if hxi64.shape != (kx,) or het64.shape != (kx,) or hzt64.shape != (kx,):
        raise ValueError("hxi, het, and hzt must have shape (kx,).")

    if np.any(hxi64 <= 0) or np.any(het64 <= 0) or np.any(hzt64 <= 0):
        raise ValueError("hxi, het, and hzt must be positive.")

    dxi_f = float(dxi)
    det_f = float(det)
    dzt_f = float(dzt)

    if dxi_f <= 0 or det_f <= 0 or dzt_f <= 0:
        raise ValueError("dxi, det, and dzt must be positive.")

    # Scale the magnetic field components.
    # From the theory: B_tilde_xi = h_eta * h_zeta * B_xi, etc.
    # hxi64 has shape (kx,) and broadcasts against bxi64 of shape (ix, jx, kx).
    b_tilde_xi = bxi64 * (het64 * hzt64)
    b_tilde_et = bet64 * (hzt64 * hxi64)
    b_tilde_zt = bzt64 * (hxi64 * het64)

    dx_arr = np.full(kx, dxi_f)
    dy_arr = np.full(kx, det_f)
    dz_arr = np.full(kx, dzt_f)

    raw = trace_field_lines(
        bx=b_tilde_xi,
        by=b_tilde_et,
        bz=b_tilde_zt,
        dx=dx_arr,
        dy=dy_arr,
        dz=dz_arr,
        seed_i=seed_i,
        seed_j=seed_j,
        seed_k=seed_k,
        line_center=line_center,
        line_length=line_length,
        margin=margin,
        n_substeps=n_substeps,
    )

    return CurvilinearFieldLineResult(
        xi=(raw.i - 1.0) * dxi_f,
        eta=(raw.j - 1.0) * det_f,
        zeta=(raw.k - 1.0) * dzt_f,
        lmin=raw.lmin,
        lmax=raw.lmax,
        line_center=raw.line_center,
        num_lines=raw.num_lines,
        line_length=raw.line_length,
    )


def _as_float64_array(name: str, value: Any, ndim: int) -> NDArray[np.float64]:
    arr = np.asarray(value, dtype=np.float64)
    if arr.ndim != ndim:
        raise ValueError(f"{name} must be a {ndim}D array.")
    return arr
