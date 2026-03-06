from __future__ import annotations

import importlib
from typing import Any, cast

import numpy as np
from numpy.typing import NDArray

from .types import CurvilinearFieldLineResult, FieldLineResult


def trace_field_lines(
    bx: NDArray[np.floating[Any]],
    by: NDArray[np.floating[Any]],
    bz: NDArray[np.floating[Any]],
    dx: NDArray[np.floating[Any]],
    dy: NDArray[np.floating[Any]],
    dz: NDArray[np.floating[Any]],
    icen_bln: NDArray[np.floating[Any]],
    jcen_bln: NDArray[np.floating[Any]],
    kcen_bln: NDArray[np.floating[Any]],
    lcen_bln: int,
    lx_bln: int,
    margin: int,
    nsubstepx: int = 3,
) -> FieldLineResult:
    """Trace magnetic field lines using the Fortran backend via f2py.

    Parameters
    ----------
    bx, by, bz:
        Magnetic field components with shape ``(ix, jx, kx)``.
    dx, dy, dz:
        Grid spacing profiles along the vertical index with shape ``(kx,)``.
    icen_bln, jcen_bln, kcen_bln:
        Per-line center coordinates with shape ``(nx_bln,)``.
    lcen_bln:
        1-based center index along the line coordinate.
    lx_bln:
        Number of points along each field line.
    margin:
        Number of ghost cells for periodic boundary handling.
    nsubstepx:
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

    ix, jx, kx = bx64.shape

    dx64 = _as_float64_array("dx", dx, ndim=1)
    dy64 = _as_float64_array("dy", dy, ndim=1)
    dz64 = _as_float64_array("dz", dz, ndim=1)

    if dx64.shape != (kx,) or dy64.shape != (kx,) or dz64.shape != (kx,):
        raise ValueError("dx, dy, and dz must have shape (kx,).")

    icen64 = _as_float64_array("icen_bln", icen_bln, ndim=1)
    jcen64 = _as_float64_array("jcen_bln", jcen_bln, ndim=1)
    kcen64 = _as_float64_array("kcen_bln", kcen_bln, ndim=1)

    if not (icen64.shape == jcen64.shape == kcen64.shape):
        raise ValueError("icen_bln, jcen_bln, and kcen_bln must have the same shape.")

    nx_bln = icen64.shape[0]

    if lx_bln <= 0:
        raise ValueError("lx_bln must be greater than 0.")
    if not (1 <= lcen_bln <= lx_bln):
        raise ValueError("lcen_bln must satisfy 1 <= lcen_bln <= lx_bln.")
    if margin < 0:
        raise ValueError("margin must be non-negative.")
    if nsubstepx <= 0:
        raise ValueError("nsubstepx must be greater than 0.")

    i_bln = np.zeros((nx_bln, lx_bln), dtype=np.float64, order="F")
    j_bln = np.zeros((nx_bln, lx_bln), dtype=np.float64, order="F")
    k_bln = np.zeros((nx_bln, lx_bln), dtype=np.float64, order="F")
    lmin_bln = np.zeros(nx_bln, dtype=np.int32)
    lmax_bln = np.zeros(nx_bln, dtype=np.int32)

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
        lcen_bln=int(lcen_bln),
        icen_bln=np.asfortranarray(icen64),
        jcen_bln=np.asfortranarray(jcen64),
        kcen_bln=np.asfortranarray(kcen64),
        bx=np.asfortranarray(bx64),
        by=np.asfortranarray(by64),
        bz=np.asfortranarray(bz64),
        dx=np.asfortranarray(dx64),
        dy=np.asfortranarray(dy64),
        dz=np.asfortranarray(dz64),
        nsubstepx=int(nsubstepx),
        margin=int(margin),
    )

    return FieldLineResult(
        i=i_bln,
        j=j_bln,
        k=k_bln,
        lmin=lmin_bln,
        lmax=lmax_bln,
        lcen=int(lcen_bln),
        nx=int(nx_bln),
        lx=int(lx_bln),
    )


def trace_field_lines_curvilinear(
    bxi: NDArray[np.floating[Any]],
    bet: NDArray[np.floating[Any]],
    bzt: NDArray[np.floating[Any]],
    dxi: float,
    det: float,
    dzt: float,
    hxi: NDArray[np.floating[Any]],
    het: NDArray[np.floating[Any]],
    hzt: NDArray[np.floating[Any]],
    icen_bln: NDArray[np.floating[Any]],
    jcen_bln: NDArray[np.floating[Any]],
    kcen_bln: NDArray[np.floating[Any]],
    lcen_bln: int,
    lx_bln: int,
    margin: int,
    nsubstepx: int = 3,
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
    icen_bln, jcen_bln, kcen_bln:
        Per-line seed coordinates as 1-based grid indices with shape
        ``(nx_bln,)``.
    lcen_bln:
        1-based center index along the line coordinate.
    lx_bln:
        Number of points along each field line.
    margin:
        Number of ghost cells for periodic boundary handling.
    nsubstepx:
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

    ix, jx, kx = bxi64.shape

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
    b_tilde_eta = bet64 * (hzt64 * hxi64)
    b_tilde_zet = bzt64 * (hxi64 * het64)

    dx_arr = np.full(kx, dxi_f)
    dy_arr = np.full(kx, det_f)
    dz_arr = np.full(kx, dzt_f)

    raw = trace_field_lines(
        bx=b_tilde_xi,
        by=b_tilde_eta,
        bz=b_tilde_zet,
        dx=dx_arr,
        dy=dy_arr,
        dz=dz_arr,
        icen_bln=icen_bln,
        jcen_bln=jcen_bln,
        kcen_bln=kcen_bln,
        lcen_bln=lcen_bln,
        lx_bln=lx_bln,
        margin=margin,
        nsubstepx=nsubstepx,
    )

    return CurvilinearFieldLineResult(
        xi=(raw.i - 1.0) * dxi_f,
        eta=(raw.j - 1.0) * det_f,
        zeta=(raw.k - 1.0) * dzt_f,
        lmin=raw.lmin,
        lmax=raw.lmax,
        lcen=raw.lcen,
        nx=raw.nx,
        lx=raw.lx,
    )


def _as_float64_array(name: str, value: Any, ndim: int) -> NDArray[np.float64]:
    arr = np.asarray(value, dtype=np.float64)
    if arr.ndim != ndim:
        raise ValueError(f"{name} must be a {ndim}D array.")
    return arr
