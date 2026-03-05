from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

from .types import FieldLineResult, ObservationMapResult, OpenFieldResult


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
        from . import _bbtobln
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


def compute_open_field_fraction(
    i_bln: NDArray[np.floating[Any]],
    j_bln: NDArray[np.floating[Any]],
    k_bln: NDArray[np.floating[Any]],
    lmin_bln: NDArray[np.integer[Any]],
    lmax_bln: NDArray[np.integer[Any]],
    bzt: NDArray[np.floating[Any]],
    k_min: int,
    k_max: int,
    margin: int,
) -> OpenFieldResult:
    """Compute open-field area filling factors from traced field lines.

    Parameters
    ----------
    i_bln, j_bln, k_bln:
        Field-line coordinates with shape ``(nx_bln, lx_bln)``.
    lmin_bln, lmax_bln:
        Valid line-index bounds for each line with shape ``(nx_bln,)``.
    bzt:
        Vertical magnetic component used for normalization, shape ``(ix, jx, kx)``.
    k_min, k_max:
        Inclusive 1-based vertical index range where the filling factor is evaluated.
    margin:
        Number of ghost cells for periodic boundary handling.

    Returns
    -------
    OpenFieldResult
        Filled factor field ``f_opn`` with shape ``(ix, jx, kx)``.

    Raises
    ------
    ValueError
        If shapes or index bounds are invalid.
    ImportError
        If the compiled Fortran extension is not available.
    """

    i64 = _as_float64_array("i_bln", i_bln, ndim=2)
    j64 = _as_float64_array("j_bln", j_bln, ndim=2)
    k64 = _as_float64_array("k_bln", k_bln, ndim=2)
    if i64.shape != j64.shape or i64.shape != k64.shape:
        raise ValueError("i_bln, j_bln, and k_bln must have the same shape.")
    nx_bln, lx_bln = i64.shape

    lmin32 = _as_int32_array("lmin_bln", lmin_bln, ndim=1)
    lmax32 = _as_int32_array("lmax_bln", lmax_bln, ndim=1)
    if lmin32.shape != (nx_bln,) or lmax32.shape != (nx_bln,):
        raise ValueError("lmin_bln and lmax_bln must have shape (nx_bln,).")

    bzt64 = _as_float64_array("bzt", bzt, ndim=3)
    ix, jx, kx = bzt64.shape
    if not (1 <= k_min <= k_max <= kx):
        raise ValueError("k_min and k_max must satisfy 1 <= k_min <= k_max <= kx.")
    if margin < 0:
        raise ValueError("margin must be non-negative.")

    try:
        from . import _bbtobln
    except ImportError as exc:
        raise ImportError(
            "Failed to import the Fortran extension. Run `pip install -e '.[dev]'`."
        ) from exc

    f_opn = np.zeros((ix, jx, kx), dtype=np.float64, order="F")
    _bbtobln.blntofopn(
        f_opn=f_opn,
        i_bln=np.asfortranarray(i64),
        j_bln=np.asfortranarray(j64),
        k_bln=np.asfortranarray(k64),
        lmin_bln=np.asfortranarray(lmin32),
        lmax_bln=np.asfortranarray(lmax32),
        bzt=np.asfortranarray(bzt64),
        k_min=int(k_min),
        k_max=int(k_max),
        margin=int(margin),
    )
    return OpenFieldResult(f_opn=f_opn)


def map_field_lines_to_height(
    i_bln: NDArray[np.floating[Any]],
    j_bln: NDArray[np.floating[Any]],
    k_bln: NDArray[np.floating[Any]],
    lmin_bln: NDArray[np.integer[Any]],
    lmax_bln: NDArray[np.integer[Any]],
    lcen_bln: int,
    k_obs: int,
) -> ObservationMapResult:
    """Map field-line center points to a target observation height.

    Parameters
    ----------
    i_bln, j_bln, k_bln:
        Field-line coordinates with shape ``(nx_bln, lx_bln)``.
    lmin_bln, lmax_bln:
        Valid line-index bounds for each line with shape ``(nx_bln,)``.
    lcen_bln:
        1-based center index along the line coordinate.
    k_obs:
        1-based target vertical index for observation.

    Returns
    -------
    ObservationMapResult
        Arrays of mapped horizontal positions and vertical mismatch ``dk_obs``.

    Raises
    ------
    ValueError
        If input shapes or center index are invalid.
    ImportError
        If the compiled Fortran extension is not available.
    """

    i64 = _as_float64_array("i_bln", i_bln, ndim=2)
    j64 = _as_float64_array("j_bln", j_bln, ndim=2)
    k64 = _as_float64_array("k_bln", k_bln, ndim=2)
    if i64.shape != j64.shape or i64.shape != k64.shape:
        raise ValueError("i_bln, j_bln, and k_bln must have the same shape.")
    nx_bln, lx_bln = i64.shape

    lmin32 = _as_int32_array("lmin_bln", lmin_bln, ndim=1)
    lmax32 = _as_int32_array("lmax_bln", lmax_bln, ndim=1)
    if lmin32.shape != (nx_bln,) or lmax32.shape != (nx_bln,):
        raise ValueError("lmin_bln and lmax_bln must have shape (nx_bln,).")

    if not (1 <= lcen_bln <= lx_bln):
        raise ValueError("lcen_bln must satisfy 1 <= lcen_bln <= lx_bln.")

    try:
        from . import _bbtobln
    except ImportError as exc:
        raise ImportError(
            "Failed to import the Fortran extension. Run `pip install -e '.[dev]'`."
        ) from exc

    i_obs = np.zeros(nx_bln, dtype=np.float64, order="F")
    j_obs = np.zeros(nx_bln, dtype=np.float64, order="F")
    dk_obs = np.zeros(nx_bln, dtype=np.float64, order="F")

    _bbtobln.blntobmap(
        i_obs=i_obs,
        j_obs=j_obs,
        dk_obs=dk_obs,
        i_bln=np.asfortranarray(i64),
        j_bln=np.asfortranarray(j64),
        k_bln=np.asfortranarray(k64),
        lmin_bln=np.asfortranarray(lmin32),
        lmax_bln=np.asfortranarray(lmax32),
        lcen_bln=int(lcen_bln),
        k_obs=int(k_obs),
    )
    return ObservationMapResult(i_obs=i_obs, j_obs=j_obs, dk_obs=dk_obs)


def _as_float64_array(name: str, value: Any, ndim: int) -> NDArray[np.float64]:
    arr = np.asarray(value, dtype=np.float64)
    if arr.ndim != ndim:
        raise ValueError(f"{name} must be a {ndim}D array.")
    return arr


def _as_int32_array(name: str, value: Any, ndim: int) -> NDArray[np.int32]:
    arr = np.asarray(value, dtype=np.int32)
    if arr.ndim != ndim:
        raise ValueError(f"{name} must be a {ndim}D array.")
    return arr
