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
    """f2py実装を用いて磁力線をトレースする。"""

    bx64 = _as_float64_array("bx", bx, ndim=3)
    by64 = _as_float64_array("by", by, ndim=3)
    bz64 = _as_float64_array("bz", bz, ndim=3)

    if bx64.shape != by64.shape or bx64.shape != bz64.shape:
        raise ValueError("bx, by, bz のshapeは一致する必要があります。")

    ix, jx, kx = bx64.shape

    dx64 = _as_float64_array("dx", dx, ndim=1)
    dy64 = _as_float64_array("dy", dy, ndim=1)
    dz64 = _as_float64_array("dz", dz, ndim=1)

    if dx64.shape != (kx,) or dy64.shape != (kx,) or dz64.shape != (kx,):
        raise ValueError("dx, dy, dz のshapeは (kx,) である必要があります。")

    icen64 = _as_float64_array("icen_bln", icen_bln, ndim=1)
    jcen64 = _as_float64_array("jcen_bln", jcen_bln, ndim=1)
    kcen64 = _as_float64_array("kcen_bln", kcen_bln, ndim=1)

    if not (icen64.shape == jcen64.shape == kcen64.shape):
        raise ValueError("icen_bln, jcen_bln, kcen_bln のshapeは一致する必要があります。")

    nx_bln = icen64.shape[0]

    if lx_bln <= 0:
        raise ValueError("lx_bln は1以上である必要があります。")
    if not (1 <= lcen_bln <= lx_bln):
        raise ValueError("lcen_bln は 1 以上 lx_bln 以下である必要があります。")
    if margin < 0:
        raise ValueError("margin は0以上である必要があります。")
    if nsubstepx <= 0:
        raise ValueError("nsubstepx は1以上である必要があります。")

    i_bln = np.zeros((nx_bln, lx_bln), dtype=np.float64, order="F")
    j_bln = np.zeros((nx_bln, lx_bln), dtype=np.float64, order="F")
    k_bln = np.zeros((nx_bln, lx_bln), dtype=np.float64, order="F")
    lmin_bln = np.zeros(nx_bln, dtype=np.int32)
    lmax_bln = np.zeros(nx_bln, dtype=np.int32)

    try:
        from . import _bbtobln
    except ImportError as exc:
        raise ImportError(
            "Fortran拡張の読み込みに失敗しました。`pip install -e .[dev]` を実行してください。"
        ) from exc

    _bbtobln.bbtobln(
        i_bln,
        j_bln,
        k_bln,
        lmin_bln,
        lmax_bln,
        int(lcen_bln),
        np.asfortranarray(icen64),
        np.asfortranarray(jcen64),
        np.asfortranarray(kcen64),
        np.asfortranarray(bx64),
        np.asfortranarray(by64),
        np.asfortranarray(bz64),
        np.asfortranarray(dx64),
        np.asfortranarray(dy64),
        np.asfortranarray(dz64),
        int(nsubstepx),
        int(margin),
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
    """磁力線配列から開放磁場の面積充填率を計算する。"""

    i64 = _as_float64_array("i_bln", i_bln, ndim=2)
    j64 = _as_float64_array("j_bln", j_bln, ndim=2)
    k64 = _as_float64_array("k_bln", k_bln, ndim=2)
    if i64.shape != j64.shape or i64.shape != k64.shape:
        raise ValueError("i_bln, j_bln, k_bln のshapeは一致する必要があります。")
    nx_bln, lx_bln = i64.shape

    lmin32 = _as_int32_array("lmin_bln", lmin_bln, ndim=1)
    lmax32 = _as_int32_array("lmax_bln", lmax_bln, ndim=1)
    if lmin32.shape != (nx_bln,) or lmax32.shape != (nx_bln,):
        raise ValueError("lmin_bln, lmax_bln のshapeは (nx_bln,) である必要があります。")

    bzt64 = _as_float64_array("bzt", bzt, ndim=3)
    ix, jx, kx = bzt64.shape
    if not (1 <= k_min <= k_max <= kx):
        raise ValueError("k_min, k_max は 1 <= k_min <= k_max <= kx を満たす必要があります。")
    if margin < 0:
        raise ValueError("margin は0以上である必要があります。")

    try:
        from . import _bbtobln
    except ImportError as exc:
        raise ImportError(
            "Fortran拡張の読み込みに失敗しました。`pip install -e .[dev]` を実行してください。"
        ) from exc

    f_opn = np.zeros((ix, jx, kx), dtype=np.float64, order="F")
    _bbtobln.blntofopn(
        f_opn,
        np.asfortranarray(i64),
        np.asfortranarray(j64),
        np.asfortranarray(k64),
        np.asfortranarray(lmin32),
        np.asfortranarray(lmax32),
        np.asfortranarray(bzt64),
        int(k_min),
        int(k_max),
        int(margin),
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
    """磁力線の中心高さから観測高さへのマッピングを計算する。"""

    i64 = _as_float64_array("i_bln", i_bln, ndim=2)
    j64 = _as_float64_array("j_bln", j_bln, ndim=2)
    k64 = _as_float64_array("k_bln", k_bln, ndim=2)
    if i64.shape != j64.shape or i64.shape != k64.shape:
        raise ValueError("i_bln, j_bln, k_bln のshapeは一致する必要があります。")
    nx_bln, lx_bln = i64.shape

    lmin32 = _as_int32_array("lmin_bln", lmin_bln, ndim=1)
    lmax32 = _as_int32_array("lmax_bln", lmax_bln, ndim=1)
    if lmin32.shape != (nx_bln,) or lmax32.shape != (nx_bln,):
        raise ValueError("lmin_bln, lmax_bln のshapeは (nx_bln,) である必要があります。")

    if not (1 <= lcen_bln <= lx_bln):
        raise ValueError("lcen_bln は 1 以上 lx_bln 以下である必要があります。")

    try:
        from . import _bbtobln
    except ImportError as exc:
        raise ImportError(
            "Fortran拡張の読み込みに失敗しました。`pip install -e .[dev]` を実行してください。"
        ) from exc

    i_obs = np.zeros(nx_bln, dtype=np.float64, order="F")
    j_obs = np.zeros(nx_bln, dtype=np.float64, order="F")
    dk_obs = np.zeros(nx_bln, dtype=np.float64, order="F")

    _bbtobln.blntobmap(
        i_obs,
        j_obs,
        dk_obs,
        np.asfortranarray(i64),
        np.asfortranarray(j64),
        np.asfortranarray(k64),
        np.asfortranarray(lmin32),
        np.asfortranarray(lmax32),
        int(lcen_bln),
        int(k_obs),
    )
    return ObservationMapResult(i_obs=i_obs, j_obs=j_obs, dk_obs=dk_obs)


def _as_float64_array(name: str, value: Any, ndim: int) -> NDArray[np.float64]:
    arr = np.asarray(value, dtype=np.float64)
    if arr.ndim != ndim:
        raise ValueError(f"{name} は{ndim}次元配列である必要があります。")
    return arr


def _as_int32_array(name: str, value: Any, ndim: int) -> NDArray[np.int32]:
    arr = np.asarray(value, dtype=np.int32)
    if arr.ndim != ndim:
        raise ValueError(f"{name} は{ndim}次元配列である必要があります。")
    return arr
