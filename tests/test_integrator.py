from __future__ import annotations

import importlib

import numpy as np
import pytest

from wind3d_field_lines import trace_field_lines


def _build_uniform_case(dtype: type[np.floating]) -> tuple[np.ndarray, ...]:
    ix, jx, kx = 4, 4, 6
    bx = np.zeros((ix, jx, kx), dtype=dtype)
    by = np.zeros((ix, jx, kx), dtype=dtype)
    bz = np.ones((ix, jx, kx), dtype=dtype)
    dx = np.ones(kx, dtype=dtype)
    dy = np.ones(kx, dtype=dtype)
    dz = np.ones(kx, dtype=dtype)
    icen = np.array([2.0, 3.0], dtype=dtype)
    jcen = np.array([2.0, 3.0], dtype=dtype)
    kcen = np.array([3.0, 3.0], dtype=dtype)
    return bx, by, bz, dx, dy, dz, icen, jcen, kcen


@pytest.mark.skipif(
    importlib.util.find_spec("wind3d_field_lines._bbtobln") is None,
    reason="Fortran extension is not built.",
)
def test_trace_uniform_vertical_field() -> None:
    bx, by, bz, dx, dy, dz, icen, jcen, kcen = _build_uniform_case(np.float64)

    result = trace_field_lines(
        bx=bx,
        by=by,
        bz=bz,
        dx=dx,
        dy=dy,
        dz=dz,
        icen_bln=icen,
        jcen_bln=jcen,
        kcen_bln=kcen,
        lcen_bln=3,
        lx_bln=6,
        margin=1,
        nsubstepx=2,
    )

    np.testing.assert_allclose(result.i[:, 2], icen)
    np.testing.assert_allclose(result.j[:, 2], jcen)
    np.testing.assert_allclose(result.k[:, 2], kcen)
    assert np.all((1 <= result.lmin) & (result.lmin <= result.lx))
    assert np.all((1 <= result.lmax) & (result.lmax <= result.lx))


@pytest.mark.skipif(
    importlib.util.find_spec("wind3d_field_lines._bbtobln") is None,
    reason="Fortran extension is not built.",
)
def test_trace_boundary_clipping() -> None:
    bx, by, bz, dx, dy, dz, icen, jcen, _ = _build_uniform_case(np.float64)
    kcen = np.array([4.8, 4.8], dtype=np.float64)

    result = trace_field_lines(
        bx=bx,
        by=by,
        bz=bz,
        dx=dx,
        dy=dy,
        dz=dz,
        icen_bln=icen,
        jcen_bln=jcen,
        kcen_bln=kcen,
        lcen_bln=3,
        lx_bln=6,
        margin=1,
        nsubstepx=2,
    )

    assert np.all(result.lmax < result.lx)


@pytest.mark.skipif(
    importlib.util.find_spec("wind3d_field_lines._bbtobln") is None,
    reason="Fortran extension is not built.",
)
def test_float32_input_is_accepted() -> None:
    bx, by, bz, dx, dy, dz, icen, jcen, kcen = _build_uniform_case(np.float32)

    result = trace_field_lines(
        bx=bx,
        by=by,
        bz=bz,
        dx=dx,
        dy=dy,
        dz=dz,
        icen_bln=icen,
        jcen_bln=jcen,
        kcen_bln=kcen,
        lcen_bln=3,
        lx_bln=6,
        margin=1,
    )

    assert result.i.dtype == np.float64
    assert result.j.dtype == np.float64
    assert result.k.dtype == np.float64


def test_invalid_shape_raises_value_error() -> None:
    bx = np.zeros((2, 2, 2), dtype=np.float64)
    by = np.zeros((2, 2, 2), dtype=np.float64)
    bz = np.zeros((2, 2, 3), dtype=np.float64)

    with pytest.raises(ValueError, match="shape"):
        trace_field_lines(
            bx=bx,
            by=by,
            bz=bz,
            dx=np.ones(2),
            dy=np.ones(2),
            dz=np.ones(2),
            icen_bln=np.array([1.0]),
            jcen_bln=np.array([1.0]),
            kcen_bln=np.array([1.0]),
            lcen_bln=1,
            lx_bln=2,
            margin=0,
        )


def test_import_extension_module() -> None:
    spec = importlib.util.find_spec("wind3d_field_lines._bbtobln")
    if spec is None:
        pytest.skip("Fortran extension is not built.")
    mod = importlib.import_module("wind3d_field_lines._bbtobln")
    assert hasattr(mod, "bbtobln")
