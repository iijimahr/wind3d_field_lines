from __future__ import annotations

import importlib

import numpy as np
import pytest

from wind3d_field_lines import CurvilinearFieldLineResult, trace_field_lines_curvilinear

_HAS_FORTRAN = importlib.util.find_spec("wind3d_field_lines._bbtobln") is not None

_skip_no_fortran = pytest.mark.skipif(
    not _HAS_FORTRAN,
    reason="Fortran extension is not built.",
)


def _build_uniform_curvilinear_case(
    dtype: type[np.floating] = np.float64,
) -> tuple[np.ndarray, ...]:
    ix, jx, kx = 4, 4, 6
    bxi = np.zeros((ix, jx, kx), dtype=dtype)
    bet = np.zeros((ix, jx, kx), dtype=dtype)
    bzt = np.ones((ix, jx, kx), dtype=dtype)
    hxi = np.ones(kx, dtype=dtype)
    het = np.ones(kx, dtype=dtype)
    hzt = np.ones(kx, dtype=dtype)
    icen = np.array([2.0, 3.0], dtype=dtype)
    jcen = np.array([2.0, 3.0], dtype=dtype)
    kcen = np.array([3.0, 3.0], dtype=dtype)
    return bxi, bet, bzt, hxi, het, hzt, icen, jcen, kcen


@_skip_no_fortran
def test_cartesian_identity() -> None:
    """With h=1 and d=1 the physical coordinates match the grid-index result."""
    bxi, bet, bzt, hxi, het, hzt, icen, jcen, kcen = _build_uniform_curvilinear_case()

    result = trace_field_lines_curvilinear(
        bxi=bxi,
        bet=bet,
        bzt=bzt,
        dxi=1.0,
        det=1.0,
        dzt=1.0,
        hxi=hxi,
        het=het,
        hzt=hzt,
        seed_i=icen,
        seed_j=jcen,
        seed_k=kcen,
        line_center=3,
        line_length=6,
        margin=1,
        n_substeps=2,
    )

    # With d=1 the physical coordinate equals (index - 1).
    np.testing.assert_allclose(result.xi[:, 2], icen - 1.0)
    np.testing.assert_allclose(result.eta[:, 2], jcen - 1.0)
    np.testing.assert_allclose(result.zeta[:, 2], kcen - 1.0)


@_skip_no_fortran
def test_uniform_scale_factor_geometry() -> None:
    """Uniform scale factors of different magnitudes yield the same geometry."""
    bxi, bet, bzt, _, _, _, icen, jcen, kcen = _build_uniform_curvilinear_case()
    kx = bxi.shape[2]

    common_kwargs = dict(
        bxi=bxi,
        bet=bet,
        bzt=bzt,
        dxi=1.0,
        det=1.0,
        dzt=1.0,
        seed_i=icen,
        seed_j=jcen,
        seed_k=kcen,
        line_center=3,
        line_length=6,
        margin=1,
        n_substeps=2,
    )

    result_h1 = trace_field_lines_curvilinear(
        hxi=np.ones(kx),
        het=np.ones(kx),
        hzt=np.ones(kx),
        **common_kwargs,
    )
    result_h2 = trace_field_lines_curvilinear(
        hxi=np.full(kx, 2.0),
        het=np.full(kx, 2.0),
        hzt=np.full(kx, 2.0),
        **common_kwargs,
    )

    np.testing.assert_allclose(result_h1.xi, result_h2.xi, rtol=1e-6)
    np.testing.assert_allclose(result_h1.eta, result_h2.eta, rtol=1e-6)
    np.testing.assert_allclose(result_h1.zeta, result_h2.zeta, rtol=1e-6)


@_skip_no_fortran
def test_output_type_and_shape() -> None:
    """Result is a CurvilinearFieldLineResult with correct array shapes."""
    bxi, bet, bzt, hxi, het, hzt, icen, jcen, kcen = _build_uniform_curvilinear_case()

    result = trace_field_lines_curvilinear(
        bxi=bxi,
        bet=bet,
        bzt=bzt,
        dxi=1.0,
        det=1.0,
        dzt=1.0,
        hxi=hxi,
        het=het,
        hzt=hzt,
        seed_i=icen,
        seed_j=jcen,
        seed_k=kcen,
        line_center=3,
        line_length=6,
        margin=1,
    )

    assert isinstance(result, CurvilinearFieldLineResult)
    assert result.xi.shape == (2, 6)
    assert result.eta.shape == (2, 6)
    assert result.zeta.shape == (2, 6)
    assert result.lmin.shape == (2,)
    assert result.lmax.shape == (2,)


@_skip_no_fortran
def test_float32_input_accepted() -> None:
    """float32 inputs are cast to float64 and the result uses float64 arrays."""
    bxi, bet, bzt, hxi, het, hzt, icen, jcen, kcen = _build_uniform_curvilinear_case(
        dtype=np.float32
    )

    result = trace_field_lines_curvilinear(
        bxi=bxi,
        bet=bet,
        bzt=bzt,
        dxi=1.0,
        det=1.0,
        dzt=1.0,
        hxi=hxi,
        het=het,
        hzt=hzt,
        seed_i=icen,
        seed_j=jcen,
        seed_k=kcen,
        line_center=3,
        line_length=6,
        margin=1,
    )

    assert result.xi.dtype == np.float64
    assert result.eta.dtype == np.float64
    assert result.zeta.dtype == np.float64


@_skip_no_fortran
def test_varying_scale_factor_in_zeta() -> None:
    """A zeta-varying scale factor produces valid, finite output."""
    ix, jx, kx = 4, 4, 8
    bxi = np.zeros((ix, jx, kx))
    bet = np.zeros((ix, jx, kx))
    bzt = np.ones((ix, jx, kx))
    hxi = np.ones(kx)
    het = np.ones(kx)
    hzt = np.linspace(1.0, 2.0, kx)

    result = trace_field_lines_curvilinear(
        bxi=bxi,
        bet=bet,
        bzt=bzt,
        dxi=1.0,
        det=1.0,
        dzt=1.0,
        hxi=hxi,
        het=het,
        hzt=hzt,
        seed_i=np.array([2.0, 3.0]),
        seed_j=np.array([2.0, 3.0]),
        seed_k=np.array([4.0, 4.0]),
        line_center=5,
        line_length=10,
        margin=1,
    )

    assert np.all(np.isfinite(result.zeta))
    assert np.all((1 <= result.lmin) & (result.lmin <= result.line_length))
    assert np.all((1 <= result.lmax) & (result.lmax <= result.line_length))


def test_invalid_field_shape_raises() -> None:
    """Mismatched bxi/bet/bzt shapes raise ValueError."""
    bxi = np.zeros((2, 2, 3))
    bet = np.zeros((2, 2, 3))
    bzt = np.zeros((2, 2, 4))

    with pytest.raises(ValueError, match="shape"):
        trace_field_lines_curvilinear(
            bxi=bxi,
            bet=bet,
            bzt=bzt,
            dxi=1.0,
            det=1.0,
            dzt=1.0,
            hxi=np.ones(3),
            het=np.ones(3),
            hzt=np.ones(3),
            seed_i=np.array([1.0]),
            seed_j=np.array([1.0]),
            seed_k=np.array([1.0]),
            line_center=1,
            line_length=2,
            margin=0,
        )


def test_invalid_scale_factor_shape_raises() -> None:
    """hxi/het/hzt with wrong size raise ValueError."""
    bxi = np.zeros((2, 2, 4))  # kx=4
    bet = np.zeros((2, 2, 4))
    bzt = np.zeros((2, 2, 4))

    with pytest.raises(ValueError, match="shape"):
        trace_field_lines_curvilinear(
            bxi=bxi,
            bet=bet,
            bzt=bzt,
            dxi=1.0,
            det=1.0,
            dzt=1.0,
            hxi=np.ones(3),  # wrong size
            het=np.ones(4),
            hzt=np.ones(4),
            seed_i=np.array([1.0]),
            seed_j=np.array([1.0]),
            seed_k=np.array([1.0]),
            line_center=1,
            line_length=2,
            margin=0,
        )


def test_non_positive_grid_spacing_raises() -> None:
    """dxi <= 0 raises ValueError."""
    bxi = np.zeros((2, 2, 3))
    bet = np.zeros((2, 2, 3))
    bzt = np.zeros((2, 2, 3))

    with pytest.raises(ValueError, match="positive"):
        trace_field_lines_curvilinear(
            bxi=bxi,
            bet=bet,
            bzt=bzt,
            dxi=0.0,
            det=1.0,
            dzt=1.0,
            hxi=np.ones(3),
            het=np.ones(3),
            hzt=np.ones(3),
            seed_i=np.array([1.0]),
            seed_j=np.array([1.0]),
            seed_k=np.array([1.0]),
            line_center=1,
            line_length=2,
            margin=0,
        )


def test_non_positive_scale_factor_raises() -> None:
    """Scale factors with non-positive elements raise ValueError."""
    bxi = np.zeros((2, 2, 3))
    bet = np.zeros((2, 2, 3))
    bzt = np.zeros((2, 2, 3))

    with pytest.raises(ValueError, match="positive"):
        trace_field_lines_curvilinear(
            bxi=bxi,
            bet=bet,
            bzt=bzt,
            dxi=1.0,
            det=1.0,
            dzt=1.0,
            hxi=np.array([-1.0, 1.0, 1.0]),  # negative element
            het=np.ones(3),
            hzt=np.ones(3),
            seed_i=np.array([1.0]),
            seed_j=np.array([1.0]),
            seed_k=np.array([1.0]),
            line_center=1,
            line_length=2,
            margin=0,
        )


def test_positional_arguments_rejected() -> None:
    bxi = np.zeros((2, 2, 3))
    bet = np.zeros((2, 2, 3))
    bzt = np.zeros((2, 2, 3))

    with pytest.raises(TypeError):
        trace_field_lines_curvilinear(  # type: ignore[misc]
            bxi,
            bet,
            bzt,
            1.0,
            1.0,
            1.0,
            np.ones(3),
            np.ones(3),
            np.ones(3),
            np.array([1.0]),
            np.array([1.0]),
            np.array([1.0]),
            1,
            2,
        )


def test_legacy_keyword_arguments_rejected() -> None:
    bxi = np.zeros((2, 2, 3))
    bet = np.zeros((2, 2, 3))
    bzt = np.zeros((2, 2, 3))

    with pytest.raises(TypeError):
        legacy_kwargs = {
            "b" + "_xi": bxi,
            "b" + "_et": bet,
            "b" + "_zt": bzt,
        }
        trace_field_lines_curvilinear(
            **legacy_kwargs,  # type: ignore[arg-type]
            dxi=1.0,
            det=1.0,
            dzt=1.0,
            hxi=np.ones(3),
            het=np.ones(3),
            hzt=np.ones(3),
            seed_i=np.array([1.0]),
            seed_j=np.array([1.0]),
            seed_k=np.array([1.0]),
            line_center=1,
            line_length=2,
        )
