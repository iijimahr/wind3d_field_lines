"""Tests for compute_potential_field and _thomas_solve."""

from __future__ import annotations

import numpy as np
import pytest

from wind3d_field_lines import compute_potential_field
from wind3d_field_lines.potential_field import _thomas_solve

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _cartesian_h(kx: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return unit scale factors (Cartesian coordinates)."""
    ones = np.ones(kx)
    return ones.copy(), ones.copy(), ones.copy()


def _make_grid(ix: int, jx: int, kx: int, lzt: float):
    """Return (dxi, det, lzt, kx, hxi, het, hzt) for Cartesian domain."""
    dxi = 1.0
    det = 1.0
    hxi, het, hzt = _cartesian_h(kx)
    return dxi, det, lzt, kx, hxi, het, hzt


# ---------------------------------------------------------------------------
# Correctness tests
# ---------------------------------------------------------------------------


class TestCartesianUniformField:
    """Constant surface field => bzt = const, bxi = bet = 0 everywhere."""

    def test_bzt_equals_boundary(self):
        ix, jx, kx = 8, 8, 16
        lzt = 10.0
        dxi, det, lzt, kx, hxi, het, hzt = _make_grid(ix, jx, kx, lzt)
        bzt_bottom = np.ones((ix, jx)) * 3.7

        bxi, bet, bzt = compute_potential_field(
            bzt_bottom=bzt_bottom,
            dxi=dxi,
            det=det,
            lzt=lzt,
            kx=kx,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        np.testing.assert_allclose(bzt, 3.7, atol=1e-10)

    def test_bxi_bet_zero(self):
        ix, jx, kx = 8, 8, 16
        lzt = 10.0
        dxi, det, lzt, kx, hxi, het, hzt = _make_grid(ix, jx, kx, lzt)
        bzt_bottom = np.ones((ix, jx)) * 2.5

        bxi, bet, bzt = compute_potential_field(
            bzt_bottom=bzt_bottom,
            dxi=dxi,
            det=det,
            lzt=lzt,
            kx=kx,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        np.testing.assert_allclose(bxi, 0.0, atol=1e-10)
        np.testing.assert_allclose(bet, 0.0, atol=1e-10)

    def test_output_shape(self):
        ix, jx, kx = 6, 10, 12
        lzt = 5.0
        dxi, det, lzt, kx, hxi, het, hzt = _make_grid(ix, jx, kx, lzt)
        bzt_bottom = np.ones((ix, jx))

        bxi, bet, bzt = compute_potential_field(
            bzt_bottom=bzt_bottom,
            dxi=dxi,
            det=det,
            lzt=lzt,
            kx=kx,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        assert bxi.shape == (ix, jx, kx)
        assert bet.shape == (ix, jx, kx)
        assert bzt.shape == (ix, jx, kx)


class TestCartesianAnalyticalSingleMode:
    """Single Fourier mode in xi1 only; compare bzt with the sinh analytical solution.

    For Cartesian coordinates with bzt_bottom = sin(k1*xi1) and all other
    wavenumbers zero:
        f(xi3) = sinh(kappa*(lzt - xi3)) / (kappa * cosh(kappa*lzt))
        bzt(xi1, xi3) = sin(k1*xi1) * cosh(kappa*(lzt-xi3)) / cosh(kappa*lzt)
    where kappa = |k1|.
    """

    @pytest.mark.parametrize("kx", [32, 64])
    def test_bzt_matches_analytical(self, kx: int):
        ix = 32
        lzt = 20.0
        dxi = 1.0
        det = 1.0
        hxi, het, hzt = _cartesian_h(kx)

        dz = lzt / (kx - 0.5)
        xi3 = np.arange(kx) * dz  # (kx,)
        k1_phys = 2.0 * np.pi / (ix * dxi)  # fundamental mode
        xi1 = np.arange(ix) * dxi  # (ix,)

        bzt_bottom = np.sin(k1_phys * xi1)[:, None]  # (ix, 1)

        bxi, bet, bzt = compute_potential_field(
            bzt_bottom=bzt_bottom,
            dxi=dxi,
            det=det,
            lzt=lzt,
            kx=kx,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        kappa = k1_phys
        bzt_analytical = (
            np.sin(k1_phys * xi1)[:, None, None]
            * (np.cosh(kappa * (lzt - xi3)) / np.cosh(kappa * lzt))[None, None, :]
        )  # (ix, 1, kx)

        # Second-order FD scheme => error ~ dz^2
        np.testing.assert_allclose(bzt, bzt_analytical, atol=5e-4)


class TestFluxConservation:
    """Horizontally integrated bzt must equal the surface integral at all heights."""

    def test_flux_conserved(self):
        rng = np.random.default_rng(42)
        ix, jx, kx = 16, 16, 20
        lzt = 8.0
        dxi, det, _, _, hxi, het, hzt = _make_grid(ix, jx, kx, lzt)

        bzt_bottom = rng.standard_normal((ix, jx))

        bxi, bet, bzt = compute_potential_field(
            bzt_bottom=bzt_bottom,
            dxi=dxi,
            det=det,
            lzt=lzt,
            kx=kx,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        flux_bottom = bzt_bottom.sum()
        for k in range(kx):
            np.testing.assert_allclose(
                bzt[:, :, k].sum(),
                flux_bottom,
                rtol=1e-10,
                err_msg=f"Flux not conserved at level k={k}",
            )


class TestUpperBoundaryCondition:
    """Horizontal field (bxi, bet) should vanish at xi3 = lzt.

    The discrete scheme enforces f = 0 at the half-integer point
    xi3 = lzt.  The last interior grid point is at xi3 = (kx-1)*dz,
    which differs from lzt by dz/2.  For large kx the horizontal field
    at the last interior level is O(kappa * dz / 2) smaller than the
    surface amplitude.
    """

    def test_bxi_bet_small_at_top(self):
        rng = np.random.default_rng(7)
        ix, jx, kx = 16, 16, 64
        lzt = 20.0
        dxi, det, _, _, hxi, het, hzt = _make_grid(ix, jx, kx, lzt)
        bzt_bottom = rng.standard_normal((ix, jx))

        bxi, bet, bzt = compute_potential_field(
            bzt_bottom=bzt_bottom,
            dxi=dxi,
            det=det,
            lzt=lzt,
            kx=kx,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        amp = np.abs(bzt_bottom).max()
        assert np.abs(bxi[:, :, -1]).max() < 0.1 * amp
        assert np.abs(bet[:, :, -1]).max() < 0.1 * amp


class TestPotentialCondition:
    """For a potential field curl B = 0 (within numerical discretisation error)."""

    def test_curl_b_small(self):
        ix, jx, kx = 16, 16, 32
        lzt = 10.0
        dxi = 1.0
        det = 1.0
        hxi, het, hzt = _cartesian_h(kx)

        rng = np.random.default_rng(99)
        bzt_bottom = rng.standard_normal((ix, jx))

        bxi, bet, bzt = compute_potential_field(
            bzt_bottom=bzt_bottom,
            dxi=dxi,
            det=det,
            lzt=lzt,
            kx=kx,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        # Estimate (curl B)_z = d(bet)/dxi1 - d(bxi)/dxi2 in interior
        # Use FFT for horizontal derivatives (exact)
        bxi_fft = np.fft.fft2(bxi, axes=(0, 1))
        bet_fft = np.fft.fft2(bet, axes=(0, 1))
        k1 = 2.0 * np.pi * np.fft.fftfreq(ix) / dxi
        k2 = 2.0 * np.pi * np.fft.fftfreq(jx) / det
        curl_z = np.real(
            np.fft.ifft2(
                1j * k1[:, None, None] * bet_fft - 1j * k2[None, :, None] * bxi_fft,
                axes=(0, 1),
            )
        )  # (ix, jx, kx)

        amp = np.abs(bzt_bottom).max()
        # Curl should be negligible compared to field amplitude
        assert np.abs(curl_z).max() < 1e-10 * amp


# ---------------------------------------------------------------------------
# Input validation tests
# ---------------------------------------------------------------------------


class TestInputValidation:
    def _default_args(self):
        ix, jx, kx = 4, 4, 8
        return dict(
            bzt_bottom=np.ones((ix, jx)),
            dxi=1.0,
            det=1.0,
            lzt=5.0,
            kx=kx,
            hxi=np.ones(kx),
            het=np.ones(kx),
            hzt=np.ones(kx),
        )

    def test_bzt_bottom_not_2d(self):
        args = self._default_args()
        args["bzt_bottom"] = np.ones(4)
        with pytest.raises(ValueError, match="2D"):
            compute_potential_field(**args)

    def test_wrong_scale_factor_shape(self):
        args = self._default_args()
        args["hxi"] = np.ones(99)
        with pytest.raises(ValueError, match="shape"):
            compute_potential_field(**args)

    def test_non_positive_scale_factor(self):
        args = self._default_args()
        args["hzt"] = np.full(8, -1.0)
        with pytest.raises(ValueError, match="positive"):
            compute_potential_field(**args)

    def test_non_positive_dxi(self):
        args = self._default_args()
        args["dxi"] = 0.0
        with pytest.raises(ValueError, match="positive"):
            compute_potential_field(**args)

    def test_non_positive_lzt(self):
        args = self._default_args()
        args["lzt"] = -1.0
        with pytest.raises(ValueError, match="positive"):
            compute_potential_field(**args)

    def test_kx_too_small(self):
        args = self._default_args()
        args["kx"] = 1
        args["hxi"] = np.ones(1)
        args["het"] = np.ones(1)
        args["hzt"] = np.ones(1)
        with pytest.raises(ValueError, match="kx"):
            compute_potential_field(**args)

    def test_float32_input_accepted(self):
        args = self._default_args()
        args["bzt_bottom"] = args["bzt_bottom"].astype(np.float32)
        args["hxi"] = args["hxi"].astype(np.float32)
        bxi, bet, bzt = compute_potential_field(**args)
        assert bxi.dtype == np.float64


# ---------------------------------------------------------------------------
# _thomas_solve unit tests
# ---------------------------------------------------------------------------


class TestThomasSolve:
    """Unit tests for the internal TDMA solver."""

    def _solve(self, lower, diag, upper, rhs):
        """Thin wrapper that converts lists to arrays."""
        return _thomas_solve(
            lower=np.array(lower, dtype=np.float64),
            diag=np.array(diag, dtype=np.float64),
            upper=np.array(upper, dtype=np.float64),
            rhs=np.array(rhs, dtype=np.float64),
        )

    # --- correctness against numpy.linalg.solve ---

    def test_n2_system(self):
        # [2 -1] [x0]   [1]
        # [-1 3] [x1] = [2]
        lower = [0.0, -1.0]
        diag = [2.0, 3.0]
        upper = [-1.0, 0.0]
        rhs = [1.0, 2.0]

        x = self._solve(lower, diag, upper, rhs)

        A = np.array([[2, -1], [-1, 3]], dtype=np.float64)
        expected = np.linalg.solve(A, [1.0, 2.0])
        np.testing.assert_allclose(x, expected, atol=1e-14)

    def test_n5_system(self):
        rng = np.random.default_rng(0)
        n = 5
        # Build a diagonally dominant tridiagonal matrix.
        lower = np.concatenate([[0.0], rng.uniform(-1, 0, n - 1)])
        upper = np.concatenate([rng.uniform(-1, 0, n - 1), [0.0]])
        diag = -(np.abs(lower) + np.abs(upper)) - 1.0  # strictly dominant
        rhs = rng.standard_normal(n)

        x = _thomas_solve(lower, diag, upper, rhs)

        A = np.diag(diag) + np.diag(lower[1:], -1) + np.diag(upper[:-1], 1)
        expected = np.linalg.solve(A, rhs)
        np.testing.assert_allclose(x, expected, atol=1e-12)

    def test_identity_like_system(self):
        # Diagonal-only system: solution is rhs / diag.
        n = 8
        diag = np.full(n, 3.0)
        lower = np.zeros(n)
        upper = np.zeros(n)
        rhs = np.arange(1.0, n + 1)

        x = _thomas_solve(lower, diag, upper, rhs)

        np.testing.assert_allclose(x, rhs / 3.0, atol=1e-14)

    # --- batched (multi-dimensional) RHS ---

    def test_batched_2d_rhs(self):
        # Same matrix, multiple RHS vectors stacked in axis 0.
        n = 4
        lower = np.array([0.0, -1.0, -1.0, -1.0])
        diag = np.array([-3.0, -3.0, -3.0, -3.0])
        upper = np.array([-1.0, -1.0, -1.0, 0.0])
        A = np.diag(diag) + np.diag(lower[1:], -1) + np.diag(upper[:-1], 1)

        rng = np.random.default_rng(1)
        batch = 5
        rhs_batch = rng.standard_normal((batch, n))  # (5, n)

        x_batch = _thomas_solve(
            lower=lower,
            diag=np.broadcast_to(diag, (batch, n)).copy(),
            upper=upper,
            rhs=rhs_batch,
        )

        for i in range(batch):
            expected = np.linalg.solve(A, rhs_batch[i])
            np.testing.assert_allclose(x_batch[i], expected, atol=1e-12)

    def test_batched_3d_rhs(self):
        # Three-dimensional batch (m, n, system_size).
        n = 6
        lower = np.array([0.0, -0.5, -0.5, -0.5, -0.5, -0.5])
        diag = np.full(n, -2.0)
        upper = np.array([-0.5, -0.5, -0.5, -0.5, -0.5, 0.0])
        A = np.diag(diag) + np.diag(lower[1:], -1) + np.diag(upper[:-1], 1)

        rng = np.random.default_rng(2)
        m1, m2 = 3, 4
        rhs_3d = rng.standard_normal((m1, m2, n))

        x_3d = _thomas_solve(
            lower=lower,
            diag=np.broadcast_to(diag, (m1, m2, n)).copy(),
            upper=upper,
            rhs=rhs_3d,
        )

        assert x_3d.shape == (m1, m2, n)
        for i in range(m1):
            for j in range(m2):
                expected = np.linalg.solve(A, rhs_3d[i, j])
                np.testing.assert_allclose(x_3d[i, j], expected, atol=1e-12)

    # --- output properties ---

    def test_output_shape_matches_rhs(self):
        n = 7
        diag = np.full((3, 5, n), -2.0)
        lower = np.zeros(n)
        upper = np.zeros(n)
        rhs = np.ones((3, 5, n))

        x = _thomas_solve(lower, diag, upper, rhs)
        assert x.shape == (3, 5, n)

    def test_output_dtype_is_float64(self):
        diag = np.array([-2.0, -2.0])
        lower = np.zeros(2)
        upper = np.zeros(2)
        rhs = np.ones(2)

        x = _thomas_solve(lower, diag, upper, rhs)
        assert x.dtype == np.float64

    def test_does_not_modify_input_rhs(self):
        diag = np.array([-2.0, -3.0, -2.0])
        lower = np.array([0.0, -0.5, -0.5])
        upper = np.array([-0.5, -0.5, 0.0])
        rhs = np.array([1.0, 2.0, 3.0])
        rhs_copy = rhs.copy()

        _thomas_solve(lower, diag, upper, rhs)

        np.testing.assert_array_equal(rhs, rhs_copy)
