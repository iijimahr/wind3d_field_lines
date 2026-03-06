"""Tests for compute_potential_field and _thomas_solve."""

from __future__ import annotations

import numpy as np
import pytest

from wind3d_field_lines import compute_potential_field
from wind3d_field_lines.potential_field import _thomas_solve

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _cartesian_h(n3: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return unit scale factors (Cartesian coordinates)."""
    ones = np.ones(n3)
    return ones.copy(), ones.copy(), ones.copy()


def _make_grid(n1: int, n2: int, n3: int, l3: float):
    """Return (dxi, det, l3, n3, hxi, het, hzt) for Cartesian domain."""
    dxi = 1.0
    det = 1.0
    hxi, het, hzt = _cartesian_h(n3)
    return dxi, det, l3, n3, hxi, het, hzt


# ---------------------------------------------------------------------------
# Correctness tests
# ---------------------------------------------------------------------------


class TestCartesianUniformField:
    """Constant surface field => B3 = const, B1 = B2 = 0 everywhere."""

    def test_b3_equals_boundary(self):
        n1, n2, n3 = 8, 8, 16
        l3 = 10.0
        dxi, det, l3, n3, hxi, het, hzt = _make_grid(n1, n2, n3, l3)
        b3_bottom = np.ones((n1, n2)) * 3.7

        b1, b2, b3 = compute_potential_field(
            b3_bottom=b3_bottom,
            dxi=dxi,
            det=det,
            l3=l3,
            n3=n3,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        np.testing.assert_allclose(b3, 3.7, atol=1e-10)

    def test_b1_b2_zero(self):
        n1, n2, n3 = 8, 8, 16
        l3 = 10.0
        dxi, det, l3, n3, hxi, het, hzt = _make_grid(n1, n2, n3, l3)
        b3_bottom = np.ones((n1, n2)) * 2.5

        b1, b2, b3 = compute_potential_field(
            b3_bottom=b3_bottom,
            dxi=dxi,
            det=det,
            l3=l3,
            n3=n3,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        np.testing.assert_allclose(b1, 0.0, atol=1e-10)
        np.testing.assert_allclose(b2, 0.0, atol=1e-10)

    def test_output_shape(self):
        n1, n2, n3 = 6, 10, 12
        l3 = 5.0
        dxi, det, l3, n3, hxi, het, hzt = _make_grid(n1, n2, n3, l3)
        b3_bottom = np.ones((n1, n2))

        b1, b2, b3 = compute_potential_field(
            b3_bottom=b3_bottom,
            dxi=dxi,
            det=det,
            l3=l3,
            n3=n3,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        assert b1.shape == (n1, n2, n3)
        assert b2.shape == (n1, n2, n3)
        assert b3.shape == (n1, n2, n3)


class TestCartesianAnalyticalSingleMode:
    """Single Fourier mode in xi1 only; compare B3 with the sinh analytical solution.

    For Cartesian coordinates with b3_bottom = sin(k1*xi1) and all other
    wavenumbers zero:
        f(xi3) = sinh(kappa*(l3 - xi3)) / (kappa * cosh(kappa*l3))
        B3(xi1, xi3) = sin(k1*xi1) * cosh(kappa*(l3-xi3)) / cosh(kappa*l3)
    where kappa = |k1|.
    """

    @pytest.mark.parametrize("n3", [32, 64])
    def test_b3_matches_analytical(self, n3: int):
        n1 = 32
        l3 = 20.0
        dxi = 1.0
        det = 1.0
        hxi, het, hzt = _cartesian_h(n3)

        dz = l3 / (n3 - 0.5)
        xi3 = np.arange(n3) * dz  # (n3,)
        k1_phys = 2.0 * np.pi / (n1 * dxi)  # fundamental mode
        xi1 = np.arange(n1) * dxi  # (n1,)

        b3_bottom = np.sin(k1_phys * xi1)[:, None]  # (n1, 1)

        b1, b2, b3 = compute_potential_field(
            b3_bottom=b3_bottom,
            dxi=dxi,
            det=det,
            l3=l3,
            n3=n3,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        kappa = k1_phys
        b3_analytical = (
            np.sin(k1_phys * xi1)[:, None, None]
            * (np.cosh(kappa * (l3 - xi3)) / np.cosh(kappa * l3))[None, None, :]
        )  # (n1, 1, n3)

        # Second-order FD scheme => error ~ dz^2
        np.testing.assert_allclose(b3, b3_analytical, atol=5e-4)


class TestFluxConservation:
    """Horizontally integrated B3 must equal the surface integral at all heights."""

    def test_flux_conserved(self):
        rng = np.random.default_rng(42)
        n1, n2, n3 = 16, 16, 20
        l3 = 8.0
        dxi, det, _, _, hxi, het, hzt = _make_grid(n1, n2, n3, l3)

        b3_bottom = rng.standard_normal((n1, n2))

        b1, b2, b3 = compute_potential_field(
            b3_bottom=b3_bottom,
            dxi=dxi,
            det=det,
            l3=l3,
            n3=n3,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        flux_bottom = b3_bottom.sum()
        for k in range(n3):
            np.testing.assert_allclose(
                b3[:, :, k].sum(),
                flux_bottom,
                rtol=1e-10,
                err_msg=f"Flux not conserved at level k={k}",
            )


class TestUpperBoundaryCondition:
    """Horizontal field (B1, B2) should vanish at xi3 = l3.

    The discrete scheme enforces f = 0 at the half-integer point
    xi3 = l3.  The last interior grid point is at xi3 = (n3-1)*dz,
    which differs from l3 by dz/2.  For large n3 the horizontal field
    at the last interior level is O(kappa * dz / 2) smaller than the
    surface amplitude.
    """

    def test_b1_b2_small_at_top(self):
        rng = np.random.default_rng(7)
        n1, n2, n3 = 16, 16, 64
        l3 = 20.0
        dxi, det, _, _, hxi, het, hzt = _make_grid(n1, n2, n3, l3)
        b3_bottom = rng.standard_normal((n1, n2))

        b1, b2, b3 = compute_potential_field(
            b3_bottom=b3_bottom,
            dxi=dxi,
            det=det,
            l3=l3,
            n3=n3,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        amp = np.abs(b3_bottom).max()
        assert np.abs(b1[:, :, -1]).max() < 0.1 * amp
        assert np.abs(b2[:, :, -1]).max() < 0.1 * amp


class TestPotentialCondition:
    """For a potential field curl B = 0 (within numerical discretisation error)."""

    def test_curl_b_small(self):
        n1, n2, n3 = 16, 16, 32
        l3 = 10.0
        dxi = 1.0
        det = 1.0
        hxi, het, hzt = _cartesian_h(n3)

        rng = np.random.default_rng(99)
        b3_bottom = rng.standard_normal((n1, n2))

        b1, b2, b3 = compute_potential_field(
            b3_bottom=b3_bottom,
            dxi=dxi,
            det=det,
            l3=l3,
            n3=n3,
            hxi=hxi,
            het=het,
            hzt=hzt,
        )

        # Estimate (curl B)_z = dB2/dxi1 - dB1/dxi2 in interior
        # Use FFT for horizontal derivatives (exact)
        B1_fft = np.fft.fft2(b1, axes=(0, 1))
        B2_fft = np.fft.fft2(b2, axes=(0, 1))
        k1 = 2.0 * np.pi * np.fft.fftfreq(n1) / dxi
        k2 = 2.0 * np.pi * np.fft.fftfreq(n2) / det
        curl_z = np.real(
            np.fft.ifft2(
                1j * k1[:, None, None] * B2_fft - 1j * k2[None, :, None] * B1_fft,
                axes=(0, 1),
            )
        )  # (n1, n2, n3)

        amp = np.abs(b3_bottom).max()
        # Curl should be negligible compared to field amplitude
        assert np.abs(curl_z).max() < 1e-10 * amp


# ---------------------------------------------------------------------------
# Input validation tests
# ---------------------------------------------------------------------------


class TestInputValidation:
    def _default_args(self):
        n1, n2, n3 = 4, 4, 8
        return dict(
            b3_bottom=np.ones((n1, n2)),
            dxi=1.0,
            det=1.0,
            l3=5.0,
            n3=n3,
            hxi=np.ones(n3),
            het=np.ones(n3),
            hzt=np.ones(n3),
        )

    def test_b3_bottom_not_2d(self):
        args = self._default_args()
        args["b3_bottom"] = np.ones(4)
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

    def test_non_positive_l3(self):
        args = self._default_args()
        args["l3"] = -1.0
        with pytest.raises(ValueError, match="positive"):
            compute_potential_field(**args)

    def test_n3_too_small(self):
        args = self._default_args()
        args["n3"] = 1
        args["hxi"] = np.ones(1)
        args["het"] = np.ones(1)
        args["hzt"] = np.ones(1)
        with pytest.raises(ValueError, match="n3"):
            compute_potential_field(**args)

    def test_float32_input_accepted(self):
        args = self._default_args()
        args["b3_bottom"] = args["b3_bottom"].astype(np.float32)
        args["hxi"] = args["hxi"].astype(np.float32)
        b1, b2, b3 = compute_potential_field(**args)
        assert b1.dtype == np.float64


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
