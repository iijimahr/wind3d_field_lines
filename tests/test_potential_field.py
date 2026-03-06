"""Tests for compute_potential_field."""

from __future__ import annotations

import numpy as np
import pytest

from wind3d_field_lines import compute_potential_field


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
            dxi=dxi, det=det, l3=l3, n3=n3,
            hxi=hxi, het=het, hzt=hzt,
        )

        np.testing.assert_allclose(b3, 3.7, atol=1e-10)

    def test_b1_b2_zero(self):
        n1, n2, n3 = 8, 8, 16
        l3 = 10.0
        dxi, det, l3, n3, hxi, het, hzt = _make_grid(n1, n2, n3, l3)
        b3_bottom = np.ones((n1, n2)) * 2.5

        b1, b2, b3 = compute_potential_field(
            b3_bottom=b3_bottom,
            dxi=dxi, det=det, l3=l3, n3=n3,
            hxi=hxi, het=het, hzt=hzt,
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
            dxi=dxi, det=det, l3=l3, n3=n3,
            hxi=hxi, het=het, hzt=hzt,
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
        n1, n2 = 32, 1
        l3 = 20.0
        dxi = 1.0
        det = 1.0
        hxi, het, hzt = _cartesian_h(n3)

        dz = l3 / (n3 - 0.5)
        xi3 = np.arange(n3) * dz                            # (n3,)
        k1_phys = 2.0 * np.pi / (n1 * dxi)                 # fundamental mode
        xi1 = np.arange(n1) * dxi                           # (n1,)

        b3_bottom = np.sin(k1_phys * xi1)[:, None]          # (n1, 1)

        b1, b2, b3 = compute_potential_field(
            b3_bottom=b3_bottom,
            dxi=dxi, det=det, l3=l3, n3=n3,
            hxi=hxi, het=het, hzt=hzt,
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
            dxi=dxi, det=det, l3=l3, n3=n3,
            hxi=hxi, het=het, hzt=hzt,
        )

        flux_bottom = b3_bottom.sum()
        for k in range(n3):
            np.testing.assert_allclose(
                b3[:, :, k].sum(), flux_bottom, rtol=1e-10,
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
            dxi=dxi, det=det, l3=l3, n3=n3,
            hxi=hxi, het=het, hzt=hzt,
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
        dz = l3 / (n3 - 0.5)

        rng = np.random.default_rng(99)
        b3_bottom = rng.standard_normal((n1, n2))

        b1, b2, b3 = compute_potential_field(
            b3_bottom=b3_bottom,
            dxi=dxi, det=det, l3=l3, n3=n3,
            hxi=hxi, het=het, hzt=hzt,
        )

        # Estimate (curl B)_z = dB2/dxi1 - dB1/dxi2 in interior
        # Use FFT for horizontal derivatives (exact)
        B1_fft = np.fft.fft2(b1, axes=(0, 1))
        B2_fft = np.fft.fft2(b2, axes=(0, 1))
        k1 = 2.0 * np.pi * np.fft.fftfreq(n1) / dxi
        k2 = 2.0 * np.pi * np.fft.fftfreq(n2) / det
        curl_z = np.real(np.fft.ifft2(
            1j * k1[:, None, None] * B2_fft - 1j * k2[None, :, None] * B1_fft,
            axes=(0, 1),
        ))  # (n1, n2, n3)

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
            dxi=1.0, det=1.0, l3=5.0, n3=n3,
            hxi=np.ones(n3), het=np.ones(n3), hzt=np.ones(n3),
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
