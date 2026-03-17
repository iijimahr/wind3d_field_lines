from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray


def compute_potential_field(
    *,
    bzt_bottom: NDArray[np.floating[Any]],
    dxi: float,
    det: float,
    lzt: float,
    kx: int,
    hxi: NDArray[np.floating[Any]],
    het: NDArray[np.floating[Any]],
    hzt: NDArray[np.floating[Any]],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Compute the potential magnetic field from the surface normal component.

    Solves the Laplace equation ``nabla^2 Psi = 0`` for the scalar magnetic
    potential ``Psi``, with ``B = -nabla Psi``.  The computation is performed
    in an orthogonal curvilinear coordinate system ``(xi, eta, zeta)``
    where the scale factors ``h`` depend only on ``zeta``.

    The algorithm uses a 2-D FFT in the horizontal ``(xi, eta)`` directions
    and solves a 1-D boundary-value ODE for each horizontal wavenumber pair
    using the Thomas (TDMA) algorithm.  See
    ``SPEC/potential_field_extrapolation.md`` for the full derivation.

    Parameters
    ----------
    bzt_bottom:
        Normal component of the magnetic field at the lower boundary
        (``k = 0``) with shape ``(ix, jx)``.
    dxi, det:
        Uniform grid spacing in the ``xi`` and ``eta`` directions.
    lzt:
        Domain length in the ``zeta`` direction.
    kx:
        Number of grid points in the ``zeta`` direction. Grid points are
        located at ``zeta_k = k * dz`` where ``dz = lzt / (kx - 0.5)``.
    hxi, het, hzt:
        Scale factors ``hxi``, ``het``, ``hzt`` sampled at the grid points
        ``zeta_k = k * dz`` (``k = 0, ..., kx - 1``), with shape ``(kx,)``.

    Returns
    -------
    bxi, bet, bzt:
        Physical components of the magnetic field with shape ``(ix, jx, kx)``.

    Raises
    ------
    ValueError
        If input shapes or values are invalid.
    """
    bztb = np.asarray(bzt_bottom, dtype=np.float64)
    if bztb.ndim != 2:
        raise ValueError("bzt_bottom must be a 2D array.")
    ix, jx = bztb.shape

    hxi64 = np.asarray(hxi, dtype=np.float64)
    het64 = np.asarray(het, dtype=np.float64)
    hzt64 = np.asarray(hzt, dtype=np.float64)

    kx_i = int(kx)
    if (
        hxi64.shape != (kx_i,)
        or het64.shape != (kx_i,)
        or hzt64.shape != (kx_i,)
    ):
        raise ValueError("hxi, het, and hzt must have shape (kx,).")
    if np.any(hxi64 <= 0) or np.any(het64 <= 0) or np.any(hzt64 <= 0):
        raise ValueError("hxi, het, and hzt must be positive.")

    dxi_f = float(dxi)
    det_f = float(det)
    lzt_f = float(lzt)

    if dxi_f <= 0 or det_f <= 0:
        raise ValueError("dxi and det must be positive.")
    if lzt_f <= 0:
        raise ValueError("lzt must be positive.")
    if kx_i < 2:
        raise ValueError("kx must be at least 2.")

    dz = lzt_f / (kx_i - 0.5)

    # Scale factors at interior half-integer points: h_{k+1/2} = (h_k + h_{k+1}) / 2
    # shape (kx-1,)
    hxi_half = 0.5 * (hxi64[:-1] + hxi64[1:])
    het_half = 0.5 * (het64[:-1] + het64[1:])
    hzt_half = 0.5 * (hzt64[:-1] + hzt64[1:])

    # Ghost half-point at k = -1/2: 0th-order extrapolation from k = 0.
    hxi_neg = hxi64[0]
    het_neg = het64[0]
    hzt_neg = hzt64[0]

    # Ghost half-point at k = KX - 1/2: 0th-order extrapolation from k = KX-1.
    hxi_top = hxi64[-1]
    het_top = het64[-1]
    hzt_top = hzt64[-1]

    # Jacobian J = 1 / (hxi * het * hzt)
    J = 1.0 / (hxi64 * het64 * hzt64)  # (kx,)
    J_half = 1.0 / (hxi_half * het_half * hzt_half)  # (kx-1,)
    J_neg = 1.0 / (hxi_neg * het_neg * hzt_neg)  # scalar
    J_top = 1.0 / (hxi_top * het_top * hzt_top)  # scalar

    # alpha_{k+1/2} = 1 / (dz^2 * J_{k+1/2} * hzt_{k+1/2}^2)
    alpha = 1.0 / (dz**2 * J_half * hzt_half**2)  # (kx-1,)
    alpha_neg = 1.0 / (dz**2 * J_neg * hzt_neg**2)  # scalar
    alpha_top = 1.0 / (dz**2 * J_top * hzt_top**2)  # scalar

    # --- Build tridiagonal coefficient arrays (geometry only, shape (kx,)) ---
    # lower[k] multiplies f_{k-1}.  lower[0] is unused (no ghost in equation).
    lower_geo = np.zeros(kx_i)
    lower_geo[1:] = J[1:] * alpha  # J_k * alpha_{k-1/2}

    # upper[k] multiplies f_{k+1}.  upper[-1] is unused.
    upper_geo = np.zeros(kx_i)
    # k=0: ghost substitution merges lower ghost relation into the f[k+1] term.
    upper_geo[0] = J[0] * (alpha[0] + alpha_neg)
    upper_geo[1:-1] = J[1:-1] * alpha[1:]  # J_k * alpha_{k+1/2}
    # upper_geo[-1] = 0  (ghost f_{KX} = -f_{KX-1} absorbed into diagonal)

    # diag[k] = -(lower[k] + upper[k]) (before subtracting kappa^2)
    diag_geo = np.zeros(kx_i)
    diag_geo[0] = -upper_geo[0]
    diag_geo[1:-1] = -(lower_geo[1:-1] + upper_geo[1:-1])
    # k=KX-1: ghost gives extra -2*J[-1]*alpha_top on diagonal.
    diag_geo[-1] = -J[-1] * (2.0 * alpha_top + alpha[-1])

    # RHS: only the lower-boundary term is nonzero.
    rhs_geo = np.zeros(kx_i)
    rhs_geo[0] = -2.0 * J[0] * alpha_neg * hzt64[0] * dz

    # --- Wavenumbers ---
    k1 = 2.0 * np.pi * np.fft.fftfreq(ix) / dxi_f  # (ix,)
    k2 = 2.0 * np.pi * np.fft.fftfreq(jx) / det_f  # (jx,)

    # kappa^2 = k1^2/hxi^2 + k2^2/het^2, shape (ix, jx, kx)
    kappa2 = (k1[:, None, None] / hxi64[None, None, :]) ** 2 + (
        k2[None, :, None] / het64[None, None, :]
    ) ** 2

    # Full diagonal: shape (ix, jx, kx)
    diag = diag_geo[None, None, :] - kappa2

    # --- Solve tridiagonal system (Thomas algorithm, vectorised over m, n) ---
    # f[m, n, k] is the normalised potential for wavenumber pair (k1[m], k2[n]).
    f = _thomas_solve(
        lower=lower_geo,
        diag=diag,
        upper=upper_geo,
        rhs=np.broadcast_to(rhs_geo, (ix, jx, kx_i)).copy(),
    )  # shape (ix, jx, kx), real

    # --- Compute df/dzeta (needed for bzt) ---
    df_dz = np.empty_like(f)
    # k=0: exact from the lower BC  partial_f / partial_zeta = -hzt(0)
    df_dz[:, :, 0] = -hzt64[0]
    # Interior: centred finite difference
    df_dz[:, :, 1:-1] = (f[:, :, 2:] - f[:, :, :-2]) / (2.0 * dz)
    # k=KX-1: use upper ghost f_{KX} = -f_{KX-1}
    df_dz[:, :, -1] = (-f[:, :, -1] - f[:, :, -2]) / (2.0 * dz)

    # --- Magnetic field in wavenumber space ---
    A = np.fft.fft2(bztb)  # (ix, jx)
    Af = A[:, :, None] * f  # (ix, jx, kx)

    bxi_fft = -1j * k1[:, None, None] * Af
    bet_fft = -1j * k2[None, :, None] * Af
    bzt_fft = -A[:, :, None] * df_dz

    # --- Inverse FFT and normalise by scale factors ---
    bxi = np.real(np.fft.ifft2(bxi_fft, axes=(0, 1))) / hxi64[None, None, :]
    bet = np.real(np.fft.ifft2(bet_fft, axes=(0, 1))) / het64[None, None, :]
    bzt = np.real(np.fft.ifft2(bzt_fft, axes=(0, 1))) / hzt64[None, None, :]

    return bxi, bet, bzt


def _thomas_solve(
    lower: NDArray[np.float64],
    diag: NDArray[np.float64],
    upper: NDArray[np.float64],
    rhs: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Solve a tridiagonal system using the Thomas (TDMA) algorithm.

    Supports batched problems: the system axis is the last axis, and all
    leading axes are solved independently.

    Parameters
    ----------
    lower:
        Sub-diagonal with shape broadcastable to ``(..., n)``.
        ``lower[..., 0]`` is unused.
    diag:
        Main diagonal with shape ``(..., n)``.
    upper:
        Super-diagonal with shape broadcastable to ``(..., n)``.
        ``upper[..., -1]`` is unused.
    rhs:
        Right-hand side with shape ``(..., n)``.

    Returns
    -------
    NDArray[np.float64]
        Solution array with the same shape as ``rhs``.
    """
    n = rhs.shape[-1]
    diag = np.array(diag, dtype=np.float64)  # working copy
    rhs = rhs.copy()

    # Forward sweep
    for k in range(1, n):
        factor = lower[..., k] / diag[..., k - 1]
        diag[..., k] -= factor * upper[..., k - 1]
        rhs[..., k] -= factor * rhs[..., k - 1]

    # Back substitution
    x = np.zeros_like(rhs)
    x[..., -1] = rhs[..., -1] / diag[..., -1]
    for k in range(n - 2, -1, -1):
        x[..., k] = (rhs[..., k] - upper[..., k] * x[..., k + 1]) / diag[..., k]

    return x
