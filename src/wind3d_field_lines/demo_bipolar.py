"""Bipolar sunspot potential field demo with matplotlib 3D visualization."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .integrator import trace_field_lines
from .potential_field import compute_potential_field


@dataclass(frozen=True)
class BipolarDemoConfig:
    """Configuration for the bipolar sunspot potential field demo."""

    # Physical parameters
    b0: float = 100.0
    sigma: float = 8.0
    spot_distance: float = 20.0

    # Grid parameters
    nx: int = 128
    ny: int = 128
    nz: int = 64
    x_min: float = -60.0
    x_max: float = 60.0
    y_min: float = -60.0
    y_max: float = 60.0
    lzt: float = 60.0

    # Seed / tracing parameters
    seed_count: int = 30
    seed_rng_seed: int = 42
    line_length: int = 601
    line_center: int = 301
    margin: int = 0
    n_substeps: int = 3

    # Output
    output: str | None = None


def build_bipolar_surface_field(
    *,
    b0: float,
    sigma: float,
    spot_distance: float,
    x: NDArray[np.float64],
    y: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Build the bipolar sunspot surface magnetic field B_3(x, y, z=0).

    The field consists of two Gaussian spots of opposite polarity placed
    symmetrically at ``x = ±spot_distance``.

    Parameters
    ----------
    b0:
        Peak field strength [G].
    sigma:
        Gaussian width [Mm].
    spot_distance:
        Distance from the origin to each spot centre [Mm].
    x:
        1-D array of x coordinates [Mm].
    y:
        1-D array of y coordinates [Mm].

    Returns
    -------
    NDArray[np.float64]
        Surface field with shape ``(nx, ny)``.

    Raises
    ------
    ValueError
        If ``b0``, ``sigma``, or ``spot_distance`` are non-positive.
    """
    if b0 <= 0.0:
        raise ValueError("b0 must be positive.")
    if sigma <= 0.0:
        raise ValueError("sigma must be positive.")
    if spot_distance <= 0.0:
        raise ValueError("spot_distance must be positive.")

    xx, yy = np.meshgrid(x, y, indexing="ij")
    r2_pos = (xx - spot_distance) ** 2 + yy**2
    r2_neg = (xx + spot_distance) ** 2 + yy**2
    two_sig2 = 2.0 * sigma**2
    return b0 * (np.exp(-r2_pos / two_sig2) - np.exp(-r2_neg / two_sig2))


def _uniform_spacing(axis: NDArray[np.float64], name: str) -> float:
    if axis.size < 2:
        raise ValueError(f"{name} must contain at least 2 grid points.")
    spacing = float(axis[1] - axis[0])
    if spacing <= 0.0:
        raise ValueError(f"{name} must be strictly increasing.")
    return spacing


def _to_index(
    values: NDArray[np.float64], vmin: float, step: float
) -> NDArray[np.float64]:
    return (values - vmin) / step + 1.0


def _to_physical(
    values: NDArray[np.float64], vmin: float, step: float
) -> NDArray[np.float64]:
    return vmin + (values - 1.0) * step


def _domain_segment(
    xv: NDArray[np.float64],
    yv: NDArray[np.float64],
    center: int,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
) -> tuple[int, int]:
    """Return ``[a, b]``: the contiguous in-domain segment containing *center*.

    Parameters
    ----------
    xv, yv:
        Physical x and y coordinates along the traced line (1-D).
    center:
        0-based index of the seed point within *xv* / *yv*.
    x_min, x_max, y_min, y_max:
        Domain bounds.

    Returns
    -------
    tuple[int, int]
        Inclusive ``[a, b]`` range, or ``(-1, -1)`` if the centre itself is
        outside the domain.
    """
    in_domain = (xv >= x_min) & (xv <= x_max) & (yv >= y_min) & (yv <= y_max)
    if not in_domain[center]:
        return -1, -1

    out_idx = np.where(~in_domain)[0]
    before = out_idx[out_idx < center]
    after = out_idx[out_idx > center]

    a = int(before[-1]) + 1 if len(before) else 0
    b = int(after[0]) - 1 if len(after) else len(xv) - 1
    return a, b


def run_demo(config: BipolarDemoConfig) -> None:
    """Extrapolate the potential field and visualize traced field lines.

    Parameters
    ----------
    config:
        Demo configuration parameters.
    """
    import matplotlib
    import matplotlib.pyplot as plt

    if config.output is not None:
        matplotlib.use("Agg")

    # --- Build horizontal grids ---
    x = np.linspace(config.x_min, config.x_max, config.nx, dtype=np.float64)
    y = np.linspace(config.y_min, config.y_max, config.ny, dtype=np.float64)
    dxi = _uniform_spacing(x, "x")
    det = _uniform_spacing(y, "y")

    # --- Surface boundary field ---
    bzt_bottom = build_bipolar_surface_field(
        b0=config.b0,
        sigma=config.sigma,
        spot_distance=config.spot_distance,
        x=x,
        y=y,
    )

    # --- Potential field extrapolation ---
    # The vertical grid uses dz = lzt / (nz - 0.5) with z_k = k * dz.
    dz_step = config.lzt / (config.nz - 0.5)
    hones = np.ones(config.nz, dtype=np.float64)

    bxi, bet, bzt = compute_potential_field(
        bzt_bottom=bzt_bottom,
        dxi=dxi,
        det=det,
        lzt=config.lzt,
        kx=config.nz,
        hxi=hones,
        het=hones,
        hzt=hones,
    )

    # --- Uniform grid spacing profiles for the tracer ---
    dx_profile = np.full(config.nz, dxi, dtype=np.float64)
    dy_profile = np.full(config.nz, det, dtype=np.float64)
    dz_profile = np.full(config.nz, dz_step, dtype=np.float64)

    # --- Seed points scattered over both polarity spots at z = 0 ---
    # Split seed_count evenly between positive (x = +spot_distance) and
    # negative (x = -spot_distance) spots, each drawn from a 2-D Gaussian
    # with standard deviation equal to the spot width sigma.
    rng = np.random.default_rng(config.seed_rng_seed)
    n_pos = config.seed_count // 2
    n_neg = config.seed_count - n_pos
    seed_x = np.concatenate(
        [
            rng.normal(config.spot_distance, config.sigma, n_pos),
            rng.normal(-config.spot_distance, config.sigma, n_neg),
        ]
    )
    seed_y = rng.normal(0.0, config.sigma, config.seed_count)
    seed_z = np.zeros(config.seed_count, dtype=np.float64)

    icen = _to_index(seed_x, config.x_min, dxi)
    jcen = _to_index(seed_y, config.y_min, det)
    kcen = _to_index(seed_z, 0.0, dz_step)

    result = trace_field_lines(
        bx=bxi,
        by=bet,
        bz=bzt,
        dx=dx_profile,
        dy=dy_profile,
        dz=dz_profile,
        seed_i=icen,
        seed_j=jcen,
        seed_k=kcen,
        line_center=config.line_center,
        line_length=config.line_length,
        margin=config.margin,
        n_substeps=config.n_substeps,
    )

    x_line = _to_physical(result.i, config.x_min, dxi)
    y_line = _to_physical(result.j, config.y_min, det)
    z_line = _to_physical(result.k, 0.0, dz_step)

    # --- Visualization ---
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Field lines — clip each traced line to the horizontal domain so that
    # periodic wrap-around artefacts are not plotted.
    valid_line_count = 0
    for n in range(result.num_lines):
        lmin = int(max(1, result.lmin[n]))
        lmax = int(min(result.line_length, result.lmax[n]))
        if lmax - lmin + 1 < 2:
            continue

        sl = slice(lmin - 1, lmax)
        xv = x_line[n, sl]
        yv = y_line[n, sl]
        zv = z_line[n, sl]

        # Centre index within this slice
        center = max(0, min(len(xv) - 1, config.line_center - lmin))
        a, b = _domain_segment(
            xv,
            yv,
            center,
            config.x_min,
            config.x_max,
            config.y_min,
            config.y_max,
        )
        if b - a + 1 < 2:
            continue

        ax.plot(
            xv[a : b + 1],
            yv[a : b + 1],
            zv[a : b + 1],
            color="royalblue",
            linewidth=1.2,
        )
        valid_line_count += 1

    # Seed points
    ax.scatter(
        seed_x,
        seed_y,
        seed_z,
        color="crimson",
        s=40,
        zorder=5,
        label="Seed points",
    )

    ax.set_xlabel("x [Mm]")
    ax.set_ylabel("y [Mm]")
    ax.set_zlabel("z [Mm]")
    ax.set_title("Bipolar sunspot — potential field lines")
    ax.set_zlim(0.0, config.lzt)

    print(
        "Demo summary: "
        f"grid=({config.nx}, {config.ny}, {config.nz}), "
        f"seeds={config.seed_count}, "
        f"valid_lines={valid_line_count}/{result.num_lines}, "
        f"l-range=[{int(result.lmin.min())}, {int(result.lmax.max())}]"
    )

    if config.output is not None:
        out_path = Path(config.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved visualization to: {out_path}")
    else:
        plt.show()

    plt.close(fig)


if __name__ == "__main__":
    run_demo(BipolarDemoConfig(output="bipolar_demo.png"))
