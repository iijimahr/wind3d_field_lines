"""Linear force-free arcade field demo with matplotlib 3D visualization."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .integrator import trace_field_lines


@dataclass(frozen=True)
class ArcadeDemoConfig:
    """Configuration for the linear force-free arcade demo."""

    ba: float = 6.0
    la: float = 12.0
    decay_a: float = 30.0
    nx: int = 65
    ny: int = 65
    nz: int = 65
    x_min: float = -12.0
    x_max: float = 12.0
    y_min: float = -50.0
    y_max: float = 50.0
    z_min: float = 0.0
    z_max: float = 65.0
    seed_count: int = 30
    seed_rng_seed: int = 42
    lx_bln: int = 301
    lcen_bln: int = 151
    margin: int = 0
    nsubstepx: int = 3
    output: str | None = None


def build_arcade_field(
    *,
    ba: float,
    la: float,
    decay_a: float,
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    z: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Build linear force-free arcade field components on a regular grid.

    Parameters
    ----------
    ba:
        Field strength scale [G].
    la:
        Arcade half-width scale [Mm].
    decay_a:
        Vertical decay length [Mm].
    x:
        1-D array of x coordinates [Mm].
    y:
        1-D array of y coordinates [Mm].
    z:
        1-D array of z coordinates [Mm].

    Returns
    -------
    tuple of (bx, by, bz)
        Magnetic field components on the 3-D grid, each with shape
        ``(nx, ny, nz)``.

    Raises
    ------
    ValueError
        If ``la`` or ``decay_a`` are non-positive, or if
        ``2*la / (pi*decay_a) > 1``.
    """
    if la <= 0.0:
        raise ValueError("la must be positive.")
    if decay_a <= 0.0:
        raise ValueError("decay_a must be positive.")

    c = 2.0 * la / (math.pi * decay_a)
    if c > 1.0:
        raise ValueError("Invalid parameters: 2*la/(pi*decay_a) must be <= 1.")

    k = math.pi / (2.0 * la)
    s = math.sqrt(max(0.0, 1.0 - c * c))

    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    del yy  # The analytic field is independent of y.

    ez = np.exp(-zz / decay_a)
    cos_kx = np.cos(k * xx)
    sin_kx = np.sin(k * xx)

    bx = -c * ba * cos_kx * ez
    by = -s * ba * cos_kx * ez
    bz = ba * sin_kx * ez
    return bx, by, bz


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


def run_demo(config: ArcadeDemoConfig) -> None:
    """Trace field lines and visualize with matplotlib 3D.

    Parameters
    ----------
    config:
        Demo configuration parameters.
    """
    import matplotlib
    import matplotlib.pyplot as plt

    if config.output is not None:
        matplotlib.use("Agg")

    x = np.linspace(config.x_min, config.x_max, config.nx, dtype=np.float64)
    y = np.linspace(config.y_min, config.y_max, config.ny, dtype=np.float64)
    z = np.linspace(config.z_min, config.z_max, config.nz, dtype=np.float64)

    dx_step = _uniform_spacing(x, "x")
    dy_step = _uniform_spacing(y, "y")
    dz_step = _uniform_spacing(z, "z")

    bx, by, bz = build_arcade_field(
        ba=config.ba,
        la=config.la,
        decay_a=config.decay_a,
        x=x,
        y=y,
        z=z,
    )

    dx_profile = np.full(config.nz, dx_step, dtype=np.float64)
    dy_profile = np.full(config.nz, dy_step, dtype=np.float64)
    dz_profile = np.full(config.nz, dz_step, dtype=np.float64)

    # Scatter seed points randomly over the photospheric surface (z = 0).
    # Avoid |x| < 0.5 Mm near the polarity inversion line where Bz ~ 0.
    rng = np.random.default_rng(config.seed_rng_seed)
    x_half = config.x_max * 0.85
    y_half = config.y_max * 0.96
    seed_x_raw = rng.uniform(-x_half, x_half, config.seed_count)
    # Nudge seeds very close to the PIL away from it so tracing is meaningful.
    pil_margin = 0.5
    seed_x = np.where(
        np.abs(seed_x_raw) < pil_margin,
        np.sign(seed_x_raw + 1e-9) * pil_margin,
        seed_x_raw,
    )
    seed_y = rng.uniform(-y_half, y_half, config.seed_count)
    seed_z = np.zeros(config.seed_count, dtype=np.float64)

    icen = _to_index(seed_x, config.x_min, dx_step)
    jcen = _to_index(seed_y, config.y_min, dy_step)
    kcen = _to_index(seed_z, config.z_min, dz_step)

    result = trace_field_lines(
        bx=bx,
        by=by,
        bz=bz,
        dx=dx_profile,
        dy=dy_profile,
        dz=dz_profile,
        icen_bln=icen,
        jcen_bln=jcen,
        kcen_bln=kcen,
        lcen_bln=config.lcen_bln,
        lx_bln=config.lx_bln,
        margin=config.margin,
        nsubstepx=config.nsubstepx,
    )

    x_line = _to_physical(result.i, config.x_min, dx_step)
    y_line = _to_physical(result.j, config.y_min, dy_step)
    z_line = _to_physical(result.k, config.z_min, dz_step)

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    valid_line_count = 0
    for n in range(result.nx):
        lmin = int(max(1, result.lmin[n]))
        lmax = int(min(result.lx, result.lmax[n]))
        if lmax - lmin + 1 < 2:
            continue

        ax.plot(
            x_line[n, lmin - 1 : lmax],
            y_line[n, lmin - 1 : lmax],
            z_line[n, lmin - 1 : lmax],
            color="royalblue",
            linewidth=1.2,
        )
        valid_line_count += 1

    ax.scatter(
        seed_x,
        seed_y,
        seed_z,
        color="crimson",
        s=40,
        zorder=5,
        label="Seed points",
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("Linear force-free arcade field lines")

    print(
        "Demo summary: "
        f"grid=({config.nx}, {config.ny}, {config.nz}), "
        f"seeds={config.seed_count}, "
        f"valid_lines={valid_line_count}/{result.nx}, "
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
