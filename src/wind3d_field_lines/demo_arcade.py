"""Linear force-free arcade field demo with matplotlib 3D visualization."""

from __future__ import annotations

import argparse
import math
from collections.abc import Sequence
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
    y_min: float = -40.0
    y_max: float = 40.0
    z_min: float = 0.0
    z_max: float = 65.0
    seed_count: int = 9
    seed_x_min: float = -8.0
    seed_x_max: float = 8.0
    seed_y: float = 0.0
    seed_z: float = 0.0
    lx_bln: int = 101
    lcen_bln: int = 51
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


def run_demo(config: ArcadeDemoConfig) -> int:
    """Trace field lines and visualize with matplotlib 3D.

    Parameters
    ----------
    config:
        Demo configuration parameters.

    Returns
    -------
    int
        Exit code (0 on success).
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

    seed_x = np.linspace(config.seed_x_min, config.seed_x_max, config.seed_count)
    seed_y = np.full(config.seed_count, config.seed_y, dtype=np.float64)
    seed_z = np.full(config.seed_count, config.seed_z, dtype=np.float64)

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

    ax.set_xlabel("x [Mm]")
    ax.set_ylabel("y [Mm]")
    ax.set_zlabel("z [Mm]")
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
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a linear force-free arcade magnetic field and save it to a .npz file."
        )
    )
    parser.add_argument(
        "--ba", type=float, default=6.0, help="Field strength scale Ba [G]."
    )
    parser.add_argument(
        "--la", type=float, default=12.0, help="Arcade half-width scale La [Mm]."
    )
    parser.add_argument(
        "--decay-a", type=float, default=30.0, help="Vertical decay length a [Mm]."
    )
    parser.add_argument("--nx", type=int, default=65, help="Number of x-grid points.")
    parser.add_argument("--ny", type=int, default=65, help="Number of y-grid points.")
    parser.add_argument("--nz", type=int, default=65, help="Number of z-grid points.")
    parser.add_argument(
        "--x-min", type=float, default=-12.0, help="Grid x minimum [Mm]."
    )
    parser.add_argument(
        "--x-max", type=float, default=12.0, help="Grid x maximum [Mm]."
    )
    parser.add_argument(
        "--y-min", type=float, default=-40.0, help="Grid y minimum [Mm]."
    )
    parser.add_argument(
        "--y-max", type=float, default=40.0, help="Grid y maximum [Mm]."
    )
    parser.add_argument(
        "--z-min", type=float, default=0.0, help="Grid z minimum [Mm]."
    )
    parser.add_argument(
        "--z-max", type=float, default=65.0, help="Grid z maximum [Mm]."
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output .npz file path for the field data.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the wind3d-build-arcade-field command.

    Builds a linear force-free arcade magnetic field on a regular grid and
    saves the result to a NumPy ``.npz`` file containing arrays ``x``, ``y``,
    ``z``, ``bx``, ``by``, and ``bz``.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    x = np.linspace(args.x_min, args.x_max, args.nx, dtype=np.float64)
    y = np.linspace(args.y_min, args.y_max, args.ny, dtype=np.float64)
    z = np.linspace(args.z_min, args.z_max, args.nz, dtype=np.float64)

    bx, by, bz = build_arcade_field(
        ba=args.ba,
        la=args.la,
        decay_a=args.decay_a,
        x=x,
        y=y,
        z=z,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out_path, x=x, y=y, z=z, bx=bx, by=by, bz=bz)
    print(f"Saved arcade field to: {out_path}")
    print(f"  grid: ({args.nx}, {args.ny}, {args.nz})")
    print(f"  ba={args.ba} G, la={args.la} Mm, decay_a={args.decay_a} Mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
