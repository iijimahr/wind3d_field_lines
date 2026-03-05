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
    off_screen: bool = False
    screenshot: str | None = None


def build_arcade_field(
    *,
    ba: float,
    la: float,
    decay_a: float,
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    z: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Build linear force-free arcade field components on a regular grid."""

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
    """Run tracing and visualize field lines with PyVista."""

    try:
        import pyvista as pv
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "pyvista is required for the demo. Install with: pip install -e '.[demo]'"
        ) from exc

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

    plotter = pv.Plotter(off_screen=config.off_screen)
    plotter.set_background("white")

    valid_line_count = 0
    for n in range(result.nx):
        lmin = int(max(1, result.lmin[n]))
        lmax = int(min(result.lx, result.lmax[n]))
        if lmax - lmin + 1 < 2:
            continue

        points = np.column_stack(
            (
                x_line[n, lmin - 1 : lmax],
                y_line[n, lmin - 1 : lmax],
                z_line[n, lmin - 1 : lmax],
            )
        )
        line = pv.lines_from_points(points, close=False)
        plotter.add_mesh(line, color="royalblue", line_width=4)
        valid_line_count += 1

    seeds = np.column_stack((seed_x, seed_y, seed_z))
    plotter.add_points(
        seeds,
        color="crimson",
        point_size=11,
        render_points_as_spheres=True,
    )

    plotter.add_axes()
    plotter.show_grid(xlabel="x [Mm]", ylabel="y [Mm]", zlabel="z [Mm]")
    plotter.add_title("Linear force-free arcade field-line demo")

    screenshot = config.screenshot
    if config.off_screen and screenshot is None:
        screenshot = "arcade_field_demo.png"

    print(
        "Demo summary: "
        f"grid=({config.nx}, {config.ny}, {config.nz}), "
        f"seeds={config.seed_count}, "
        f"valid_lines={valid_line_count}/{result.nx}, "
        f"l-range=[{int(result.lmin.min())}, {int(result.lmax.max())}]"
    )

    if screenshot is not None:
        plotter.show(screenshot=screenshot)
        print(f"Saved visualization to: {screenshot}")
    else:
        plotter.show()

    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Visualize traced arcade magnetic field lines with PyVista."
    )
    parser.add_argument(
        "--ba", type=float, default=6.0, help="Field strength scale Ba."
    )
    parser.add_argument(
        "--la", type=float, default=12.0, help="Arcade half-width scale La."
    )
    parser.add_argument(
        "--decay-a", type=float, default=30.0, help="Vertical decay length a."
    )
    parser.add_argument("--nx", type=int, default=65, help="Number of x-grid points.")
    parser.add_argument("--ny", type=int, default=65, help="Number of y-grid points.")
    parser.add_argument("--nz", type=int, default=65, help="Number of z-grid points.")
    parser.add_argument(
        "--seed-count", type=int, default=9, help="Number of seed points."
    )
    parser.add_argument("--seed-xmin", type=float, default=-8.0, help="Min seed x.")
    parser.add_argument("--seed-xmax", type=float, default=8.0, help="Max seed x.")
    parser.add_argument("--seed-y", type=float, default=0.0, help="Seed y coordinate.")
    parser.add_argument("--seed-z", type=float, default=0.0, help="Seed z coordinate.")
    parser.add_argument(
        "--lx-bln",
        type=int,
        default=101,
        help="Number of points along each field line.",
    )
    parser.add_argument(
        "--lcen-bln", type=int, default=51, help="Center index along each field line."
    )
    parser.add_argument("--margin", type=int, default=0, help="Ghost-cell margin.")
    parser.add_argument(
        "--nsubstepx", type=int, default=3, help="Integration substeps per segment."
    )
    parser.add_argument(
        "--off-screen",
        action="store_true",
        help="Enable off-screen rendering (useful for headless environments).",
    )
    parser.add_argument(
        "--screenshot",
        type=str,
        default=None,
        help=(
            "Screenshot output path. "
            "If omitted with --off-screen, defaults to arcade_field_demo.png."
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    config = ArcadeDemoConfig(
        ba=args.ba,
        la=args.la,
        decay_a=args.decay_a,
        nx=args.nx,
        ny=args.ny,
        nz=args.nz,
        seed_count=args.seed_count,
        seed_x_min=args.seed_xmin,
        seed_x_max=args.seed_xmax,
        seed_y=args.seed_y,
        seed_z=args.seed_z,
        lx_bln=args.lx_bln,
        lcen_bln=args.lcen_bln,
        margin=args.margin,
        nsubstepx=args.nsubstepx,
        off_screen=args.off_screen,
        screenshot=args.screenshot,
    )

    if config.seed_count <= 0:
        parser.error("--seed-count must be greater than 0.")
    if config.lx_bln <= 0:
        parser.error("--lx-bln must be greater than 0.")
    if not (1 <= config.lcen_bln <= config.lx_bln):
        parser.error("--lcen-bln must satisfy 1 <= lcen-bln <= lx-bln.")

    if config.screenshot is not None:
        out_path = Path(config.screenshot)
        if out_path.parent != Path("."):
            out_path.parent.mkdir(parents=True, exist_ok=True)

    return run_demo(config)


if __name__ == "__main__":
    raise SystemExit(main())
