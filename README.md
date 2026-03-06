# wind3d_field_lines

Magnetic field-line tracer for wind3d data.
Supports both Cartesian grids (`trace_field_lines`) and orthogonal curvilinear
coordinate systems (`trace_field_lines_curvilinear`).

Published docs: <https://iijimahr.github.io/wind3d_field_lines/>

## Quick start

Cartesian grid:

```python
from wind3d_field_lines import trace_field_lines

result = trace_field_lines(
    bx=bx, by=by, bz=bz,
    dx=dx_profile, dy=dy_profile, dz=dz_profile,
    icen_bln=icen, jcen_bln=jcen, kcen_bln=kcen,
    lcen_bln=151, lx_bln=301, margin=0,
)
# result.i/j/k : traced grid indices, shape (n_seeds, lx_bln)
```

Orthogonal curvilinear coordinates:

```python
from wind3d_field_lines import trace_field_lines_curvilinear

result = trace_field_lines_curvilinear(
    bxi=bxi, bet=bet, bzt=bzt,
    dxi=dxi, det=det, dzt=dzt,
    hxi=hxi, het=het, hzt=hzt,
    icen_bln=icen, jcen_bln=jcen, kcen_bln=kcen,
    lcen_bln=151, lx_bln=301, margin=0,
)
# result.xi/eta/zeta : traced physical coordinates, shape (n_seeds, lx_bln)
```

## Demo

See the [Arcade Field Demo](https://iijimahr.github.io/wind3d_field_lines/demo_arcade.html)
for a step-by-step example of tracing and visualizing magnetic field lines.

```python
from wind3d_field_lines.demo_arcade import ArcadeDemoConfig, run_demo

run_demo(ArcadeDemoConfig())                          # interactive
run_demo(ArcadeDemoConfig(output="arcade_demo.png"))  # save to file
```

## For developers

### Installation

```shell
git clone https://github.com/iijimahr/wind3d_field_lines.git
cd wind3d_field_lines
python -m venv venv
. venv/bin/activate
pip install -U pip && pip install -e ".[dev,docs]"
```

### Task automation with Makefile

```shell
make test         # Run lint, doctest, and unit tests
make pytest       # Run unit tests
make typecheck    # Run static type checking
make docs         # Build documentation
make clean        # Clean build artifacts
```

### CI/CD

GitHub Actions runs validation on pushes and pull requests:

- `ruff check` (without auto-fix)
- `ruff format --check`
- `mypy src`
- `pytest`
- `sphinx-build -b doctest`
- `sphinx-build -b html`

Documentation is deployed to GitHub Pages only when changes are merged into
`main` (or manually via workflow dispatch), not on regular branch pushes.
