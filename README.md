# wind3d_field_lines

Magnetic field-line tracer for RAMENS wind3d data.
Supports both Cartesian grids (`trace_field_lines`) and orthogonal curvilinear
coordinate systems (`trace_field_lines_curvilinear`).
Also provides potential magnetic field extrapolation (`compute_potential_field`).

Published docs: <https://iijimahr.github.io/wind3d_field_lines/>

## Quick start

Cartesian grid:

```python
from wind3d_field_lines import trace_field_lines

result = trace_field_lines(
    bx=bx, by=by, bz=bz,
    dx=dx_profile, dy=dy_profile, dz=dz_profile,
    seed_i=icen, seed_j=jcen, seed_k=kcen,
    line_center=151, line_length=301, margin=0,
)
# result.i/j/k : traced grid indices, shape (n_seeds, line_length)
```

Orthogonal curvilinear coordinates:

```python
from wind3d_field_lines import trace_field_lines_curvilinear

result = trace_field_lines_curvilinear(
    bxi=bxi, bet=bet, bzt=bzt,
    dxi=dxi, det=det, dzt=dzt,
    hxi=hxi, het=het, hzt=hzt,
    seed_i=icen, seed_j=jcen, seed_k=kcen,
    line_center=151, line_length=301, margin=0,
)
# result.xi/eta/zeta : traced physical coordinates, shape (n_seeds, line_length)
```

Breaking changes in 0.2.0:

- Tracing APIs are keyword-only.
- Legacy `*_bln` argument names were removed.
- Result metadata fields are now `line_center`, `num_lines`, and `line_length`.

Potential field extrapolation:

```python
from wind3d_field_lines import compute_potential_field
import numpy as np

ix, jx, kx = 64, 64, 32
lzt = 20.0  # domain height [Mm]
dzt = lzt / (kx - 0.5)

# Surface normal field (lower boundary)
bzt_bottom = ...  # shape (ix, jx)

bxi, bet, bzt = compute_potential_field(
    bzt_bottom=bzt_bottom,
    dxi=1.0, det=1.0, dzt=dzt,
    hxi=np.ones(kx), het=np.ones(kx), hzt=np.ones(kx),  # Cartesian
)
# bxi/bet/bzt : magnetic field components, shape (ix, jx, kx)
```

## Demo

See the [Arcade Field Demo](https://iijimahr.github.io/wind3d_field_lines/demo_arcade.html)
for a step-by-step example of tracing and visualizing magnetic field lines.

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
