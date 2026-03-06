# wind3d_field_lines

Initial implementation of a magnetic field-line tracer for wind3d data.

## Documentation

Local HTML docs can be built with:

```shell
make docs
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

To use the demo visualization (requires matplotlib):

```shell
pip install -e ".[demo]"
```

### Task automation with Makefile

```shell
make test         # Run lint, doctest, and unit tests
make pytest       # Run unit tests
make typecheck    # Run static type checking
make docs         # Build documentation
make clean        # Clean build artifacts
```

## Demo

See the [Arcade Field Demo](https://iijimahr.github.io/wind3d_field_lines/demo_arcade.html)
page in the documentation for a step-by-step example of tracing and
visualizing magnetic field lines.

Quick start:

```python
from wind3d_field_lines.demo_arcade import ArcadeDemoConfig, run_demo

run_demo(ArcadeDemoConfig())                          # interactive
run_demo(ArcadeDemoConfig(output="arcade_demo.png"))  # save to file
```

## CI/CD

GitHub Actions runs validation on pushes and pull requests:

- `ruff check` (without auto-fix)
- `ruff format --check`
- `mypy src`
- `pytest`
- `sphinx-build -b doctest`
- `sphinx-build -b html`

Documentation is deployed to GitHub Pages only when changes are merged into
`main` (or manually via workflow dispatch), not on regular branch pushes.

Published docs URL:
<https://iijimahr.github.io/wind3d_field_lines/>
