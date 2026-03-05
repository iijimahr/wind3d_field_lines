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

### Task automation with Makefile

```shell
make test         # Run lint, doctest, and unit tests
make pytest       # Run unit tests
make docs         # Build documentation
make clean        # Clean build artifacts
```
