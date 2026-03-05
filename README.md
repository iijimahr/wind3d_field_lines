# wind3d_field_lines

Initial implementation of a magnetic field-line tracer for wind3d data.

## Public API (Current Stage)

- `trace_field_lines`: Trace magnetic field lines

## Development Setup

```bash
pip install -e .[dev]
```

## Test

```bash
pytest
```

## Documentation

Install docs dependencies:

```bash
pip install -e .[docs]
```

Build HTML docs:

```bash
make -C docs html
```
