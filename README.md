# wind3d_field_lines

Initial implementation of a magnetic field-line tracer for wind3d data.

## Public API (Current Stage)

- `trace_field_lines`: Trace magnetic field lines
- `compute_open_field_fraction`: Compute open-field area filling factors
- `map_field_lines_to_height`: Map field-line footpoints to an observation height

## Development Setup

```bash
pip install -e .[dev]
```

## Test

```bash
pytest
```
