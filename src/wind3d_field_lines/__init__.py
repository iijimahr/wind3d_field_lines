"""Public API for the wind3d_field_lines package."""

from .integrator import trace_field_lines, trace_field_lines_curvilinear
from .potential_field import compute_potential_field
from .types import CurvilinearFieldLineResult, FieldLineResult

__all__ = [
    "trace_field_lines",
    "trace_field_lines_curvilinear",
    "compute_potential_field",
    "FieldLineResult",
    "CurvilinearFieldLineResult",
]
