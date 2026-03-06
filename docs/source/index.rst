wind3d_field_lines
==================

Magnetic field-line tracer for wind3d data, backed by a Fortran extension.

Two tracing functions are provided:

- :func:`~wind3d_field_lines.trace_field_lines` — Cartesian grids with
  non-uniform spacing along the vertical axis.
- :func:`~wind3d_field_lines.trace_field_lines_curvilinear` — orthogonal
  curvilinear coordinate systems, using scale-factor rescaling of the field
  components (see :doc:`theory_curvilinear`).

Both functions return a frozen dataclass containing the traced coordinates and
valid index ranges for each field line.

Status
======

This package is under active development for research and learning purposes.

Installation
============

.. code-block:: bash

   pip install -e .

API Reference
=============

.. toctree::
   :maxdepth: 2

   api
   demo_arcade
   theory_curvilinear
