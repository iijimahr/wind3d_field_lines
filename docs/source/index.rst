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

Quick start
===========

Cartesian grid:

.. code-block:: python

   from wind3d_field_lines import trace_field_lines

   result = trace_field_lines(
       bx=bx, by=by, bz=bz,
       dx=dx_profile, dy=dy_profile, dz=dz_profile,
       icen_bln=icen, jcen_bln=jcen, kcen_bln=kcen,
       lcen_bln=151, lx_bln=301, margin=0,
   )
   # result.i/j/k : traced grid indices, shape (n_seeds, lx_bln)

Orthogonal curvilinear coordinates:

.. code-block:: python

   from wind3d_field_lines import trace_field_lines_curvilinear

   result = trace_field_lines_curvilinear(
       bxi=bxi, bet=bet, bzt=bzt,
       dxi=dxi, det=det, dzt=dzt,
       hxi=hxi, het=het, hzt=hzt,
       icen_bln=icen, jcen_bln=jcen, kcen_bln=kcen,
       lcen_bln=151, lx_bln=301, margin=0,
   )
   # result.xi/eta/zeta : traced physical coordinates, shape (n_seeds, lx_bln)

API Reference
=============

.. toctree::
   :maxdepth: 2

   api
   demo_arcade
   theory_curvilinear
