wind3d_field_lines
==================

Magnetic field-line tracer for RAMENS wind3d data.

Two tracing functions and one field extrapolation function are provided:

- :func:`~wind3d_field_lines.trace_field_lines` — Cartesian grids with
  non-uniform spacing along the vertical axis.
- :func:`~wind3d_field_lines.trace_field_lines_curvilinear` — orthogonal
  curvilinear coordinate systems, using scale-factor rescaling of the field
  components (see :doc:`theory_curvilinear`).
- :func:`~wind3d_field_lines.compute_potential_field` — potential magnetic
  field extrapolation from a surface boundary condition, independent of the
  tracing functions (see :doc:`theory_potential_field`).

The tracing functions return a frozen dataclass containing the traced
coordinates and valid index ranges for each field line.

Status
======

This package is under active development for research and learning purposes.

Installation
============

.. code-block:: bash

   pip install git+https://github.com/iijimahr/wind3d_field_lines.git

Quick start
===========

Cartesian grid:

.. code-block:: python

   from wind3d_field_lines import trace_field_lines

   result = trace_field_lines(
       bx=bx, by=by, bz=bz,
       dx=dx_profile, dy=dy_profile, dz=dz_profile,
       seed_i=icen, seed_j=jcen, seed_k=kcen,
       line_center=151, line_length=301, margin=0,
   )
   # result.i/j/k : traced grid indices, shape (n_seeds, line_length)

Orthogonal curvilinear coordinates:

.. code-block:: python

   from wind3d_field_lines import trace_field_lines_curvilinear

   result = trace_field_lines_curvilinear(
       bxi=bxi, bet=bet, bzt=bzt,
       dxi=dxi, det=det, dzt=dzt,
       hxi=hxi, het=het, hzt=hzt,
       seed_i=icen, seed_j=jcen, seed_k=kcen,
       line_center=151, line_length=301, margin=0,
   )
   # result.xi/eta/zeta : traced physical coordinates, shape (n_seeds, line_length)

Potential field extrapolation:

.. code-block:: python

   import numpy as np
   from wind3d_field_lines import compute_potential_field

   ix, jx, kx = 64, 64, 32
   bxi, bet, bzt = compute_potential_field(
       bzt_bottom=bzt_bottom,       # surface normal field, shape (ix, jx)
       dxi=1.0, det=1.0,          # horizontal grid spacing
       lzt=20.0, kx=kx,            # vertical domain length and grid points
       hxi=np.ones(kx),           # scale factors hxi, het, hzt at each level
       het=np.ones(kx),
       hzt=np.ones(kx),
   )
   # bxi/bet/bzt : magnetic field components, shape (ix, jx, kx)

API Reference
=============

.. toctree::
   :maxdepth: 2

   api
   demo_arcade
   demo_bipolar
   theory_curvilinear
   theory_potential_field
