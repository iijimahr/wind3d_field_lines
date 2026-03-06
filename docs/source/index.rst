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

   pip install -e .

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
       b_xi=b_xi, b_eta=b_eta, b_zeta=b_zeta,
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

   n1, n2, n3 = 64, 64, 32
   b1, b2, b3 = compute_potential_field(
       b3_bottom=b3_bottom,       # surface normal field, shape (n1, n2)
       dxi=1.0, det=1.0,          # horizontal grid spacing
       l3=20.0, n3=n3,            # vertical domain length and grid points
       hxi=np.ones(n3),           # scale factors h1, h2, h3 at each level
       het=np.ones(n3),
       hzt=np.ones(n3),
   )
   # b1/b2/b3 : magnetic field components, shape (n1, n2, n3)

API Reference
=============

.. toctree::
   :maxdepth: 2

   api
   demo_arcade
   demo_bipolar
   theory_curvilinear
   theory_potential_field
