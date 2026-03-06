Arcade Field Demo
=================

This page demonstrates how to use ``wind3d_field_lines`` to trace and visualize
magnetic field lines in a linear force-free arcade field.

Physics model
=============

The demo uses the analytic field described in ``arcade_field_demo.md``.  Given
arcade half-width :math:`L_a`, vertical decay length :math:`a`, and field
strength scale :math:`B_a`, the dimensionless parameters are:

.. math::

   k = \frac{\pi}{2L_a},\quad
   C = \frac{2L_a}{\pi a},\quad
   S = \sqrt{1 - C^2}

The three magnetic-field components on a right-handed Cartesian grid
:math:`(x, y, z)` — with :math:`z` pointing upward and the polarity inversion
line along :math:`x = 0` — are:

.. math::

   B_x = -C B_a \cos(kx)\, e^{-z/a},\quad
   B_y = -S B_a \cos(kx)\, e^{-z/a},\quad
   B_z =    B_a \sin(kx)\, e^{-z/a}

Default parameters (coordinates in Mm, field amplitude in G):

- :math:`B_a = 6\,\text{G}`
- :math:`L_a = 12\,\text{Mm}`
- :math:`a   = 30\,\text{Mm}`

The model is valid only when :math:`C = 2L_a/(\pi a) \leq 1`.

Usage example
=============

The following snippet reproduces the figure below.  No command-line call is
needed — you can embed this directly in a script or notebook:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from wind3d_field_lines.demo_arcade import (
       ArcadeDemoConfig,
       build_arcade_field,
       run_demo,
   )

   # Run with default parameters and display interactively
   config = ArcadeDemoConfig()
   run_demo(config)

To save the figure to a file instead:

.. code-block:: python

   config = ArcadeDemoConfig(output="arcade_demo.png")
   run_demo(config)

You can also call :func:`build_arcade_field` directly to obtain field
components for your own analysis:

.. code-block:: python

   import numpy as np
   from wind3d_field_lines.demo_arcade import build_arcade_field

   x = np.linspace(-12.0, 12.0, 65)
   y = np.linspace(-40.0, 40.0, 65)
   z = np.linspace(0.0, 65.0, 65)

   bx, by, bz = build_arcade_field(ba=6.0, la=12.0, decay_a=30.0, x=x, y=y, z=z)
   print(bx.shape)  # (65, 65, 65)

Visualization
=============

.. figure:: _static/arcade_demo.png
   :alt: Traced field lines of a linear force-free arcade field
   :width: 80%
   :align: center

   Magnetic field lines of the linear force-free arcade (blue) traced from
   nine seed points at the photospheric base (red dots).  Coordinates are
   in Mm and the field amplitude scale is 6 G.

Building field data from the command line
=========================================

The package also provides a command-line tool that saves the arcade field
to a NumPy ``.npz`` file for offline analysis:

.. code-block:: bash

   wind3d-build-arcade-field --output arcade_field.npz

   # Custom parameters
   wind3d-build-arcade-field --ba 8.0 --la 15.0 --decay-a 25.0 \
       --nx 33 --ny 33 --nz 33 \
       --output arcade_field_custom.npz

The output file contains arrays ``x``, ``y``, ``z``, ``bx``, ``by``, ``bz``
on the requested grid.  Run ``wind3d-build-arcade-field --help`` for the full
list of options.

Notes
=====

- If :math:`2L_a / (\pi a) > 1`, the model is invalid and both
  :func:`build_arcade_field` and the CLI exit with an error.
- The field is independent of :math:`y` by construction; the
  :math:`y` dimension is retained so that the full 3-D tracer can be
  used without modification.
