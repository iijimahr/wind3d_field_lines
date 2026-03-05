Arcade Field Demo (PyVista)
===========================

This demo visualizes traced magnetic field lines for a linear force-free arcade
field using ``pyvista``.

Physics model
=============

The demo uses the analytic field in ``arcade_field_demo.md``:

.. math::

   k = \frac{\pi}{2L_a},\quad
   C = \frac{2L_a}{\pi a},\quad
   S = \sqrt{1 - C^2}

.. math::

   B_x = -C B_a \cos(kx) e^{-z/a},\quad
   B_y = -S B_a \cos(kx) e^{-z/a},\quad
   B_z = B_a \sin(kx) e^{-z/a}

The default parameters are:

- ``B_a = 6``
- ``L_a = 12``
- ``a = 30``

with coordinates interpreted in Mm and field amplitude in G for demo readability.

How to run
==========

Install with demo dependencies:

.. code-block:: bash

   pip install -e ".[demo]"

Run the interactive demo:

.. code-block:: bash

   wind3d-arcade-demo

Run in off-screen mode and save a screenshot:

.. code-block:: bash

   wind3d-arcade-demo --off-screen --screenshot docs/build/arcade_demo.png

You can also run the module directly:

.. code-block:: bash

   python -m wind3d_field_lines.demo_arcade_pyvista --help

Notes
=====

- If ``2L_a / (\pi a) > 1``, the model is invalid and the demo exits with an error.
- This demo is intended as a usage example and exploration tool.
