Potential Magnetic Field Extrapolation
=======================================

This page describes the algorithm used by
:func:`~wind3d_field_lines.compute_potential_field` to extrapolate a
potential magnetic field from a surface boundary condition in an orthogonal
curvilinear coordinate system.

Problem statement
-----------------

The magnetic field is derived from a scalar potential :math:`\Psi`:

.. math::

   \mathbf{B} = -\nabla\Psi.

The divergence-free condition :math:`\nabla\cdot\mathbf{B}=0` then requires
:math:`\Psi` to satisfy the Laplace equation:

.. math::

   \nabla^2\Psi = 0.

Boundary conditions:

- **Lower boundary** (:math:`\xi_3 = 0`): the normal component of the
  surface magnetic field is prescribed, :math:`B_3 = B_0(\xi_1, \xi_2)`.
- **Upper boundary** (:math:`\xi_3 = L_3`): the horizontal field vanishes,
  :math:`\mathbf{B}_h = 0`.

Coordinate system
-----------------

An orthogonal curvilinear coordinate system :math:`(\xi_1, \xi_2, \xi_3)` is
assumed with scale factors :math:`h_j(\xi_3)` that depend only on
:math:`\xi_3`.  Both horizontal directions are periodic.  Special cases:

- **Cartesian**: :math:`h_j = 1`.
- **Local spherical** (Iijima et al. 2023):
  :math:`h_1 = h_2 = r\sqrt{f}`, :math:`h_3 = 1`.

The Jacobian is :math:`J = 1/(h_1 h_2 h_3)`, and the Laplacian is

.. math::

   \nabla^2\Psi
   = J \sum_j \frac{\partial}{\partial \xi_j}
     \left(\frac{1}{J h_j^2} \frac{\partial \Psi}{\partial \xi_j}\right).

Solution method
---------------

Because the horizontal directions are periodic, we expand in horizontal
Fourier modes.  Writing

.. math::

   \Psi(\xi_1,\xi_2,\xi_3)
   = \mathcal{F}^{-1}\!\left\{A(k_1,k_2)\,f(k_1,k_2,\xi_3)\right\},

where :math:`A(k_1,k_2) = \mathcal{F}\{B_3(\xi_1,\xi_2,0)\}` is the FFT of
the surface field, the Laplace equation reduces to an independent 1-D ODE for
each wavenumber pair :math:`(k_1,k_2)`:

.. math::

   J \frac{\partial}{\partial \xi_3}
   \!\left(\frac{1}{J h_3^2} \frac{\partial f}{\partial \xi_3}\right)
   - \left(\frac{k_1^2}{h_1^2} + \frac{k_2^2}{h_2^2}\right) f = 0.

Boundary conditions for :math:`f`:

.. math::

   \frac{\partial f}{\partial \xi_3} = -h_3 \quad (\xi_3 = 0), \qquad
   f = 0 \quad (\xi_3 = L_3).

The magnetic field components are then recovered as

.. math::

   B_1 = \mathcal{F}^{-1}\!\left\{-ik_1 A f\right\} / h_1, \qquad
   B_2 = \mathcal{F}^{-1}\!\left\{-ik_2 A f\right\} / h_2, \qquad
   B_3 = \mathcal{F}^{-1}\!\left\{-A\,\partial f/\partial\xi_3\right\} / h_3.

Discretisation
--------------

Grid points are placed at :math:`\xi_{3,k} = k\,\Delta\xi_3`
(:math:`k = 0, \ldots, N_3-1`) with

.. math::

   \Delta\xi_3 = \frac{L_3}{N_3 - 1/2}.

The ODE is discretised with a finite-volume scheme at interior points:

.. math::

   \frac{J_k}{\Delta\xi_3}
   \left(
     \frac{f_{k+1}-f_k}{J_{k+1/2}\,h_{3,k+1/2}^2\,\Delta\xi_3}
     - \frac{f_k-f_{k-1}}{J_{k-1/2}\,h_{3,k-1/2}^2\,\Delta\xi_3}
   \right)
   - \left(\frac{k_1^2}{h_{1,k}^2} + \frac{k_2^2}{h_{2,k}^2}\right) f_k = 0.

Quantities at half-integer points are linearly interpolated from the adjacent
grid points.  Ghost-cell values encode the boundary conditions:

- **Lower** (:math:`k = -1`):
  :math:`f_{-1} = f_1 + 2\,h_{3,0}\,\Delta\xi_3`
- **Upper** (:math:`k = N_3`):
  :math:`f_{N_3} = -f_{N_3-1}`

The resulting tridiagonal system is solved with the Thomas (TDMA) algorithm,
vectorised simultaneously over all :math:`(k_1,k_2)` wavenumber pairs.

Analytical solution (Cartesian)
---------------------------------

For Cartesian coordinates (:math:`h_j = 1`) and a single horizontal mode
with :math:`\kappa = \sqrt{k_1^2 + k_2^2}`:

.. math::

   f(\xi_3) =
   \begin{cases}
     L_3 - \xi_3 & (\kappa = 0), \\[4pt]
     \dfrac{\sinh\!\left(\kappa(L_3-\xi_3)\right)}{\kappa\cosh(\kappa L_3)}
     & (\kappa \neq 0).
   \end{cases}

API
---

.. autofunction:: wind3d_field_lines.compute_potential_field
   :noindex:
