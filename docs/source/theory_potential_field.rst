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

- **Lower boundary** (:math:`\zeta = 0`): the normal component of the
  surface magnetic field is prescribed, :math:`B_\zeta = B_0(\xi, \eta)`.
- **Upper boundary** (:math:`\zeta = L_\zeta`): the horizontal field vanishes,
  :math:`\mathbf{B}_h = 0`.

Coordinate system
-----------------

An orthogonal curvilinear coordinate system :math:`(\xi, \eta, \zeta)` is
assumed with scale factors :math:`h_j(\zeta)` that depend only on
:math:`\zeta`.  Both horizontal directions are periodic.  Special cases:

- **Cartesian**: :math:`h_j = 1`.
- **Local spherical** (Iijima et al. 2023):
  :math:`h_\xi = h_\eta = r\sqrt{f}`, :math:`h_\zeta = 1`.

The Jacobian is :math:`J = 1/(h_\xi h_\eta h_\zeta)`, and the Laplacian is

.. math::

   \nabla^2\Psi
   = J \sum_j \frac{\partial}{\partial \xi_j}
     \left(\frac{1}{J h_j^2} \frac{\partial \Psi}{\partial \xi_j}\right).

Solution method
---------------

Because the horizontal directions are periodic, we expand in horizontal
Fourier modes.  Writing

.. math::

   \Psi(\xi,\eta,\zeta)
   = \mathcal{F}^{-1}\!\left\{A(k_\xi,k_\eta)\,f(k_\xi,k_\eta,\zeta)\right\},

where :math:`A(k_\xi,k_\eta) = \mathcal{F}\{B_\zeta(\xi,\eta,0)\}` is the FFT of
the surface field, the Laplace equation reduces to an independent 1-D ODE for
each wavenumber pair :math:`(k_\xi,k_\eta)`:

.. math::

   J \frac{\partial}{\partial \zeta}
   \!\left(\frac{1}{J h_\zeta^2} \frac{\partial f}{\partial \zeta}\right)
   - \left(\frac{k_\xi^2}{h_\xi^2} + \frac{k_\eta^2}{h_\eta^2}\right) f = 0.

Boundary conditions for :math:`f`:

.. math::

   \frac{\partial f}{\partial \zeta} = -h_\zeta \quad (\zeta = 0), \qquad
   f = 0 \quad (\zeta = L_\zeta).

The magnetic field components are then recovered as

.. math::

   B_\xi = \mathcal{F}^{-1}\!\left\{-ik_\xi A f\right\} / h_\xi, \qquad
   B_\eta = \mathcal{F}^{-1}\!\left\{-ik_\eta A f\right\} / h_\eta, \qquad
   B_\zeta = \mathcal{F}^{-1}\!\left\{-A\,\partial f/\partial\zeta\right\} / h_\zeta.

Discretisation
--------------

Grid points are placed at :math:`\zeta_k = k\,\Delta\zeta`
(:math:`k = 0, \ldots, N_\zeta-1`), where :math:`\Delta\zeta` is the prescribed
uniform grid spacing. The upper boundary is therefore located at

.. math::

   L_\zeta = (N_\zeta - 1/2)\,\Delta\zeta.

The ODE is discretised with a finite-volume scheme at interior points:

.. math::

   \frac{J_k}{\Delta\zeta}
   \left(
     \frac{f_{k+1}-f_k}{J_{k+1/2}\,h_{\zeta,k+1/2}^2\,\Delta\zeta}
     - \frac{f_k-f_{k-1}}{J_{k-1/2}\,h_{\zeta,k-1/2}^2\,\Delta\zeta}
   \right)
   - \left(\frac{k_\xi^2}{h_{\xi,k}^2} + \frac{k_\eta^2}{h_{\eta,k}^2}\right) f_k = 0.

Quantities at half-integer points are linearly interpolated from the adjacent
grid points.  Ghost-cell values encode the boundary conditions:

- **Lower** (:math:`k = -1`):
  :math:`f_{-1} = f_1 + 2\,h_{\zeta,0}\,\Delta\zeta`
- **Upper** (:math:`k = N_\zeta`):
  :math:`f_{N_\zeta} = -f_{N_\zeta-1}`

The resulting tridiagonal system is solved with the Thomas (TDMA) algorithm,
vectorised simultaneously over all :math:`(k_\xi,k_\eta)` wavenumber pairs.

Analytical solution (Cartesian)
---------------------------------

For Cartesian coordinates (:math:`h_j = 1`) and a single horizontal mode
with :math:`\kappa = \sqrt{k_\xi^2 + k_\eta^2}`:

.. math::

   f(\zeta) =
   \begin{cases}
     L_\zeta - \zeta & (\kappa = 0), \\[4pt]
     \dfrac{\sinh\!\left(\kappa(L_\zeta-\zeta)\right)}{\kappa\cosh(\kappa L_\zeta)}
     & (\kappa \neq 0).
   \end{cases}

API
---

.. autofunction:: wind3d_field_lines.compute_potential_field
   :noindex:
