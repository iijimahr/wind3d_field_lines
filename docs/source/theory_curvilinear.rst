Magnetic Field Line Tracing in Orthogonal Curvilinear Coordinates
=================================================================

This page explains why field-line tracing in an orthogonal curvilinear
coordinate system is mathematically equivalent to the standard Cartesian
streamline problem after a simple rescaling of the field components.

Coordinate system
-----------------

Consider an orthogonal curvilinear coordinate system

.. math::

   (\xi, \eta, \zeta)

with scale factors

.. math::

   (h_\xi, h_\eta, h_\zeta).

The magnetic field is expressed through its physical components

.. math::

   \mathbf{B} =
   B_\xi \hat{\mathbf{e}}_\xi
   + B_\eta \hat{\mathbf{e}}_\eta
   + B_\zeta \hat{\mathbf{e}}_\zeta.

Definition of magnetic field lines
------------------------------------

A magnetic field line is a curve whose tangent vector is parallel to
the magnetic field:

.. math::

   \frac{d\mathbf{x}}{ds} \parallel \mathbf{B}.

Using the coordinate basis, the infinitesimal displacement is

.. math::

   d\mathbf{x}
   = h_\xi\,d\xi\,\hat{\mathbf{e}}_\xi
   + h_\eta\,d\eta\,\hat{\mathbf{e}}_\eta
   + h_\zeta\,d\zeta\,\hat{\mathbf{e}}_\zeta.

The field-line condition therefore gives

.. math::

   \frac{h_\xi\,d\xi}{B_\xi}
   = \frac{h_\eta\,d\eta}{B_\eta}
   = \frac{h_\zeta\,d\zeta}{B_\zeta}.

Reduction to a Cartesian-type form
------------------------------------

Define the Jacobian of the coordinate transformation

.. math::

   J = h_\xi h_\eta h_\zeta.

Introduce the scaled magnetic field components

.. math::

   \tilde{B}_\xi   = h_\eta h_\zeta B_\xi, \qquad
   \tilde{B}_\eta  = h_\zeta h_\xi  B_\eta, \qquad
   \tilde{B}_\zeta = h_\xi  h_\eta  B_\zeta.

With this definition, the divergence-free condition
:math:`\nabla\cdot\mathbf{B}=0` becomes

.. math::

   \frac{\partial \tilde{B}_\xi}{\partial \xi}
   + \frac{\partial \tilde{B}_\eta}{\partial \eta}
   + \frac{\partial \tilde{B}_\zeta}{\partial \zeta}
   = 0,

which is identical in form to the divergence condition in Cartesian
coordinates.

Field line equations
---------------------

Magnetic field lines can therefore be described by

.. math::

   \frac{d\xi}{d\lambda}   = \tilde{B}_\xi, \qquad
   \frac{d\eta}{d\lambda}  = \tilde{B}_\eta, \qquad
   \frac{d\zeta}{d\lambda} = \tilde{B}_\zeta.

The vector :math:`\tilde{\mathbf{B}} = (\tilde{B}_\xi, \tilde{B}_\eta, \tilde{B}_\zeta)`
is divergence-free in the computational coordinates :math:`(\xi, \eta, \zeta)`,
so magnetic field lines correspond to streamlines of this vector field exactly
as in Cartesian coordinates.

Summary
-------

For a magnetic field given by physical components
:math:`(B_\xi, B_\eta, B_\zeta)` in an orthogonal curvilinear coordinate
system:

1. Compute the scaled components

   .. math::

      \tilde{B}_\xi   = h_\eta h_\zeta B_\xi, \qquad
      \tilde{B}_\eta  = h_\zeta h_\xi  B_\eta, \qquad
      \tilde{B}_\zeta = h_\xi  h_\eta  B_\zeta.

2. The divergence constraint becomes the standard Cartesian form

   .. math::

      \partial_\xi \tilde{B}_\xi
      + \partial_\eta \tilde{B}_\eta
      + \partial_\zeta \tilde{B}_\zeta = 0.

3. Magnetic field lines are obtained from

   .. math::

      \frac{d\boldsymbol{\xi}}{d\lambda} = \tilde{\mathbf{B}}.

Thus the field-line tracing problem in orthogonal curvilinear coordinates
reduces analytically to the standard Cartesian streamline problem in the
coordinate space :math:`(\xi, \eta, \zeta)`.

API
---

.. autofunction:: wind3d_field_lines.trace_field_lines_curvilinear
   :noindex:

.. autoclass:: wind3d_field_lines.CurvilinearFieldLineResult
   :members:
   :noindex:
