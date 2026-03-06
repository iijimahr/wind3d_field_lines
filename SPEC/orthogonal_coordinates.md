# Magnetic Field Line Tracing in Orthogonal Curvilinear Coordinates

Consider an orthogonal curvilinear coordinate system

\[
(\xi,\eta,\zeta)
\]

with scale factors

\[
(h_\xi, h_\eta, h_\zeta).
\]

The magnetic field is expressed using its physical components

\[
\mathbf{B} =
B_\xi \hat{\mathbf e}_\xi
+
B_\eta \hat{\mathbf e}_\eta
+
B_\zeta \hat{\mathbf e}_\zeta .
\]

---

## Definition of magnetic field lines

A magnetic field line is a curve whose tangent vector is parallel to the magnetic field:

\[
\frac{d\mathbf{x}}{ds} \parallel \mathbf{B}.
\]

Using the coordinate basis, the infinitesimal displacement is

\[
d\mathbf{x} = h_\xi d\xi \hat{\mathbf e}_\xi + h_\eta d\eta \hat{\mathbf e}_\eta + h_\zeta d\zeta \hat{\mathbf e}_\zeta.
\]

The field line condition therefore gives

\[
\frac{h_\xi d\xi}{B_\xi} = \frac{h_\eta d\eta}{B_\eta} = \frac{h_\zeta d\zeta}{B_\zeta}.
\]

This relation determines the geometry of the magnetic field lines.

---

## Reduction to a Cartesian-type form

Define the Jacobian of the coordinate transformation

\[
J = h_\xi h_\eta h_\zeta .
\]

Introduce the scaled magnetic field components

\[
\tilde B_\xi = h_\eta h_\zeta B_\xi,
\]

\[
\tilde B_\eta = h_\zeta h_\xi B_\eta,
\]

\[
\tilde B_\zeta = h_\xi h_\eta B_\zeta .
\]

With this definition, the divergence-free condition

\[
\nabla\cdot\mathbf{B}=0
\]

becomes

\[
\frac{\partial \tilde B_\xi}{\partial \xi}
+
\frac{\partial \tilde B_\eta}{\partial \eta}
+
\frac{\partial \tilde B_\zeta}{\partial \zeta}
= 0 .
\]

This is identical in form to the divergence condition in Cartesian coordinates.

---

## Field line equations

Magnetic field lines can therefore be described by

\[
\frac{d\xi}{d\lambda}=\tilde B_\xi,
\]

\[
\frac{d\eta}{d\lambda}=\tilde B_\eta,
\]

\[
\frac{d\zeta}{d\lambda}=\tilde B_\zeta .
\]

The vector

\[
\tilde{\mathbf B}
= (\tilde B_\xi,\tilde B_\eta,\tilde B_\zeta)
\]

is divergence-free in the computational coordinates

\[
(\xi,\eta,\zeta),
\]

so magnetic field lines correspond to streamlines of this vector field exactly as in Cartesian coordinates.

---

## Summary

For a magnetic field given by physical components \((B_\xi,B_\eta,B_\zeta)\) in an orthogonal curvilinear coordinate system:

1. Define the scaled components

\[
\tilde B_\xi = h_\eta h_\zeta B_\xi,\quad
\tilde B_\eta = h_\zeta h_\xi B_\eta,\quad
\tilde B_\zeta = h_\xi h_\eta B_\zeta.
\]

1. The divergence constraint becomes

\[
\partial_\xi \tilde B_\xi +
\partial_\eta \tilde B_\eta +
\partial_\zeta \tilde B_\zeta = 0 .
\]

1. Magnetic field lines are obtained from

\[
\frac{d\boldsymbol{\xi}}{d\lambda}
= \tilde{\mathbf B}.
\]

Thus the field-line tracing problem in orthogonal curvilinear coordinates reduces analytically to the standard Cartesian streamline problem in the coordinate space \((\xi,\eta,\zeta)\).
