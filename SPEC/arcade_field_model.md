# Linear force-free arcade field (demo data for field-line tracing)

This file provides a simple analytic magnetic field that is convenient for testing a **magnetic field-line tracing** module.

It is adapted from the initial arcade field used in Kaneko & Yokoyama (2017, *ApJ* 845:12), but **rewritten so that `z` is the height (vertical) direction** while keeping a right-handed `(x, y, z)` coordinate system.

---

## Coordinate system

- Right-handed coordinates: **`x` (horizontal), `y` (horizontal), `z` (height)**.
- The horizontal plane is the **`x–y`** plane.
- The polarity inversion line (PIL) is located at **`x = 0`**.
- The field decays exponentially with height `z`.

---

## Magnetic field

Let

\[
k \equiv \frac{\pi}{2L_a},
\qquad
C \equiv \frac{2L_a}{\pi a},
\qquad
S \equiv \sqrt{1 - C^2},
\]

where \(B_a\) is the field-strength scale, \(L_a\) controls the horizontal periodicity/arcade width, and \(a\) is the magnetic decay length (also a typical coronal scale height in the original setup).

Then the field is

\[
B_x(x,y,z)= - C\, B_a \cos(kx)\, e^{-z/a},
\]

\[
B_y(x,y,z)= - S\, B_a \cos(kx)\, e^{-z/a},
\]

\[
B_z(x,y,z)= \;\;\;\; B_a \sin(kx)\, e^{-z/a}.
\]

Notes:

- This is a **linear force-free** field (constant-\(\alpha\)) and is divergence-free.
- Condition for real-valued \(S\): \(C \le 1\) i.e. \( \dfrac{2L_a}{\pi a} \le 1\).

---

## Typical parameter values (recommended demo defaults)

These are representative values used in the paper and work well for demos:

| Symbol | Meaning | Typical value | Notes |
|---:|---|---:|---|
| \(B_a\) | field strength scale | **6 G** | 1 G = 1e-4 T |
| \(L_a\) | arcade half-width scale | **12 Mm** | 1 Mm = 1e6 m |
| \(a\) | decay length | **30 Mm** | must satisfy \(2L_a/(\pi a) < 1\) |

With these values:

- \(C = 2L_a/(\pi a) \approx 0.2546\)
- \(S = \sqrt{1-C^2} \approx 0.9670\)
- \(k = \pi/(2L_a) \approx 0.1309\ \text{Mm}^{-1}\)

---

## Suggested demo domain & seed points

A compact domain that shows clear arcade structure:

- \(x \in [-12,\, 12]\) Mm  
- \(y \in [-40,\, 40]\) Mm  *(field is independent of `y`, but include it to test 3D code paths)*  
- \(z \in [0,\, 65]\) Mm

Typical seed points for tracing:

- Footpoint line at the base: `z = 0`, `y = 0`, with `x` sampled in `[-10, 10]` Mm.
- Example: `(x, y, z) = (-8, 0, 0), (-6, 0, 0), ..., (8, 0, 0)`.

---

## Units & scaling tips

- If your tracer expects SI units, convert:
  - \(B_a = 6\ \text{G} = 6\times 10^{-4}\ \text{T}\)
  - \(L_a = 12\ \text{Mm} = 1.2\times 10^{7}\ \text{m}\)
  - \(a = 30\ \text{Mm} = 3.0\times 10^{7}\ \text{m}\)

- You can freely rescale:
  - **amplitude** via \(B_a\),
  - **horizontal scale** via \(L_a\),
  - **vertical decay** via \(a\).

---

## Quick sanity checks (recommended)

1. **Right-handedness / vertical direction**
   - Confirm that increasing `z` decreases \(|B|\) via \(e^{-z/a}\).

2. **Divergence-free**
   - Numerically verify \(\nabla\cdot B \approx 0\) on a grid.

3. **Symmetry**
   - At `x=0`, \(B_z=0\) and \((B_x,B_y)\propto -\cos(0)=-1\).
   - At `x=±L_a`, \(kx=±\pi/2\), so \(\cos=0\) and the field is purely vertical: \(B_z=\pm B_a e^{-z/a}\).

---

## Minimal reference implementation (optional)

Below is a small function you can paste into a test harness (units consistent with your choice).

```python
import math

def arcade_B(x, y, z, Ba=6.0, La=12.0, a=30.0):
    # x,y,z in Mm, Ba in G by default
    k = math.pi / (2.0 * La)
    C = 2.0 * La / (math.pi * a)
    S = math.sqrt(max(0.0, 1.0 - C*C))
    ez = math.exp(-z / a)
    cx = math.cos(k * x)
    sx = math.sin(k * x)
    Bx = -C * Ba * cx * ez
    By = -S * Ba * cx * ez
    Bz =  Ba * sx * ez
    return (Bx, By, Bz)
```

(If your code uses SI, convert inputs/outputs accordingly.)
