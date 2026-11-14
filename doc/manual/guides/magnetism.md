---
title: Magnetism
---

# Magnetism

Collinear spin polarisation is selected with `[scf] spin = 2`. This page collects
what is known about its numerics. Non-collinear magnetism is not available.

## Starting a polarised calculation

* `zeta0` is the initial polarisation of the atomic density.
* `spin_sign = -1` on a `[[species]]` starts that species antiparallel
  (ferrimagnets, antiferromagnets on `cell` supercells).
* For restarts from an existing potential:

```toml
[scf]
spin = 2
start_potential = "run1_pot.dat"   # written earlier by save_potential
```

## Default seeds

| `method` | `zeta0` default | Rationale |
|---|---|---|
| `kspace` | 0.4 | Ordered route is stable; use ≈ 0.4 for Fe |
| `asr` | **0.1** | Larger seeds run away |

!!! danger "Runaway with a large seed"
    In the alloy route the seed is applied to the whole atomic density
    *including the core*, splitting the core levels by about 13 Ry at the
    nucleus and starting at ≈ 3 μ<sub>B</sub>. For bcc Fe, `zeta0 = 0.4` **runs
    away** (moment → 8, \(E_F\) → +1.5 Ry) for every mixing strength and contour
    resolution tried, while `zeta0 = 0.1` converges in 35 iterations. A seed
    that polarises only the valence region is planned.

## What to expect

| System | Result | Note |
|---|---|---|
| bcc Fe, ordered route | 2.2 μ<sub>B</sub> | |
| bcc Fe, alloy route (A = B), depth 8/10/12/16 | 2.262 / 2.362 / 2.344 / 2.351 μ<sub>B</sub> | depth-converged ≈ 2.35; the ≈ 5–6 % excess over the ordered first-order ASA route is a *representation effect* (orthogonal second-order moments vs. generalised first-order), not a depth effect — **not yet explained** |
| Ni, ordered | 0.64 μ<sub>B</sub> | |
| B2 FeCo (k-space) | Fe 2.760, Co 1.768 μ<sub>B</sub> | vs. Elk LSDA 2.75 / ≈ 1.75 at 2.856 Å |
| B2 FeCo (real-space recursion, depth 24) | Fe 2.745, Co 1.716 μ<sub>B</sub> | |

## Limitations

* The moment of Ni-based random alloys needs more levels than exact generation
  affords without symmetry (see [Convergence](convergence.html)).
* Absolute moments differ by ≈ 4 % between ASA codes; LDA moments are not
  comparable with GGA references at tight tolerance (see
  [Validation](../validation.html)).
