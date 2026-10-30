---
title: Roadmap and open items
---

# Roadmap and open items

The authoritative, current lists are `docs/ROADMAP.md`, `docs/TODO.md`,
`docs/STATUS.md` and `docs/MAGNETISM_TODO.md` in the repository. This page
states the direction in one place.

## Goal

A robust, modern Fortran ASR code for problems CPA cannot handle: high-entropy
alloys, disordered altermagnets, short-range order, and off-diagonal /
environment disorder.

## Fixed decisions

| Item | Choice |
|---|---|
| Build | `fpm` after sourcing `env.sh` |
| Compiler | Intel `ifx` primary; gfortran is a bonus, not a gate |
| Linear algebra | MKL |
| xc | libxc ≥ 7 through the C ABI |
| Parallelism | OpenMP; MPI / coarrays only if measurement shows a need |
| Production method | Energy-dependent ASR in the AMscf form (`asr_ed`); the energy-independent route and KKR-ASA are cross-checks |

## Planned and open

* Non-collinear magnetism → needs block / complex-Hermitian recursion.
* Plug-in Hamiltonians: Wannier (`*_hr.dat`) and generic tight binding through
  [[asr_operator]].
* Short-range order (Warren–Cowley) and correlated disorder.
* Explaining the ≈ 5 % excess Fe moment in the alloy route versus the ordered ASA.
* A seed that polarises only the valence region.
* Verifying the Luchini–Nex joiner against the original paper.
* Better handling of the bcc single-band terminator problem.
* Deeper recursion at fixed memory for magnetic d-band alloys.
