---
title: Roadmap and open items
---

# Roadmap and open items

This page states the direction in one place.

## Goal

A robust, modern Fortran ASR code for problems CPA cannot handle: high-entropy
alloys, partially disordered compounds, short-range order, and off-diagonal /
environment disorder.

## Fixed decisions

| Item | Choice |
|---|---|
| Build | `fpm` after sourcing `env.sh` |
| Compiler | Intel `ifx` primary; gfortran is a bonus, not a gate |
| Linear algebra | MKL |
| xc | libxc ≥ 7 through the C ABI |
| Parallelism | OpenMP + MPI (`-DUSE_MPI`, [[asr_mpi]]) |
| Production method | Energy-dependent ASR (`route = "edc"` default, [[asr_edc]], and `route = "ed3"`, [[asr_ed3]]); the energy-independent route and KKR-ASA are reference cross-checks |

## Planned and open

* Non-collinear magnetism → needs block / complex-Hermitian recursion.
* Plug-in Hamiltonians: Wannier (`*_hr.dat`) and generic tight binding through
  [[asr_operator]].
* Multi-shell short-range order and correlated disorder beyond the first shell (`[[sro]]` and `[epi]` handle first-shell SRO and pair interactions).
* Explaining the ≈ 5 % excess Fe moment in the alloy route versus the ordered ASA.
* A seed that polarises only the valence region.
* Verifying the Luchini–Nex joiner against the original paper.
* Better handling of the bcc single-band terminator problem.
* Deeper recursion at fixed memory for magnetic d-band alloys.
