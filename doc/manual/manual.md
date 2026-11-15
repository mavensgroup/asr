---
title: Manual
author: Rudra Banerjee
ordered_subpage: installation.md
                 quickstart.md
                 input
                 output.md
                 hamiltonian
                 theory
                 guides
                 architecture.md
                 validation.md
                 developing.md
                 appendix
---

# asr User and Developer Manual

[![GitHub](https://img.shields.io/badge/GitHub-mavensgroup%2Fasr-181717?logo=github&logoColor=white)](https://github.com/mavensgroup/asr)
[![Release](https://img.shields.io/badge/release-v0.1.0-blue.svg)](https://github.com/mavensgroup/asr/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Fortran](https://img.shields.io/badge/Fortran-2008%20%7C%202018-734f96?logo=fortran&logoColor=white)](https://fortran-lang.org)
[![Build with fpm](https://img.shields.io/badge/build-fpm-008080?logo=fortran&logoColor=white)](https://fpm.fortran-lang.org)
[![Compiler](https://img.shields.io/badge/compiler-Intel%20oneAPI-0071C5?logo=intel&logoColor=white)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
[![libxc](https://img.shields.io/badge/libxc-≥%207.0-2e7d32)](https://libxc.gitlab.io)
[![OpenMP](https://img.shields.io/badge/parallel-OpenMP-f15a24)](https://www.openmp.org)

**asr** is a from-scratch modern-Fortran code for the electronic structure of
*random alloys beyond the coherent-potential approximation (CPA)*: high-entropy
alloys, partially disordered Heusler compounds, disordered magnets, and
short-range-order problems. Disorder is averaged **exactly** in an *augmented
space* and the averaged Green function is obtained by the Haydock (Lanczos)
recursion with a terminated continued fraction.

The method is **Hamiltonian-agnostic**: it needs an operator on a graph with
species-dependent site parameters and species-independent hopping, and nothing
more ([Hamiltonians](hamiltonian/index.html)). At present one Hamiltonian is
implemented, a self-consistent **TB-LMTO-ASA** that the code builds itself:
radial scalar-relativistic solver, libxc exchange–correlation, screened
structure constants, Ewald/Madelung sums, Anderson mixing. No potential
parameters are tabulated or imported.

## Key Capabilities

- **Beyond the CPA**: Averages the Green function exactly over disorder in the Mookerjee augmented space via Haydock Lanczos recursion, preserving local environment fluctuations and cluster effects without single-site approximations.
- **Hamiltonian-agnostic method**: The disorder average and the recursion see only an operator ([[asr_operator]]) and class parameters ([[asr_augcls]]); the Hamiltonian is a separate, replaceable layer.
- **From-scratch TB-LMTO-ASA (the implemented Hamiltonian)**: Self-contained scalar-relativistic Koelling–Harmon radial solver, libxc exchange–correlation, screened structure constants, Madelung sums, and Anderson mixing without tabulated empirical parameters.
- **Sublattices & HEAs**: Handles ordered crystals, site-selective disorder, multicomponent alloys, and partially disordered Heuslers ($L2_1$) on decorated lattices with species-resolved magnetic moments and SIM electrostatics.
- **Rigorously Validated**: Every claim is backed by automated tests: hydrogenic levels to $10^{-10}$ Ry, NIST atoms to $3\times 10^{-7}$ Ha, exact $2^N$ configuration averages to $10^{-14}$, and Stuttgart LMTO benchmarks.

!!! tip "Source code repository"
    The source code is hosted on GitHub at **<https://github.com/mavensgroup/asr>**.
    To get started, clone the repository:

        git clone https://github.com/mavensgroup/asr.git
        cd asr

    See the [Installation](installation.html) guide for prerequisites and build instructions.

!!! note "How this manual and the API reference fit together"
    The pages you are reading are the *guide*: how to run, configure, interpret
    and extend the code, and the physics behind it. The **Modules**,
    **Procedures** and **Derived Types** menus in the top bar are the *reference*,
    generated from the `!>` comments in `src/`. Pages link into the reference
    with `[[module_name]]`, e.g. [[asr_recursion]].

## What can I compute?

One input format serves every system. The `[calculation] method` key selects
the solver; the structure decides whether the system is ordered or disordered.

| `method` | System | Solver | Self-consistent? | Example |
|---|---|---|---|---|
| `asr` (default) | Any decorated sc/bcc/fcc structure: pure metal, ordered compound (B2, Heusler), partial disorder, binary to high-entropy alloys, first-shell short-range order | Energy-dependent augmented-space recursion | yes (spin 1 or 2) | `example/cuzn_alloy.toml`, `fe66ni14cr20_fcc.toml`, `co2fe04cr06al_l21.toml`, `fe66ni14cr20_fcc_sro.toml` |
| `kspace` | Ordered structures of any cell, SQS supercells | Brillouin zone, tetrahedra: the exact reference | yes | `example/fe_bcc.toml`, `feco_b2_kspace.toml` |

## Reading guide

* **New user** — [Installation](installation.html) → [Quick start](quickstart.html) →
  [Convergence guide](guides/convergence.html).
* **Setting up a calculation** — [Input reference](input/index.html), [Hamiltonians](hamiltonian/index.html) and
  [Output files](output.html).
* **Want to know what the numbers mean** — [Theory](theory/index.html), [Hamiltonians](hamiltonian/index.html) and
  [Validation](validation.html).
* **Contributor** — [Architecture](architecture.html) and
  [Developing](developing.html).

## Status at a glance

!!! warning "Research code"
    Every claim in the *Validation* chapter is backed by an automated test or a
    documented measurement. Items listed under *Known limitations* there are
    open problems, not bugs to be worked around silently. In particular, the
    required recursion depth for magnetic d-band alloys is still a limiting factor.

## Conventions used in this manual

* Energies in the code and in all output files are **Rydberg**; lengths are
  **Bohr** internally. `alat` in the input is in **Å**. See
  [Units and conventions](appendix/conventions.html).
* `monospace` is a literal input key, file or module name.
* Section names such as `[scf]` refer to TOML tables; `[[species]]` is an array
  of tables.
