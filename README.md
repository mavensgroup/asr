# asr — augmented-space recursion for disordered alloys

[![GitHub](https://img.shields.io/badge/GitHub-mavensgroup%2Fasr-181717?logo=github&logoColor=white)](https://github.com/mavensgroup/asr)
[![Release](https://img.shields.io/badge/release-v0.1.0-blue.svg)](https://github.com/mavensgroup/asr/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Fortran](https://img.shields.io/badge/Fortran-2008%20%7C%202018-734f96?logo=fortran&logoColor=white)](https://fortran-lang.org)
[![Build with fpm](https://img.shields.io/badge/build-fpm-008080?logo=fortran&logoColor=white)](https://fpm.fortran-lang.org)
[![Compiler](https://img.shields.io/badge/compiler-Intel%20oneAPI%20(ifx%20%2F%20MKL)-0071C5?logo=intel&logoColor=white)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
[![libxc](https://img.shields.io/badge/libxc-≥%207.0-2e7d32)](https://libxc.gitlab.io)
[![OpenMP](https://img.shields.io/badge/parallel-OpenMP-f15a24)](https://www.openmp.org)
[![Documentation](https://img.shields.io/badge/docs-FORD-4051b5)](https://mavensgroup.github.io/asr/)

A from-scratch modern-Fortran code for electronic structure of random alloys beyond CPA: high-entropy alloys,
disordered altermagnets, short-range order.

## Method and Hamiltonian

The code has two separate layers.

* **The method: augmented-space recursion (ASR).** It takes the configuration average over the disorder and gives the
  Green function by recursion. It is Hamiltonian-agnostic. It needs an operator on a graph whose site-diagonal part
  depends on the local species and whose hopping does not. It never asks where that operator comes from.
* **The Hamiltonian: TB-LMTO-ASA.** It supplies that operator: the species parameters, the screened hopping blocks
  and their energy dependence. It is built in from scratch (own radial solver, libxc, own structure constants, own
  SCF). At present it is the only Hamiltonian. Others (Wannier, tight binding) would plug in through the same operator
  interface.

The input follows the same split: the species block and a few `[scf]` keys are Hamiltonian-specific, the rest
(structure, disorder, recursion, output) is not. See [`doc/manual/hamiltonian/index.md`](doc/manual/hamiltonian/index.md).
The boundary is not yet perfectly clean: the SCF loop calls the ASA density and Madelung steps, and the site-class
operator assumes the LMTO form h = c + d S d.

Repository: **<https://github.com/mavensgroup/asr>**

Everything is validated against independent references (exact enumeration, NIST atoms, Stuttgart LMTO
output, KPM on large supercells). `docs/STATUS.md` is the status dashboard and lists only verified claims.
`docs/ROADMAP.md` gives the plan, and `docs/FINDINGS_LOG.md` holds the detailed results.

## Dependencies

* **Intel oneAPI (`ifx` + oneMKL + OpenMP)**: Primary Fortran compiler and BLAS/LAPACK/threading runtime.
* **[fpm](https://fpm.fortran-lang.org)** (≥ 0.10): Build system and test runner (`conda install -c conda-forge fpm`, `pip install --user fpm`, or download the standalone binary from [fpm releases](https://github.com/fortran-lang/fpm/releases)).
* **[libxc](https://libxc.gitlab.io)** (≥ 7.0, C library): Exchange–correlation (`vbh`, `vwn`, `pw`, `pwmod`, `pbe`), called via its C ABI (`iso_c_binding`; no Fortran `.mod` files required). Install via `sudo dnf install libxc-devel`, `conda install -c conda-forge "libxc>=7.0"`, or CMake from source.
* **[toml-f](https://github.com/toml-f/toml-f)** (`v0.5.2`): Fetched automatically by `fpm` on first build and cached under `build/dependencies/toml-f`.

See [`doc/manual/installation.md`](doc/manual/installation.md) for platform-specific dependency installation commands (Fedora/RHEL, Ubuntu/Debian, Conda, Spack, and building `libxc` from source).

## Getting the code, building, and testing
```bash
git clone https://github.com/mavensgroup/asr.git
cd asr
source env.sh                 # Intel oneAPI (ifx, MKL) + fpm flags
asr_build --tests             # fpm build (fetches toml-f on first run)
asr_test                      # run all tests (about 40 min on an 18-core laptop; the SCF tests dominate)
```
After a module *signature* change do `rm -rf build/ifx_*` — fpm's incremental build can keep stale objects.  The typical
symptom is a segmentation fault at the entry of a routine whose source did not change (seen in `asa_scf`).

An MPI build (`-DUSE_MPI`, module `asr_mpi`) splits the recursion chains over ranks with no input change.  Without the flag
the same source gives the serial program.

## Run
```bash
asr example/cuzn_b2.toml             # ordered B2 CuZn by the energy-dependent recursion
asr example/feco_b2_asr.toml         # ferromagnetic B2 FeCo by the energy-dependent recursion
asr example/cuzn_alloy.toml          # random binary bcc Cu50Zn50 via augmented-space recursion
asr example/fe66ni14cr20_fcc.toml    # random fcc Fe66Ni14Cr20, multicomponent (Alam & Mookerjee 2009, Fig. 2)
asr example/fe66ni14cr20_fcc_sro.toml  # the same with first-shell Warren-Cowley SRO, triad (1, 1, 0) (Fig. 4)
asr example/co2fe04cr06al_l21.toml   # sublattice-disordered L2_1 Heusler Co2Fe0.4Cr0.6Al
```
Outputs are written next to the input file. Input is one validated TOML record; unknown keys are errors.
One input format for every system, solved by the energy-dependent augmented-space recursion. See [`doc/manual/quickstart.md`](doc/manual/quickstart.md); the Hamiltonian-specific input is separated in [`doc/manual/input/hamiltonian.md`](doc/manual/input/hamiltonian.md).

## What is inside (src/)
| layer | modules |
|---|---|
| method: recursion | `asr_operator` (abstract `linop`, CSR), `asr_recursion` (Lanczos/Haydock), `asr_greens` (continued fraction, tails), `asr_terminator` (sqrt, Beer–Pettifor, Luchini–Nex join, period-p), `asr_contour`, `asr_partial` (species-resolved G), `asr_kpm` (independent DOS reference) |
| method: disorder | `asr_augspace` (binary augmented space, hashed states), `asr_augn`, `asr_augcls` (multi-sublattice/multicomponent augmented space), `asr_augblk` (block kernels), `asr_augsym`, `asr_symmetry` (O_h/D_4h reduction), `asr_sublattice`, `asr_sro` (first-shell Warren–Cowley SRO), `asr_lattices`, `asr_lattice_toy` (model clusters for tests) |
| Hamiltonian: TB-LMTO-ASA | `asr_radial` (scalar-relativistic Koelling–Harmon, exact E-derivatives), `asr_xc` (libxc), `asr_atom`, `asr_potpar`, `asr_pfun` (exact ASA potential function), `asr_harmonics`, `asr_strconst` (screened S^α), `asr_lmto_ham`, `asr_lmto_alloy`, `asr_lmto_realization`, `asr_augorb`, energy-dependent forms `asr_ed`, `asr_edn`, `asr_endep`, `asr_edc` (route `edc`), `asr_ed3` (route `ed3`, full ASA), `asr_input_lmto` |
| SCF and electrostatics (joins the two) | `asr_class_scf`, `asr_alloy_scf`, `asr_mixing` (Anderson), `asr_ewald` (Madelung), `asr_screening` (screened electrostatics of the charge deviations) |
| base and driver | `asr_kinds`, `asr_testing`, `asr_tomlutil`, `asr_report`, `asr_mpi` (MPI layer, serial stubs without `-DUSE_MPI`), `asr_driver` (the one input format and driver of `asr <input.toml>`) |

## Validation highlights
* radial solver: hydrogenic levels 1e-10 Ry; free atoms He/Ne vs NIST LDA to 3e-7 Ha
* augmented-space recursion equals the exact average over all 2^N configurations (1e-14), also for the multi-orbital TB-LMTO form and multi-sublattice/multicomponent classes
* serial and 1, 2, 4 MPI ranks agree bit for bit (CuZn B2, random alloys on routes `edc` and `ed3`)

## Conventions
Rydberg atomic units; `alat` in Å in inputs, Bohr internally; orbital order s, p(y,z,x), d(xy,yz,3z²−r²,xz,x²−y²).
Stuttgart (AMscf) and Questaal are reference *code* only — no numbers from them are used as inputs.

## Machines and code sync
Heavy runs go to the desktop (`ssh desktop`, needs the SRM VPN: `nmcli con up id "SRM"`). Code moves through GitHub (<https://github.com/mavensgroup/asr>):
commit locally, then `asr_sync` (pushes `main` to GitHub and fast-forwards the desktop working copy `~/asr_remote4`) or `asr_pull` (fetches commits made on the other side).
The desktop working copy must stay clean (commit there before syncing); run outputs live in the ignored `trial/`.

