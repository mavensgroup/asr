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
disordered altermagnets, short-range order. The LMTO part is built in (own radial solver, libxc, own structure
constants, own SCF); other Hamiltonians (Wannier, tight binding) plug in through the same operator interface.

Repository: **<https://github.com/mavensgroup/asr>**

Everything is validated against independent references (exact enumeration, k-space, NIST atoms, Stuttgart LMTO
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
asr_test                      # run all tests (≈ 3 min on 16 cores)
```
After a module *signature* change do `rm -rf build/ifx_*` — fpm's incremental build can keep stale objects.

## Run
```bash
asr example/fcc_ordered.toml         # model tight-binding s-band via Haydock recursion (1 s)
asr example/ag_fcc.toml              # self-consistent fcc Ag: SCF → parameters → l-resolved DOS (9 s)
asr example/fe_bcc.toml              # ferromagnetic bcc Fe
asr example/fecu_b2.toml             # ordered B2 FeCu (two classes, Madelung, charge transfer)
asr example/cuzn_alloy.toml          # random binary bcc Cu50Zn50 via augmented-space recursion
asr example/co2fe04cr06al_l21.toml   # sublattice-disordered L2_1 Heusler Co2Fe0.4Cr0.6Al
```
Outputs are written next to the input file. Input is one validated TOML record; unknown keys are errors.
Supported `[calculation] kind` values: `"model"`, `"lmto"`, `"alloy"`, and `"sublattice"` (see [`doc/manual/quickstart.md`](doc/manual/quickstart.md)).

## What is inside (src/)
| layer | modules |
|---|---|
| base | `asr_kinds`, `asr_operator` (abstract `linop`, CSR), `asr_testing`, `asr_tomlutil`, `asr_input`, `asr_input_lmto`, `asr_report` |
| recursion | `asr_recursion` (Lanczos/Haydock), `asr_greens` (continued fraction, tails), `asr_terminator` (sqrt, Beer–Pettifor, Luchini–Nex join, period-p), `asr_contour`, `asr_partial` (species-resolved G) |
| LMTO | `asr_radial` (scalar-relativistic Koelling–Harmon, exact E-derivatives), `asr_xc` (libxc), `asr_atom`, `asr_potpar`, `asr_harmonics`, `asr_strconst` (screened S^α), `asr_lmto_ham`, `asr_lmto_cell`, `asr_ewald`, `asr_tetra`, `asr_mixing` (Anderson), `asr_asa` (ordered SCF, spin, multi-class), `asr_lmto_driver` |
| disorder | `asr_augspace` (binary augmented space, hashed states), `asr_augsym`, `asr_augn`, `asr_augcls` (multi-sublattice/multicomponent augmented space), `asr_symmetry` (O_h/D_4h reduction), `asr_augorb`, `asr_ed`, `asr_endep`, `asr_edc`, `asr_sublattice`, `asr_lmto_alloy`, `asr_lmto_realization`, `asr_alloy_scf`, `asr_class_scf`, `asr_alloy_driver`, `asr_class_driver`, `asr_kpm`, `asr_lattices`, `asr_kspace` |

## Validation highlights
* radial solver: hydrogenic levels 1e-10 Ry; free atoms He/Ne vs NIST LDA to 3e-7 Ha
* own ASA SCF reproduces Stuttgart Ag parameters (√Δ, γ, p) to ~1%; Fe 2.2 μB, Ni 0.64 μB, Cu 0; B2 FeCu occupations to 0.01–0.02 e
* augmented-space recursion equals the exact average over all 2^N configurations (1e-14), also for the multi-orbital TB-LMTO form and multi-sublattice/multicomponent classes
* real-space augmented-space recursion equals the ordered k-space recursion coefficients (1e-14)
* alloy SCF in the trivial-alloy limit reproduces the ordered SCF

## Conventions
Rydberg atomic units; `alat` in Å in inputs, Bohr internally; orbital order s, p(y,z,x), d(xy,yz,3z²−r²,xz,x²−y²).
Stuttgart (AMscf) and Questaal are reference *code* only — no numbers from them are used as inputs.

## Machines and code sync
Heavy runs go to the desktop (`ssh desktop`, needs the SRM VPN: `nmcli con up id "SRM"`). Code moves through GitHub (<https://github.com/mavensgroup/asr>):
commit locally, then `asr_sync` (pushes `main` to GitHub and fast-forwards the desktop working copy `~/asr_remote4`) or `asr_pull` (fetches commits made on the other side).
The desktop working copy must stay clean (commit there before syncing); run outputs live in the ignored `trial/`.

