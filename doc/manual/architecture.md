---
title: Architecture
---

# Architecture

## Design principles

1. **Hamiltonian-agnostic core.** The recursion needs only an operator with
   `apply` ([[asr_operator]]). LMTO, tight-binding and augmented-space operators
   all extend the abstract `linop`; Wannier models can plug in the same way.
2. **Nothing tabulated.** Radial solver → potential parameters → structure
   constants → SCF are all in the code; reference codes are used only as
   validation targets.
3. **Strict input.** One validated TOML record; unknown keys are errors.
4. **Every claim has a test.** See [Validation](validation.html).
5. **Matrix-free.** Augmented-space operators are applied as a sequence of
   elementary operators (random diagonal, hop) without forming a matrix.

## Layers

The layering below is the intended structure (the ROADMAP's "no upward dependencies" principle); it is conceptual and not enforced by the build.

```text
 drivers      asr_lmto_driver  asr_alloy_driver  asr_class_driver  asr_report      app/main.f90
                     │                 │                 │
 SCF          asr_asa          asr_alloy_scf      asr_class_scf
                     │                 │                 │
 operators    asr_lmto_cell    asr_augorb/augsym/augn/augcls   asr_endep/ed/edn/edc
              asr_lmto_ham     asr_lmto_alloy  asr_lmto_realization  asr_sublattice
                     │                 │
 physics      asr_radial  asr_potpar  asr_atom  asr_xc  asr_strconst  asr_harmonics
              asr_ewald   asr_tetra   asr_mixing  asr_pfun
                     │
 recursion    asr_recursion  asr_greens  asr_terminator  asr_contour  asr_partial
                     │
 base         asr_kinds  asr_operator  asr_symmetry  asr_lattices  asr_tomlutil  asr_testing
```

## Program flow

`app/main.f90` ([[asr_main]]) reads `[calculation] kind` with `peek_kind` and dispatches:

```text
asr input.toml
  └─ peek_kind(file)
       ├─ "lmto"       read_lmto_input   → run_lmto   → asa_scf → asa_dos  → files + report
       ├─ "alloy"      read_alloy_input  → run_alloy  → alloy_scf → alloy_ldos → files + report
       ├─ "sublattice" read_class_input  → run_class  → class_scf → files + report
       └─ otherwise    read_input (model) → cubic_cluster → lanczos → cf_green → dos file
```

Exit codes: 0 success; 1 input error; 2 usage; 3 SCF not converged.

## The augmented-space pipeline in one picture

```text
species + [asr]
    │
    ▼
asr_alloy_scf / asr_class_scf (per SCF iteration)
    │  1. asr_radial → asr_potpar : C, √Δ, p, γ, E_ν per species and spin
    │  2. asr_lmto_alloy / asr_sublattice : graph + bond blocks S^α (common screening)
    │  3. asr_augsym / asr_augcls : node set (hashed), D4h reduction, params → random operators
    │  4. recursion: Lanczos (energy-independent: stored basis; ED: one chain per z)
    │  5. asr_partial / asr_endep : species-resolved G(z) → asr_contour moments M_k
    │  6. density, core, Hartree, xc, SIM/Madelung → asr_mixing
    ▼
converged potentials → final DOS on the real axis
```

## Module index

All modules link to their generated API pages.

### Base

| Module | Purpose |
|---|---|
| [[asr_kinds]] | Kind parameters (`dp`, `i4`, `i8`), \(\pi\), `ci` |
| [[asr_operator]] | Abstract Hermitian operator `linop` (`apply`, `dim`, `mapply`), CSR matrix |
| [[asr_testing]] | Minimal test harness: `check`, `check_close`, `report` |
| [[asr_tomlutil]] | Sub-table access and strict key checking for the TOML readers |
| [[asr_input]] | Reader for `kind = "model"` |
| [[asr_input_lmto]] | Reader for `kind = "lmto"`, `peek_kind` |

### Recursion core

| Module | Purpose |
|---|---|
| [[asr_recursion]] | Lanczos/Haydock recursion with normalised vectors |
| [[asr_greens]] | Continued-fraction Green function, square-root tails, LDOS |
| [[asr_terminator]] | Asymptote estimators (last, avg2, Beer–Pettifor, period 2) and joiners (butt, Luchini–Nex) |
| [[asr_contour]] | Complex-energy contour (closed, open, Fermi-smeared) and occupied moments |
| [[asr_partial]] | Species-resolved Green functions from the augmented-space recursion |

### Ordered LMTO-ASA

| Module | Purpose |
|---|---|
| [[asr_radial]] | Scalar-relativistic Koelling–Harmon solver on a shifted log mesh, exact energy derivatives |
| [[asr_xc]] | LDA/GGA exchange–correlation via the libxc C ABI, spin-polarised |
| [[asr_atom]] | All-electron free-atom LDA; cumulative integrals and Poisson–Hartree |
| [[asr_potpar]] | ASA potential parameters \(E_\nu, C, \sqrt\Delta, p, \gamma\) from the radial solutions |
| [[asr_harmonics]] | Real spherical harmonics, Gaunt coefficients by exact quadrature |
| [[asr_strconst]] | Bare and screened TB-LMTO structure constants |
| [[asr_lmto_ham]] | Second-order TB-LMTO-ASA Hamiltonian for an ordered one-atom cell |
| [[asr_lmto_cell]] | General ordered cell: clusters, screening, \(H(\mathbf k)\), DOS |
| [[asr_ewald]] | Ewald lattice sums for point charges (Madelung) |
| [[asr_tetra]] | Linear tetrahedron method (Blöchl–Jepsen–Andersen) |
| [[asr_mixing]] | Anderson mixing with safeguard |
| [[asr_asa]] | Self-consistent ordered TB-LMTO-ASA: multi-class, collinear spin |
| [[asr_lmto_driver]] | Driver for `kind = "lmto"` |

### Disorder and augmented space

| Module | Purpose |
|---|---|
| [[asr_augspace]] | Binary augmented space (single orbital), hashed states |
| [[asr_symmetry]] | O\(_h\) site permutations on periodic sc/bcc/fcc supercells; D\(_{4h}\) orbital representations |
| [[asr_augorb]] | Multi-orbital augmented-space operator \(H=E_\nu+h-hoh\) |
| [[asr_augsym]] | D\(_{4h}\)-reduced multi-orbital augmented-space operator |
| [[asr_augn]] | n-component (multicomponent/HEA) augmented space |
| [[asr_augcls]] | Site-class (sublattice) resolved, point-group reduced augmented space |
| [[asr_lattices]] | Periodic sc/bcc/fcc graphs and random alloy realisations |
| [[asr_sublattice]] | Decoration cells: sublattices from Wyckoff positions and centering |
| [[asr_lmto_alloy]] | Binary-alloy description from own LMTO output (common screening) |
| [[asr_lmto_realization]] | One numerical alloy realisation, order-2/3 Hamiltonian (KPM/ED ground truth) |

### Energy-dependent routes

| Module | Purpose |
|---|---|
| [[asr_pfun]] | Exact ASA potential function at complex energy (KKR-ASA) |
| [[asr_endep]] | Energy-dependent recursion \(P(z)-S\) (complex-symmetric Lanczos per energy) |
| [[asr_ed]] | AMscf-form energy-dependent ASR (binary) |
| [[asr_edn]] | AMscf-form energy-dependent ASR, n species |
| [[asr_edc]] | AMscf-form energy-dependent ASR on the site-class operator |

### Self-consistent alloy drivers and reporting

| Module | Purpose |
|---|---|
| [[asr_alloy_scf]] | SCF for binary alloys (SIM electrostatics, three recursion routes, B2 real-space) |
| [[asr_alloy_driver]] | Reader and driver for `kind = "alloy"` |
| [[asr_class_scf]] | SCF for sublattice/multicomponent alloys |
| [[asr_class_driver]] | Reader and driver for `kind = "sublattice"` |
| [[asr_report]] | `asr.out` writer: banner, input echo, SCF history, results, timing |

### Reference methods and toys

| Module | Purpose |
|---|---|
| [[asr_kpm]] | Chebyshev kernel-polynomial DOS with stochastic trace and Jackson kernel |
| [[asr_kspace]] | Exact k-space Green function of ordered cubic NN bands |
| [[asr_lattice_toy]] | NN single-orbital clusters for the `model` kind and tests |

## Programs

| Program | Purpose |
|---|---|
| [[asr_main]] | The `asr` executable |
| [[term_bench]] | Terminator benchmark: for every (asymptote, join, depth) reports max \|LDOS\(_\text{rec}\) − LDOS\(_\text{k-space}\)\| over the band; kinds `sc`, `bcc`, `fcc`, `b2`. Output on stdout; redirect into `trial/`. |

## Other directories

| Directory | Content |
|---|---|
| `src/` | The library (all modules above; `fpm.toml` sets `library = true`) |
| `app/` | Executables |
| `test/` | `test_*` (automatic tests, run by `fpm test`) and `probe_*` (exploratory studies, not assertions) |
| `example/` | Ready-to-run input files |
| `tools/` | `sqs.py`, a SQS generator |
| `data/reference/` | Validation data: Stuttgart Ag parameters (`ag_fcc_stuttgart`), digitised CuZn/MnCr DOS |
| `docs/` | Working notes: `STATUS.md` (verified claims), `TODO.md`, `ROADMAP.md`, `VALIDATION_TARGETS.md`, findings |
| `trial/` | Scratch (git-ignored except `README.md`) |
| `doc/` | This manual (FORD pages in `doc/manual`) |
