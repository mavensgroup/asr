---
title: Architecture
---

# Architecture

## Design principles

1. **Hamiltonian-agnostic core.** The recursion needs only an operator with
   `apply` ([[asr_operator]]) and, for disorder, class parameters
   ([[asr_augcls]]). The method does not know which Hamiltonian produced them.
   TB-LMTO-ASA is the one Hamiltonian implemented; a Wannier or tight-binding
   model would plug in the same way. See [Hamiltonians](hamiltonian/index.html).
2. **Nothing tabulated.** Radial solver → potential parameters → structure
   constants → SCF are all in the code; reference codes are used only as
   validation targets.
3. **Strict input.** One validated TOML record; unknown keys are errors.
4. **Every claim has a test.** See [Validation](validation.html).
5. **Matrix-free.** Augmented-space operators are applied as a sequence of
   elementary operators (random diagonal, hop) without forming a matrix.

## Layers

The code has two layers that meet only in the SCF: the **method** (augmented-space
recursion, Hamiltonian-agnostic) and the **Hamiltonian** (TB-LMTO-ASA, the one
implemented). The diagram shows the modules of the program (`app/main.f90`) only.

```text
 driver        asr_driver  asr_inputdump  asr_report  asr_runstats   app/main.f90
                     │
 SCF           asr_class_scf  asr_mixing  asr_ewald  asr_screening   (joins the two layers)
                     │                                  │
 ──── METHOD (Hamiltonian-agnostic) ────────┬──── HAMILTONIAN: TB-LMTO-ASA ───────────────────
 disorder      asr_augcls  asr_augblk        │  operator    asr_edc  asr_ed3  asr_lmto_alloy
               asr_augspace  asr_sublattice  │              asr_lmto_ham  asr_lmto_cell
               asr_sro  asr_lattices         │  physics     asr_radial  asr_potpar  asr_atom
 recursion     asr_chain  asr_contour        │              asr_xc  asr_strconst  asr_harmonics
                                             │  input       asr_input_lmto
 ──────────────────────────────────────────────────────────────────────────────────────────────
 base          asr_kinds  asr_operator  asr_symmetry  asr_tomlutil  asr_lapack  asr_mpi
```

The method sees only the abstract operator ([[asr_operator]]) and the class
parameters ([[asr_augcls]]); the Hamiltonian supplies them
([Hamiltonians](hamiltonian/index.html)). The boundary is not yet perfectly
clean: the SCF loop still calls the ASA density and Madelung steps, and
`asr_augcls` assumes the \(h=c+dSd\) structure of the LMTO operator. The ordered
k-space path ([[asr_lmto_driver]], [[asr_asa]], [[asr_tetra]]) is an internal
check and sits on the Hamiltonian side.

### The reference layer

`src/reference/` holds the earlier generations of the method and independent
numerics: the binary operators and routes, the n-species operator on one
sublattice, the binary SCF, plain Lanczos, terminators, kernel polynomials. The
program uses none of them. The tests and the probes use them to cross-check the
production code, for example the class operator with `edc` against the binary
route, or the recursion against an exact configuration average. A reference
module may use production modules; a production module never uses a reference
module. `tools/check_layers.py` checks the rule and prints the size of the
program (38 modules, about 10 700 lines; reference layer 16 modules).

## Program flow

`app/main.f90` ([[asr_main]]) reads the one input format with `read_run_input` ([[asr_driver]]) and dispatches on `[calculation] method`:

```text
asr input.toml
  └─ read_run_input(file)        strict TOML, one format
       └─ run_input
            ├─ method = "asr"     (default)  class_scf → ldos on the real axis → files + report
            └─ method = "kspace"             run_lmto → asa_scf → asa_dos → files + report
```

Exit codes: 0 success; 1 input error; 2 usage; 3 SCF not converged.

## The augmented-space pipeline in one picture

```text
species + [asr]
    │
    ▼
asr_class_scf (per SCF iteration)
    │  1. asr_radial → asr_potpar : C, √Δ, p, γ, E_ν per species and spin
    │  2. asr_lmto_alloy / asr_sublattice : graph + bond blocks S^α (common screening)
    │  3. asr_augcls : node set (hashed), point-group reduction, params → random operators
    │  4. asr_edc / asr_ed3 : one complex-symmetric Lanczos chain per energy (asr_augblk: blocks of energies)
    │  5. asr_chain : terminated continued fraction → G(z) per component → asr_contour moments M_k
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
| [[asr_lapack]] | Explicit interfaces to the LAPACK routines used |
| [[asr_mpi]] | MPI wrappers (serial stubs without `USE_MPI`) |
| [[asr_testing]] | Minimal test harness: `check`, `check_close`, `report` |
| [[asr_tomlutil]] | Sub-table access and strict key checking for the TOML readers |
| [[asr_driver]] | The one input reader and the dispatch on `method` |
| [[asr_inputdump]] | The resolved input (every key with the value used) for the report |
| [[asr_input_lmto]] | Input record of the k-space method (`lmto_input_t`, species) |

### TB-LMTO-ASA Hamiltonian and the k-space method

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
| [[asr_lmto_driver]] | Driver for `method = "kspace"` |

### Disorder, augmented space and the energy-dependent recursion

| Module | Purpose |
|---|---|
| [[asr_augspace]] | Graph of the lattice and hashed augmented-space states |
| [[asr_symmetry]] | O\(_h\) site permutations on periodic sc/bcc/fcc supercells; D\(_{4h}\) orbital representations |
| [[asr_augcls]] | Site-class (sublattice) resolved, point-group reduced augmented space, n species per class |
| [[asr_augblk]] | Kernels on blocks of contour energies |
| [[asr_lattices]] | Periodic sc/bcc/fcc graphs and random alloy realisations |
| [[asr_sublattice]] | Decoration cells: sublattices from Wyckoff positions and centering |
| [[asr_sro]] | First-shell short-range order: conditional probabilities |
| [[asr_screening]] | Screened charge-deviation electrostatics (SIM, SRO closure) |
| [[asr_lmto_alloy]] | Graph and bond blocks \(S^\alpha\) from own LMTO output (common screening) |
| [[asr_edc]] | AMscf-form energy-dependent ASR on the site-class operator (`route = "edc"`, the default) |
| [[asr_ed3]] | Full generalised ASA problem on the site-class operator (`route = "ed3"`) |
| [[asr_chain]] | Terminated continued fraction of a chain, square-root tails, adaptive depth, chain statistics |
| [[asr_contour]] | Complex-energy contour (closed, open, Fermi-smeared) and occupied moments |

### Self-consistency and output

| Module | Purpose |
|---|---|
| [[asr_class_scf]] | SCF of `method = "asr"`: any decorated structure, ordered or disordered |
| [[asr_runstats]] | Timers and the names of the total-energy terms |
| [[asr_report]] | `asr.out`, `asr.json` and the DOS file header ([Output files](output.html)) |

### Reference layer (`src/reference/`, tests and probes only)

| Module | Purpose | Production counterpart |
|---|---|---|
| [[asr_alloy_scf]] | SCF for binary alloys (SIM electrostatics, three recursion routes, B2 real-space) | [[asr_class_scf]] |
| [[asr_augorb]] | Binary multi-orbital augmented-space operator \(H=E_\nu+h-hoh\) | [[asr_augcls]] |
| [[asr_augsym]] | D\(_{4h}\)-reduced binary operator | [[asr_augcls]] |
| [[asr_augn]] | n-component augmented space on one sublattice | [[asr_augcls]] |
| [[asr_ed]] | AMscf-form energy-dependent ASR, binary | [[asr_edc]] |
| [[asr_edn]] | AMscf-form energy-dependent ASR, n species | [[asr_edc]] |
| [[asr_endep]] | Potential-function routes \(P(z)-S\), binary; re-exports [[asr_chain]] | [[asr_chain]] |
| [[asr_lmto_realization]] | One numerical alloy realisation, order-2/3 Hamiltonian (KPM/ED ground truth) | |
| [[asr_recursion]] | Lanczos/Haydock recursion with normalised vectors | [[asr_edc]] |
| [[asr_greens]] | Continued-fraction Green function, square-root tails, LDOS | [[asr_chain]] |
| [[asr_terminator]] | Asymptote estimators and joiners | [[asr_chain]] |
| [[asr_partial]] | Species-resolved Green functions, partial disorder | [[asr_edc]] |
| [[asr_kpm]] | Chebyshev kernel-polynomial DOS with stochastic trace and Jackson kernel | |
| [[asr_pfun]] | Exact ASA potential function at complex energy (KKR-ASA) | |
| [[asr_kspace]] | Exact k-space Green function of ordered cubic NN bands | |
| [[asr_lattice_toy]] | NN single-orbital clusters for tests | |

## Programs

| Program | Purpose |
|---|---|
| [[asr_main]] | The `asr` executable |
| [[term_bench]] | Terminator benchmark: for every (asymptote, join, depth) reports max \|LDOS\(_\text{rec}\) − LDOS\(_\text{k-space}\)\| over the band; kinds `sc`, `bcc`, `fcc`, `b2`. Output on stdout; redirect into `trial/`. |

## Other directories

| Directory | Content |
|---|---|
| `src/` | The library (all modules above; `fpm.toml` sets `library = true`); `src/reference/` the reference layer |
| `app/` | Executables |
| `test/` | `test_*`: the fast unit tests, run by `fpm test` |
| `test_slow/` | `test_*`: slow self-consistent and enumeration tests (fpm examples, `asr_test_slow` in `env.sh`) |
| `probe/` | `probe_*`: exploratory studies that print results and take arguments (fpm examples, `fpm run --example NAME -- ARGS`) |
| `example/` | Ready-to-run input files |
| `tools/` | `sqs.py` (SQS generator), `check_layers.py` (the layer rule), `gen_build_info.sh` (build provenance) |
| `data/reference/` | Validation data: Stuttgart Ag parameters (`ag_fcc_stuttgart`), digitised CuZn/MnCr DOS |
| `trial/` | Scratch (git-ignored except `README.md`) |
| `doc/` | This manual (mdBook in `doc/manual`; the FORD API reference is generated into `doc/ford_api`) |
