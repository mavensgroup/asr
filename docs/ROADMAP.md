# ASR (new) — Roadmap

Goal: a robust, modern Fortran augmented-space-recursion code for problems CPA
cannot handle — high-entropy alloys, disordered altermagnets, short-range
order (SRO), off-diagonal / environment disorder.

## Decisions (fixed)

| item | choice |
|---|---|
| build | `fpm`, `source /opt/intel/oneapi/setvars.sh` first |
| compiler | Intel `ifx` (primary). Keep code standard-conforming; gfortran CI is a bonus, not a gate |
| linear algebra | MKL (BLAS/LAPACK/FFT if needed) |
| XC | libxc 7.1.2 (C ABI via `iso_c_binding`) |
| parallelism | OpenMP from day one; MPI / coarrays only if measurement shows a need (decision gate at M7) |
| LMTO | take only the required parts from `../AMscf/lmto` (Stuttgart) and `../questall`; generate our own potential parameters |
| references | `../ASR_CORRECT` (algorithm), `../references/*.tex,*.pdf` (spec + gates), `../asr2` (validation targets and pitfalls only — no code copied) |
| legacy bug lessons | `../asr2/docs/STATUS.md`: weak spots were magnetic moments, absolute band positions, disordered-site electrostatics, E_F gauge drift |

## Architecture principle: Hamiltonian-agnostic core

ASR only needs three things from a model; everything physical lives behind
this interface so LMTO, Wannier and plain TB are interchangeable:

```
HamiltonianProvider
  - sites & orbital basis      : N_site, n_orb(site), lattice, neighbour list
  - per-species local data     : on-site (matrix, may be non-diagonal in L), and
                                 any disorder-dependent quantity
  - hopping blocks             : H_{RL,R'L'}(species_R, species_R'), possibly
                                 species-pair dependent (off-diagonal disorder)
```

Layers (each is an fpm library module group, no upward dependencies):

1. `kinds/consts/io/units` — base
2. `lattice` — Bravais + basis, neighbour shells, symmetry ops
3. `ham` — provider interface + implementations: `ham_tb` (toy/Slater–Koster),
   `ham_lmto` (TB-LMTO second order, h − h o h), `ham_wannier` (reads
   `*_hr.dat`)
4. `disorder` — species, concentrations, SRO (Warren–Cowley), multi-component
   (HEA) generalisation of the binary M operator
5. `augspace` — configuration states, map generation, symmetry reduction
6. `recursion` — Haydock/block Haydock, terminators (sqrt, Beer–Pettifor,
   Lucas/Lanczos-error), continued fractions
7. `gf` — Green function, DOS/LDOS, spectral function, contour integration
8. `scf` — XC, Poisson/Madelung, Fermi level, mixing (needs `ham_lmto`)
9. `drivers` — `asr` executable, TOML input

## Milestones and deliverables

Every milestone ends with: runnable example, automated test(s) that must pass,
and a short entry in `docs/STATUS.md` stating only what is *verified*.

### M0 — Skeleton  (≈ days)
- `fpm.toml`, ifx + MKL + OpenMP flags, test harness, CI-style `make check`
- **Deliverable:** `fpm test` runs on ifx; one trivial passing test.

### M1 — Recursion core on a plain TB model
- Block-capable Haydock on a generic sparse operator; sqrt terminator;
  continued fraction G(z); LDOS.
- **Test:** pure (ordered) s-band tight-binding on sc/bcc/fcc: DOS vs
  analytic/k-space BZ integration; recursion == dense diagonalisation of a
  finite cluster to 1e-12 on exhausted Krylov space.
- **Deliverable (first user-visible):** `asr dos` for an *ordered* bcc and fcc
  single-band model, plot of DOS.

### M2 — Augmented space for binary disorder (single orbital)
- Binary M operator, configuration-state map, BFS generation with a depth
  cutoff, neighbour-vector tables, OpenMP-parallel H·ψ.
- **Test:** exact configuration average by brute-force enumeration on small
  clusters (2^N configs) == ASR to machine precision when the Krylov space is
  exhausted; CPA-comparable limit for diagonal disorder.
- **Deliverable:** DOS of a binary alloy on bcc/fcc/sc TB model with diagonal
  and off-diagonal disorder, any concentration.

### M3 — Symmetry reduction and performance
- Point-group fold of the augmented space, orbital-degeneracy reduction,
  memory-lean state storage, OpenMP scaling study.
- **Test:** folded == unfolded results; scaling table recorded.
- **Deliverable:** benchmark of states/s vs threads; depth reachable on target
  hardware stated in numbers.

### T — Terminator workstream (starts at M1, gates M5; runs in parallel)
The tail is the main accuracy limit: in ASR the depth is tiny (augmented space
grows exponentially, ~10-20 levels vs hundreds for ordered systems), so
termination error is not a footnote. Legacy `ASR_CORRECT/cgreen.f90` uses
Beer-Pettifor: band edges from the finite chain with last beta halved,
iterated to a symmetric spectrum, then a single-band sqrt tail. Known
weaknesses to attack: one band only (d-band + s-p, gaps, sharp edges are
mis-modelled), edges set by the finite chain, and the oscillating
alpha_n/beta_n of disordered systems never "converging" to a constant.

Candidates, all implemented behind one `terminator` interface and compared on
the same coefficients:
1. sqrt, constant tail (baseline)
2. Beer-Pettifor (reproduce legacy first; must match `ASR_CORRECT`)
3. Nex (user wrote "nycx": assumed Nex, to be confirmed) / Luchini-Nex
4. Turchi-Ducastelle-Treglia-style multiband / gap-aware (analytic
   continuation of the coefficients with a tail that has several bands)
5. Chebyshev/KPM-with-kernel as an independent non-continued-fraction control
6. newer ideas (e.g. rational/Pade or Stieltjes-constrained tails, to be
   surveyed from the literature before implementing)

Test bench (independent truth, built before any terminator is trusted):
- ordered lattices: coefficients to depth N vs exact k-space DOS, error vs N
- 1D / small-cluster binary alloys: exhausted-Krylov exact average, then
  truncated at depth N, measure error vs N
- large-supercell ensemble exact diagonalisation (many configs, MKL) as the
  truth for 3D disordered cases
- invariants for every terminator: positivity of DOS, normalisation (sum rule,
  total weight = 1), first 2N moments reproduced, E_F and integrated charge
  stable vs depth, Kramers-Kronig consistency
- Fermi level / charge integrated on the contour must be insensitive to the
  tail (asr2 saw E_F drift with chain length)

### M4 — LMTO ingredients (own generation)
- Port only what is needed: radial solver / potential parameters (C, Δ, o, p,
  E_ν) from atomic-sphere calculation (from ATOM in `AMscf/lmto`), screened
  structure constants (TBSTR), libxc for XC.
- **Test:** potential parameters vs Stuttgart ATOM and Questaal `lm` for Cu,
  Fe, Ni, Cr; ordered-crystal band structure/DOS vs `lm`/`lmf-asa` (same
  input, same parameters).
- **Deliverable:** `asr dos` for **pure fcc Cu / bcc Fe** with real s,p,d
  orbitals, matching Questaal ASA DOS.

### M5 — Real alloys, non-self-consistent
- `ham_lmto` with species-dependent parameters (random C, √Δ, o) in the
  augmented space; full spd hop; Fermi level via energy contour.
- **Test / gate:** CuZn bcc (Tarafder 2006, Fig. 3) DOS from fixed, tabulated
  potential parameters; compare against ASR_CORRECT run on identical inputs
  (must agree to numerical noise) and against the paper.
- **Deliverable:** DOS and E_F for fcc/bcc binary alloys (CuNi, CuZn).

### M6 — Self-consistency (the hard part; legacy asr2 failed here)
- Charge per species, **correct disordered-site electrostatics** (species-
  resolved Madelung / screened-impurity model, with explicit validation vs
  supercell LDA/ASA for an ordered ↔ disordered limit), moment/charge mixing,
  Fermi-level determination without gauge drift.
- **Test / gate:** ordered limit reproduces `lm` total energy and moments;
  FeNi/FeCr moments vs literature and vs a large-supercell reference.
- **Deliverable:** self-consistent binary alloy, charges, moments, DOS.

**SPRKKR notes (`/home/rudra/Software/sprkkr`, read 2026-10):** Madelung
potential belongs to the *site*, not the atom type; the charge multipoles
`CMNTQ` of a site are the concentration-weighted sum over the types on it
(`scfcmnt.f`, `scfmad_pot.f`). Charge-transfer beyond single-site is handled
by the Screened Impurity Model: a type-dependent shift
`V_SIM = -SIM * 2 sqrt(4pi) * Q_net(type) / R_nn` with `SIM` in [0,1]
(0.5 default, `fpscfnewpot.f`, `scfnewpot.f`). NL-CPA exists (`nlcpa.f`,
Rowlands 2003) and can return a correction to the Madelung potential, with an
`alpha` short-range-order parameter on cluster configuration probabilities.
For us ASR can do better than SIM: charge per site *and per local environment*
is available from the augmented-space states, so electrostatics can be built
from environment-resolved charge rather than a screened-impurity fudge. This is
to be tested against SIM and supercell references, not assumed.

### M7 — Multi-component, SRO, HEA
- n-component generalisation of the augmented-space operator (non-binary
  local space), Warren–Cowley SRO (first shells, then general), sublattice
  disorder (B2, L1₂, etc.), HEA examples (CoCrFeNi, Cantor alloy).
- **Test:** vs enumerated supercells for small N; SRO = 0 limit identical to
  the uncorrelated map; SRO checks vs Alam 2009 FeNiCr.
- **Decision gate:** is memory/time OpenMP-bound? Decide MPI / coarrays only
  from measured data.
- **Deliverable:** HEA DOS/spectral function with and without SRO.

### M8 — Spectral functions, transport, magnetism
- Bloch spectral function A(k,E), optical conductivity (Tarafder–Mookerjee),
  exchange couplings J_ij, collinear then spin-polarised magnetism with
  disorder.
- **Deliverable:** A(k,E) and σ(ω) for a binary alloy; J(R) for FeCo etc.

### M9 — Altermagnets, non-collinear, SOC
- Spin-split bands from sublattice + spin structure under disorder, noncollinear
  (Kaphle 2012 formalism), optional SOC.
- **Deliverable:** disordered altermagnet (e.g. MnTe-like or RuO₂-type model)
  spin-split spectral function with SRO.

### M10 — Wannier front-end and release
- `ham_wannier` for arbitrary materials; user docs; examples; regression suite.

## Verification policy (to avoid repeating asr2)
- Nothing is called "working" until a test against an independent reference
  (exact diagonalisation, brute-force average, Questaal, paper figure) passes.
- `docs/STATUS.md` records only verified claims, each with its test name.
- One milestone at a time; no new feature while the previous gate is red.
- Legacy code is read as a reference, never pasted wholesale.

## Real-disorder roadmap: energy-dependent (KKR-ASA) route  (set 2026-10-07)
The production route for every system is the energy-dependent recursion on P(z) - S (user decision, see STATUS.md).  The energy-independent
H^(2)/H^(3) recursion and the k-space ASA stay as cross-checks.  Each stage has a gate; no stage starts while the previous gate is red.

### D0 - Exact potential function (no linearisation)
P_l(z), dP/dz, d2P/dz2 and the normalised phi_l(z, r) from the scalar-relativistic radial equation at each complex contour energy;
P^alpha = P^0/(1 - alpha P^0), G = lambda + mu g mu with mu^2 = dP^alpha/dz, lambda = -P''/(2P'); valence density
rho_l(r) = -(1/pi) Im sum_j w_j G_l(z_j) phi_l(z_j, r)^2 (no E_nu, no C/Delta/gamma/p in the band problem).
Gates: (a) unit: exact P agrees with the two-parameter form to second order at E_nu and with the real-axis radial P (probe_pfun) on the real axis;
mu^2, lambda vs complex finite differences; (b) ordered: bcc Fe, B2 CuZn, B2 FeCo moments / E_F / charge transfer vs k-space ASA within the
linearisation error of the k-space ASA (expect a few mRy, ~0.01-0.02 muB); (c) sum rule (all states) and charge neutrality to 1e-4 e.

### D1 - Representation and pole control for disorder
With exact P, choose screening alpha_l per l so that the real poles of P^alpha_l lie outside the occupied window of the contour for every species.
Gates: fixed-potential disordered CuZn at depth 3, 4, 5: species-resolved d occupations <= 10, total charge = nelec, M0 sum rule per orbital = 1,
and the DOS converges with depth (the current two-parameter ED over-counts at depth >= 4: Cu d 11.7 e).

### D2 - Non-magnetic disordered SCF
CuZn bcc at x_Cu = 0.25, 0.5, 0.75: SCF convergence (all three), adaptive depth (`depth_tol`), rank and node-budget study, SIM beta.
Gates: Cu-d / Zn-d peak positions and their concentration trend vs Banerjee & Mookerjee 2010 (digitised data in data/reference/); qnet trend;
results stable under depth +4 and rank +1.  Second system: CuNi or AgPd (Chakrabarti & Mookerjee 2005 moments/charges).

### D3 - Electrostatics in disorder
SIM beta from the ordered limit and from supercell averages (own k-space ASA on SQS-like supercells: same ASA model, exact average).
Gate: species charges and E_F within 0.01 e / 5 mRy of the supercell average for CuZn x = 0.5.

### D4 - Magnetic disorder
Spin-polarised ED SCF: FeCo bcc disordered (x = 0.5), FeV 50-50 (Chakrabarti & Mookerjee 2005 Table 1: Fe 1.12, V -0.11 muB), NiFe (Invar side
excluded).  Needs: robust spin mixing (moment-aware Anderson, spin-flip start control), SCF stall criterion (stop when dV, E_F and moments are
stationary).  Gates: FeCo disordered ~ FeCo ordered within literature spread (KKR-CPA, SQS); FeV signs/magnitudes; Slater-Pauling trend for FeCo
x = 0.2..0.8.

### D5 - Independent validation of the disorder average
Brute-force supercell averages (k-space ASA over many random 16-128 atom cells, same ASA model) and published KKR-CPA / SQS data.
Gate: species DOS and moments agree with the supercell average within the statistical error; documented in VALIDATION_TARGETS.md.

### D6 - Cost and scale
Node generation only to the depth actually used (adaptive), symmetry reduction kept, OpenMP efficiency on the desktop; profile per SCF iteration.
Gate: CuZn x = 0.5 depth 5 SCF in < 1 h on the desktop.

### D7 - Beyond binary
Multi-component (ternary, HEA), short-range order (Warren-Cowley), off-diagonal (size) disorder; non-collinear magnetism (M8/M9).

## Progress snapshot (kept current; details in `docs/STATUS.md`, open items in `docs/TODO.md`)
| milestone | state |
|---|---|
| M0 skeleton | done |
| M1 recursion core + ordered tight-binding DOS | done; bcc single-band tail error is a known terminator limitation |
| T terminator workstream | stage 1 done (5 tails, Luchini-Nex from description only); disordered-case measurement waits for M2 |
| M4 own LMTO | **mostly done for non-magnetic single-atom metals**: radial solver, libxc, atom SCF, potential parameters, ASA SCF reproduce Stuttgart Ag to ~1%. Next: spin polarisation, multi-atom cells + Madelung, Fermi-level integration, TOML wiring |
| M2 augmented space | binary single-orbital done and exact vs brute force (`test_augspace`); point-group reduction (M3 part) done, ~40x; first disordered terminator study measured. Next: multi-orbital / LMTO augmented Hamiltonian (M5), multi-component + SRO (M7) |
| M5/M6 alloy SCF | done for the EI route (exact augmented space, D4h reduction, SIM, ordered B2 real-space); ED route: ordered systems validated (FeCo, CuZn, Fe), disordered only at depth 3 |
| D0-D7 real disorder (ED) | D0 started 2026-10-07 |

## Immediate next steps
See `docs/TODO.md` (ordered list, updated every session).
