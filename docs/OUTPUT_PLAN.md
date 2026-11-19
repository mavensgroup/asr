# Output architecture for ASR

> **Status**: plan, 2026-11-17.  Output code is not changed yet.
> **Scope**: how ASR writes what it has computed.  The physics modules already hold the data.  The data that was missing
> (energy terms, timers, SCF statistics, DOS at E_F) was added on branch `output-data` (2d998ac).  Nothing in this plan
> changes a computed number.
> **Modules**: `src/asr_report.f90` (new layout), `src/asr_driver.f90` and `src/asr_lmto_driver.f90` (fill the report record),
> one new generated module for build provenance.

---

## 1. Principles

1. **Every printed number is a stored number or a stated formula of stored numbers.**  The report never computes physics.
   Derived quantities (gamma_e, chi_P, Warren-Cowley alpha) print their formula and constants once.
2. **One unit system.**  Energy Ry, length Bohr, charge e, moment mu_B, DOS states/Ry.  Every column header carries its unit.
   The banner gives the conversions once.  eV and Angstrom appear only as a second value in parentheses where a reader
   compares with experiment (lattice constant, N(E_F)).
3. **Spin convention stated, not assumed.**  The SCF results hold several quantities per spin channel.  The report prints
   electrons and states for both spins unless the column says `per spin`.  Section 3 lists the convention of every field.
4. **Width at most 90 columns, no merged numbers.**  Every numeric field has at least one blank before it.  Wide tables
   (many components) are split into blocks of columns, never wrapped.  A unit test checks the width of every line.
5. **Greppable.**  Sections start with `=== name ===`.  The run ends with exactly one status line (Section 4.13).
6. **Health checks are a table with status words.**  `OK`, `WARN`, `FAIL`, `info`.  A collected warning list at the end
   repeats every `WARN` and `FAIL`.
7. **A machine-readable copy.**  Every run writes `<out>.json` with the same final numbers (Section 5).  Scripts and paper
   tables read this file, never the text report.
8. **Provenance.**  Git revision and dirty flag, compiler and flags, host, ranks and threads, and an input fingerprint.
9. **The recursion route comes first.**  It is the production route.  The k-space route uses the same layout and marks
   itself as an internal check.

---

## 2. What other codes do well, and what ASR takes from them

| Code | Feature worth copying | ASR |
| :--- | :--- | :--- |
| VASP | `vasprun.xml`: a machine-readable copy of the run beside the text output | `<out>.json` (Section 5) |
| VASP `OUTCAR` | "General timing and accounting" block at the end | timing table, Section 4.11 |
| Quantum ESPRESSO | timing per routine with calls, CPU and wall time | timing table with calls, CPU, wall, CPU/wall |
| Quantum ESPRESSO | energy split into one-electron, Hartree, xc, Ewald terms that add to the total | energy card, Section 4.7 |
| ABINIT | counts of comments and warnings at the end; YAML documents inside the main output | warning list, JSON copy |
| Gaussian | one fixed final line (`Normal termination`) that scripts grep | one status line |
| SPR-KKR | per-component (type) charges, moments, DOS in CPA | everything per component, weighted sums per site |
| KKR-CPA codes generally | thin diagnostics of the effective medium | causality, contour bottom, E_nu window, depth: Section 4.10 |

ASR adds what these codes cannot report: diagnostics of the recursion itself (non-causal contour values, weight below the
contour bottom, recursion depth used, terminator fallbacks) and the SRO screening modes.

---

## 3. Data inventory

Status: **old** = in the SCF results before this plan; **new** = added on `output-data` (2d998ac); **input** = from the
TOML input; **report** = computed in the report from the fields named.  Recursion = `class_scf_result_t`
(`asr_class_scf`); k-space = `asa_result_t` (`asr_asa`).

| Item | Recursion field | k-space field | Status | Convention |
| :--- | :--- | :--- | :--- | :--- |
| E_tot | `etot`, `hist_etot` | `etot`, `hist_etot` | old | per site / per atom, Ry |
| energy terms | `eterm(n_et)`, `ecomp(n_et, comp)` | `eterm`, `ecomp(n_et, class)` | new | sum = etot to rounding |
| term names | `asr_runstats:et_name` | same | new | |
| dV, E_F per iteration | `hist_dv`, `hist_ef` | `hist_dv`, `hist_ef` | old | |
| moment, qnet per iteration | `hist_moment`, `hist_qnet` (iter, comp) | (iter, class) | old | per atom of the component |
| E_F passes, max dE_nu | `hist_npass`, `hist_denu` | none | new | |
| wall time per iteration | `hist_wall` | `hist_wall` | new | s, this rank |
| timers | `timer` (`run_timer_t`) | `timer` | new | exclusive buckets, `tb_name` |
| chain evaluations | `timer%work(tb_*)` | none | new | all ranks |
| potential parameters | `comp(p)%pp(l, s)`, `comp(p)%enu` | `cls(c)%pp`, `cls(c)%enu` | old | gamma representation |
| moments q0 q1 q2 | `comp(p)%q(k, l, s)` | `cls(c)%q(k, l, s)` | old | **recursion: per spin channel; k-space: both spins when nspin = 1** |
| channel charge Q_l | `qch(comp, s, l)` | = q0 + p q2 (report) | old | per spin channel (ed3: q0 + p q2) |
| sphere qnet, moment | `comp(p)%qnet`, `%moment` | `cls(c)%qnet`, `%moment` | old | per atom |
| DOS at E_F | `dosef(comp, s, 0:2)` | DOS grid (report interpolates) | new | per spin channel, at E_F + i ed_eta |
| DOS curve | `edos_e`, `edos(e, l, comp, s)` | `dos` of `asa_dos` | old | per spin channel |
| cell charge residual | `qcell` | report: sum natc qnet | new | e per site |
| E_nu window clamps | `nclamp`, `nclamp_run` | `nclamp`, `nclamp_run` | old | |
| weight below contour | `wbelow1`, `wbelow`, `below_ok`; `wbelow_tol` | none | old | e per site |
| causality | `nnode`, `nnoncausal` | none | old | |
| recursion depth, bp fallbacks | `depth_mean`, `nfall` | none | new | last iteration |
| Anderson restarts | `mix_nrestart`, `mix_nh` | none | new | |
| band energy check | none | `eband`, `sumev` | old | must be equal |
| screening model | `screen` (a, g00, g1, w00, lambda, beta_max) | none | old | |
| SRO modes | `mode_lam`, `mode_gam`, `fscr`, `udrive` | none | old | |
| SRO probabilities | `inp%ci%sro_cond`, `inp%ci%conc` | none | input | alpha = 1 - P(j given i) / c_j |
| nodes per sublattice | `nodes` | none | old | |
| components | `comp_cls`, `comp_spec`, `comp_w`, `comp_x` | `natc` (driver) | old | |
| provenance | `asr_build_info` (generated), `compiler_version()` | same | **to add** (output side) | Section 6 |

The report converts per-spin-channel values to electrons and states for both spins with fs = 2 / nspin.  This
conversion lives in one adapter per route (Section 7), so the report routines never see the difference between the routes.

---

## 4. The text report (`asr.out`)

Numbers in the mock-ups show the layout.  The energy card (4.7) has real numbers: fcc Ag, k-space route, from
`test_asa_ag`.  The other values are illustrative.  The CuZn values are not results.

### 4.1 Banner and provenance

```text
==========================================================================================
 ASR  augmented-space recursion for disordered alloys  *  TB-LMTO-ASA
==========================================================================================
 program            : asr 0.1.0   git 2d998ac   (working tree clean)
 compiler           : Intel(R) Fortran Compiler 2025.x
 flags              : -O2 -g -traceback -qopenmp -qmkl=parallel -free -standard-semantics
                      -warn all -fp-model=precise -heap-arrays 1024
 started            : 2026-11-17 12:00:14   host desktop   pid 41233
 parallel           : MPI ranks 1, OpenMP threads 18 per rank
 input              : example/cuzn_alloy.toml   fingerprint 9c41d0a7e2b35f18
 working directory  : /home/user/asr_remote4
 units              : energy Ry, length Bohr, charge e, moment mu_B, DOS states/Ry
 conversions        : 1 Ry = 13.605693123 eV = 0.5 Ha;  1 Bohr = 0.529177211 Angstrom
==========================================================================================
```

The fingerprint is FNV-1a 64 of the input bytes.  It identifies the input; it is not a security hash.  `dirty` replaces
`clean` when the build had uncommitted changes.

### 4.2 Input (verbatim)

As now (`rep_echo_input`): every input line behind `   | `.  It stays: a report must reproduce its run on its own.

### 4.3 Method

```text
=== method ===============================================================================
 Hamiltonian        : TB-LMTO-ASA, energy-dependent ASR, route edc (H^(2) form)
 basis              : s p d, scalar relativistic; parameters in the gamma representation
 XC functional      : vbh (libxc), LDA;  spin polarised: no
 components         : 2 (sublattice, species) pairs, one potential each
 augmented space    : levels 4, recursion steps 4, max rank 1, supercell 8
 nodes / sublattice : A 3517
 S^alpha            : kept range 1.200 alat, inversion radius 1.200 alat
 screening alpha    : s 0.34850   p 0.05303   d 0.01071
 terminator         : hull;  depth tolerance 1.0E-05 (relative)
 occupation contour : 64 nodes (Fermi function, kT 1.0 mRy, 8 Matsubara poles)
                      bottom 1.200 Ry below the lowest E_nu;  E_F probe eta 0.010 Ry
 electrostatics     : Madelung of the sublattice averages (Ewald, decoration cell)
                      + screened charge deviations (SIM, sim_beta 0.700)
 screening model    : hardness a 1.2345 Ry/e   G00 0.4567   G1 -0.0123   W00/G00 0.8123
 mixing             : Anderson, beta 0.300, history 6
                      step limits 0.15 Ry (V), 0.08 Ry (E_nu)
 convergence        : max over r of r/(1+r) |V_out - V_in| < 1.0E-05 Ry
                      (with SRO also a |u_out - u_in|);  E_nu inside its band window
```

The k-space route prints `route kspace (internal check)`, the k mesh, the BZ method and the smearing instead of the
augmented-space lines.

### 4.4 Structure

```text
=== structure ============================================================================
 graph lattice           : bcc   a = 5.385719 Bohr (2.850000 Angstrom)
 decoration cell         : 1^3 conventional cubes, 2 sites
 sphere radius S         : 2.651780 Bohr (equal spheres, space filling)
 volume per site         : 78.11163 Bohr^3 (11.57494 Angstrom^3)

 sublattice  site fraction  species and concentrations
 ---------------------------------------------------------------------
 A                  1.0000  Cu 0.5000   Zn 0.5000
```

Ordered compounds list their sublattices the same way (one species, concentration 1).

### 4.5 Short-range order (only with `[[sro]]`)

```text
=== short-range order ====================================================================
 first shell: centre on A, neighbours on A, 8 sites
   centre    P(Cu|centre)  P(Zn|centre)   alpha(Cu)   alpha(Zn)
   --------------------------------------------------------------
   Cu              0.6000        0.4000     -0.2000      0.2000
   Zn              0.4000        0.6000      0.2000     -0.2000
 alpha_ij = 1 - P(j|i) / c_j.  For i /= j: alpha > 0 means fewer unlike neighbours
 (clustering), alpha < 0 more (ordering).

 screening modes of the SRO closure
   mode    lambda(Ahat)    Gamma = dq/db    G00 + G1 lambda    F [Ry/e]
   -------------------------------------------------------------------
      1          0.0000          0.41230            0.45670      1.4310
      2          1.6000          0.20112            0.43702      3.9876
 drive u [e] (last iteration):  Cu:A  0.01234   Zn:A -0.01234
```

A mode with Gamma = 0 prints `infinite (full clustering: this deviation vanishes)` in the F column, as now.

### 4.6 SCF history

Main table, one line per iteration:

```text
=== SCF history ==========================================================================
 iter       E_tot [Ry]      dE [Ry]     dV [Ry]  dE_nu [Ry]      E_F [Ry]  pass  wall [s]
 ---------------------------------------------------------------------------------------
    1   -3329.12345678          ---   1.044E+00   8.000E-02   -0.27936800     5     12.31
    2   -3329.13001234   -6.556E-03   3.211E-01   8.000E-02   -0.10123400     3      8.02
  ...
   19   -3329.13456789   -2.100E-08   8.100E-06   1.200E-06   -0.03543366     1      6.75
 ---------------------------------------------------------------------------------------
 converged after 19 iterations (dV < 1.0E-05 Ry, no E_nu held at its window)
```

`pass` is the number of contour evaluations of the E_F search (`hist_npass`).  `dE_nu` is the largest change of E_nu
proposed by the output moments, before mixing (`hist_denu`).

Second table, per component, in blocks of at most four components (spin polarised: moment and qnet; non-magnetic:
qnet only, eight per block).  The last column of the first block is the moment per site.

```text
 per component:  m [mu_B per atom]  qnet [e per atom]
 iter    Fe:A m  Fe:A qnet     Ni:A m  Ni:A qnet     Cr:A m  Cr:A qnet   site m
 ---------------------------------------------------------------------------------
    1    2.1034   -0.01234     0.6012    0.02345    -0.8123    0.03456   1.2341
```

### 4.7 Total energy

```text
=== total energy =========================================================================
 Kohn-Sham form: output density in the input potential; no Fermi smearing entropy
 term                                                            per atom [Ry]
 ------------------------------------------------------------------------------
 valence band energy        sum_l (q1 + E_nu Q_l)                  -5.34589914
 core eigenvalue sum                                            -6386.65436378
 double counting            -int n (v_in + 2Z/r)                -8267.89945822
 sphere Hartree             1/2 int n V_H                        4326.11556616
 exchange-correlation       int n eps_xc                         -290.29543954
 Madelung                   1/2 q V_M                              -0.00000000
 screened deviations        1/2 q V_F                               0.00000000
 ------------------------------------------------------------------------------
 total energy E_tot                                            -10624.07959452
```

These are the real terms of fcc Ag (k-space, `test_asa_ag`).  Their sum agrees with `etot` to 0 in that test.  The nuclear
terms -2Z int n/r of the kinetic energy and of the sphere electrostatics cancel, so the Hartree term holds the electrons
only.  The large Hartree and double-counting terms nearly cancel.  This is expected for an all-electron sphere.

For alloys a second table gives the same rows per component (per atom of the component), in blocks of four components,
with the last column `per site = sum_p w_p`.  The report prints no mixing or formation energy.  Those need reference runs
that ASR does not have in one run.

### 4.8 Converged results per component

```text
=== converged results ====================================================================
 E_F                       :  -0.03543366 Ry
 cell charge residual      :  -1.775E-05 e per site   (sum_p w_p qnet_p; 0 = neutral)
 spin moment per site      :   0.00000 mu_B

 component Cu:A   Z 29  valence 11   c 0.5000   weight per site 0.5000   S 2.65178 Bohr
   electrons in sphere 28.98765 e   qnet -0.01235 e   moment 0.00000 mu_B

   potential parameters, gamma representation [Ry]
   spin  l        E_nu           C   sqrt(Delta)          p       gamma
   -----------------------------------------------------------------
      1  s   -0.439626   -0.342362    -0.405310    0.057631    0.421712
      1  p   -0.392049    0.673705     0.398889    0.027060    0.112814
      1  d   -0.420481   -0.325525     0.092836    2.280993    0.004528

   moments and DOS at E_F, both spins
   spin  l    q0 [e]   q1 [e Ry]   q2 [e Ry^2]    Q_l [e]   N_l(E_F) [1/Ry]
   ---------------------------------------------------------------------
      1  s   0.68472     0.00000       0.01988    0.68472          0.3576
      1  p   0.66416    -0.00000       0.01603    0.66416          0.4770
      1  d   9.63877     0.00001       0.06749    9.63877          3.7234
   sum                                            10.98765          4.5580
```

The two tables replace the present 100-column table.  Q_l is the sphere charge of the channel (edc: q0; ed3:
q0 + p q2).  With nspin = 2 the rows repeat for spin 2, and a third line per l gives up minus down.

### 4.9 Fermi level card

```text
=== Fermi level ==========================================================================
 E_F                                       :  -0.03543366 Ry
 N(E_F) per site, both spins               :   3.3002 states/Ry   (0.24256 states/eV)
 broadening of N(E_F)                      :   0.0100 Ry  (G at E_F + i eta)
 Sommerfeld coefficient, bare band          :   0.5718 mJ/(mol K^2)
 Pauli susceptibility, bare (no Stoner)     :   7.841E-06 emu/mol
 gamma_e = (pi^2/3) k_B^2 N(E_F) = 0.173247 mJ/(mol K^2) per state/Ry;
 chi_P = mu_B^2 N(E_F) = 2.37605E-06 emu/mol per state/Ry;  mol = mole of sites

 component      s [1/Ry]   p [1/Ry]   d [1/Ry]   total [1/Ry]   weight
 ---------------------------------------------------------------------
 Cu:A             0.3576     0.4770     3.7234         4.5580   0.5000
 Zn:A             0.4084     1.5090     0.1250         2.0424   0.5000
 ---------------------------------------------------------------------
 per site (sum_p w_p N_p)                              3.3002
```

The recursion route takes N(E_F) from `dosef` (the same G as the occupation contour).  The k-space route takes it from the
DOS grid and prints the smearing.  Check of the constants: fcc Ag has N(E_F) = 2 x 1.7861 states/Ry.  This gives gamma_e =
0.619 mJ/(mol K^2) (measured 0.646) and chi_P = 8.49E-06 emu/mol.  The bare values omit the electron-phonon and Stoner
enhancements; the card says so.

### 4.10 Health checks

```text
=== health checks ========================================================================
 check                                          value            limit         status
 ------------------------------------------------------------------------------------
 SCF converged (dV)                             8.1E-06 Ry       1.0E-05 Ry    OK
 E_nu held at its window, last / whole run      0 / 0            0             OK
 non-causal contour values (-Im G <= 0)         0 of 4480        0             OK
 valence weight below contour, first / final    2E-09 / 1E-09    1.0E-03 e     OK
 cell charge residual                           -1.8E-05 e       1.0E-03 e     OK
 recursion depth used, mean                     4.00 of 4        -             info
 bp tail fallbacks, last iteration              0                -             info
 Anderson restarts in the run                   2                -             info
 E_F passes, last iteration                     1                -             info
```

Limits come from the code (`tol`, `wbelow_tol`) and are printed, not hidden.  A non-causal count above zero is `WARN`.
The SCF and the E_nu window are `FAIL` when not met, as `converged` is.  The k-space route adds `band energy check
|eband - sumev|` and drops the recursion rows.

### 4.11 Timing

```text
=== timing ===============================================================================
 part                                   calls    wall [s]     cpu [s]  cpu/wall   share
 ------------------------------------------------------------------------------------
 setup (graph, operators, screening)        2        2.10        2.05       1.0    1.1%
 potential parameters                      22        0.41        0.41       1.0    0.2%
 recursion chains: SCF contour             48      171.20     3051.33      17.8   91.3%
 recursion chains: bottom check             2        3.52       62.10      17.6    1.9%
 recursion chains: DOS                      1        6.80      120.40      17.7    3.6%
 density, core, Hartree, xc, energy        19        1.92        1.92       1.0    1.0%
 mixing and E_nu windows                   19        0.05        0.05       1.0    0.0%
 restart file                              19        0.09        0.09       1.0    0.0%
 other (outside the parts above)                     1.43        1.40               0.8%
 ------------------------------------------------------------------------------------
 total                                             187.52     3239.75      17.3  100.0%
 chain evaluations, SCF contour: 1843200   (10766 per wall second)
 rank 0 of 1; CPU time is the process time of this rank (all threads)
```

`cpu/wall` is not a parallel efficiency: busy-waiting threads count as CPU.  The report does not print an efficiency.  With
MPI the table is that of rank 0; chain counts are totals over all ranks.

### 4.12 Warnings and files

```text
=== warnings =============================================================================
 none
=== files written ========================================================================
 cuzn_alloy_dos.dat          DOS, 401 energies
 cuzn_alloy_pot.dat          restart potentials (last iteration)
 asr.json                    machine-readable summary
```

Every `WARN` and `FAIL` of 4.10 and every warning printed during the run (E_F search not converged, states below the
contour) appears here with its section.

### 4.13 End of run

```text
=== end of run ===========================================================================
 finished           : 2026-11-17 12:03:22   total wall time 187.52 s
 status             : NORMAL TERMINATION   SCF CONVERGED
==========================================================================================
```

The status line is one of `SCF CONVERGED`, `SCF NOT CONVERGED`, `SCF CONVERGED WITH WARNINGS`.  It is the last non-rule
line, so `tail -2` finds it.

---

## 5. Machine-readable summary (`<out>.json`)

Plain JSON, written by rank 0 next to the report.  Every language reads JSON without extra packages (Python `json`, `jq`).
Arrays keep the order of the report.  Units are in the key names.

```json
{
  "schema": "asr-summary-1",
  "provenance": {"version": "0.1.0", "git": "2d998ac", "dirty": false, "compiler": "...",
                 "host": "desktop", "mpi_ranks": 1, "omp_threads": 18,
                 "input": "example/cuzn_alloy.toml", "input_fnv64": "9c41d0a7e2b35f18"},
  "method": {"route": "edc", "levels": 4, "nrec": 4, "rank": 1, "nspin": 1, "xc": "vbh"},
  "status": {"converged": true, "iterations": 19, "warnings": []},
  "fermi_energy_ry": -0.03543366,
  "total_energy_ry_per_site": -3329.13456789,
  "energy_terms_ry_per_site": {"band": -0.0, "core": -0.0, "double_counting": -0.0,
                               "hartree": 0.0, "xc": -0.0, "madelung": -0.0, "screened": 0.0},
  "components": [
    {"label": "Cu:A", "species": "Cu", "sublattice": "A", "z": 29, "concentration": 0.5,
     "weight_per_site": 0.5, "qnet_e": -0.01235, "moment_mub": 0.0,
     "occupation_e": {"s": 0.68472, "p": 0.66416, "d": 9.5098},
     "dos_ef_per_ry": {"s": 0.3576, "p": 0.477, "d": 3.7234},
     "energy_terms_ry": {"band": 0.0}}
  ],
  "dos_ef_per_site_per_ry": 3.3002,
  "history": {"etot_ry": [], "dv_ry": [], "ef_ry": [], "npass": [], "wall_s": []},
  "checks": {"nnoncausal": 0, "nnode": 4480, "wbelow_first_e": 2e-9, "wbelow_final_e": 1e-9,
             "nclamp_last": 0, "nclamp_run": 0, "qcell_e": -1.8e-5},
  "timing_s": {"total_wall": 187.52, "total_cpu": 3239.75, "parts": {}}
}
```

Numbers are written with `es24.16` (the full double).  The JSON is the source for paper tables; `tools/` scripts move to it.

---

## 6. Provenance at build time

fpm has no configure step.  Plan:

* `src/asr_build_info.f90` is committed with `git_rev = 'unknown'`, `git_dirty = .false.`.
* `tools/gen_build_info.sh` writes the file from `git describe --always --dirty` and the flags.  It rewrites the file only
  when the content changes, so fpm recompiles one module and its users.  It sets `git update-index --skip-worktree` on the
  file, so the generated version never shows as a change.
* `asr_build`, `asr_test` and `asr_run` in `env.sh` call the script first.  A plain `fpm build` still works and reports
  `unknown`.  The smriti build script calls the script too.
* Compiler and flags come from `compiler_version()` and `compiler_options()` (`iso_fortran_env`, Fortran 2008).  Host from
  `/proc/sys/kernel/hostname`, falling back to `HOSTNAME`.

---

## 7. Code structure

* `asr_report` stops using `asr_asa` types.  It gets one record, `rep_run_t`: labels, Z, valence, weights, per-component
  parameters and moments **already in both-spin electrons**, energy terms, history arrays, checks, timer and provenance.
* Two adapters fill it: `rep_from_class_scf(inp, res, rep)` in `asr_driver` and `rep_from_asa(inp, cell, res, rep)` in
  `asr_lmto_driver`.  The spin conversion (fs = 2 / nspin for the recursion fields) happens there and only there.
* Report routines, one per section: `rep_banner`, `rep_echo_input` (kept), `rep_method`, `rep_structure`, `rep_sro`,
  `rep_scf_history`, `rep_energy`, `rep_components`, `rep_fermi`, `rep_checks`, `rep_timing`, `rep_warnings`,
  `rep_close`; `rep_json` writes the summary from the same record.
* One helper writes a table in column blocks so that no table exceeds 90 columns for any number of components.
* Only rank 0 writes, as now.

---

## 8. DOS file (`*_dos.dat`)

The data columns stay as they are (plotting scripts read them).  Only the header changes: one line per column.

```text
# ASR density of states
# run        : CuZn bcc alloy   (example/cuzn_alloy.toml, git 2d998ac)
# route      : recursion edc, levels 4;  broadening eta 0.020000 Ry (G at E + i eta)
# E_F        : -0.03543366 Ry
# mesh       : 401 energies, -0.900000 .. 0.700000 Ry
# units      : E in Ry; DOS in states/Ry per atom of a component, or per site
# spin       : nspin 1; DOS columns are sums over both spins
# column  1  : E [Ry]
# column  2  : total per site, sum_p w_p N_p
# column  3  : Cu:A s
# column  4  : Cu:A p
# column  5  : Cu:A d
# column  6  : Zn:A s
# column  7  : Zn:A p
# column  8  : Zn:A d
# column  9  : total per site, spin 1
# column 10  : total per site, spin 2
```

`numpy.loadtxt` and gnuplot skip the `#` lines.  The k-space route writes the same header with `route kspace` and the
smearing.

---

## 9. Restart file and large data

* **Restart potentials** already check Z, mesh size and sphere radius on reading, and carry the SRO drive u.  No change
  in this plan.  Saving the Anderson history is mixer state, not output; it is out of scope.
* **Spectral function A(k,E)** (branch `spectral`, not merged): a k path of a few hundred points times a thousand energies
  times components and spins is about 10^6 to 10^7 numbers.  Write it as `.npy` (NumPy format: a short text header and
  raw little-endian doubles).  A Fortran stream writer is about 30 lines.  NumPy reads it directly and can memory-map
  it.  No library, no compiler-specific `.mod` files.  HDF5 waits until a dataset is too large for this.
  The `ef_scan_dump` chain coefficients stay a diagnostic text file.

---

## 10. Implementation and verification

1. `rep_run_t` and the two adapters; the section routines; the column-block helper.  The old routines go.
2. Provenance module and script; JSON writer; DOS header.
3. Fast unit test `test/test_report.f90` (laptop):
   * a synthetic `rep_run_t` with 1, 3 and 9 components, nspin 1 and 2;
   * every line of the report at most 90 columns; every section header present, in order;
   * the energy rows add up to the total line as printed;
   * gamma_e and chi_P of N = 3.5722 states/Ry are 0.6189 mJ/(mol K^2) and 8.488E-06 emu/mol;
   * per-spin conversion: a recursion record and a k-space record of the same state print the same q0.
4. Desktop runs (not on the laptop): `data/kspace/ag_fcc.toml`, `example/cuzn_b2.toml`, `example/cuzn_alloy.toml`,
   `example/fe66ni14cr20_fcc_sro.toml`, `example/co2fe04cr06al_l21.toml` (spin polarised).  For each run check: the
   numbers in the new report equal the old report and the JSON; the slow test `test_scf_guards` passes (energy terms,
   timers, DOS at E_F).
5. Merge after the desktop runs pass.  `output-data` (2d998ac) merges cleanly with `spectral` and `epi-sro`.

---

## 11. Not in this plan

* Anderson history in the restart file (mixer state).
* Mixing or formation energies (need reference runs).
* Electron-phonon or Stoner enhancement of gamma_e and chi_P.
* HDF5.
