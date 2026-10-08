# asr — augmented-space recursion for disordered alloys

A from-scratch modern-Fortran code for electronic structure of random alloys beyond CPA: high-entropy alloys,
disordered altermagnets, short-range order. The LMTO part is built in (own radial solver, libxc, own structure
constants, own SCF); other Hamiltonians (Wannier, tight binding) plug in through the same operator interface.

Everything is validated against independent references (exact enumeration, k-space, NIST atoms, Stuttgart LMTO
output, KPM on large supercells). `docs/STATUS.md` lists only verified claims; `docs/TODO.md` and `docs/ROADMAP.md`
are kept current.

## Build and test
```
source env.sh                 # Intel oneAPI (ifx, MKL) + fpm flags
asr_build --tests             # fpm build
asr_test                      # all tests (≈ 3 min on 16 cores)
```
Needs: ifx, MKL, fpm, libxc ≥ 7 (C library), toml-f (fetched by fpm; cached under `build/dependencies`).
After a module *signature* change do `rm -rf build/ifx_*` — fpm's incremental build can keep stale objects.

## Run
```
asr example/ag_fcc.toml       # self-consistent fcc Ag: SCF → parameters → l-resolved DOS (5 s)
asr example/fe_bcc.toml       # ferromagnetic bcc Fe
asr example/fecu_b2.toml      # ordered B2 FeCu (two classes, Madelung, charge transfer)
```
Outputs are written next to the input file. Input is one validated TOML record; unknown keys are errors.
`[calculation] kind = "lmto"`: `[structure]`, `[[species]]`, `[scf]`, `[dos]`, `[output]` (see `src/asr_input_lmto.f90`).
`kind = "model"`: ordered tight-binding DOS via recursion (`example/*_ordered.toml`).

## What is inside (src/)
| layer | modules |
|---|---|
| base | `asr_kinds`, `asr_operator` (abstract `linop`, CSR), `asr_testing`, `asr_tomlutil`, `asr_input`, `asr_input_lmto` |
| recursion | `asr_recursion` (Lanczos/Haydock), `asr_greens` (continued fraction, tails), `asr_terminator` (sqrt, Beer–Pettifor, Luchini–Nex join, period-p), `asr_contour`, `asr_partial` (species-resolved G) |
| LMTO | `asr_radial` (scalar-relativistic Koelling–Harmon, exact E-derivatives), `asr_xc` (libxc), `asr_atom`, `asr_potpar`, `asr_harmonics`, `asr_strconst` (screened S^α), `asr_lmto_ham`, `asr_lmto_cell`, `asr_ewald`, `asr_tetra`, `asr_mixing` (Anderson), `asr_asa` (ordered SCF, spin, multi-class), `asr_lmto_driver` |
| disorder | `asr_augspace` (binary augmented space, hashed states), `asr_symmetry` (O_h reduction, ~40×), `asr_augorb` (multi-orbital H = E + h − hoh), `asr_lmto_alloy`, `asr_lmto_realization`, `asr_alloy_scf`, `asr_kpm`, `asr_lattices`, `asr_kspace` |

## Validation highlights
* radial solver: hydrogenic levels 1e-10 Ry; free atoms He/Ne vs NIST LDA to 3e-7 Ha
* own ASA SCF reproduces Stuttgart Ag parameters (√Δ, γ, p) to ~1%; Fe 2.2 μB, Ni 0.64 μB, Cu 0; B2 FeCu occupations to 0.01–0.02 e
* augmented-space recursion equals the exact average over all 2^N configurations (1e-14), also for the multi-orbital TB-LMTO form
* real-space augmented-space recursion equals the ordered k-space recursion coefficients (1e-14)
* alloy SCF in the trivial-alloy limit reproduces the ordered SCF

## Conventions
Rydberg atomic units; `alat` in Å in inputs, Bohr internally; orbital order s, p(y,z,x), d(xy,yz,3z²−r²,xz,x²−y²).
Stuttgart (AMscf) and Questaal are reference *code* only — no numbers from them are used as inputs.

## Machines and code sync
Heavy runs go to the desktop (`ssh desktop`, needs the SRM VPN: `nmcli con up id "SRM"`).  Code moves only through git:
the bare hub `desktop:~/asr.git` (remote `hub`) is shared by this checkout and the desktop working copy `~/asr_remote4`.
Commit, then `asr_sync` (laptop -> hub -> desktop fast-forward) or `asr_pull` (fetch commits made on the other side).
The desktop working copy must stay clean (commit there before syncing); run outputs live in the ignored `trial/`.

