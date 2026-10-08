# STATUS — verified claims only

Build: `source env.sh; asr_test` (fpm + ifx 2026.1 + MKL + OpenMP).
Last updated 2026-10-02.

## M0 skeleton — done
fpm project, `env.sh` helpers, in-repo test harness (`asr_testing`).

## M1 recursion core — partly done
| claim | test | status |
|---|---|---|
| Lanczos/Haydock recursion with exhausted Krylov space == dense diagonalisation (MKL dsyev), G(z) to 1e-14, 2nd moment to 1e-15 | `test_recursion_dense` | pass |
| ordered sc/bcc/fcc NN band: a1 = 0, b1^2 = z (6, 8, 12) | `test_ordered_dos` | pass |
| sc, fcc LDOS vs k-space (nk=128, eta=0.1): error falls with depth, 2.6e-4 / 3.9e-4 at depth 40 with sqrt tail | `test_ordered_dos` | pass |
| `asr <sc\|bcc\|fcc>` writes DOS, integral 0.995-0.998 (eta=0.05, window ±1.2 z) | manual | ok |

## Input layer — done (model kind)
`asr <input.toml>` from any directory; one validated TOML record
(asr2 section layout: `[calculation] [structure] [model] [recursion] [dos]
[output]`), output written next to the input file. Unknown sections/keys and
bad values are errors (`test_input`), never silent defaults. Examples in
`example/`. Sections for `[[species]]`, `[scf]`, `[augmented_space]`, `[sro]`
will be added as the corresponding milestones land.

## M4 start: ordered TB-LMTO-ASA (fcc Ag, Stuttgart converged parameters) — verified
Own code (`asr_harmonics`, `asr_strconst`, `asr_lmto_ham`), nothing taken from asr2.
| claim | test | status |
|---|---|---|
| real Y_lm orthonormal, addition theorem; Gaunt coefficients by exact quadrature | `test_harmonics` | pass |
| bare S0_ss = -2w/d, S0_{LL'}(d)=S0_{L'L}(-d) | `test_strconst` | pass |
| screened (sigma=0.7) fcc on-site 9x9 block == legacy INFO_FCC (s 3.2326, p 3.0223, t2g 2.7422, eg 1.9492) to 6e-4; fixes the orbital order and the alpha_l = sigma^(2l+1)/(2(2l+1)) set | `test_strconst` | pass |
| H(k) Hermitian; cubic degeneracies at Gamma (p3, t2g3, eg2) | `test_lmto_ham` | pass |
| DOS vs Stuttgart LMASA Ag (ctrl parameters in `data/reference/ag_fcc_stuttgart`): d weight 9.58 vs 9.50, d centroid within 0.6 mRy, occupations s/p/d 0.687/0.633/9.635 vs 0.683/0.648/9.565, N(E_F) 3.165 vs 3.177 /Ry/atom | `test_lmto_ham`, `probe_ag_dos` | pass |

Hard-won conventions (each was a bug once): (1) the Stuttgart ctrl stores c, sqrdel, p in the
**gamma representation** and they must be transformed to the structure-constant alpha with
`pptrns` before use (skipping it put the p band at +3 Ry); (2) S enters with sign +; sign of
sqrdel(p) is a pure gauge; (3) the finite-cluster S^alpha row is Hermitian only to truncation
error — the Bloch sum is Hermitianised.

**Known, understood-not-fixed difference:** above ~+0.4 Ry (high s,p) our H^(2) departs from
Stuttgart (p band top 1.27 vs 1.58 Ry; p weight too high). The Stuttgart run has the combined
correction (CCOR=T) and downfolded f; ours has neither. Irrelevant near E_F, relevant for
anything that needs the full sp band. Not yet tested with CCOR off (that needs a Stuttgart run).

## Own LMTO atom stack (no Stuttgart numbers as inputs) — verified
`asr_radial` (scalar-relativistic Koelling-Harmon on a shifted log mesh, RK4, analytic E-derivatives to 2nd
order, bound states by node-count + outward/inward matching at the turning point), `asr_xc` (libxc C ABI),
`asr_atom` (Poisson + free-atom SCF), `asr_potpar` (C, sqrt(Delta), p, gamma from phi, phidot at S).
| claim | test | status |
|---|---|---|
| hydrogenic levels E=-Z^2/n^2 (nonrel) to 2e-8 Ry (Z=10), 2e-10 (Z=1) | `test_radial` | pass |
| scalar-relativistic shift = first-order Pauli (mass-velocity + Darwin) to O((Z alpha)^2) | `test_radial` | pass |
| dP/dE, d2P/dE2 from variational equations == finite differences (5e-8, 8e-8) | `test_radial` | pass |
| free atom LDA (VWN): He E=-2.8348356 Ha (NIST -2.834836), eps_1s -0.570425; Ne E=-128.233481 Ha (NIST), eps -30.305855/-1.322809/-0.498034 | `test_atom` | pass, 3e-7 Ha |
| potential parameter C == exact zero of W{K,phi}; error collapses as E_nu -> C | `test_potpar` | pass |
Not yet verified: sqrt(Delta), gamma, p normalisations (taken from the Stuttgart definitions; confirmed
only when the SCF reproduces the Ag parameters).
Hard-won: (a) mesh derivative is dr/ds = a(r+b) per index step — RK4/Simpson step is 1, not a; (b) the
outward solution blows up past the turning point, so densities need outward-to-match + inward; (c) form the
matching function from per-solution-normalised (P,Q) — squaring 1e167 overflows.

## Self-consistent ASA fcc Ag, parameters generated entirely by our own code — verified
`asr_asa` (+ `asr_lmto_cell`, `asr_potpar`, `asr_radial`, `asr_atom`, `asr_xc`): scalar-relativistic radial
solver, relaxed core in the sphere (Hankel tail at S), neutral-sphere Hartree, libxc von Barth-Hedin, first-order
ASA bands (generalised problem), moments -> charge (rho = q0 phi^2 + 2 q1 phi phidot + q2 phidot^2, gamma rep),
E_nu to band centre, potential mixing.  Input: Z, valence count, lattice constant, valence configuration only.
Starting potential: superposed free-atom (a screened-Coulomb start converges to a *wrong* delocalised-d fixed
point: E_nu,d=+0.36 Ry — keep the atomic start).
| quantity (ours vs Stuttgart oracle) | s | p | d |
|---|---|---|---|
| sqrt(Delta) | 0.36493 / 0.36457 | 0.36540 / 0.36518 | 0.10446 / 0.10341 |
| gamma | 0.42493 / 0.42478 | 0.11544 / 0.11538 | 0.00519 / 0.00550 |
| p | 0.08173 / 0.08230 | 0.03609 / 0.03623 | 2.0604 / 2.0814 |
E_F - V(S): 0.5080 vs 0.5139 Ry; V(S) -0.6846 vs -0.7003; E_nu,d - E_F: 0.308 vs 0.327.
`test_asa_ag`: 13 assertions pass. This also **confirms** the sqrt(Delta), gamma, p normalisations in `asr_potpar`.
Not explained: 10-30 mRy absolute offsets of E_nu/C/V(S) (Stuttgart xcloc1 is the same VBH functional; numerics
converged: nk and smearing change E_F by <1 mRy; PW instead of VBH moves everything ~37 mRy). Candidates: snapshot
not fully converged, Hartree reference convention, core treatment details. Not tested: spin polarisation,
multi-atom cells (needs Madelung), tetrahedron integration.

## Collinear magnetism and multi-class SCF — verified
`asr_asa` now handles any number of atoms/classes and nspin = 1 or 2 (libxc spin-polarised, per-spin core and
E_nu, common Fermi level). Gates are literature LDA/ASA values (no files): `test_asa_mag` (9 assertions).
| system | moment (muB) | gate |
|---|---|---|
| bcc Fe (a=2.8665 A) | 2.181 (nk=12, 10 mRy) -> 2.233 (nk=24, 3 mRy) | 2.15-2.30 |
| fcc Ni (a=3.524 A) | 0.639 | 0.58-0.68 |
| fcc Cu, seeded zeta0=0.4 | 0.0000 | |m|<0.002 |
`test_xc`: polarised == unpolarised at zeta=0; majority potential deeper. Bug fixed on the way: out-of-bounds write of
the second spin column in a non-magnetic run (heap corruption, SIGSEGV at exit) — assumed-shape + -check bounds did not flag it.

## Multi-class SCF with Madelung (B2 FeCu from scratch) — verified
New: `asr_ewald` (point-charge Ewald with neutralising background; CsCl 1.762675 and NaCl 1.747565 to 2-4e-7,
cell-independent), `asr_mixing` (Anderson/Pulay mixing; converges a stiff map plain linear mixing cannot),
multi-class `asa_scf` with Madelung shifts from net sphere charges. `test_asa_b2` (9 assertions), `test_ewald` (4),
`test_mixing` (3).
Two instabilities found and fixed: (1) charge sloshing — linear mixing is unstable when beta x DOS x Madelung
coupling > 2 (initial transfer ~3 e, self-potential ~2 Ry); (2) mixing the raw potential is wrong because
V = -2Z/r + ... is ~1e6 Ry near the nucleus, so clamps/interpolation errors there never relax — mix the smooth part
u = V + 2Z/r. With both fixed Ag converges in 14 iterations (was 37) and FeCu in 17.
B2 FeCu (nonmagnetic, a=3.0 A, Fe:Cu radius ratio from the oracle geometry): occupations Fe 0.657/0.789/6.744, Cu
0.707/0.687/9.417 vs Stuttgart 0.651/0.777/6.763, 0.706/0.677/9.426; E_F -0.1935 vs -0.1871 Ry; charge transfer 0.19 e to Fe.

## BZ integration — linear tetrahedra (verified)
`asr_tetra`: Bloechl-Jepsen-Andersen weights on a periodic regular grid in fractional coordinates (any lattice),
closed-form occupied volume for the Fermi level (tetrahedra sorted once). `test_tetra` (7): free-electron sphere
converges at second order (error x4 per grid doubling, 1.4% at nk=32 for a sphere only ~4 spacings wide), exact
half filling of the sc band, sum rules, closed-form count == weight sum. SCF default is now `bz = "tetra"`.
Gaussian smearing biased N(E_F) up on the Fe d peak; tetrahedra remove that bias: Fe 2.204 -> 2.214 muB (nk 12 -> 16),
Ni 0.603 -> 0.634 (nk 12 -> 28) approaching the small-sigma Gaussian 0.640. B2 FeCu unchanged (E_F -0.1932).
Bug fixed: `bound_state` no longer aborts when a transient SCF potential has an unbound core level (returns best estimate + flag).

## M2 augmented space (binary, single orbital) — verified; first disordered terminator study
`asr_augspace` (hashed configuration states, BFS generator of the reachable (site, flipped-set) states to a given level,
exactly symmetric truncated operator), `asr_lattices` (periodic sc/bcc/fcc graphs, random realisations), `asr_kpm`
(Chebyshev/Jackson DOS, stochastic trace; independent truth, `test_kpm`).
| claim | test | status |
|---|---|---|
| all N 2^N augmented states reached on a 2x2x2 cluster; operator exactly symmetric | `test_augspace` | pass |
| augmented-space recursion == exact average over ALL 256 configurations (diagonal + off-diagonal disorder, x=0.3), 1.8e-14 | `test_augspace` | pass |
| KPM vs exact spectrum of the same supercell, 1.65% with 60 vectors (falls as 1/sqrt(R)) | `test_kpm` | pass |
Disordered sc alloy (x=0.5, eps_A-eps_B = 3t, diagonal disorder), eta = 0.15, KPM truth on 131 072 sites (noise 0.9%),
max|dDOS|/peak at depth 2/3/4/5/6 (augmented states 230 / 1.8e3 / 1.4e4 / 1.1e5 / 8.1e5):
  sqrt-last 11.4/6.0/2.4/2.6/2.1 %; sqrt-avg2 10.3/8.6/4.3/2.2/2.4 %; Beer-Pettifor 156/50/15/8.0/6.0 %;
  BP + Luchini-Nex join 102/27/11/4.4/1.9 %; period-2 234/-/47/-/11.5 %.
Reading: (i) depth limits accuracy, ~2% at depth 6 for strong disorder; (ii) the LN join is what makes Beer-Pettifor usable
at short depth; (iii) the period-2 tail is wrong for this alloy; (iv) augmented states grow ~7.5x per level, so depth 8
needs ~45M states without symmetry -> M3 (point-group reduction) is the lever. (`test/probe_asr_depth.f90`)

## Disordered terminator study to depth 12 (run on `desktop`, 28 cores, 8 min)
sc, x=0.5, eps_A-eps_B = 3t, diagonal disorder, eta=0.15, point-group reduced augmented space, KPM truth 786k sites (noise 0.9%).
max|dDOS|/peak by depth 2,4,6,8,10,11,12 (reduced states 7 / 47 / 375 / 3.6e3 / 4.0e4 / 1.4e5 / 4.7e5):
  sqrt-last  11.4 / 2.4 / 1.9 / 1.3 / 1.0 / 0.63 / 0.43 %     sqrt-avg2 10.3 / 4.2 / 2.2 / 1.4 / 1.1 / 0.82 / 0.52 %
  BP         156 / 15 / 5.8 / 2.5 / 1.1 / 1.6 / 0.74 %         BP+LN     102 / 11 / 1.9 / 1.1 / 1.0 / 0.72 / 0.56 %
  period-2 234 / 47 / 11.5 / 7.9 / 6.0 / - / 4.4 %
Reading: ASR converges to the truth (at its noise floor by depth ~10); with only diagonal disorder the state count is small
(zero-coupling states are never created); a constant sqrt tail is as good as anything at depth >= 10; BP needs the LN join
below depth ~8; the period-2 tail is never competitive for this alloy.
Off-diagonal disorder (t_AA=-1.3, t_AB=-1, t_BB=-0.7, delta=2), truth noise 0.27%, rank <= 2: sqrt-last 16.5 / 5.2 / 4.6 / 3.0 / 2.4 %
at depth 2/4/6/8/9 (states 22 ... 3.2e4) — slower than the diagonal case (ranks 3, 4, unlimited still running).

## M6 alloy SCF from the recursion — ordered-limit gate verified (`test_alloy_scf`, 14 assertions)
New: `asr_partial` (species-resolved Green functions, transition element; exact vs enumeration to 2e-15, `test_partial`),
`asr_contour` (semicircle contour moments; M0 within 3-5e-4 of tetrahedra, exactly 0.5 at half filling, `test_contour`),
`asr_alloy_scf` (four recursions s,p_z,t2g,eg with stored Lanczos bases, species-resolved moments, E_F from the electron
count, orthogonal-representation density rho_l = M0 phi^2 + 2 M1 phi phidot + M2 (phidot^2 + phi phiddot), relaxed core,
Hartree + xc + SIM Madelung, Anderson mixing), `augorb` flip-radius option (flips only within a graph distance of the start; the
medium outside is the virtual average) alongside the rank cutoff, `asr_lmto_realization` (numerical realisation; truth operator).
Gate: A = B = Ag (no disorder), 8 recursion levels, fcc: converges in 13 iterations to the ordered first-order-ASA SCF:
sqrt(Delta) s/p/d 0.36494/0.36515/0.10455 vs 0.36491/0.36517/0.10442; gamma 0.42481/0.11532/0.00523 vs 0.42492/0.11534/0.00520;
p 0.0820/0.0363/2.055 vs 0.0818/0.0363/2.061; E_nu,d within 1.2 mRy, C_d-E_nu within 0.8 mRy; E_F -0.1603 vs -0.1737 (the
13 mRy is the orthogonal vs first-order representation). Not yet: a genuinely disordered alloy (CuNi etc.), spin polarisation,
convergence of the ASR approximation (levels / rank / flip radius) for real LMTO alloys (flip_scan running on desktop).

## First disordered alloy SCF (Cu0.5 Ni0.5, fcc, a=3.57 A, non-magnetic) and the flip-radius result (2026-10-03, desktop)
`probe_cuni_scf 5 20 2 1`: 5 levels, rank <= 2, flip radius 1, 650 297 nodes; ~40 iterations (final dV 2.3e-5), E_F = -0.1362 Ry,
charge transfer 0.012 e from Cu to Ni, species parameters differ as expected (d: Cu enu -0.317, Ni -0.273 Ry). A working end-to-end
disordered-alloy SCF; at this truncation it is a demonstration, not a quantitative result.
**Flip-radius truncation is a poor control** (sc model, x=0.5, delta=2t, d1=0.3, rank <= 4, KPM truth noise 0.6%): the error plateaus at
~11-15% for flip radius 1..4 (rf=1: 11.0% at depth 9; rf=2: 11.1% at depth 6; rf=3: 14.6% at depth 5; rf=4: 13.7% at depth 4), whereas exact
generation (no flip radius) with rank <= 2 reached 2.4% at depth 9. Outside the flip radius the medium is the virtual crystal and the
walker spends most of its time there at depth. Conclusion: use exact breadth-first generation with a rank cutoff and make it affordable
with symmetry reduction (the orbital-aware point-group reduction is now the next engineering step).

## LMTO alloy accuracy of the integrated moments (Cu0.5Ni0.5, exact generation, rank <= 2; desktop 2026-10-03)
The right benchmark for the SCF is the integrated moment, not the fine DOS. s-orbital moments up to E_F = -0.10 Ry about E_nu = -0.5 vs KPM
truth (16 000 sites x 4 realisations): depth 2 (97 357 nodes): M0 -11.5%, M2 +146%; depth 3 (2 046 337 nodes): M0 **-1.1%**, M2 **-0.4%**
(M1 ~ -0.008 is a cancellation, abs. error 0.002). DOS at eta = 0.15 Ry: 37% -> 12% (sqrt-avg2). A DOS comparison at eta = 0.02 Ry is not
meaningful at 2-3 levels. Node growth ~21x per level without symmetry (depth 4 ~ 4e7), so the orbital-aware point-group reduction
is how to go deeper, but depth 3 already gives ~1% occupations.

## Cu-Ni SCF at depth 3 (exact generation, rank <= 2) and the limit for magnetic alloys (desktop 2026-10-03)
Cu0.5Ni0.5 non-magnetic: converged in 29 iterations, 2 046 337 nodes, E_F = -0.1357 Ry, charge transfer 0.047 e Cu -> Ni (0.012 e at 5 levels
with flip radius 1, so the transfer is not depth-converged). Cu0.2Ni0.8 spin-polarised (zeta0 = 0.1): converged in 43 iterations but the
moment decays to ~0. **Pure Ni (A = B = Ni, no disorder) through the alloy route at depth 3 also loses its moment** (starts 0.5 muB, ends
1e-4; ordered ASA route: 0.64), so the zero moment is a recursion-depth limitation (3 levels cannot resolve the d-band DOS at E_F that a
ferromagnetic instability needs), not a disorder effect. Fe (E_F in the d band too) was fine at 8-16 levels. zeta0 = 0.3 for Cu0.2Ni0.8 runs
away (core splitting, +-3 e charge transfer), same instability as Fe with a large seed.
Consequence: magnetic Ni-based alloys need more levels than exact generation affords without symmetry (~2e6 nodes at depth 3, ~4e7 at depth 4).
The orbital-aware point-group reduction is now the gating item for magnetic disordered alloys and for HEA-scale depth.

## Spin-polarised alloy route: Fe gate — instability found and understood (2026-10-02)
Trivial alloy A=B=bcc Fe through the recursion route. Non-magnetic: converges (32 it, E_F -0.0943). Spin-polarised with initial
polarisation zeta0 = 0.4 (applied to the whole atomic density incl. the core, splitting core levels by ~13 Ry at the nucleus and
starting at ~3 muB) **runs away** (moment -> 8, then wanders; E_F -> +1.5 Ry) for every mixing strength (beta 0.3 / 0.1) and
contour resolution (60 / 200 points) tried; zeta0 = 0.1 converges in 35 iterations to |m| = 2.36 muB (ordered ASA route 2.22; the
6% excess is not yet explained: orthogonal-vs-first-order representation or ASR depth). Default zeta0 is now 0.1.
Also fixed: the alloy-SCF convergence measure was unweighted (dominated by r ~ 1e-5 Bohr at the nucleus); now weighted by
r/(1+r) like the ordered code. Sum rule (E_F far above the bands: M0 = 1, 3, 5 for s, p, d) verified exact for Fe at 8 and 14 levels.
Fe moment vs recursion depth with the gentle seed (desktop): 8 / 10 / 12 / 16 levels -> 2.262 / 2.362 / 2.344 / 2.351 muB, E_F -0.082...-0.093 Ry,
all converged in 29-35 iterations. Depth-converged at ~2.35 muB, so the ~5-6% excess over the ordered first-order-ASA route
(2.22; Gaussian-limit 2.23) is a representation effect (orthogonal second-order moments vs generalised first-order problem), not a
recursion-depth effect — still to be explained. `test_alloy_scf_spin` (5 assertions) passes with tolerance 0.2 muB.
TODO: a seed that polarises only the valence region.

## M5 core: multi-orbital TB-LMTO-form augmented space (verified)
`asr_augorb`: H = E_nu + h - h o h with h = c + d S d and per-site, per-orbital binary-random c, d=sqrt(Delta), o
(X = X0 + X1 n; n -> M, a projector). Matrix-free application as a sequence of the elementary augmented operators (random
diagonal D, hop S via a bond-block table); node set = closure of the support of H (paths R->R1->R2, any flip subset of
{R,R1,R2}), which contains every intermediate vector; optional rank truncation. `test_augorb`: ring of 6 sites, 2 orbitals,
all 64 configurations enumerated: all N 2^N nodes reached, operator symmetric (2e-13), recursion == exact configuration
average to 4e-15. `test_alloy_ordered`: with our own Ag SCF parameters, fcc, 18 bond slots, the real-space augmented-space recursion
reproduces the ordered k-space (`hk_grid`) recursion coefficients to 2.5e-14 (s) and 1.6e-14 (e_g) over 6 levels — validates
the graph, bond-block table (incl. the on-site S^alpha block), parameter mapping and H^(2) apply. `test_lmto_real`: realisation operator.
Not yet: real structure constants (multi-shell graph + bond blocks), orbital/point-group reduction,
allocation-free apply, rank-truncation accuracy study (running on `desktop`).

## Ordered B2 FeCu (bcc lattice, two sublattices) from the Stuttgart `AMscf/lmto/fe` parameters — verified
General cell code (`asr_lmto_cell`: any lattice/basis, per-species parameters, per-sphere hard-core alpha).
Note: `AMscf/lmto/fe` is **ordered FeCu in a CsCl cell**, not elemental bcc Fe. Nonmagnetic, spd, no downfolding.
| claim | test | status |
|---|---|---|
| cell code == single-atom Ag code (spectrum at generic k, 3e-15) | `test_lmto_cell` | pass |
| screened-S cluster converged: occupations move <0.003 from 27 to 137 sites | `test_lmto_cell` (earlier run) | ok |
| N(E_F=-0.18714 Ry) = 19.11 (nk=64, sigma=1.5 mRy) vs 19 electrons | `test_lmto_cell` | pass, within 0.2 |
| d occupations: Fe 6.766 / Cu 9.429 vs Stuttgart moments 6.763 / 9.426 (nk=64); 6.79 / 9.43 at the test's nk=48 | `test_lmto_cell` | pass |
| Fe d DOS vs Stuttgart `fe/data.dos` (identified as Fe d): peak 173.5 vs 176.2 at the same energy -0.1817 Ry, correlation 0.95 | manual | ok |

Open: s,p occupations run 3-5% above the Stuttgart moments (Fe s 0.678 vs 0.651, p 0.798 vs 0.777;
Cu s 0.740 vs 0.706, p 0.698 vs 0.677) = the whole 0.1 e surplus. Suspected cause: Stuttgart
moments include the o and p (phi-dot) terms, ours are plain orthonormal |u|^2. Untested.
Numerical lesson: E_F sits on the sharp Fe d peak, so Gaussian smearing biases N(E_F) up
(19.35 at 6 mRy -> 19.11 at 1.5 mRy); use a tetrahedron/Methfessel-Paxton integrator or tiny sigma.

## Recursion on the real LMTO Hamiltonian (ordered fcc Ag) — measured
`kgrid_op` (H(k) on an nk^3 grid, re/im doubled) + Lanczos from |0,L>, compared with the exact
orbital Green function on the *same* grid, so the only difference is the tail
(`test/probe_lmto_term.f90`, output `trial/ag_term_bench_nk24_eta0.02.txt`; nk=24, eta=0.02 Ry,
window -0.9..0.4 Ry; error = max|dLDOS|/peak LDOS):
- depth 10 / 20 / 40 / 60:  s ~1.8e-1 / 9e-2 / 2.2e-2 / 7.5e-3;  t2g, eg ~1.6e-1 / 5e-2 / 1e-2 / 2.5e-3.
- **No terminator rescues shallow depth.** All one-band tails agree within ~2x at every depth
  (Beer-Pettifor with an abrupt join is best at depth 10 for eg: 1.05e-1 vs 1.6e-1; the
  Luchini-Nex join erases that gain). The period-2 tail is *worse* here (Ag coefficients are not
  period-2). So on a real d band, depth is the lever: ~1% LDOS error needs depth ~40.
- Consequence for ASR: the augmented space reaches only ~10-20 levels, where the ordered-case error
  is already 5-18% whatever the tail. Disordered-case terminator behaviour is still unmeasured (M2).

## Terminator workstream T — first results (`app/term_bench`, `test_terminators`)
Implemented behind one `terminator_t`: asymptote = last | avg2 | Beer-Pettifor
(own derivation: denominator vanishes at both band edges; BP on a uniform chain
recovers (0,1), error ~1/n^2) | periodic-2; join = butt | Luchini-Nex.
Luchini-Nex is implemented from the published *description* only (linear
interpolation of coefficients from the last computed level to the asymptote);
the original paper has not been read — **unverified**.

Measured against exact k-space truth (max |dLDOS|, eta=0.1, depths 8/16/32/60):
- sc, fcc: all single-band terminators agree within a few % once depth >= 16;
  depth, not the tail, limits accuracy. Differences only at depth 8
  (BP + butt join worst, 1e-2 vs 7e-3 on sc; LN join repairs it).
- bcc: nothing helps; ~5e-3 at depth 60 for every one-band tail (van Hove).
- **B2 ordered compound (gapped, period-2 coefficients), delta=3t:**
  one-band tails 5.7e-1/1.9e-1/6.2e-2/1.9e-2 vs period-2 tail
  1.1e-2/3.0e-3/1.4e-4/2.7e-5 — 1 to 3 orders of magnitude better.
  This is the case that matters for alloys (gaps, d-band + s-p, sublattices).

Not yet shown: that period-2 (or any) tail helps in the *disordered* case, where
coefficients converge differently. That needs M2.

## Open problems (do not quote around them)
- **bcc terminator (single-band tails):** LDOS error vs k-space is 1.4e-2 at depth 40, falling only
  slowly with depth. The k-space reference is converged (nk=128 vs 256 differ
  ~1e-3), so this is the sqrt tail, not the recursion: bcc NN has a log van Hove
  singularity at band centre and b_n^2 oscillate even/odd (≈16 ± 0.2 at depth
  20). Averaging the last two levels does not fix it. First benchmark for the
  terminator workstream T.
- Recursion is real-symmetric, scalar only; block/complex-Hermitian versions
  needed for degenerate orbitals and non-collinear magnetism.

## Orbital-aware D4h reduction of the multi-orbital augmented space (verified)
- `asr_augsym` (16 D4h ops, canonical nodes, per-link (rep, op), normalised coordinates) reproduces the unreduced
  `asr_augorb` Lanczos coefficients to ~1e-10 for the s, p_z, t2g(xy), e_g(3z^2-r^2) starts (`test_augsym_orb`,
  fcc, depth 2, rank<=2, random c/d/o/E_nu).  Node reduction 13x at depth 2, ~68x at depth 3 (28,211 vs 1,918,600).
- Lesson: `d4h_ops(coord, m, ...)` takes m = 2*l for supercell of l cells (doubled coordinates); wrong m = silent
  out-of-bounds (found with -check bounds; -check uninit/MSan gives unrelated false positives).
- CuZn bcc (a=2.85 A, rank<=2, desktop): depth 2 E_F=0.0095 Ry, qnet=0.188; depth 3 E_F=0.0176 Ry, qnet=0.218
  (949k nodes): depth 2 not converged. DOS comparison with Tarafder 2006 still to do.

## Alloy driver, SIM sign, and the depth problem (2026-10-02)
- `asr_augsym` is the default operator in `alloy_scf` (`symmetry = true`); ordered-limit gates unchanged
  (`test_alloy_scf` 14/14, `test_alloy_scf_spin` 5/5: Fe 2.36 muB alloy route vs 2.22 ordered).
- TOML `kind = "alloy"` driver (`asr_alloy_driver`, `example/cuzn_alloy.toml`, `test_alloy_input` 7/7): SCF + species/l DOS in one
  input; run end-to-end on the desktop (CuZn bcc, depth 3, 63k nodes, 32 SCF iterations).
- BUG FIXED: SIM shift had the wrong sign (V = -beta 2 alpha_M qnet/S with qnet = N - Z gave positive feedback: charge transfer grew
  with beta, 0.05 -> 0.19 e).  Now V = +beta 2 alpha_M qnet/S (screening): CuZn qnet 0.024 e at depth 3.  NOT yet validated against a
  supercell reference; the ordered B2 ASA (own Ewald) gives qnet 0.135 and C_d(Cu)-C_d(Zn) = 0.32 Ry.
- Finding (ordered Ni, no disorder, alloy route): the Ni moment needs many Lanczos levels:
  depth 6/8/12 -> 0.0, depth 20 -> 0.53, depth 30 -> 0.58 muB (ordered ASA ~0.6).  Magnetic Ni-based disordered alloys are therefore
  limited by terminator quality / depth, not by the symmetry reduction.
- CuZn vs Tarafder 2006 (digitised): at depth 2-3 the Cu d and Zn d peaks merge (separation of C_d 0.07 Ry vs 0.32 Ry in the ordered B2 ASA
  and ~0.37 Ry in the published alloy DOS) -> two well separated d bands need more levels than 2-3.  Depth 4/5 runs on the desktop.
- Process: heavy runs go to the desktop (user instruction); probes run with `OMP_NUM_THREADS` set; never `pkill -f`.
- CuZn bcc depth 4 (desktop, 587k reduced nodes, rank<=2, SIM beta 0.7): converged in 25 iterations, qnet = +/-0.029 e,
  C_d(Cu) = -0.282, C_d(Zn) = -0.639 Ry (separation 0.357 Ry vs 0.32 in ordered B2 ASA, ~0.37 in Tarafder 2006): two separate d
  bands now appear.  Peaks relative to E_F: Zn d -0.68, Cu d -0.33 Ry vs published -0.49 / -0.12: both ~0.2 Ry too deep, and
  DOS(E_F) 1.5 vs 2.4 states/Ry/spin -> E_F (sp-band filling) is probably too high at this depth.  Depth 5 pending.

## Review of docs/ENERGY_DEPENDENT_ASR_AND_CUZN_FINDINGS.md (checked 2026-10-03)
Checked independently against code and data (not taken on trust):
- CONFIRMED by test: the contour arc stopped at Im z = eta (smooth Fermi step) made N(E_F) depend on eta (Ag, fixed potentials,
  8 levels: E_F crossing -0.167/-0.185/-0.215/-0.223 Ry for eta = 0.02/0.01/0.005/0.0025).  With the arc closed by the vertical
  segment E_F+i eta -> E_F (`contour_nodes_closed`) the crossing is eta-independent (-0.2257 for all eta), as Cauchy requires.
  NOT confirmed: the quoted mechanism "50% charge loss => +44 mRy": a Lorentzian step also adds weight above E_F; the net effect for a
  smooth DOS is small.  The eta-dependence here comes from the spiky model DOS.
- CONFIRMED: the staircase N(E) is a property of the finite-depth model DOS (real tail when 0 lies outside [a-2b, a+2b]); with the
  closed contour it shrinks with depth: crossing -0.226 (8 levels), -0.190 (12), -0.248 (16), -0.245 (24) at fixed atomic-start
  potentials, jump heights 0.07 -> 0.5 -> 0.004 -> 0.02.  Ordered crystals need >= 16 levels in this form.
- Legacy check (AMscf/ASR/main.f90, hop_endep.f90, newfit.f90): AMscf also uses recursion on E - D^(1/2)(P(E) - S)D^(1/2) (diagonal
  terms (E - C)/Delta, random parts ap3/ap4), but at only ienum+1 real seed energies, with cubic-spline interpolation of the
  coefficients to a real-axis mesh, a Beer-Pettifor terminator, and real-axis integration.  So our P(z)-S recursion is the same
  object; differences: complex contour vs real mesh, terminator, seeds+interpolation (cheap).  The review's statement that AMscf
  has a "branch cut because E is inside the sqrt" is true of any standard form; our zero-argument form is equivalent.
- CuZn Zn d position: Banerjee & Mookerjee 2010 Fig. 3 (digitised data peaks at -0.682 Ry; also checked by eye on the rendered PDF page: sharp Zn d peak at about -0.7 and a broad Cu d feature near -0.3 Ry), our depth-5 -0.684 Ry: agrees.
  Tarafder 2006 (-0.49) is a non-self-consistent figure; d-band splitting 0.37 Ry is common.  The earlier "Zn d 0.2 Ry too deep" was
  a wrong target (Tarafder), not a code error.

## Energy-dependent route status (2026-10-03, end of session)
- Ordered Ag, ED route, 8 levels: `test_alloy_scf` 18/18; E_F -0.1702 vs -0.1737 Ry (first-order ASA); params within ~1%.
  Ingredients that made it work: complex-symmetric Lanczos on P(z)-S (`asr_endep`), closed contour (arc + vertical segment),
  Beer-Pettifor-type tail from the chain's extreme eigenvalues (not last coefficients), secant E_F search on N(E).
- CuZn bcc ED depth 3 (rank 2, SIM beta 0.7, mix 0.3): converged (37 its), qnet +-0.048, E_F 0.045; Zn d -0.709, Cu d -0.175 vs E_F;
  C_d separation 0.495 Ry (energy-independent depth 5: ~0.38; Banerjee 2010 / Tarafder ~0.37): not yet matching; depth-4 pending.
- Unattended queue on the desktop: `~/asr_remote3/trial/queue.sh` (log `trial/queue.log`, per-run files and a one-line summary per run
  in `trial/queue/SUMMARY.txt`): ED depth 4 (ncont 32), ED depth 3 with ncont 32 (check vs 60), SIM beta 0.3 and 1.0, x = 0.25.

## Overnight queue results (read 2026-10-04): energy-dependent CuZn
| run (rank 2, mix 0.3) | conv. | E_F | Cu d | Zn d | C_d-sep | DOS(E_F)/spin |
|---|---|---|---|---|---|---|
| ED d3, 32-pt contour | yes (35 its) | 0.0455 | -0.174 | -0.709 | 0.535 | 1.72 |
| ED d3, 60-pt contour (earlier) | yes (37) | 0.0452 | -0.175 | -0.709 | 0.534 | 1.72 |
| ED d3, SIM beta 0.3 | dV 3e-5 at 50 its | 0.0467 | -0.177 | -0.705 | 0.528 | 1.68 |
| ED d3, SIM beta 1.0 | no (dV 0.5) | -0.197 | -0.339 | -0.703 | 0.364 | 3.25 |
| ED d3, x_Cu = 0.25 | no (dV 0.35) | 0.330 | -0.658 | -0.892 | 0.235 | 2.29 |
| ED d4, 32-pt contour | no (dV 0.025 at 50 its) | -0.016 | -0.116 | -0.663 | 0.547 | 2.33 |
| energy-independent d5 (reference run) | yes | 0.006 | -0.263 | -0.684 | 0.42 | 1.50 |
| Banerjee 2010 Fig. 3 (digitised) | | | -0.309 | -0.682 / -0.718 | ~0.37 | 2.47 (units uncertain) |
- 32-point contour == 60-point (same E_F, peaks): adequate.  SIM beta 0.3 vs 0.7 changes the separation by 0.007 Ry: electrostatics is
  NOT the cause of the ED separation.
- The ED route does NOT reproduce Banerjee for CuZn: Cu d too shallow by 0.14-0.19 Ry, separation ~0.54 vs ~0.37; depth 3 -> 4 does
  not improve it (0.535 -> 0.548).  The energy-independent route (depth 5) is within ~0.05 Ry on Cu d and 0.002 on Zn d.
  The ED SCF is also fragile (beta 1.0 and x = 0.25 do not converge).  Driver default therefore reverted to energy-independent.
- Next diagnostic: apples-to-apples DOS of both routes at FIXED (energy-independent, converged) potentials; if the DOS agree the
  difference sits in the ED moments/density path, if not in the Hamiltonian form.

## Fixed-potential comparison of the two routes (probe_route_dos, converged energy-independent potentials)
| depth | route | Cu d | Zn d | sep | electrons below E_F (tot / Cu d / Zn d) |
|---|---|---|---|---|---|
| 3 | EI | -0.468 | -0.496 | 0.03 | 11.18 / 9.15 / 9.76 |
| 3 | ED | -0.346 | -0.588 | 0.24 | 11.17 / 9.35 / 10.12 |
| 4 | EI | -0.334 | -0.682 | 0.35 | 11.14 / 9.19 / 9.59 |
| 4 | ED | -0.276 | -0.648 | 0.37 | 12.47 / 11.74 / 9.53 |
- Peak positions of the two routes converge towards each other (depth 4: within 0.06 Ry, separation 0.35 vs 0.37), BUT the ED species
  DOS is not normalised: at depth 4 the Cu d band integrates to 11.7 electrons (> 10) with a peak 3.4x higher than the EI one, and
  Zn d exceeds 10 at depth 3.
- Likely cause (hypothesis, not yet proven): the potential function P(z) = (z-C)/(Delta + (gamma-alpha)(z-C)) has a real pole at
  z_p = C + Delta/|gamma-alpha| (Zn d: about E_F - 0.26 Ry, i.e. inside the Cu d band; Cu d: about E_F + 0.2).  G = lambda + mu g mu
  has lambda, mu^2 ~ 1/den but the poles cancel only for the exact g; for a truncated recursion the cancellation is incomplete and a
  spurious spike appears at z_p (the spike at E_F - 0.276 matches the Zn z_p).  The formulation that avoids inverse potential
  functions is the second-order Hamiltonian re-linearised at the energy: H(E_nu = E) = E + h - h o h with parameters at E_nu = E
  (seed energies + interpolation, as AMscf does).

## Queue 2 results (read 2026-10-04): energy-independent CuZn checks, depth 4
| run | converged | E_F | Cu d | Zn d | sep | DOS(E_F)/spin | nodes |
|---|---|---|---|---|---|---|---|
| rank 1, x=0.5 | yes (25) | 0.0196 | -0.335 | -0.683 | 0.348 | 1.51 | 17,630 |
| rank 2, x=0.5 (earlier) | yes | 0.0196 | -0.335 | -0.683 | 0.348 | 1.51 | 587k |
| rank 3, x=0.5 | yes (25) | 0.0196 | -0.335 | -0.683 | 0.348 | 1.51 | 6.57M |
| rank 2, x_Cu=0.24 | yes (26) | 0.0922 | -0.454 | -0.750 | 0.297 | 1.28 | |
| rank 2, x_Cu=0.76 | yes (19) | -0.0480 | -0.258 | -0.601 | 0.343 | 1.94 | |
- Rank is converged: rank 1 = rank 2 = rank 3 to 3 digits at depth 4 (a 370x cheaper operator).  Use rank 1 for scans.
- Concentration trend DISAGREES with Banerjee 2010 Fig. 3: published Zn-d / Cu-d peaks get deeper with more Cu (Zn d -0.579, -0.682, -0.790 for
  x_Cu = 0.24, 0.50, 0.75; Cu d -0.224, -0.309, -0.387); ours get shallower (Zn d -0.750, -0.683, -0.601; Cu d -0.454, -0.335, -0.258).  Ours is
  the rigid-band direction (more Zn -> higher E_F -> d bands deeper below E_F).  Unresolved: check the published energy zero and the
  figure's panel labels before treating this as a code error.
- Fixed-potential route comparison at depth 5 (2.88M nodes): EI Cu d -0.260 / Zn d -0.684; ED -0.314 / -0.658; electrons below E_F
  EI 11.12, ED 12.81 -> the ED normalisation problem persists with depth.

## External-target proposal checked (2026-10-04) and first B2 FeCo (FM) result
Checked a pasted validation-target proposal against the repo and the cited sources:
- Repo numbers in the proposal are right (FeCu B2 E_F -0.1932 vs Stuttgart -0.1871; Fe 2.20-2.23 ordered / ~2.35 alloy route; Ni 0.53 / 0.58 at depth 20 / 30).
  Corrections: the CuZn Cu d -0.260 / Zn d -0.684 are the self-consistent energy-independent depth-5 peaks (not "fixed potentials");
  ordered Ni is 0.639 muB (tetrahedra, nk 28), not 0.60; FeCu B2 and its E_F are Stuttgart-oracle regression numbers, not an external published
  benchmark; the proposal omits the unresolved CuZn concentration-trend disagreement and the ED-route failure.
- Cited external numbers, verified in the sources: B2 FeCo DFT 2.80 muB (Fe) / 1.81 muB (Co) = VASP-PBE values quoted in a bond-order-potential paper
  (arXiv:2208.12973, via Park et al. 2004), i.e. GGA, secondary source; B2 FeAl ~0.7 muB per f.u. FM = VASP-PBE (arXiv:1605.06967) vs paramagnetic
  experiment.  Our code is LDA/ASA: moments will differ by a few % (GGA > LDA); FeAl is functional-sensitive and should not be a pass/fail gate.
- First run (not yet a test): B2 FeCo FM, a = 2.857 A, equal spheres, tetrahedra nk 16, `example/feco_b2.toml`: converged in 43 its,
  Fe 2.760 muB, Co 1.768 muB (total 4.53 per f.u.), qnet = -/+0.103 e (Co loses to Fe), E_F -0.0593 Ry.  vs VASP-PBE 2.80 / 1.81 / 4.61:
  -1.4% / -2.3% / -1.8%.  Multi-class spin SCF in the ordered code works without changes.

## Lattice-constant sensitivity, CuZn (EI route, depth 4, rank 1) and FeCo (2026-10-04)
- CuZn: a = 2.85 / 2.9206 (mp-987) / 2.95 A -> see table below (E_F rel. peaks); our default 2.85 A is 3.5% below experiment (2.95-2.96) and its
  provenance ("published a = 2.85") is not confirmed by the reference PDFs.
- B2 FeCo FM at a = 2.844 A: Fe 2.743, Co 1.764, total 4.507 muB (a = 2.857: 2.760 / 1.768 / 4.528); literature totals 4.53 (a = 2.844) and 4.61.

## Run report `asr.out` (2026-10-05)
`asr_report` module; written by `kind = "lmto"` and `kind = "alloy"` next to the input (`[output] out_file`, default `asr.out`).
Sections (grep-able `=== name ===`): banner (version, date, threads), verbatim input, method, structure (volume, S, filling),
SCF history (dV, E_F, moments, qnet per iteration), converged per-class results (qnet, E_nu, C, sqrt(Delta), p, gamma, q0/q1/q2,
s/p/d occupations), N(E_F) per species and l, timing, status.  Checked on B2 CuZn (a = 2.92055 A, converged in 16 its) and the
bcc CuZn alloy (depth 2): occupations sum to nval + qnet exactly.  Not yet covered by an automated test.
Gotcha: alloy_scf stores q per spin channel, so the alloy driver doubles q for spin = 1 before reporting.
Gotcha: `build/ifx_*/app/asr` can be a stale binary; run via `asr_run --target asr -- input.toml`.
Planned: HDF5 for DOS/spectral/history (needs hdf5-devel headers and h5py; neither is installed).

## Ordered B2 in the real-space (recursion) driver (2026-10-05)
`kind = "alloy"`, `[structure] lattice = "b2"` (no concentration): A on the even / B on the odd sublattice of the bcc graph; real-space
recursion on the periodic supercell (`lmto_real_op`, new `build_lmto_ordered`), no augmented space, no flips; Lanczos started on an A
and on a B site (4 orbitals x spin x sublattice); exact two-site Madelung (`madelung_phi`) instead of SIM; EI route only.
Test `test_alloy_b2` (4 assertions): at the SAME converged parameters, k-space H^(2) = E + h - h o h (`ham_k_cell`, nk 24) gives E_F within
5.9 / 4.9 / 2.6 mRy of the recursion at depth 16 / 20 / 24 (supercell 12-16; supercell size irrelevant).  => the real-space recursion and
terminator are accurate to a few mRy for B2 CuZn.
B2 CuZn a = 2.92055 A (mp-987 CIF), depth 32, supercell 16: E_F -0.0503 Ry, Cu +-0.137 e; Cu d / Zn d peaks relative to E_F -0.238 / -0.619 Ry.
First-order k-space ASA SCF (`kind = "lmto"`): E_F -0.0809, Cu d / Zn d -0.206 / -0.555 Ry, Cu +0.135 e.  The ~30 mRy (E_F) and ~60 mRy (Zn d)
difference is the Hamiltonian (second-order H^(2) vs first-order generalised problem), not the recursion: the real-space route reproduces H^(2),
and the first-order k-space route is closer to Materials Project (Cu d -2.8/-3.4 eV: k-space -2.80/-3.47, real-space about -3.2/-3.3).
B2 FeCo FM (a = 2.857 A, depth 24, supercell 12, zeta0 0.1, mixing 0.2): converged only after 124 iterations (wanders for ~100: E_F up to +0.6 Ry,
qnet up to +-1.4 e; 120 iterations were not enough): Fe 2.745, Co 1.716 muB (total 4.461) vs Elk LSDA 2.75/~1.75 and first-order k-space 2.760/1.768.
SCF stability of the magnetic ordered case is poor and should be improved (mixing / clamps) before using it for alloys; depth convergence of the moments is untested.

## Energy-dependent recursion for ordered B2 (2026-10-06)
Existing before: energy-independent (H^(2)) for ordered B2 and disordered alloys; energy-dependent (P(z)-S, `asr_endep`) for disordered
alloys only (augmented space, needs symmetry = true) and ordered fcc/bcc test cases.  New: `ed_ordered_green` (asr_endep) = complex-symmetric
Lanczos on A(z) = P_site(z) - S for one site of the ordered real-space operator (`lmto_s_apply`, new, in asr_lmto_realization), wired into
`alloy_scf` (ordered + energy_dependent: closed contour, secant E_F search, species LDOS) and the TOML driver (`lattice = "b2"` with
`hamiltonian = "energy-dependent"`).  Examples: `example/cuzn_b2_ed.toml` (ordered), `example/cuzn_alloy_ed.toml` (disordered, depth 3).
- Test `test_endep_ordered` (new): recursion (10 levels) vs dense inverse of P(z)-S on the same 128-site periodic B2 supercell, 3 complex
  energies x 2 sublattices x 4 orbitals: max relative deviation 5.4e-3 (< 2% gate).
- B2 CuZn a = 2.92055 A, depth 12, supercell 6, 24-pt contour (smoke run, local): ED converged in 53 its, Cu +0.154 / Zn -0.155 e, E_F -0.0188 Ry,
  Cu d / Zn d peaks vs E_F -0.231 / -0.656 Ry, occupations sum to nval + qnet, d occupations < 10 (the normalisation problem of the
  disordered ED route is not visible in the ordered case).  Same depth, EI: converged, Cu +0.138, E_F -0.0541, peaks -0.246 / -0.608.
  So ED vs EI differ by ~35 mRy (E_F) and ~50 mRy (Zn d) at depth 12 -> a depth/convergence study (24, 32) is running on the desktop
  (~/asr_remote4/trial/ord, trial/dis).
- Desktop production (2026-10-06): B2 CuZn ED depth 24, supercell 12, 32-pt contour: converged in 48 its (32 min), E_F -0.0790 Ry,
  Cu +0.120 e, Cu d / Zn d peaks vs E_F -0.221 / -0.591 Ry.  First-order k-space ASA at the same a: E_F -0.0809, -0.206 / -0.555, Cu +0.135;
  EI H^(2) recursion depth 32: E_F -0.0503, -0.238 / -0.619, Cu +0.137.  => the ED recursion reproduces the k-space ASA E_F to 2 mRy, the
  EI (second-order H) route is 29 mRy off: the EI/ASA difference is the H^(2) truncation, as suspected.
  Disordered CuZn ED depth 3 (rank 2, 32-pt contour): converged in 44 its (20 min), E_F 0.0456, Cu +0.048 e, Cu d -0.178 / Zn d -0.706 Ry
  vs E_F (same as the earlier depth-3 queue run: the ordered-route changes did not disturb the disordered ED path).
- Disordered ED remains as documented above (depth-3 CuZn only; depth >= 4 not normalised).  Cost: ~25 s per E_F evaluation at depth 3 (63k nodes, 18 threads).

## Why the recursion (EI) route over-estimates the Fe moment and needs depth for Ni (2026-10-06, desktop, `probe_moment_routes`)
Pure ferromagnets through the alloy route (A = B, rank 0, supercell 36, 60-pt contour) vs our k-space ASA (generalised problem, tetrahedra
nk 24, same vBH LSDA): Fe 2.226 muB (E_F -0.0926), Ni 0.632 (E_F -0.1014).  Legacy TB-LMTO-ASR (Chakrabarti & Mookerjee, EPJB 44, 21 (2005),
Table 1, "Pure Fe"): 2.25 muB (provenance of that row not stated; for a pure metal the paper's effective Hamiltonian reduces to H^(2)
because the disorder operators B, F are proportional to Z^A - Z^B).
| Fe, EI route | moment | E_F |
|---|---|---|
| depth 8 / 12 / 16 / 24 (sigma 0.7, rcl 1.2) | 2.262 / 2.344 / 2.351 / 2.344 | -0.082 / -0.091 / -0.093 / -0.094 |
| depth 8 / 12 / 16, Beer-Pettifor tail | 2.270 / 2.319 / 2.352 | |
| depth 12 / 16, Stuttgart density n = (M0 - p M2) phi^2 + 2 M1 phi phidot + M2 phidot^2 | 2.349 / 2.357 | |
| depth 12, sigma 0.65 / 0.70 / 0.75 (rcl 1.2) | 2.280 / 2.344 / 2.485 | |
| k-space ASA, sigma 0.65 / 0.70 / 0.75 | 2.2259 / 2.2261 / 2.2260 | |
| depth 12, sigma 0.7, kept range rcl 1.2 / 1.5 / 2.0 | 2.344 / 2.334 / 2.294 | |
| depth 12, sigma 0.7, rcl 1.2, inversion cluster 2.5 | 2.335 | |
| depth 12, sigma 0.64 (alpha_d 0.0107 = standard TB-LMTO), rcl 1.2 / rcl 1.2 + inversion 2.5 / rcl 2.0 + inversion 2.5 | 2.278 / 2.274 / **2.240** | -0.092 / -0.093 / -0.090 |
Conclusions (each tested, not assumed):
- NOT the cause: recursion depth (Fe converged by 12 levels), terminator (BP == avg2 to 0.03), density formula (Stuttgart vs phi(E) expansion: 0.005).
- CAUSE: the second-order H = E_nu + h - h o h is representation dependent; with our hard-core sigma = 0.7 (alpha = 0.350, 0.0572, 0.0168) the Fe moment is
  0.12 too large and moves by 0.2 muB between sigma 0.65 and 0.75, while the exact ASA is invariant to 1e-4.  The references use the standard TB-LMTO
  constants (alpha_s 0.3485, alpha_p 0.05303, alpha_d 0.010714; Andersen-Jepsen; the thesis prints alpha_d as "0.00107", a typo); alpha_d = 0.0107 is
  sigma ~ 0.64 for d.  Second cause: the kept S^alpha range rcl = 1.2 alat (bcc: 2 shells, 14 neighbours) is short; 2.0 alat (64 neighbours) lowers
  the moment by 0.04-0.05.  The inversion-cluster size matters little (<= 0.01).  sigma 0.64 + rcl 2.0 gives 2.240 vs ASA 2.226 (0.6%).
- Ni: the moment collapses to 0 at depth <= 12 with either terminator (avg2 or BP), 0.28 at 16 (BP), 0.53 at 20, 0.58 at 30 (sigma 0.7); with sigma 0.64 +
  inversion 2.5: 0.566 at 20, ~0.615 at 30 (E_F -0.107 vs ASA -0.101).  Ni needs both the better representation and depth >= 30 (its E_F sits on the
  d-band-top peak, which 12 levels smear below the Stoner threshold).  The legacy "11 steps" claim (Kaphle et al. 2012) is for Mn-based moments.
- k-space H^(2) (new diagnostic `asa_input_t%second_order`): converged at nk 12 to 2.29 (E_F -0.138), DIVERGED at nk 24 -> not used as evidence.
- Energy-dependent route for Fe: augmented-space ED (A = B) at depth 8 does not converge (E_F wanders -0.53..+0.60 Ry); ordered real-space ED (route `edo`)
  depth 16/24 running (slow; E_F still wandering after 4 iterations).
New options (defaults unchanged): `alloy_input_t%terminator` ('avg2' | 'bp'), `%density` ('expand' | 'stuttgart'), `%rinv` (inversion cluster radius).

## Default changed: standard TB-LMTO screening in the recursion route (2026-10-06, all runs on the desktop)
`alloy_input_t%screening = 'auto'` (TOML `[asr] screening`, default): the energy-independent H^(2) recursion (alloy / ordered real-space) uses the
standard constants alpha = 0.3485, 0.05303, 0.010714 (Andersen & Jepsen 1984; `standard_alpha` in asr_strconst) instead of hard-core sigma = 0.7;
the energy-dependent P(z) - S route keeps hard-core sigma (see ED note below).  'standard' / 'hardcore' force either choice.
The k-space ASA (`asa_scf`) still uses sigma (it is representation independent: 2.2259/2.2261/2.2260 for sigma 0.65/0.70/0.75).
Kept S range with standard alpha (Fe depth 12 / Ni depth 20, inversion radius 3.0): rcl 1.2 / 1.5 / 1.75 / 2.0 / 2.3 / 2.5 -> Fe 2.258 / 2.261 / 2.255 /
2.254 / 2.254 / 2.254; Ni rcl 1.2 / 1.45 / 1.6 / 1.75 -> 0.539 / 0.541 / 0.540 / 0.540.  Inversion radius 2.5 / 3.0 / 3.5: identical.  => the default
rcl = 1.2 is converged to 0.004 muB with standard alpha (the strong rcl dependence was a sigma = 0.7 artefact); no cost increase.
Full suite on the desktop after the change: 35/35 test programs, 241 assertions pass.  Changes visible in tests: test_alloy_scf_spin Fe alloy route
2.285 vs ordered ASA 2.216 (was ~2.36); test_alloy_b2 recursion vs k-space H^(2) E_F 0.2 mRy at depth 16 (was 5.9 mRy).
Re-baselined (EI unless noted; peaks are relative to E_F, Ry):
| system | before (sigma 0.7) | now (standard alpha) | reference |
|---|---|---|---|
| bcc Fe, depth 12 / 16 / 24 | 2.344 / 2.351 / 2.344 | 2.258 / 2.280 / 2.261 | k-space ASA 2.226; legacy ASR 2.25 |
| fcc Ni, depth 20 / 30 / 40 | 0.534 / 0.578 / - | 0.540 / 0.594 / 0.612 (E_F -0.1033 at 40) | k-space ASA 0.632 (E_F -0.1014) |
| B2 FeCo real-space depth 24 (Fe / Co) | 2.745 / 1.716, 124 its | 2.743 / 1.710, 51 its | k-space ASA 2.760 / 1.768 |
| B2 CuZn real-space depth 24: E_F, Cu d, Zn d, Cu qnet | -0.0503 (d32), -0.238, -0.619, 0.137 | -0.0447, -0.244, -0.619, 0.137 | k-space ASA -0.0809, -0.206, -0.555, 0.135 |
| CuZn bcc disordered depth 4 rank 1: E_F, Cu d, Zn d, qnet | 0.0196, -0.335, -0.683, 0.029 | 0.0220, -0.330, -0.686, 0.029 | Banerjee 2010: Cu d -0.309, Zn d -0.682 |
| CuZn bcc disordered depth 3 rank 2: E_F, qnet | ~0.02, 0.024 | 0.0171, 0.024 (d peaks merged at depth 3, as before) | |
Remaining gaps: (1) the H^(2) truncation still puts B2 CuZn E_F ~36 mRy above the k-space ASA (sp states far from E_nu); the ED route gets it to 2 mRy;
(2) Ni needs depth >= 40 for 0.61 (ASA 0.632); (3) B2 FeCo Co moment 3% low, unchanged by the screening (not yet explained).
Energy-dependent route and the screening (desktop, 2026-10-06):
- Ordered real-space ED, bcc Fe (sigma 0.7): depth 16 converged 2.254 muB, E_F -0.0917; depth 24 2.252, E_F -0.0922 (k-space ASA 2.226, -0.0926).
  The ED route (representation independent in principle) agrees with the corrected EI route (2.26) and the legacy ASR (2.25); E_F within 1 mRy of the ASA.
- CuZn ED with the standard alpha: neither the ordered (depth 24) nor the disordered (depth 3) SCF converged in 17+ iterations (E_F swinging 0.0..0.24 Ry,
  charge swings up to 1 e), while sigma 0.7 converges (48 / 44 its).  Likely mechanism: the real pole of P^alpha(z) at z_p = C - Delta/(gamma - alpha);
  for Zn d it moves from ~ -0.33 Ry (inside the Cu d band) to ~ -0.13 Ry (just below E_F) with the smaller alpha_d; the truncated recursion does not cancel
  it in lambda + mu g mu.  Hence 'auto' keeps sigma 0.7 for ED.  With 'auto' the ED CuZn runs reproduce the earlier converged runs digit for digit.

## B2 FeCo Co moment: third-order Hamiltonian + depth (2026-10-06/07, desktop)
External targets (LSDA, same sphere-free quantities differ by method): Elk FLAPW Fe 2.75 / Co ~1.75 (a = 2.856, arXiv:1907.06613); ASW-LSDA
(Schwarz et al. 1984) Co 1.73-1.80 across bcc Fe-Co; ecalj LDA 4.44 muB/cell, experiment 4.70/cell (arXiv:1904.09093); GGA 2.80/1.81.
Our k-space ASA (nk 24): Fe 2.770, Co 1.774, E_F -0.0600.  Recursion (order 2, depth 24): 2.743 / 1.710 (Co 0.064 low).
Diagnosis (`probe_feco_terms`, one-shot at fixed k-space ASA parameters, orthonormal expansion of the exact ASA, K = h(1+oh)^-1):
| H | Fe | Co | E_F |
|---|---|---|---|
| full (1+KpK)^-1/2 K (1+KpK)^-1/2 | 2.784 | 1.792 | -0.0602 |
| K (all orders in o, no p = ED route) | 2.750 | 1.758 | -0.0536 |
| H^(2) = h - hoh | 2.743 | 1.751 | -0.0552 |
| H^(3) = H^(2) + hohoh - (hphh + hhph)/2 | 2.775 | 1.785 | -0.0638 |
=> ~80% of the H^(2) deficit is the missing overlap term h p h.  Same term explains the B2 CuZn E_F gap (H^(2) -0.0420, H^(3) -0.0741, full -0.0764).
New: `order = 3` ([asr] order; alloy_input_t%order) in all three recursion operators (lmto_real_op, augorb, augsym; augmented-space closure over
three-hop paths with flips on up to 4 sites).  Gates: test_lmto_real3 (matrix-free == dense, 2e-15), test_augorb ring (order 3 == brute-force
configuration average, 8e-15), test_augsym_orb (reduced == unreduced at order 3), test_alloy_b2 (order 3 E_F -0.08073 vs k-space ASA -0.08071 at depth 16;
order 2 is 33 mRy off).  Order 3 is unstable from the atomic start (cubic terms with E_nu far from band centres): staged SCF, order 2 until
dV < order_switch (1e-3), then order 3 (Anderson history reset).
Depth (`probe_feco_depth`, fixed k-space parameters, order 3): E_F / Co = -0.074 / 1.756 (16), -0.073 / 1.755 (24), -0.069 / 1.747 (32), -0.063 / 1.758 (48),
-0.063 / 1.781 (64) vs k-space H^(3) -0.0638 / 1.778: FeCo needs depth >= 48-64.
Self-consistent B2 FeCo, order 3: depth 24 2.752 / 1.741; depth 48 2.805 / 1.764; depth 64 **2.810 / 1.779** (E_F -0.0724).  Co now matches the k-space ASA
(1.774) and is within the LSDA literature range; Fe is 0.04 above the k-space ASA (2.770) and Elk (2.75) -- open.
Tried and rejected as default: sphere-charge start vectors (1 - hph/2)|L> (`charges = "sphere"`): fixes the s/d partition at fixed parameters (Lowdin weights
shift ~0.05 e from Co d to s) but loses ~0.03 e/cell at this order (E_F +12 mRy in CuZn) and did not move the SCF Co moment (1.740 vs 1.741).
Process note: fpm did not relink test executables after module type changes -> 5 tests segfaulted on stale objects; after `rm -rf build/ifx_*`
36/36 test programs pass (desktop).

## B2 FeCo on the energy-dependent route (2026-10-07, desktop; ED is the target production route)
Ordered real-space ED (P(z) - S, hard-core sigma 0.7 via screening 'auto', 32-pt closed contour, supercell 12):
depth 24: converged in 75 its, Fe 2.793, Co 1.763 muB, E_F -0.0582 Ry; depth 48: dV stalls at ~2e-4 after 132 its, Fe 2.786, Co 1.760, E_F -0.0585
(converged into the spin-reversed FM state, same magnitudes).  vs k-space ASA 2.770 / 1.774 / -0.0600; Elk FLAPW 2.75 / ~1.75.
=> ED: E_F within 2 mRy of the ASA, Co within 0.012, Fe +0.02, depth-converged at 24 (EI needed order 3 and depth 64 for the Co moment).
Remaining ED gap is of the size expected from the missing p (overlap) term in the 2-parameter P(z).  Cost ~1.5 min per SCF iteration at depth 24 (14 threads).

## Energy-dependent route: adaptive depth works; third-order P(z) fails (2026-10-07, desktop)
- Adaptive depth (`depth_tol`, default 1e-5; `levels` = maximum): each ED chain stops when the terminated G changes < depth_tol between checks every 4 levels
  from level 8.  B2 FeCo ED, max 48 levels: mean depth used 27.4, converged in 82 its, Fe 2.786 / Co 1.760 / E_F -0.0583 = the fixed-depth-48 result
  (2.786 / 1.760 / -0.0585) and depth 24 (2.793 / 1.763).
- Exact potential function (probe_pfun, bcc Fe from the radial solver): two-parameter P_d is off by 6-8% within +-0.2 Ry of E_nu, third-order
  (z - C + p eps^3)/Delta by ~1%; s, p third-order essentially exact; exact P_d has real poles ~ -0.85 / +0.65 Ry from E_nu.
- Third-order P(z) in the ED SCF (staged after dV < 1e-3): diverges for B2 CuZn (qnet -> 10 e), B2 FeCo and bcc Fe (moments -> 10 muB) within a few
  iterations of the switch.  Disabled (`ed_third_order = .false.`); the ED route keeps the two-parameter P.  Proper fix: exact P_l(z) and phi_l(z, r)
  from the radial equation at each complex contour energy (true KKR-ASA Green function, no linearisation).

## D0: exact potential function on the energy-dependent route (2026-10-07/08, desktop)
New `asr_pfun` (P_l(z), dP/dz from the variational equations, phi_l(z,r)), `test_pfun` 5/5.  Lessons (each verified):
lambda must be C'/C (log-derivative of P's denominator), NOT -P''/(2P') (spurious upper-half-plane poles for a meromorphic P: Cu/Zn d counts -> 0);
the contour bottom must be away from the real poles of exact P (`choose_ebot`); exact P needs ~64 contour points (32: noisy N(E));
the density -1/pi Im int G phi(z)^2 alone is incomplete (needs the single-site irregular-solution term; dV ~ 6-8 Ry) -> density from the
energy moments of the exact G for now (`ed_exact_density = .false.`).
Ordered results (exact P, moment density, adaptive depth <= 48, 64-pt contour):
| system | exact P | parametrised P (ED) | k-space ASA |
|---|---|---|---|
| B2 CuZn: E_F, Cu qnet | -0.1011, +0.130 | -0.0790, +0.120 | -0.0807, +0.136 |
| B2 FeCo: Fe, Co, E_F | 2.801, 1.802, -0.0591 | 2.793, 1.763, -0.0582 | 2.770, 1.774, -0.0600 |
| bcc Fe: moment, E_F | 2.261, -0.0906 | 2.252, ~-0.092 | 2.226, -0.0921 |
D0 gate (b) NOT passed: exact P does not bring the ED route closer to the k-space ASA (CuZn E_F 20 mRy low; moments +0.03).  The remaining
linearisation is the moment-built density (phi, phidot, phiddot at E_nu) -> next: complete exact density (irregular solution).
D0 continued (2026-10-08): complete exact density.  rho_l(r) = -1/pi Im int [ sum_m G_m(z) phi_l(z,r)^2 + (2l+1) T_l(z,r) ] with the single-site
term T = f - (int f) phi^2, f = (P_phi P_eta + Q_phi Q_eta/c^2)/W (G = (z - H)^-1 convention), eta the irregular solution with log-derivative
D* = (l + alpha K (l+1))/(1 - alpha K) at S.  T is required: the normalised phi^2 = phi_raw^2/N(z) has poles at the complex zeros of the analytic
norm N(z), and the lambda phi^2 part picks them up (Zn density negative in the core region, right total charge); T cancels them.
Two errors found on the way: T omitted (dV 6 Ry), T with the (H - z)^-1 sign (dV 12 Ry).  test_pfun 7/7 (T real on the real axis, no charge).
B2 CuZn smoke (depth 12, 64-pt contour): dV 0.157 -> 0.096 -> 0.067 (moment density: 0.157 -> 0.133 -> 0.140).  Desktop gate runs: trial/d0x.
D0 gate result (2026-10-08, desktop, trial/d0x): exact P + exact density, adaptive depth <= 48 (tol 1e-5), 64-pt contour, and the new
same-model check `kkr_check = 24` (k-space KKR-ASA, BZ sum of [P(z) - S(k)]^-1 with the same exact P, S^alpha, contour, at the final potentials):
| system | E_F k-space / recursion | occupations | moments k-space / recursion |
|---|---|---|---|
| bcc Fe | -0.08915 / -0.08936 | all within 0.005 e | 2.268 / 2.260 |
| B2 FeCo | -0.05681 / -0.05770 | within 0.016 e (Fe minority d) | Fe 2.806 / 2.801, Co 1.803 / 1.804 |
| B2 CuZn | -0.09098 / -0.09669 | d within 0.005, p recursion +0.016-0.017 e | - |
=> the ED recursion reproduces the exact KKR-ASA of the same model (magnetic: <1 mRy, 0.01 muB; CuZn: 5.7 mRy from the slowly converging p channel).
=> the earlier "gaps" to our k-space LMTO-ASA (CuZn E_F 10-16 mRy, Fe moment +0.03-0.04) are the LMTO linearisation of the k-space code, not ED errors.
Exact KKR-ASA values: bcc Fe 2.26 muB; B2 FeCo 2.80 / 1.80; CuZn E_F -0.091 (all-electron LSDA literature: Fe ~2.2, FeCo 2.75 / ~1.75 -- ASA itself
is ~0.03-0.05 muB high here).  Open: CuZn p-channel convergence (depth_tol / terminator for the wide sp band) -> D1 setup.

## D1 first pass (2026-10-08/09, desktop trial/d1): disordered CuZn, exact ED (exact P + density), sigma 0.7 screening
| depth/rank | E_F | Cu d / Zn d occ. | Cu d / Zn d peak vs E_F | min real-axis DOS |
|---|---|---|---|---|
| 3 / 2 | -0.0236 | 9.64 / 10.06 | -0.183 / -0.599 | -1.9 |
| 4 / 1 | -0.0425 (not conv., 80 its) | 9.69 / 9.97 | -0.218 / -0.623 | -11.5 |
| 5 / 1 | -0.0577 | 9.56 / 10.04 | -0.149 / -0.576 | -15.0 |
Gross over-count gone (old route: Cu d 11.7 e at depth 4), BUT gate failed: negative real-axis DOS growing with depth, real-axis vs contour
occupations disagree (Cu d 3.06 vs 4.78 per spin at depth 5), E_F not depth-converged (34 mRy), depth-4 SCF unconverged.
Diagnosis: real pole of the exact P^alpha_d of Zn inside the valence window (~ -0.54 Ry with sigma-0.7 alpha_d).  In augmented space the start
node carries <P> = x P_A + y P_B, so a truncated chain cannot cancel a single species' pole in lambda_A + mu_A^2 <n g>/x (in the ordered case the
start site carries the exact P of its species and the cancellation holds).  Next (D1 pole control): choose alpha_l per l (alpha_d <= 0 moves
the d pole to just above the lower pole of P^0, ~ -2.1 Ry) so that no P^alpha pole of any species lies in [ebot, E_F + 0.3].

## Decision 2026-10-09: ASR as in AMscf (not KKR-ASA)
AMscf's energy-dependent ASR (ASR/main.f90 ap1..ap13, hop_endep.f90): H(E) = E - T(E - E_nu - h)T - T hoh T, T = Delta_c^1/2 Delta^-1/2 (linear
P(E) = (E - C)/Delta scaled by the central species' Delta, central site fixed to one species, seeds + spline, BP tail).  Since
E - H(E) = T (E - H^(2)) T, its central-site G equals the H^(2) Green function exactly (only the Lanczos basis differs).  => production disorder
runs use the energy-independent augmented-space route (exact Delta disorder, species via transition elements); the exact KKR-ASA ED stays an
ordered-system reference.  D1 pole control is dropped.

## D2 (2026-10-09, desktop trial/d2): disordered bcc CuZn, EI route, depth 4 rank 1, a = 2.85 A
| x_Cu | E_F | Zn d / Cu d peak vs E_F (ours) | Banerjee & Mookerjee 2010 Fig. 3 (digitised) | min DOS |
|---|---|---|---|---|
| 0.24 | 0.0922 | -0.750 / -0.444 | -0.579 / -0.224 | 0.505 |
| 0.50 | 0.0220 | -0.686 / -0.330 | -0.682 / -0.309 | 0.406 |
| 0.75 | -0.0413 | -0.621 / -0.261 | -0.790 / -0.387 | 0.271 |
Panel labels checked (caption: left Cu0.24Zn0.76, middle Cu0.5Zn0.5, right Cu0.75Zn0.25; the Cu0.24 panel has the big deep Zn-d peak).  But the
published DOS integrated to E = 0 gives 11.57 / 11.68 / 12.29 electrons vs 11.76 / 11.50 / 11.25: the Cu0.75 panel is +1.04 e off, so its energy
zero is not E_F and its peaks are too deep.  Our trend (Zn-rich -> d bands deeper below E_F) is the rigid-band expectation; at 50:50 (reference
electron count +0.18) we agree to 4 / 21 mRy.  The concentration trend is not a gate against this figure; needs a second reference (KKR-CPA).
Order-3 runs (same inputs) pending: node generation over three-hop paths is slow (single-threaded).
