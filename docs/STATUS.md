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
- [WRONG FIX, corrected 2026-10-14, see "SIM sign and magnitude error"] SIM shift had the wrong sign (V = -beta 2 alpha_M qnet/S with qnet = N - Z gave positive feedback: charge transfer grew
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
E - H(E) = T (E - H^(2)) T, its central-site G equals the H^(2) Green function exactly (only the Lanczos basis differs).  => the user's decision:
production = AMscf's energy-dependent ASR in its literal form (to be implemented); the EI augmented-space route and the exact KKR-ASA ED are
cross-checks/references.  D1 pole control (of the KKR-ASA route) is dropped.  [Correction 2026-10-09: an earlier draft of this entry said
"production uses the EI route"; that contradicted the user's instruction and was withdrawn.]

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

## AMscf-form energy-dependent ASR implemented (2026-10-09): asr_ed, `hamiltonian = "asr_ed"`
Operator (AMscf ap1..ap13 in our augmented space): for central species c (start site frozen, non-flippable, M_R -> 1 or 0 there; one operator and one
recursion per species and start orbital) H(E) = [E - D_c (E - E_nu,R)/D_R] + h' - h' q h', h' = g + S', S' = D_c^1/2 S D_c^1/2 (non-random),
g_R = (D_c/D_R)(C_R - E_nu,R), q_R = o_R D_R/D_c (D = Delta^alpha).  It is the existing H^(2)-type operator with c -> g, d -> sqrt(D_c), o -> q and the E_nu term replaced
by the complex E-dependent diagonal; the real part is applied to Re/Im of the Lanczos vector, the diagonal through M (`linop%mapply`);
complex-symmetric chain per contour energy, eigenvalue-hull tail (`cf_species`).  G_c(E) = <1|(E - H(E))^-1|1> is used directly (no lambda/mu) in the
EI moment formulas.  Default screening for `asr_ed` = standard TB-LMTO alpha (the form has no P^alpha poles).  Code: src/asr_ed.f90,
`augorb/augsym_freeze_site`, `alloy_input_t%ed_asr`; example/cuzn_alloy_asr_ed.toml.
Gates (all pass): test_asr_ed (ring, 9 orbitals, exact configuration average of [E-H(E)]^-1_00 for both central species == dense H(E) == dense
[E-H^(2)]^-1_00 to 1e-15, i.e. E - H(E) = T(E - H^(2))T with T_c = 1); test_asr_ed_sym (D4h-reduced == unreduced chain 2e-11, fcc 94702 -> 6843
nodes; ordered limit A = B: ED chain == H^(2) chain 5e-14); sum rule M0(E_F = 6 Ry) = 1.0000 for all l, both species.  Desktop: clean rebuild, all 39 test
programs pass.
First production point (desktop trial/amscf/x50_l4r1; bcc Cu50Zn50, a = 2.85 A, levels 4, rank 1, supercell 8, 32 contour points; 14.5 min, 39 its,
17515 nodes) vs the EI cross-check at the same depth (trial/d2/x50_o2, order 2):
| | E_F (Ry) | qnet Cu | Cu d / Zn d C (Ry) | d-peak vs E_F Zn / Cu | min DOS |
|---|---|---|---|---|---|
| AMscf ED | 0.0145 | +0.044 | -0.254 / -0.677 | -0.728 / -0.278 | 0.48 |
| EI order 2 | 0.0220 | +0.029 | -0.279 / -0.641 | -0.686 / -0.330 | 0.41 |
DOS positive everywhere.  Banerjee-Mookerjee 2010 (50:50): Zn d -0.682, Cu d -0.309 vs E_F; Zn d C -0.677 vs their -0.682.  The two routes are equal at
infinite depth, not at depth 4: differences here are depth/rank-1 truncation (the ED peaks are broader/lower: Zn d 47 vs 70 states/Ry).
Open: depth convergence (levels 5-6, rank 2) of both routes at x = 0.5; then x = 0.24 / 0.75.

## asr_ed vs EI depth/rank convergence, CuZn x = 0.5 (2026-10-09/10, desktop trial/conv_x50, trial/amscf, trial/d2)
bcc Cu50Zn50, a = 2.85 A, supercell 8, 32 contour points (asr_ed), EI = order 2.  Peaks relative to E_F; C_d = gamma-representation band centre (Ry).
| run | E_F | qnet Cu | Cu C_d | Zn C_d | Zn d peak (height) | Cu d peak | min DOS | nodes | wall |
|---|---|---|---|---|---|---|---|---|---|
| EI  l4 r1 | 0.0220 | 0.0286 | -0.2791 | -0.6405 | -0.686 (70) | -0.330 | 0.41 | 17630 | min |
| EI  l4 r2 | 0.0220 | 0.0286 | -0.2790 | -0.6406 | -0.686 (70) | -0.330 | 0.41 | 554471 | 2 min |
| EI  l5 r1 | 0.0071 | 0.0334 | -0.2636 | -0.6462 | -0.685 (49) | -0.259 | 0.48 | 40817 | 20 s |
| ED  l4 r1 | 0.0145 | 0.0444 | -0.2544 | -0.6769 | -0.728 (47) | -0.278 | 0.48 | 17515 | 14.5 min |
| ED  l4 r2 | 0.0146 | 0.0445 | -0.2543 | -0.6769 | -0.729 (47) | -0.279 | 0.48 | 536956 | 9 h 43 min |
| ED  l5 r1 | 0.0170 | 0.0367 | -0.2566 | -0.6705 | -0.733 (36) | -0.287 | 0.49 | 40702 | 35 min |
Findings: (1) rank 1 -> 2 changes nothing at level 4 on either route (E_F < 0.1 mRy, C_d < 0.1 mRy, peaks < 1 mRy): rank 1 is enough here; rank 2 costs
31x the nodes (ED: 40x the time) for no gain.  (2) Depth 4 -> 5 (rank 1): EI moves E_F by 15 mRy, Cu d C by 15 mRy, Cu d peak by 70 mRy, Zn d peak
height 70 -> 49; asr_ed moves E_F by 2.6 mRy, C_d by 2-6 mRy, peaks by 5-9 mRy: asr_ed is much better converged in depth, as expected for a
recursion whose non-random hops start from the right site.  (3) At level 5 the routes approach each other (qnet 0.0334 vs 0.0367, Cu C_d -0.264 vs
-0.257) but Zn C_d differs by 25 mRy and E_F by 10 mRy; they are equal only at infinite depth.  Level 6 not yet run (cost per level ~2.4x: ~85 min).
(4) The disordered Zn d peak stays at -0.73 Ry (-9.9 eV) vs E_F in asr_ed (EI -0.69), i.e. 0.13-0.18 Ry deeper than the ordered B2 values (-0.555 k-space ASA,
-0.619 real space, LAPW-PBE Keast 2015 Fig. 2c about -7.3 eV).  Not a depth/rank artefact.  Open: those ordered numbers are at a = 2.92 A, these runs at
2.85 A => run ordered B2 at 2.85 A to separate lattice constant from disorder.

## Lattice constant vs disorder for the CuZn Zn d position (2026-10-10, desktop trial/q2): ordered B2, real-space EI, levels 24, supercell 16
| a (A) | E_F (Ry) | Cu qnet | Cu C_d | Zn C_d | Zn d peak vs E_F | Cu d peak |
|---|---|---|---|---|---|---|
| 2.85 | 0.0024 | +0.138 | -0.277 | -0.634 | -0.636 (-8.7 eV) | -0.259 |
| 2.92055 | -0.0447 | +0.137 | -0.307 | -0.664 | -0.619 (-8.4 eV) | -0.244 |
The lattice constant moves the Zn d peak (vs E_F) by only 17 mRy, so it does not explain the shift to the disordered alloys at 2.85 A: EI x = 0.5 depth 4-5: -0.686,
asr_ed depth 4-5: -0.73 Ry.  Disorder deepens Zn d (vs E_F) by 50 mRy (EI) to ~95 mRy (asr_ed) relative to the ordered EI value at the same a (depths differ:
24 vs 4-5 levels).  LAPW-PBE ordered (Keast et al. 2015, Fig. 2c, read by eye) -7.3 eV.

## n-component augmented space (2026-10-10): asr_augn, asr_edn
`asr_augn`: random species s_R in {1..ns} (probabilities x_s, uncorrelated).  Local basis phi_0 = 1, phi_1..phi_{ns-1} (Gram-Schmidt on species indicators),
indicator-of-s matrix M^(s)_lk = x_s phi_l(s) phi_k(s) (rank one), site parameter X_R -> W_lk = sum_s X_s M^(s)_lk; nodes = (site, flipped sites with states 1..ns-1),
closure over paths R -> R1 -> R2 with every visited site in any of its ns states (ns^3 combinations), rank limit, non-flippable sites carry the average W_00, a
site may be frozen to one species.  `asr_edn` = the asr_ed form for n species (T-scaled parameters per central species, per-species complex diagonal).
Gates (all pass): test_augn: ring N = 4, nl = 2, ns = 3: nodes N ns^N, recursion == exact configuration average of the resolvent 2e-15; frozen start site (species 2)
2e-15; ns = 2 == asr_augorb (same node count, coefficients 3e-16); complex-parameter application.  test_asr_ed_n: ring, 9 orbitals, three species (0.3/0.5/0.2),
each central species: ED recursion == dense average of [E - H(E)]^-1 == dense average of [E - H^(2)]^-1 to 3e-15.
Not yet: D4h point-group reduction of the n-component space (asr_augn is unreduced), alloy_scf with ns species (driver is two-species), short-range order.

## asr_ed depth series to level 7, concentration series, magnetic validation (2026-10-10/11, desktop trial/q2, trial/q3)
CORRECTION of the convergence entry above: it concluded from levels 4-5 that asr_ed is converged to ~6 mRy in depth.  Levels 6 and 7 show otherwise for the DOS peaks.
CuZn x = 0.5, a = 2.85 A, rank 1, asr_ed (Zn d / Cu d peaks vs E_F, Ry):
| levels | E_F | qnet Cu | Cu C_d | Zn C_d | Zn d peak (height) | Cu d peak | wall |
|---|---|---|---|---|---|---|---|
| 4 | 0.0145 | 0.0444 | -0.2544 | -0.6769 | -0.728 (47) | -0.278 | 14.5 min |
| 5 | 0.0170 | 0.0367 | -0.2566 | -0.6705 | -0.733 (36) | -0.287 | 35 min |
| 6 | 0.0149 | 0.0412 | -0.2389 | -0.6771 | -0.691 (42) | -0.267 | ~1 h |
| 7 | 0.0014 | 0.0386 | -0.2464 | -0.6732 | -0.675 (44) | -0.241 | 1 h 29 min |
The band centres C_d (Zn -0.673 +- 0.004) and qnet (0.037-0.044) are stable, but E_F still moves by 13 mRy at level 7 and the DOS peak positions drift by 40-60 mRy
(Zn d -0.733 -> -0.675, towards the ordered -0.636): not converged at level 7; the earlier statement "disordered Zn d is 0.1 Ry deeper than ordered" is therefore
only ~0.04 Ry at level 7 and still moving.  Level 8 (~3.5 h) not run.
Concentration series, asr_ed level 5 rank 1 (NOT depth converged, peaks uncertain by ~50 mRy at this depth):
| x_Cu | E_F | qnet Cu | Cu C_d | Zn C_d | Zn d peak | Cu d peak |
|---|---|---|---|---|---|---|
| 0.24 | 0.0983 | 0.0589 | -0.2106 | -0.6391 | -0.724 | -0.322 |
| 0.50 | 0.0170 | 0.0367 | -0.2566 | -0.6705 | -0.733 | -0.287 |
| 0.75 | -0.0536 | 0.0171 | -0.2973 | -0.6961 | -0.668 | -0.266 |
The band centres C_d move down with Cu content (Zn C_d -0.639 -> -0.696), E_F falls (rigid-band direction), qnet Cu rises with Zn content; peak positions are not monotonic at this
depth, so no trend statement against Banerjee-Mookerjee 2010 yet.
Magnetic binary validation, asr_ed (spin 2, zeta0 0.1, a as stated):
| system | levels / rank | result | reference | verdict |
|---|---|---|---|---|
| bcc Fe (A = B), a = 2.8665 | 12 / 0 | 2.199 muB, E_F -0.0907 | k-space ASA 2.226, -0.0926 | -1.2 %: pass |
| fcc Ni (A = B), a = 3.524 | 24 / 0 | 0.554 muB (sign = seed orientation), E_F -0.1069 | k-space ASA 0.632, -0.1014 | -12 %: needs depth >= 30 (as EI): inconclusive |
| bcc FeCo 50:50, a = 2.857 | 4 / 1 | Fe 2.624, Co 1.584 (cell mean 2.104), qnet Fe -0.019 | ordered B2: Fe 2.76 / Co 1.77 (ASA), Elk 2.75 / 1.75 | Fe -5 %, Co -10 % at depth 4: depth study needed |
| bcc FeCo, EI order 2, same input | 4 / 1 | collapsed to the NON-magnetic solution (moments 0.0004) | | EI cross-check failed (seed/stability), no comparison available |
| bcc FeV 50:50, a = 2.95 (Vegard guess, not the paper's) | 4 / 1 | Fe -2.682, V -0.080 (V parallel to Fe) | Chakrabarti-Mookerjee 2005: Fe 1.12, V -0.11 (antiparallel) | FAILS: Fe far too large, V sign wrong |
| bcc FeV, EI order 2 | 4 / 1 | non-magnetic collapse | | cross-check failed |
Status: magnetic binary systems are NOT yet validated beyond pure Fe.  Open: depth study of FeCo (levels 5-7), FeV (lattice constant scan 2.85-3.0 A, spin seeds, depth), Ni at 30-40 levels,
a working magnetic EI cross-check (larger seed / order 3).

## q4 (2026-10-11, desktop trial/q4, same binary as q2/q3): CuZn level 8, Ni depth, FeCo depth, FeV scan, magnetic EI cross-checks
- CuZn x = 0.5, asr_ed level 8: E_F 0.0044, qnet Cu 0.040, Cu C_d -0.243, Zn C_d -0.677, Zn d peak -0.680, Cu d peak -0.232 (level 7: 0.0014 / -0.675 / -0.241).
  Levels 7 -> 8 change E_F by 3 mRy and the peaks by 5-9 mRy: approaching convergence.  Disordered Zn d peak -0.68 vs ordered -0.636 (same a): ~0.04 Ry.
- fcc Ni (A = B) asr_ed level 36: 0.609 muB (level 24: 0.554) vs k-space ASA 0.632: -3.7 %, converging with depth as on the EI route.
- bcc FeCo 50:50 asr_ed, rank 1: level 4 Fe 2.624 / Co 1.584; level 5 2.151 / 1.297 (E_F -0.105); level 6 2.493 / 1.730 (E_F -0.079).  NOT converged, oscillates by
  ~0.4 muB with depth.  EI order 2 level 4 with seed 0.25: Fe 2.049 / Co 1.573 (seed 0.1 collapsed to non-magnetic).  Ordered B2 reference 2.76 / 1.77.
- bcc FeV 50:50 asr_ed level 4: Fe 2.477 / 2.588 / 2.683 / 2.764 and V +0.008 / +0.036 / +0.080 / +0.146 at a = 2.85 / 2.90 / 2.95 / 3.00 A; seed 0.25 at 2.95 gives the same
  solution (2.683 / 0.080).  EI order 2 stays non-magnetic with seed 0.25.  Paper (Chakrabarti-Mookerjee 2005): Fe 1.12, V -0.11.  The two routes (equal at infinite depth)
  disagree qualitatively at level 4 (2.7 vs 0) and neither matches the paper.
Status: pure Fe and Ni pass (Ni with depth 36); disordered magnetic binaries (FeCo, FeV) are NOT validated: depth 4-6 is far from converged for them (moment oscillation),
and FeV contradicts the reference.  Needs investigation in the development thread (terminator / depth for spin-split d bands near E_F, multiple SCF solutions,
EI-vs-ED disagreement at fixed potentials).

## Why disordered magnetic binaries fail at levels 4-6 (2026-10-11, desktop trial/q5-q9; code 9204e4e)
Question: q3/q4 left FeCo (moment oscillating ~0.4 muB with depth) and FeV (Fe 2.5-2.8 vs paper 1.12, V parallel) unexplained.  Tests, each isolating one cause:
1. Pure bcc Fe (A = B, no disorder; a = 2.8665, supercell 16).  In this limit the asr_ed operator IS H^(2) (gate test_asr_ed_sym 5e-14), so asr_ed vs EI differ only
   in the tail and contour.  k-space ASA 2.226 muB:
| levels | 4 | 5 | 6 | 8 | 12 |
|---|---|---|---|---|---|
| asr_ed, hull tail (default) | 2.64 | 2.01 | 2.27 | 2.13 | 2.20 |
| asr_ed, `ed_terminator = "bp"` | 2.44 | 2.15 | 2.33 | 2.16 | |
| asr_ed, `ed_terminator = "avg2"` | 2.29 | runaway (-4.78, E_F +0.67) | 2.32 | 2.16 | |
| EI (avg2 tail) | 2.31 | 2.21 | 2.34 | 2.19 | |
   => The FeCo oscillation (Fe 2.62 / 2.15 / 2.49 at levels 4/5/6) is the pure-Fe behaviour of the hull tail (2.64 / 2.01 / 2.27), not disorder.  Hull = band edges
   at the extreme Ritz values of the m x m chain, well inside the band at m = 4-6 -> too narrow -> too magnetic.  New option `ed_terminator` (asr_endep%ed_tail):
   'bp' = complex Beer-Pettifor fixed point (edges where the terminated fraction has its poles; == bp_asymptote on real chains, test_terminators) halves the
   spread; 'avg2' tracks EI but ran away once (complex avg2 is not guaranteed analytic).  Default still 'hull'.  bp costs ~15 small zgeev per fraction: with
   threaded MKL this dominated tiny runs (pure Fe 50x slower); run with MKL_NUM_THREADS=1.
2. Ordered B2 FeV, our k-space LMTO-ASA (independent of the recursion), a = 2.85 / 2.91 / 2.95 / 3.00 A: Fe 1.007 / 1.201 / 1.353 / 1.575, V -0.112 / -0.292 /
   -0.434 / -0.644 (antiparallel), qnet(Fe) +0.331 / 0.333 / 0.332 / 0.329, E_F 0.029 / -0.018 / -0.047 / -0.079.  (Moruzzi-Marcus PRB 45, 2934 (1992), FLAPW
   ordered FeV: Fe ~0.66, V ~0, alignment volume dependent.)  The Chakrabarti-Mookerjee 2005 disordered numbers (Fe 1.12, V -0.11, 0.31 e V -> Fe) are close to
   our ordered B2 values; that paper gives neither the lattice constant nor the recursion depth (linear mixing, 250 iterations) -> weak gate.
3. Ordered B2 FeV through the real-space EI recursion (same code path as the alloys, no disorder), a = 2.95: levels 4 / 6 / 8 / 12 / 24 -> Fe 1.85 / 0.98 / 0 / 0 /
   1.42, V -0.46 / -0.38 / 0 / 0 / 0.62 (opposite sign = antiparallel), qnet 0.336 at every depth (= k-space).  FeV is a marginal magnet (like Ni): its moment
   needs >= 24 levels even without disorder; at 8-12 levels the magnetic solution is lost.
4. SIM strength (disordered FeV, asr_ed bp level 4, a = 2.95): sim_beta 0.7 -> qnet 0.04; 0 -> qnet 0.084, Fe -2.47, V +0.18 (antiparallel).  The SIM is not what
   keeps Fe at ~2.5.  But the SIM coefficient is suspicious: V = beta 2 alpha_M q/S with alpha_M = 1.79 gives 1.25 x 2q/S, stronger than the exact B2 shift
   (2 (1.7627/d_NN) q = 1.00 x 2q/S; ordered transfer 0.33 e).  Ruban-Skriver write V = -alpha e^2 q/S with alpha ~ 0.6-0.8 (no alpha_M factor) -> our SIM
   over-screens by ~1.8x (disordered transfer 0.04 vs ordered 0.33 in FeV; 0.04 vs 0.137 in CuZn).  Not yet changed (affects CuZn results; needs a supercell gate).
Conclusions: (i) levels 4-6 cannot give moments of disordered magnetic alloys: Fe-type needs ~8-12 levels with a good tail, FeV/Ni-type >= 24; exact rank-1
generation stops at ~8 (2.4x nodes per level).  (ii) The hull tail adds a +-0.3 muB depth oscillation at 4-6 levels; bp is the better ED tail.  (iii) The FeV
reference is weak; our own ordered B2 k-space values and SQS supercells (k-space code, multi-atom cells) are the right gates.  (iv) Every magnetic SCF passes through
wild early iterations (E_F +-0.4 Ry, qnet +-1.5 e, moment sign flips, also for pure Fe at 12 levels): which magnetic solution is reached is not controlled.
Legacy note: AMscf and ASR_CORRECT hard-code `read(11)(map(i,j),j=1,15)`; the BCC map has 17 integers per record (FCC 15), so with PARAM_BCC the flip flag (col 16)
and the real-space site index used as the central-site test (col 17) are never read: bcc runs never freeze the central site and drop the (y - x) diagonal terms
(ap9/ap13/p3; irrelevant at x = 0.5).  AMscf's map is exact to 6 augmented-space hops (73109 nodes) and it runs 12 recursion steps on it.

## How AMscf differs (re-checked 2026-10-11) and the hard-wall recursion (`recursion_steps`, commit 0973dd0; desktop trial/q10)
AMscf (ASR/main.f90, hop_endep.f90, newfit.f90, cgreen.f90, pardos.f90; files/BCC_MAP decoded):
- Augmented space = a FIXED map: 1105 real-space sites (6 hops of the 1st+2nd-shell bcc graph), configurations with <= 3 flipped sites, 73109 nodes;
  links leaving the map are dropped (hard wall).  On it AMscf runs 12 Lanczos steps (maxrec), real-symmetric, at 8 real seed energies (-0.9 ... 3.0 Ry),
  spline-interpolates a_n(E), b_n(E) to 900 real energies, Beer-Pettifor tail (eigen of the modified chain, a_inf fixed point), real-axis DOS and E_F.
  Its map is SMALLER than our level-4 rank-1 space (6 hops vs 8; each H^(2) step moves 2 hops), but it does not stop the chain at the edge of the exact
  space: the extra steps see the finite cluster.  Our code tied the Lanczos depth to the generation depth (`levels`), so it terminated after 4-6 steps.
- Madelung: vm = ruban * x q / w with ruban = 0.10 (pardos.f90), i.e. essentially no SIM (ours: beta 2 alpha_M q/S = 2.5 q/S).
- bcc reading bug (15 of 17 map columns, see above).
New option `[asr] recursion_steps` (default = levels): Lanczos steps run on the node set generated to `levels` (links outside the set are dropped by
build_augsym_operator, i.e. a hard wall as in AMscf).  Pure Fe (A = B), levels 4 (215 reduced nodes), asr_ed, k-space ASA 2.226 muB:
| recursion steps | 4 | 12 | 16 | 24 |
|---|---|---|---|---|
| hull tail | 2.643 | 2.200 | 2.299 | 2.246 |
| bp tail | 2.438 | 2.234 | 2.288 | 2.249 |
=> generating the augmented space to a modest depth and running many steps on it (AMscf's recipe) removes the depth swing for Fe (levels-only: 2.64 / 2.01 /
2.27 / 2.13 at 4/5/6/8).
Disordered, asr_ed, levels 4, rank 1 (17515 nodes), recursion_steps 12, bp tail, SIM beta 0.7 (q10, ~5.5 min per SCF iteration with two runs sharing 28 cores):
| system | result (Fe positive) | E_F | qnet Fe | before (levels 4 only, hull) | reference |
|---|---|---|---|---|---|
| FeCo 50:50, a = 2.857 | Fe 2.512 / Co 1.767 (parallel) | -0.0806 | -0.013 | 2.624 / 1.584 | ordered B2 2.76 / 1.77 |
| FeV 50:50, a = 2.95 | Fe 1.726 / V -0.361 (antiparallel) | -0.0668 | +0.051 | 2.682 / +0.080 | ordered B2 k-space 1.353 / -0.434 (qnet 0.332); paper 1.12 / -0.11 |
The hard-wall recursion fixes both qualitatively (FeV: V antiparallel, Fe moment down by 1 muB).  Open: convergence in recursion_steps (16/24), the SIM
(FeV transfer 0.05 e vs 0.33 ordered; more charge on Fe lowers its moment), and an independent disordered reference (SQS).
Also: FeCo levels 4 (no hard-wall extension) with bp tail: Fe 2.426 / Co 1.404, E_F -0.0769 (hull: 2.624 / 1.584).

## FeRh vs SPR-KKR-CPA (2026-10-12, desktop trial/q11)
Reference: Greeshma & Banerjee, Pramana 97, 95 (2023) (references/FeRh_SPRKKR.pdf): SPR-KKR, full potential, PBE, lmax 2; Table 1 at the optimised a
(fig. 1b minimum ~5.78, axis labelled A but read as Bohr -> 3.059 A): ordered B2 Fe 3.271 / Rh 0.942, disordered A2 (CPA) Fe 3.070 / Rh 0.725.
Ours (LDA vbh, ASA equal spheres, a = 3.059 A):
| | Fe | Rh | qnet Fe | E_F |
|---|---|---|---|---|
| ordered B2, k-space LMTO-ASA | 3.229 | 1.064 | +0.126 | -0.1132 |
| disordered A2, asr_ed levels 4 rank 1, recursion_steps 12, bp, SIM 0.7 | 2.867 | 0.751 | +0.022 | -0.1272 |
| disordered A2, same, recursion_steps 24 (110 its) | 2.901 | 0.710 | +0.021 | -0.1237 |
Ordered: within 0.04 (Fe) / 0.12 (Rh) of SPR-KKR despite LDA vs PBE.  Disordered: Rh matches CPA (0.75 vs 0.73); Fe 0.20 low (ordered -> disordered
drop 0.36 vs CPA 0.20).  12 -> 24 steps: Fe +0.03, Rh -0.04, i.e. converged in steps; the remaining Fe deficit (0.17) is not recursion depth.
Candidates: SIM (Fe gains 0.02 e vs 0.13 ordered; SPR-KKR-CPA has no SIM), the hard-wall node set (levels 4), LDA vs PBE, or a genuine ASR-vs-CPA difference.
SIM strength and functional (2026-10-13/14, desktop trial/q12, q13; disordered A2, levels 4, recursion_steps 12, bp; PBE added in 7c64127):
| run | Fe | Rh | qnet Fe | E_F |
|---|---|---|---|---|
| LDA, sim_beta 0.7 | 2.867 | 0.751 | 0.022 | -0.1272 |
| LDA, sim_beta 0.39 (mixing 0.1) | 2.863 | 0.755 | 0.028 | -0.1270 |
| LDA, sim_beta 0.2 (mixing 0.1) | 2.859 | 0.758 | 0.033 | -0.1271 |
| LDA, sim_beta 0 (mixing 0.2) | SCF fell onto the non-magnetic solution (early swings) | | 0.058 | -0.1474 |
| PBE, sim_beta 0.7 (mixing 0.1) | 2.951 | 0.905 | 0.020 | -0.0575 |
| ordered B2 k-space, PBE | 3.305 | 1.068 | 0.107 | -0.0397 |
SPR-KKR (PBE): ordered 3.271 / 0.942, CPA 3.070 / 0.725.
=> The SIM does not matter (moments within 0.01 for beta 0.2-0.7; the transfer stays 0.02-0.03 e: small transfer is the disordered environment, not the SIM).
   [WITHDRAWN 2026-10-14: these runs used the wrong-sign SIM (see "SIM sign and magnitude error"); with beta 0 the vacuum-boundary self term +2q/S alone
   still over-screens, so the beta scan did not test physical screening.]
=> PBE vs PBE: Rh +0.13 (ordered) / +0.18 (disordered) high in both: ASA (equal spheres) offset.  Fe: ordered +0.03, disordered -0.12; the
   ordered -> disordered Fe drop is 0.35 (ours, LDA and PBE alike) vs 0.20 (CPA).  Not depth, node set, SIM or functional: either ASR (local
   environments) vs single-site CPA, or ASA vs full potential.  Decisive check: SQS supercell with our own k-space ASA (same ASA and functional).

## MnCr vs Banerjee & Mookerjee arXiv:0812.0216 Fig. 1 (2026-10-14, desktop trial/q14)
Reference (data/reference/mncr/banerjee2009_fig1_50-50.dat, oracle only): disordered bcc Mn50Cr50, legacy TB-LMTO-ASR (VBH, a = 5.4802 Bohr = 2.900 A,
12 recursion steps on the legacy bcc map); peaks (E - E_F, Ry): majority -0.360 / -0.227 / +0.022, minority -0.377 / -0.239 / +0.014, N(E_F) 15.6 / 18.4
per spin; nearly non-magnetic.  Ours, same a and xc, asr_ed levels 4, bp:
| run | moments Mn / Cr | E_F | peaks (E - E_F) | N(E_F) per spin |
|---|---|---|---|---|
| disordered, recursion_steps 12 | 0.003 / -0.001 | -0.0349 | -0.163, +0.072 | 8.5 |
| disordered, recursion_steps 24 | -0.986 / +0.181 | -0.0340 | spin 1 -0.128, +0.074; spin 2 -0.161, +0.009, +0.087 | 6.0 / 13.7 |
| ordered B2, k-space, non-magnetic | - | -0.0269 | -0.271, -0.178 / -0.146, +0.032 (0.02 Ry Gaussian) | 14.8 |
| ordered B2, k-space, spin | 0.072 / 0.429 | -0.0267 | | |
| pure bcc Cr, k-space | - | -0.0217 | -0.248, -0.158 / -0.111, +0.087 | 4.8 |
| pure bcc Cr, asr_ed route (A = B), 12 / 24 / 48 steps | - | -0.023 / -0.017 / -0.015 | 48: -0.253, -0.183 / -0.113, +0.067 / +0.105 | 5.6 / 4.8 / 5.0 |
Findings: (i) the recursion route reproduces the bcc d-band shape of pure Cr vs our k-space LMTO (all features at 48 steps, N(E_F) within 0.2-0.8);
(ii) Fig. 1 is NOT reproduced: its bonding peaks (-0.36 / -0.23) are 0.08-0.09 Ry deeper than in both our ordered B2 MnCr k-space (-0.27 / -0.18) and
pure Cr k-space (-0.25 / -0.16) at the same a and xc, so the 2009 band is wider than any of our independent calculations; candidate legacy causes: the bcc
map-reading bug (15 of 17 columns: central site not frozen), legacy CTRL sphere radii (Mn R = 2.14 vs Cr 2.70 Bohr, not space filling), 14-neighbour S.
(iii) Our disordered DOS has N(E_F) 8.5 (12 steps) vs 14.8 in ordered B2 and the top peak 40 mRy higher than B2; the disordered magnetic state is not
converged in steps (non-magnetic at 12, Mn -0.99 / Cr +0.18 at 24): MnCr is a marginal magnet like FeV.  Open: 48-step disordered run, antiparallel seed.

## SIM sign and magnitude error (found 2026-10-14 in the legacy/reference audit; fixed in 326759b)
The sphere Hartree potential is solved with the vacuum boundary condition (asr_atom%poisson_hartree: v_H(S) = 2 N/S, nucleus -2Z/r), so at S the
sphere already carries its own +2 qnet/S (qnet = N - Z), exactly as in the k-space ASA (asr_asa: v = vacuum Hartree + 2 sum_j Phi_ij qnet_j; for B2
2 (Phi_AA - Phi_AB) = -2 x 1.00 / S, i.e. the external term is NEGATIVE on the electron-rich site).  The SIM (Ruban & Skriver PRB 66, 024201 (2002); thesis
sec. 2.6.6: V_i = -e^2 Q_i / R_1, screening charge on the first shell) must therefore be V = -beta 2 qnet / R1 (bcc: R1 = S / 0.569).
The code had V = +beta 2 alpha_M qnet / S (alpha_M = 1.79, beta 0.7 -> +2.5 q/S).  History: the first version had the minus sign but also the alpha_M factor,
-2.5 q/S over-cancels the self term (net -0.5 q/S: runaway); on 2026-10-02 the sign was flipped instead of removing alpha_M.
Net restoring term at S per unit transfer (bcc): old +4.5 q/S; corrected (beta 1) +0.86 q/S; ordered B2 exact Ewald ~0; AMscf (pardos.f90 replaces the whole
ves by 0.1 x q/S) +0.05 q/S.  => every disordered result before 326759b had its charge transfer suppressed (CuZn 0.03 vs ordered 0.137, FeRh 0.02 vs 0.13,
FeV 0.05 vs 0.33), and may have stiffened the early SCF swings.  Ordered-limit runs (exact Ewald, `vmo`) and all k-space results are unaffected.
Reruns with beta 1 (desktop trial/q15): FeRh r12 bp, FeV r12 bp, CuZn x50 l4 r1, MnCr r24 bp.  The 48-step MnCr run (q14) was stopped at it 29
(E_F drifted to +0.19 Ry, Cr 2.2 muB; adaptive depth 12-48 per energy).
First q15 result, CuZn x = 0.5 (asr_ed levels 4 rank 1, same input as trial/amscf/x50_l4r1 except the SIM; 22 its vs 39):
| SIM | qnet Cu | E_F | Zn d peak vs E_F | Cu d peak | DOS(E_F) |
|---|---|---|---|---|---|
| old (+1.25 x 2q/S) | +0.044 | 0.0145 | -0.728 | -0.278 | 1.74 |
| corrected (beta 1) | +0.118 | 0.0130 | -0.709 | -0.291 | 1.67 |
Ordered B2 qnet 0.137; Banerjee-Mookerjee 2010 (50:50) Zn d -0.682, Cu d -0.309: both peaks move toward the reference (Zn 46 -> 27 mRy, Cu 31 -> 18 mRy off).
Magnetic q15 reruns (corrected SIM, otherwise the q10/q11 inputs: asr_ed levels 4 rank 1, recursion_steps 12, bp):
| system | SIM | Fe | Rh / V | qnet Fe | E_F | its |
|---|---|---|---|---|---|---|
| FeRh A2, a = 3.059 | old | 2.867 | +0.751 | +0.022 | -0.1272 | 70 |
| FeRh A2 | corrected | 2.833 | +0.781 | +0.067 | -0.1267 | 63 |
| FeV 50:50, a = 2.95 | old | 1.726 | -0.361 | +0.051 | -0.0668 | 51 |
| FeV 50:50 | corrected | 1.647 | -0.351 | +0.157 | -0.0650 | 40 |
References: FeRh SPR-KKR-CPA 3.070 / 0.725, ordered B2 (ours) 3.229 / 1.064, qnet 0.126; FeV ordered B2 k-space 1.353 / -0.434, qnet 0.332.
=> transfers x3 toward the ordered values; moments change little (FeRh Fe -0.03: the ordered -> disordered Fe drop, 0.40 vs CPA 0.20, is NOT the SIM).
FeRh passed through an antiparallel-Rh state (Rh -0.74 at it 26) before settling on the ferromagnetic solution: the start still decides the path.
MnCr r24 (corrected SIM): not converging (it 74: dV 0.3, Mn 1.41 / Cr -0.57, E_F -0.090; at it 25 qnet(Mn) + qnet(Cr) = -0.14, i.e. species charges
did not balance); same pathology as the stopped 48-step run.  Suspect: complex-symmetric chains beyond 12 steps with adaptive depth (12-24 per energy) on MnCr.

## Our LMTO-ASA vs Questaal `lm` (2026-10-15; Questaal ASA built from ../questall, gcc+MKL, no libxc; scratch runs, oracle only)
Same model in both: a = 5.48021 Bohr, equal spheres S = 2.6983, spd, hard-core 0.7 S screening, scalar relativistic, no combined correction,
three-centre ASA (ours O = 1 + (o_i+o_j) h + h p h; Questaal TWOC=0), tetrahedra.  Energies relative to E_F (the codes' zeros differ).
- B2 MnCr, non-magnetic (Questaal XCFUN=2): qnet Mn +0.113391 vs ours +0.11339; d parameters C-E_F within 0.3 mRy (Mn) / 0.1 mRy (Cr),
  sqrt(Delta), gamma, p within 0.2%; s, p C within 1.4 mRy; q2 identical; DOS (both broadened 0.005 Ry) peaks within 2-6 mRy, N(E_F) 59 vs 64
  states/Ry/cell (E_F on a steep peak).
- bcc Cr (non-magnetic): C_d - E_F -0.0114 vs -0.0137, E_nu,d - E_F -0.161 vs -0.164, Delta and gamma identical to 0.1%.
- XC: Questaal XCFUN=2 "Barth-Hedin" uses the ASW/Hedin-Lundqvist constants (c_P/c_F = 0.045/0.0225, r_P/r_F = 21/52.92); libxc VBH (ours) and
  the legacy AMscf LMTO (lmto/ATOM/xcloc1.f, cp = -0.0504) use the original von Barth-Hedin ones (0.0504/0.0254, 30/75) -> compare moments with VWN.
- bcc Fe (VWN, a = 2.8665 A, 16^3): ours 2.252, Questaal 2.275 muB (1%).
- B2 MnCr spin-polarised (VWN): Questaal Mn 0.033 / Cr 0.612 (16^3), 0.033 / 0.585 (24^3), 0.027 / 0.566 (32^3) (non-magnetic state unstable: the
  moment grows smoothly from 0.24); ours 0.035 / 0.531 (16^3, overall sign arbitrary), NON-magnetic at 24^3 and 32^3 (moment decays from 0.15).
  Questaal's structure-constant range (default 2.7 avw, 15-site clusters) is not the cause: with our 4 avw (59 sites) Cr 0.610.
=> Our LMTO-ASA reproduces Questaal for the band structure and charge transfer; magnetism is ~1% weaker (Fe), which tips the marginal B2 MnCr
   (Cr ~0.6 muB in Questaal) to non-magnetic at converged k.  Open: source of the 1% (E_nu choice ~3 mRy lower in Questaal, k mesh shifted vs
   Gamma-centred, Bloechl correction, core polarisation).

## MnCr structure and magnetic states (2026-10-15; Questaal built on the desktop ~/questaal_build, desktop trial/q17, q18)
Materials Project 50:50 MnCr (Cmmm, a = b = 2.518, c = 3.936 A, gamma = 109.1 deg; PBE, hypothetical ordered cell): a bcc superstructure, 8 + 6
neighbours at 2.45-2.52 / 2.84-2.92 A (ideal bcc at the same volume 2.483 / 2.868), 11.79 A^3/atom (a_bcc 2.868); its NN shell is 4 unlike + 4 like
(B2: 8 unlike).  => the bcc parent of the disordered alloy is right; only the ordered B2 reference is not the PBE ground-state order.
Questaal, a = 2.900 A, VWN, 24^3 (ASA equal spheres, same model as ours), energies per two-atom cell:
| state | Mn | Cr | E - E(NM) (mRy) |
|---|---|---|---|
| bcc Cr AF (CsCl cell) | +0.571 | -0.571 | -0.085 (FM start collapses to NM) |
| B2 MnCr NM | 0 | 0 | 0 |
| B2 MnCr weak FM (parallel or high-spin start) | 0.033 | 0.585 | -0.26 |
| B2 MnCr ferrimagnet (antiparallel start +-1.5 or +-3.8) | +1.985 | -1.080 | -2.12 |
Ours (k-space, new [[species]] spin_sign = -1 for Cr, 24^3): VWN zeta0 0.2 -> Mn +1.944 / Cr -1.057 (Questaal 1.985 / -1.080); VBH zeta0 0.2 -> 1.662 / -0.823;
zeta0 0.1 -> weak state (Cr -0.54 VWN / -0.43 VBH); zeta0 0.4 -> RUNAWAY: E_F climbs ~0.15 Ry per iteration with growing charge transfer (1.9 e at
it 10) and ends at a fixed point Mn 4.9 / Cr 3.4 (parallel), E_F +0.21 Ry, which Questaal never reaches from the same start (it converges to the
ferrimagnet).  => (i) the ground state of B2 MnCr is the ferrimagnet, missed by parallel starts (ours, the 2009 paper); (ii) our SCF has a spurious
runaway mode at large initial polarisation -- the same E_F drift (+0.2 Ry) seen in the disordered MnCr / FeRh / FeV runs.  Next: find it.

## SCF runaway found and fixed (2026-10-15; desktop trial/q17-q20b; commits 7bc52d2, 460f2f2, 4ba5ca0, ca48346)
Diagnosed with the new `[scf] trace = true` (k-space ASA: E_nu, C, sqrt(Delta)_d, q0, V(S), Madelung per iteration/class/spin; V_in/V_out at S).
1. Mixer (the cause of the early E_F excursions in every magnetic SCF, also pure Fe): unweighted Anderson on the log mesh, whose many points near
   the nucleus dominate the least-squares fit; near S the step went against the residual every iteration (bcc Fe: V_out(S) = -0.71 fixed while
   V_in(S) climbed by the 0.15 Ry clamp per iteration; E_F to +1.5 Ry by iteration 12, moment flipped, convergence only at iteration ~60-70).
   Fix: Anderson inner products weighted by (r/(1+r))^2 (the dV measure; E_nu entries weight 1) + restart of the history when the step opposes the
   residual or the residual doubles.  (A linear start while dV > 0.1 was tried and dropped: B2 MnCr charge-sloshed at +-1.4 e.)
2. E_nu could leave its band: past phi_l(S) = 0 the parametrisation jumps to the next radial branch (Mn minority d: E_nu 1.06, C 4.86 Ry,
   sqrt(Delta) -1.12) and the SCF converged to a spurious Mn 4.9 / Cr 3.4 muB state with E_F +0.21 Ry.  Fix: band_window / clamp_enu keep each E_nu in
   [E(D_l = l), E(phi_l(S) = 0) - 5%] of the branch with n-l-1 nodes (Questaal's pnu floor / fixed principal quantum number); test_band_window.
   (A first version with the floor at D_l = 0 pinned the p E_nu above the centre of gravity of the occupied p weight: no convergence.)
Results (k-space, 24^3 unless stated; max E_F reached during the SCF in brackets):
| run | before | after |
|---|---|---|
| bcc Fe VWN / VBH | 2.252 / 2.214, 58-72 its (E_F max +1.5) | 2.252 / 2.214, 48-49 its (E_F max +0.07) |
| B2 MnCr VWN, Cr antiparallel, zeta0 0.4 | spurious Mn 4.9 / Cr 3.4, E_F +0.21 | ferrimagnet 1.944 / -1.057, 57 its (E_F max 0.17) |
| B2 MnCr VWN, zeta0 0.2 antiparallel / parallel start | 1.944 / -1.057 (66 its) / weak state | 1.944 / -1.057 in 40 / 56 its (Questaal 1.985 / -1.080) |
| B2 MnCr VBH, antiparallel (zeta0 0.2 or 0.4) | 1.662 / -0.823 (111-136 its) or spurious | 1.662 / -0.823 in 41-57 its |
| B2 MnCr VBH, parallel start | weak state | not converged in 150 its (weak states, dV 1e-4 level) |
Earlier magnetic alloy runs (q5-q16) used the old mixer; their wild early iterations and the E_F +0.2 Ry drifts are this bug.  Converged results
are fixed points of the same equations and stand, but which solution was reached was decided by the broken mixer.

## Adaptive depth and recursion depth at fixed potential; E_F search fix (2026-10-15; desktop trial/q16, commit f6ebe11)
Disordered bcc Mn50Cr50 (VBH, a = 2.900, asr_ed levels 4, bp, SIM beta 1), non-magnetic SCF at 12 steps (44 its, qnet Mn +0.054; old SIM 0.018),
potential saved; then ONE iteration from it at 12 / 24 / 48 steps with depth_tol 0 and 1e-5:
| steps | E_F (tol 0) | E_F (tol 1e-5) | qnet Mn / Cr | dV vs input |
|---|---|---|---|---|
| 12 | -0.035031 | -0.035031 | +0.05396 / -0.05394 | 5.7e-5 (restart check) |
| 24 | -0.029774 | -0.029774 | +0.10999 / -0.11000 | 3.8e-2 |
| 48 | -0.029176 | -0.029176 | +0.10857 / -0.10911 | 3.9e-2 |
=> the adaptive depth (depth_tol) changes nothing at any depth; charges balance.  An external review's claim that depth_tol breaks the contour
integral (and caused the unbalanced species charges) is refuted.  MnCr needs 24 steps (12 -> 24: charge transfer x2, E_F 5 mRy; 24 -> 48: 0.6 mRy).
The unbalanced charges (MnCr -0.14 e, CuZn +0.34 e in early iterations) came from the E_F search: 6 contour passes, then a first-order moment
correction with a slope different from the DOS used in it.  Fix f6ebe11: 12 passes, correction with DOS(E_F), remaining N - N_elec reported.

## Disordered MnCr with all fixes, ferrimagnetic start (2026-10-15; desktop trial/q23; code 17416b2)
bcc Mn50Cr50, VBH, a = 2.900 A, asr_ed levels 4 rank 1, bp, SIM beta 1, weighted Anderson, E_nu window, bracketed E_F search; start Cr spin_sign -1,
zeta0 0.2; potentials saved (pot.dat).
| steps | its | Mn | Cr | qnet Mn | E_F | N(E_F) spin 1 / 2 | peaks vs E_F, spin 1 ; spin 2 |
|---|---|---|---|---|---|---|---|
| 12 | 66 | 0.000 | 0.000 | +0.054 | -0.0350 | 8.5 / 8.5 | -0.162, +0.075 |
| 24 | 58 | +0.969 | -0.178 | +0.055 | -0.0341 | 13.6 / 6.0 | -0.163, +0.009, +0.087 ; -0.131, +0.074 |
Both converged cleanly (one E_F-search warning each, in the first iterations; DOS positive everywhere).  24 steps is the converged depth (q16).
The 24-step state is the same as the old q14 r24 run (Mn -0.986 / Cr +0.181, reached through wild iterations): disordered MnCr is a weak
ferrimagnet, ~0.40 muB per atom net, far below ordered B2 (VBH 1.66 / -0.82).  12 steps loses the moment even from an antiparallel start.
vs Banerjee & Mookerjee 2009 Fig. 1 (nearly non-magnetic, N(E_F) 15.6 / 18.4, bonding peaks -0.36 / -0.23): still not reproduced; the legacy band
is 0.1-0.2 Ry wider than ours, which agrees with Questaal for the ordered cell (see "Our LMTO-ASA vs Questaal").

## Clean reruns of all disordered binaries, VBH and PBE (2026-10-15/16; desktop trial/q24, code a28a0d1)
All fixes in (corrected SIM beta 1, weighted Anderson, E_nu window, bracketed E_F search); asr_ed levels 4 rank 1, 24 recursion steps, bp tail;
starts: FeRh / FeCo parallel (zeta0 0.1), FeV / MnCr antiparallel (spin_sign -1 on V / Cr, zeta0 0.2).  Moments A / B (muB), qnet of A.
| system (a, A) | VBH | PBE | reference |
|---|---|---|---|
| FeCo (2.857) | 2.578 / 1.746, qnet -0.040, 45 its | 2.665 / 1.793, -0.041, 53 its | ordered B2 (ours) 2.76 / 1.77 |
| FeV (2.95) | 1.848 / -0.469, +0.157, 43 its | 2.08 / -0.65, +0.15 (stationary from it ~18; dV floor 1e-3, stopped at 101) | ordered B2 k-space (VBH) 1.353 / -0.434; Chakrabarti-Mookerjee 2005 1.12 / -0.11 |
| FeRh (3.059) | 2.867 / 0.742, +0.063 (stationary from it ~30; dV floor 1e-4, stopped at 103) | 2.956 / 0.991, +0.058, 39 its | SPR-KKR-CPA (PBE) 3.070 / 0.725 |
| CuZn x = 0.5 (2.85) | NM, qnet Cu +0.110; Zn d -0.652, Cu d -0.222 vs E_F | +0.110; -0.650 / -0.220 | Banerjee-Mookerjee 2010: -0.682 / -0.309 |
| MnCr (2.900), 24 steps | 0.969 / -0.178, +0.055 (q23) | 1.395 / -0.382, +0.054, 58 its | 2009 Fig. 1 nearly non-magnetic |
| MnCr, 48 steps (restart) | 1.132 / -0.205, 34 its | 1.497 / -0.400, 21 its | |
Findings: (i) every run converged or became stationary; max E_F during the SCF stayed physical (FeCo overshot to +0.31 Ry at it 7 from the
atomic start, then converged).  (ii) Depth: FeV 12 -> 24 steps raises Fe 1.65 -> 1.85 and V -0.35 -> -0.47; MnCr moments are NOT converged at 24
steps (24 -> 48: +17% VBH, +7% PBE) although its non-magnetic E_F / charges are (q16); FeRh 12 -> 24 changes < 0.04.  (iii) CuZn at 24 steps: Cu d
moves 70 mRy towards E_F vs the depth-4 runs (-0.291 -> -0.222), close to ordered B2 (-0.26), away from the 2010 figure; VBH vs PBE identical.
(iv) PBE vs VBH: moments +0.05..0.4 muB; FeRh PBE Fe 2.96 (CPA 3.07), Rh 0.99 (CPA 0.73; ordered PBE ours 1.07 vs SPR-KKR 0.94: equal-sphere ASA).
(v) Two runs (FeRh VBH, FeV PBE) stalled at a dV noise floor: N(E) steps by ~0.01 e at E_F and the sharp-E_F search cannot resolve it ->
Fermi-smeared contour ([asr] fermi_kt, cc2d513/b2367c2); test reruns with kT = 3 mRy in trial/q25.

## Fermi smearing, depth 48, SQS references, first fcc alloys (2026-10-16/17; desktop trial/q25-q28, code faf2c78 / cbd4250)
**Fermi-smeared contour** ([asr] fermi_kt): both runs that stalled with a sharp E_F now converge.
| run | kT (mRy) | its | m A / m B | note |
|---|---|---|---|---|
| FeRh VBH r24 | 3 | 42 | 2.866 / 0.729 | sharp: 2.867 / 0.742, stalled dV 6e-5 |
| FeV PBE r24 | 3 | 39 | 2.068 / -0.621 | sharp: 2.08 / -0.67, stalled dV 7e-3 |
| FeV PBE r24 | 1 | 21 (restart) | 2.074 / -0.639 | dV 1.5e-6 |
The shift scales ~kT^2: kT = 1 mRy is within ~0.003 muB of sharp, kT = 3 mRy pulls V in by ~0.02.  Use kT = 1 mRy (default still 0 in the code).
**Depth 48** (restart from the 24-step potentials, kT 3 mRy): FeRh Fe 2.877 / Rh 0.632 (24: 2.866 / 0.729); FeV VBH Fe 1.890 / V -0.479 (24: 1.848 / -0.469).
**SQS references** (tools/sqs.py, lattice = "cell" + site_classes, our own k-space LMTO-ASA, same spheres / functional, VBH, tetrahedra; species mean
(spread over sites)).  Gate: B2 FeCo as a 2-site cell == lattice "b2" to all digits.  sqs.py fix cbd4250: smallest triplet = least-perimeter triangle
(fcc equilateral R1).  bcc 16-atom cells exact through shell 4 (#1) / 5 (#2), 32-atom essentially random to shell 6; fcc 50:50 exact; fcc 75:25
objective 0.045 / 0.054 (16) / 0.0035 (32).
| alloy (a, A) | ASR r24 | ASR r36/48 | SQS 16 #1 | SQS 16 #2 | SQS 32 | verdict |
|---|---|---|---|---|---|---|
| bcc FeCo (2.857) | 2.578 / 1.746 | - | 2.592 (0.07) / 1.778 (0.01) | 2.594 / 1.785 | not converged (q28 rerun) | agree; Co 0.035 low |
| bcc FeV (2.95) | 1.848 / -0.469 | r48 1.890 / -0.479 | 1.756 (0.47) / -0.364 (0.22) | 1.718 / -0.353 | not converged (q28 rerun) | ASR ~0.1 too magnetic, 48 steps worse |
| bcc FeRh (3.059) | 2.866 / 0.729 | r48 2.877 / 0.632 | 2.977 (0.06) / 0.630 (0.18) | 2.997 / 0.744 | 2.991 (0.05) / 0.709 (0.15) | **Fe 0.11 low at 24 and 48**; Rh agrees |
| fcc Fe0.5Ni0.5 (3.588) | 2.512 / 0.593 | r36 running (it 1: 2.546 / 0.619) | 2.535 (0.11) / 0.662 (0.01) | 2.530 / 0.647 | running | Fe agrees, Ni 9% low |
| fcc Ni0.75Cu0.25 (3.547) | 0.319 / 0.004 | r36 running | 0.517 (0.04) / -0.003 | 0.507 / -0.000 | not converged | **Ni 38% low** |
Charge transfer agrees everywhere (FeRh 0.063 vs 0.064; FeNi 0.086 vs 0.096 e; NiCu ~0).  E_F agrees (NiCu -0.1189 vs -0.1172).
Local environment (SQS): FeCo almost none (Fe 2.50-2.68); FeRh smooth (Rh 0.92 -> 0.43 as like-neighbours go 2 -> 6); FeV strong (Fe 0.87 with 2
Fe neighbours -> 2.27 with 6; V -0.74 -> -0.10).  The like-NN = 4 sites: FeCo 2.600 / 1.787, FeV 1.885 / -0.275, FeRh 3.006 / 0.742.
Reading: the weak / environment-sensitive moments carry the ASR error (Ni in NiCu, Ni in FeNi, Fe in FeRh, FeV); strong Fe / Co moments agree.
Three error sources, separable: (a) augmented-space truncation (levels 4, **rank 1** = at most one flipped site; FeRh Fe does not move with
steps -> rank is the suspect); (b) energy resolution of a 24-step chain whose asymptotic width is the whole spd band (~2 Ry -> ~50-100 mRy,
coarser than the Ni d peak at E_F; pure fcc Ni through ASR 0.554 / 0.609 at 24 / 36 vs 0.632); (c) termination artefacts beyond the exact
region (hard wall -> finite-cluster oscillation of a_n, b_n; the sqrt tail is attached wherever step m lands: MnCr 0.97 -> 1.13, FeRh Rh 0.74 -> 0.63).
**k-space SCF weakness**: 3 of 4 32-atom SQS cells (FeCo, FeV, NiCu; 32 classes, mixing 0.2) reach dV 1e-2..1e-3 and then burst (dV 1-5, E_F
jumps); all 16-atom cells and FeRh 32 converge.  q28 reruns FeCo / FeV 32 with mixing 0.1, 600 its (after q27).

## External audit `suggestion.md` (2026-10-17): checked against the source, mostly wrong
- asr_ed "inconsistent scaling" (C - E_nu unscaled, random hops): false; the quoted code does not exist.  asr_ed has g_R = (Delta_c/Delta_R)(C_R - E_nu,R),
  non-random hops sqrt(Delta_c) S sqrt(Delta_c), o -> o Delta_R/Delta_c: exactly T (E - H^(2)) T.  Its "Option A" (switch to EI H^(2)) is excluded.
- "AMscf / Questaal use spin-averaged E_nu": false (AMscf enutcg loops over isp; Questaal idmod 0 per spin, idmod 4 is the averaged option).
- Shift moments before rho / missing inner atomsc loop / mixpq: path only, same fixed point (rho uses moments about E_nu,in with phi(E_nu,in)).
- levels 4 = "<= 3 flips", bp = "average of the last 4 steps": false (levels = generation depth, rank 1; bp = Beer-Pettifor fixed point).
- 1% vs Questaal: its numbers are wrong (ours 2.252 < 2.275, VWN, no CCOR in either).  Questaal's "vBH" = ASW constants is true but known;
  AMscf xcloc1 uses the original vBH = our libxc vbh.  keep_p exists as density = 'stuttgart' (0.005 muB).  **Real but tiny**: the scalar-
  relativistic factor gfac = 1 + l(l+1)/(tmc r)^2 in Questaal's densities (rhocor.f:298, atomsr.f:2012), absent in ours (~1e-4 muB expected).

## Tail scan at fixed potential (2026-10-18, desktop trial/q29; MAGNETISM_TODO item 1)
`[asr] tail_scan = m_max` (66acbf1): start potential, production E_F search (24 steps, bp), then one chain per contour node to m_max on the contours
of 5 Fermi levels (7.5 mRy apart); every truncation m = 8..m_max x tail (hull, bp, avg2) gets its own E_F (cubic interpolation of N) and moments.
Self-check: (m = 24, bp) reproduces the production E_F and moments to all printed digits in every case.  A = B: same-model k-space reference
(E_nu + h - h o h, same alpha, same S cluster, same potential, nk 40 tetrahedra) = the exact m -> infinity answer of the recursion.
Pure fcc Ni (A = B) uses the Ni potential of q27 nicu_r24 (a = 3.547; not self-consistent for pure Ni, irrelevant for a fixed-potential test).
Moments (bp; hull and avg2 within 0.005 of bp for m >= 26, within 0.01-0.03 below):
| m | Ni A=B l4 | l5 | l6 | l8 | FeRh Fe / Rh (l4) | NiCu Ni (l4) | MnCr Mn / Cr (l4) |
|---|---|---|---|---|---|---|---|
| 12 | 0.291 | 0.291 | 0.291 | 0.291 | 2.825 / 0.760 | 0.259 | 0.804 / -0.191 |
| 24 | 0.379 | 0.384 | 0.391 | 0.387 | 2.866 / 0.729 | 0.319 | 0.969 / -0.178 |
| 32 | 0.436 | 0.427 | 0.431 | 0.425 | 2.896 / 0.704 | 0.350 | 1.025 / -0.159 |
| 40 | 0.461 | 0.457 | 0.458 | 0.456 | 2.882 / 0.691 | 0.387 | 1.025 / -0.160 |
| 48 | 0.447 | 0.473 | 0.471 | 0.476 | 2.876 / 0.685 | 0.395 | 1.013 / -0.168 |
| 64 | - | 0.486 | 0.486 | 0.486 | - | - | - |
| exact / SQS | k-space 0.511 | | | | SQS 2.99 / 0.71-0.74 | SQS 0.51 (SCF) | |
Reading:
- (c) the tail is irrelevant: hull / bp / avg2 agree, their edges converge to the same spd band (Ni: -0.748 .. 1.767 Ry, width 2.5 Ry).  No
  oscillation about a stable mean (only +-0.01 jitter from m to m) -> termination averaging (item 3) and a period-2 terminator (item 4) cannot
  remove the error.  b_n (Ni d) noisy about 0.62 = width / 4; eg alternates only for n < 10.
- (b) energy resolution dominates for Ni: the moment rises with m like m_inf - c/m (Ni l8: 48 -> 0.476, 64 -> 0.486 gives m_inf ~ 0.516 vs exact
  0.511).  24 steps give 74% of the exact Ni moment at fixed potential (SCF amplifies, cf. 0.554 vs 0.632).  NiCu Ni: 0.319 -> 0.395 (24 -> 48).
- levels: identical to m ~ 24; levels 4 (175 nodes for A = B) saturates at m ~ 36 (0.46) and drifts down, levels 5 (506) = 8 (1675) to 0.003
  up to m = 64.  Production for Ni-type moments needs levels 5 and m >= 64.
- FeRh Fe is NOT resolution: flat 2.88 +- 0.01 from m = 26 to 48 (SQS 2.99); Rh drifts 0.73 -> 0.69.  Remaining suspects: rank 1 (item 5) and the
  levels-4 node set.  MnCr: Mn +5% from 24 to 32, then flat 1.01-1.03.

## Contour path independence (2026-10-18, desktop trial/q30; MAGNETISM_TODO item 2) and q27 r36 / q28 SQS 32 results
Iteration-1 output moments at the q29 fixed potentials (24 steps, kT 3 mRy, depth_tol 0); `[asr] fermi_poles` (7343afb) moves the Fermi
contour (Gamma = 2 n pi kT; auto n = 5 -> 0.094 Ry):
| case | auto, 32 pts | auto, 64 | n = 2 (0.038 Ry), 64 | n = 10, 64 | n = 20, 64 |
|---|---|---|---|---|---|
| Ni A=B (bp) | 0.3790 | 0.3790 | 0.3790 | 0.3790 | 0.3790 (hull 0.3802 / 0.3801 at n = 2 / 10) |
| FeRh Fe / Rh, qnet | 2.8657 / 0.7291, 0.0635 | 2.8658 / 0.7302, 0.0632 | 2.8615 / 0.7015, 0.0707 | 2.8658 / 0.7301, 0.0632 | - |
| NiCu Ni | 0.3187 | - | 0.3187 | running | - |
Path independent to 1e-4 (Ni, NiCu) and 1e-3 (FeRh, = the 32/64-point quadrature difference) for Gamma >= 0.09 Ry.  Only the lowest path
(FeRh, Gamma 0.038 Ry) deviates (Rh -0.03, qnet +0.007): quadrature of the near-axis leg or the per-z edge choice; production paths are not
affected.  => the per-z tail edges are not a source of the weak-moment errors; no edge freezing needed.
SCF results that finished meanwhile (q27 r36, q28 32-atom SQS with mixing 0.1):
| alloy | ASR r24 | ASR r36 | SQS (16 / 32) |
|---|---|---|---|
| fcc Ni0.75Cu0.25, Ni | 0.319 | **0.493** | 0.517 / 0.507 (16) |
| fcc Fe0.5Ni0.5, Fe / Ni | 2.512 / 0.593 | 2.555 / **0.641** | 2.535 / 0.662 (16) |
| bcc FeCo, Fe / Co | 2.578 / 1.746 | - | 32: 2.583 (sd 0.06) / 1.780 (sd 0.01) |
| bcc FeV, Fe / V | 1.848 / -0.469 (r48 1.890 / -0.479) | - | 32: 1.808 (sd 0.32) / -0.382 (sd 0.19) |
The Ni moments close most of the gap to SQS between 24 and 36 steps (NiCu -38% -> -4%, FeNi -9% -> -3%): energy resolution of the chain
(item 1) is the cause for Ni; the SCF amplifies the fixed-potential deficit.  FeV stays too magnetic (and gets worse with steps), FeRh Fe flat.
Physics of the resolution limit (see MAGNETISM_TODO): a d start orbital hybridises with s, p, so the chain's asymptotic width is the whole spd
band (Ni: b_1 = 0.22, b_n ~ 0.62 for n >= 2, edges -0.75 .. 1.77 Ry); m steps resolve ~2.5 Ry / m.  A weak moment is set by a d-band edge at E_F
(Ni minority d top); smearing it across E_F moves holes below E_F -> moment systematically low.  Strong moments (both d edges away from E_F) are
immune.  AMscf (12 steps, NN map, real-axis seed energies + spline, BP tail) has no remedy; Ganguly et al. PRB 83, 094407 (2011, 11 steps)
report Mn moments 0.57-0.67 muB below CPA - consistent with this limit (A = B Ni at 12 steps: 57% of the exact moment).

## Average-medium (VCA) correction in alloys (2026-10-18, desktop trial/q31, code 75e5732): negative
SCF from the r24 potentials, `medium = "vca"`, medium_nk 60, kT 1 mRy, 24 steps.  Iteration 1 = fixed-potential effect of the correction alone.
| alloy | r24 plain | it 1 (fixed pot.) | converged | r36 plain | SQS |
|---|---|---|---|---|---|
| NiCu Ni / Cu | 0.319 / 0.004 | 0.3185 / -0.051 | killed it 6: 0.311 / -0.076 (Cu growing) | 0.493 | 0.51 / 0.00 |
| FeNi Fe / Ni | 2.512 / 0.593 | 2.538 / 0.594 | 2.531 / 0.616 (12 its) | 2.555 / 0.641 | 2.535 / 0.662 |
| FeRh Fe / Rh | 2.866 / 0.729 | - | 2.879 / 0.737 (21 its) | - | 2.99 / 0.71-0.74 |
Exact for A = B (gate), useless for alloys: the VCA medium averages the d levels and exchange splittings, so its d edge (and its truncation
error near E_F) is not where an alloy site has it; it adds nothing to the Ni moment and puts spurious spin on Cu.  FeCo / FeV q31 runs still
going (for the record).  Plain 36 steps (q27) already reach NiCu -4%, FeNi -3% vs SQS.

## One potential per species SQS (2026-10-18, desktop trial/q32): how much of ASR - SQS is intrinsic
Same SQS cells as q26 / q27 / q28 with `site_classes = false` (one ASA class = one self-consistent potential per species, the single-site
constraint ASR / CPA have) instead of one class per site.  All ASR numbers are energy-dependent (asr_ed, levels 4, rank 1).
Moments A / B (muB); SQS per-site = species mean of the site-resolved run; 2c = one class per species (16 #1 / #2; 32 #1).
| alloy | ASR | SQS per-site | SQS 2c | ASR vs 2c |
|---|---|---|---|---|
| bcc FeCo | r24 2.578 / 1.746 | 2.592 / 1.778 (16), 2.583 / 1.780 (32) | 2.582 / 1.793, 2.587 / 1.796 | -0.2% / -2.7% |
| fcc Fe0.5Ni0.5 | r36 2.555 / 0.641 (r24 2.512 / 0.593) | 2.535 / 0.662 | 2.516 / 0.680, 2.509 / 0.664 | +1.7% / -4.6% |
| fcc Ni0.75Cu0.25 | r36 Ni 0.493 (r24 0.319) | 0.517 / 0.507 | 0.517, 0.508 | -3.9% |
| bcc FeV | r24 1.848 / -0.469 | 1.756 / -0.364 (16), 1.808 / -0.382 (32) | 1.842 / -0.414, 1.795 / -0.393; 32: 1.852 / -0.399 | Fe -0.2%; V 0.07 muB |
| bcc FeRh | r24 2.866 / 0.729 (r48 2.877 / 0.632) | 2.977 / 0.630, 2.997 / 0.744; 32: 2.991 / 0.709 | 3.006 / 0.646, 3.016 / 0.741; 32: 3.009 / 0.715 | **Fe -4.8%**; Rh +2% |
Reading (user's bar: 3-4% between formalisms is agreement, do not chase exact numbers):
- FeV: the single-potential constraint raises SQS Fe by ~0.08 and V by ~0.03; ASR Fe then agrees (-0.2% vs 2c 32).  The "too magnetic" Fe was
  intrinsic (one potential per species), not a recursion error.  V stays 0.07 muB more negative (small induced moment).
- FeRh: the constraint changes nothing (Fe 3.01); ASR Fe -0.14 (-4.8%), flat in steps (24 -> 48).  The one genuine outlier.  Rank is unlikely
  (rank 1 = 2 = 3 to 3 digits in CuZn, "Queue 2 results"); not pursued further for now.
- Ni-type moments need 36 steps (r24: NiCu -38%, FeNi Ni -12%); r36 within ~4-5%.  FeCo agrees at 24.
Production recommendation: asr_ed, levels 4, rank 1, bp, kT 1 mRy; 24 steps for strong moments, 36 for weak (Ni-type) moments.

## Sublattices: one site-class operator for ordered / partially disordered / multicomponent alloys (2026-10-19, branch pd, desktop ~/asr_pd/trial/q33)
Goal (user): partial disorder -> multicomponent -> HEA with ONE machinery, no binary or Heusler shortcuts.
New code (binary route untouched, kept as regression reference):
- `asr_augcls`: augmented space where every site class (sublattice) has its own species list and concentrations (1 = ordered, 2, ..., N);
  local basis per class as asr_augn (W_lk = sum_k' X_k' x_k' phi_l phi_k), canonical-representative reduction as asr_augsym, the group
  supplied by the caller (class-preserving subset of D4h; identity = unreduced); order 2 (all asr_ed needs).
- `asr_edc`: asr_ed (AMscf form) on it, one recursion per (sublattice, species) component; `class_ops` (asr_symmetry): class-preserving
  D4h subset + cubic-site check; `asr_sublattice`: decoration of the sc/bcc/fcc graph by Wyckoff positions + centering (P/I/F),
  per-sublattice origin shift, site fractions, Madelung matrix of sublattice-average charges (Ewald over the decoration cell).
- `asr_class_scf` + `asr_class_driver` (`kind = "sublattice"`, `[[sublattice]] positions / species / concentrations`): one potential per
  component, V_p = sum Phi(k,k') (qbar_k' - qbar_cell) - sim_beta 2 (q_p - qbar_k)/R1 (Madelung of the sublattice averages + SIM of the
  deviation; = the binary SIM for one sublattice, = exact Madelung for ordered sublattices).
- `tools/sqs.py --decorate`: SQS on the disordered sublattice, ordered atoms at fixed offsets (L2_1 4b disorder).
- `atom_shells` leaves out empty valence shells (Al: an empty 3d is unbound in the free atom and stopped the start).
Gates (verified):
| gate | result |
|---|---|
| ring of 4, three classes (ordered {1} / binary {2,3} / ternary {1,3,4}, shared species), frozen start, vs exact configuration average | 2.7e-15 |
| decorated bcc (L2_1 pattern, ternary 4b), class-reduced vs unreduced, origin on 4a / 4b (D4h) and 8c (T_d -> D2d) | 4e-13; nodes 7288 -> 704, 8000 -> 1135 |
| one binary class vs asr_augsym / asr_ed | 3.4e-15, same node count |
| one ternary class (reduced) vs unreduced asr_augn / asr_edn | 1.8e-12; nodes 188625 -> 13600 (14x) |
| SCF binary limit: `kind = "sublattice"` (one sublattice, I-centred) vs `kind = "alloy"` asr_ed, CuZn and FeCo, levels 3 / 12 steps | iteration 1 identical to all printed digits (sum rules, E_F search, E_F, q, m); from it 2 ~1e-6 drift = neutral-Madelung convention (E_F-search residual charge) |
Same-model k-space references (VBH, equal spheres, nk 16):
| system | Co | Y | Al | total / f.u. | paper ASR (cond-mat/0606737 Table 3) |
|---|---|---|---|---|---|
| Co2FeAl a = 5.727 | 1.161 | Fe 2.758 | -0.102 | 4.98 (Slater-Pauling 5) | Co 1.119, Fe 2.709, Al -0.101 |
| Co2CrAl a = 5.727 | 0.902 | Cr 1.269 | -0.072 | 3.00 (Slater-Pauling 3) | Co 0.760, Cr 1.548, Al -0.071 |
| B2 FeCo a = 2.857 | Co 1.768 | Fe 2.760 | | | |
The Co/Cr split of Co2CrAl differs from the paper (total equal): sphere partition; our gate is the same-model k-space / SQS.
ASR results (q33: levels 4, 24 steps, rank 1, bp, kT 1 mRy, VBH) vs the same-model references (muB):
| system | site | ASR r24 | reference | diff |
|---|---|---|---|---|
| B2 FeCo | Fe / Co | 2.767 / 1.711 | k-space 2.760 / 1.768 | +0.3% / -3.2% |
| Co2FeAl | Co / Fe / Al | 1.106 / 2.694 / -0.107 | k-space 1.161 / 2.758 / -0.102 | -4.7% / -2.3% / ~0 |
| Co2CrAl | Co / Cr / Al | 0.860 / 1.189 / -0.070 | k-space 0.902 / 1.269 / -0.072 | -4.7% / -6.3% / ~0 |
| Co2Fe0.4Cr0.6Al | Co / Fe / Cr / Al | 0.936 / 2.751 / 1.237 / -0.082 | SQS 80 atoms x 2 (one potential per species) 1.017 / 2.849 / 1.186 / -0.086 | -8.0% / -3.4% / +4.3% / ~0 |
(The two SQS cells agree to 1e-3.  Co2CrAl: moments constant to 1e-3 from it 30, but dV stalled at ~1e-4 for 120 iterations -> FIXED, below.)
SCF did not stop at convergence (Co2CrAl) -- cause and fix (2026-10-19, desktop trial/q35): the E_F search exits at |N - N_elec| < 2e-4 and the
residual went into the sphere charges (self-Hartree 2 q/S ~1.5e-4 Ry): a dV floor above tol = 1e-5.  In metals false position usually lands far
inside the window; with E_F in a gap N(E) is flat and the residual is random in [0, 2e-4] (E_F jitter +-3e-4 Ry at constant charges; Anderson
amplified it to spikes of 1e-2).  Now the count is made exact after the search (moments interpolated linearly between the two bracketing
evaluations, N linear in the moments; else a first-order DOS step).  Co2CrAl from the stalled potential: converged in 12 iterations
(dV 3e-6, cell charge 6.5e-7 instead of ~1e-3), moments unchanged (Co 0.8595, Cr 1.189, Al -0.070).
Ported to asr_alloy_scf for the asr_ed route only.  The older P(z) - S route keeps the plain exit: its A = B test SCF (test_alloy_scf, Ag, 8 levels,
hard-core screening) has an antisymmetric A/B charge mode with gain ~ -9 per iteration (physical estimate ~ -0.6) that diverges as soon as the
count is accurate -- with the exact-count step, with a first-order step only, and with the OLD code at a 1e-6 tolerance (diagnostic clones,
desktop ~/asr_v1, ~/asr_v2); the 2e-4 residual noise happened to keep it damped.  Open issue for that (cross-check) route, not for asr_ed.
36 steps (q34): Co2FeAl Co 1.154 / Fe 2.721 / Al -0.104 (k-space 1.161 / 2.758 / -0.102: -0.6% / -1.3%); Co2CrAl Co 0.816 / Cr 1.273 /
Al -0.064 (k-space 0.902 / 1.269: Co -9.5%, Cr +0.3%); Co2Fe0.4Cr0.6Al running.
Total per formula unit: ASR 4.80 (Co2FeAl) and 3.63 (Co2Fe0.4Cr0.6Al) vs 5.00 / 3.80 (k-space / SQS = Slater-Pauling): the 24-step chain
smears the minority gap.  Co (minority d edge at E_F in these half-metals) is the low one, the weak-moment / chain-resolution case (as Ni);
strong moments (Fe, B2 FeCo) agree.  Paper ASR (Table 3): Co 1.075, Fe 2.617, Cr 0.972, Al -0.117, total 3.66 (same smearing).
Running (q34): the three Heuslers at 36 steps from the r24 potentials.
Co2Fe0.4Cr0.6Al at 36 steps did not stop either (fixed binary, q35/cfca_r36_fix): dV 2.5e-5 at it 15, then 2-3e-4; the Cr charge
alternates in sign between iterations (+-5e-4 e), Co / Fe / Al stable to 1e-4.  Ruled out (q35, same start potential): the adaptive depth
(depth_tol = 0: same floor), the BP tail (edges and depths frozen per chain once dV < 2e-3: same floor; the avg2 tail: no floor reached
either, slower, Co 0.930).  Probe (q36): one SCF evaluation at the r24 potential with the Cr potential shifted by t 5e-4 Ry, t = -1, 0, 1:
the output is linear to 2e-5 in the charges (Cr qnet -0.28250 / -0.29387 / -0.30522; production bp + adaptive depth) -> no noise; the Cr
response is dq/dV = -23 e/Ry, gain of the Cr potential ~ -12 (ordinary ASA charge sloshing, strong).  Cause found in asr_mixing: after a
history restart the new pairs were written at the old ring position while the fit uses columns 1..nh, so the newest pair was left out and
stale pairs from before the restart were fitted; restarts (|r| grew by 2) fire in exactly the near-converged regime.  Fixed (history refilled
from column 1; test_mixing checks the secant step after a restart).  Second change: in class_scf the two restart safeguards (|r| grew
by 2; step against the residual) act only while dV > 1e-3 -- near convergence they fired every 2-4 iterations and reduced Anderson to
linear mixing.  q37, same start potential, production settings (36 steps, bp, adaptive depth, mixing 0.1):
  old mixer: no stop (dV 1e-4 .. 3e-4 for 21 iterations);  ring fix only (m1): converged in 33 iterations (dV 9.2e-6);
  ring fix + no safeguards below 1e-3 (m3, committed): converged in 26 iterations (dV 1.9e-6).  Cell charge -7.5e-7.
  (m2, ring fix without the growth test only: still restarted by the residual test, no better; stopped.)
Converged r36 Co2Fe0.4Cr0.6Al (m3): Co 0.942, Fe 2.773, Cr 1.222, Al -0.082 muB (m1 agrees to 3e-4); vs SQS (one potential per
species) 1.017 / 2.849 / 1.186 / -0.086: -7.4% / -2.7% / +3.0% / ~0; total 3.645 muB / f.u. (Slater-Pauling 3.8; chain smearing as above).

## One input format, one program (2026-10-19)
`asr <input.toml>` reads one format (asr_driver): structure (decorated sc/bcc/fcc + [[sublattice]], or an explicit cell + [[site]]) and
species; `[calculation] method = "asr"` (default, every structure: pure, ordered compounds, partial disorder, binary .. HEA) or
`"kspace"` (ordered only; the decoration is reduced to its primitive cell by the centering translations that map every sublattice onto itself,
one ASA class per sublattice; B2 unequal spheres via [kspace] sphere_ratio).  Removed: kind = "lmto" / "alloy" / "model" readers
(asr_input_lmto keeps the types; asr_alloy_driver, asr_input deleted; asr_alloy_scf stays as a library route for the tests).  An old file gets
"kind is no longer read" with a pointer to method.  Examples converted (test_run_input reads all of them).  Gate: bcc Fe, B2 FeCo, B2 FeCu
(sphere ratio 1.0968) in k-space, old reader vs new: identical iterations, dV, E_F, charges and moments (desktop ~/asr_pd/trial/one).
