# Validation targets — what we compare against, and what each number is

Rules
1. A number from another code is a **target**, never a result of asr_new. Every row states source, code, functional, lattice constant, and whether *we* verified it.
2. asr_new is **LDA** (von Barth-Hedin `vbh`, libxc), ASA (space-filling spheres, spd, scalar-relativistic, no non-spherical terms).
   Do not gate on GGA/PBE numbers at tight tolerance: PBE enlarges lattice constants and magnetic moments relative to LDA
   (Fe moment, Fe-Co, Ni).  Prefer an LDA reference at the SAME lattice constant; otherwise state the expected sign/size of the bias.
3. Compare at the same lattice constant.  Compare relative energies (E_d - E_F, d-band splitting), integrated quantities (electron count,
   moments) and N(E_F); absolute E_F in Ry is only an internal same-representation regression number.
4. Oracle outputs from the Stuttgart code (AMscf/lmto) are same-representation regression targets, not physics references.
5. Digitised published figures are shape/position checks (peak positions), not precision numbers; peak heights depend on broadening.

| system | quantity | target value | source / code / functional | lattice | verified by us? | our value (LDA/ASA) |
|---|---|---|---|---|---|---|
| B2 FeCo (FM) | Fe, Co local moments, **LSDA** | **2.75, ~1.75 muB** (total ~4.50) | LSDA, FLAPW (Elk), experimental structure, arXiv:1907.06613 (Nature of the magnetic moment of cobalt in ordered FeCo) | **2.856 A** | text read by us (numbers 2.75 and "~1.75" quoted from the paper; not from its tables) | **2.760, 1.768 (a = 2.857 A): +0.4% / +1.0%** |
| B2 FeCo (FM), **real-space recursion** | Fe, Co moments | 2.75, ~1.75 muB (Elk LSDA, row above) | same | 2.857 A | see row above | `kind="alloy"`, `lattice="b2"`, depth 24, supercell 12, converged in 124 its: **2.745, 1.716** (total 4.461): Fe -0.2% / Co -1.9% vs Elk; first-order k-space route 2.760/1.768.  Depth/supercell convergence of the moments not yet checked. |
| B2 FeCo (FM) | Fe moment, LSDA, KKR-ASA | 2.66 muB | same paper (KKR in the atomic sphere approximation, quoted for comparison with DLM) | 2.856 | text read | 2.760 (+3.8%): ASA-to-ASA spread of ~4% between codes |
| B2 FeCo | LDA+DMFT, T = 1000 K | 2.69 / 1.82 muB | same paper | | text read | not comparable (finite T, correlations) |
| B2 FeCo (FM) | Fe, Co local moments | 2.80, 1.81 muB | VASP-PBE quoted in arXiv:2208.12973 (BOP paper; via Park et al. 2004) | not stated in what I read | source text checked | 2.760, 1.768 (a = 2.857 A); 2.743, 1.764 (a = 2.844 A) |
| B2 FeCo (FM) | total moment | 4.61 (= 2.80+1.81) | same | | derived | 4.528 / 4.507 |
| B2 FeCo (FM) | total moment | 4.530 muB/f.u., a = 2.844 A | pasted table, "DFT calculation" - **functional/code/source unknown** | 2.844 | NO | 4.507 (a = 2.844) |
| B2 FeAl | FM moment | ~0.7 muB/f.u. | VASP-PBE, arXiv:1605.06967; experiment paramagnetic | | source text checked | not run; LDA likely smaller - not a pass/fail gate |
| B2 FeAl | lattice | 2.872 A | pasted table, "DFT theoretical prediction" - functional unknown | | NO | |
| B2 CuZn | lattice | 2.966 A (DFT) vs 2.959 A (exp) | pasted table - functional unknown; Materials Project mp-987: 2.9206 A | | NO | we use 2.85 A (provenance unconfirmed) |
| B2 CuZn | Zn d / Cu d vs E_F, d-splitting | -0.682 / -0.309 / 0.373 Ry | Banerjee & Mookerjee 2010 Fig. 3, digitised (TB-LMTO-ASR, LDA, 50:50 disordered) | not stated | figure checked by eye | EI depth 4-5: Zn -0.65..-0.68, Cu -0.26..-0.34, sep 0.35-0.42 (lattice-dependent, see STATUS) |
| B2 CuZn (ordered) | Cu d peaks vs E_F (main / second) | **-2.8 / -3.4 to -3.5 eV** (Zn d below the -5 eV plot cut, not read) | Materials Project mp-987 DOS plot, **PBE**, digitised by eye from a screenshot (+-0.2 eV) | **2.92055 A** (CIF `Cu50Zn50_2.cif`, Pm-3m) | figure read by us; not the underlying data | **-2.80 / -3.47 eV** (E_F -0.0809 Ry; Cu d -0.206/-0.255 Ry, Zn d -0.555 Ry = -7.6 eV; Cu/Zn net charge +-0.135 e; nonmagnetic). Agrees within reading error. Zn d depth NOT checked. |
| B2 CuZn (ordered), **real-space recursion** | same Cu d / Zn d peaks vs E_F | as the MP row above (Cu d -2.8 / -3.4 eV) | same | 2.92055 A | see MP row | recursion depth 32, supercell 16 (`kind="alloy"`, `lattice="b2"`): E_F -0.0503 Ry; Cu d **-0.238 Ry = -3.2 eV**, Zn d -0.619 Ry; Cu +0.137 e.  At the same parameters it agrees with k-space H^(2) in E_F to 2.6 mRy (depth 24, `test_alloy_b2`).  It is ~0.4 eV deeper in Cu d than MP / the first-order k-space SCF because the real-space route uses H^(2)=E+h-hoh. |
| bcc CuZn | same, Tarafder 2006 | -0.488 / -0.122 / 0.366 | digitised; non-self-consistent starting potentials (per review doc; not verified by us) | | partly | not a target |
| fcc Cu-Ni | moment onset ~46 at.% Ni; ~linear above ~25 at.% | | pasted table, source unknown | | NO | not testable yet (needs depth ~25-30 for Ni) |
| fcc Ag | d bands ~0.2 Ry below E_F | | pasted table, (L)APW; "about" | | NO | Ag matches Stuttgart oracle within ~1% on parameters |
| bcc Fe | moment | 2.15-2.30 | literature range (LDA/GGA mixed) | 2.8665 A | range only | ordered 2.20-2.23; alloy route 2.35 |
| fcc Ni | moment | 0.58-0.68 | literature range | 3.524 A | range only | ordered 0.639 |
| B2 FeCu | E_F, occupations | -0.1871 Ry etc. | Stuttgart oracle (same representation) | 3.0 A | oracle | -0.1932 |

Open action: a VASP or QE **LDA** run (same lattice constants) for B2 FeCo and B2 CuZn would separate functional effects from the ASA
representation; not available in this environment.

## Manuscript to-do (added 2026-10-05)
- DOS figures needed for B2 FeCo (spin-resolved, Fe/Co) and B2 CuZn (Cu/Zn d) at the lattice constants above; overlay the MP/Elk
  reference positions as markers, not as digitised curves.  Need the underlying MP DOS data (or an uncropped plot covering -10..-5 eV)
  to check the CuZn Zn d position.
- Structure files: `Cu50Zn50_2.cif` = ordered B2 CuZn, mp-987 (usable).  `Cu50Zn50.cif` is mislabelled: 16 Zn + 10 Cu, I-43m, 26 atoms
  = Cu5Zn8 gamma-brass, not 50:50; open structure (12.95 A^3/atom) -> ASA without empty spheres is expected to be unreliable; not attempted.
- FeCo moments (Fe 2.760 / Co 1.768 at 2.857 A vs Elk LSDA 2.75 / ~1.75): within +0.4% / +1.0%; this is the magnetic validation.

## What the published literature does and does not give (searched 2026-10-04)
- **B2 FeCo:** a proper LDA-family reference exists: Elk FLAPW LSDA at a = 2.856 A, Fe 2.75 / Co ~1.75 muB (arXiv:1907.06613).  Our LDA/ASA at the
  same lattice constant agrees within 0.4% (Fe) and 1% (Co); a KKR-ASA LSDA value in the same paper (Fe 2.66) is 4% lower than ours, so ASA-vs-ASA
  differences of a few % exist between codes.  The earlier PBE numbers (2.80/1.81) are 1.5-3.5% above the LSDA ones, the expected PBE bias.
- **B2 CuZn:** I found NO first-principles (LDA/GGA) paper with numerical d-band positions or a lattice constant that I could read.
  Both TB-LMTO-ASR papers we have (Tarafder 2006 = arXiv:cond-mat/0606623 / 0607767, Banerjee & Mookerjee 2010) give no lattice constant in the text
  I could extract, and their curves differ by ~0.2 Ry in absolute d-peak positions relative to E_F (Zn d -0.49 vs -0.68 Ry).  Tarafder's text says the
  "density functional self-consistent TB-LMTO takes care of the charge transfer", which does NOT support the earlier review-note claim that the 2006
  curves came from non-self-consistent elemental potentials (not verified either way).  Reported for context: LDA places Zn 3d too shallow compared to
  photoemission (6.8-8.5 eV in ZnO; XPS for brass ~9.5-10 eV), so the experimental Zn-d position is not an LDA target.
  => There is no trustworthy external CuZn d-band target at the moment; our own checks are internal (EI vs ED routes, depth/rank convergence).
