---
title: Validation
---

# Validation

The code's policy: **only verified claims are stated**. Every row below is
either an automated test (`fpm test`, harness in [[asr_testing]]) or a measurement recorded in
`docs/STATUS.md`. A number taken from another code is a *target*, never a result of
asr; targets are listed with source, functional and lattice constant in
`docs/VALIDATION_TARGETS.md`.

!!! note "What asr is, for validation purposes"
    LDA (von Barth–Hedin `vbh` through libxc), ASA (space-filling spheres,
    spd, scalar-relativistic, no non-spherical terms). Do not compare to GGA
    numbers at tight tolerance (PBE enlarges lattice constants and moments).
    Compare at the same lattice constant; relative energies (\(E_d-E_F\), splittings), electron
    counts and moments are meaningful, absolute \(E_F\) is only an internal
    regression number.

## Test suite

```bash
source env.sh
asr_test                  # everything, about 3 minutes on 16 cores
asr_test test_radial      # a single test (fpm test <name>)
```

| Area | Tests | Claim |
|---|---|---|
| Recursion | `test_recursion_dense` | Exhausted Krylov space equals dense diagonalisation (MKL `dsyev`): \(G(z)\) to 10⁻¹⁴, second moment to 10⁻¹⁵ |
| | `test_ordered_dos` | NN band on sc/bcc/fcc: \(a_1=0\), \(b_1^2=z\) (6, 8, 12); LDOS vs k-space (nk = 128, η = 0.1) error falls with depth, 2.6·10⁻⁴ / 3.9·10⁻⁴ at depth 40 (sc, fcc) |
| Terminators | `test_terminators` | unit tests of every asymptote/join, plus the gapped-band (B2) stress test |
| Contour | `test_contour` | Contour moments from a terminated continued fraction equal the tetrahedron method (independent) for the ordered sc band (\(M_0\) = occupied states, \(M_1\)); the Fermi-smeared contour reproduces \(f(\varepsilon)\) occupation |
| Radial / atom / xc | `test_radial`, `test_atom`, `test_xc`, `test_band_window` | hydrogenic levels to 10⁻¹⁰ Ry; free atoms He (−2.834836 Ha), Ne (−128.233481 Ha) vs NIST LDA (VWN) to 3·10⁻⁷ Ha; libxc binding (unpolarised = polarised at ζ = 0, spin-splitting sign); the valence window of the radial solution contains \(E_\nu\) and \(C\) |
| Parameters | `test_potpar`, `test_harmonics`, `test_strconst`, `test_orbital_rep` | \(C\) agrees with the exact zero of \(W\{K,\phi\}\) (error second order in \(E_\nu-C\)); Gaunt coefficients exact; orbital representation of D\(_{4h}\) is a group homomorphism; screened fcc on-site block equals the legacy `INFO_FCC` to 6·10⁻⁴ (s 3.2326, p 3.0223, t2g 2.7422, eg 1.9492) |
| Ordered LMTO | `test_lmto_ham`, `test_lmto_cell`, `test_asa_ag`, `test_asa_mag`, `test_asa_b2` | \(H(\\mathbf k)\) Hermitian with the correct cubic degeneracies; Ag vs Stuttgart; cell code equals the single-atom code at generic k to 3·10⁻¹⁵; moments of Fe, Ni, B2 systems |
| BZ / electrostatics | `test_tetra`, `test_ewald`, `test_mixing` | tetrahedra: free-electron sphere volume, half-filled tight binding, weight sum rule, grid convergence; Ewald: CsCl 1.762675, NaCl 1.747565, independent of the cell choice; Anderson mixing converges on a stiff map where linear mixing diverges |
| Augmented space | `test_augspace`, `test_augorb`, `test_augn`, `test_augcls`, `test_augsym`, `test_augsym_orb`, `test_partial` | Recursion equals the exact average over **all** \(2^N\) configurations to 10⁻¹⁴ (also multi-orbital, n-species, class-resolved); the point-group-reduced operators give identical Lanczos coefficients to the unreduced ones; species-resolved \(\langle\langle n_RG_{RR}\rangle\rangle\) equals the occupation-weighted enumeration |
| Alloy LMTO | `test_alloy_ordered`, `test_lmto_real`, `test_lmto_real3`, `test_alloy_scf`, `test_alloy_scf_spin`, `test_alloy_b2` | Real-space augmented-space recursion equals the ordered k-space recursion coefficients to 2.5·10⁻¹⁴ (s) / 1.6·10⁻¹⁴ (e\(_g\)) over 6 levels; the trivial-alloy limit reproduces the ordered SCF; third-order real-space Hamiltonian equals dense matrix algebra; ordered B2 CuZn real-space SCF against the k-space \(H^{(2)}\) at the same parameters (separates recursion error from the representation difference — \(H^{(2)}\) misses the \(hph\) overlap term and puts \(E_F\) ≈ 35 mRy too high, which \(H^{(3)}\) must close) |
| Energy-dependent | `test_endep`, `test_endep_ordered`, `test_pfun`, `test_asr_ed`, `test_asr_ed_n`, `test_asr_ed_sym` | Ordered limit equals the exact k-space inverse of \(P(z)-S(\mathbf k)\); disordered species-resolved \(G\) equals ensemble averages of dense inverses (ring, brute force) for `asr_endep`, `asr_ed` (binary), `asr_edn` (3 species); reduced = unreduced; exact potential function equals the real-axis radial solver on the real axis |
| KPM | `test_kpm` | KPM on a periodic sc lattice reproduces the exact k-space DOS (independent truth for disordered tests) |
| Input | `test_input`, `test_lmto_input`, `test_alloy_input` | Unknown keys and bad values are errors, never silent defaults |

`probe_*` programs in `test/` are exploratory studies (e.g. terminator depth
scans, route comparisons) that print results rather than assert them.

## Headline results

| Claim | Evidence |
|---|---|
| Own ASA SCF reproduces Stuttgart Ag parameters \((\\sqrt\\Delta,\\gamma,p)\) | ≈ 1 % |
| d weight 9.58 vs 9.50; d centroid within 0.6 mRy; occupations s/p/d 0.687/0.633/9.635 vs 0.683/0.648/9.565; \(N(E_F)\) 3.165 vs 3.177 /Ry/atom | `test_lmto_ham`, `probe_ag_dos` (parameters from `data/reference/ag_fcc_stuttgart`) |
| Fe 2.2 μ\(_B\), Ni 0.64 μ\(_B\), Cu 0 | ordered ASA SCF |
| B2 FeCu occupations | 0.01–0.02 e vs Stuttgart (d: Fe 6.766 / Cu 9.429 vs 6.763 / 9.426) |
| B2 FeCo (FM) moments | Fe 2.760, Co 1.768 μ\(_B\) (k-space) vs Elk LSDA 2.75 / ≈ 1.75 (2.856 Å); real-space recursion, depth 24: 2.745 / 1.716 |
| Augmented-space recursion = exact configuration average | 10⁻¹⁴ |
| Real-space augmented-space recursion = ordered k-space recursion | 10⁻¹⁴ |
| Alloy SCF in the trivial-alloy limit = ordered SCF | `test_alloy_scf` (14 assertions) |
| Energy-dependent exact KKR-ASA = same-model k-space KKR | ordered limit passed |
| Sum rule: for \(E_F\) far above the bands \(M_0=1,3,5\) (s, p, d) | exact for Fe at 8 and 14 levels |

## Known limitations

!!! warning "Open problems — do not quote around them"
    * **bcc terminator.** Single-band tails give 1.4·10⁻² LDOS error at depth 40,
      falling only slowly; the k-space reference is converged (nk 128 vs 256
      differ ≈ 10⁻³). The cause is the van Hove singularity at the band centre and
      even/odd oscillation of \(b_n^2\).
    * **Depth.** On a real d band ≈ 1 % LDOS error needs depth ≈ 40 for an *ordered*
      Ag test, while the augmented space reaches only ≈ 10–20 levels. No
      terminator removes this.
    * **Magnetic d-band alloys** need more levels than exact generation affords
      without symmetry (Ni loses its moment at depth 3).
    * **Fe in the alloy route** gives ≈ 5–6 % more moment than the ordered
      first-order ASA route; depth-converged and not yet explained.
    * **Luchini–Nex joiner** is implemented from a secondary description and not
      checked against the original paper.
    * **s,p occupations** are 3–5 % above the Stuttgart moments for B2 FeCu
      (≈ 0.1 e surplus); hypothesis (o and p terms) untested.
    * Recursion is **real-symmetric and scalar**; no block / complex-Hermitian
      recursion (needed for non-collinear magnetism).
    * The CuZn comparison with Banerjee & Mookerjee (2010): the 50:50 point
      agrees (Zn d 4 mRy, Cu d 21 mRy at depth 4, rank 1), the concentration trend
      is opposite to that paper, whose Fig. 3 is internally inconsistent
      (Cu\(_{0.75}\) panel holds +1.04 e below its "E\(_F\)").

## Reference data

| Path | Content |
|---|---|
| `data/reference/ag_fcc_stuttgart/` | Stuttgart `ctrl` and s, p, d, f DOS of Ag (regression target only) |
| `data/reference/mncr/banerjee2009_fig1_50-50.dat` | Digitised published DOS (shape and peak positions only) |
