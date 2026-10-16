# Disordered magnetism: what is left (2026-10-16)

Production route: asr_ed (AMscf form), levels 4, recursion_steps 24, bp tail, SIM beta 1, VBH or PBE.
Clean results so far: docs/STATUS.md "Clean reruns of all disordered binaries, VBH and PBE" (desktop trial/q24, q23).

## Done (do not redo)
- [x] SIM sign and magnitude (326759b): V = -beta 2 qnet / R1.
- [x] SCF runaway: Anderson metric weighted by (r/(1+r))^2 + restart safeguard (ca48346); E_nu kept inside its band window (7bc52d2 ff.).
- [x] E_F search: bracketed false position (17416b2); species charges balance.
- [x] Antiparallel starts: [[species]] spin_sign (da1a24b); restarts: save_potential / start_potential / flip_spin (6692e95).
- [x] depth_tol (adaptive chain length) shown harmless at 12 / 24 / 48 steps (q16).
- [x] Our LMTO-ASA vs Questaal lm: band structure and charge transfer agree; Fe moment 1% low (STATUS "Our LMTO-ASA vs Questaal").

## Open, in order
1. [ ] E_F noise floor (dV stalls at 1e-4..1e-3 in FeRh VBH, FeV PBE): Fermi-smeared contour implemented ([asr] fermi_kt, cc2d513);
       test running (trial/q25, kT = 3 mRy).  Gate: dV < 1e-5 with moments within 0.01 muB of the stalled runs.  Then decide the default kT.
2. [ ] Independent disordered reference: SQS supercells with our own k-space LMTO-ASA (same spheres, functional, parameters; no recursion).
       FeCo, FeV, FeRh at x = 0.5 (16 and 32 atoms, two SQS each; per-site classes -> local-environment moments).  Gate: ASR species moments
       within the SQS spread.  Settles whether the FeRh Fe deficit vs CPA (0.20 muB) is ASR truncation or the ASA model.
3. [ ] Recursion depth protocol for moments: MnCr changes 7-17% from 24 to 48 steps, FeV changed a lot from 12 to 24.  Options: always
       report 24 / 48; extrapolate in steps; better tail (bp vs hull vs multi-band).  Gate: moments stable to 0.02 muB in steps, matching SQS (2).
4. [ ] Generated depth (levels) at fixed steps: the hard-wall node set is levels 4; test levels 5 at 24 / 48 steps for FeV and MnCr.
5. [ ] The 1% weaker magnetism vs Questaal (bcc Fe 2.252 vs 2.275, VWN): trace E_nu choice (~3 mRy lower in Questaal), core polarisation,
       Bloechl correction.  Matters only for marginal magnets (B2 MnCr flips to non-magnetic at converged k).
6. [ ] SCF start: FeCo still overshoots E_F to +0.31 Ry from the atomic start; dV sat at the 0.15 Ry clamp for several iterations in the alloy
       SCF.  Add the k-space `trace` diagnostics to the alloy SCF and check V_in / V_out there.  Option: start alloys from converged ordered
       or pure-element potentials (start_potential).
7. [ ] Magnetic state selection: always run parallel and antiparallel starts for systems with antiferromagnetic tendencies (Mn, Cr, V);
       compare total energies (needs the alloy total energy, not implemented) instead of picking by start.
8. [ ] Model limits (document, do not fix here): one collinear moment per species cannot describe sublattice AF, local-environment moment
       flips (Mn in FeMn), or non-collinear order (Kaphle 2012, Ganguly 2011).  Non-collinear ASR is a separate project.
9. [ ] References to re-check with clean runs: FeV vs Chakrabarti-Mookerjee 2005 (theirs 1.12 / -0.11; ours 1.85 / -0.47 VBH);
       MnCr vs Banerjee-Mookerjee 2009 Fig. 1 (not reproduced; legacy band 0.1-0.2 Ry wider than ours and Questaal).

## After this (separate projects)
- Partial (sublattice) disorder: general lattice + basis graph, per-site configuration spaces (ordered sites never enter the cardinality
  sequence), automatic site symmetry, multi-class SCF with Ewald + SIM.  Target: Co2Fe0.4Cr0.6Al (Chakraborty 2006; experiment).
- HEA: n-component augmented space exists (asr_augn / asr_edn); needs the multi-class SCF and symmetry reduction from the item above.
- Kubo-Greenwood optical conductivity (Sadhukhan 2017): LMTO current operator + two-Green-function average (KPM route).
