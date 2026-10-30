# J_ij plan: exchange parameters in energy-dependent ASR  (2026-10-21)

Status: deferred, nothing implemented.  Starts after SRO (HANDOFF_hea_sro.md Part B).  Literature marked *(unverified)* is from memory.
## 1. Formula
Magnetic force theorem (Lichtenstein-Katsnelson-Antropov-Gubanov) in the TB-LMTO form of Turek, Kudrnovsky, Drchal *(unverified
citation details)*, for site i occupied by Q and site j by Q':

    J_ij^{QQ'} = (1/4 pi) Im Int^{E_F} dz tr_L [ dP_i^Q(z) g_ij^{QQ',up}(z) dP_j^{Q'}(z) g_ji^{Q'Q,dn}(z) ],
    dP = P^up - P^dn,   g = (P - S)^-1 (auxiliary Green function),

with g the conditional average with i and j frozen.  Sign and normalisation convention (Heisenberg H = -sum J e_i.e_j, per pair or
per bond) to be fixed against an ordered-metal gate.

## 2. Why no orbital peeling for J_ij
The Mookerjee-era codes had only the on-site G_ii from scalar recursion, and pair energies as differences of large, nearly equal E_ij
terms lost the difference in numerical error (Chowdhuri thesis section 2.7.1).  Orbital peeling (Burke 1976) fixed that cancellation.
With an inter-site g_ij available the formula above is evaluated directly; peeling is not needed.

## 3. Implementation on the production route
- **ED route only** (production route, STATUS "AMscf-form energy-dependent ASR"): `asr_ed` already works with P(E) - S in the Delta_c scaling, so g is what its chains
  produce up to the similarity T = Delta_c^1/2 Delta^-1/2.  With two frozen sites, undo the scaling at both ends (T_i, T_j); derive and
  test it, do not substitute the EI H^(2) route.
- **Two frozen sites**: i -> Q, j -> Q' (generalise the single frozen start site); remaining sites averaged (with SRO if active).
  Symmetry: subgroup fixing the pair (bond little group); more nodes than the single-site chain.
- **Off-diagonal block**: block Lanczos / block recursion with the start block = orbitals of i and j (2 x 9 spd, or 2 x 5 d) gives the
  whole i-j block of g in one chain per energy.  Alternative (more chains): scalar recursions from |i> +- |j>, |i> +- i|j>.  Block
  recursion does not exist in the code yet; it is the new numerical piece (block terminator too).
- **Contour**: same closed contour as the moments; J from the converged potential (post-processing, no SCF).
- **Cost**: one block chain per spin x contour point x species pair (Fe-Fe, Fe-Cr, Cr-Cr) x neighbour shell; parallel over all of them
  (same OpenMP-over-chains scheme as class_scf).

## 4. Orbital peeling: where it is still the right tool
Pair energies for a finite perturbation (species exchange: effective pair interactions V^(2) for ordering / SRO energetics) by the
Lloyd-formula log det difference, evaluated from diagonal recursions of successively peeled Hamiltonians.  Cost: one recursion per
peeled orbital per configuration.  Use only if the Fe-Cr paper needs ordering energies (Fe-Cr ordering energies, docs/PUBLICATION_PLAN.md section 2 stretch).

## 5. Approximation to state in any paper
The up and down propagators see the same configuration; averaging them separately drops that correlation (vertex correction).  CPA
work makes the same approximation and argues it is small *(unverified)*.  Conditioning on both sites keeps environment effects CPA
misses, but J_ij is still not exact.

## 6. Gates
1. Ordered bcc Fe and fcc Ni (frozen-pair chain on the pure metal == k-space J_ij from the same LMTO parameters, our own k-space code):
   first shells and the sum rule J_0 = sum_j J_0j vs the on-site expression.
2. Block recursion vs scalar |i> +- |j> recursions on a disordered binary (same g_ij).
3. Ordered B2 FeCo vs a published LSDA J_ij *(reference to find)*.
4. Disordered bcc Fe-Co or Fe-Cr vs KKR-CPA J_ij *(reference to find)*: trends, qualitatively and at the 3-4% bar
   used for the moments.
