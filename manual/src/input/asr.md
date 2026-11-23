---
title: Recursion and disorder (method = "asr")
---

# Input: the recursion

`method = "asr"` is the default. The keys below are read only for it.

## `[asr]`

### Augmented space

| Key | Default | Meaning |
|---|---|---|
| `levels` | `4` | Recursion depth, ≥ 2 |
| `recursion_steps` | `0` | Lanczos steps; `0` = `levels`; otherwise ≥ `levels` |
| `rank` | `1` | Maximum number of simultaneously flipped sites in an augmented-space state (≥ 0). Sites beyond are an average medium. Rank 2 is a refinement; for soft-moment alloys it changes moments by several percent ([Convergence](../guides/convergence.html)) |
| `supercell` | `8` | Graph size in conventional cells; a multiple of `period` |
| `flip_radius` | `99` | Only sites within this distance may flip |
| `cluster_radius`, `inversion_radius` | `1.2`, `0.0` | Real-space cluster for the structure constants and for their inversion (in `alat`; `0.0` keeps `cluster_radius`) |

### Energy integration

| Key | Default | Meaning |
|---|---|---|
| `contour_points` | `32` | Gauss–Legendre nodes on the complex-energy contour |
| `fermi_kt` | `0.001` | Electronic temperature (Ry); gives a smooth \(N(E_F)\) |
| `fermi_poles` | `0` | Matsubara poles enclosed (0 = automatic, range 0 to 200) |
| `ed_terminator` | `hull` | Tail of each chain: `hull`, `bp`, `avg2`. Keep `hull` ([Recursion and terminators](../theory/recursion.html)) |
| `depth_tol` | `1e-5` | A chain stops early when \(G\) changes by less than this; `levels` is the maximum |

### Route (Hamiltonian-dependent)

| Key | Default | Meaning |
|---|---|---|
| `route` | `edc` | `edc`: AMscf energy-dependent form of \(H^{(2)}\). `ed3`: the full generalised ASA problem including the \(hph\) term. See [Energy-dependent ASR](../theory/energy_dependent.html) |

### Speed

| Key | Default | Meaning |
|---|---|---|
| `energy_block` | `0` | Contour energies per chain block. `0` runs one chain per energy. `8` runs eight in one block, 1.35× faster for CoCrFeNi at about 6× the memory. Same result within the chain's length error. The chains of `[epi]` use the same setting |
| `ef_linear_max` | `0` | Ry. \(E_F\) steps up to this size are taken to first order without a new contour (exact at self-consistency); `0` = off |

An MPI build (`-DUSE_MPI`) splits the chain tasks over ranks with no input change.  Build it with the flags of your MPI
library added to the usual ones, for OpenMPI `--flag "... -DUSE_MPI -I<mpi>/include" --link-flag "... -L<mpi>/lib -lmpi_mpifh
-lmpi"` (MPICH and Intel MPI: `-lmpifort -lmpi`).  Run one rank per node with all its cores as OpenMP threads.  A serial build
started by `mpirun` stops at once: its copies would all write the same files.

Threads.  A chain task is one start orbital, species, spin and contour energy, so a run has
\(4\,N_{\text{species}}\,N_{\text{spin}}\,N_E\) of them (random CuZn, 32 energies: 256).  While every thread has several, each
thread runs whole chains.  When a rank holds fewer chains than threads (many nodes), teams of threads share the
matrix-vector products of one chain (nested OpenMP); the log says `chains: ... teams of ... threads`.  The results do not
depend on the split.  With the Intel runtime set `KMP_HOT_TEAMS_MAX_LEVEL=2`, so that the inner teams are kept between
products (gfortran's runtime keeps them anyway).  Keep the maths library single-threaded (`MKL_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`): its calls come from inside the threads.  Thread binding is optional.  If you bind, use
`OMP_PLACES=cores` with `OMP_PROC_BIND=spread,close`: `spread` gives each outer thread its own share of the cores, and its
team stays inside that share.  `close` at the outer level does not split the cores between the teams, so teams can land on
the same cores.

With several ranks, the timing section of `asr.out` (and `timing_s` in `asr.json`) gives each rank's own chain time:
min, mean, max and max/mean.  A max/mean well above 1 means some ranks wait for others at the sum over ranks.

### Diagnostics

| Key | Default | Meaning |
|---|---|---|
| `ef_scan` | `0` | Diagnostic: scan \(N(E)\) on `2 * ef_scan + 1` energies around `ef_scan_center` at the start potential, then stop |
| `ef_scan_step` | `2e-5` | Energy step (Ry) for `ef_scan` |
| `ef_scan_center` | start \(E_F\) | Centre energy (Ry) for `ef_scan` |
| `ef_scan_dump` | — | File: dump \(N(E)\) at the start potential, then stop |
| `alpha_scale` | `1` | Multiply the screening constants by this factor |

### Constraints enforced by the reader

* `levels` ≥ 2; `recursion_steps` is 0 or ≥ `levels`.
* `supercell` is a multiple of `period`.
* `fermi_poles` in 0 to 200; `route` is `edc` or `ed3`.
* `energy_block` ≥ 0; `ef_linear_max` ≥ 0; `alpha_scale` > 0.

## `[[sro]]`

First-shell Warren–Cowley short-range order (decorated lattice, `method = "asr"`).
The first shell of the frozen centre carries the conditional concentrations
\(P(s\mid c)\) ([[asr_sro]]); shell sites are independent of each other, and sites
beyond the first shell are uncorrelated.

| Key | Meaning |
|---|---|
| `sublattices` | `["A", "B"]`; optional with one sublattice |
| `shell` | `1` |
| `pairs` | `[["Fe", "Ni", alpha], …]`. \(\alpha>0\): the pair avoids contact (clustering); \(\alpha=1\): never neighbours; \(\alpha<0\): ordering |

With one sublattice the unlike pairs are symmetric and the diagonal is derived.
With two sublattices, `[a, b, alpha]` is species `a` on the first and `b` on the
second; unset pairs are 0 and both marginals must vanish. Example:
`example/fe66ni14cr20_fcc_sro.toml`.

The electrostatics of the SRO follow the same model: the first-shell structure factor
is used wherever it is non-negative, and beyond that range it is completed with as
little further-shell order as possible ([SCF and electrostatics](../theory/scf_electrostatics.html)).
Every \(\alpha\) up to complete clustering is allowed. A species with \(\alpha=1\) to
all others (only its own kind as neighbours) stays neutral. `asr.out` lists the SRO
modes and their screened potentials. With `sim_beta` above about 1 the result for strong
clustering depends on the order beyond the first shell, which this input does not fix.

## `[epi]`

Effective pair interactions of the converged alloy ([[asr_epi]]). They are the input
of a Monte Carlo or cluster-variation treatment of the short-range order, or can be
compared with interactions fitted to experiment. The section needs
`method = "asr"`. Both routes work. `edc` uses \(H^{(2)}\), and `ed3` uses the full
generalised ASA problem of the SCF.

| Key | Default | Meaning |
|---|---|---|
| `shells` | `0` | Coordination shells per sublattice. `0` switches the calculation off |
| `lambda_points` | `4` | Gauss–Legendre nodes of the species path |
| `symmetric` | `false` | `true`: on one sublattice only \(U(c,s)\) with \(s\ge c\) is computed, and \(U(s,c)=U(c,s)\) is copied |
| `file` | `epi.toml` | Pair tables (TOML) |

The run computes the pair table \(U(c,s)\) for every pair type after the SCF.
A pair type is one shell of one pair of sublattices. Species \(c\) sits on the
origin and species \(s\) on the partner site. Every other site is averaged with
the concentrations of its sublattice.

\[
U(c,s)=\Omega(c,s)-\Omega(c_0,s)-\Omega(c,s_0)+\Omega(c_0,s_0)
\]

\(\Omega\) is the grand potential at frozen potentials. \(c_0\) and \(s_0\) are
the first species of the two sublattices. The pair model has the energy
\(E=\sum_{i<j}U(s_i,s_j)\). On one sublattice the exchange interaction of an unlike
pair is

\[
V(p,q)=U(p,p)+U(q,q)-U(p,q)-U(q,p).
\]

\(V>0\) favours unlike neighbours. The SRO energy per site is
\(\sum_n \tfrac{z_n}{2}\sum_{p<q}x_p x_q V_n(p,q)\,\alpha_n(p,q)\).

The table has two parts.

* **Band part.** The species of the origin moves from \(c_0\) to \(c\) on a straight
  path in the parameters of \(H^{(2)}\). The derivative of \(\Omega\) along the path
  is \(\mathrm{Tr}[f(H)\,\partial H/\partial\lambda]\). It needs only the diagonal
  blocks \(G\), \(hG\) and \(hGh\) of the origin. The expression is linear in \(G\).
  So the configuration average is exact in augmented space. The partner site is
  frozen. The path integral uses `lambda_points` nodes.
* **Screened Coulomb part.** \(W(R)\,\delta b_c\,\delta b_s\) with the pair kernel
  \(W\) of the charge-excess model and its drive deviations \(\delta b\)
  ([SCF and electrostatics](../theory/scf_electrostatics.html)).

`asr.out` lists both parts and their sum. \(U(p,q)\) and \(U(q,p)\) come from
different chains. Their difference measures the error of the chains and of the path
quadrature. The chains run the full `recursion_steps`. The medium around the pair is
uncorrelated, also in a run with `[[sro]]`.

Cost. One pair type needs
\((n_1-1)\times\) `lambda_points` \(\times\,n_2\times 9\times 2\) chains per
contour energy and spin. \(n_1\) and \(n_2\) are the numbers of species on the two
sublattices. One Fermi-level pass of the SCF needs 4 chains per component. The pair
group has order 2 to 4, so the node set is 4 to 8 times that of the SCF. A quaternary
alloy with 4 path points costs roughly 200 to 400 Fermi-level passes of the SCF per
shell. Use MPI for alloys with four or more species. Three settings reduce the cost.

* **Path points.** Two points are enough for species that are close, such as neighbours
  in a transition series. The unit test `test_epi` finds an error of 0.01 % there.
  One point is the midpoint rule. It costs a quarter of the default. It is not safe.
  Its error is 2 % of the largest entry and can reach a factor of two on a small entry.
  Use it for a first survey only. Species that differ strongly need four points (0.1 %).
* **`symmetric = true`.** The exact table of one sublattice is symmetric. The option
  computes one half. For four species the chains drop by a quarter, for five species
  by 30 %. The error measure \(U(p,q)-U(q,p)\) is then not available. Check one shell
  without the option first.
* **`energy_block` in `[asr]`.** The pair chains run in blocks of contour energies, as
  in the SCF. The block size is capped so that one thread needs less than 400 MB. The
  table changes at rounding level only. The gain for the pair chains is not yet
  measured. The SCF gains 1.35× with `energy_block = 8`.

`epi.toml` holds one `[[pair]]` table per pair type, energies in Ry. Any Monte Carlo
code of the pair model \(E=\sum_{i<j}U(s_i,s_j)\) can read it.
