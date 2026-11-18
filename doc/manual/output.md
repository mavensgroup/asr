---
title: Output files
---

# Output files

All files are written to the directory of the input file. Every
`method` prints a short summary to the terminal. The full information is in the
report (`out_file`, default `asr.out`), in its machine-readable copy (`asr.json`) and in
plain-text data files. The data files load directly into gnuplot, matplotlib or
`numpy.loadtxt`; lines starting with `#` are comments.

## The report `asr.out`

Plain text, at most 90 columns ([[asr_report]]). Every number in the report is a
number the run computed, or a stated formula of such numbers. Every table header
carries its unit. Section headers have a fixed form `=== name ===…`, so the file can
be grepped:

```bash
grep -A12 "=== health checks" asr.out
sed -n '/=== total energy/,/=== converged results/p' asr.out
tail -2 asr.out          # the status line
```

| Section | Content |
|---|---|
| *(banner)* | Program version, git revision and whether the working tree had changes, compiler and flags, start time, host, process id, MPI ranks and OpenMP threads, input file and its fingerprint, working directory, units and conversions |
| `input` | Every input key the run used, with its value, in TOML form. Keys that are not in the input file carry the mark `# default`. Sections that the chosen method does not read are left out. The block is valid TOML: saved to a file, it reruns the same calculation. |
| `method` | Hamiltonian and route, basis, XC functional; for the recursion: levels, recursion steps, rank, supercell, augmented-space nodes per sublattice, screening \(\alpha\), terminator, the occupation contour (nodes, Matsubara poles, \(kT\), bottom, \(E_F\) probe \(\eta\)), electrostatics and the screening model; for k-space: screening \(\sigma\) and the Brillouin-zone integration; mixing and the convergence criterion |
| `structure` | Lattice, lattice constant (Bohr and Å), decoration cell or primitive vectors, sphere radius, volume per site, the sublattices with their species and concentrations, or the atoms with positions (k-space) |
| `short-range order` | Only with `[[sro]]`: the conditional probabilities \(P(j\mid i)\), the Warren-Cowley \(\alpha_{ij}\), the screening modes (eigenvalue of \(\hat A\), \(\Gamma=\delta q/\delta b\), the first-shell value, \(F\) in Ry/e, "infinite" for a fully clustered mode) and the drive \(u\) |
| `SCF history` | One line per iteration: \(E_{tot}\), its change, \(dV\), the largest change of \(E_\nu\), \(E_F\), the number of \(E_F\) passes and the wall time. Then per component the moment and the net charge per iteration, in blocks of columns. |
| `total energy` | The terms of \(E_{tot}\) (see below) per site, and per component for alloys |
| `converged results` | \(E_F\), the cell charge residual, the moment per site; per component \(Z\), valence, concentration, weight, \(S\), electrons in the sphere, net charge, moment; the potential parameters and the moments with the DOS at \(E_F\) for each spin and \(l\) |
| `Fermi level` | \(N(E_F)\) per site (both spins), its broadening, the bare Sommerfeld coefficient \(\gamma_e\) and the bare Pauli susceptibility \(\chi_P\), and \(N(E_F)\) per component and \(l\) |
| `health checks` | Each check with its value, its limit and a status: `OK`, `WARN`, `FAIL` or `info` (see below) |
| `timing` | Wall and CPU time per part of the run, calls, CPU/wall and share; the time outside the parts; chain evaluations per second |
| `warnings` | Every `WARN` and `FAIL` of the health checks and every warning of the run; `none` if there are none |
| `files written` | The files of this run |
| `end of run` | Finish time, total wall time and one status line |

The status line is one of

```text
 status             : NORMAL TERMINATION   SCF CONVERGED
 status             : NORMAL TERMINATION   SCF CONVERGED WITH WARNINGS
 status             : NORMAL TERMINATION   SCF NOT CONVERGED
```

!!! note "Units and spin in the report"
    Energies in Rydberg, lengths in Bohr (Å in parentheses where stated), charges in
    electrons, moments in μ<sub>B</sub>, DOS in states/Ry. The banner repeats this
    with the conversions. With `spin = 1` the charges and the DOS are sums over both
    spins. With `spin = 2` the tables list each spin channel, and a row `1-2` gives
    spin 1 minus spin 2.

### The total energy

The total energy has the Kohn-Sham form: the output density in the input potential.
The Fermi smearing entropy is not included. The report lists its terms; their sum is
\(E_{tot}\) to rounding.

| Term | Formula |
|---|---|
| valence band energy | \(\sum_l (q_{1,l} + E_{\nu,l} Q_l)\), \(Q_l\) the sphere charge of the channel |
| core eigenvalue sum | sum of the relaxed core eigenvalues |
| double counting | \(-\int n\,(v_{in} + 2Z/r)\) |
| sphere Hartree | \(\tfrac12 \int n\, V_H[n]\) of the electrons |
| exchange-correlation | \(\int n\, \varepsilon_{xc}[n]\) |
| Madelung | \(\tfrac12 q V_M\) (sublattice-average charges; k-space: the cell charges) |
| screened deviations | \(\tfrac12 q V_F\) (screened charge deviations, SIM; with SRO the drive form) |

The nuclear terms \(-2Z\int n/r\) of the kinetic energy and of the sphere
electrostatics cancel. So neither the double counting nor the Hartree term holds them,
and these two large terms nearly cancel each other. For alloys a second table gives
each term per atom of each component and per site
(\(\sum_p w_p E_p\)).

### Reading `converged results`

```text
   potential parameters, gamma representation [Ry]
   spin  l         E_nu            C  sqrt(Delta)            p        gamma
      1  d    -0.260398    -0.252370     0.097030     2.500828    -0.003822

   moments and DOS at E_F, both spins
   spin  l      q0 [e]   q1 [e Ry]  q2 [e Ry^2]     Q_l [e]   N_l(E_F) [1/Ry]
      1  d     9.55236     0.00015      0.07210     9.55236            1.3913
```

| Column | Meaning |
|---|---|
| `E_nu` | Linearisation energy (Ry) |
| `C` | Band-centre potential parameter (Ry) |
| `sqrt(Delta)` | Band-width parameter. Its sign is a pure gauge (observed e.g. for p). |
| `p` | \(\langle\dot\phi|\dot\phi\rangle\), the overlap parameter |
| `gamma` | \(\gamma\) of the gamma representation |
| `q0`, `q1`, `q2` | Energy moments of the \(l\) channel |
| `Q_l` | Sphere charge of the channel: \(q_0\) for the recursion route `edc`, \(q_0 + p\,q_2\) for `ed3` and k-space |
| `N_l(E_F)` | DOS of the channel at \(E_F\) (states/Ry) |

The row `sum` gives the valence electrons in the sphere (\(\sum Q_l\)) and
\(N(E_F)\) of the component. \(C\), \(\sqrt\Delta\), \(p\), \(\gamma\), \(q_0..q_2\)
are in the **gamma representation**. This matters if you compare to other LMTO codes:
a code that stores parameters in a different representation must transform them first
(see [Conventions](appendix/conventions.html#gamma-representation)).

### Derived quantities at the Fermi level

\(\gamma_e = (\pi^2/3)\,k_B^2\,N(E_F)\) and \(\chi_P = \mu_B^2\,N(E_F)\), with
\(N(E_F)\) per site and both spins, and one mole of sites. The report prints the
constants: 0.173247 mJ/(mol K²) and 2.37605×10⁻⁶ emu/mol per state/Ry. Both are bare
band values: they leave out the electron-phonon and the Stoner enhancement.

### Health checks

| Check | Limit | Status when not met |
|---|---|---|
| SCF converged (\(dV\)) | `tolerance` | `FAIL` |
| \(E_\nu\) held at its band window, last iteration / whole run | 0 in the last iteration | `FAIL` |
| non-causal contour values (\(-\mathrm{Im}\,G \le 0\)), recursion | 0 | `WARN` |
| valence weight below the contour bottom, first / final potential, recursion | \(10^{-3}\) e | `WARN` |
| cell charge residual, recursion | \(10^{-3}\) e per site | `WARN` |
| band energy from the moments against the eigenvalue sum, k-space | \(10^{-5}\) Ry | `WARN` |
| recursion depth used, terminator fallbacks, Anderson restarts, \(E_F\) passes | none | `info` |

## The summary `asr.json`

The same final numbers as the report, in JSON, written next to it (the name of
`out_file` with the extension `.json`). Scripts and tables should read this file, not
the text report. Numbers carry the full double precision; units are in the key names.

| Key | Content |
|---|---|
| `schema` | `"asr-summary-1"` |
| `provenance` | version, git revision, dirty flag, compiler, host, MPI ranks, OpenMP threads, input file, input fingerprint |
| `method` | route, levels, recursion steps, rank, spin, XC |
| `status` | `converged`, `iterations`, `warnings` (the list of the report) |
| `fermi_energy_ry`, `total_energy_ry_per_site` | \(E_F\), \(E_{tot}\) |
| `energy_terms_ry_per_site` | the terms of \(E_{tot}\) |
| `components` | per component: label, species, sublattice, \(Z\), concentration, weight, net charge, moment, s/p/d occupations, s/p/d DOS at \(E_F\), energy terms |
| `dos_ef_per_site_per_ry` | \(N(E_F)\) per site |
| `history` | per iteration: `etot_ry`, `dv_ry`, `ef_ry`, `npass`, `wall_s` |
| `checks` | non-causal values, contour values, weight below the contour (first, final), \(E_\nu\) clamps (last, run), cell charge residual |
| `timing_s` | total wall and CPU time, wall time per part |

```python
import json
s = json.load(open("asr.json"))
print(s["status"]["converged"], s["total_energy_ry_per_site"], s["components"][0]["moment_mub"])
```

## `dos_file`

States/Ry on the energy mesh of `[dos]`. The header records the run, the route, the
broadening, \(E_F\), the mesh, the units and the spin convention. Two aligned comment
rows close the header: the column numbers, then a label over each column.

```text
#          1             2             3             4             5             6             7             8
#     E [Ry]         total      Cu@bcc s      Cu@bcc p      Cu@bcc d      Zn@bcc s      Zn@bcc p      Zn@bcc d
   -0.900000  9.890731E-01  1.288636E-01  4.872272E-02  1.680909E-01  3.151507E-01  6.290973E-02  1.254408E+00
```

**`method = "asr"`.** `total` is per site, \(\sum_p w_p N_p\). Then for each component
(`species@sublattice`) the s, p, d DOS per atom, summed over both spins. With
`spin = 2` follow the per-site totals of spin up and spin down, and for each component
the s, p, d DOS of each spin channel.

**`method = "kspace"`.** For each atom and each spin channel the s, p, d DOS
(states/Ry per atom and spin channel), atom by atom.

```python
import numpy as np
d = np.loadtxt("dos.dat")
E, total = d[:, 0], d[:, 1]
```

## `params_file` (`method = "kspace"` only)

Potential parameters in the gamma representation (Ry):

```text
# fcc Ag  (generated by asr lmto; gamma representation, Ry)
# E_F =  -0.17371795
# class symbol  spin l      enu            c          sqrt(Delta)        p            gamma
  1  Ag         1   0    -0.53399071    -0.44941372     0.36490979     0.08176447     0.42492068
```

## Other files

| File | When | Content |
|---|---|---|
| `save_potential` | `method = "asr"` | Species potentials and \(E_\nu\) of the last iteration (written every iteration); read back with `start_potential` |

## Terminal output of SCF iterations

The recursion prints one line per iteration in the terminal:

```text
it 30  dV= 4.021E-06  EF=  0.017089  qnet(A,B)=  0.06709 -0.06709  m(A,B)=  0.0000  0.0000
```

The convergence criterion is the maximum change of \(rV/(1+r)\) between
iterations (below `tolerance`, in Ry): the weight \(r/(1+r)\) keeps the
nucleus region \(r\sim10^{-5}\) Bohr from dominating the measure.
