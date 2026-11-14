---
title: SQS cross-checks
---

# SQS cross-checks

A **special quasirandom structure** reproduces the correlation functions of a
random alloy in a small periodic supercell. It is useful as an *independent*
check of an augmented-space result (a different approximation to the same
physics) and for local environment effects.

[`tools/sqs.py`](https://github.com/mavensgroup/asr/blob/main/tools/sqs.py)
searches supercells (Hermite normal forms of determinant \(N\) of the primitive
cell, compact shapes only) and, in each, the occupation by simulated annealing
with pair swaps that best matches

\[
\Pi_k=\langle\sigma_i\sigma_j\rangle_{\text{shell }k}=(2x-1)^2 ,\qquad
\Pi_3=(2x-1)^3\ \text{(smallest triplet)},
\]

with \(\sigma=+1\) for A, \(-1\) for B. It writes one TOML fragment per SQS
with `[structure] lattice = "cell"` and the `[[site]]` list for
`method = "kspace"`.

```bash
python3 -I tools/sqs.py --lattice bcc --n 16 --x 0.5 --species Fe Co --nsqs 2 --out sqs_dir
```

**Partial (sublattice) disorder.** Build the SQS on the disordered sublattice
only and decorate every site with the ordered atoms at fixed offsets (in
cubic-lattice-constant units), e.g. L2\(_1\) Co\(_2\)(Fe,Cr)Al with Fe/Cr as the fcc SQS
lattice:

```bash
python3 -I tools/sqs.py --lattice fcc --n 10 --x 0.4 --species Fe Cr \
    --decorate Al:0.5,0.5,0.5 Co:0.25,0.25,0.25 Co:-0.25,-0.25,-0.25 --out d
```

Combine the generated fragment with `[[species]]`, `[scf]`, `[dos]`, `[output]`
and set `site_classes = true` to get local moments and charges per site.
Requires Python 3 with NumPy.
