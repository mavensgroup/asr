#!/usr/bin/env python3
"""Gate B4.4: Warren-Cowley triads of an A-B-C alloy against Alam & Mookerjee, JPCM 21, 195503 (2009), figures 4 and 5.

Reads <dir>/rand/dos.dat and <dir>/tXYZ/dos.dat, where XYZ = (alpha_FeNi, alpha_FeCr, alpha_NiCr) in {0, 1} (our run names).
Alam labels a triad (alpha_AB, alpha_AC, alpha_BC); the paper's text and table 1 give (A, B, C) = (Fe, Ni, Cr).  The energy moment
functions are the paper's, M_n(E) = int_{-inf}^E dE' E'^n n(E') of the total DOS, energies absolute (Ry).  Prints
  - the L1 change of the total and component d DOS against the random alloy (which pair controls the DOS shape),
  - for every assignment of (A, B, C) to the species: the fraction of energies in [-0.6, 0.3] Ry where the paper's orderings
    M2(110) < M2(101) < M2(011) and M2(100) < M2(010) < M2(001) hold, and whether the paper's figure 4 families (alpha_BC = 1 vs 0)
    separate in the size of the DOS change,
and writes a figure (total DOS of each triad over the random one; M2 and M4 differences to the random alloy).

  python3 -I tools/sro_triads.py <dir> [--perm Fe Ni Cr] [--out fig.png]
"""
import argparse
import itertools

import numpy as np

SP = ("Fe", "Ni", "Cr")
PAIRS = (("Fe", "Ni"), ("Fe", "Cr"), ("Ni", "Cr"))            # order of our run names tXYZ
CLAIMS = (("110", "101", "011"), ("100", "010", "001"))        # paper: M2 increasing along each triple
FAM1, FAM0 = ("101", "011", "001", "111"), ("110", "100", "010")  # paper figure 4: alpha_BC = 1 / 0 shape families
RUNS = ["t110", "t101", "t011", "t100", "t010", "t001", "t111"]


def ours(lab, perm):
    """Our run name for the paper's label (alpha_AB, alpha_AC, alpha_BC) when (A, B, C) = perm."""
    a = dict(zip([frozenset(perm[:2]), frozenset((perm[0], perm[2])), frozenset(perm[1:])], lab))
    return "t" + "".join(a[frozenset(p)] for p in PAIRS)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dir")
    ap.add_argument("--perm", nargs=3, default=list(SP), help="(A, B, C) used for the figure labels")
    ap.add_argument("--out", default="sro_triads.png")
    arg = ap.parse_args()
    dat = {r: np.loadtxt(f"{arg.dir}/{r}/dos.dat") for r in ["rand"] + RUNS}
    E = dat["rand"][:, 0]
    dE = E[1] - E[0]
    grid = (E >= -0.6) & (E <= 0.3)

    def cum(y):
        return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * dE)])

    def mom(r, n):
        y = dat[r][:, 1] / cum(dat[r][:, 1])[-1]                # normalised to 1 over the DOS window
        return cum(E**n * y)

    print("L1 change against the random alloy (states/atom): total DOS, then d DOS of Fe, Ni, Cr")
    l1 = {}
    for r in RUNS:
        l1[r] = np.abs(dat[r][:, 1] - dat["rand"][:, 1]).sum() * dE
        dd = [np.abs(dat[r][:, k] - dat["rand"][:, k]).sum() * dE for k in (4, 7, 10)]
        print(f"  {r}  {l1[r]:6.3f}   " + " ".join(f"{x:6.3f}" for x in dd))
    print("Paper's M2 orderings (fraction of energies in [-0.6, 0.3] Ry) and figure 4 families, per (A, B, C):")
    for perm in itertools.permutations(SP):
        fr = []
        for c in CLAIMS:
            m = [mom(ours(l, perm), 2) for l in c]
            fr.append(np.mean(((m[0] < m[1]) & (m[1] < m[2]))[grid]))
        sep = min(l1[ours(l, perm)] for l in FAM1) > max(l1[ours(l, perm)] for l in FAM0)
        print(f"  ({', '.join(perm)}): {fr[0]:.2f} {fr[1]:.2f}   families separate: {sep}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    blue, orange, aqua, gray = "#2a78d6", "#eb6834", "#1baf7a", "#8a8a85"
    ef = float(open(f"{arg.dir}/rand/dos.dat").readlines()[1].split()[-1])
    perm = tuple(arg.perm)
    inv = {ours(l, perm): l for l in ("".join(t) for t in itertools.product("01", repeat=3)) if l != "000"}
    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(3, 4, hspace=0.45, wspace=0.28)
    for i, r in enumerate(RUNS):
        ax = fig.add_subplot(gs[i // 4, i % 4])
        ax.plot(E, dat["rand"][:, 1], color=gray, lw=1.5, label="random")
        ax.plot(E, dat[r][:, 1], color=blue, lw=2, label="SRO")
        ax.axvline(ef, color=gray, lw=0.8, ls=":")
        ax.set_title(f"ours {r[1:]}   paper {inv[r][0]}-{inv[r][1]}-{inv[r][2]}", fontsize=10)
        ax.set_xlim(-0.65, 0.3)
        ax.tick_params(labelsize=8)
        if i % 4 == 0:
            ax.set_ylabel("total DOS (states/Ry/atom)", fontsize=9)
        if i == 0:
            ax.legend(fontsize=8, frameon=False)
    ax = fig.add_subplot(gs[1, 3])
    ax.axis("off")
    ax.text(0, 0.9, f"(A, B, C) = ({', '.join(perm)})\nlabels (alpha_AB, alpha_AC, alpha_BC)\n"
            "ours tXYZ = (FeNi, FeCr, NiCr)\ndotted: E_F of the random alloy", fontsize=9, va="top")
    for j, (n, nm) in enumerate(((2, "M2"), (4, "M4"))):
        for k, c in enumerate(CLAIMS):
            ax = fig.add_subplot(gs[2, 2 * j + k])
            base = mom("rand", n)
            for l, col in zip(c, (blue, orange, aqua)):
                ax.plot(E, (mom(ours(l, perm), n) - base) * 1e3, color=col, lw=2, label=f"{l[0]}-{l[1]}-{l[2]}")
            ax.axhline(0, color=gray, lw=0.8)
            ax.set_xlim(-0.6, 0.3)
            ax.set_title(f"{nm}(E) - {nm} random  (x 1e-3)", fontsize=10)
            ax.tick_params(labelsize=8)
            ax.legend(fontsize=8, frameon=False)
            ax.set_xlabel("E (Ry)", fontsize=9)
    fig.savefig(arg.out, dpi=110, bbox_inches="tight")
    print("wrote", arg.out)


if __name__ == "__main__":
    main()
