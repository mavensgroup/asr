"""Weyl-disc criterion: W = -Im [[g11, g1m],[g1m, gmm]] of Ktilde_m^-1 positive definite  <=>  every tail with Im t <= 0 gives -Im G > 0."""
import numpy as np, sys
from f4lab import *
rng = np.random.default_rng(7)
def corner(al, be, m):
    Ki = np.linalg.inv(tri(al, be, m))
    return Ki[0, 0], Ki[0, m - 1], Ki[m - 1, m - 1]
tot = dict(n=0, wbad=0, neg0=0, negh=0, negr=0, viol_suff=0, viol_nec=0, negh_wok=0, negr_wok=0, thull_noncausal=0)
for tv, seed in (((1.0, 0.25, 3.0), 1), ((1.0, 0.25, 3.0), 2), ((1.0, 0.4, 2.0), 3), ((1.0, 0.6, 1.5), 2)):
    A, M, H, T = make_pencil(L=7, seed=seed, tvals=tv, eps=(0.0, 0.15, -0.2))
    for y in (0.1, 0.02, 0.005, 0.001):
        for m in (24,):
            loc = dict(n=0, wbad=0, neg0=0, negh=0, negr=0)
            for E in np.linspace(-0.9, 0.9, 401):
                z = E + 1j * y
                al, be, Q = cs_lanczos(z * M - A, m)
                g11, g1m, gmm = corner(al, be, m)
                W = -np.array([[g11.imag, g1m.imag], [g1m.imag, gmm.imag]])
                wok = np.linalg.eigvalsh(W).min() > 0
                th = hull_tail_K(al, be, m, z); tr = hull_tail_K(al, be, m, z, real_edges=True)
                G0 = g11; Gh = g_chain(al, be, m, th); Gr = g_chain(al, be, m, tr)
                # Moebius form check
                assert abs(Gh - (g11 + g1m ** 2 * th / (1 - gmm * th))) < 1e-8 * (1 + abs(Gh))
                for k in ("n",): tot[k] += 1; loc[k] += 1
                if not wok: tot["wbad"] += 1; loc["wbad"] += 1
                if G0.imag >= 0: tot["neg0"] += 1; loc["neg0"] += 1
                if Gh.imag >= 0: tot["negh"] += 1; loc["negh"] += 1; tot["negh_wok"] += wok
                if Gr.imag >= 0: tot["negr"] += 1; loc["negr"] += 1; tot["negr_wok"] += (wok and tr.imag <= 0)
                if th.imag > 0: tot["thull_noncausal"] += 1
                # sufficiency: random causal tails
                ts = (rng.normal(size=200) - 1j * np.abs(rng.normal(size=200))) * 10.0 ** rng.uniform(-4, 2, size=200)
                Gt = g11 + g1m ** 2 * ts / (1 - gmm * ts)
                if wok and (Gt.imag >= 0).any(): tot["viol_suff"] += 1
                # necessity: if W not PD there is a causal tail with -Im G <= 0 (boundary circle: real t)
                if not wok:
                    tt = np.concatenate([np.linspace(-50, 50, 20001) / max(abs(gmm), 1e-30), [0.0]])
                    tt = np.concatenate([tt, np.tan(np.linspace(-1.5707, 1.5707, 20001)) / max(abs(gmm), 1e-30)])
                    Gt = g11 + g1m ** 2 * tt / (1 - gmm * tt)
                    if not (Gt.imag >= -1e-12 * abs(Gt)).any(): tot["viol_nec"] += 1
            print("T %s seed %d Im z %.3f m %d: W not PD %3d, -ImG<=0 tail-free %2d hull %2d real-edge %2d  of %d" % (tv, seed, y, m, loc["wbad"], loc["neg0"], loc["negh"], loc["negr"], loc["n"]))
            sys.stdout.flush()
print(tot)
