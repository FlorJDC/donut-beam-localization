# -*- coding: utf-8 -*-
"""r06 verifier MC: fixed SBR vs background constant per exposure (fig 5 points, iterative),
plus a large-N MC check that the naive-MLE bias is the N-independent pseudo-true shift.
Own seed 20260930; no donutloc."""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcore import tcp, intens, probs, crb, mle_batch, lms_tcp, mlms_tcp, sigma

HERE = os.path.dirname(os.path.abspath(__file__))
FW = 300.0
SEED = 20260930
which = sys.argv[1] if len(sys.argv) > 1 else "all"
out = {}


def boot_se(err, nb=500, seed=1):
    rng = np.random.default_rng(seed)
    n = err.shape[0]
    return np.std([sigma(err[rng.integers(0, n, n)]) for _ in range(nb)], ddof=1)


if which in ("all", "fig5"):
    L, N, REP = 50.0, 100, 20000
    cen = tcp(L)
    bg = intens(np.zeros(2), cen, FW).sum() / 40.0
    models = {"fixed": lambda q: probs(q, cen, FW, sbr=10.0),
              "phys": lambda q: probs(q, cen, FW, bg=bg)}
    for x in (25.0, 50.0):
        r = np.array([x, 0.0])
        for mk, pf in models.items():
            rng = np.random.default_rng(SEED + int(x) + (0 if mk == "fixed" else 1000))
            p = pf(r)
            cnt = rng.multinomial(N, p / p.sum(), size=REP).astype(float)
            c = crb(pf, r, N)
            t0 = time.time()
            ests = {"mle": mle_batch(cnt, pf, 2 * L, h0=1.0),
                    "lms": lms_tcp(cnt, L, FW, 10.0),
                    "mlms": mlms_tcp(cnt, L, FW, 10.0)}
            for ek, e in ests.items():
                err = e - r
                k = "%s_%s_x%d" % (mk, ek, x)
                out[k + "_bias_x"] = float(err[:, 0].mean())
                out[k + "_bias_x_se"] = float(err[:, 0].std(ddof=1) / np.sqrt(REP))
                out[k + "_bias_y"] = float(err[:, 1].mean())
                out[k + "_bias_y_se"] = float(err[:, 1].std(ddof=1) / np.sqrt(REP))
                s = sigma(err)
                out[k + "_sigma_over_crb"] = float(s / c)
                out[k + "_sigma_over_crb_se"] = float(boot_se(err) / c)
            print(x, mk, "done %.1fs" % (time.time() - t0), flush=True)

if which in ("all", "iter"):
    NT, REP = 1000, 10000
    Ls = 150.0 * (25.0 / 150.0) ** (np.arange(4) / 3.0)
    Nk = [250] * 4

    def S0(Lk):
        return intens(np.zeros(2), tcp(Lk), FW).sum()

    runs = {"matchL150": dict(bg=S0(150.0) / 40.0), "matchL25": dict(bg=S0(25.0) / 40.0),
            "fixedSBR10": dict(sbr=10.0), "nobg": dict()}
    for irun, (name, kw) in enumerate(runs.items()):
        rng = np.random.default_rng(SEED + 17 * (irun + 1))
        rad = 37.5 * np.sqrt(rng.uniform(size=REP))
        ph = rng.uniform(0, 2 * np.pi, REP)
        rt = np.stack([rad * np.cos(ph), rad * np.sin(ph)], -1)
        c = np.zeros_like(rt)
        t0 = time.time()
        for k, Lk in enumerate(Ls):
            cen = tcp(Lk)
            pf = lambda q, cen=cen: probs(q, cen, FW, **kw)
            P = pf(rt - c)
            P = P / P.sum(1, keepdims=True)
            cnt = np.array([rng.multinomial(Nk[k], pp) for pp in P], float)
            rel = mle_batch(cnt, pf, 0.75 * Lk, h0=0.75 * Lk / 80)
            c = c + rel
        err = c - rt
        out["iter_%s_sigma" % name] = float(sigma(err))
        out["iter_%s_sigma_se" % name] = float(boot_se(err))
        out["iter_%s_rmse" % name] = float(np.sqrt(np.mean(np.sum(err ** 2, 1)) / 2))
        print(name, out["iter_%s_sigma" % name], "%.1fs" % (time.time() - t0), flush=True)

if which in ("all", "largeN"):
    # naive MLE bias at large N approaches the noise-free (pseudo-true) value
    cen = tcp(50.0)
    r = np.array([20.0, 0.0])
    pt = lambda q: probs(q, cen, FW, eps=0.01, sbr=10.0)
    pn = lambda q: probs(q, cen, FW, sbr=10.0)
    for N in (2000, 20000):
        rng = np.random.default_rng(SEED + N)
        cnt = rng.multinomial(N, pt(r), size=4000).astype(float)
        e = mle_batch(cnt, pn, 37.5, h0=0.25) - r
        out["largeN_naive_eps0p01_x20_N%d_bias" % N] = [float(e[:, 0].mean()), float(e[:, 1].mean())]
        out["largeN_naive_eps0p01_x20_N%d_se" % N] = [float(e[:, 0].std() / 63.2), float(e[:, 1].std() / 63.2)]
        print(N, out["largeN_naive_eps0p01_x20_N%d_bias" % N], flush=True)

json.dump(out, open(os.path.join(HERE, "mc_%s.json" % which), "w"), indent=1)
for k, v in out.items():
    print(k, v)
