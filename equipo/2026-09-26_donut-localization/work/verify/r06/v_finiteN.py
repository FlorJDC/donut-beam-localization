# -*- coding: utf-8 -*-
"""Finite-N MC mean bias of naive vs honest MLE (eps=0.01, x=20, L=50, SBR=10), vs noise-free shift."""
import numpy as np, sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcore import tcp, probs, mle_batch
cen = tcp(50.0); FW = 300.0; r = np.array([20.0, 0.0])
pt = lambda q: probs(q, cen, FW, eps=0.01, sbr=10.0)
pn = lambda q: probs(q, cen, FW, sbr=10.0)
out = {}
for N in (100, 500, 2000, 20000):
    rng = np.random.default_rng(777 + N)
    cnt = rng.multinomial(N, pt(r), size=8000).astype(float)
    for nm, pf in (("naive", pn), ("honest", pt)):
        e = mle_batch(cnt, pf, 37.5, h0=0.25) - r
        m = e.mean(0); se = e.std(0, ddof=1) / np.sqrt(len(e))
        out["N%d_%s" % (N, nm)] = [float(m[0]), float(m[1]), float(se[0]), float(se[1]), float(np.hypot(*m))]
        print(N, nm, np.round(m, 3), np.round(se, 3), "|b|=%.3f" % np.hypot(*m), flush=True)
json.dump(out, open("finiteN.json", "w"), indent=1)
