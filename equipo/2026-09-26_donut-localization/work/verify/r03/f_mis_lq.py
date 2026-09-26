# -*- coding: utf-8 -*-
"""Follow-up on the (L/4,0) naive-bias discrepancy: (a) noise-free naive MLE on expected counts over
many patterns (population mean); (b) second independent 400x200 MC draw, naive only, at (L/4,0)."""
import os, sys, time
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from mymle import tcp, pvec, mle
L, N, SBR = 100., 500, 10.; C0 = tcp(L)
for pn, r0 in (('center', np.zeros(2)), ('Lq', np.array([25., 0.]))):
    rng = np.random.default_rng(8080); dirs = 2*np.pi*rng.random((4000, 4))
    line = []
    for d in (2., 5., 10.):
        n = np.array([N*pvec(r0, C0 + d*np.stack([np.cos(q), np.sin(q)], 1), SBR) for q in dirs])
        e = mle(n, C0, SBR, 0.75*L) - r0; b = np.linalg.norm(e, axis=1)
        line.append('d=%g: <|b0|>/d=%.4f±%.4f' % (d, b.mean()/d, b.std()/np.sqrt(len(b))/d))
    print('noise-free naive bias, %s, 4000 patterns: %s' % (pn, '  '.join(line)), flush=True)
P, R = 400, 200
dirs = 2*np.pi*np.random.default_rng(271828).random((P, 4)); rng = np.random.default_rng(161803)
r0 = np.array([25., 0.]); out = []
for d in (2., 5., 10.):
    b = []
    for j in range(P):
        C = C0 + d*np.stack([np.cos(dirs[j]), np.sin(dirs[j])], 1)
        n = rng.multinomial(N, pvec(r0, C, SBR), size=R)
        b.append(np.linalg.norm((mle(n, C0, SBR, 0.75*L) - r0).mean(0)))
    b = np.array(b); out.append((d, b.mean(), b.std(ddof=1)/np.sqrt(P)))
    print('MC draw 2, Lq naive d=%g: |bias|=%.3f±%.3f  ratio %.4f' % (d, b.mean(), out[-1][2], b.mean()/d), flush=True)
dd = np.array([o[0] for o in out]); bb = np.array([o[1] for o in out]); ss = np.array([o[2] for o in out]); w = 1/ss**2
print('draw 2 slope Lq = %.4f ± %.4f' % ((w*dd*bb).sum()/(w*dd*dd).sum(), 1/np.sqrt((w*dd*dd).sum())))
