# -*- coding: utf-8 -*-
import numpy as np
from w13 import pv, mle, C0, N
rng = np.random.default_rng(7)
for d in (2., 5., 10.):
    for which in ('all4', 'outer3'):
        B = []
        for s in range(100):
            phi = 2*np.pi*rng.random((20, 4)); D = d*np.stack([np.cos(phi), np.sin(phi)], -1)
            if which == 'outer3': D[:, 3] = 0
            ns = np.array([N*pv(np.zeros(2), C0+D[k]) for k in range(20)])
            e = mle(ns, C0)          # asymptotic (noise-free) bias per pattern
            B.append(np.linalg.norm(e, axis=1).mean())
        B = np.array(B)
        print('delta', d, which, 'mean-over-20 bias: mean %.2f sd %.2f  [5%%,95%%]=[%.2f,%.2f]' % (B.mean(), B.std(), *np.percentile(B, [5, 95])))
