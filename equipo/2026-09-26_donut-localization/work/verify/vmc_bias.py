# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import minimize
from scipy.special import xlogy
from vmc import run, pgrid, summ
# larger run for bias at (2,0)
ds = []
for seed in [11, 12, 13]:
    est, n = run([2, 0], 100, np.inf, 8000, seed)
    ds.append(est[:, 0]-2)
d = np.concatenate(ds); print('bias dx (24000 reps)', d.mean(), '+-', d.std()/np.sqrt(len(d)))
# cross-check the grid MLE against Nelder-Mead from multiple starts on 300 samples
est, n = run([2, 0], 100, np.inf, 300, 5)
worse = 0; maxdiff = 0
for k in range(300):
    f = lambda r: -xlogy(n[k], pgrid(np.array([r]), np.inf)[0]).sum()
    best = min((minimize(f, s0, method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-10)) for s0 in [est[k], [0, 0], [5, 5], [-5, 3], [2, -6]]), key=lambda o: o.fun)
    if best.fun < f(est[k])-1e-7: worse += 1
    maxdiff = max(maxdiff, np.abs(best.x-est[k]).max())
print('grid worse than NM in', worse, 'of 300; max |diff|', maxdiff)
