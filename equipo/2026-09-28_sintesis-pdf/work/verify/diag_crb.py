# -*- coding: utf-8 -*-
import sys, numpy as np
sys.path.insert(0, '/home/user/donut-beam-localization/src')
exec(open('verify_pminflux_B.py').read().split('# --- MLE propio')[0].split('print("== 2. CRB ratios")')[0])
from donutloc import fisher, beams, patterns, photons
import inspect
src = inspect.getsource(fisher.fisher_matrix); print(src[-900:])
pts = np.c_[xd, yd]
ideal = lambda r: p_of(np.asarray(r)[..., 0], np.asarray(r)[..., 1])
theirs = fisher.crb(ideal, pts, 1.0, zero_policy="limit")
d = theirs - c0; i = np.argsort(-np.abs(d))[:8]
print(theirs.mean(), c0.mean()); print(np.c_[xd[i], yd[i], c0[i], theirs[i]])
