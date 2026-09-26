# -*- coding: utf-8 -*-
import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), *[".."] * 5))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from donutloc import fisher, photons, patterns, beams, estimators as est
p = photons.make_model(patterns.tcp_centers(50.), beams.make_beam("donut", fwhm=300.))
N = np.array([100., 400.])
print("point, r=(2,), N=(2,):", fisher.crb(p, [0., 0.], N))
try:
    print("limit, r=(2,), N=(2,):", fisher.crb(p, [0., 0.], N, zero_policy="limit"))
except Exception as e:
    print("limit, r=(2,), N=(2,): EXC", type(e).__name__, e)
try:
    print("crb_map N array (1,5):", fisher.crb_map(p, np.linspace(-2, 2, 5), [0.], np.full((1, 5), 100.)))
except Exception as e:
    print("crb_map N (1,5): EXC", type(e).__name__, e)
# mle with integer-like float counts from Poisson mode OK; center not origin with p_fn check
print("mle center offset:", est.mle([30, 30, 30, 2], p, 20., center=(1., 1.)))
# chunk float
print("mle chunk=2.5:", est.mle([[30, 30, 30, 2]]*3, p, 20., chunk=2.5))
