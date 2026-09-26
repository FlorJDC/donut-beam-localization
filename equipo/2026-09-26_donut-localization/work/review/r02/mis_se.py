# -*- coding: utf-8 -*-
import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), *[".."] * 5))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from donutloc import experiments as ex, photons, patterns, beams, montecarlo as mc, estimators as est
beam = beams.make_beam("donut", fwhm=300.)
ideal = patterns.tcp_centers(100.)
pn = photons.make_model(ideal, beam, sbr=10)
for d in (10.,):
    rng = np.random.default_rng(42)
    sig = {"h": [], "n": []}
    for j in range(20):
        tc = patterns.perturb_centers(ideal, d, rng=rng)
        pt = photons.make_model(tc, beam, sbr=10)
        for key, pf in (("h", pt), ("n", pn)):
            o = mc.run_mc(lambda C: est.mle(C, pf, search_radius=75.), pt, [0., 0.], 500, 200, seed=1000 + j)
            sig[key].append(o["sigma"])
    for key in sig:
        s = np.array(sig[key])
        print("delta %g %s: sigma_mean %.4f formula SE %.4f  between-pattern SE %.4f (std %.4f)" % (
            d, key, s.mean(), s.mean() / (2 * np.sqrt(200 * 20)), s.std(ddof=1) / np.sqrt(20), s.std(ddof=1)))
