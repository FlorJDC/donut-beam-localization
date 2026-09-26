# -*- coding: utf-8 -*-
import os, sys, time
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), *[".."] * 5))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from donutloc import experiments as ex, photons, patterns, beams, fisher, estimators as est

t = time.time()
r = ex.iterative_minflux(1000, n_rep=2000, seed=7)
print("iter fixed noBG n_rep=2000 seed7: sigma %.4f±%.4f rmse %.4f bias %s cam %.4f (%.1fs)" % (r["sigma"], r["sigma_se"], r["rmse"], r["bias"], r["camera_sigma"], time.time() - t))
print("  sigma_iter", r["sigma_iter"], "crb_iter", r["crb_center_iter"])
# outliers per iteration: max |error| relative to L_k
rt = r["r_true"]
# replicate internally to inspect iteration-0 error distribution
beam = beams.make_beam("donut", fwhm=300.)
rng = np.random.default_rng(7)
rad = 37.5 * np.sqrt(rng.uniform(size=2000)); phi = rng.uniform(0, 2*np.pi, 2000)
rt2 = np.stack([rad*np.cos(phi), rad*np.sin(phi)], 1)
print("  r_true reproduced:", np.allclose(rt2, rt))
p0 = photons.make_model(patterns.tcp_centers(150.), beam)
C = rng.multinomial(250, p0(rt2) / p0(rt2).sum(1, keepdims=True))
e0 = est.mle(C, p0, search_radius=112.5) - rt2
d = np.hypot(*e0.T)
print("  iter0 |err| quantiles 50/99/max: %.2f %.2f %.2f ; frac>12.5nm (outside next L/2): %.4f" % (np.median(d), np.quantile(d, .99), d.max(), np.mean(d > 12.5)))
# expected sigma iter0 from CRB averaged over disk
pts = rt2[:500]
cr = fisher.crb(p0, pts, 250)
print("  rms CRB over emitter disk (iter0): %.3f" % np.sqrt(np.mean(cr**2)))
# final-iteration efficiency: CRB at actual relative positions in last iteration
# SBR=10 and adaptive
for kw in ({"sbr": 10.}, {"rule": "adaptive"}, {"rule": "adaptive", "sbr": 10.}):
    r = ex.iterative_minflux(1000, n_rep=2000, seed=7, **kw)
    print(kw, "sigma %.4f±%.4f rmse %.4f L %s" % (r["sigma"], r["sigma_se"], r["rmse"], np.round(r["L"], 2)))
# SE check: bootstrap of sigma
r = ex.iterative_minflux(1000, n_rep=2000, seed=7)
e = r["estimates"] - r["r_true"]
bs = []
g = np.random.default_rng(3)
for _ in range(400):
    i = g.integers(0, len(e), len(e)); bs.append(np.sqrt(0.5 * e[i].var(0, ddof=1).sum()))
print("  sigma_se formula %.4f vs bootstrap %.4f; kurtosis-x %.2f" % (r["sigma_se"], np.std(bs), np.mean((e[:,0]-e[:,0].mean())**4)/e[:,0].var()**2))
# N_total small edge: split remainder
print("split 1003/4:", ex._split_photons(1003, 4, "equal"), " weights:", ex._split_photons(1000, 4, [1, 1, 1, 5]))
# misalignment: honest uses perturbed centers? replicate one pattern manually
t = time.time()
m = ex.misalignment_study([0., 5.], n_patterns=3, n_rep=100)
print("misalign small (%.1fs): honest bias" % (time.time() - t), m["honest"]["bias_abs_mean"], "naive", m["naive"]["bias_abs_mean"], "identical", m["identical"])
rng = np.random.default_rng(42)
ideal = patterns.tcp_centers(100.)
tc = patterns.perturb_centers(ideal, 5., rng=rng)
print("  perturbed centres pattern0:", np.round(tc - ideal, 3))
# perturbed centre-exposure too? (4 zeros all moved) -> hypot:
print("  displacement norms:", np.round(np.hypot(*(tc - ideal).T), 6))
# time for figure-scale misalignment: 4 deltas x 20 patterns x 2 pos x 200 reps -> reported 11 s. For <2% SE on sigma need R>=625
t = time.time(); ex.misalignment_study([5.], n_patterns=5, n_rep=625); print("  1 delta x 5 patterns x 625 reps: %.1fs" % (time.time() - t))
