# -*- coding: utf-8 -*-
"""Round-1 code-review integration checks: make_model p_fn -> fisher / estimators / montecarlo /
closed_forms. Run from project root: python equipo/2026-09-26_donut-localization/work/review/r01_integration.py
"""
import os, sys, time
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "tests"))
import numpy as np
from donutloc import beams, patterns, photons, fisher, estimators as est, montecarlo as mc, closed_forms as cf
import test_acceptance as acc   # read-only import of the independent reference

def model(L=50.0, fwhm=300.0, **kw):
    bkw = {k: kw.pop(k) for k in ("eps", "zero_model") if k in kw}
    return photons.make_model(patterns.tcp_centers(L), beams.make_beam("donut", fwhm=fwhm, **bkw), **kw)

print("A. crb_limit via make_model vs acceptance _crb_center vs closed form")
for L, f in ((50., 300.), (100., 300.), (50., 360.), (150., 250.)):
    p = model(L, f)
    a = fisher.crb_limit(p, 100); b = acc._crb_center(L, 100, f); c = cf.crb_tcp_center_limit(L, 100, f)
    pt = float(fisher.crb(p, [0., 0.], 100)); s27 = cf.crb_tcp_center_point(L, 100, f)
    print("  L=%g f=%g lim pkg=%.7f acc=%.7f cf=%.7f | pt pkg=%.7f S27=%.7f  rel(pkg/acc-1)=%.1e rel(pkg/cf-1)=%.1e"
          % (L, f, a, b, c, pt, s27, a / b - 1, a / c - 1))

print("B. eps (both models) + sbr through photons (beam convention) vs cf.crb_tcp_center_eps")
for zm in ("constant", "gaussian"):
    for eps in (0.002, 0.05):
        for sbr in (None, 5.0):
            p = model(100., 300., eps=eps, zero_model=zm, sbr=sbr)
            got = float(fisher.crb(p, [0., 0.], 100))
            ref = cf.crb_tcp_center_eps(100., 100, 300., eps, np.inf if sbr is None else sbr, zm, "beam")
            lim = fisher.crb_limit(p, 100)
            print("  %s eps=%g sbr=%s pt=%.6f cf=%.6f rel=%.1e  lim/pt-1=%.1e" % (zm, eps, sbr, got, ref, got / ref - 1, lim / got - 1))

print("C. lms (general, make_model) vs lms_tcp, expected counts at r=(0.5,-0.3)")
for sbr in (None, 10.0):
    p = model(50., 300., sbr=sbr)
    C = p(np.array([0.5, -0.3])) * 1000
    print("  sbr=%s lms=%s lms_tcp=%s" % (sbr, est.lms(C, p), est.lms_tcp(C, 50., 300., sbr=sbr)))

print("D. MLE Monte Carlo at centre through make_model (L=50, N=100, fwhm=300)")
for sbr in (10.0, None):
    p = model(50., 300., sbr=sbr)
    t = time.time()
    res = mc.run_mc(lambda Cc: est.mle(Cc, p, search_radius=50.), p, [0., 0.], 100, 5000)
    dt = time.time() - t
    crb = fisher.crb_limit(p, 100)
    print("  sbr=%s sigma=%.4f +- %.4f crb_lim=%.5f eff=%.4f +- %.4f bias=%s  (%.2f s)"
          % (sbr, res["sigma"], res["sigma_err"], crb, res["sigma"] / crb, res["sigma_err"] / crb, res["bias"], dt))
    # empirical SE of sigma by bootstrap vs reported sigma_err
    e = res["estimates"]; rng = np.random.default_rng(1)
    bs = [np.sqrt(0.5 * e[rng.integers(0, len(e), len(e))].var(0, ddof=1).sum()) for _ in range(300)]
    print("     bootstrap SE(sigma)=%.4f vs reported sigma_err=%.4f" % (np.std(bs), res["sigma_err"]))

print("E. poisson mode through make_model, all-zero counts at low N")
p = model(50., 300.)
res = mc.run_mc(lambda Cc: est.mle(Cc, p, search_radius=50.), p, [5., 0.], 2, 2000, mode="poisson")
print("  N=2 poisson: n_valid=%d, zero-total reps=%d, bias=%s" % (res["n_valid"], int((res["counts"].sum(1) == 0).sum()), res["bias"]))

print("F. fisher.crb default p_min near the zero (pointwise discontinuity radius)")
for r in (1e-2, 1e-3, 1e-4, 5e-5, 1e-5):
    print("  r=%g crb=%.6f" % (r, float(fisher.crb(p, [r, 0.], 100))))

print("G. far field with gaussian beam (underflow)")
pg = photons.make_model(patterns.tcp_centers(50.), beams.make_beam("gaussian", fwhm=300.))
print("  crb at r=6000 nm:", fisher.crb(pg, [6000., 0.], 100), " p:", pg([6000., 0.]))

print("H. mle grid size for fine grid_step (memory of chunk @ grid)")
for gs in (2.0, 0.5, 0.2):
    G = est._disk_grid((0., 0.), 50., gs).shape[0]
    print("  grid_step=%g G=%d -> chunk(4096) x G float64 = %.2f GB" % (gs, G, 4096 * G * 8 / 1e9))

print("I. mle with off-origin search centre (pattern search disk)")
p = model(100., 300., sbr=10.)
truth = np.array([30., 10.])
print("  ", est.mle(np.round(1e7 * p(truth)), p, search_radius=10., center=(25., 5.)), "(truth", truth, ", disk r=10 around (25,5))")
