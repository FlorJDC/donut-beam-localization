# -*- coding: utf-8 -*-
"""R2 review: re-run the four refuted R1 scenarios + unclear items."""
import os, sys, time, warnings, tracemalloc
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), *[".."] * 5))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from donutloc import beams, patterns, photons, fisher, estimators as est, closed_forms as cf

def model(L=50., fwhm=300., sbr=None, kind="donut"):
    return photons.make_model(patterns.tcp_centers(L), beams.make_beam(kind, fwhm=fwhm), sbr=sbr)

p = model()
print("1) crb_map default:", fisher.crb_map(p, np.linspace(-2, 2, 5), [0.0], 100))
print("   crb_map point  :", fisher.crb_map(p, np.linspace(-2, 2, 5), [0.0], 100, zero_policy="point"))
# fine scan near origin: continuity of 'limit' map
xs = np.array([-1e-3, -1e-5, -1e-7, 0, 1e-7, 1e-5, 1e-4, 1e-3, 1e-2])
print("   near-origin map:", fisher.crb_map(p, xs, [0.0], 100))
# 2D map with origin + peripheral zero pixel exactly
m = fisher.crb_map(p, np.arange(-30, 31, 5.0), np.arange(-30, 31, 5.0), 100)
print("   2D map finite:", np.all(np.isfinite(m)), "min", m.min(), "at", np.unravel_index(m.argmin(), m.shape), "val(0,25)", m[list(np.arange(-30,31,5.)).index(25.), 6])
# 2) mle tol<=0
for kw in ({"tol": 0.0}, {"tol": -1.}, {"grid_step": 0.}):
    try:
        est.mle([10, 10, 10, 1], p, 50., **kw); print("2) no error!", kw)
    except ValueError as e:
        print("2) ValueError", kw, str(e)[:60])
# tiny tol (positive) - does it terminate quickly?
t = time.time(); est.mle([10, 10, 10, 1], p, 50., tol=1e-300); print("   tol=1e-300 time %.2fs" % (time.time() - t))
# 3) memory
pb = model(sbr=10.)
C = np.random.default_rng(1).multinomial(500, pb(np.array([3., 1.])), size=2000)
tracemalloc.start(); t = time.time()
r = est.mle(C, pb, search_radius=50., grid_step=0.2)
pk = tracemalloc.get_traced_memory()[1]; tracemalloc.stop()
print("3) M=2000 G~196k: peak %.0f MiB, %.1f s" % (pk / 2**20, time.time() - t))
tracemalloc.start(); t = time.time()
r = est.mle(C[:50], pb, search_radius=200., grid_step=0.2)
pk = tracemalloc.get_traced_memory()[1]; tracemalloc.stop()
print("   radius 200 step 0.2 (G~3.1M), M=50: peak %.0f MiB, %.1f s" % (pk / 2**20, time.time() - t))
# 4) c in (1,2)
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    for c in (1.1, 1.5, 2.0):
        pc = photons.make_model(patterns.tcp_centers(50.), beams.make_beam("donut", fwhm=300., power=c))
        a = cf.crb_tcp_center_limit(50., 100, 300., power=c); b = fisher.crb_limit(pc, 100)
        print("4) c=%.1f closed %.5f numeric %.5f rel %.2e" % (c, a, b, b / a - 1))
    print("   warnings:", [str(x.message)[:40] for x in w])
# unclear: far field NaN; negative counts
pg = model(kind="gaussian")
print("5) far field p:", pg(np.array([6000., 0])), "crb:", fisher.crb(pg, [6000., 0], 100))
print("   far field with sbr=10:", model(kind="gaussian", sbr=10.)(np.array([6000., 0])))
try:
    est.mle([-5, 10, 10, 0], p, 50.); print("   neg counts accepted!")
except ValueError as e:
    print("   neg counts ValueError")
# far field with 'limit' map around NaN
print("   crb_map gaussian far:", fisher.crb_map(pg, [5990., 6000.], [0.], 100))
