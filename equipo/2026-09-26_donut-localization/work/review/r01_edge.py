# -*- coding: utf-8 -*-
"""Edge cases for estimators.mle / lms (round-1 review). Run from project root."""
import os, sys, subprocess
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from donutloc import beams, patterns, photons, estimators as est, fisher

p = photons.make_model(patterns.tcp_centers(50.), beams.make_beam(fwhm=300.))
code = ("import sys; sys.path.insert(0, r'%s');" % os.path.join(ROOT, "src") +
        "from donutloc import beams, patterns, photons, estimators as est;"
        "p = photons.make_model(patterns.tcp_centers(50.), beams.make_beam(fwhm=300.));"
        "print(est.mle([10,10,10,1], p, 50., tol=0.0))")
try:
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=20)
    print("tol=0:", out.stdout.strip(), out.stderr.strip()[-200:])
except subprocess.TimeoutExpired:
    print("tol=0: HANGS (killed after 20 s)")

# grid_step > search_radius: grid is only the centre
print("grid_step>R:", est.mle([30, 30, 30, 10], p, 5., grid_step=10.))
# empty batch
try:
    print("empty batch:", est.mle(np.zeros((0, 4)), p, 50.).shape)
except Exception as e:
    print("empty batch raises:", type(e).__name__, e)
# lms at r_lin where J^T J singular (flat model)
flat = lambda r: np.full(np.asarray(r).shape[:-1] + (4,), 0.25)
try:
    est.lms([1, 2, 3, 4], flat)
except Exception as e:
    print("lms singular raises:", type(e).__name__)
# counts of wrong K
try:
    print("mle K mismatch:", est.mle([1, 2, 3], p, 50.))
except Exception as e:
    print("mle K mismatch raises:", type(e).__name__, str(e)[:80])
# negative counts accepted silently?
print("mle negative counts:", est.mle([-5, 10, 10, 0], p, 50.))
# crb_map grid containing the exact origin
xs = np.linspace(-2, 2, 5)
print("crb_map row through origin:", fisher.crb_map(p, xs, [0.0], 100))
