# -*- coding: utf-8 -*-
import os, sys, time
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), *[".."] * 5))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from donutloc import vectorial as v, photons, patterns, fisher, beams

# n_theta convergence of zero depths / curvature
for nt in (201, 401, 801, 3201):
    t = time.time()
    zo = v.zero_depth(handedness=-1, n_theta=nt)
    zl = v.zero_depth(polarization="linear", n_theta=nt)
    c = v.zero_curvature(n_theta=nt)
    D = v.peak_to_peak_diameter(n_theta=nt)
    print("n_theta %5d: opp %.7f lin %.7f curv %.6e Dpp %.4f (%.1fs)" % (nt, zo, zl, c, D, time.time() - t))
print("note-B params n=1.5: opp %.5f lin %.5f" % (v.zero_depth(handedness=-1, n=1.5), v.zero_depth(polarization="linear", n=1.5)))
print("correct hand zero_depth:", v.zero_depth(), " at z=300:", v.zero_depth(z=300.))
print("charge=-1,s=-1 (ls=+1):", v.zero_depth(charge=-1, handedness=-1), " charge=-1,s=+1:", v.zero_depth(charge=-1, handedness=1))
print("charge=0 (no vortex) zero_depth (expect 1):", v.zero_depth(charge=0))
print("linear pol_angle 0 vs 0.7:", v.zero_depth(polarization="linear"), v.zero_depth(polarization="linear", pol_angle=0.7))
# independent 2D check at generic points, several cases
rho = np.array([0., 37., 150., 260.]); phi = np.array([0., 0.4, 1.9, -2.5])
for opt in ({}, {"handedness": -1}, {"polarization": "linear", "pol_angle": 0.3}, {"charge": 2, "handedness": 1}, {"z": 200.}):
    a = v.focal_field(rho, phi, **opt); b = v.focal_field(rho, phi, method="2d", n_phi=128, **opt)
    Ia = sum(abs(c) ** 2 for c in a); Ib = sum(abs(c) ** 2 for c in b)
    err = max(np.max(abs(x - y)) for x, y in zip(a, b)) / np.sqrt(Ia.max())
    print("2d vs bessel", opt, "max comp err rel %.1e" % err)
# linear superposition: does the 2D linear field match the bessel superposition component-wise? (above)
# beam interpolation near zero (circular)
b = v.make_vectorial_beam()
be = v.make_vectorial_beam(mode="exact")
r = np.array([0.1, 0.5, 1.0, 1.5, 3.3, 25.0, 50.3, 191.0])
print("interp/exact ratio:", b(r, 0 * r) / be(r, 0 * r))
cc = v.zero_curvature()
print("interp I/(c r^2) small r:", b(r[:4], 0 * r[:4]) / (cc * r[:4] ** 2))
# CRB with interp vs exact beam
for L in (50., 100.):
    ce = fisher.crb_limit(photons.make_model(patterns.tcp_centers(L), be), 100)
    ci = fisher.crb_limit(photons.make_model(patterns.tcp_centers(L), b), 100)
    cg = float(fisher.crb(photons.make_model(patterns.tcp_centers(L), b), [7., 3.], 100))
    cge = float(fisher.crb(photons.make_model(patterns.tcp_centers(L), be), [7., 3.], 100))
    print("L=%g crb_limit exact %.6f interp %.6f (rel %.1e) ; at (7,3) exact %.6f interp %.6f (rel %.1e)" % (L, ce, ci, ci/ce-1, cge, cg, cg/cge-1))
# linear beam grid interpolation near zero: I ~ quadratic but bilinear interpolation on 5 nm grid
bl = v.make_vectorial_beam(polarization="linear", rho_max=400.)
ble = v.make_vectorial_beam(polarization="linear", mode="exact")
pts = np.array([[1., 0.5], [2.5, 2.5], [0.3, -0.2], [12., 7.], [31., -3.]])
print("linear interp/exact:", bl(pts[:, 0], pts[:, 1]) / ble(pts[:, 0], pts[:, 1]))
pl = photons.make_model(patterns.tcp_centers(50.), bl); ple = photons.make_model(patterns.tcp_centers(50.), ble)
print("linear crb at (7,3): interp %.5f exact %.5f; crb_limit interp %.5f exact %.5f" % (
    fisher.crb(pl, [7., 3.], 100), fisher.crb(ple, [7., 3.], 100), fisher.crb_limit(pl, 100), fisher.crb_limit(ple, 100)))
# timing of exact beam in an MLE-like call
t = time.time(); g = np.random.default_rng(0).uniform(-40, 40, (2000, 2)); pe = photons.make_model(patterns.tcp_centers(50.), be); pe(g); print("exact model on 2000 pts: %.2fs" % (time.time() - t))
t = time.time(); pi_ = photons.make_model(patterns.tcp_centers(50.), b); pi_(g); print("interp model on 2000 pts: %.4fs" % (time.time() - t))
