# -*- coding: utf-8 -*-
"""Verifier r03: package-behaviour claims (package = system under test; references are my own)."""
import os, sys, json
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H, '..'))
ROOT = os.path.abspath(os.path.join(H, *['..']*5)); sys.path.insert(0, os.path.join(ROOT, 'src'))
from vfisher import crb as mycrb
from mymle import sig_c, boot_se
from donutloc import fisher, photons, patterns, beams, montecarlo, experiments
p = photons.make_model(patterns.tcp_centers(50.), beams.make_beam('donut', fwhm=300.))
out = fisher.crb(p, [0, 0], np.array([100, 400]), zero_policy='limit')
print('crb limit array N ->', np.shape(out), out, ' my limits', [mycrb([1e-7, 0], 50., 300., N) for N in (100, 400)])
print('  crb_limit each N:', [fisher.crb_limit(p, N) for N in (100, 400)])
out2 = fisher.crb(p, np.array([[0, 0], [3, 1]]), np.array([100, 400]), zero_policy='limit')
print('  r (2,2) & N (2,):', out2, ' mine', mycrb([1e-7, 0], 50., 300., 100), mycrb([3, 1], 50., 300., 400))
try:
    print('  N (3,) vs r (2,2):', fisher.crb(p, np.array([[0, 0], [3, 1]]), np.array([100, 400, 900]), zero_policy='limit'))
except Exception as ex: print('  N(3,) vs r(2,2) raises', type(ex).__name__, ex)
print('  point policy array N:', fisher.crb(p, [0, 0], np.array([100, 400])))
# bootstrap
rng = np.random.default_rng(3)
for R in (500, 5000):
    e = rng.normal(size=(R, 2))*1.7
    a = montecarlo.bootstrap_sigma_se(e, n_boot=2000, seed=1)
    a = a[0] if isinstance(a, tuple) else a
    print('bootstrap gaussian R=%d: pkg %.5f  mine %.5f  sigma/(2sqrtR) %.5f  ratio pkg/gauss %.3f' % (R, a, boot_se(e, 2000, 7), sig_c(e)/(2*np.sqrt(R)), a/(sig_c(e)/(2*np.sqrt(R)))))
e = rng.standard_t(3, size=(5000, 2))
a = montecarlo.bootstrap_sigma_se(e, n_boot=2000, seed=1); a = a[0] if isinstance(a, tuple) else a
print('bootstrap t3 R=5000: pkg %.5f mine %.5f gauss %.5f' % (a, boot_se(e, 2000, 7), sig_c(e)/(2*np.sqrt(5000))))
# misalignment sigma_se: package 20x200, delta=10, centre
m = experiments.misalignment_study([10.0], n_patterns=20, n_rep=200, seed=5)
for k in ('honest', 'naive'):
    d = m[k]
    print('mis pkg %s: sigma_mean %s sigma_se %s sigma_se_within %s ratio %s' % (k, np.round(d['sigma_mean'], 3), np.round(d['sigma_se'], 4), np.round(d['sigma_se_within'], 4), np.round(np.asarray(d['sigma_se'])/np.asarray(d['sigma_se_within']), 2)))
