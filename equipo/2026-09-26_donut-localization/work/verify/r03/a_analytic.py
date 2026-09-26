# -*- coding: utf-8 -*-
"""Verifier r03: analytic sample of paper_numbers (own Fisher from ../vfisher.py, no donutloc)."""
import json, sys, os
import numpy as np
from scipy.optimize import minimize_scalar, brentq
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'r02'))
from vfisher import crb, fisher, probs, LN2
from cam import cam_crb
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), *['..']*5))
PN = json.load(open(os.path.join(ROOT, 'data', 'paper_numbers.json')))
out = {}
def chk(key, mine, tol=1e-6):
    v = PN[key]['value']; rel = abs(mine - v)/max(abs(v), 1e-300)
    ok = rel < tol
    out[key] = dict(mine=mine, theirs=v, rel=rel, ok=bool(ok))
    print('%-50s mine=%.10g theirs=%.10g rel=%.2e %s' % (key, mine, v, rel, 'OK' if ok else 'MISMATCH'))

F = 300.
def lim(L, N, fw=F, **kw):   # r->0 limit, average over directions (trace is direction-independent)
    vals = [crb(1e-6*np.array([np.cos(t), np.sin(t)]), L, fw, N, **kw) for t in (0.2, 1.4, 2.7)]
    return float(np.mean(vals))
# exponents
LF = np.array([5., 10., 20., 40.]); NF = np.array([100, 300, 1000, 3000, 10000.])
cL = [lim(L, 100) for L in LF]; cN = [lim(50., N) for N in NF]
chk('crb_exponent_L', np.polyfit(np.log(LF), np.log(cL), 1)[0], 1e-4)
chk('crb_exponent_N', np.polyfit(np.log(NF), np.log(cN), 1)[0], 1e-6)
# local exponents (robustness: range dependence)
print('  local L slopes', np.round(np.diff(np.log(cL))/np.diff(np.log(LF)), 5))
for rng_ in ([10, 50], [25, 150], [5, 40], [1, 10]):
    Ls = np.linspace(*rng_, 6); print('  exponent over', rng_, np.polyfit(np.log(Ls), np.log([lim(L, 100) for L in Ls]), 1)[0])
chk('crb_center_lg_L50_N100_nm', lim(50., 100), 1e-6)
chk('crb_center_point_S27_L50_N100_nm', crb([0., 0.], 50., F, 100))
chk('crb_center_sbr10_L50_N100_nm', crb([0., 0.], 50., F, 100, sbr=10))
chk('crb_center_sbr5_N500_fwhm360_L150_nm', crb([0., 0.], 150., 360., 500, sbr=5))
chk('crb_center_lg300_L150_N100_nm', lim(150., 100))
chk('crb_center_multiphoton_c2_L100_N100_quadratic_nm', crb([0., 0.], 100., F, 100, c=2.0, quad=True))
chk('limit_to_point_ratio_L150', lim(150., 100)/crb([0., 0.], 150., F, 100))
chk('limit_to_point_ratio_L5', lim(5., 100)/crb([0., 0.], 5., F, 100))
# limit ellipse axes L=50 (numeric Fisher incl. central term at r = 1e-6 x-hat)
Fm = fisher(np.array([1e-6, 0.]), 50., F, 100, drop_zero=False); C = np.linalg.inv(Fm)
chk('crb_limit_axis_par_L50_N100_nm', np.sqrt(C[0, 0]), 1e-5)
chk('crb_limit_axis_perp_L50_N100_nm', np.sqrt(C[1, 1]), 1e-5)
Fm = fisher(np.array([1e-6, 0.]), 100., F, 100, drop_zero=False, quad=True); C = np.linalg.inv(Fm)
chk('crb_limit_axis_ratio_quadratic', np.sqrt(C[0, 0]/C[1, 1]), 1e-5)
# eps degradation (V9 path, own)
chk('crb_eps_degradation_L100_eps0p05_gaussian', crb([0., 0.], 100., F, 100, eps=0.05, ped='gauss')/crb([0., 0.], 100., F, 100))
chk('crb_eps_degradation_L50_eps0p15_constant', crb([0., 0.], 50., F, 100, eps=0.15, ped='const')/crb([0., 0.], 50., F, 100))
# L_opt
for eps, tag in ((0.002, '0p002'), (0.05, '0p05'), (0.15, '0p15')):
    f = lambda L: crb([0., 0.], L, F, 100, eps=eps, ped='gauss')
    LL = np.geomspace(1, 450, 400); v = [f(L) for L in LL]; i = int(np.argmin(v))
    r = minimize_scalar(f, bounds=(LL[i-1], LL[i+1]), method='bounded', options={'xatol': 1e-9})
    chk('L_opt_eps%s_nm' % tag, r.x, 1e-5); chk('crb_opt_eps%s_nm' % tag, r.fun, 1e-7)
    chk('L_opt_over_fwhm_sqrt_eps_eps%s' % tag, r.x/(F*np.sqrt(eps)), 1e-5)
# divergence at ring diameter
chk('eps0_crb_divergence_L_nm', F/np.sqrt(LN2))
g = lambda L: 1/lim(L, 100)
print('  1/CRB at L=359, 360.3367, 361:', g(359.), g(360.3367), g(361.))
a = 4*LN2/F**2
chk('zero_depth_transition_scale_nm', np.sqrt(0.002/(np.e*a)))
# photons for 5 nm
chk('photons_for_5nm_L50_sbr20', (crb([0., 0.], 50., F, 1, sbr=20)/5)**2)
chk('photons_for_5nm_L50_limit_nobg', (lim(50., 1)/5)**2, 1e-5)
chk('iterative_crb_all_photons_L25_N1000_nm', lim(25., 1000))
chk('adaptive_iter0_crb_center_nm', lim(150., 250))
# camera
chk('camera_pixelated_9x9_N1000_nm', cam_crb(1000))
chk('camera_pixelated_9x9_N400_nm', cam_crb(400))
chk('camera_perpixel_sbr500_N600_nm', cam_crb(600, sbr=500/81.))
chk('camera_photons_for_5nm_perpixel_sbr500', (cam_crb(1, sbr=500/81.)/5)**2)
chk('camera_photons_for_5nm_total_sbr500', (cam_crb(1, sbr=500.)/5)**2)
chk('camera_sigma_nm', 100/np.sqrt(1000))
json.dump(out, open(os.path.join(os.path.dirname(__file__), 'a_analytic.json'), 'w'), indent=1)
print('n mismatches:', sum(not v['ok'] for v in out.values()), 'of', len(out))
