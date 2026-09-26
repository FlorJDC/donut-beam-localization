# -*- coding: utf-8 -*-
"""Verifier r03: vectorial numbers with own direct-2D RW quadrature (../r02/rw.py); then test the
package's new linear harmonic table against my exact CRB (package called only as system under test)."""
import json, os, sys, time
import numpy as np
from scipy.optimize import minimize_scalar
H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(H, '..')); sys.path.insert(0, os.path.join(H, '..', 'r02'))
from rw import RW
from vfisher import centers
ROOT = os.path.abspath(os.path.join(H, *['..']*5))
PN = json.load(open(os.path.join(ROOT, 'data', 'paper_numbers.json')))
LN2 = np.log(2)
def rep(k, mine, tol):
    v = PN[k]['value']; rel = abs(mine-v)/max(abs(v), 1e-300) if v != 0 else abs(mine)
    print('%-44s mine=%.8g theirs=%.8g rel=%.2e %s' % (k, mine, v, rel, 'OK' if rel < tol else 'MISMATCH'))

def vcrb(b, L, r, N=100):
    C = centers(L); P = np.zeros((4, 3)); P[:, :2] = np.asarray(r, float)[None]-C
    E, Ex, Ey = b.grad(P); I = (np.abs(E)**2).sum(1)
    g = np.stack([2*np.real((np.conj(E)*Ex).sum(1)), 2*np.real((np.conj(E)*Ey).sum(1))], 1)
    S = I.sum(); gS = g.sum(0); p = I/S; gp = g/S-np.outer(I, gS)/S**2
    Fm = N*sum(np.outer(gp[i], gp[i])/p[i] for i in range(4) if p[i] > 0)
    return np.sqrt(np.trace(np.linalg.inv(Fm))/2)
def vlim(b, L, N=100):
    return np.mean([vcrb(b, L, 1e-3*np.array([np.cos(t), np.sin(t)]), N) for t in (0.1, 0.9, 1.7, 2.5, 3.3, 4.1)])

bc = RW(pol='right', Nt=300, Np=64)
# ring peak & diameter (circular symmetric -> radial along x)
f = lambda r: -bc.I(np.array([[r, 0, 0.]]))[0]
res = minimize_scalar(f, bounds=(120, 260), method='bounded', options={'xatol': 1e-7})
Imax = -res.fun; D = 2*res.x
rep('vectorial_peak_to_peak_diameter_nm', D, 1e-5)
E, Ex, Ey = bc.grad(np.zeros((1, 3))); c = (np.abs(Ex[0])**2).sum()/Imax
rep('vectorial_zero_curvature_nm2', c, 1e-5)
rep('vectorial_lg_fwhm_curvature_nm', np.sqrt(4*np.e*LN2/c), 1e-5)
rep('vectorial_lg_fwhm_diameter_nm', D*np.sqrt(LN2), 1e-5)
print('  zero depth correct I(0)/Imax =', bc.I(np.zeros((1, 3)))[0]/Imax)
bl = RW(pol='x', Nt=300, Np=64)
g = np.linspace(-320, 320, 161); X, Y = np.meshgrid(g, g)
pts = np.stack([X.ravel(), Y.ravel(), 0*X.ravel()], 1)
I = np.concatenate([bl.I(pts[i:i+2000]) for i in range(0, len(pts), 2000)])
j = np.argmax(I); from scipy.optimize import minimize
r2 = minimize(lambda q: -bl.I(np.array([[q[0], q[1], 0.]]))[0], pts[j, :2], method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-14})
rep('vectorial_zero_depth_linear', bl.I(np.zeros((1, 3)))[0]/(-r2.fun), 1e-5)
bo = RW(pol='left', Nt=300, Np=64)
Io = -minimize_scalar(lambda r: -bo.I(np.array([[r, 0, 0.]]))[0], bounds=(0, 300), method='bounded').fun
Iog = bo.I(pts[::7]).max()
print('  opposite: I(0)/max(line search)=%.7f, grid max check %.4f' % (bo.I(np.zeros((1, 3)))[0]/Io, bo.I(np.zeros((1, 3)))[0]/max(Io, Iog)))
rep('vectorial_zero_depth_wrong_handedness', bo.I(np.zeros((1, 3)))[0]/max(Io, Iog), 1e-4)
# CRBs
for L in (50., 100., 150.):
    rep('crb_center_vec_linear_L%d_N100_nm' % L, vlim(bl, L), 2e-4)
rep('crb_center_vec_correct_L150_N100_nm', vlim(bc, 150.), 2e-4)
rep('crb_center_vec_opposite_L100_N100_nm', vcrb(bo, 100., [0, 0]), 2e-4)
from vfisher import crb as lgcrb
rep('crb_vec_over_lgcurv_L150_N100', vlim(bc, 150.)/lgcrb([1e-7, 0], 150., PN['vectorial_lg_fwhm_curvature_nm']['value'], 100), 2e-4)
ex73 = vcrb(bl, 50., [7., 3.]); exlim = vlim(bl, 50.)
print('MY exact linear L=50,N=100: CRB(7,3)=%.6f  limit=%.6f' % (ex73, exlim))
# ---- package under test (fix claim) ----
sys.path.insert(0, os.path.join(ROOT, 'src'))
import warnings
from donutloc import vectorial as V, photons, patterns, fisher
for rho_max, nth in ((120.0, 401), (1500.0, 801)):
    t = time.time(); bh = V.make_vectorial_beam(rho_max=rho_max, polarization='linear', n_theta=nth)
    pm = photons.make_model(patterns.tcp_centers(50.0), bh)
    a = float(fisher.crb(pm, np.array([7., 3.]), 100)); la = float(fisher.crb_limit(pm, 100))
    print('PKG harmonic rho_max=%g n_theta=%d: CRB(7,3)=%.6f (rel %.2e vs mine)  limit=%.6f (rel %.2e)  [%.1fs] method=%s' %
          (rho_max, nth, a, a/ex73-1, la, la/exlim-1, time.time()-t, bh.linear_method))
    # direct intensity comparison vs my RW at random points (normalise by my Imax)
    rng = np.random.default_rng(5); q = rng.uniform(-80, 80, (30, 2))
    mine = bl.I(np.c_[q, np.zeros(30)])/(-r2.fun); pk = bh(q[:, 0], q[:, 1])
    print('   intensity max rel diff vs my RW:', np.max(np.abs(pk/mine-1)))
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    bcart = V.make_vectorial_beam(rho_max=120.0, grid_step=5.0, polarization='linear', n_theta=401, linear_method='cartesian')
    print('cartesian warns:', [str(x.category.__name__) for x in w])
pm = photons.make_model(patterns.tcp_centers(50.0), bcart)
print('PKG cartesian CRB(7,3)=%.4f (%+.2f%%)  limit=%.4f (%+.2f%%)' % (float(fisher.crb(pm, np.array([7., 3.]), 100)), 100*(float(fisher.crb(pm, np.array([7., 3.]), 100))/ex73-1), float(fisher.crb_limit(pm, 100)), 100*(float(fisher.crb_limit(pm, 100))/exlim-1)))
