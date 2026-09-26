# -*- coding: utf-8 -*-
import numpy as np, sys
from scipy.special import erf
sys.path.insert(0, '..')
from vfisher import crb

def cam_crb(N, s=100., a=100., K=9, sbr=np.inf, x0=0., y0=0., norm=True):
    e = (np.arange(K+1)-K/2)*a
    def cdf(t, m): return 0.5*(1+erf((t-m)/(np.sqrt(2)*s)))
    def pdf(t, m): return np.exp(-(t-m)**2/(2*s*s))/(np.sqrt(2*np.pi)*s)
    qx = np.diff(cdf(e, x0)); qy = np.diff(cdf(e, y0))
    dqx = -np.diff(pdf(e, x0)); dqy = -np.diff(pdf(e, y0))
    q = np.outer(qx, qy); gx = np.outer(dqx, qy); gy = np.outer(qx, dqy)
    if norm:
        S = q.sum(); gxs = gx.sum(); gys = gy.sum()
        gx = gx/S - q*gxs/S**2; gy = gy/S - q*gys/S**2; q = q/S
    t = 1.0 if np.isinf(sbr) else sbr/(sbr+1)
    p = t*q + (1-t)/K**2
    gx *= t; gy *= t
    F = N*np.array([[(gx*gx/p).sum(), (gx*gy/p).sum()], [(gx*gy/p).sum(), (gy*gy/p).sum()]])
    return np.sqrt(np.trace(np.linalg.inv(F))/2)

print('W7 9x9 a=100 N=400:', cam_crb(400), ' unnorm:', cam_crb(400, norm=False))
print('   emitter at pixel corner:', cam_crb(400, x0=50, y0=50))
print('W7 a=10 K=121:', cam_crb(400, a=10., K=121)/5.0, cam_crb(400, a=10., K=121, norm=False)/5.0)
print('   a->1, K=1201:', cam_crb(400, a=1., K=1201)/5.0)
for sbr in (500,):
    print('W8 SBR=500 (total) N=600:', cam_crb(600, sbr=sbr))
    # per-pixel interpretation: signal peak pixel / bkg per pixel = 500?
    c1 = cam_crb(1, sbr=sbr); print('   N for 5 nm:', (c1/5)**2)
# alternative: b per pixel with SBR defined vs brightest pixel
e = (np.arange(10)-4.5)*100
from scipy.special import erf as E_
qx = np.diff(0.5*(1+E_(e/(np.sqrt(2)*100)))); q = np.outer(qx, qx); q /= q.sum()
print('   peak pixel fraction', q.max(), ' mean pixel frac', 1/81)
print('W9 S31 N needed for 5 nm, L=50 fwhm=300:')
for sbr in (np.inf, 50, 20, 10, 5):
    c1 = crb([0., 0.], 50., 300., 1, sbr=sbr)
    print('  SBR', sbr, (c1/5)**2)
print('  (limit r->0 at inf:', (crb([1e-7, 0.], 50., 300., 1)/5)**2, ')')
