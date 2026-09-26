# -*- coding: utf-8 -*-
"""Verifier r03: own TCP model + grid/refine MLE (independent of donutloc). nm units."""
import numpy as np
from scipy.special import xlogy
LN2 = np.log(2); FW = 300.

def tcp(L):
    ang = np.pi/2 + 2*np.pi*np.arange(3)/3
    c = np.zeros((4, 2)); c[:3, 0] = L/2*np.cos(ang); c[:3, 1] = L/2*np.sin(ang)
    return c

def pvec(pts, C, sbr=np.inf, fw=FW):
    """pts (...,2) absolute; C (4,2) zero positions -> p (...,4). SBR against total beam signal."""
    a = 4*LN2/fw**2
    d2 = ((np.asarray(pts)[..., None, :]-C)**2).sum(-1)
    I = np.e*a*d2*np.exp(-a*d2)
    q = I/I.sum(-1, keepdims=True)
    s = 1.0 if np.isinf(sbr) else sbr/(sbr+1.)
    return s*q + (1-s)/4

def mle(n, C, sbr, R, center=(0., 0.), ncoarse=60, chunk=1000):
    """n (M,4) counts; maximise log-lik over disk |r-center|<=R. coarse grid (2R/ncoarse) + 5 refinements."""
    n = np.asarray(n, float); M = n.shape[0]; center = np.asarray(center, float)
    h0 = 2*R/ncoarse
    g = np.arange(-R, R+h0/2, h0)
    G = np.stack(np.meshgrid(g, g, indexing='ij'), -1).reshape(-1, 2)
    G = G[(G**2).sum(1) <= R*R*(1+1e-12)] + center
    lp = np.log(np.clip(pvec(G, C, sbr), 1e-300, None))
    est = np.empty((M, 2))
    for k in range(0, M, chunk):
        est[k:k+chunk] = G[np.argmax(n[k:k+chunk]@lp.T, 1)]
    h = h0
    loc = np.arange(-2, 2+1e-9, 0.2)
    off0 = np.stack(np.meshgrid(loc, loc, indexing='ij'), -1).reshape(-1, 2)
    for lev in range(5):
        off = off0*h
        for k in range(0, M, chunk):
            pts = est[k:k+chunk, None, :]+off[None]
            ll = xlogy(n[k:k+chunk, None, :], pvec(pts, C, sbr)).sum(-1)
            ll[((pts-center)**2).sum(-1) > R*R] = -np.inf
            est[k:k+chunk] = pts[np.arange(pts.shape[0]), np.argmax(ll, 1)]
        h /= 5
    return est

def multinom(rng, N, P):
    """row-wise multinomial for P (M,4)."""
    return rng.multinomial(N, P)

def sig_c(e):   # centred per-axis sigma sqrt((var_x+var_y)/2)
    return np.sqrt((e.var(0, ddof=1)).sum()/2)

def boot_se(e, nb=1000, seed=0, fn=sig_c):
    rng = np.random.default_rng(seed); R = len(e)
    return np.std([fn(e[rng.integers(0, R, R)]) for _ in range(nb)], ddof=1)
