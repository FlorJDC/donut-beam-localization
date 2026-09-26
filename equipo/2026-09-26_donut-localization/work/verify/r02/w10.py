# -*- coding: utf-8 -*-
"""Verifier r02: own iterative MINFLUX loop. seed chosen by verifier (not 42 unless asked)."""
import numpy as np, sys
from scipy.special import xlogy
LN2 = np.log(2); FW = 300.

def pvec(pts, L, sbr):
    """pts (...,2) relative to TCP centre -> (...,4)"""
    a = 4*LN2/FW**2
    ang = np.pi/2 + 2*np.pi*np.arange(3)/3
    c = np.zeros((4, 2)); c[:3, 0] = L/2*np.cos(ang); c[:3, 1] = L/2*np.sin(ang)
    d2 = ((pts[..., None, :]-c)**2).sum(-1)
    I = np.e*a*d2*np.exp(-a*d2)
    q = I/I.sum(-1, keepdims=True)
    s = 1.0 if np.isinf(sbr) else sbr/(sbr+1)
    return s*q + (1-s)/4

def mle(n, L, sbr, R):
    """n (M,4) counts, TCP at origin; maximise over disk |r|<=R. returns (M,2)."""
    M = n.shape[0]
    h = R/40
    g = np.arange(-R, R+h/2, h)
    G = np.stack(np.meshgrid(g, g, indexing='ij'), -1).reshape(-1, 2)
    G = G[(G**2).sum(1) <= R*R*(1+1e-12)]
    lp = np.log(np.clip(pvec(G, L, sbr), 1e-300, None))
    ll = n.astype(float)@lp.T
    est = G[np.argmax(ll, 1)]
    for lev in range(4):
        loc = np.arange(-2*h, 2*h+1e-12, h/5)
        off = np.stack(np.meshgrid(loc, loc, indexing='ij'), -1).reshape(-1, 2)
        pts = est[:, None, :]+off[None]
        P = pvec(pts, L, sbr)
        ll = xlogy(n[:, None, :].astype(float), P).sum(-1)
        ll[(pts**2).sum(-1) > R*R] = -np.inf
        est = pts[np.arange(M), np.argmax(ll, 1)]
        h /= 5
    return est

def run(nrep, Ntot=1000, K=4, L0=150., Lf=25., sbr=np.inf, recenter=True, seed=1, rdisk=37.5, frac=0.75):
    rng = np.random.default_rng(seed)
    Ls = L0*(Lf/L0)**(np.arange(K)/(K-1))
    rr = rdisk*np.sqrt(rng.random(nrep)); th = 2*np.pi*rng.random(nrep)
    x = np.stack([rr*np.cos(th), rr*np.sin(th)], 1)
    c = np.zeros((nrep, 2)); Nk = Ntot//K
    errs = []
    for k, L in enumerate(Ls):
        p = pvec(x-c, L, sbr)
        n = np.array([rng.multinomial(Nk, pi) for pi in p])
        est = c + mle(n, L, sbr, frac*L)
        errs.append(est-x)
        if recenter: c = est
    return Ls, errs

def sig(e):
    s = np.sqrt((e**2).sum(1).mean()/2)
    # SE of RMS via delta method on mean of squared radial errors
    se = s*np.std((e**2).sum(1))/((e**2).sum(1).mean()*2*np.sqrt(len(e)))
    return s, se

if __name__ == '__main__':
    nrep = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    for sbr in (np.inf, 10.):
        for seed in (1, 2):
            Ls, errs = run(nrep, sbr=sbr, seed=seed)
            print('SBR', sbr, 'seed', seed, 'L', np.round(Ls, 2), 'sigma per iter', ['%.3f±%.3f' % sig(e) for e in errs],
                  ' bias final', np.round(errs[-1].mean(0), 3))
    Ls, errs = run(nrep, recenter=False, seed=3)
    print('no recenter, sigma per iter', ['%.2f±%.2f' % sig(e) for e in errs])
    Ls, errs = run(nrep, recenter=False, seed=4)
    print('no recenter, sigma per iter', ['%.2f±%.2f' % sig(e) for e in errs])
