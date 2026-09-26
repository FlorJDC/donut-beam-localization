# -*- coding: utf-8 -*-
"""Independent MC: multinomial sampling + MLE via vectorized coarse-to-fine grids. Verifier seed 20260926."""
import numpy as np, sys
from scipy.special import xlogy
from vfisher import probs, crb, LN2, E

L, F = 50.0, 300.0

def pgrid(pts, sbr):
    # pts (M,2) -> (M,4) vectorized p
    a = 4*LN2/F**2
    ang = np.pi/2 + 2*np.pi*np.arange(3)/3
    c = np.zeros((4, 2)); c[:3, 0] = L/2*np.cos(ang); c[:3, 1] = L/2*np.sin(ang)
    d2 = ((pts[:, None, :]-c[None])**2).sum(-1)
    I = E*a*d2*np.exp(-a*d2)
    q = I/I.sum(1, keepdims=True)
    s = 1.0 if np.isinf(sbr) else sbr/(sbr+1)
    return s*q+(1-s)/4

def mle(n, sbr, center, half, step0):
    """n (R,4). coarse grid around center (+-half) then 4 refinements x10."""
    R = n.shape[0]
    g = np.arange(-half, half+step0/2, step0)
    G = np.stack(np.meshgrid(g, g, indexing='ij'), -1).reshape(-1, 2)+center
    lp = np.log(np.clip(pgrid(G, sbr), 1e-300, None))
    est = np.empty((R, 2))
    for k in range(0, R, 200):
        nn = n[k:k+200].astype(float)
        ll = nn@lp.T
        # xlogy-like: where p=0 and n=0 contribution handled by clip (0*log(1e-300)=0)
        est[k:k+200] = G[np.argmax(ll, 1)]
    step = step0
    for lev in range(4):
        loc = np.arange(-1.5*step, 1.5*step+1e-12, step/10)
        off = np.stack(np.meshgrid(loc, loc, indexing='ij'), -1).reshape(-1, 2)
        for k in range(0, R, 200):
            pts = est[k:k+200, None, :]+off[None]
            P = pgrid(pts.reshape(-1, 2), sbr).reshape(pts.shape[0], -1, 4)
            ll = xlogy(n[k:k+200, None, :].astype(float), P).sum(-1)
            est[k:k+200] = pts[np.arange(pts.shape[0]), np.argmax(ll, 1)]
        step /= 10
    return est

def run(r0, N, sbr, R, seed, half=None):
    rng = np.random.default_rng(seed)
    p = pgrid(np.array([r0], float), sbr)[0]
    n = rng.multinomial(N, p, size=R)
    s27 = crb([0, 0], L, F, N)
    if half is None: half = 10*s27
    est = mle(n, sbr, np.array(r0, float), half, s27/8)
    return est, n

def summ(est, r0):
    d = est-np.array(r0)
    sig = np.sqrt(((d-d.mean(0))**2).sum(1).mean()/2)  # sqrt((var x + var y)/2)
    R = len(d)
    # SE of sigma via delta method on mean of per-rep squared deviations
    se = np.std(((d-d.mean(0))**2).sum(1)/2)/np.sqrt(R)/(2*sig)
    return sig, se, d.mean(0)

if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'V13':
        est, _ = run([0, 0], 100, 10, 8000, 20260926)
        sig, se, b = summ(est, [0, 0]); c = crb([0, 0], L, F, 100, sbr=10)
        print('V13 sigma', sig, '+-', se, 'crb', c, 'eff', sig/c, '+-', se/c, 'bias', b)
        print('edge hits', np.mean(np.abs(est).max(1) > 10*crb([0,0],L,F,100)-1))
    if which == 'V14':
        for N in [100, 1000, 10000]:
            est, n = run([0, 0], N, np.inf, 6000, 1000+N)
            sig, se, b = summ(est, [0, 0])
            lim = crb([1e-7, 0], L, F, N); s27 = crb([0, 0], L, F, N)
            print('V14 N', N, 'sigma', sig, '+-', se, 'sig/lim', sig/lim, '+-', se/lim, 'sig/S27', sig/s27, 'bias', b, 'n3>0 frac', np.mean(n[:, 3] > 0))
        est, n = run([2, 0], 100, np.inf, 8000, 777)
        sig, se, b = summ(est, [2, 0]); d = est[:, 0]-2
        print('V14 bias at (2,0): mean dx', d.mean(), '+-', d.std()/np.sqrt(len(d)), 'dy', b[1], 'sigma', sig, 'crb here', crb([2, 0], L, F, 100))
    if which == 'V15':
        est, _ = run([10, 0], 1000, np.inf, 6000, 4242)
        sig, se, b = summ(est, [10, 0]); c = crb([10, 0], L, F, 1000)
        print('V15 sigma', sig, '+-', se, 'crb', c, 'eff', sig/c, '+-', se/c, 'bias', b)
    if which == 'V16':
        rng = np.random.default_rng(99); N = 100; R = 40000
        for sbr in [np.inf, 10]:
            p = pgrid(np.zeros((1, 2)), sbr)[0]
            n = rng.multinomial(N, p, size=R)/N
            ang = np.pi/2+2*np.pi*np.arange(3)/3
            ri = L/2*np.stack([np.cos(ang), np.sin(ang)], 1)
            x = L**2*LN2/F**2; g = 1-x
            est0 = -(n[:, :3]@ri)/g
            sig = np.sqrt(est0.var(0).sum()/2)
            print('LMS sbr', sbr, 'sigma', sig, 'SE~', sig/np.sqrt(2*R), 'sig/S27', sig/crb([0,0],L,F,N), 'sig/lim', sig/crb([1e-7,0],L,F,N))
        # V17: contraction without 1/s, at true r off-center, large N, SBR=10
        for r0 in [[2, 0], [0, 3], [4, -2]]:
            for sbr in [10, np.inf]:
                p = pgrid(np.array([r0], float), sbr)[0]
                est_mean = -(p[:3]@ri)/g   # expectation of LMS = LMS on expected p
                print('V17 r0', r0, 'sbr', sbr, 'E[LMS] (no 1/s)', est_mean, 'ratio', est_mean/np.array(r0, float))
