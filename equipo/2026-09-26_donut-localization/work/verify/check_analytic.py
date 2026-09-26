# -*- coding: utf-8 -*-
import numpy as np
from vfisher import *

def xg(L, f): x = L**2*LN2/f**2; return x, 1-x
# gradient check vs central FD
for kw in [dict(), dict(sbr=5), dict(eps=0.05), dict(eps=0.05, ped='gauss'), dict(c=2.0), dict(sbr=3, c=3.0)]:
    r = np.array([7.3, -4.1]); h = 1e-5
    p, gp = probs(r, 80, 300, **kw)
    fd = np.array([(probs(r+h*e, 80, 300, **kw)[0]-probs(r-h*e, 80, 300, **kw)[0])/(2*h) for e in np.eye(2)]).T
    print('gradcheck', kw, np.abs(fd-gp).max()/np.abs(gp).max())

print('--- V1 point value')
for L, N, f in [(50, 100, 300), (100, 100, 300), (150, 1000, 360)]:
    x, g = xg(L, f)
    print(L, N, f, crb([0, 0], L, f, N), L/(2*np.sqrt(2*N))/g)

print('--- V2 limit r->0 (numeric at |r|=1e-6 several directions) vs formula')
for L, N, f in [(50, 100, 300), (5, 100, 300), (100, 100, 300), (150, 100, 300), (150, 100, 200)]:
    x, g = xg(L, f)
    form = np.sqrt(L**2/(8*N*g**2)*(3*g**2+np.exp(x))/(3*g**2+2*np.exp(x)))
    nums = [crb(1e-6*np.array([np.cos(t), np.sin(t)]), L, f, N) for t in [0, 0.3, 1.1, 2.5]]
    pt = crb([0, 0], L, f, N)
    print(L, f, 'formula', round(form, 6), 'numeric', np.round(nums, 6), 'rho', round(form/pt, 5))
print('quadratic rho', crb([1e-6, 0], 50, 300, 100, quad=True)/crb([0, 0], 50, 300, 100, quad=True), 2/np.sqrt(5))
print('quadratic sigma^2*10N/L^2', crb([1e-6, 0], 50, 300, 100, quad=True)**2*1000/2500)

print('--- V4 ellipse L=100 N=100 quadratic')
F = fisher([1e-7, 0], 100, 300, 100, quad=True); ev = np.linalg.eigvalsh(F)
print('sigma par, perp', 1/np.sqrt(ev[::-1]))
F = fisher(1e-7*np.array([np.cos(1), np.sin(1)]), 100, 300, 100, quad=True); w, v = np.linalg.eigh(F)
print('eigvec of largest (should be r-hat=(0.540,0.841))', v[:, 1])
x, g = xg(50, 300); N = 100; L = 50
al = 8*N*g**2/L**2; be = 16*N*np.exp(x)/(3*L**2)
F = fisher([1e-7, 0], L, 300, N); print('donut L50 F', F.diagonal(), 'alpha+beta, alpha', al+be, al)
# center term alone at small r
for rr in [1e-2, 1e-4, 1e-6]:
    p, gp = probs([rr, 0], L, 300); print('center term', rr, gp[3, 0]**2/p[3]*N, be)

print('--- V6 finite SBR: point vs limit vs S31')
for sbr in [5, 10]:
    for L in [50, 100]:
        x, g = xg(L, 300)
        s31 = L/(2*np.sqrt(2*N))/g*np.sqrt((1+1/sbr)*(1+3/(4*sbr)))
        print(sbr, L, crb([0, 0], L, 300, N, sbr=sbr), crb([1e-6, 0], L, 300, N, sbr=sbr), s31)
# transition radius: CRB vs r at large SBR
L = 50; x, g = xg(L, 300)
for sbr in [1e4, 1e6]:
    rc = np.sqrt(3)*L/4*np.exp(-x/2)/np.sqrt(sbr)
    p, gp = probs([rc, 0], L, 300, sbr=sbr); s = sbr/(sbr+1)
    # at r_c the center 'signal' part of p3 should equal background part
    q3 = (p[3]-(1-s)/4)/s
    print('rc', sbr, rc, 'signal part', s*q3, 'bg part', (1-s)/4,
          'center term / beta', gp[3, 0]**2/p[3]*N/be, '(expect ~1/2)')

print('--- V7 constant pedestal = S31 with SBR_eps')
for L in [50, 100]:
    x, g = xg(L, 300)
    for eps in [0.002, 0.01, 0.05, 0.15]:
        se = 3*E*x*np.exp(-x)/(4*eps)
        a = crb([0, 0], L, 300, N, eps=eps); b = crb([0, 0], L, 300, N, sbr=se)
        lim = crb([1e-6, 0], L, 300, N, eps=eps)
        print(L, eps, a, lim, b)
# combination
L = 50; x, g = xg(L, 300); eps = 0.01; se = 3*E*x*np.exp(-x)/(4*eps)
S = 3*E*x*np.exp(-x)
for sbr in [5, 20]:
    # (a) SBR vs pure LG signal: background per exposure b = S/(4 sbr) ; emulate via eps_tot
    b = S/(4*sbr); a1 = crb([0, 0], L, 300, N, eps=eps+b)
    eff = 1/(1/sbr+1/se); print('combo LG', a1, crb([0, 0], L, 300, N, sbr=eff))
    b = (S+4*eps)/(4*sbr); a2 = crb([0, 0], L, 300, N, eps=eps+b)
    eff = 1/((1+1/sbr)*(1+1/se)-1); print('combo total', a2, crb([0, 0], L, 300, N, sbr=eff))

print('--- V8 gaussian pedestal ratio')
for L in [50, 100]:
    x, g = xg(L, 300)
    for eps in [0.01, 0.15, 0.5]:
        F = fisher([0, 0], L, 300, N, eps=eps, ped='gauss')[0, 0]; F0 = 8*N*g**2/L**2
        form = 3*x**2*(E*g-eps)**2/(g**2*(E*x+eps)*(3*(E*x+eps)+eps*np.exp(x)))
        Fl = fisher([1e-6, 0], L, 300, N, eps=eps, ped='gauss')
        # equivalent background: S31 with SBR matching p3 share
        print(L, eps, F/F0, form, 'limit diag', Fl.diagonal()/F0)
# not equivalent to background: compare with S30 at the same center-share
L = 50; eps = 0.05; p, _ = probs([0, 0], L, 300, eps=eps, ped='gauss')
s = 1-4*p[3]; sbr = s/(1-s)
print('gauss ped CRB', crb([0, 0], L, 300, N, eps=eps, ped='gauss'), 'bg with same p', crb([0, 0], L, 300, N, sbr=sbr))

print('--- V9 degradation ratios')
for L in [50, 100]:
    x, g = xg(L, 300); s27 = L/(2*np.sqrt(2*N))/g
    print(L, [(eps, round(crb([0, 0], L, 300, N, eps=eps)/s27, 3), round(crb([0, 0], L, 300, N, eps=eps, ped='gauss')/s27, 3)) for eps in [0.002, 0.01, 0.05, 0.15]])

print('--- V10 eps=0.002 L=100: CRB vs r')
for r in [0, 1, 2, 3, 4, 4.9, 5, 6, 8, 10, 15]:
    print(r, round(crb([r, 0], 100, 300, N, eps=0.002), 4), round(crb([0, r], 100, 300, N, eps=0.002), 4), round(crb([r/np.sqrt(2), r/np.sqrt(2)], 100, 300, N, eps=0.002), 4))
rr = np.linspace(0, 20, 2001); cc = [crb([q, 0], 100, 300, N, eps=0.002) for q in rr]; print('min along x', rr[np.argmin(cc)], min(cc))
for L in [50, 100, 150]:
    cc = [crb([q, 0], L, 300, N, eps=0.002) for q in rr]; print('L', L, 'argmin r', rr[np.argmin(cc)], 'r_t', np.sqrt(0.002/(E*4*LN2/300**2)))

print('--- V11 multiphoton')
for cexp in [1, 2, 3]:
    print(cexp, crb([0, 0], 100, 300, 100, c=cexp, quad=True), crb([1e-4, 0], 100, 300, 100, c=cexp, quad=True), crb([1e-6, 0], 100, 300, 100, c=cexp, quad=True))
x, g = xg(80, 300)
for cexp in [1.5, 2, 3]:
    print('donut L80', cexp, crb([0, 0], 80, 300, 100, c=cexp), 80/(2*cexp*np.sqrt(200))/g, crb([1e-6, 0], 80, 300, 100, c=cexp))

print('--- V12 Masullo')
for f in [360, 300]:
    print(f, [round(crb([0, 0], L, f, 500, sbr=5), 3) for L in [50, 100, 150]])
print('1D', 50/(4*10))
# limiting cases
print('fwhm->inf', crb([0,0],50,1e7,100), 50/(2*np.sqrt(200)))
