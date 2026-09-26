# -*- coding: utf-8 -*-
"""Independent Fisher/CRB for TCP donut (verifier, round 1). Analytic gradients."""
import numpy as np
E = np.e
LN2 = np.log(2)

def centers(L):
    ang = np.pi/2 + 2*np.pi*np.arange(3)/3
    c = np.zeros((4, 2))
    c[:3, 0] = L/2*np.cos(ang); c[:3, 1] = L/2*np.sin(ang)
    return c

def intens(r, L, fwhm, eps=0.0, ped='const', quad=False):
    """returns I (4,), grad (4,2). LG donut peak-normalized to 1 (Eq S17)."""
    a = 4*LN2/fwhm**2
    c = centers(L)
    d = np.asarray(r, float)[None, :] - c
    d2 = (d**2).sum(1)
    if quad:   # quadratic limit: I = e a d2 (so fwhm cancels in p)
        I = E*a*d2; dId2 = E*a*np.ones(4)
        if eps: I = I + eps
    else:
        ex = np.exp(-a*d2)
        I = E*a*d2*ex; dId2 = E*a*ex*(1 - a*d2)
        if eps:
            if ped == 'const':
                I = I + eps
            else:
                I = I + eps*ex; dId2 = dId2 - a*eps*ex
    g = (dId2*2)[:, None]*d
    return I, g

def probs(r, L, fwhm, sbr=np.inf, eps=0.0, ped='const', c=1.0, quad=False):
    I, gI = intens(r, L, fwhm, eps, ped, quad)
    J = I**c
    gJ = (c*I**(c-1))[:, None]*gI if c != 1 else gI
    S = J.sum(); gS = gJ.sum(0)
    q = J/S; gq = gJ/S - np.outer(J, gS)/S**2
    s = 1.0 if np.isinf(sbr) else sbr/(sbr+1)
    p = s*q + (1-s)/4
    return p, s*gq

def fisher(r, L, fwhm, N=1, drop_zero=True, **kw):
    p, gp = probs(r, L, fwhm, **kw)
    F = np.zeros((2, 2))
    for i in range(4):
        if p[i] <= 0:
            if drop_zero: continue
        F += np.outer(gp[i], gp[i])/p[i]
    return N*F

def crb(r, L, fwhm, N=1, **kw):
    F = fisher(r, L, fwhm, N, **kw)
    return np.sqrt(np.trace(np.linalg.inv(F))/2)
