# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from rw import RW
def rad(b, rr):
    P = np.zeros((len(rr), 3)); P[:, 0] = rr; return b.I(P)
def peak_radial(b):
    rr = np.linspace(0, 500, 251); I = rad(b, rr); i = np.argmax(I)
    r = minimize_scalar(lambda x: -rad(b, [x])[0], bounds=(rr[max(i-2, 0)], rr[i+2]), method='bounded', options={'xatol': 1e-6})
    return r.x, -r.fun
def peak_2d(b):
    best = (0, None)
    for ang in np.linspace(0, np.pi, 37):
        rr = np.linspace(0, 500, 126); P = np.zeros((len(rr), 3)); P[:, 0] = rr*np.cos(ang); P[:, 1] = rr*np.sin(ang)
        I = b.I(P); i = np.argmax(I)
        if I[i] > best[0]: best = (I[i], P[i, :2])
    r = minimize(lambda q: -b.I([[q[0], q[1], 0]])[0], best[1], method='Nelder-Mead', options={'xatol': 1e-7, 'fatol': 1e-15})
    return -r.fun, r.x
n = 1.5
for pol in ('left', 'x', 'y'):
    b = RW(n, pol=pol, Nt=200, Np=96); Im, q = peak_2d(b); I0 = b.I([[0, 0, 0]])[0]
    print(pol, 'depth n=1.5:', I0/max(Im, I0), q)
for F in (0.5, 1., 2., 100., 1e4):
    b = RW(n, F=F, Nt=200, Np=96); rp, _ = peak_radial(b); print('F', F, 'Dpp n=1.5 = %.3f' % (2*rp))
for F in (0.5, 1., 2., 100.):
    b = RW(1.518, F=F, Nt=200, Np=96); rp, _ = peak_radial(b); print('F', F, 'Dpp n=1.518 = %.3f' % (2*rp))
