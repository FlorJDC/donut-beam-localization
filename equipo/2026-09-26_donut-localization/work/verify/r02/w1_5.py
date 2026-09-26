# -*- coding: utf-8 -*-
import numpy as np, sys
from scipy.optimize import minimize_scalar, minimize
from rw import RW

def rad_prof(b, rr, z=0.):
    P = np.zeros((len(rr), 3)); P[:, 0] = rr; P[:, 2] = z
    return b.I(P)

def peak_radial(b):
    rr = np.linspace(0, 600, 1201); I = rad_prof(b, rr)
    i = np.argmax(I)
    res = minimize_scalar(lambda r: -rad_prof(b, [r])[0], bounds=(rr[max(i-2,0)], rr[i+2]), method='bounded', options={'xatol': 1e-6})
    return res.x, -res.fun

def peak_2d(b):
    best = (0, None)
    for ang in np.linspace(0, np.pi, 73):
        rr = np.linspace(0, 600, 601)
        P = np.zeros((len(rr), 3)); P[:, 0] = rr*np.cos(ang); P[:, 1] = rr*np.sin(ang)
        I = b.I(P); i = np.argmax(I)
        if I[i] > best[0]: best = (I[i], P[i, :2])
    res = minimize(lambda q: -b.I([[q[0], q[1], 0]])[0], best[1], method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-14})
    return -res.fun, res.x

for n in (1.518, 1.5):
    for Nt, Np in ((200, 96), (400, 160)):
        print('n=', n, 'grid', Nt, Np)
        bR = RW(n, pol='right', Nt=Nt, Np=Np)
        rp, Ip = peak_radial(bR)
        I0 = bR.I([[0, 0, 0], [0, 0, 300.], [0, 0, -500.]])
        print('  W1 right: I(0)/Imax at z=0,300,-500:', I0/Ip)
        # check rotational symmetry
        a = np.linspace(0, 2*np.pi, 7)[:-1]
        Pr = np.stack([rp*np.cos(a), rp*np.sin(a), 0*a], 1)
        print('  rot sym spread:', np.ptp(bR.I(Pr))/Ip)
        print('  W3 D_pp = %.4f nm' % (2*rp))
        E, Ex, Ey = bR.grad([[0, 0, 0]])
        c = (np.abs(Ex)**2).sum()/Ip; cy = (np.abs(Ey)**2).sum()/Ip
        print('  W4 c = %.5e (y: %.5e)  fwhm_curv = %.3f  fwhm_diam = %.3f' % (c, cy, np.sqrt(4*np.e*np.log(2)/c), 2*rp*np.sqrt(np.log(2))))
        # numeric check of c
        for rho in (0.1, 1.0):
            print('    I(rho)/Imax/rho^2 at', rho, rad_prof(bR, [rho])[0]/Ip/rho**2)
        # W5 z-flux fraction via Parseval: int |E_j|^2 dA ∝ sum |a_j|^2/(W * cos th) ... use direct weights
        # a_j contains w*sin th; energy density in (kx,ky): |a_j/(w sin th)|^2 / cos th * sin th dth dphi (up to const)
        bL = RW(n, pol='left', Nt=Nt, Np=Np)
        rpL, IpL = peak_2d(bL)[0], None
        ImL = rpL
        IL0 = bL.I([[0, 0, 0]])[0]
        print('  W2 opposite depth I0/Imax =', IL0/max(ImL, IL0), ' (Imax off-axis', ImL, 'I0', IL0, ')')
        bX = RW(n, pol='x', Nt=Nt, Np=Np)
        ImX, qX = peak_2d(bX); IX0 = bX.I([[0, 0, 0]])[0]
        print('  W2 linear depth =', IX0/ImX, ' peak at', qX)
