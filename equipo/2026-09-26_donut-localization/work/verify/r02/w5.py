# -*- coding: utf-8 -*-
import numpy as np
from rw import RW
n = 1.518; NA = 1.4; tmax = np.arcsin(NA/n); F = 5/3.
bL = RW(n, pol='left', Nt=200, Np=96)
E0 = bL.field([[0, 0, 0]])[0]
print('left, on-axis field components (abs):', np.abs(E0))
x, w = np.polynomial.legendre.leggauss(400); th = 0.5*tmax*(x+1); wt = 0.5*tmax*w
g = np.exp(-(np.sin(th)/(F*np.sin(tmax)))**2)*np.sqrt(np.cos(th))*np.sin(th)
Ez_pred = -2*np.pi*(g*np.sin(th)*wt).sum()/np.sqrt(2)   # J0(0)=1, Er=e^{-i phi}/sqrt2 times e^{i phi}
print('Ez(0) numeric', E0[2], ' predicted -2pi/sqrt2 int g sin =', Ez_pred)
# z-dependence: E_z(0,z) = -2pi/sqrt2 int g sin e^{ikz cos}
k = 2*np.pi*n/640
for z in (200., 500.):
    print(' z', z, bL.field([[0, 0, z]])[0][2], -2*np.pi/np.sqrt(2)*(g*np.sin(th)*np.exp(1j*k*z*np.cos(th))*wt).sum())
# W5: z-energy fraction, right-handed, Parseval
for FF in (5/3., 1., 100.):
    b = RW(n, F=FF, pol='right', Nt=400, Np=128)
    A = b.a  # includes sin th * weights; energy weight: |A/(sin w)|^2 / cos * sin * w  = |A|^2/(cos sin w)
    Wt = np.abs(b.a[0])**2*0  # placeholder
    TH = b.TH
    # recover weights: sin*w = |kv|... recompute
    xx, ww = np.polynomial.legendre.leggauss(400); tt = 0.5*tmax*(xx+1); wtt = 0.5*tmax*ww
    sw = np.repeat(np.sin(tt)*wtt, 128)*(2*np.pi/128)
    e = (np.abs(A)**2)/(np.cos(TH)*sw)
    print('F', FF, 'Parseval z fraction:', e[2].sum()/e.sum())
# brute force in focal plane (radially symmetric components) for F=5/3
b = RW(n, pol='right', Nt=300, Np=128)
rr = np.concatenate([np.linspace(0, 3000, 3001), np.linspace(3001, 20000, 4000)])
P = np.zeros((len(rr), 3)); P[:, 0] = rr
Es = np.concatenate([b.field(P[i:i+500]) for i in range(0, len(rr), 500)])
Ij = np.abs(Es)**2
tot = np.trapz(Ij.sum(1)*rr, rr); zz = np.trapz(Ij[:, 2]*rr, rr)
print('focal-plane integral z fraction (to 20 um):', zz/tot, ' (to 3 um:', np.trapz(Ij[:3001, 2]*rr[:3001], rr[:3001])/np.trapz(Ij[:3001].sum(1)*rr[:3001], rr[:3001]), ')')
