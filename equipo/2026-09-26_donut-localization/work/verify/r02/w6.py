# -*- coding: utf-8 -*-
import numpy as np, sys
sys.path.insert(0, '..')
from vfisher import crb, centers
from rw import RW
n = 1.518; Nph = 100
def vec_crb(b, L, r, drop=True):
    C = centers(L)
    P = np.zeros((4, 3)); P[:, :2] = np.asarray(r)[None]-C
    E, Ex, Ey = b.grad(P)
    I = (np.abs(E)**2).sum(1)
    g = np.stack([2*np.real((np.conj(E)*Ex).sum(1)), 2*np.real((np.conj(E)*Ey).sum(1))], 1)
    S = I.sum(); gS = g.sum(0); p = I/S; gp = g/S - np.outer(I, gS)/S**2
    Fm = sum(np.outer(gp[i], gp[i])/p[i] for i in range(4) if p[i] > 0)*Nph
    return np.sqrt(np.trace(np.linalg.inv(Fm))/2)
def formula(b, L, Imax):
    # f, f' at d=L/2 ; c curvature
    P = np.array([[L/2, 0, 0], [0, 0, 0]])
    E, Ex, Ey = b.grad(P)
    f = (np.abs(E[0])**2).sum()/Imax; dfdr = 2*np.real((np.conj(E[0])*Ex[0]).sum())/Imax
    fp = dfdr/(2*(L/2)); c = (np.abs(Ex[1])**2).sum()/Imax
    al = (fp*L)**2/(2*f*f); be = 4*c/(3*f)
    return np.sqrt(0.5*(1/al+1/(al+be))/Nph)
for pol in ('right', 'left', 'x'):
    b = RW(n, pol=pol, Nt=300, Np=128)
    # Imax irrelevant for CRB; set 1
    out = []
    for L in (50., 100., 150.):
        v = [vec_crb(b, L, 1e-3*np.array([np.cos(t), np.sin(t)])) for t in (0.1, 1.3, 2.9)]
        s = '%s L=%g  CRB(|r|=1e-3, 3 dirs)=%s' % (pol, L, np.round(v, 5))
        if pol == 'right': s += '  formula=%.5f' % formula(b, L, 1.0)
        else: s += '  point r=0: %.5f' % vec_crb(b, L, [0., 0.])
        print(s)
for fw in (327.14, 300.):
    print('LG fwhm', fw, [round(crb([1e-7, 0.], L, fw, Nph), 5) for L in (50., 100., 150.)])
