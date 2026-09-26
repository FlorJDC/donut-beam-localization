# -*- coding: utf-8 -*-
"""Verifier r02: Richards-Wolf focal field by direct 2D pupil quadrature (plane-wave sum).
No Bessel reduction. theta: Gauss-Legendre; phi: uniform (spectrally exact for trig polys)."""
import numpy as np
LAM = 640.0; NA = 1.4

class RW:
    def __init__(self, n=1.518, F=5/3., pol='right', l=1, Nt=300, Np=128, NA=NA, lam=LAM):
        self.k = 2*np.pi*n/lam
        tmax = np.arcsin(NA/n)
        x, w = np.polynomial.legendre.leggauss(Nt)
        th = 0.5*tmax*(x+1); wt = 0.5*tmax*w
        ph = 2*np.pi*np.arange(Np)/Np; wp = 2*np.pi/Np
        TH, PH = np.meshgrid(th, ph, indexing='ij')
        W = (wt[:, None]*wp)*np.ones_like(TH)
        # pupil radius rho' = f sin(theta), aperture h = f sin(tmax); amplitude Gaussian exp(-(rho'/w0)^2), F=w0/h
        amp = np.exp(-(np.sin(TH)/(F*np.sin(tmax)))**2) * np.sqrt(np.cos(TH)) * np.exp(1j*l*PH)
        if pol == 'right': Ei = np.array([1, 1j])/np.sqrt(2)      # sigma=+1
        elif pol == 'left': Ei = np.array([1, -1j])/np.sqrt(2)    # sigma=-1
        elif pol == 'x': Ei = np.array([1, 0j])
        elif pol == 'y': Ei = np.array([0j, 1])
        c, s = np.cos(PH), np.sin(PH)
        Er = Ei[0]*c + Ei[1]*s            # E . rho_hat
        Ep = -Ei[0]*s + Ei[1]*c           # E . phi_hat
        ct, st = np.cos(TH), np.sin(TH)
        a = np.empty((3,)+TH.shape, complex)
        a[0] = Ep*(-s) + Er*ct*c
        a[1] = Ep*c + Er*ct*s
        a[2] = Er*(-st)
        self.a = (a*amp*st*W).reshape(3, -1)       # includes sin(theta) dtheta dphi
        self.kv = self.k*np.stack([st*c, st*s, ct]).reshape(3, -1)
        self.TH = TH.ravel(); self.amp2 = (np.abs(amp)**2).ravel()
    def field(self, pts):
        pts = np.atleast_2d(pts)            # (M,3)
        ph = np.exp(1j*pts@self.kv)         # (M,K)
        return ph@self.a.T                   # (M,3)
    def grad(self, pts):
        """returns E (M,3), dE/dx (M,3), dE/dy (M,3)"""
        pts = np.atleast_2d(pts); ph = np.exp(1j*pts@self.kv)
        E = ph@self.a.T
        Ex = ph@(self.a*1j*self.kv[0]).T; Ey = ph@(self.a*1j*self.kv[1]).T
        return E, Ex, Ey
    def I(self, pts):
        E = self.field(pts); return (np.abs(E)**2).sum(1)
