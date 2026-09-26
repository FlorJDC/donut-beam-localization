# -*- coding: utf-8 -*-
import numpy as np
from scipy.special import xlogy
LN2 = np.log(2); FW = 300.; L = 100.; N = 500; SBR = 10.
a = 4*LN2/FW**2
ang = np.pi/2 + 2*np.pi*np.arange(3)/3
C0 = np.zeros((4, 2)); C0[:3, 0] = L/2*np.cos(ang); C0[:3, 1] = L/2*np.sin(ang)
def pv(pts, C):
    d2 = ((pts[..., None, :]-C)**2).sum(-1)
    I = np.e*a*d2*np.exp(-a*d2); q = I/I.sum(-1, keepdims=True)
    s = SBR/(SBR+1); return s*q + (1-s)/4
def mle(n, C, R=30.):
    M = n.shape[0]; h = 1.0
    g = np.arange(-R, R+h/2, h); G = np.stack(np.meshgrid(g, g, indexing='ij'), -1).reshape(-1, 2)
    est = G[np.argmax(n.astype(float)@np.log(pv(G, C)).T, 1)]
    for lev in range(4):
        loc = np.arange(-2*h, 2*h+1e-12, h/5); off = np.stack(np.meshgrid(loc, loc, indexing='ij'), -1).reshape(-1, 2)
        pts = est[:, None, :]+off[None]
        ll = xlogy(n[:, None, :].astype(float), pv(pts, C)).sum(-1)
        est = pts[np.arange(M), np.argmax(ll, 1)]; h /= 5
    return est
if __name__ == "__main__":
  rng = np.random.default_rng(2026)
  for d in (0., 2., 5., 10.):
      bn, bh, sn, sh = [], [], [], []
      for k in range(20):
          phi = 2*np.pi*rng.random(4)
          C = C0 + d*np.stack([np.cos(phi), np.sin(phi)], 1)
          p = pv(np.zeros(2), C)
          n = rng.multinomial(N, p, size=200)
          en = mle(n, C0); eh = mle(n, C)
          bn.append(np.linalg.norm(en.mean(0))); bh.append(np.linalg.norm(eh.mean(0)))
          sn.append(np.sqrt(en.var(0).sum()/2)); sh.append(np.sqrt(eh.var(0).sum()/2))
      print('delta=%g naive bias %.2f  honest bias %.2f  | sigma naive %.2f honest %.2f' % (d, np.mean(bn), np.mean(bh), np.mean(sn), np.mean(sh)))
