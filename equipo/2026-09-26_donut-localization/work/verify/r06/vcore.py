# -*- coding: utf-8 -*-
"""r06 verifier: own minimal model (no donutloc import).

LG donut I = e a d^2 exp(-a d^2) + eps, a = 4 ln2 / fwhm^2 (ring peak 1).
TCP: zeros at (L/2)(cos phi_k, sin phi_k), phi_k = pi/2 + 2 pi k/3, centre last.
Probabilities: fixed SBR (Balzarotti S30) or constant background per exposure (S28).
Fisher by complex-step derivatives; batched MLE by multilevel grid search in a disk.
"""
import numpy as np

E = np.e
LN2 = np.log(2.0)


def tcp(L, n_ring=3, centre=True, phi0=np.pi / 2):
    ang = phi0 + 2 * np.pi * np.arange(n_ring) / n_ring
    c = 0.5 * L * np.stack([np.cos(ang), np.sin(ang)], -1)
    if centre:
        c = np.vstack([c, [0.0, 0.0]])
    return c


def intens(r, cen, fwhm, eps=0.0):
    """r (...,2) possibly complex -> (...,K)."""
    a = 4 * LN2 / fwhm ** 2
    d = r[..., None, :] - cen
    d2 = d[..., 0] * d[..., 0] + d[..., 1] * d[..., 1]
    return E * a * d2 * np.exp(-a * d2) + eps


def probs(r, cen, fwhm, eps=0.0, sbr=None, bg=None):
    I = intens(r, cen, fwhm, eps)
    K = cen.shape[0]
    S = I.sum(-1)[..., None]
    if bg is not None:
        return (I + bg) / (S + K * bg)
    p0 = I / S
    if sbr is None or np.isinf(sbr):
        return p0
    s = sbr / (sbr + 1.0)
    return s * p0 + (1 - s) / K


def fisher(pf, r, N, h=1e-20, drop=()):
    r = np.asarray(r, float)
    p = np.real(pf(r.astype(complex)))
    g = []
    for j in range(2):
        rr = r.astype(complex)
        rr[j] += 1j * h
        g.append(np.imag(pf(rr)) / h)
    g = np.array(g)  # (2,K)
    F = np.zeros((2, 2))
    for i in range(p.size):
        if i in drop:
            continue
        F += np.outer(g[:, i], g[:, i]) / p[i]
    return N * F


def crb(pf, r, N, **kw):
    F = fisher(pf, r, N, **kw)
    return np.sqrt(0.5 * np.trace(np.linalg.inv(F)))


def _ll(counts, P):
    with np.errstate(divide="ignore"):
        lp = np.log(P)
    lp = np.where(np.isfinite(lp), lp, -1e300)
    return np.einsum("rk,rgk->rg", counts, lp) if P.ndim == 3 else counts @ lp.T


def mle_batch(counts, pf, R, h0=None, levels=7, shrink=4.0, n_loc=11, chunk=256):
    """Maximize sum n_i log p_i(r) over |r|<=R. Global grid (spacing h0) + local grids."""
    counts = np.atleast_2d(np.asarray(counts, float))
    if h0 is None:
        h0 = R / 60.0
    g = np.arange(-R, R + h0 / 2, h0)
    X, Y = np.meshgrid(g, g)
    pts = np.stack([X.ravel(), Y.ravel()], -1)
    pts = pts[np.hypot(pts[:, 0], pts[:, 1]) <= R]
    with np.errstate(divide="ignore"):
        LP = np.log(pf(pts))
    LP = np.where(np.isfinite(LP), LP, -1e300)
    out = np.empty((counts.shape[0], 2))
    for s0 in range(0, counts.shape[0], chunk):
        c = counts[s0:s0 + chunk]
        best = pts[np.argmax(c @ LP.T, axis=1)]
        hp = h0
        off = np.linspace(-1.25, 1.25, n_loc)
        for _ in range(levels):
            ox, oy = np.meshgrid(off * hp, off * hp)
            hp = 2.5 * hp / (n_loc - 1)
            o = np.stack([ox.ravel(), oy.ravel()], -1)
            cand = best[:, None, :] + o[None]
            P = pf(cand)
            ll = _ll(c, P)
            ll[np.hypot(cand[..., 0], cand[..., 1]) > R] = -np.inf
            best = cand[np.arange(c.shape[0]), np.argmax(ll, axis=1)]
        out[s0:s0 + chunk] = best
    return out


def lms_tcp(counts, L, fwhm, sbr=None, phi0=np.pi / 2):
    """Balzarotti S49-S50: r = -(1 - L^2 ln2/fwhm^2)^-1 (1/s) sum_periph p_i b_i."""
    counts = np.atleast_2d(np.asarray(counts, float))
    P = counts / counts.sum(1, keepdims=True)
    b = tcp(L, centre=False, phi0=phi0)
    g = 1.0 / (1.0 - L ** 2 * LN2 / fwhm ** 2)
    if sbr is not None and np.isfinite(sbr):
        g *= (sbr + 1.0) / sbr
    return -g * (P[:, :3] @ b)


def mlms_tcp(counts, L, fwhm, sbr=None, beta=(1.27, 3.8)):
    counts = np.atleast_2d(np.asarray(counts, float))
    p0 = counts[:, 3] / counts.sum(1)
    return (beta[0] + beta[1] * p0)[:, None] * lms_tcp(counts, L, fwhm, sbr)


def sigma(err):
    return np.sqrt(0.5 * (err[:, 0].var(ddof=1) + err[:, 1].var(ddof=1)))
