# -*- coding: utf-8 -*-
"""Verificacion independiente (worker C, sesgo por fondo en TCP fijo).
No importa donutloc: modelo, Fisher (derivadas numericas) y MLE (grilla vectorizada + Nelder-Mead) propios.
"""
import numpy as np
from scipy.optimize import minimize

FWHM, L, N = 300.0, 100.0, 500
A = 4 * np.log(2) / FWHM ** 2
ang = np.pi / 2 + 2 * np.pi * np.arange(3) / 3
C = np.vstack([np.c_[L / 2 * np.cos(ang), L / 2 * np.sin(ang)], [0, 0]])


def lam(x, y):
    x = np.asarray(x, float)[..., None]; y = np.asarray(y, float)[..., None]
    r2 = (x - C[:, 0]) ** 2 + (y - C[:, 1]) ** 2
    return np.e * A * r2 * np.exp(-A * r2)


def prob(x, y, b):
    l = lam(x, y) + b
    return l / l.sum(-1, keepdims=True)


b5 = lam(0, 0).sum() / (4 * 5.0)
b20 = lam(0, 0).sum() / (4 * 20.0)
print("b(SBR5)=%.6g b(SBR20)=%.6g" % (b5, b20))
for b in (b5, b20):
    print("  SBR(x=50)=%.4f" % (lam(50, 0).sum() / (4 * b)))


def fisher(x, y, b, h=1e-4):
    th = np.array([x, y, b])
    p0 = prob(x, y, b)
    D = []
    for k in range(3):
        e = np.zeros(3); e[k] = h if k < 2 else h * 1e-3
        D.append((prob(*(th + e)) - prob(*(th - e))) / (2 * e[k]))
    D = np.array(D)
    return N * (D / p0) @ D.T


def crb(x, b):
    F = fisher(x, 0.0, b)
    ck = np.linalg.inv(F[:2, :2]); cf = np.linalg.inv(F)[:2, :2]
    return np.sqrt(np.trace(ck) / 2), np.sqrt(np.trace(cf) / 2), F[0, 2], F[1, 2]


for name, b in (("SBR5", b5), ("SBR20", b20)):
    for x in (0, 5, 10, 20, 50):
        k, f, fxb, fyb = crb(x, b)
        print("%s x=%2d CRBknown=%.4f CRBfree=%.4f Fxb=%.3g Fyb=%.3g" % (name, x, k, f, fxb, fyb))

# ---- MLE -------------------------------------------------------------------
R = 150.0
g = np.arange(-R, R + 1e-9, 2.0)
GX, GY = np.meshgrid(g, g, indexing="ij")
mask = GX ** 2 + GY ** 2 <= R ** 2
GX, GY = GX[mask], GY[mask]


def fit(counts, bfit, free=False):
    """counts (R,4). bfit fixed value (free=False) or start. returns (R,2) or (R,3)."""
    if not free:
        LP = np.log(np.clip(prob(GX, GY, bfit), 1e-300, None))  # (G,4)
        ll = counts @ LP.T
        i0 = ll.argmax(1)
        out = []
        for c, i in zip(counts, i0):
            f = lambda t: -(c * np.log(np.clip(prob(t[0], t[1], bfit), 1e-300, None))).sum()
            best = None
            for s in ([GX[i], GY[i]],):
                r = minimize(f, s, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-9, maxiter=2000))
                if best is None or r.fun < best.fun:
                    best = r
            out.append(best.x)
        return np.array(out)
    bgrid = np.concatenate([[0.0], np.geomspace(1e-4, 0.3, 25)])
    out = []
    LPs = [np.log(np.clip(prob(GX, GY, bb), 1e-300, None)) for bb in bgrid]
    lls = np.stack([counts @ LP.T for LP in LPs], 1)  # (R,B,G)
    for c, llc in zip(counts, lls):
        j, i = np.unravel_index(llc.argmax(), llc.shape)
        f = lambda t: np.inf if (t[2] < 0 or t[2] > 0.3) else -(c * np.log(np.clip(prob(t[0], t[1], t[2]), 1e-300, None))).sum()
        s = [GX[i], GY[i], max(bgrid[j], 1e-5)]
        r = minimize(f, s, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9, maxiter=4000))
        # también con b=0 exacto (frontera)
        f0 = lambda t: -(c * np.log(np.clip(prob(t[0], t[1], 0.0), 1e-300, None))).sum()
        r0 = minimize(f0, s[:2], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-9))
        out.append(r.x if r.fun <= r0.fun else np.r_[r0.x, 0.0])
    return np.array(out)


def pop_fit(x, btrue, bfit):
    """MLE sobre cuentas esperadas (sesgo poblacional); varios arranques."""
    c = N * prob(x, 0.0, btrue)
    f = lambda t: -(c * np.log(np.clip(prob(t[0], t[1], bfit), 1e-300, None))).sum()
    LP = np.log(np.clip(prob(GX, GY, bfit), 1e-300, None))
    ll = LP @ c
    order = np.argsort(-ll)[:200]
    res = []
    for i in order[::20]:
        r = minimize(f, [GX[i], GY[i]], method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-12, maxiter=4000))
        res.append((r.fun, r.x))
    res.sort(key=lambda z: z[0])
    return res


print("\n--- sesgo poblacional (modelo sin fondo) ---")
for name, b in (("SBR5", b5), ("SBR20", b20)):
    for x in (0, 5, 10, 50):
        res = pop_fit(x, b, 0.0)
        best = res[0]
        near = [z for z in res if z[0] - best[0] < 1e-6]
        print("%s x=%2d best=(%.3f,%.3f) nll=%.6f  n_equiv_minima=%d %s" % (
            name, x, best[1][0] - x, best[1][1], best[0], len(near),
            np.round([z[1] for z in near], 2).tolist()[:4]))
    for rel in (0.7, 1.3):
        m = 0
        for x in (0, 5, 10, 20, 30, 40, 50):
            bx = pop_fit(x, b, rel * b)[0][1]
            m = max(m, abs(bx[0] - x), abs(bx[1]))
        print("%s fondo x%.1f: sesgo pobl. max por eje = %.3f nm" % (name, rel, m))

print("\n--- Monte Carlo ---")
NREP = 1500
for name, b, xs, free in (("SBR5", b5, (0, 10, 20, 50), True), ("SBR20", b20, (0, 5), False)):
    for x in xs:
        rng = np.random.default_rng([42, 7, int(x)])
        cnt = rng.multinomial(N, prob(x, 0.0, b), size=NREP).astype(float)
        k, fcrb, _, _ = crb(x, b)
        for est, bfit in (("sin_fondo", 0.0), ("fondo_exacto", b)):
            e = fit(cnt, bfit)
            d = e - [x, 0]
            bias = d.mean(0); se = d.std(0) / np.sqrt(NREP)
            sig = np.sqrt((d.var(0).sum()) / 2); rmse = np.sqrt((d ** 2).sum(1).mean() / 2)
            print("%s x=%2d %-12s bias=(%.2f±%.2f, %.2f±%.2f) sigma=%.3f rmse=%.3f CRBk=%.3f sig/CRB=%.3f" % (
                name, x, est, bias[0], se[0], bias[1], se[1], sig, rmse, k, sig / k))
        if free and x <= 20:
            e = fit(cnt[:600], b, free=True)
            d = e[:, :2] - [x, 0]
            sig = np.sqrt(d.var(0).sum() / 2)
            print("%s x=%2d fondo_libre(600) sigma=%.3f CRBfree=%.3f frac_b0=%.3f bmean=%.4f" % (
                name, x, sig, fcrb, np.mean(e[:, 2] < 1e-6), e[:, 2].mean()))
