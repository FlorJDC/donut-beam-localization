# -*- coding: utf-8 -*-
"""Verificacion independiente (worker B, p-MINFLUX timing). No importa donutloc."""
import numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(12345)
FWHM, L, T, K = 300.0, 100.0, 12.5, 4
a = 4*np.log(2)/FWHM**2
ang = np.pi/2 + 2*np.pi*np.arange(3)/3
C = np.vstack([0.5*L*np.c_[np.cos(ang), np.sin(ang)], [[0, 0]]])

def I(x, y):  # (..., 4)
    r2 = (x[..., None]-C[:, 0])**2 + (y[..., None]-C[:, 1])**2
    return np.e*a*r2*np.exp(-a*r2)
def p_of(x, y):
    i = I(np.asarray(x, float), np.asarray(y, float)); return i/i.sum(-1, keepdims=True)
def dp(x, y, h=1e-4):
    return np.stack([(p_of(x+h, y)-p_of(x-h, y))/(2*h), (p_of(x, y+h)-p_of(x, y-h))/(2*h)], -1)

# --- 1. cross-talk: forma cerrada vs MC de fotones
def Mmat(tau, order=(0, 1, 2, 3)):
    q = np.exp(-T/tau); slot = np.empty(4, int); slot[list(order)] = np.arange(4)
    M = np.empty((4, 4))
    for i in range(4):
        for j in range(4):
            M[i, j] = (1-q)*q**((slot[i]-slot[j]) % 4)/(1-q**4)
    return M
def M_mc(tau, n=2_000_000, order=(0, 1, 2, 3)):
    # pulso en slot s a tiempo s*T + k*P; ventana del foton = floor(((sT + t) mod P)/T)
    P = 4*T; M = np.zeros((4, 4))
    for j in range(4):
        s = order.index(j); t = rng.exponential(tau, n)
        w = np.floor(((s*T + t) % P)/T).astype(int)
        cnt = np.bincount(w, minlength=4)/n
        for sw in range(4): M[order[sw], j] = cnt[sw]
    return M
print("== 1. cross-talk tau=5")
M5 = Mmat(5.0); print(np.round(M5[:, 0], 4)); print("MC", np.round(M_mc(5.0)[:, 0], 4))
print("col sums", M5.sum(0), "row sums", M5.sum(1))

# --- 2. CRB
def crb(x, y, M=None):
    p = p_of(x, y); d = dp(x, y)
    if M is not None: p = p@M.T; d = np.einsum('ij,...jk->...ik', M, d)
    with np.errstate(divide='ignore', invalid='ignore'):
        w = np.where(p > 1e-14, 1/p, 0.0)
    F = np.einsum('...i,...ia,...ib->...ab', w, d, d)
    Fi = np.linalg.inv(F); return np.sqrt((Fi[..., 0, 0]+Fi[..., 1, 1])/2)
g = np.arange(-50, 50.001, 2.5); X, Y = np.meshgrid(g, g); m = X**2+Y**2 <= 50**2+1e-9
xd, yd = X[m], Y[m]; print("disk pts", xd.size)
c0 = crb(xd, yd); c0c = crb(np.array(0.), np.array(0.))
# limite direccional en el centro (promedio sobre direcciones de r pequeno)
print("CRB ideal disco medio %.4f, centro punto %.4f, centro r=1e-3 dirs:" % (c0.mean(), c0c),
      [round(float(crb(np.array(1e-3*np.cos(t)), np.array(1e-3*np.sin(t)))), 3) for t in (0, .5, 1.)])
print("== 2. CRB ratios")
for order in [(0, 1, 2, 3), (0, 2, 1, 3), (3, 0, 1, 2), (1, 0, 3, 2)]:
    rr = []
    for tau in (1, 2, 3, 4, 5):
        M = Mmat(tau, order); cm = crb(xd, yd, M)
        rr.append((round(cm.mean()/c0.mean(), 4), round(np.mean(cm/c0), 4),
                   round(float(crb(np.array(0.), np.array(0.), M))/c0c, 4)))
    print(order, rr)

# --- MLE propio: grilla en disco R=75 paso 1 + Nelder-Mead
R = 75.0; gg = np.arange(-R, R+1e-9, 1.0); GX, GY = np.meshgrid(gg, gg); gm = GX**2+GY**2 <= R*R
GX, GY = GX[gm], GY[gm]
def mle(f, model):
    # f: (n,4) conteos o fracciones; model(x,y)->(...,4)
    LP = np.log(np.clip(model(GX, GY), 1e-300, None))
    ll = f@LP.T; k = ll.argmax(1); out = np.empty((len(f), 2))
    for n in range(len(f)):
        x0 = np.array([GX[k[n]], GY[k[n]]])
        def nll(v):
            if v@v > R*R: return 1e300
            return -(f[n]*np.log(np.clip(model(np.array(v[0]), np.array(v[1])), 1e-300, None))).sum()
        r = minimize(nll, x0, method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-12, maxiter=2000))
        out[n] = r.x
    return out
print("== 3. sesgo MLE ingenuo")
xs = np.linspace(0, 50, 51)
for tau, order in [(3, (0, 1, 2, 3)), (3, (0, 2, 1, 3)), (5, (0, 1, 2, 3)), (5, (0, 2, 1, 3))]:
    M = Mmat(tau, order); f = p_of(xs, 0*xs)@M.T
    e = mle(f, p_of) - np.c_[xs, 0*xs]
    em = mle(f, lambda x, y: p_of(x, y)@M.T) - np.c_[xs, 0*xs]
    b = np.hypot(e[:, 0], e[:, 1])
    print("tau", tau, order, "max|b| %.3f at x=%g  (max|bx| %.3f) ; withM max %.2e" % (
        b.max(), xs[b.argmax()], np.abs(e[:, 0]).max(), np.abs(em).max()))

# --- 4. flickering
def telegraph(n, Ttot=400.0, ton=100.0, toff=100.0):
    ivs = []
    for _ in range(n):
        on = rng.random() < ton/(ton+toff); t = 0.0; cur = []
        while t < Ttot:
            d = rng.exponential(ton if on else toff)
            if on: cur.append((t, min(t+d, Ttot)))
            t += d; on = not on
        ivs.append(cur)
    return ivs
def w_seq(iv, r, Ttot=400.0):
    b = Ttot/(4*r); w = np.zeros(4)
    for k in range(4*r):
        lo, hi = k*b, (k+1)*b
        for s, e in iv: w[k % 4] += max(0.0, min(hi, e)-max(lo, s))
    return w
def w_int(iv, Ttot=400.0):
    P = 0.05; Tu = 0.0125; w = np.zeros(4)
    for s, e in iv:
        for j in range(4):
            w[j] += np.ceil((e-j*Tu)/P - 1e-12) - np.ceil((s-j*Tu)/P - 1e-12)
    return w
NR = 2000; ivs = telegraph(NR)
p0 = p_of(np.array(0.), np.array(0.)); print("p centro", p0)
res = {}
for key in ["seq1", "seq5", "seq25", "int"]:
    W = np.array([w_seq(iv, int(key[3:])) if key != "int" else w_int(iv) for iv in ivs])
    lam = W*p0; ok = lam.sum(1) > 0
    f = lam[ok]/lam[ok].sum(1, keepdims=True)
    e = mle(f, p_of)
    sfl = np.sqrt(0.5*np.mean((e**2).sum(1)))
    res[key] = (W, ok, sfl); print("flicker", key, "n_ok", ok.sum(), "sigma_fl pop %.3f nm" % sfl,
                                  "(std-based %.3f)" % np.sqrt(0.5*(e.var(0).sum())))
# --- 5. SimuFLUX: N=100 multinomial
print("== 5. N=100")
for key in ["seq1", "seq5", "int"]:
    W, ok, sfl = res[key]; lam = W[ok]*p0; qq = lam/lam.sum(1, keepdims=True)
    c = np.array([rng.multinomial(100, q) for q in qq]).astype(float)
    c0_ = rng.multinomial(100, p0, size=len(qq)).astype(float)
    e = mle(c, p_of); e0 = mle(c0_, p_of)
    v = 0.5*(e.var(0).sum()); v0 = 0.5*(e0.var(0).sum())
    print(key, "STD %.3f ctrl %.3f sqrt(STD^2-ctrl^2) %.3f  vs pop sigma_fl %.3f  pred STD %.3f" % (
        np.sqrt(v), np.sqrt(v0), np.sqrt(max(v-v0, 0)), sfl, np.sqrt(v0+sfl**2)))
