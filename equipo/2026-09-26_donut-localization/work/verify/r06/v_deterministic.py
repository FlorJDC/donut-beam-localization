# -*- coding: utf-8 -*-
"""r06 verifier: deterministic numbers (CRBs, SBRs, SimuFLUX conventions, noise-free naive MLE/LMS)."""
import json, os, sys
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcore import tcp, intens, probs, crb, fisher, mle_batch, lms_tcp, E, LN2

FW = 300.0
out = {}

# ---------- bgphys: SBR and CRB at x=25,50 ----------
L, N = 50.0, 100
cen = tcp(L)
bg = intens(np.zeros(2), cen, FW).sum() / (4 * 10.0)
for x in (25.0, 50.0):
    r = np.array([x, 0.0])
    out["bgphys_sbr_x%d" % x] = intens(r, cen, FW).sum() / (4 * bg)
    out["bgphys_fixed_crb_x%d_nm" % x] = crb(lambda q: probs(q, cen, FW, sbr=10.0), r, N)
    out["bgphys_phys_crb_x%d_nm" % x] = crb(lambda q: probs(q, cen, FW, bg=bg), r, N)
# closed form for the SBR ratio: S(r)/S(0)
# ---------- iterative centre SBR (Eq. S32) ----------
def S0(Lk):
    return intens(np.zeros(2), tcp(Lk), FW).sum()
Ls = 150.0 * (25.0 / 150.0) ** (np.arange(4) / 3.0)
for Lm in (150.0, 25.0):
    b = S0(Lm) / 40.0
    sb = [S0(Lk) / (4 * b) for Lk in Ls]
    out["iter_L_schedule"] = list(Ls)
    out["bgphys_iter_match_L%d_sbr_all" % Lm] = sb
    out["S32_closed_form_ratio_25_over_150"] = (25 / 150.) ** 2 * np.exp(-LN2 * (25 ** 2 - 150 ** 2) / FW ** 2)

# ---------- SimuFLUX convention ----------
L2, fw2 = 75.0, 310.0
c2 = tcp(L2)
pf = lambda q: probs(q, c2, fw2)
out["sf_point_value_drop_centre"] = crb(pf, np.zeros(2), 100, drop=(3,))
for r0 in (1e-2, 1e-3, 1e-4):
    for ang in (0.0, 0.7, np.pi / 2):
        out["sf_limit_r%g_ang%.2f" % (r0, ang)] = crb(pf, r0 * np.array([np.cos(ang), np.sin(ang)]), 100)
# closed form: point value; peripheral p_i=1/3 at centre; grad p_i for LG with envelope
a = 4 * LN2 / fw2 ** 2
# own closed form: I_i(0) = e a (L/2)^2 exp(-a L^2/4) = I0; grad I_i at 0 = -2 e a exp(-q)(1-q) b_i, q = a L^2/4
q = a * L2 ** 2 / 4
b = tcp(L2, centre=False)
I0 = E * a * (L2 / 2) ** 2 * np.exp(-q)
gI = -2 * E * a * np.exp(-q) * (1 - q) * b          # (3,2)
S = 3 * I0; gS = gI.sum(0)                            # = 0
gp = gI / S                                           # grad p_i = gI/S (gS=0)
F = 100 * sum(np.outer(gp[i], gp[i]) / (1 / 3.) for i in range(3))
out["sf_point_value_closed_form"] = np.sqrt(0.5 * np.trace(np.linalg.inv(F)))
# L/sqrt(8N) (quadratic limit) for scale
out["sf_L_over_sqrt8N"] = L2 / np.sqrt(800)
# hexagon without centre (point value is ring-pattern independent?)
ch = tcp(L2, n_ring=6, centre=False)
out["sf_hexagon_point_value"] = crb(lambda qq: probs(qq, ch, fw2), np.zeros(2), 100)
# SimuFLUX-style finite difference: central difference with step eps, p+1e-4, per-axis sqrt(Finv_xx)
for step in (1.0, 0.1):
    def pfr(qq):
        return probs(qq, c2, fw2)
    g = []
    for j in range(2):
        d = np.zeros(2); d[j] = step / 2
        g.append((pfr(d) - pfr(-d)) / step)
    g = np.array(g)
    p = pfr(np.zeros(2)) + 1e-4
    Fs = 100 * sum(np.outer(g[:, i], g[:, i]) / p[i] for i in range(4))
    Fi = np.linalg.inv(Fs)
    out["sf_style_fd_step%g_xaxis" % step] = np.sqrt(Fi[0, 0])
    out["sf_style_fd_step%g_N1000_xaxis" % step] = np.sqrt(Fi[0, 0] * 100 / 1000)
# what L would give 2.8 exactly? scan
for LL in (60.0, 70.0, 75.0, 80.0, 90.0):
    cc = tcp(LL)
    out["sf_point_value_L%d" % LL] = crb(lambda qq: probs(qq, cc, fw2), np.zeros(2), 100, drop=(3,))
# e factor: SimuFLUX-form donut 4 ln2 (r/fwhm)^2 exp(-4 ln2 r^2/fwhm^2) max and power
rr = np.linspace(0, 1000, 200001)
sfd = 4 * LN2 * (rr / fw2) ** 2 * np.exp(-4 * LN2 * rr ** 2 / fw2 ** 2)
gau = np.exp(-4 * LN2 * rr ** 2 / fw2 ** 2)
out["sf_donut_max"] = sfd.max(); out["one_over_e"] = 1 / E
out["sf_power_ratio_donut_over_gauss"] = np.trapz(sfd * rr, rr) / np.trapz(gau * rr, rr)
ours = intens(np.stack([rr, 0 * rr], -1), np.zeros((1, 2)), fw2)[:, 0]
out["ours_over_sf_ratio"] = float(np.median(ours[1:] / sfd[1:]))

# ---------- naive / mis-specified MLE, noise-free ----------
cen = tcp(50.0)
R = 0.75 * 50.0

def nf_mle(p_true_counts, pf_model, x0):
    # global: my batched grid, then scipy polish; also start from truth
    r1 = mle_batch(p_true_counts[None], pf_model, R, h0=0.25, levels=10)[0]
    f = lambda v: -np.sum(p_true_counts * np.log(pf_model(np.asarray(v))))
    best = None
    for st in (r1, np.array([x0, 0.0])):
        res = minimize(f, st, method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-14, maxiter=20000))
        if np.hypot(*res.x) <= R and (best is None or res.fun < best.fun):
            best = res
    return r1, (r1 if best is None else best.x)

honest = 0.0
for eps in (0.002, 0.01):
    for x in (10.0, 20.0):
        r = np.array([x, 0.0])
        ptrue = lambda q, e=eps: probs(q, cen, FW, eps=e, sbr=10.0)
        pnaive = lambda q: probs(q, cen, FW, sbr=10.0)
        cnt = 100 * ptrue(r)
        g1, est = nf_mle(cnt, pnaive, x)
        b = np.hypot(*(est - r))
        c = crb(ptrue, r, 100)
        key = "naive_eps%s_x%d" % (str(eps).replace(".", "p"), x)
        out[key + "_mle_bias_abs"] = b
        out[key + "_mle_bias_vec"] = list(est - r)
        out[key + "_grid_vs_polish"] = float(np.hypot(*(g1 - est)))
        out[key + "_crb"] = c
        out[key + "_Neq"] = 100 * (c / b) ** 2
        out[key + "_lms_extra"] = float(lms_tcp(cnt, 50, FW, 10.0)[0, 0] - lms_tcp(100 * pnaive(r), 50, FW, 10.0)[0, 0])
        g2, h = nf_mle(cnt, ptrue, x)
        honest = max(honest, np.hypot(*(h - r)), np.hypot(*(g2 - r)))
out["naive_honest_max"] = honest
for sa in (20.0, np.inf):
    for x in (10.0, 20.0):
        r = np.array([x, 0.0])
        cnt = 100 * probs(r, cen, FW, sbr=10.0)
        pa = lambda q, s=sa: probs(q, cen, FW, sbr=s)
        g1, est = nf_mle(cnt, pa, x)
        out["naive_sbr%s_x%d_mle_bias_abs" % (sa, x)] = np.hypot(*(est - r))
        out["naive_sbr%s_x%d_mle_bias_vec" % (sa, x)] = list(est - r)
# eps=0.05 (omitted in the claim: naive MLE hits the disk edge?)
for x in (10.0, 20.0):
    r = np.array([x, 0.0])
    cnt = 100 * probs(r, cen, FW, eps=0.05, sbr=10.0)
    g1, est = nf_mle(cnt, lambda q: probs(q, cen, FW, sbr=10.0), x)
    out["naive_eps0p05_x%d_est" % x] = list(est)

json.dump({k: (float(v) if np.isscalar(v) else [float(t) for t in v]) for k, v in out.items()},
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "deterministic.json"), "w"), indent=1)
for k, v in out.items():
    print(k, v)
