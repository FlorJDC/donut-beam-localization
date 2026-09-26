# -*- coding: utf-8 -*-
"""Verify every closed form of donutloc.closed_forms against an independent numerical Fisher.

Usage (from the project root):  python scripts/verify_crb_closed_forms.py

Prints "case | closed form | numerical | rel. error [| package]" and exits with code 1 if any
relative error exceeds 1e-6 (point values, exact formulas) or 1e-3 (r -> 0 limits, evaluated
numerically at r0 = 1e-3 nm averaged over 12 directions, as in tests/test_acceptance.py).
Rows marked "approx" (leading-order crossover formula) are informative and never fail.

If donutloc.fisher / photons / patterns / beams import, an extra column recomputes each
exact case through the package (make_model + fisher.crb / fisher.crb_limit).
Derivation: docs/derivations/crb_tcp_center.md
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from donutloc import closed_forms as cf  # noqa: E402

try:
    from donutloc import beams as _beams, fisher as _fisher, patterns as _patterns, photons as _photons
    HAVE_PKG = True
except ImportError:
    HAVE_PKG = False

TOL = {"point": 1e-6, "limit": 1e-3, "approx": None}


# ---------------- own numerical reference (independent of the package) --------------------

def _centers(L):
    ang = np.pi / 2 + 2 * np.pi * np.arange(3) / 3
    return np.vstack([np.stack([0.5 * L * np.cos(ang), 0.5 * L * np.sin(ang)], axis=1), [[0.0, 0.0]]])


def probs(r, L, fwhm, eps=0.0, zero_model="constant", sbr=np.inf, power=1, sbr_ref="beam"):
    d2 = ((np.asarray(r, float)[None, :] - _centers(L)) ** 2).sum(axis=1)
    if np.isinf(fwhm):
        lg, g = d2.copy(), np.ones(4)
    else:
        a = 4 * np.log(2) / fwhm ** 2
        g = np.exp(-a * d2)
        lg = np.e * a * d2 * g
    lam = (lg + (eps if zero_model == "constant" else eps * g)) ** power
    if not np.isinf(sbr):
        lam = lam + (lam.sum() if sbr_ref == "beam" else lg.sum()) / (4.0 * sbr)
    return lam / lam.sum()


def crb_num(r, L, N, fwhm, h=1e-3, p_min=1e-12, **kw):
    r = np.asarray(r, float)
    p = probs(r, L, fwhm, **kw)
    grad = np.array([(probs(r + h * e, L, fwhm, **kw) - probs(r - h * e, L, fwhm, **kw)) / (2 * h)
                     for e in np.eye(2)])
    m = p > p_min
    F = N * np.einsum("ik,jk->ij", grad[:, m] / p[m], grad[:, m])
    return float(np.sqrt(0.5 * np.trace(np.linalg.inv(F))))


def crb_ring(L, N, fwhm, r0, n_dir=12, **kw):
    return float(np.mean([crb_num(r0 * np.array([np.cos(t), np.sin(t)]), L, N, fwhm,
                                  h=min(r0 * 1e-2, 1e-3), p_min=0.0, **kw)
                          for t in np.linspace(0, 2 * np.pi, n_dir, endpoint=False)]))


def crb_lim_num(L, N, fwhm, **kw):
    return crb_ring(L, N, fwhm, 1e-3, **kw)


# ---------------- optional package column ------------------------------------------------

def pkg_value(kind, L, N, fwhm, eps=0.0, zero_model="constant", sbr=np.inf, power=1):
    if not HAVE_PKG:
        return None
    try:
        if np.isinf(fwhm):
            beam = _beams.make_beam("quadratic", fwhm=300.0, power=power)
        else:
            beam = _beams.make_beam("donut", fwhm=fwhm, eps=eps, zero_model=zero_model, power=power)
        p_fn = _photons.make_model(_patterns.tcp_centers(L), beam,
                                   sbr=None if np.isinf(sbr) else sbr)
        if kind == "point":
            return float(_fisher.crb(p_fn, np.zeros(2), N))
        return float(_fisher.crb_limit(p_fn, N))
    except Exception as exc:  # the package is written in parallel; never block on it
        return "err:%s" % type(exc).__name__


# ---------------- cases -------------------------------------------------------------------

rows = []


def add(name, kind, closed, numerical, pkg=None):
    rel = abs(numerical / closed - 1.0)
    rows.append((name, kind, closed, numerical, rel, pkg))


N = 100
for L in (5.0, 50.0, 100.0, 150.0):
    for f in (300.0, 360.0, np.inf):
        if L == 5.0 and f == 360.0:
            continue
        fs = "inf" if np.isinf(f) else "%g" % f
        add("point  S27  L=%g fwhm=%s" % (L, fs), "point", cf.crb_tcp_center_point(L, N, f),
            crb_num([0, 0], L, N, f), pkg_value("point", L, N, f))
        add("limit  r->0 L=%g fwhm=%s" % (L, fs), "limit", cf.crb_tcp_center_limit(L, N, f),
            crb_lim_num(L, N, f), pkg_value("limit", L, N, f))

for L in (50.0, 100.0):
    for f in (300.0, 360.0):
        for sbr in (5.0, 10.0):
            cl = cf.crb_tcp_center_point(L, N, f, sbr=sbr)
            add("point  S31  L=%g fwhm=%g SBR=%g" % (L, f, sbr), "point", cl,
                crb_num([0, 0], L, N, f, sbr=sbr), pkg_value("point", L, N, f, sbr=sbr))
            add("limit  S31  L=%g fwhm=%g SBR=%g" % (L, f, sbr), "limit", cl,
                crb_lim_num(L, N, f, sbr=sbr), pkg_value("limit", L, N, f, sbr=sbr))

for zm in ("constant", "gaussian"):
    for eps in (0.002, 0.01, 0.1):
        for L in (50.0, 100.0):
            for sbr in (np.inf, 5.0):
                for ref in ("beam", "lg"):
                    if np.isinf(sbr) and ref == "lg":
                        continue
                    cl = cf.crb_tcp_center_eps(L, N, 300.0, eps, sbr, zm, ref)
                    pk = pkg_value("point", L, N, 300.0, eps, zm, sbr) if ref == "beam" else None
                    add("eps    %-8s e=%g L=%g SBR=%g ref=%s" % (zm, eps, L, sbr, ref), "point", cl,
                        crb_num([0, 0], L, N, 300.0, eps=eps, zero_model=zm, sbr=sbr, sbr_ref=ref), pk)
            if eps == 0.01:
                cl = cf.crb_tcp_center_eps(L, N, 300.0, eps, np.inf, zm)
                add("eps    %-8s e=%g L=%g limit" % (zm, eps, L), "limit", cl,
                    crb_lim_num(L, N, 300.0, eps=eps, zero_model=zm),
                    pkg_value("limit", L, N, 300.0, eps, zm))

for c in (2, 3):
    for f in (300.0, np.inf):
        fs = "inf" if np.isinf(f) else "%g" % f
        add("multiphoton c=%d L=100 fwhm=%s point" % (c, fs), "point",
            cf.crb_tcp_center_point(100.0, N, f, power=c), crb_num([0, 0], 100.0, N, f, power=c),
            pkg_value("point", 100.0, N, f, power=c))
        add("multiphoton c=%d L=100 fwhm=%s limit" % (c, fs), "limit",
            cf.crb_tcp_center_limit(100.0, N, f, power=c), crb_lim_num(100.0, N, f, power=c),
            pkg_value("limit", 100.0, N, f, power=c))

# 1D forms (Eq. S21 numerically)
for kind, f, L in (("donut", 300.0, 100.0), ("quadratic", np.inf, 50.0), ("gaussian", 300.0, 300.0)):
    def p0(x, kind=kind, f=f, L=L):
        def I(u):
            if kind == "quadratic":
                return u ** 2
            a = 4 * np.log(2) / f ** 2
            return (np.e * a * u ** 2 if kind == "donut" else 1.0) * np.exp(-a * u ** 2)
        return I(x + L / 2) / (I(x + L / 2) + I(x - L / 2))
    hh = 1e-3
    d = (p0(hh) - p0(-hh)) / (2 * hh)
    add("1D     %s L=%g" % (kind, L), "point", cf.crb_1d_center(L, N, f, kind),
        np.sqrt(p0(0.0) * (1 - p0(0.0))) / abs(d) / np.sqrt(N))


# leading-order crossover (informative)
def approx_near_center(r, L, Nn, f, eps=0.0, sbr=np.inf):
    x = cf.x_param(L, f)
    a = 4 * np.log(2) / f ** 2
    S0 = 3 * (np.e * x * np.exp(-x) + eps) + eps
    b = 0.0 if np.isinf(sbr) else S0 / (4 * sbr)
    alpha = 1.0 / cf.crb_tcp_center_eps(L, Nn, f, eps, sbr, "constant") ** 2
    gamma = Nn * 4 * np.e * a / (S0 + 4 * b)
    rc = cf.crossover_radius(L, f, eps, sbr)
    fr = r ** 2 / (r ** 2 + rc ** 2)
    return float(np.sqrt(0.5 * (1 / alpha + 1 / (alpha + gamma * fr)))), rc


for eps, sbr in ((0.0, 1e4), (0.002, np.inf), (0.0, 10.0)):
    rc = cf.crossover_radius(100.0, 300.0, eps, sbr)
    for k in (1.0, 3.0):
        ap, _ = approx_near_center(k * rc, 100.0, N, 300.0, eps, sbr)
        add("approx crossover eps=%g SBR=%g r=%g r_c (r_c=%.3g nm)" % (eps, sbr, k, rc), "approx",
            ap, crb_ring(100.0, N, 300.0, k * rc, eps=eps, sbr=sbr))

# ---------------- print --------------------------------------------------------------------

print("donutloc closed forms vs own numerical Fisher (N=%d; tol point 1e-6, limit 1e-3)" % N)
print("package column: %s" % ("donutloc.fisher/photons/patterns/beams" if HAVE_PKG else "not importable"))
hdr = "%-58s %-6s %12s %12s %10s" % ("case", "kind", "closed form", "numerical", "rel. error")
if HAVE_PKG:
    hdr += " %12s" % "package"
print(hdr)
print("-" * len(hdr))
fail = []
for name, kind, cl, nu, rel, pk in rows:
    line = "%-58s %-6s %12.6f %12.6f %10.2e" % (name, kind, cl, nu, rel)
    if HAVE_PKG:
        line += " %12s" % ("" if pk is None else (pk if isinstance(pk, str) else "%.6f" % pk))
    tol = TOL[kind]
    bad = tol is not None and rel > tol
    if HAVE_PKG and isinstance(pk, float) and tol is not None and abs(pk / cl - 1) > max(tol, 1e-6):
        bad = True
    if bad:
        fail.append(name)
        line += "  <-- FAIL"
    print(line)

print()
print("Ratio sigma_lim/sigma_point (no background), rho^2=(3g^2+e^x)/(3g^2+2e^x):")
for L in (5.0, 50.0, 100.0, 150.0):
    print("  fwhm=300 L=%5g : %.5f   (quadratic 2/sqrt5 = %.5f)" % (L, cf.limit_to_point_ratio(L, 300.0),
                                                                   2 / np.sqrt(5)))

print()
print("Degradation CRB(eps)/CRB_S27(eps=0) at r=0, SBR=inf  [constant | gaussian] ; "
      "(ratio to the r->0 limit value in brackets, constant model)")
eps_list = (0.002, 0.01, 0.03, 0.05, 0.1, 0.15)
for f in (300.0, 360.0):
    for L in (50.0, 100.0, 150.0):
        p = cf.crb_tcp_center_point(L, N, f)
        lim = cf.crb_tcp_center_limit(L, N, f)
        cells = []
        for e in eps_list:
            c = cf.crb_tcp_center_eps(L, N, f, e, zero_model="constant")
            g = cf.crb_tcp_center_eps(L, N, f, e, zero_model="gaussian")
            cells.append("%.3f|%.3f[%.3f]" % (c / p, g / p, c / lim))
        print("  fwhm=%g L=%3g: " % (f, L) + "  ".join("e=%g:%s" % (e, s) for e, s in zip(eps_list, cells)))

print()
acc = cf.crb_tcp_center_limit(50.0, 100, 300.0)
print("acceptance: crb_center_lg_L50_N100_nm (fwhm=300, r->0 limit) = %.6f nm" % acc)
print("            (S27 point value would be %.6f nm)" % cf.crb_tcp_center_point(50.0, 100, 300.0))

print()
if fail:
    print("FAILED: %d case(s): %s" % (len(fail), "; ".join(fail)))
    sys.exit(1)
n_exact = sum(1 for r in rows if r[1] != "approx")
print("all %d exact checks pass (max rel. error point %.2e, limit %.2e)" % (
    n_exact, max(r[4] for r in rows if r[1] == "point"), max(r[4] for r in rows if r[1] == "limit")))
