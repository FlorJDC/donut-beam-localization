# -*- coding: utf-8 -*-
"""Closed-form Cramer-Rao bounds (CRB) at the centre of the MINFLUX TCP and in 1D.

Derivation: ``docs/derivations/crb_tcp_center.md`` (every formula below is derived there
step by step and cross-checked numerically by ``scripts/verify_crb_closed_forms.py``).

Conventions (see the round-1 plan):
  * TCP: 3 donuts on a circle of DIAMETER ``L`` (radius R = L/2) plus one at the centre.
  * LG donut (Balzarotti 2017, Eq. S17), peak 1:
        I(r) = 4 e ln2 r^2/fwhm^2 exp(-4 ln2 r^2/fwhm^2);   fwhm = inf -> quadratic zero.
  * CRB metric: sigma = sqrt(tr(F^-1)/2), F = N sum_i grad p_i grad p_i^T / p_i (Eq. S11, S13).
  * Residual zero ``eps`` (intensity at r=0 relative to the ring peak):
        "constant": I_LG + eps ;  "gaussian": I_LG + eps exp(-4 ln2 r^2/fwhm^2).
  * Background: Eq. S30 with a fixed SBR.

Key dimensionless number: x = L^2 ln2 / fwhm^2 (x = a R^2 with a = 4 ln2/fwhm^2).

Units: nm. All functions broadcast over numpy arrays; scalars in -> floats out.
Only numpy is used.
"""

import numpy as np

__all__ = [
    "x_param",
    "crb_tcp_center_point",
    "crb_tcp_center_limit",
    "crb_tcp_center_limit_axes",
    "limit_to_point_ratio",
    "crb_tcp_center_eps",
    "sbr_eff_constant_pedestal",
    "crossover_radius",
    "sbr_center_vs_L",
    "crb_1d_center",
]

_LN2 = np.log(2.0)
_E = np.e


def _out(v):
    v = np.asarray(v, dtype=float)
    return float(v) if v.ndim == 0 else v


def x_param(L, fwhm=np.inf):
    """x = L^2 ln2 / fwhm^2 (0 for the quadratic limit fwhm = inf)."""
    L = np.asarray(L, dtype=float)
    fwhm = np.asarray(fwhm, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.where(np.isinf(fwhm), 0.0, L ** 2 * _LN2 / np.where(np.isinf(fwhm), 1.0, fwhm) ** 2)
    return x


def _sbr_factor(sbr, K=4):
    """sqrt((1+1/SBR)(1+(K-1)/(K SBR))) -- Eq. S31 factor; 1 for SBR = inf."""
    sbr = np.asarray(sbr, dtype=float)
    inv = np.where(np.isinf(sbr), 0.0, 1.0 / np.where(np.isinf(sbr), 1.0, sbr))
    return np.sqrt((1.0 + inv) * (1.0 + (K - 1.0) / K * inv))


def _norm_ref(sbr_ref):
    if sbr_ref == "total":
        return "beam"
    if sbr_ref not in ("beam", "lg"):
        raise ValueError("sbr_ref must be 'beam' (= 'total') or 'lg'")
    return sbr_ref


def _check_x(x):
    if np.any(np.asarray(x) >= 1.0):
        raise ValueError("L^2 ln2/fwhm^2 must be < 1 (peripheral zeros inside the ring maxima)")


def crb_tcp_center_point(L, N, fwhm=np.inf, sbr=np.inf, power=1):
    """POINT value of the CRB exactly at r = 0 (Eq. S27; Eq. S31 with background).

        sigma = L / (2 c sqrt(2N)) * (1 - L^2 ln2/fwhm^2)^-1 * sqrt((1+1/SBR)(1+3/(4 SBR)))

    Without background the central exposure has p = 0 at r = 0 and its 0/0 Fisher term is
    EXCLUDED (the Eq. S27 convention); for power c = 1 this differs from the r -> 0 limit
    (see ``crb_tcp_center_limit``). ``power`` = c is the multiphoton exponent (lambda ~ I^c);
    c = 1 reproduces Eq. S27/S31 and the 1/c law is derived in the .md, section (e).
    """
    x = x_param(L, fwhm)
    _check_x(x)
    if np.any(np.asarray(power) < 1):
        raise ValueError("power must be >= 1")
    s = np.asarray(L, dtype=float) / (2.0 * np.asarray(power, float) * np.sqrt(2.0 * np.asarray(N, float)))
    return _out(s / (1.0 - x) * _sbr_factor(sbr))


def limit_to_point_ratio(L, fwhm=np.inf):
    """rho = sigma_lim / sigma_point without background, power 1:

        rho^2 = (3 g^2 + e^x) / (3 g^2 + 2 e^x),   g = 1 - x,  x = L^2 ln2/fwhm^2.

    Quadratic limit (x -> 0): rho = 2/sqrt(5) = 0.894427...
    """
    x = x_param(L, fwhm)
    _check_x(x)
    g2 = (1.0 - x) ** 2
    ex = np.exp(x)
    return _out(np.sqrt((3.0 * g2 + ex) / (3.0 * g2 + 2.0 * ex)))


def crb_tcp_center_limit(L, N, fwhm=np.inf, power=1):
    """LIMIT r -> 0 of the CRB at the TCP centre, no background (definition of the acceptance
    test). For c = 1 it is the same for every direction of approach:

        sigma_lim^2 = L^2 / (8 N g^2) * (3 g^2 + e^x) / (3 g^2 + 2 e^x)

    (= L^2/(10 N) for a quadratic zero, i.e. sigma_lim/sigma_S27 = 2/sqrt(5)).
    For c >= 2 the central term vanishes as r^(2c-2) and the limit equals the point value.
    """
    p = crb_tcp_center_point(L, N, fwhm, np.inf, power)
    if np.all(np.asarray(power) > 1):
        return p
    rho = limit_to_point_ratio(L, fwhm)
    rho = np.where(np.asarray(power) > 1, 1.0, rho)
    return _out(p * rho)


def crb_tcp_center_limit_axes(L, N, fwhm=np.inf):
    """(sigma_parallel, sigma_perp) of the r -> 0 limit covariance (no background, c = 1):
    parallel/perpendicular to the direction of approach. The ellipse is anisotropic,
    sigma_par/sigma_perp = sqrt(3 g^2/(3 g^2 + 2 e^x)) (= sqrt(3/5) quadratic), although
    sqrt(tr/2) is the same for every direction.
    """
    x = x_param(L, fwhm)
    _check_x(x)
    L = np.asarray(L, float)
    N = np.asarray(N, float)
    alpha = 8.0 * N * (1.0 - x) ** 2 / L ** 2
    beta = 16.0 * N * np.exp(x) / (3.0 * L ** 2)
    return _out(1.0 / np.sqrt(alpha + beta)), _out(1.0 / np.sqrt(alpha))


def _beam_center_values(x, eps, zero_model):
    """Normalised beam profile values used by the master formula (units: ring peak = 1,
    lengths in units of R): returns (hR, dhR_times_R2, h0), with
      hR  = h(R^2),  dhR_times_R2 = R^2 h'(R^2),  h0 = h(0),   h(u) as a function of u = r^2.
    LG part: h = e a u exp(-a u)  ->  hR = e x e^-x,  R^2 h' = e x (1 - x) e^-x.
    """
    ex = np.exp(-x)
    hR = _E * x * ex
    dh = _E * x * (1.0 - x) * ex
    if zero_model == "constant":
        h0 = eps
        hR = hR + eps
    elif zero_model == "gaussian":
        h0 = eps
        hR = hR + eps * ex
        dh = dh - eps * x * ex          # R^2 d/du [eps exp(-a u)] = -eps a R^2 e^-x
    else:
        raise ValueError("zero_model must be 'constant' or 'gaussian'")
    return hR, dh, h0


def crb_tcp_center_eps(L, N, fwhm, eps, sbr=np.inf, zero_model="constant", sbr_ref="beam"):
    """EXACT CRB at r = 0 with a residual zero ``eps`` (both models) and optional background.

    Master formula (.md, section d): for any radial beam h(u), u = r^2, the TCP Fisher matrix at
    the centre is isotropic,
        F = 6 N R^2 h'(R^2)^2 / [ (h(R^2) + b) (3 h(R^2) + h(0) + 4 b) ] * Id,
    b = background per exposure in the same units. ``sbr_ref`` fixes what "signal" means:
      * "beam" (default; alias "total"): SBR of Eq. S29/S30 measured against the FULL beam
        signal INCLUDING the pedestal (b = sum_j I_j(0) / (4 SBR)). This is the convention of
        donutloc.photons.probabilities(sbr=...) with a beam that has eps (Eq. S30 composed on
        top of the eps-beam); for the constant model
        1/SBR_eff = 1/SBR + (4 eps/sum_j I_LG,j(0)) (1 + 1/SBR).
      * "lg": SBR measured against the pure LG signal (b = sum_j I_LG,j(0) / (4 SBR)); for the
        constant model this gives 1/SBR_eff = 1/SBR + 4 eps / sum_j I_LG,j(0).
    The "gaussian" pedestal is NOT exactly an effective background (its pedestal differs per
    exposure: eps at the centre, eps e^-x at the periphery, and it also changes the gradient);
    the master formula is nevertheless exact for it.
    For eps > 0 or finite SBR the central p is > 0, so point value = r -> 0 limit (continuous).
    For eps = 0, SBR = inf it returns the POINT value (Eq. S27).
    """
    fwhm_a = np.asarray(fwhm, float)
    if np.any(np.isinf(fwhm_a)):
        raise ValueError("eps is defined relative to the ring peak: fwhm must be finite")
    if np.any(np.asarray(eps) < 0):
        raise ValueError("eps must be >= 0")
    sbr_ref = _norm_ref(sbr_ref)
    x = x_param(L, fwhm)
    _check_x(x)
    eps = np.asarray(eps, float)
    hR, dh, h0 = _beam_center_values(x, eps, zero_model)
    S0 = 3.0 * hR + h0
    sbr = np.asarray(sbr, float)
    inv = np.where(np.isinf(sbr), 0.0, 1.0 / np.where(np.isinf(sbr), 1.0, sbr))
    if sbr_ref == "beam":
        b = S0 * inv / 4.0
    else:
        b = 3.0 * _E * x * np.exp(-x) * inv / 4.0
    L = np.asarray(L, float)
    R2 = (L / 2.0) ** 2
    # F = 6 N R^2 h'^2 / (...)  with h' = dh / R^2
    F = 6.0 * np.asarray(N, float) * dh ** 2 / (R2 * (hR + b) * (S0 + 4.0 * b))
    return _out(1.0 / np.sqrt(F))


def sbr_eff_constant_pedestal(L, fwhm, eps, sbr=np.inf, sbr_ref="beam"):
    """Effective SBR at r = 0 that makes the constant-pedestal model identical to Eq. S30:
      pedestal alone:  SBR_eps = sum_j I_LG,j(0) / (4 eps) = 3 e x e^-x / (4 eps)
      "lg":            1/SBR_eff = 1/SBR + 1/SBR_eps
      "beam"/"total":  1 + 1/SBR_eff = (1 + 1/SBR)(1 + 1/SBR_eps)
    Plugging SBR_eff into Eq. S31 reproduces crb_tcp_center_eps(..., "constant") exactly.
    """
    x = x_param(L, fwhm)
    eps = np.asarray(eps, float)
    sig = 3.0 * _E * x * np.exp(-x)
    with np.errstate(divide="ignore"):
        inv_e = 4.0 * eps / sig
    sbr = np.asarray(sbr, float)
    inv_b = np.where(np.isinf(sbr), 0.0, 1.0 / np.where(np.isinf(sbr), 1.0, sbr))
    sbr_ref = _norm_ref(sbr_ref)
    if sbr_ref == "lg":
        inv = inv_b + inv_e
    else:
        inv = (1.0 + inv_b) * (1.0 + inv_e) - 1.0
    with np.errstate(divide="ignore"):
        return _out(np.where(inv == 0, np.inf, 1.0 / np.where(inv == 0, 1.0, inv)))


def crossover_radius(L, fwhm, eps=0.0, sbr=np.inf):
    """Leading-order radius r_c over which the central-exposure Fisher term switches on:
        term(r) ~= (4 e a / S0') * r^2 / (r^2 + r_c^2) * rhat rhat^T,
        r_c^2 = (eps + b) / (e a),  a = 4 ln2/fwhm^2,
    with b = S0/(4 SBR) the background per exposure ("beam" convention, constant pedestal).
    r_c = 0 only for eps = 0 and SBR = inf (the discontinuity). Valid for r_c << L.
    """
    x = x_param(L, fwhm)
    hR, dh, h0 = _beam_center_values(x, np.asarray(eps, float), "constant")
    S0 = 3.0 * hR + h0
    sbr = np.asarray(sbr, float)
    inv = np.where(np.isinf(sbr), 0.0, 1.0 / np.where(np.isinf(sbr), 1.0, sbr))
    b = S0 * inv / 4.0
    a = 4.0 * _LN2 / np.asarray(fwhm, float) ** 2
    return _out(np.sqrt((np.asarray(eps, float) + b) / (_E * a)))


def sbr_center_vs_L(L, L0, sbr0, fwhm):
    """Eq. S32 (Balzarotti 2017, p. 20): SBR at the TCP centre for fixed background,
        SBR(0, L) = (L^2/L0^2) exp(ln2 (L0^2 - L^2)/fwhm^2) SBR(0, L0).
    (Follows from sum_j I_j(0) = 3 e x e^-x ∝ L^2 exp(-L^2 ln2/fwhm^2).)
    """
    L = np.asarray(L, float)
    L0 = np.asarray(L0, float)
    fwhm = np.asarray(fwhm, float)
    if np.any(np.isinf(fwhm)):
        return _out(L ** 2 / L0 ** 2 * np.asarray(sbr0, float))
    return _out(L ** 2 / L0 ** 2 * np.exp(_LN2 * (L0 ** 2 - L ** 2) / fwhm ** 2) * sbr0)


def crb_1d_center(L, N, fwhm=np.inf, kind="donut"):
    """1D, two exposures separated by L, CRB at x = 0 (Balzarotti 2017, pp. 15-16):
      "quadratic": L/(4 sqrt N)                              (Eq. S22c)
      "donut":     L/(4 sqrt N) (1 - ln2 L^2/fwhm^2)^-1      (Eq. S22f; fwhm=inf -> S22c)
      "gaussian":  fwhm^2 / (4 ln2 L sqrt N)                 (Eq. S23c; needs finite fwhm)
    """
    L = np.asarray(L, float)
    N = np.asarray(N, float)
    if kind == "quadratic":
        return _out(L / (4.0 * np.sqrt(N)))
    if kind == "donut":
        x = x_param(L, fwhm)
        _check_x(x)
        return _out(L / (4.0 * np.sqrt(N)) / (1.0 - x))
    if kind == "gaussian":
        fwhm = np.asarray(fwhm, float)
        if np.any(np.isinf(fwhm)):
            raise ValueError("gaussian beam needs a finite fwhm")
        return _out(fwhm ** 2 / (4.0 * _LN2 * L * np.sqrt(N)))
    raise ValueError("kind must be 'donut', 'quadratic' or 'gaussian'")
