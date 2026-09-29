# -*- coding: utf-8 -*-
"""Background as a nuisance parameter: model, free-background MLE and augmented Fisher/CRB.

Model (Balzarotti2017 Eq. S28, constant background per exposure): with beam intensities
``lambda_i(r)`` (``photons.intensities``) and a background ``b >= 0`` per exposure in the same
units, conditioning on the total number of detected photons N (signal + background),

    p_i(r, b) = (lambda_i(r) + b) / sum_j (lambda_j(r) + b) .

For a fixed ``b`` this is exactly ``photons.probabilities(..., bg_per_exposure=b)``.  Here ``b`` is
also allowed to be an *unknown* parameter: the parameter vector is ``theta = (x, y, b)``.

With K = 4 exposures (TCP) the multinomial has K - 1 = 3 free probabilities, so (x, y, b) is
exactly determined: the free-background MLE solves ``p(theta_hat) = n / N`` whenever that point is
reachable with b >= 0, and it pays for the extra parameter with a larger variance.  The CRB of
(x, y) with b unknown is the (x, y) block of the inverse of the 3 x 3 Fisher matrix
(:func:`crb_free_bg`); it is never smaller than the known-b CRB of ``fisher.crb``.

Functions
---------
* :func:`bg_from_sbr` -- background per exposure giving a chosen SBR at a reference point
  (inverse of ``photons.sbr_at``, Balzarotti2017 Eq. S29).
* :func:`probabilities_bg` -- ``p(x, y, b)`` for ``theta`` of shape ``(..., 3)``.
* :func:`mle_free_bg` -- vectorized MLE of (x, y, b), b >= 0.
* :func:`fisher_matrix_bg`, :func:`crb_free_bg` -- 3 x 3 Fisher information and marginal CRB.
"""

import numpy as np

from . import photons

__all__ = ["bg_from_sbr", "probabilities_bg", "mle_free_bg", "fisher_matrix_bg", "crb_free_bg"]

_P_FLOOR = 1e-300
_MEM_BUDGET = 64 * 2 ** 20


def _pos(name, v):
    try:
        ok = np.isfinite(float(v)) and float(v) > 0
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise ValueError("%s must be a finite number > 0, got %r" % (name, v))


def bg_from_sbr(centers, beam, sbr, r_ref=(0.0, 0.0)):
    """Background per exposure ``b`` such that ``photons.sbr_at(r_ref, centers, beam, b) == sbr``:
    ``b = sum_j lambda_j(r_ref) / (K * sbr)`` (Balzarotti2017 Eq. S29 inverted)."""
    _pos("sbr", sbr)
    lam = photons.intensities(np.asarray(r_ref, float), centers, beam)
    return float(np.sum(lam) / (lam.shape[-1] * float(sbr)))


def probabilities_bg(theta, centers, beam):
    """``p_i(x, y, b) = (lambda_i + b) / sum_j (lambda_j + b)`` for ``theta`` of shape ``(..., 3)``.

    Returns shape ``(..., K)``.  ``b`` must be >= 0 (``ValueError`` otherwise); with ``b = 0`` and
    all ``lambda_i = 0`` the result is NaN (as in ``photons.probabilities``).
    """
    th = np.asarray(theta, dtype=float)
    if th.shape[-1:] != (3,):
        raise ValueError("theta must have shape (..., 3), got %r" % (th.shape,))
    b = th[..., 2]
    if np.any(~(b >= 0)) or np.any(~np.isfinite(b)):
        raise ValueError("background b must be finite and >= 0")
    lam = photons.intensities(th[..., :2], centers, beam) + b[..., None]
    with np.errstate(invalid="ignore", divide="ignore"):
        return lam / lam.sum(axis=-1, keepdims=True)


def _disk_grid(center, radius, step):
    n = int(np.floor(radius / step + 1e-9))
    g = step * np.arange(-n, n + 1)
    X, Y = np.meshgrid(g, g, indexing="xy")
    pts = np.stack([X.ravel(), Y.ravel()], axis=-1)
    d2 = np.sum(pts ** 2, axis=-1)
    inside = d2 <= radius ** 2 * (1 + 1e-12)
    pts, d2 = pts[inside], d2[inside]
    order = np.argsort(d2, kind="stable")
    return pts[order] + np.asarray(center, dtype=float)


def _logp_u(cand, centers, beam):
    """Clipped ln p_i for candidates cand (..., 3) = (x, y, u), b = u^2."""
    th = np.concatenate([cand[..., :2], cand[..., 2:3] ** 2], axis=-1)
    p = probabilities_bg(th, centers, beam)
    return np.log(np.maximum(p, _P_FLOOR))


def mle_free_bg(counts, centers, beam, search_radius, b_max, center=(0.0, 0.0), grid_step=None,
                n_b=41, tol=1e-4, mem_budget=_MEM_BUDGET):
    """Maximum-likelihood estimate of (x, y, b) with the background b >= 0 as a free parameter.

    Algorithm (vectorized over the batch, same spirit as ``estimators.mle``):

    1. Global grid: (x, y) on the square lattice of step ``grid_step`` (default
       ``search_radius / 25``) inside the disk ``|r - center| <= search_radius``, times ``n_b``
       values of ``u = sqrt(b)`` uniformly spaced in ``[0, sqrt(b_max)]`` (the square-root grid
       is finer at small b, where high-SBR data live).
    2. Refinement: 3-D pattern search in (x, y, u) (5 x 5 x 5 local grid, step halved each
       pass) until the spatial step is below ``tol`` nm; the u-step starts at the u-grid spacing
       and is halved together with the spatial step.  ``b = u^2`` so b >= 0 automatically and
       b = 0 is an interior point for u.  Spatial candidates outside the disk are rejected; u is
       not bounded above during the refinement.

    Ties are broken towards the grid point closest to ``center`` and the smallest b.  An
    all-zero count vector returns ``(center, 0)``.

    Parameters
    ----------
    counts : ``(K,)`` or ``(M, K)`` non-negative counts.
    centers, beam : pattern and beam (as in ``photons.make_model``).
    search_radius : float > 0 (nm).
    b_max : float > 0, upper end of the initial b grid (beam intensity units).
    n_b : int >= 2, number of u-grid values.
    tol : float > 0, final spatial step (nm).
    mem_budget : bytes for one (rows, G) block of the global grid.

    Returns ``(3,)`` or ``(M, 3)``: ``(x_hat, y_hat, b_hat)``.
    """
    C = np.asarray(counts, dtype=float)
    if C.ndim not in (1, 2):
        raise ValueError("counts must have shape (K,) or (M, K), got %r" % (C.shape,))
    single = C.ndim == 1
    C = np.atleast_2d(C)
    if not np.all(np.isfinite(C)) or np.any(C < 0):
        raise ValueError("counts must be finite and >= 0")
    c = np.asarray(centers, dtype=float)
    if c.ndim != 2 or c.shape[1] != 2:
        raise ValueError("centers must have shape (K, 2), got %r" % (c.shape,))
    if C.shape[1] != c.shape[0]:
        raise ValueError("counts have K = %d exposures but the pattern has K = %d"
                         % (C.shape[1], c.shape[0]))
    center = np.asarray(center, dtype=float)
    if center.shape != (2,):
        raise ValueError("center must have shape (2,)")
    _pos("search_radius", search_radius)
    _pos("b_max", b_max)
    _pos("tol", tol)
    _pos("mem_budget", mem_budget)
    if grid_step is not None:
        _pos("grid_step", grid_step)
    if int(n_b) < 2:
        raise ValueError("n_b must be >= 2, got %r" % (n_b,))
    R = float(search_radius)
    step = R / 25.0 if grid_step is None else float(grid_step)

    rg = _disk_grid(center, R, step)                                   # (Gr, 2)
    ug = np.linspace(0.0, np.sqrt(float(b_max)), int(n_b))             # (nb,)
    du = ug[1] - ug[0]
    # order: spatial (closest to centre first) outer, u inner (smallest b first) -> tie-breaking
    grid = np.concatenate([np.repeat(rg, ug.size, axis=0),
                           np.tile(ug, rg.shape[0])[:, None]], axis=1)  # (G, 3)
    logp = _logp_u(grid, c, beam)                             # (G, K)
    G = grid.shape[0]
    M = C.shape[0]
    chunk = max(1, int(np.floor(float(mem_budget) / (8.0 * G))))
    best = np.empty((M, 3))
    for s in range(0, M, chunk):
        ll = C[s:s + chunk] @ logp.T
        best[s:s + chunk] = grid[np.argmax(ll, axis=1)]

    off = np.arange(-2, 3, dtype=float)
    OX, OY, OU = np.meshgrid(off, off, off, indexing="ij")
    offs = np.stack([OX.ravel(), OY.ravel(), OU.ravel()], axis=-1)     # (125, 3)
    offs = offs[np.argsort(np.sum(offs ** 2, axis=-1), kind="stable")]  # zero offset first
    s_r, s_u = step / 2.0, du / 2.0
    r2max = R ** 2 * (1 + 1e-12)
    cur = best
    while True:
        cand = cur[:, None, :] + offs[None] * np.array([s_r, s_r, s_u])  # (M, 125, 3)
        ll = np.sum(C[:, None, :] * _logp_u(cand, c, beam), axis=-1)
        out = np.sum((cand[..., :2] - center) ** 2, axis=-1) > r2max
        ll[out] = -np.inf
        cur = cand[np.arange(M), np.argmax(ll, axis=1)]
        if s_r < tol:
            break
        s_r /= 2.0
        s_u /= 2.0
    res = np.column_stack([cur[:, :2], cur[:, 2] ** 2])
    zero = C.sum(axis=1) <= 0
    res[zero] = [center[0], center[1], 0.0]
    return res[0] if single else res


def fisher_matrix_bg(centers, beam, r, b, N, h=1e-3, hb=None):
    """3 x 3 Fisher information of (x, y, b) for N multinomial photons (Balzarotti2017 Eq. S11
    with theta = (x, y, b)):  F_kl = N sum_i d_k p_i d_l p_i / p_i.

    Centred finite differences with step ``h`` (nm) in x, y and ``hb`` in b (default
    ``1e-4 * b``, or 1e-8 if b = 0; for b = 0 a one-sided forward difference is used).
    ``r`` is ``(2,)``; returns ``(3, 3)``.
    """
    r = np.asarray(r, dtype=float)
    if r.shape != (2,):
        raise ValueError("r must have shape (2,)")
    b = float(b)
    if b < 0:
        raise ValueError("b must be >= 0")
    if hb is None:
        hb = 1e-4 * b if b > 0 else 1e-8
    th0 = np.array([r[0], r[1], b])
    E = np.diag([h, h, hb])
    p = probabilities_bg(th0, centers, beam)
    grads = []
    for k in range(3):
        if k == 2 and b - hb < 0:
            d = (probabilities_bg(th0 + E[k], centers, beam) - p) / hb
        else:
            d = (probabilities_bg(th0 + E[k], centers, beam)
                 - probabilities_bg(th0 - E[k], centers, beam)) / (2.0 * E[k, k])
        grads.append(d)
    D = np.stack(grads)                                                # (3, K)
    return float(N) * (D / p) @ D.T


def crb_free_bg(centers, beam, r, b, N, **kw):
    """CRB of the position when the background per exposure is an unknown parameter.

    Returns a dict: ``sigma`` = sqrt((S_xx + S_yy)/2) from the (x, y) block of the inverse of the
    3 x 3 Fisher matrix (same scalar convention as ``fisher.crb``), ``sigma_x``, ``sigma_y``,
    ``sigma_b`` (= sqrt(S_bb)), and ``sigma_known_b`` = the same scalar from the 2 x 2 (x, y) block
    of F itself inverted (b known; equals ``fisher.crb`` of the bg_per_exposure model).
    """
    F = fisher_matrix_bg(centers, beam, r, b, N, **kw)
    S = np.linalg.inv(F)
    Sk = np.linalg.inv(F[:2, :2])
    return {
        "sigma": float(np.sqrt(0.5 * (S[0, 0] + S[1, 1]))),
        "sigma_x": float(np.sqrt(S[0, 0])),
        "sigma_y": float(np.sqrt(S[1, 1])),
        "sigma_b": float(np.sqrt(S[2, 2])),
        "sigma_known_b": float(np.sqrt(0.5 * (Sk[0, 0] + Sk[1, 1]))),
    }
