# -*- coding: utf-8 -*-
"""Numerical Fisher information and Cramer-Rao bound for arbitrary excitation patterns.

Everything here works on a *model callable* ``p_fn(r)``:

* ``r`` has shape ``(..., 2)`` (positions in nm),
* the return value has shape ``(..., K)`` and sums to 1 over the last axis,
* it must be vectorized (any leading shape).

The Fisher information of a multinomial with ``N`` photons and probabilities ``p_i(r)`` is
(Balzarotti2017 Eq. S11)

    F(r) = N * sum_i  grad p_i  grad p_i^T / p_i ,

and the scalar CRB used throughout the project is (Balzarotti2017 Eq. S13, d = 2)

    sigma_CRB = sqrt( tr(F^-1) / 2 ) = sqrt( (Sigma_xx + Sigma_yy) / 2 ).

Gradients are computed with centred finite differences of step ``h`` (nm).

Point value vs. limit at a perfect zero
---------------------------------------
Without background, a perfect zero gives ``p_i = 0`` and ``grad p_i = 0`` at the zero, so the term
``(grad p_i)^2 / p_i`` is 0/0.  :func:`fisher_matrix` *excludes* every term with ``p_i <= p_min``
(default ``1e-12``).  At the exact centre of the TCP this reproduces Balzarotti2017 Eq. S27 (the
"point value").  The limit r -> 0 is different (the central term tends to a non-zero,
direction-dependent matrix); it is computed by :func:`crb_limit`, which averages the CRB over
``n_dir`` directions at distance ``r0`` -- the same definition used by ``tests/test_acceptance.py``.

``zero_policy`` (keyword of :func:`crb` and :func:`crb_map`) selects which of the two values is
returned at points where some ``p_i <= p_min``:

* ``"point"`` -- the point value above (Balzarotti2017 Eq. S27 at the TCP centre); default of
  :func:`crb`, so ``crb(p_fn, 0, N)`` is unchanged;
* ``"limit"`` -- the value is replaced by :func:`crb_limit` around that point (``r0 = 1e-3`` nm,
  ``n_dir = 12``), which makes maps continuous across a perfect zero (no isolated S27 pixel at the
  origin of a symmetric grid); default of :func:`crb_map`.

With a finite background or a residual zero ``eps > 0`` no ``p_i`` vanishes and both policies give
the same number (point value = limit, Balzarotti2017 Eq. S31).

Far field: if every exposure intensity underflows to 0 and there is no background,
``photons.probabilities`` returns NaN and :func:`crb` returns **NaN** there (not ``inf``); ``inf``
is reserved for a genuinely singular (finite) Fisher matrix.
"""

import numpy as np

__all__ = ["fisher_matrix", "crb", "crb_axes", "crb_limit", "crb_map"]


def _as_points(r):
    r = np.asarray(r, dtype=float)
    if r.shape[-1] != 2:
        raise ValueError("r must have shape (..., 2), got %r" % (r.shape,))
    return r


def fisher_matrix(p_fn, r, N, h=1e-3, p_min=1e-12):
    """Fisher information matrix F(r) = N sum_i grad p_i grad p_i^T / p_i (Balzarotti2017 Eq. S11).

    Parameters
    ----------
    p_fn : callable, ``(..., 2) -> (..., K)``, vectorized, rows sum to 1.
    r : array_like ``(..., 2)``, positions in nm.
    N : float or array broadcastable to ``r.shape[:-1]``, total detected photons.
    h : float, centred finite-difference step in nm.
    p_min : float. Terms with ``p_i <= p_min`` are **excluded** from the sum (point-value
        convention of Balzarotti2017 Eq. S27 at a perfect zero).  With ``p_min = 0`` only exact
        zeros are excluded.

    Returns
    -------
    F : ndarray ``(..., 2, 2)``.
    """
    r = _as_points(r)
    ex = np.array([h, 0.0])
    ey = np.array([0.0, h])
    # one vectorized call: stack (5, ..., 2)
    pts = np.stack([r, r + ex, r - ex, r + ey, r - ey], axis=0)
    P = np.asarray(p_fn(pts), dtype=float)
    p, pxp, pxm, pyp, pym = P[0], P[1], P[2], P[3], P[4]
    dx = (pxp - pxm) / (2.0 * h)
    dy = (pyp - pym) / (2.0 * h)
    keep = p > p_min
    inv_p = np.where(keep, 1.0 / np.where(keep, p, 1.0), 0.0)
    fxx = np.sum(dx * dx * inv_p, axis=-1)
    fxy = np.sum(dx * dy * inv_p, axis=-1)
    fyy = np.sum(dy * dy * inv_p, axis=-1)
    F = np.stack([np.stack([fxx, fxy], axis=-1), np.stack([fxy, fyy], axis=-1)], axis=-2)
    N = np.asarray(N, dtype=float)
    return F * N[..., None, None]


def _cov_from_F(F):
    """Closed-form inverse of 2x2 Fisher matrices; singular -> inf entries."""
    a = F[..., 0, 0]
    b = F[..., 0, 1]
    d = F[..., 1, 1]
    det = a * d - b * b
    with np.errstate(divide="ignore", invalid="ignore"):
        good = det > 0
        inv = np.where(good, 1.0 / np.where(good, det, 1.0), np.inf)
        sxx = np.where(good, d * inv, np.inf)
        syy = np.where(good, a * inv, np.inf)
        sxy = np.where(good, -b * inv, 0.0)
        bad = np.isnan(det)                 # undefined probabilities (e.g. far-field underflow)
        if np.any(bad):
            sxx = np.where(bad, np.nan, sxx)
            syy = np.where(bad, np.nan, syy)
            sxy = np.where(bad, np.nan, sxy)
    return sxx, sxy, syy


_ZERO_POLICIES = ("point", "limit")


def crb(p_fn, r, N, zero_policy="point", limit_r0=1e-3, limit_n_dir=12, **kw):
    """Scalar CRB sigma = sqrt(tr(F^-1)/2) in nm (Balzarotti2017 Eq. S11-S13).

    Parameters
    ----------
    p_fn, r, N : as in :func:`fisher_matrix`.
    zero_policy : {"point", "limit"}
        What to return at points where some ``p_i <= p_min`` (a perfect zero without background):
        ``"point"`` (default) keeps the point value with those terms excluded (Balzarotti2017
        Eq. S27 at the TCP centre); ``"limit"`` replaces it by :func:`crb_limit` centred on that
        point, with ``r0 = limit_r0`` and ``n_dir = limit_n_dir`` (finite-difference step
        ``limit_r0 * 1e-2`` and ``p_min = 0`` there, the acceptance-test definition).
    limit_r0, limit_n_dir : parameters of the ``"limit"`` replacement.
    **kw : forwarded to :func:`fisher_matrix` (``h``, ``p_min``); ``p_min`` (default 1e-12) is also
        the threshold that flags a zero for ``"limit"``.

    Returns shape ``r.shape[:-1]``. A singular Fisher matrix gives ``inf``; undefined
    probabilities (NaN from ``p_fn``, e.g. far-field underflow without background) give NaN.
    """
    if zero_policy not in _ZERO_POLICIES:
        raise ValueError("zero_policy must be one of %s, got %r" % (_ZERO_POLICIES, zero_policy))
    sxx, _, syy = _cov_from_F(fisher_matrix(p_fn, r, N, **kw))
    out = np.sqrt(0.5 * (sxx + syy))
    if zero_policy == "point":
        return out
    r = _as_points(r)
    p_min = kw.get("p_min", 1e-12)
    P = np.asarray(p_fn(r), dtype=float)
    zero = np.any(P <= p_min, axis=-1)                    # NaN compares False -> stays NaN
    if not np.any(zero):
        return out
    shape = r.shape[:-1]
    flat = np.array(out, dtype=float).reshape(-1)
    zf = np.broadcast_to(zero, shape).reshape(-1)
    pts = r.reshape(-1, 2)[zf]                                       # (M, 2)
    Nz = np.broadcast_to(np.asarray(N, dtype=float), shape).reshape(-1)[zf]
    ang = 2.0 * np.pi * np.arange(int(limit_n_dir)) / int(limit_n_dir)
    ring = pts[:, None, :] + float(limit_r0) * np.stack([np.cos(ang), np.sin(ang)], axis=-1)
    vals = crb(p_fn, ring, Nz[:, None], zero_policy="point", h=float(limit_r0) * 1e-2, p_min=0.0)
    flat[zf] = np.mean(vals, axis=-1)
    return flat.reshape(shape)


def crb_axes(p_fn, r, N, **kw):
    """Per-axis CRB and isotropy.

    Returns ``(sigma_x, sigma_y, isotropy)`` where ``sigma_x = sqrt(Sigma_xx)``,
    ``sigma_y = sqrt(Sigma_yy)`` and ``isotropy = sqrt(lambda_min / lambda_max)`` from the
    eigenvalues of Sigma = F^-1, i.e. min sigma_i / max sigma_i along the principal axes
    (Balzarotti2017 Eq. S14).
    """
    sxx, sxy, syy = _cov_from_F(fisher_matrix(p_fn, r, N, **kw))
    with np.errstate(invalid="ignore"):
        tr = sxx + syy
        disc = np.sqrt(np.maximum((0.5 * (sxx - syy)) ** 2 + sxy ** 2, 0.0))
        lmax = 0.5 * tr + disc
        lmin = 0.5 * tr - disc
        iso = np.where(np.isfinite(lmax) & (lmax > 0), np.sqrt(np.maximum(lmin, 0.0) / lmax), np.nan)
    return np.sqrt(sxx), np.sqrt(syy), iso


def crb_limit(p_fn, N, r_center=(0.0, 0.0), r0=1e-3, n_dir=12, p_min=0.0):
    """CRB in the limit r -> r_center, averaged over directions of approach.

    Mean of :func:`crb` over ``n_dir`` points ``r_center + r0 (cos a, sin a)``,
    ``a = 2 pi k / n_dir``, with finite-difference step ``h = r0 * 1e-2``.  This is exactly the
    definition of ``_crb_center`` in ``tests/test_acceptance.py``.  ``p_min`` defaults to 0 (only
    exact zeros excluded) so that the small-but-finite central probability at distance ``r0``
    is always kept, whatever ``r0``.
    """
    ang = 2.0 * np.pi * np.arange(n_dir) / n_dir
    rc = np.asarray(r_center, dtype=float)
    pts = rc + r0 * np.stack([np.cos(ang), np.sin(ang)], axis=-1)
    return float(np.mean(crb(p_fn, pts, N, h=r0 * 1e-2, p_min=p_min)))


def crb_map(p_fn, xs, ys, N, zero_policy="limit", **kw):
    """CRB on the grid ``xs x ys`` (nm).  Returns an array of shape ``(len(ys), len(xs))``
    (``meshgrid`` 'xy' indexing: row = y, column = x).

    ``zero_policy`` defaults to ``"limit"``: pixels on a perfect zero (e.g. the origin of a
    symmetric grid over a background-free TCP) get the r -> 0 limit instead of the Eq. S27 point
    value, so the map has no spurious isolated peak.  Pass ``zero_policy="point"`` for the R1
    behaviour.  Other keywords as in :func:`crb`."""
    X, Y = np.meshgrid(np.asarray(xs, float), np.asarray(ys, float), indexing="xy")
    return crb(p_fn, np.stack([X, Y], axis=-1), N, zero_policy=zero_policy, **kw)
