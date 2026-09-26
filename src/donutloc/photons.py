# -*- coding: utf-8 -*-
"""Multiplexed photon model: exposure intensities, probabilities with background, sampling.

A *beam* is any vectorized callable ``beam(x, y) -> intensity`` (e.g. from
:func:`donutloc.beams.make_beam`). ``centers`` is a ``(K, 2)`` array of beam-zero positions (nm).
Emitter positions ``r`` have shape ``(..., 2)``; outputs have shape ``(..., K)``.

Background (docs/literature/A_minflux_theory.md section 1.3):

- no background (``sbr=None`` or ``inf``): ``p_i = lambda_i / sum_j lambda_j``
  (Balzarotti2017 Eq. S4);
- fixed SBR (conditioning on N; default alternative): Balzarotti2017 Eq. S30,
  ``p_i = SBR/(SBR+1) p_i^(0) + 1/(SBR+1) * 1/K``;
- fixed background per exposure ``lambda_b`` (same units as the beam intensity, i.e. relative
  to the ring peak times the beam amplitude): Balzarotti2017 Eq. S28,
  ``p_i = (lambda_i + lambda_b)/sum_j (lambda_j + lambda_b)``; the SBR is then
  position-dependent, ``SBR = sum_j lambda_j/(K lambda_b)`` (Eq. S29).
"""

import numpy as np


def _centers(centers):
    c = np.asarray(centers, dtype=float)
    if c.ndim != 2 or c.shape[1] != 2:
        raise ValueError("centers must have shape (K, 2), got %r" % (c.shape,))
    return c


def _check_bg(sbr, bg_per_exposure):
    if sbr is not None and bg_per_exposure is not None:
        raise ValueError("pass at most one of sbr and bg_per_exposure")
    if sbr is not None and not float(sbr) >= 0:
        raise ValueError("sbr must be >= 0 (or None/inf), got %r" % (sbr,))
    if bg_per_exposure is not None and not float(bg_per_exposure) >= 0:
        raise ValueError("bg_per_exposure must be >= 0, got %r" % (bg_per_exposure,))


def intensities(r, centers, beam):
    """Excitation intensity of each exposure at emitter position(s) ``r``.

    ``lambda_i(r) = beam(x - cx_i, y - cy_i)``; returns shape ``(..., K)``.
    """
    r = np.asarray(r, dtype=float)
    if r.shape[-1:] != (2,):
        raise ValueError("r must have shape (..., 2), got %r" % (r.shape,))
    c = _centers(centers)
    return beam(r[..., 0, None] - c[:, 0], r[..., 1, None] - c[:, 1])


def probabilities(r, centers, beam, sbr=None, bg_per_exposure=None):
    """Probability that a detected photon belongs to each exposure, shape ``(..., K)``.

    Parameters
    ----------
    sbr : float or None
        Fixed signal-to-background ratio (total signal / total background, Eq. S29), applied
        with Balzarotti2017 Eq. S30. ``None`` or ``inf``: no background. ``0``: pure background.
        Note: the signal used is the full beam intensity, including any residual zero ``eps``.
    bg_per_exposure : float or None
        Fixed background per exposure, Balzarotti2017 Eq. S28 (position-dependent SBR).

    Passing both raises ``ValueError``. If all ``lambda_i`` vanish and there is no background
    the result is NaN (undefined).
    """
    _check_bg(sbr, bg_per_exposure)
    lam = intensities(r, centers, beam)
    K = lam.shape[-1]
    if bg_per_exposure is not None:
        lam = lam + float(bg_per_exposure)
        return lam / lam.sum(axis=-1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        p0 = lam / lam.sum(axis=-1, keepdims=True)
    if sbr is None or np.isinf(sbr):
        return p0
    s = float(sbr)
    if s == 0.0:
        return np.full(lam.shape, 1.0 / K)
    return (s / (s + 1.0)) * p0 + (1.0 / (s + 1.0)) / K


def sbr_at(r, centers, beam, bg_per_exposure):
    """Position-dependent SBR for a fixed background per exposure, Balzarotti2017 Eq. S29:
    ``SBR(r) = sum_j lambda_j(r) / (K lambda_b)``. Returns shape ``(...)`` (inf if
    ``bg_per_exposure == 0``)."""
    lam = intensities(r, centers, beam)
    b = float(bg_per_exposure)
    if b < 0:
        raise ValueError("bg_per_exposure must be >= 0")
    with np.errstate(divide="ignore"):
        return lam.sum(axis=-1) / (lam.shape[-1] * b)


def make_model(centers, beam, sbr=None, bg_per_exposure=None):
    """Return ``p_fn(r) -> (..., K)`` (vectorized, sums to 1 on the last axis), the glue
    contract consumed by ``fisher``, ``estimators`` and ``montecarlo``.

    Arguments are validated immediately. ``p_fn`` carries attributes ``centers, beam, sbr,
    bg_per_exposure, K``.
    """
    _check_bg(sbr, bg_per_exposure)
    c = _centers(centers).copy()

    def p_fn(r):
        return probabilities(r, c, beam, sbr=sbr, bg_per_exposure=bg_per_exposure)

    p_fn.centers = c
    p_fn.beam = beam
    p_fn.sbr = sbr
    p_fn.bg_per_exposure = bg_per_exposure
    p_fn.K = c.shape[0]
    return p_fn


def sample_counts(p, N, size=None, rng=None, mode="multinomial"):
    """Sample photon counts per exposure.

    Parameters
    ----------
    p : array_like, shape (K,)
        Probabilities (non-negative, summing to 1 within 1e-8; renormalized exactly).
    N : int or float
        ``"multinomial"``: exact total photon number (non-negative integer).
        ``"poisson"``: mean total photon number; counts are independent Poisson(N p_i).
    size : int or None
        Number of repetitions. ``None`` returns shape ``(K,)``, else ``(size, K)``.
    rng : numpy.random.Generator, int, or None
        Generator or seed; ``None`` uses a fresh unseeded ``default_rng()``.
    mode : {"multinomial", "poisson"}

    Returns integer array (int64).
    """
    p = np.asarray(p, dtype=float)
    if p.ndim != 1:
        raise ValueError("p must have shape (K,), got %r" % (p.shape,))
    if np.any(p < 0) or not np.all(np.isfinite(p)):
        raise ValueError("p must be finite and non-negative")
    tot = p.sum()
    if abs(tot - 1.0) > 1e-8:
        raise ValueError("p must sum to 1, got %r" % (tot,))
    p = p / tot
    if not isinstance(rng, np.random.Generator):
        rng = np.random.default_rng(rng)
    shape = p.shape if size is None else (int(size),) + p.shape
    if mode == "multinomial":
        if float(N) != int(N) or int(N) < 0:
            raise ValueError("N must be a non-negative integer for multinomial, got %r" % (N,))
        return rng.multinomial(int(N), p, size=size).astype(np.int64)
    if mode == "poisson":
        if not float(N) >= 0:
            raise ValueError("N must be >= 0, got %r" % (N,))
        return rng.poisson(float(N) * np.broadcast_to(p, shape)).astype(np.int64)
    raise ValueError("mode must be 'multinomial' or 'poisson', got %r" % (mode,))
