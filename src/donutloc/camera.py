# -*- coding: utf-8 -*-
"""Ideal camera localization as a ``p_fn`` model (reference for MINFLUX comparisons).

Model (Balzarotti2017 Eq. S59-S63; docs/literature/A_minflux_theory.md section 4.4):

* Gaussian PSF of standard deviation ``sigma_psf`` centred on the emitter position r
  (Eq. S59).
* A square window of ``n_pix x n_pix`` pixels of side ``pixel`` (nm), centred on (0, 0).  The
  signal probability of pixel i is the PSF integrated over the pixel, a product of two erf
  differences (Eq. S60):

      q_i(r) = 1/4 [erf((x_hi - x)/(sqrt2 s)) - erf((x_lo - x)/(sqrt2 s))]
                   [erf((y_hi - y)/(sqrt2 s)) - erf((y_lo - y)/(sqrt2 s))].

* **Renormalisation inside the window (declared choice):** ``p_i^(0) = q_i / sum_j q_j``, i.e.
  ``N`` counts the photons detected *inside the window* and the signal falling outside is
  discarded (conditioning on the ROI).  For the default 9 x 9 window of 100 nm pixels and
  sigma_psf = 100 nm the discarded fraction at r = 0 is ~1e-5.  Far outside the window all
  ``q_i`` underflow and the model returns NaN (as ``photons.probabilities`` does).
* Uniform background per pixel (Eq. S61-S63).  Two conventions for the SBR:

  - ``sbr_convention="total"`` (default; Masullo's total-signal / total-background, the same
    convention as ``photons`` for MINFLUX, Balzarotti2017 Eq. S29-S30):
    ``p_i = SBR/(SBR+1) p_i^(0) + 1/(SBR+1) * 1/K``;
  - ``sbr_convention="per_pixel"``: Balzarotti's SBR_c = total signal / background *per pixel*
    (Eq. S61), ``p_i = 1/(K+SBR_c) + SBR_c/(K+SBR_c) p_i^(0)`` (Eq. S62-S63).

  They are related by ``SBR_c = K * SBR_total`` (note A section 9): e.g. SBR_c = 500 with
  K = 81 is SBR_total = 6.17.

The output has shape ``(..., K)`` with ``K = n_pix**2``, pixels in row-major order (y index
outer, x index inner), so ``fisher.crb`` applies unchanged (Eq. S11, S13).
"""

import numpy as np
from scipy.special import erf, erfc

from . import fisher

__all__ = ["pixel_edges", "make_camera_model", "crb_camera_ideal", "crb_camera"]

_CONVENTIONS = ("total", "per_pixel")


def pixel_edges(pixel=100.0, n_pix=9):
    """Pixel edges (nm) of a window of ``n_pix`` pixels of side ``pixel`` centred on 0."""
    pixel = float(pixel)
    n_pix = int(n_pix)
    if not pixel > 0:
        raise ValueError("pixel must be > 0, got %r" % (pixel,))
    if n_pix < 1:
        raise ValueError("n_pix must be >= 1, got %r" % (n_pix,))
    return pixel * (np.arange(n_pix + 1) - 0.5 * n_pix)


def _segment(lo, hi):
    """0.5*(erf(hi) - erf(lo)) for lo <= hi, computed with erfc in the tails (no cancellation)."""
    both_pos = lo >= 0
    both_neg = hi <= 0
    with np.errstate(invalid="ignore"):
        mid = 0.5 * (erf(hi) - erf(lo))
        pos = 0.5 * (erfc(lo) - erfc(hi))
        neg = 0.5 * (erfc(-hi) - erfc(-lo))
    return np.where(both_pos, pos, np.where(both_neg, neg, mid))


def _signal_fraction(r, edges, sigma_psf):
    """Unnormalised pixel fractions q_i(r), shape (..., n_pix**2) (Balzarotti2017 Eq. S60)."""
    r = np.asarray(r, dtype=float)
    if r.shape[-1:] != (2,):
        raise ValueError("r must have shape (..., 2), got %r" % (r.shape,))
    k = 1.0 / (np.sqrt(2.0) * sigma_psf)
    x = r[..., 0, None]
    y = r[..., 1, None]
    qx = _segment((edges[:-1] - x) * k, (edges[1:] - x) * k)       # (..., n)
    qy = _segment((edges[:-1] - y) * k, (edges[1:] - y) * k)       # (..., n)
    q = qy[..., :, None] * qx[..., None, :]                          # (..., n_y, n_x)
    return q.reshape(q.shape[:-2] + (-1,))


def _sbr_total(sbr, K, sbr_convention):
    """SBR in the total/total convention (None -> no background)."""
    if sbr_convention not in _CONVENTIONS:
        raise ValueError("sbr_convention must be one of %s, got %r" % (_CONVENTIONS, sbr_convention))
    if sbr is None or np.isinf(float(sbr)):
        return None
    s = float(sbr)
    if not s >= 0:
        raise ValueError("sbr must be >= 0 (or None/inf), got %r" % (sbr,))
    return s if sbr_convention == "total" else s / K


def make_camera_model(sigma_psf=100.0, pixel=100.0, n_pix=9, sbr=None, sbr_convention="total"):
    """Ideal-camera model ``p_fn(r) -> (..., n_pix**2)`` (Balzarotti2017 Eq. S59-S63).

    Parameters
    ----------
    sigma_psf : float, PSF standard deviation (nm), Eq. S59.
    pixel : float, pixel side (nm).
    n_pix : int, the window is ``n_pix x n_pix`` pixels centred on (0, 0).
    sbr : float or None, signal-to-background ratio (None/inf: no background; 0: pure background).
    sbr_convention : {"total", "per_pixel"}; see the module docstring (SBR_c = K SBR_total).

    The returned callable carries attributes ``sigma_psf, pixel, n_pix, K, edges, sbr,
    sbr_convention, sbr_total``.
    """
    sigma_psf = float(sigma_psf)
    if not sigma_psf > 0:
        raise ValueError("sigma_psf must be > 0, got %r" % (sigma_psf,))
    edges = pixel_edges(pixel, n_pix)
    K = int(n_pix) ** 2
    s_tot = _sbr_total(sbr, K, sbr_convention)

    def p_fn(r):
        q = _signal_fraction(r, edges, sigma_psf)
        with np.errstate(invalid="ignore", divide="ignore"):
            p0 = q / q.sum(axis=-1, keepdims=True)
        if s_tot is None:
            return p0
        if s_tot == 0.0:
            return np.full(q.shape, 1.0 / K)
        return (s_tot / (s_tot + 1.0)) * p0 + (1.0 / (s_tot + 1.0)) / K

    p_fn.sigma_psf = sigma_psf
    p_fn.pixel = float(pixel)
    p_fn.n_pix = int(n_pix)
    p_fn.K = K
    p_fn.edges = edges
    p_fn.sbr = sbr
    p_fn.sbr_convention = sbr_convention
    p_fn.sbr_total = np.inf if s_tot is None else s_tot
    return p_fn


def crb_camera_ideal(sigma_psf, N):
    """Continuum, background-free camera limit sigma_PSF / sqrt(N) (per axis, nm)
    (Balzarotti2017 p. 1: sigma_PSF = 100 nm, N = 400 -> 5 nm)."""
    return np.asarray(sigma_psf, dtype=float) / np.sqrt(np.asarray(N, dtype=float))


def crb_camera(sigma_psf, N, pixel=100.0, n_pix=9, sbr=None, r=(0.0, 0.0), sbr_convention="total"):
    """CRB sqrt(tr(F^-1)/2) (nm) of the ideal pixelated camera at emitter position ``r``,
    computed with ``fisher.crb`` (Eq. S11, S13) on :func:`make_camera_model`.

    Returns a float for a single position ``r`` of shape ``(2,)``, else an array of shape
    ``r.shape[:-1]``."""
    p = make_camera_model(sigma_psf, pixel, n_pix, sbr=sbr, sbr_convention=sbr_convention)
    out = fisher.crb(p, np.asarray(r, dtype=float), N)
    out = np.asarray(out, dtype=float)
    return float(out) if out.ndim == 0 else out
