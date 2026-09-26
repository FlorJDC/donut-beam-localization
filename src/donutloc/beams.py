# -*- coding: utf-8 -*-
"""Scalar beam intensity profiles (paraxial), in nm.

All profiles take Cartesian coordinates ``x, y`` (scalars or broadcastable arrays, nm) measured
from the beam axis and return the intensity with the same broadcast shape.

Conventions (Balzarotti et al. 2017, SI):

- LG01 donut, Balzarotti2017 Eq. S17 with A0 = 1::

      I(r) = 4 e ln2 r^2/fwhm^2 * exp(-4 ln2 r^2/fwhm^2)

  ``fwhm`` is a size parameter, not the FWHM of any curve of the donut. The ring maximum is 1,
  at ``r_ring = fwhm / (2 sqrt(ln2))``; the peak-to-peak ring diameter is
  ``fwhm/sqrt(ln2) ~= 1.201 fwhm``.
- Gaussian, Balzarotti2017 Eq. S19: ``exp(-4 ln2 r^2/fwhm^2)`` (peak 1 at r = 0).
- Quadratic zero, Balzarotti2017 Eq. S16 with the curvature of the donut at its zero
  (Eq. S20): ``4 e ln2 r^2/fwhm^2``. It matches ``lg_donut`` to first order for r << fwhm.

Residual zero ``eps`` (finite zero depth) is the intensity at r = 0 relative to the ring peak:

- ``zero_model="gaussian"`` (default): ``I_LG + eps * exp(-4 ln2 r^2/fwhm^2)``;
- ``zero_model="constant"``:           ``I_LG + eps``.

There is no renormalization (the global scale cancels in the multiplexed probabilities), so
"x % of the peak at the zero" means ``eps = x/100``.
"""

import numpy as np

_LN2 = np.log(2.0)
_ZERO_MODELS = ("gaussian", "constant")
_KINDS = ("donut", "gaussian", "quadratic")


def _a(fwhm):
    fwhm = float(fwhm)
    if not fwhm > 0:
        raise ValueError("fwhm must be > 0, got %r" % (fwhm,))
    return 4.0 * _LN2 / fwhm ** 2


def _r2(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return x * x + y * y


def ring_radius(fwhm):
    """Radius of the donut ring maximum, ``fwhm / (2 sqrt(ln2))`` (from Balzarotti2017 Eq. S17).

    The peak-to-peak ring diameter is ``2 * ring_radius(fwhm) ~= 1.201 * fwhm``.
    """
    _a(fwhm)
    return float(fwhm) / (2.0 * np.sqrt(_LN2))


def gaussian(x, y, fwhm=300.0):
    """Gaussian beam, Balzarotti2017 Eq. S19 with A0 = 1: ``exp(-4 ln2 r^2/fwhm^2)``."""
    return np.exp(-_a(fwhm) * _r2(x, y))


def quadratic(x, y, fwhm=300.0):
    """Ideal quadratic zero, Balzarotti2017 Eq. S16, with the donut curvature of Eq. S20:
    ``4 e ln2 r^2/fwhm^2`` (so ``quadratic ~ lg_donut`` for 4 ln2 r^2 << fwhm^2)."""
    return np.e * _a(fwhm) * _r2(x, y)


def lg_donut(x, y, fwhm=300.0, eps=0.0, zero_model="gaussian"):
    """LG01 donut, Balzarotti2017 Eq. S17 (ring peak = 1), plus an optional residual zero.

    Parameters
    ----------
    x, y : array_like
        Coordinates relative to the beam axis (nm), broadcastable.
    fwhm : float
        Size parameter of Eq. S17 (nm). Default 300.
    eps : float
        Residual intensity at r = 0 relative to the ring peak (>= 0). ``lg_donut(0, 0, eps=e)
        == e`` for both zero models.
    zero_model : {"gaussian", "constant"}
        ``"gaussian"``: add ``eps * exp(-4 ln2 r^2/fwhm^2)``; ``"constant"``: add ``eps``.
    """
    a = _a(fwhm)
    eps = float(eps)
    if eps < 0:
        raise ValueError("eps must be >= 0, got %r" % (eps,))
    if zero_model not in _ZERO_MODELS:
        raise ValueError("zero_model must be one of %s, got %r" % (_ZERO_MODELS, zero_model))
    r2 = _r2(x, y)
    g = np.exp(-a * r2)
    out = np.e * a * r2 * g
    if eps != 0.0:
        if zero_model == "gaussian":
            out = out + eps * g
        else:
            out = out + eps
    return out


def make_beam(kind="donut", fwhm=300.0, eps=0.0, zero_model="gaussian", power=1):
    """Return a vectorized beam function ``f(x, y) -> intensity``.

    Parameters
    ----------
    kind : {"donut", "gaussian", "quadratic"}
    fwhm : float
        Size parameter (nm).
    eps, zero_model :
        Residual zero, only for ``kind="donut"`` (see :func:`lg_donut`). A nonzero ``eps`` with
        another kind raises ``ValueError``.
    power : float
        Multiphoton / nonlinear exponent c: the returned function is ``I(x, y) ** power``
        (lambda_i proportional to I^c; see docs/literature/A_minflux_theory.md, section 6).
        ``eps`` is applied to the intensity before the power.

    The returned function carries attributes ``kind, fwhm, eps, zero_model, power``.
    """
    if kind not in _KINDS:
        raise ValueError("kind must be one of %s, got %r" % (_KINDS, kind))
    _a(fwhm)
    power = float(power)
    if not power > 0:
        raise ValueError("power must be > 0, got %r" % (power,))
    if kind != "donut" and float(eps) != 0.0:
        raise ValueError("eps is only supported for kind='donut'")
    if kind == "donut":
        # validate eps / zero_model now rather than at first call
        lg_donut(0.0, 0.0, fwhm=fwhm, eps=eps, zero_model=zero_model)

        def base(x, y):
            return lg_donut(x, y, fwhm=fwhm, eps=eps, zero_model=zero_model)
    elif kind == "gaussian":
        def base(x, y):
            return gaussian(x, y, fwhm=fwhm)
    else:
        def base(x, y):
            return quadratic(x, y, fwhm=fwhm)

    if power == 1.0:
        f = base
    else:
        def f(x, y):
            return base(x, y) ** power

    f.kind = kind
    f.fwhm = float(fwhm)
    f.eps = float(eps)
    f.zero_model = zero_model
    f.power = power
    return f
