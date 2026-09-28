"""donutloc -- simulation of fluorescent-emitter localization with structured donut beams
(MINFLUX and variants).

Units: nm. Positions are arrays of shape ``(..., 2)``.

Submodules (import them explicitly, e.g. ``from donutloc import photons``):

- ``beams``     -- scalar beam intensity profiles (LG donut, Gaussian, quadratic zero).
- ``patterns``  -- exposure-pattern geometry (TCP, regular polygons, misalignment).
- ``photons``   -- multiplexed photon model: probabilities with background, sampling.
- ``fisher``, ``estimators``, ``montecarlo``, ``closed_forms`` -- inference and CRB.
- ``background`` -- background as a free parameter: (x, y, b) MLE and 3x3 Fisher / marginal CRB.
"""

__version__ = "0.1.0"
