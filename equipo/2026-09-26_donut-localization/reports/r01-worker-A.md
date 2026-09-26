# Round 1 — Worker A (physics core: beams, patterns, photons)

## What I did
Wrote the package core as Task 1 in `r01-pi.md` specifies, with tests alongside.
Files (all new):
- `src/donutloc/__init__.py`: docstring plus `__version__ = "0.1.0"`. It imports no submodules.
- `src/donutloc/beams.py`, `src/donutloc/patterns.py`, `src/donutloc/photons.py`
- `tests/test_beams.py` (10 tests), `tests/test_patterns.py` (7), `tests/test_photons.py` (11)

`pyproject.toml` already existed and I did **not** change it. It has a src layout with `packages.find where=["src"]`, which is enough.
Nothing from `docs/private/` was used. The only dependencies are numpy and scipy (scipy only in one test). Python 3.8 syntax.

## Public API (as implemented)
```
beams.ring_radius(fwhm) -> float                       # fwhm/(2 sqrt ln2)
beams.gaussian(x, y, fwhm=300.0)                       # Eq. S19
beams.quadratic(x, y, fwhm=300.0)                      # Eq. S16 with Eq. S20 curvature: 4e ln2 r^2/fwhm^2
beams.lg_donut(x, y, fwhm=300.0, eps=0.0, zero_model="gaussian")   # Eq. S17 (+ eps*exp(-a r^2) | + eps)
beams.make_beam(kind="donut", fwhm=300.0, eps=0.0, zero_model="gaussian", power=1) -> f(x, y)
    # kind in {"donut","gaussian","quadratic"}; returns I**power; attrs f.kind,.fwhm,.eps,.zero_model,.power
    # eps != 0 with a non-donut kind -> ValueError; eps is applied before the power

patterns.polygon_centers(L, M, rotation=np.pi/2, center=True) -> (M+1,2) | (M,2)
patterns.tcp_centers(L, rotation=np.pi/2, center=True) -> (4,2) center last | (3,2)
patterns.perturb_centers(centers, displacement, rng=None, mode="random_direction") -> (K,2)
    # displacement: scalar or (K,) magnitudes -> uniform random direction per center;
    # (K,2) -> explicit vectors. rng: Generator | int seed | None (unseeded default_rng()).

photons.intensities(r, centers, beam) -> (...,K)
photons.probabilities(r, centers, beam, sbr=None, bg_per_exposure=None) -> (...,K)
    # sbr None/inf: Eq. S4; sbr finite: Eq. S30 (sbr=0 -> uniform 1/K); bg_per_exposure: Eq. S28.
    # Both given -> ValueError; negative -> ValueError.
photons.sbr_at(r, centers, beam, bg_per_exposure) -> (...)   # Eq. S29, sum(lambda)/(K b); inf if b=0
photons.make_model(centers, beam, sbr=None, bg_per_exposure=None) -> p_fn(r) -> (...,K)
    # validates args at creation; copies centers; attrs p_fn.centers,.beam,.sbr,.bg_per_exposure,.K
photons.sample_counts(p, N, size=None, rng=None, mode="multinomial") -> int64 (K,) | (size,K)
    # "multinomial": N exact integer total; "poisson": N = mean total, independent Poisson(N p_i).
    # p must be 1-D, >=0, sum to 1 within 1e-8.
```
Conventions not pinned down by the plan:
- `bg_per_exposure` uses the beam-intensity units (ring peak = 1 × amplitude).
- With `sbr`, the "signal" is the full beam intensity, *including* `eps`.
- `rng=None` means an unseeded generator, following numpy semantics. Tests pass seed 42 explicitly.

## Tests (real output)
```
python -m unittest discover -s tests -p "test_beams.py"    -> Ran 10 tests ... OK
python -m unittest discover -s tests -p "test_patterns.py" -> Ran 7 tests ... OK
python -m unittest discover -s tests -p "test_photons.py"  -> Ran 11 tests ... OK
python -m unittest discover -s tests                        -> Ran 28 tests ... OK (skipped=1)  [acceptance SKIP; B/C files not present at run time]
```

## Findings
1. **Pedestal ≡ background (the plan's formula needs a qualifier).** A constant pedestal `I_LG + eps` gives exactly the same `p` as `bg_per_exposure = eps` (Eq. S28). On top of a fixed SBR, the result depends on how that SBR is defined:
   - **SBR on the pure-LG signal** (total background = ΣI_LG/SBR): `1/SBR_eff = 1/SBR + K·eps/ΣI_LG` holds exactly. This is the formula in the plan.
   - **SBR on the full beam signal** (what `photons.probabilities(..., sbr=...)` does when the beam carries `eps`): there is a cross term, `1/SBR_eff = 1/SBR + (K·eps/ΣI_LG)(1 + 1/SBR)`.
   - All three identities are tested to rtol 1e-12, for eps ∈ {0.01, 0.1}, SBR ∈ {5, 10} and 20 positions including r=0.
   - Worker C should state which SBR convention `crb_tcp_center_eps` uses.
   - For `zero_model="gaussian"` the pedestal differs between exposures, so the equivalence is not exact. I did not test it.
2. **Glue check against the acceptance test** (scratch script, not versioned). Setup: `make_model(tcp_centers(50), make_beam(fwhm=300))` with a throwaway finite-difference Fisher.
   - Point CRB at r=0, excluding p<1e-12: 1.80247191 nm. Eq. S27 gives 1.80247191 nm (rel. diff 2e-11).
   - r→0 limit (12 directions, r0=1e-3): 1.6050955713 nm. This matches `test_acceptance._crb_center(50,100,300)` to all printed digits.
   - Limit/S27 ratio at fwhm=300: 0.89050. The quadratic-case value 2/√5 is 0.89443, so at finite fwhm the ratio is *not* exactly 2/√5, which is consistent with the PI's open question. Deriving the correction belongs to Worker C. I only report the number.

## Not done / open
- `perturb_centers` implements only `mode="random_direction"` (plus explicit vectors). No Gaussian-jitter mode.
- `sample_counts` accepts only 1-D `p` (one position). It does not do per-position batches.

```claims
[
 {"text": "beams.lg_donut(ring_radius(fwhm),0,fwhm) == 1 to 1e-14 and the numerical radial maximum sits at ring_radius(fwhm)=fwhm/(2 sqrt ln2) (within 1e-4 nm) for fwhm in {200,300,360}"},
 {"text": "Peak-to-peak ring diameter 2*ring_radius(fwhm) = fwhm/sqrt(ln2) = 1.2011*fwhm"},
 {"text": "lg_donut(0,0,eps=e,zero_model=m) == e exactly for m in {gaussian,constant}, e in {0,0.002,0.05,0.15}"},
 {"text": "lg_donut/quadratic = exp(-4 ln2 r^2/fwhm^2) exactly (to 1e-14), so quadratic matches lg_donut with relative error ~4 ln2 r^2/fwhm^2 for r<<fwhm"},
 {"text": "patterns.tcp_centers(L) is bit-identical (assert_array_equal) to tests/test_acceptance.py::_tcp(L) for L in {50,100,150,37.3}"},
 {"text": "photons.probabilities sums to 1 (rtol 1e-13) with no background, sbr=5, sbr=inf, and bg_per_exposure=0.01; at r=0 without background p=(1/3,1/3,1/3,0) with p_center exactly 0"},
 {"text": "With sbr, p_i(0) equals Balzarotti2017 Eq. S30: peripheral SBR/(SBR+1)/3+1/(4(SBR+1)), center 1/(4(SBR+1)), for SBR in {1,5,10,13.6} (rtol 1e-14)"},
 {"text": "Fixed bg_per_exposure b gives Eq. S28 and equals Eq. S30 evaluated with the local SBR of Eq. S29, sum(lambda)/(K b) (rtol 1e-12)"},
 {"text": "A constant pedestal I_LG+eps gives p identical to bg_per_exposure=eps (rtol 1e-13), i.e. Eq. S30 with 1/SBR_eff = K eps/sum I_LG"},
 {"text": "With an additional fixed SBR defined on the pure-LG signal, 1/SBR_eff = 1/SBR + K eps/sum I_LG exactly; with SBR defined on the full beam signal (photons.probabilities(sbr=...) on a pedestal beam), 1/SBR_eff = 1/SBR + (K eps/sum I_LG)(1+1/SBR) exactly (rtol 1e-12, eps in {0.01,0.1}, SBR in {5,10})"},
 {"text": "sample_counts multinomial (N=100, 20000 reps, seed 42) has all rows summing to N and per-exposure means within 5 standard errors of N p; the poisson mode has means within 5 SE of N p and variance/mean within 5% of 1"},
 {"text": "make_model(tcp_centers(50), make_beam(fwhm=300)) with a finite-difference Fisher reproduces Eq. S27 at r=0 (1.802472 nm, rel. 2e-11) and the acceptance _crb_center(50,100,300)=1.6050956 nm in the r->0 limit; limit/S27 = 0.89050 at fwhm=300 versus 2/sqrt5 = 0.89443 in the quadratic limit"}
]
```
