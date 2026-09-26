# Round 3 — Verifier report

All my code is in `work/verify/r03/`. It is my own Fisher/MLE code: `mymle.py`, plus `../vfisher.py`, `../r02/rw.py` and `../r02/cam.py` from earlier verifier rounds. **Nothing imports `donutloc` for recomputation.** The package is called only in `c_fixes.py` and the tail of `b_vectorial.py`, and there it is the system under test: I compare its output against my own references.

Seeds are my own, never 42. When I compare Monte Carlo (MC) numbers, z = (mine − theirs)/√(SE_mine² + SE_theirs²).

| Script | What it checks |
|---|---|
| `a_analytic.py` | 36 analytic keys |
| `b_vectorial.py` | Vectorial numbers, with a direct 2-D Richards-Wolf quadrature, and the linear-polarization fix |
| `c_fixes.py` | The package fixes: `crb` limit with an N array, bootstrap, misalignment `sigma_se` |
| `d_mle.py` | MLE and LMS MC |
| `e_iter.py` | Iterative MINFLUX sweep, 10⁴ repetitions per N (log `e_iter.log`) |
| `f_mis.py` | Misalignment, 400 patterns × 200 repetitions (log `f_mis.log`) |
| `f_mis_lq.py` | Follow-up on the (L/4, 0) discrepancy (log `f_mis_lq.log`) |
| `g_break.py` | Attacks on `check_provenance` in scratch roots |
| `png/` | Rendered figures |

## X1 — Independent recomputation

### Analytic numbers (`a_analytic.py`, `b_vectorial.py`): 49 keys, 0 mismatches

- **Exponents.**
  - `crb_exponent_L` = 1.0043479 (theirs 1.0043478) over L = 5, 10, 20, 40 nm.
  - `crb_exponent_N` = −0.5 exactly.
  - Robustness: the local slopes are 1.0007, 1.0026 and 1.0104, and the fit depends on the range: [1, 10] gives 1.0002, [10, 50] gives 1.008, [25, 150] gives 1.073. So "≈1" holds only for L ≪ fwhm, which the claim states.
- **Centre CRBs.**
  - S27 point value, the r→0 limit, S31 at SBR = 10.
  - Masullo check with fwhm = 360 at L = 150.
  - LG 300 at L = 150.
  - Multiphoton c = 2.
- **Limit-to-point ratios** at L = 5 and L = 150.
- **Limit-ellipse axes** at L = 50 (1.37977 / 1.80247) and the quadratic ratio √(3/5).
- **Degradation with eps:** V9, two entries.
- **L_opt, CRB_opt and L_opt/(fwhm√eps)** for eps = 0.002, 0.05 and 0.15 (rel ≤ 3.4e-7).
- **Other analytic keys:**
  - divergence at L = 360.337 nm (1/CRB → 0 there);
  - transition scale 4.887 nm;
  - photons for 5 nm at SBR = 20 and with the r→0 limit (10.305);
  - `iterative_crb_all_photons_L25_N1000_nm`;
  - `adaptive_iter0_crb_center_nm` = 3.46995.
- **Camera:**
  - pixelated 9×9 CRB at N = 400 and N = 1000;
  - per-pixel SBR convention (4.9642 nm, 591.43 photons);
  - total-SBR alternative (437.12 photons);
  - ideal camera 3.1623 nm.
- **Vectorial:**
  - D_pp = 384.666;
  - curvature 7.04223e-5 (rel 1.5e-6);
  - LG fwhm 327.141 (matched curvature) and 320.256 (matched diameter);
  - zero depths 0 (1.6e-32), 0.371648 and 0.845241;
  - linear CRB at L = 50, 100, 150: 60.7742, 38.0613, 35.7231 (rel ≤ 2.3e-6);
  - correct-hand CRB at L = 150, opposite-hand CRB at L = 100;
  - vectorial/LG-curvature ratio at L = 150.

### Package-behaviour claims (`c_fixes.py`, `b_vectorial.py`)

- **`vectorial_linear_beam_harmonic_table`: verified.** My exact direct-quadrature values are CRB(7,3) = 62.35500 and limit = 60.77419.
  - The harmonic table with defaults (ρ_max = 1500, n_θ = 801) agrees to 1.4e-6 and 2.7e-7.
  - The test settings (ρ_max = 120, n_θ = 401) agree to 5.7e-6 and 4.5e-6.
  - `linear_method='cartesian'` raises a UserWarning and reproduces +7.43 % and −3.30 %.
  - The "~1e-9" in the claim is agreement with the package's own `mode='exact'`. Against an independent quadrature the agreement is ~1e-6, which is my quadrature floor.
- **`crb_limit_array_N_fix`: verified.**
  - `crb(p, [0,0], array([100,400]), zero_policy='limit')` returns shape (2,), [1.605096, 0.802548], equal to my limits and to `crb_limit` at each N.
  - With r of shape (2,2) and N of shape (2,) it pairs them correctly.
  - Mismatched shapes raise ValueError.
- **`bootstrap_se`: verified qualitatively.** The package bootstrap on Gaussian data gives 1.037× and 0.993× σ/(2√R) (R = 500 and 5000), matching my own bootstrap. My iterative sample at N = 250 gives 0.0278 (bootstrap) vs 0.0059 (Gaussian).
  - The N = 250 distribution is extremely heavy-tailed (kurtosis 43 in my sample), so the bootstrap SE there is itself unstable. Their 0.0137 and my 0.0278 differ by 2×.
  - The claimed "0.0137 vs 0.0058" is one realization, not a stable number.
- **`misalignment_sigma_se_fix`: verified.** With 20×200 at δ = 10, the package's `sigma_se`/`sigma_se_within` ratio is 3.1–8.3, so the between-pattern SE is now used.
  - My independent 400×200 between-pattern SEs match theirs closely, e.g. naive centre δ = 10: 0.0108 vs 0.0109; naive L/4 δ = 10: 0.184 vs 0.160.

### MLE / LMS Monte Carlo (`d_mle.py`, 10⁴ repetitions, 4×10⁴ for the bias)

| key | mine | theirs | z |
|---|---|---|---|
| mle_efficiency_center (SBR 10) | 0.9867 ± 0.0051 | 0.9912 ± 0.0051 | −0.63 |
| mle_sigma_center_sbr10_nm | 1.9339 ± 0.0100 | 1.9429 ± 0.0101 | −0.63 |
| mle_nobg_sigma_over_crb_limit_center | 0.8416 ± 0.0042 | 0.8397 ± 0.0041 | +0.33 |
| mle_nobg_sigma_over_s27_center | 0.7495 ± 0.0038 | 0.7478 ± 0.0037 | +0.33 |
| mle_nobg_bias_x_r2_nm | −0.3267 ± 0.0077 | −0.3445 ± 0.0077 | +1.63 |
| mle_sbr10_sigma_over_crb_x15_N100 | 1.520 ± 0.031 | 1.512 ± 0.029 | +0.20 |
| mle_sbr10_sigma_over_crb_x15_N1000 | 1.0067 ± 0.0058 | 1.0042 ± 0.0057 | +0.31 |
| lms_nobg_sigma_over_s27_center | 0.9969 ± 0.0049 | 1.0000 ± 0.0049 | −0.45 |
| lms_nobg_sigma_over_crb_limit_center | 1.1195 ± 0.0056 | 1.1229 ± 0.0055 | −0.45 |

The x = 15 tail is real.
- At N = 100, 2.8 % of the errors exceed 5×CRB (max 58 nm). The robust σ (from the MAD) is 0.98×CRB.
- The value depends on the protocol: the search disk has radius 2L = 100 nm. With a disk of 50 nm I get 1.456 ± 0.024. The paper must state the disk radius, and it does in the description.
- At N = 1000 the tail is gone.

### Iterative MINFLUX (`e_iter.py`: my own loop from the docstring protocol, 10⁴ repetitions per N, a different seed per N)

| N | mine | theirs | z |
|---|---|---|---|
| 250 | 1.178 ± 0.028 | 1.158 ± 0.014 | +0.64 |
| 500 | 0.7418 ± 0.0042 | 0.7359 ± 0.0042 | +0.99 |
| 1000 | 0.5110 ± 0.0029 | 0.5077 ± 0.0027 | +0.83 |
| 2000 | 0.3544 ± 0.0019 | 0.3544 ± 0.0019 | 0.00 |
| 4000 | 0.2485 ± 0.0013 | 0.2495 ± 0.0013 | −0.54 |
| 8000 | 0.1756 ± 0.0010 | 0.1763 ± 0.0010 | −0.47 |

Slopes:
- **Full range:** −0.5426 ± 0.0051, CI95 [−0.5527, −0.5335]. Theirs: −0.5365 [−0.5424, −0.5312]; z = −1.04.
- **N ≥ 500:** −0.5197 ± 0.0025, CI95 [−0.5244, −0.5147]. Theirs: −0.5148 [−0.5196, −0.5096]; z = −1.35.
- Both runs exclude −0.5 for N ≥ 500, so "slightly steeper than −1/2" is confirmed. The full-range value depends on the heavy-tailed N = 250 point.
- The ratio to the camera is 0.157–0.186 (theirs 0.158–0.183).
- This resolves W11 (r2 unclear). The best statement is: slope −0.515…−0.520 for N ≥ 500 and −0.537…−0.543 over the full range.

Extras at N = 1000:

| case | mine | theirs | z |
|---|---|---|---|
| SBR = 10 | 0.6102 ± 0.0032 | 0.6165 ± 0.0033 | −1.37 |
| adaptive [150, 25, 25, 25] | 0.4755 ± 0.0023 | 0.4795 ± 0.0027 | −1.15 |
| recenter=False (artefact) | 8.256 ± 0.044 | 8.145 ± 0.045 | +1.76 |

Iteration 0 (L = 150, 250 photons):
- σ = 4.229 (theirs 4.264 ± 0.027).
- Fraction with |err| > 12.5 nm: 0.0258 ± 0.0016 (theirs 0.0265 ± 0.0016).
- Only 94.1 % lie within 3×3.47 nm. So "~97 %, not 99 %" is correct.

### Misalignment (`f_mis.py`: my own model, 400 patterns × 200 repetitions)

**Centre: verified.**
- Naive |bias| at δ = 0/2/5/10: 0.173 / 1.477 / 3.748 / 7.786 nm. Theirs: 0.173 / 1.492 / 3.731 / 7.708; all |z| < 0.4.
- Honest |bias|: 0.173 / 0.168 / 0.182 / 0.165. Theirs: 0.173 / 0.169 / 0.177 / 0.180, which is at the MC floor of 0.164.
- Slope: 0.7559 ± 0.0095 vs 0.7546 ± 0.0093 (z = 0.09).
- Noise-free population, 4000 patterns: ⟨|b|⟩/δ = 0.746 / 0.755 / 0.784.
- σ and rmse agree (|z| ≤ 1.5).
- δ = 0 gives bit-identical honest and naive estimates.

**(L/4, 0): partly unclear.**
- Honest estimates and all σ values agree.
- Naive |bias| in my draw 1 is 1.571 / 3.980 / 9.059, against their 1.691 / 4.202 / 9.304 (z = −2.2 / −1.6 / −0.65). My slope is 0.8135 ± 0.0123 vs 0.8619 ± 0.0120 (z = −2.8).
- A second independent draw (400×200, new patterns) gives a slope of 0.836 ± 0.012.
- The noise-free population ratios are 0.794 / 0.802 / 0.871 per δ.
- Three draws give 0.814, 0.836 and 0.862. The spread is larger than the quoted SE of ~0.012, because at L/4 with δ = 10 the |bias| distribution across patterns is heavy-tailed. The population value is ≈0.83.
- So "≈0.86δ at (L/4, 0)" is within ~2.5σ, and the direction of the claim (larger than at the centre) is confirmed. But 0.862 ± 0.012 is not reproduced at its quoted precision.
- At L/4 the ratio grows with δ (0.79 → 0.87–0.93), so a straight line through the origin hides some nonlinearity.

Caveats on the other misalignment claims:
- **Honest bias at (L/4, 0), δ = 10:** 0.268–0.276 against a floor of 0.231, about 5σ above it. This is a small real finite-N MLE bias. "Honest at the MC floor" holds at the centre, which is what the claim's numbers cover. It does not hold at (L/4, 0) for δ ≥ 10.
- **Honest σ at (L/4, 0):** 1.1–2.4 % above the honest CRB for δ ≤ 10, and about 8 % above at δ = 15 in fig8. "Honest σ follows the honest CRB" is approximate off-centre.

## X2 — Consistency and figures

**Summaries against the registry.** Every key in `data/fig{1..8}_summary.json` that is also in `paper_numbers.json` is identical: 149 shared keys. There is one last-digit SE difference, for `mle_nobg_bias_x_r2_nm` in fig5 (rel 1.7e-15, harmless).
- Renamed duplicates also match, e.g. `fig5_mle_sigma_center_sbr10_nm` and `fig6_iterative_sigma_N*`.
- `fig5_mle_sbr10_sigma_over_crb_x0` and `fig5_crb_center_sbr10_nm` use the numeric limit CRB instead of the S31 point value (rel 3e-9).
- fig8 was checked after its regeneration: 07:26, `fig8_n_patterns` = 400, `_paperconfig` `n_patterns` = 400. All 65 shared misalignment keys are identical.
- fig8 also defines `misalignment_naive_bias_slope_{center,Lq}` (0.7546 and 0.8619, the same as `*_over_delta_*`). It also defines `fig8_naive_bias_slope_all_deltas_*` = 0.762 and 0.870, fitted over 7 deltas including δ = 15.

**All 8 figures exist.** I rendered each PDF to PNG and looked at them. What I checked in each:

- **fig1:** ring radii 180.2 and 192.3 nm; I/ρ² plateaus at 8.37e-5 (LG 300) and 7.04e-5; the curvature-matched LG overlaps the vectorial curve.
- **fig2:** zero depths 0 / 0.845 / 0.372; CRB against L; the ≤0.4 % agreement up to L = 150 is visible.
- **fig3:** 1.802 and 1.605 markers, 1.960 with SBR; the ratio tends to 2/√5.
- **fig4:** slope 1.004 over the shaded range 5–40 nm; divergence at 360 nm; photons for 5 nm at 400, 591, 10.3 and 15.4 (the 13.0 marker is not labelled); the axes for per-pixel against total SBR are labelled correctly.
- **fig7:** L_opt stars; the 0.78 rule holds only for eps ≲ 0.01; the 4.9 nm scale line.

Caption and figure problems:
- **fig8 caption is stale and contradicts the figure and the inbox decision (refuted).**
  - The caption says "50 patterns × 200 repetitions"; the figure and data use 400.
  - The caption says "line ~0.75 δ"; the figure draws "0.76 δ", fitted over all 7 deltas.
  - The inbox r3 decision requires ≈0.75δ at the centre and ≈0.86δ at (L/4, 0). The (L/4, 0) slope is neither drawn nor mentioned.
  - The caveat in `claims.json` `misalignment_naive_bias` still says "0.79 must be re-evaluated with the rerun".
- **fig5 caption labels the wrong panels.** It describes (a) to (c) and puts the bias at r = (2, 0) nm in (c). In the figure that bias is panel (d), and (c) is σ/bound against N. The caption also does not mention the outlier bump at x₀ = 10–20 nm (σ/CRB ≈ 1.5), which is visible in (b).
- **fig6 and fig7 captions are still worded as a plan.** fig6 says "inset or side panel …, where shown …" and fig7 says "Optional: …". The fig6 caption does not describe panel (b), σ/σ_camera from 0.183 to 0.158, or panel (c), σ per iteration against L_k with the centre CRB.
  - Nothing in them contradicts the numbers, and the artefact point in fig6 is grey and labelled as the caveat requires. They need rewriting in R4.
- **fig7(c)** has a y-axis label "σ_CRB" with no unit (nm). This is minor.
- **Backlog "Texto:" caveats against the captions:** the ones in scope are respected — per-pixel SBR 500/81, 10.3 photons in the r→0 limit, monotonic only for L < 360, 0.78 only for eps ≲ 0.01, recenter=False as an artefact, 4.9 nm as a scale. The exception is naive bias ≈0.75δ, handled above under fig8. No figure claims agreement with Caprile at F = 1.

## X3 — Checks

**`python scripts/check_provenance.py`** passes, but vacuously:

```
provenance: 0 tag(s), 0 entry(ies), 0 error(s), 0 warning(s)
claims: 37, numbers: 204, provenance entries: 0
all checks pass
```

It exits 0, but the stub `paper/` has zero tags.

**`python -m unittest tests.test_acceptance -v`:** 8 tests, OK. The hash of `test_acceptance.py` is still 7d1988…2303. The full suite `python -m unittest discover -s tests` ran 177 tests, OK.

**Attacks in scratch roots (`g_break.py`, via `--root`; the real `paper/` was not touched).**

Correctly rejected (exit 1):
- an unknown `\pnum`;
- an unresolved `\src`;
- `\pnumse` for a key without an SE;
- a bad `number_keys` entry;
- a missing reproduce file;
- `numbers.tex` out of sync;
- `\%` followed by a bad `\pnum`;
- a missing `\input` file;
- a claim with an unknown number;
- nested inputs;
- a circular `\input`.

Correctly accepted: the good case, and a commented-out bad `\pnum`.

Holes, which exit 0 when they should fail:
- **`\input sections/hidden` without braces.** TeX accepts this, but the file is not flattened, so its `\pnum` and `\src` are never checked.
- **`\src {key}` with a space** is not recognized as a tag by the agent-team regex, so the unresolved key passes.

A `\pnum` hidden inside a macro (`\newcommand{\nn}[1]{\pnum{#1}}`) gives a false error on `#1`. This errs on the safe side, but the real key escapes the check. These holes are low-risk while R4 writes in the plain `\pnum{k}\src{k}` style, but they can be closed with a regex that accepts whitespace and brace-less `\input`.

## Not resolved

- The precise value of `misalignment_naive_bias_over_delta_Lq` (0.862 ± 0.012). My two draws give 0.814 and 0.836, and the population value is ≈0.83. Pool several pattern draws, or quote "≈0.83–0.86" with an SE that reflects the between-draw spread.
- SBR = 10 for `mle_efficiency_center` still needs explicit confirmation from the author. This is carried over and does not block.

```claims
[{"status": "verified", "text": "X1-exponents: crb_exponent_L = 1.00435 (fit of the r->0 centre CRB over L=5,10,20,40 nm, N=100, fwhm=300; own Fisher 1.0043479) and crb_exponent_N = -0.5 exactly (N=100..1e4, L=50). Caveat: range-dependent (1.008 over 10-50, 1.073 over 25-150 nm); 'linear' only for L<<fwhm"},
 {"status": "verified", "text": "X1-analytic sample: 49 keys recomputed independently with 0 mismatches (rel <= 3.4e-7 for optimizations, <=2.3e-6 for vectorial quadrature): S27 1.80247, limit 1.60510, S31 1.96006, limit-ellipse axes 1.37977/1.80247 (L=50), sqrt(3/5), Masullo fwhm360 L150 3.16727, multiphoton c=2 1.76777, limit/point ratios L5/L150, V9 degradation, L_opt/crb_opt/ratio for eps=0.002/0.05/0.15, divergence 360.337, scale 4.887, photons-for-5nm (SBR20 14.157, limit 10.305), camera pixelated 5.2045/3.2916, per-pixel SBR 4.9642/591.43, total-SBR 437.12, ideal 3.1623, adaptive_iter0_crb 3.46995, crb_all_photons_L25 0.25094, vectorial D_pp 384.666, curvature 7.0422e-5, fwhm 327.141/320.256, zero depths 0/0.371648/0.845241, linear CRB 60.774/38.061/35.723, correct L150 5.33104, opposite L100 89.588, vec/LGcurv L150 0.99619"},
 {"status": "verified", "text": "X1-vectorial_linear_beam_harmonic_table: package harmonic linear beam gives CRB(7,3)=62.35509 and crb_limit=60.77420 (L=50,N=100) vs my independent direct 2D Richards-Wolf 62.35500/60.77419 (rel 1.4e-6/2.7e-7; test settings rho_max=120,n_theta=401: 5.7e-6/4.5e-6); linear_method='cartesian' warns (UserWarning) and is +7.43 %/-3.30 %. The '~1e-9' is vs the package's own mode='exact'; vs an independent quadrature ~1e-6"},
 {"status": "verified", "text": "X1-crb_limit_array_N_fix: fisher.crb(p,[0,0],array([100,400]),zero_policy='limit') -> shape (2,) [1.605096,0.802548] = my independent limits and crb_limit at each N; r (2,2) with N (2,) pairs correctly; mismatched shapes raise ValueError"},
 {"status": "verified", "text": "X1-bootstrap_se: montecarlo.bootstrap_sigma_se equals sigma/(2sqrtR) within 4 % for Gaussian errors (R=500, 5000) and agrees with my own bootstrap; for iterative MINFLUX at N=250 the bootstrap SE greatly exceeds the Gaussian one (mine 0.0278 vs 0.0059; theirs 0.0137 vs 0.0058). Caveat: at N=250 the distribution is extremely heavy-tailed (kurtosis 43 in my sample), so the bootstrap SE there is itself unstable (factor 2 between draws)"},
 {"status": "verified", "text": "X1-misalignment_sigma_se_fix: package sigma_se is now the between-pattern SE (20x200, delta=10: sigma_se/sigma_se_within = 3.1-8.3); my independent 400x200 between-pattern SEs agree with the registry's (e.g. naive centre d10 0.0108 vs 0.0109)"},
 {"status": "verified", "text": "X1-mle_efficient_center_sbr10: own MC (10^4 reps, own seed, disk radius L) sigma=1.9339+-0.0100 nm, efficiency 0.9867+-0.0051 vs registry 0.9912+-0.0051 (z=-0.63); within the 0.9-1.2 band"},
 {"status": "verified", "text": "X1-mle_superefficient_nobg: own MC sigma/CRB_limit=0.8416+-0.0042 (reg. 0.8397), sigma/S27=0.7495+-0.0038 (reg. 0.7478), sigma 1.3509 (reg. 1.3478); bias_x at r=(2,0) = -0.3267+-0.0077 nm with 4e4 reps (reg. -0.3445+-0.0077, z=+1.6)"},
 {"status": "verified", "text": "X1-mle_offcenter_sbr10_tail: sigma/CRB at r=(15,0), SBR=10, disk radius 2L: N=100 1.520+-0.031 (reg. 1.512+-0.029), N=1000 1.0067+-0.0058 (reg. 1.0042); outlier tail confirmed (2.8 % of errors >5 CRB, robust MAD sigma 0.98 CRB). Caveat: protocol-dependent (disk radius 50 nm gives 1.456+-0.024); the disk radius must be stated"},
 {"status": "verified", "text": "X1-lms_center_equals_s27: own MC LMS no background at centre sigma/S27=0.9969+-0.0049 (reg. 1.0000), sigma/CRB_limit=1.1195+-0.0056 (reg. 1.1229)"},
 {"status": "verified", "text": "X1-iterative_beats_camera: own iterative loop (docstring protocol, 10^4 reps, own seeds) N=1000 sigma=0.5110+-0.0029 (reg. 0.5077+-0.0027), ratio to camera 0.1616 (reg. 0.1606); SBR=10 0.6102+-0.0032 (reg. 0.6165, z=-1.4); camera 3.1623 and CRB all photons at L=25 0.25094 exact"},
 {"status": "verified", "text": "X1-iterative_photon_scaling (resolves W11): own sigmas at N=250..8000 agree (|z|<=1.0); slope N>=500 = -0.5197+-0.0025, CI95 [-0.5244,-0.5147] (reg. -0.5148 [-0.5196,-0.5096], z=-1.35); full range -0.5426+-0.0051 [-0.5527,-0.5335] (reg. -0.5365, z=-1.04); both exclude -0.5; ratio to camera 0.157-0.186 (reg. 0.158-0.183). Caveat: the full-range slope and the N=250 SE are driven by a heavy tail at N=250"},
 {"status": "verified", "text": "X1-iterative_adaptive_rule: schedule [150,25,25,25], sigma=0.4755+-0.0023 (reg. 0.4795+-0.0027, z=-1.15); iteration-0 real error 4.229 (reg. 4.264) vs centre CRB 3.470; fraction outside 12.5 nm 0.0258+-0.0016 (reg. 0.0265), only 94.1 % within 3*3.47 nm, so ~97 % not 99 %"},
 {"status": "verified", "text": "X1-iterative_no_recenter: recenter=False final sigma 8.256+-0.044 (reg. 8.145+-0.045) with per-iteration 4.25/3.37/3.19/8.26, i.e. it worsens in the last iteration (L=25, disk 18.75 nm < spread 37.5 nm): a search-disk artefact, as the caveat says"},
 {"status": "verified", "text": "X1-misalignment centre (400 patterns x 200 reps, own model): naive |bias| 0.173/1.477/3.748/7.786 nm for delta=0/2/5/10 (reg. 0.173/1.492/3.731/7.708, |z|<0.4); honest 0.173/0.168/0.182/0.165 at the MC floor 0.164; naive bias/delta slope 0.7559+-0.0095 (reg. 0.7546+-0.0093); noise-free population 0.746/0.755/0.784 per delta; sigma and rmse agree; delta=0 gives identical estimates"},
 {"status": "unclear", "text": "X1-misalignment_naive_bias_over_delta_Lq = 0.862+-0.012: two independent 400x200 draws give 0.8135+-0.0123 and 0.836+-0.012, and the noise-free population over 4000 patterns gives 0.794/0.802/0.871 per delta (population ~0.83). The quoted SE understates the spread between pattern draws (heavy-tailed |bias| at delta=10), and the ratio grows with delta (nonlinear). 'Larger than at the centre, ~0.83-0.86' holds; 0.862+-0.012 is not reproduced at its quoted precision. Pool draws or widen the SE"},
 {"status": "verified", "text": "X1-misalignment_loss_is_bias: naive sigma centre 1.851/1.861/1.893/2.021 vs rmse 1.853/2.182/3.457/6.300 (loss dominated by bias); honest sigma and all L/4 sigmas agree with the registry (|z|<=1.5). Caveat: at (L/4,0) the honest sigma is 1-2.4 % above the honest CRB for delta<=10 (8 % at 15), and the honest |bias| at delta=10 (0.27) is ~5 sigma above the MC floor (0.23): 'honest at floor / follows CRB' is exact only at the centre"},
 {"status": "verified", "text": "X2-summaries: all 149 keys shared between data/fig1..8_summary.json and data/paper_numbers.json are identical (one fig5 SE differs at 1.7e-15 rel); renamed duplicates (fig5_mle_sigma_center_sbr10_nm, fig5_mle_nobg_*_N100, fig6_iterative_sigma_N*, fig6_ratio_to_camera_N1000) match; fig8 regenerated at 07:26 with n_patterns=400, and its 65 shared keys match"},
 {"status": "verified", "text": "X2-figures 1-4 and 7: PDFs exist and were inspected visually; the axes, legends and marked values (1.802/1.605/1.960 nm, zero depths 0/0.845/0.372, slope 1.004 over 5-40 nm, divergence at 360 nm, 400/591/10.3/15.4 photons, L_opt stars, 0.78 only for eps<~0.01, 4.9 nm scale) match captions, numbers and backlog caveats. Minor: fig7(c) y label lacks (nm); the fig4(d) 13.0 marker is unlabelled"},
 {"status": "refuted", "text": "X2-fig8 caption (structure/figures.json): says '50 patterns x 200 repetitions' and 'line ~0.75 delta', but the figure and data use 400 patterns and the drawn line is '0.76 delta' (fit over 7 deltas incl. 15); the (L/4,0) slope (~0.86 per the r3 inbox decision) is neither drawn nor mentioned; the caveat of claims.json misalignment_naive_bias still refers to the superseded 0.79"},
 {"status": "refuted", "text": "X2-fig5 caption: describes panels (a)-(c) and puts the r=(2,0) bias (-0.34 nm) in (c), but the figure has 4 panels and that bias is panel (d); the caption also omits the MLE outlier bump (sigma/CRB ~1.5 at x0=10-20 nm) visible in (b)"},
 {"status": "unclear", "text": "X2-fig6 and fig7 captions are still worded as a plan ('inset or side panel', 'where shown', 'Optional:') and do not describe fig6 (b) sigma/sigma_camera 0.183->0.158 or (c) sigma per iteration vs L_k; no numerical contradiction (the artefact point is grey and labelled); they need rewriting in R4"},
 {"status": "verified", "text": "X3-checks: python scripts/check_provenance.py -> 'all checks pass', exit 0 (vacuous: stub paper, 0 tags, 0 entries); python -m unittest tests.test_acceptance -v -> 8 tests OK; unittest discover -> 177 tests OK; sha256 of tests/test_acceptance.py = 7d1988...2303 intact"},
 {"status": "verified", "text": "X3-check_provenance correctly fails (scratch --root copies) on: unknown \\pnum, unresolved \\src, \\pnumse without se, bad number_keys, missing reproduce file, numbers.tex out of sync, a bad \\pnum after an escaped \\%, a missing \\input file, a claim with an unknown number, nested inputs, a circular \\input; correctly ignores a commented-out \\pnum"},
 {"status": "refuted", "text": "X3-check_provenance holes: '\\input sections/hidden' (no braces, valid TeX) is not flattened, so its bad \\pnum and \\src pass (exit 0); '\\src {key}' with a space is not recognized as a tag, so an unresolved key passes; \\pnum inside a user macro gives a false error on '#1' and the real key is not checked. Low risk for plain \\pnum{k}\\src{k} writing, but should be closed (whitespace-tolerant regexes, brace-less \\input)"}]
```
