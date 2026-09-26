# Round 5 — Verifier (manuscript text, new r4 keys, literature)

Read: intent.md, state.json (all 125 claims, backlog incl. "Caveats para el texto"), inbox.jsonl
(7 decisions), r04-pi / r04-writer, docs/derivations/crb_tcp_center.md, docs/literature/A_*.md and
B_*.md, all paper/sections/*.tex (expanded with the actual \pnum values by my
`work/verify/r05/expand.py`), and the rendered PDF: all 13 pages at 110 dpi in
`work/verify/r05/png/p01..p13.png`, every figure inspected against its caption.
main.pdf (11:58) is newer than every section file, numbers.tex and paper_numbers.json (11:52),
so the PDF is in sync. `python scripts/check_provenance.py` -> "all checks pass" (178/178).
I did not read docs/private/. Nothing in the paper looks like unpublished private data: no
mention of p-MINFLUX, own setups or measurements. The zero depths are generic (0.002–0.15), and
the only experimental value is Balzarotti's published "<0.2 %".

Own code (none of it imports donutloc): `work/verify/r05/vpop.py` (population misalignment),
`vdemo.py` (demo keys, pedestal CRBs, LMS shrink, text arithmetic), `expand.py` (reading aid).
Outputs: `vpop_out.txt`, `vdemo_out.txt`.

---

## M1 — Text versus numbers and derivations

### Equations re-derived by hand (all correct)
I re-derived each of these myself. All match the appendix and the main text:
- master formula (A2);
- S27 point value: F = 8Ng²/L²;
- central term: 4ea/S0 = 16eˣ/(3L²);
- α, β and σ_lim² = L²/(8Ng²)·(3g²+eˣ)/(3g²+2eˣ);
- quadratic limit: β/α = 2/3, ρ² = 4/5, so ρ = 2/√5 and σ_lim² = L²/(10N);
- small-x expansion: ρ ≈ (2/√5)(1−9x/40);
- ellipse axes and their ratio √(3/5);
- background: b = S/(4 SBR) reproduces Eq. S30 and gives S31;
- transition radius r_c = (√3L/4)e^{−x/2}/√SBR;
- constant pedestal: SBR_ε = 3exe^{−x}/(4ε), and the "beam" convention combination
  1+1/SBR_eff = (1+1/SBR)(1+1/SBR_ε);
- Gaussian pedestal: h_R, R²h'_R and F_ε/F_pt;
- first-order coefficients, whose ratio is 1 + 3x²/7;
- multiphoton: F = 2Nc²g²/R².

Numeric spot checks (vdemo.py):
- L = 50: σ∥ = 1.3798, σ⊥ = 1.8025, σ_lim = 1.60510. A numeric 12-direction limit gives 1.60510;
  the point value is 1.80247.
- Masullo L = 150, fwhm = 360: 3.1673.
- Photons for 5 nm: 13.00, 13.45, 14.16, 15.37 and 17.93; limit 10.305.
- Iterative ratios: 0.1831 at N = 250 and 0.1577 at N = 8000.
- The LMS covariance at the centre is L²/(8Ng²): I derived it analytically.

Conventions respected:
- fwhm is the size parameter of Eq. S17 (ring diameter fwhm/√ln2 = 360.3 nm).
- L is the diameter.
- The SBR is total/total. The camera uses SBR_c per pixel (500/81 = 6.173).
- The CRB is √(tr F⁻¹/2) = Eq. S13.
- "Centre CRB" means the limit unless stated otherwise. Table II is explicitly relative to the
  point value.

### Ledger caveats: all present
- L exponent only for L ≪ fwhm: crb_center:76–80, open_points:25.
- 2L search radius for the x0 = 15 tail: estimators:41–47, fig 5 caption.
- "Honest at the floor" only at the centre: nonidealities:107–110. See issue 8 for a nuance.
- N = 250 bootstrap unstable: iterative:39–40, open_points:22.
- Adaptive rule ≈ 97 %; 3.47 vs 4.26 nm: iterative:47–52.
- Iterative slope excludes −0.5: iterative:31–34.
- Caprile F = 1: open_points:18.
- 0.78 only for ε ≲ 0.01: nonidealities:53–54, fig 7.
- The superseded 0.862 and 0.79 are not cited.
- The population slope is the main number (inbox r4).
- SBR = 10 for the MLE efficiency is declared a convention (inbox r1).
- The MLE superefficiency without background is reported, not hidden (inbox r1).

### Problems found (location, what is wrong, suggested wording)

1. **vectorial.tex:32–33 (overclaim, contradicted by Table I).** The text says "with the wrong
   hand the centre CRB is two orders of magnitude larger, and with linear polarization more than
   one order." Table I gives these factors:
   - opposite hand / correct: 92, 27 and 16 at L = 50, 100 and 150 nm;
   - linear / correct: 38, 11.5 and 6.7 at the same L.

   The factor is L-dependent. At L = 150 the linear case is less than one order of magnitude.

   Suggested: "with the wrong hand the centre CRB is 16–90 times larger (L = 150–50 nm), and
   with linear polarization 7–40 times larger."

2. **vectorial.tex:26–28.** "the default LG donut ... overestimates it by a few percent at the
   largest L [Fig. 2(e,f)]". Fig. 2(f) runs to L = 250 nm, where the difference is about 11 %.
   Suggested: "by 2.8 % at L = 150 nm and more at larger L (about 11 % at 250 nm, Fig. 2f)". Or
   restrict the sentence to the table range.

3. **vectorial.tex:28–31 (wrong justification).** "This is expected: the centre CRB depends only
   on the beam profile within a distance L/2 of each zero ..., where the vectorial profile is
   close to quadratic with curvature c." The overall curvature cancels in p_i: any exactly
   quadratic profile, whatever its c, gives the same CRB (g = 1). What sets the CRB is how the
   profile departs from quadratic at distance L/2 (h'_R/h_R, the g factor) together with the
   ratio of the central curvature to h_R. Matching c at a fixed peak of 1 reproduces that
   departure for these beams empirically, but curvature alone does not force it.

   Suggested: "The overall curvature cancels in p_i. The centre CRB is set by how the profile
   departs from quadratic at a distance L/2 (Eq. A2). With the peak normalized to one, the
   curvature-matched LG happens to reproduce this departure for the vectorial beam to a fraction
   of a percent."

4. **crb_center.tex:20–21.** "a well-defined isotropic 'centre CRB'". The scalar is
   direction-independent, but line 45 says correctly that the limiting ellipse is anisotropic, so
   the two sentences contradict each other. Suggested: "a well-defined, direction-independent
   'centre CRB' (the error ellipse itself is anisotropic, see below)".
   - Minor: α and β are used at lines 41 and 45 but defined only in App. A4. Add "(α, β defined in
     Eq. A4)".

5. **model.tex:50–52 (convention not stated).** "The SBR is ... at the emitter position (Eq. S29),
   so it depends on the position and on the pattern." Every simulation, however, holds the SBR
   fixed at every position (Eq. S30 with constant SBR). This is what `photons.probabilities(sbr=)`
   and the fig 3 maps do. It implies a position-dependent background: Fig. 3(b,d), the off-centre
   MLE and the misalignment study all use constant SBR, not constant background.

   Suggested: add "In all simulations we hold the SBR fixed at the value quoted, at every
   position (Eq. S30); with a fixed background instead, the SBR would fall off-centre and with L
   (Eq. S32)." Optionally add this to the scope item in open_points:7–10.

6. **estimators.tex:35–38 and fig 5.** "Away from the centre the linearized estimators are
   strongly biased toward the centre ..., and their σ falls below the CRB [Fig. 5(b)]; this is a
   consequence of the shrinkage".
   - (a) In Fig. 5(b) the mLMS has σ/CRB ≈ 1.3 > 1 for x0 ≲ 12 nm. The σ < CRB region starts at
     x0 ≳ 15 nm.
   - (b) The mLMS bias is slightly positive (away from the centre) for x0 ≈ 5–15 nm (Fig. 5a).
   - (c) The plotted LMS and mLMS include the 1/s correction (`fig_5_estimators.py` calls
     `lms_tcp(..., sbr=SBR)`). Because the sentence directly follows the uncorrected-LMS
     sentence, "the shrinkage" reads as the 1/s shrinkage.

   Suggested: "In Fig. 5(a,b) the LMS and mLMS include the 1/s correction. Beyond x0 ≈ 15 nm they
   are strongly biased toward the centre (the linearization compresses the estimate), and their σ
   falls below the CRB; this reflects the compression, not a gain in precision."

7. **estimators.tex:66 (fig 5 caption) and structure/figures.json fig5 (d).** "(d) Bias of the
   background-free MLE near the centre versus x0 ..., toward the centre". The panel shows a
   negative bias only for x0 ≲ 8 nm. It crosses zero near 9 nm and is +0.18 nm (away from the
   centre) at x0 = 12 nm.

   Suggested: "(d) ... versus x0 (N = 100): negative (toward the centre) for x0 ≲ 8 nm, with
   −0.3445 ± 0.0077 nm at r = (2,0) nm, and changing sign near x0 ≈ 9 nm." The text at
   estimators:22 (bias at (2,0)) is correct.

8. **nonidealities.tex:103–106.** "At the centre its mean |bias| at δ = 10 nm, 0.1797 ± 0.0051 nm,
   is at the MC noise floor ... 0.1643 nm." The value is 9 % (3.0 SE) above the floor. At (L/4,0)
   it is 19 % (5.4 SE) above. The r3 verifier's independent draw gave 0.165 at the centre, and
   the other δ values sit at the floor. So "at the floor" is fair as a description but
   overstated as a registry statement.

   Suggested: "is close to the MC noise floor ... (0.1643 nm; 9 % above, and at the floor in an
   independent draw), whereas at (L/4,0) it is 19 % above (0.2759 ± 0.0083 against 0.2309 nm)".

9. **nonidealities.tex:97–101.** The MC slope "differs from the population estimate in its
   weighting of the δ values and includes the MC noise floor". The noise floor would push the MC
   slope up, but the MC slope is lower (0.755 vs 0.777), and the floor is negligible at these δ.
   The difference is weighting plus the sampling of 400 patterns. The inverse-variance weights
   favour δ = 2, where |bias|/δ is smallest. Applying the MC weights to the population ratios
   gives 0.762 (my calculation); the remaining 0.007 is within the 400-pattern sampling.

   Suggested: "it differs mainly through its inverse-variance weighting, which favours δ = 2 nm
   where |bias|/δ is smallest, and through the sampling of 400 patterns."

10. **abstract.tex:23–25 and discussion.tex:21.** "a bias of about 0.7768 times the displacement".
    - Four significant figures sit badly with "about".
    - The unweighted LSQ slope is 78 % weighted by δ = 10 (weights δ²/Σδ² = 3/19/78 %), so it is
      essentially the δ = 10 ratio.
    - The ratio runs from 0.747 to 0.784 over δ = 2–10 nm.
    - The conditions are missing: mean over random patterns, L = 100 nm, SBR = 10.

    Suggested for the abstract: "produces a mean bias of 0.75–0.78 times the displacement of the
    zeros at the centre (δ = 2–10 nm, L = 100 nm)". In the discussion: "of order 0.75–0.78 δ at
    the centre and larger at (L/4,0)". "Larger off-centre" is verified only at (L/4,0).

11. **abstract.tex:17–18 and 19–22 (conditions missing).**
    - The efficiency claim holds only at L = 50 nm, N = 100, with a disk of radius L. Add
      "(L = 50 nm, N = 100)".
    - The iterative 0.5077 nm is protocol-specific. Add "with our protocol (L_k = 150 → 25 nm, no
      background)".
    - "Every number quoted is produced by a script": the 0.2 % value is literature. Say "Every
      computed number".

12. **discussion.tex:56–59.** "Iterative MINFLUX reaches a small fraction of the camera
    precision" is ambiguous: it is the error that is a fraction. And "the effects above, all
    absent from that simulation" is not quite right, because an SBR = 10 run exists (0.6165 nm).

    Suggested: "... reaches an error 0.16–0.18 times that of the ideal camera ...; the zero
    depth, misalignment and (except for one SBR = 10 run) background effects above were absent
    from that simulation ..."

13. **methods.tex:9–11.** "All MC runs use ... enough repetitions for a statistical error below
    2 %; where the dispersion is intrinsic (between random misalignment patterns) ..." Two keys
    break this and are not misalignment keys:
    - mle_nobg_bias_x_r2_nm, SE 2.2 %;
    - adaptive_frac_outside_next_radius, SE 6 %.

    Both are correctly quoted with \pnumse, but the sentence is false as written.

    Suggested: "... below 2 % for the central results; the others (misalignment quantities, the
    bias at r = (2,0) nm, the adaptive fraction) are quoted with their SE."

14. **open_points.tex:28–30.** "Population misalignment slopes ... pending an independent
    re-computation". I re-computed them this round (M2, verified). Delete the bullet, or change it
    to "verified in the final round by an independent computation (60 000 patterns, three seeds)".

15. **open_points.tex:18–21 (slight overstatement).** "do not agree with the values we
    transcribed ..., notably at filling factor 1". Per W3 (verified r2) and note B §8:
    - F = 0.5, 2 and uniform agree within 1–2.5 % (515.5/382.9/377.6 vs 508/392/380 nm);
    - F = 1 differs by 27 nm (401.4 vs 428 nm).

    Suggested: "agree within about 2 % at filling factors 0.5, 2 and uniform, but differ by 27 nm
    (401 vs 428 nm) at filling factor 1; we do not claim agreement there."

16. **crb_center.tex:67–69 (minor).** "At realistic SBR this radius is a sizeable fraction of L".
    The leading-order r_c/L is 0.19, 0.13 and 0.09 for SBR = 5, 10 and 20. The derivation notes
    that the leading-order formula fails at SBR = 10. Suggested: "(about 0.1–0.2 L for SBR = 5–20
    by the leading-order estimate)".

17. **crb_center.tex:99–101 (minor).** The fwhm = 360 values 0.9413/1.962/3.167 "reproduce" the
    Masullo column 0.94/1.96/3.16. The L = 150 value rounds to 3.17, not 3.16 (0.2 % off). The
    table legend (SI p. 9, column 2) states no SBR or N for that column; SBR = 5, N = 500 and
    fwhm = 360 are inferred.

    Suggested: "reproduces, to ≤ 0.3 %, the centre column (0.94/1.96/3.16 nm) of Supplementary
    Table 1 of Masullo et al., whose SBR = 5 and N = 500 we infer from the other columns".

18. **Formatting (numbers_tex, W1's file).**
    - SEs print with 4 significant figures, more digits than the value: "0.9912 ± 0.005131",
      "0.1606 ± 8.484×10⁻⁴", "1.943 ± 0.01006".
    - Trailing zeros are dropped in Tables I–II ("1.6", "1.06", "1.3", "2.13" next to
      four-figure entries).

    Round each SE to 1–2 significant figures and print the value to the same decimal. Use fixed
    significant figures in the tables.

19. **Attribution wording (low).**
    - model.tex:23–24: "product of the topological charge and the photon spin is +1
      [Lopez2023, Caprile2022]". This formalization is the project's own (note B §7.8). Lopez
      says "circular and along the same direction of the phase ramp". Suggested: "i.e. circular
      polarization rotating in the same sense as the phase ramp [7,8] (topological charge ×
      spin = +1 in our convention)".
    - model.tex:28: "Following Balzarotti Eq. S24": the pattern is rotated (φ₀ = π/2, whereas
      S24 has zeros at i·2π/3) and the centre is labelled 3 instead of 0. Add "up to a global
      rotation and relabelling (centre last)". This matters because "(L/4,0)" refers to this
      orientation.
    - introduction.tex:15–17: "Most analyses of MINFLUX use ... a perfect zero, an exact knowledge
      of the pattern ... [1]". Ref. [1]'s experiments used a measured PSF (SM Eq. S72–S73).
      Suggested: "Theoretical analyses of MINFLUX, starting with the Supplementary Material of
      Ref. [1], use ..."

20. **LMS shrink wording (low), estimators.tex:34–35.** "shrinks the estimate by exactly s at
    every position" is exact for the expected estimate. For a single realization it is
    tautological (it is exactly s times the corrected LMS). Suggested: "its expectation is
    exactly s = SBR/(SBR+1) times that of the background-free LMS at every position (the
    constant term cancels because Σb_k = 0)." This is also the definition of the
    lms_sbr10_shrink_factor key, which is not cited in the text.

Items checked and correct:
- intro statements;
- model Eqs. 1–5;
- the 360.3 nm ring diameter;
- vectorial zero and handedness physics;
- crb_center numbers and limits;
- LMS = S27 at the centre;
- MLE superefficiency "not a small-N effect" (Fig. 5c: 0.84/0.83/0.83/0.83/0.83);
- camera section;
- iterative section (the slope-mechanism speculation is labelled as such);
- the direction statements in nonidealities: loss larger at small L, first-order agreement O(x²),
  L_opt ∝ fwhm√ε independent of fwhm (exact by scale invariance), off-centre minimum with
  4.887 nm as a scale;
- discussion items 1, 3, 5 and 6;
- appendix.

Figure checks:
- Fig. 1–4, 6, 7 and 8 match their captions.
- In Fig. 4(d) the 13-photon point-value marker is not visible or labelled (cosmetic, known from
  r3).
- The Fig. 8 caption (tex) and figures.json agree.

---

## M2 — Independent verification of the round-4 keys

**Population misalignment (vpop.py).**

Method:
- Own LG and TCP (rotation π/2), with Eq. S30 and SBR = 10 in both the true and the nominal
  model.
- All four zeros are displaced by δ in independent uniform random directions.
- Naive MLE = argmax Σ q_i ln p_nom,i(r), with q = p_true(r0) (expected counts). I solve it by
  vectorized Fisher scoring started at r0.
- A brute-force check on a 0.5 nm grid over the 0.75L disk found the same maximum for 200
  patterns per case (0 mismatches).
- Seeds 101, 202 and 303 are independent of the worker's seed 42, with 20 000 patterns each.

| | δ=2 | δ=5 | δ=10 | LSQ slope (3 seeds, ±0.0022 / ±0.0032 each) | registry |
|---|---|---|---|---|---|
| centre | 0.7448–0.7478 | 0.7534–0.7563 | 0.7839–0.7867 | 0.7796 / 0.7768 / 0.7770 | 0.7768 ± 0.0050; ratios 0.7471/0.7549/0.7835 |
| (L/4,0) | 0.8035–0.8072 | 0.8093–0.8141 | 0.8614–0.8686 | 0.8552 / 0.8561 / 0.8495 | 0.8507 ± 0.0074; ratios 0.8023/0.8073/0.8635 |

Pooled over my 60 000 patterns the slopes are 0.7778 ± 0.0013 (centre) and 0.8536 ± 0.0018
(L/4,0). Both registry values agree within 0.2–0.4 of their SE.

Independent δ→0 check: the exact linear-response map gives |bias|/δ = 0.7471 ± 0.0005 (centre)
and 0.8038 ± 0.0006 (L/4,0). This matches the δ = 2 ratios and confirms that |bias|/δ grows
with δ (nonlinear).

Convergence (was 4000 patterns enough?):
- At (L/4,0) the per-pattern slope is right-skewed: median 0.805, p99 2.0, max 4.4.
- Running means over 400 patterns therefore fluctuate by about ±0.02 (1σ). My running means were
  0.853 / 0.853 / 0.863 / 0.855 at 400 / 1000 / 4000 / 20 000 patterns.
- The worker's sequence 0.896 / 0.876 / 0.851 is consistent with this: about 2σ at 400 and 1000,
  on nested samples.
- 4000 patterns give ±0.0074 (0.9 %). The quoted SE (std of per-pattern slopes / √P) is
  meaningful: my between-seed spread at 20 000 patterns (0.8495–0.8561) matches ±0.0032.

So 4000 is converged at the stated SE. The bias is bounded, so the variance is finite despite the
tail.

Also resolved: the writer's concern that "0.777 differs from MC 0.755 by 2.4 SE". The MC weights
applied to the population ratios give 0.762; the rest is 400-pattern sampling. It is not a
contradiction, but see issue 9 for the wording.

Not reproduced: the worker's claim that "with the same first 400 patterns ... matches the MC
within 1 MC SE" needs their seed-42 pattern stream, which I did not use. It stays unclear.

**Demo keys (vdemo.py; own finite-difference Fisher at the centre, L = 100, N = 100, fwhm = 300):**
- eps_constant_pedestal_sbr_eff_L100_eps0p05: 2.907508. I get the same value from the closed form
  and from ΣI_LG(0)/(4ε).
- eps_constant_pedestal_over_s31_L100_eps0p05: numeric CRB 4.980611 / S31(SBR_ε) 4.980611 =
  1.0000000.
- eps_gaussian_pedestal_over_s31_L100_eps0p05: 5.007113 / 4.980611 = 1.0053210.
- Both pedestal CRBs are continuous at r = 10⁻³ and 10⁻² nm. The Table II entries 1.30023 and
  1.30715 are reproduced.
- The constant pedestal plus SBR = 10 ("beam" convention): numeric 5.428956 = S31(SBR_eff)
  5.428956.
- lms_sbr10_shrink_factor: with my own LMS (Eq. S50, rotation π/2) the ratio of uncorrected to
  corrected is 0.9090909 = 10/11. The ratio of the uncorrected LMS to the background-free LMS on
  its own expected counts is also exactly 10/11 at (5,0), (4,−2) and (20,7).

Note: these four demo keys are not cited in the manuscript. They are registry-only.

---

## M3 — Literature attributions

I checked every \cite against notes A and B. I spot-checked the PDFs with PyMuPDF (read-only):
- Balzarotti p. 6: "<0.2 % of the doughnut crest".
- Balzarotti p. 13: Eq. S13 "arithmetic mean of the eigenvalues σ_i²" = √(tr/d).
- Balzarotti pp. 35 and 44: β0 = 1.27, β1 = 3.8 for the live mLMS.
- Masullo SI p. 9: Supplementary Table 1, MINFLUX 0.94/1.96/3.16 in column 2 ("precision at the
  center"; its legend gives no SBR or N).
- Stefani 2023: "for parabolic minima p(x) is scale invariant ... key to zooming-in".
- Lopez 2023: circular polarization "along the same direction of the phase ramp".

All attributions are correct:
- Balzarotti: Eqs. S13, S17, S22c, S24, S27, S29–S31, S49–S51 and S61; σ_PSF = 100 nm; iterative
  scheme only qualitative; MLE needs larger N off-centre.
- Gwosch: 3D.
- Masullo2022 and rastmin: framework.
- MasulloStefani2022: λ ∝ I^c.
- Stefani2023: scale invariance and tracking.
- Lopez2023 and Caprile2022: handedness, Richards–Wolf, Fig. 10.

No invented claims about others' work. The three wording nuances (handedness formalization, S24
rotation, "Most analyses ... [1]") and the inferred Masullo conditions are listed in issues 17
and 19.

---

## What stays open
- The worker's "same first 400 patterns vs MC within 1 SE" (needs their stream; low priority).
- The earlier unclear items (rmse vs Masullo Eq. 4.2; SBR = 10 confirmation by the author) are
  correctly declared in Sec. X.

```claims
[{"status": "verified", "text": "M2-pop: misalignment_naive_bias_over_delta_pop_center = 0.7768 +- 0.0050 and _pop_Lq = 0.8507 +- 0.0074 (4000 patterns) reproduced independently (own LG/TCP/Eq. S30 model, naive MLE on expected counts by Fisher scoring, grid-checked): three independent 20000-pattern seeds give 0.7796/0.7768/0.7770 (+-0.0022) and 0.8552/0.8561/0.8495 (+-0.0032); pooled 0.7778+-0.0013 and 0.8536+-0.0018"},
 {"status": "verified", "text": "M2-pop-ratios: per-delta population ratios |bias|/delta = 0.7471/0.7549/0.7835 (centre) and 0.8023/0.8073/0.8635 (L/4,0) for delta=2/5/10 reproduced (own 0.745-0.748/0.753-0.756/0.784-0.787 and 0.804-0.807/0.809-0.814/0.861-0.869); exact linear-response limit 0.7471+-0.0005 (centre) and 0.8038+-0.0006 (L/4,0); |bias|/delta grows with delta"},
 {"status": "verified", "text": "M2-pop-convergence: at (L/4,0) the per-pattern slope is right-skewed (median 0.805, p99 2.0, max 4.4), so 400-pattern running means fluctuate by ~+-0.02; the worker's 0.896/0.876/0.851 (400/1000/4000) is consistent with sampling; 4000 patterns are converged at the quoted SE (+-0.9 %), and the SE std(per-pattern slope)/sqrt(P) is meaningful (between-seed spread at 20000 matches +-0.0032)"},
 {"status": "verified", "text": "M2-pop-vs-MC: the centre MC weighted slope 0.7546+-0.0093 vs population 0.7768 is explained by inverse-variance weighting (MC weights applied to the population ratios give 0.762) plus 400-pattern sampling; not a contradiction"},
 {"status": "unclear", "text": "M2-worker-400: 'with the same first 400 patterns as the 400x200 MC the noise-free naive |bias| matches the MC within 1 MC SE' - not reproduced (requires the worker's seed-42 pattern stream)"},
 {"status": "verified", "text": "M2-demo: eps_constant_pedestal_sbr_eff_L100_eps0p05 = 2.907508; eps_constant_pedestal_over_s31_L100_eps0p05 = 1.0000000 (numeric CRB 4.980611 = S31(SBR_eps)); eps_gaussian_pedestal_over_s31_L100_eps0p05 = 1.0053210 (5.007113/4.980611); constant pedestal + SBR=10 (beam convention) numeric 5.428956 = S31(SBR_eff); own finite-difference Fisher"},
 {"status": "verified", "text": "M2-lms: lms_sbr10_shrink_factor = 0.9090909 = 10/11; the expected uncorrected LMS equals s times the background-free LMS at (5,0), (4,-2), (20,7) (own Eq. S50 implementation)"},
 {"status": "verified", "text": "M1-appendix: master formula A2, S27 point value, central term 16e^x/(3L^2), alpha/beta, sigma_lim^2 = L^2/(8Ng^2)(3g^2+e^x)/(3g^2+2e^x), rho=2/sqrt5 and L^2/(10N) for x->0, rho ~ (2/sqrt5)(1-9x/40), ellipse ratio sqrt(3/5), S31 via b=S/(4SBR), r_c, SBR_eps, SBR_eff combination, gaussian-pedestal F ratio, first-order ratio 1+3x^2/7, multiphoton F=2Nc^2g^2/R^2 - all re-derived by hand and correct"},
 {"status": "verified", "text": "M1-conventions: manuscript uses fwhm = Eq. S17 size parameter (ring diameter fwhm/sqrt(ln2)=360.3 nm), L = diameter, SBR total/total (camera SBR_c per pixel, 500/81=6.173), CRB = sqrt(tr F^-1/2) = Eq. S13, 'centre CRB' = r->0 limit unless stated; all backlog caveats (L exponent range, 2L search radius, honest-at-floor only at centre, N=250 bootstrap, adaptive 97 %, slope excludes -0.5, Caprile F=1, 0.78 only eps<~0.01, SBR=10 convention, superefficiency reported) are present; 0.862 and 0.79 not cited"},
 {"status": "verified", "text": "M1-pdf: paper/main.pdf (13 pp.) is newer than all sections, numbers.tex and paper_numbers.json; check_provenance 'all checks pass' (178/178); figures 1-4, 6-8 match their captions; no private/unpublished data in the paper"},
 {"status": "refuted", "text": "M1-vectorial.tex:32-33 'wrong hand CRB two orders of magnitude larger, linear more than one order' - Table I gives factors 92/27/16 (opposite) and 38/11.5/6.7 (linear) for L=50/100/150; fix: 'by factors 16-90 and 7-40 (L=150-50 nm)'"},
 {"status": "refuted", "text": "M1-fig5-caption (estimators.tex:66 and figures.json fig5 d) 'bias ... toward the centre' - panel (d) is negative only for x0<~8 nm and +0.18 nm (away from the centre) at x0=12 nm; restrict to x0<~8 nm and note the sign change"},
 {"status": "refuted", "text": "M1-estimators.tex:35-38 'the linearized estimators ... their sigma falls below the CRB' - in Fig. 5(b) the mLMS has sigma/CRB ~1.3 >1 for x0<~12 nm and a small outward bias at 5-15 nm; the plotted LMS/mLMS include the 1/s correction, so 'shrinkage' must not read as the 1/s shrinkage; restrict to x0>~15 nm and say 'linearization compression'"},
 {"status": "refuted", "text": "M1-methods.tex:9-11 'enough repetitions for a statistical error below 2 %' except misalignment - mle_nobg_bias_x_r2_nm (2.2 %) and adaptive_frac_outside_next_radius (6 %) exceed 2 % and are not misalignment keys (they are correctly quoted with \\pnumse); reword"},
 {"status": "refuted", "text": "M1-nonidealities.tex:97-101 the MC-vs-population slope difference attributed to weighting 'and the MC noise floor' - the floor would raise the MC slope, which is lower (0.755 vs 0.777); the difference is weighting (0.762 with MC weights) plus 400-pattern sampling"},
 {"status": "unclear", "text": "M1-wording (not wrong numbers, but overstated or ambiguous; fixes in r05-verifier.md): vectorial:26-28 'few percent at the largest L' (Fig. 2f reaches ~11 % at 250 nm); vectorial:28-31 curvature justification (overall curvature cancels); crb_center:20 'isotropic'; model:50-52 SBR held fixed at every position not stated; nonidealities:103-106 honest centre bias 9 % (3 SE) above floor called 'at the floor'; abstract:23 and discussion:21 'about 0.7768' (ratio 0.747-0.784, LSQ 78 % weighted by delta=10; conditions missing); abstract (iii)/(iv) conditions missing; discussion:56-59 'all absent'; open_points:28-30 now verified; open_points:18-21 overstated (agree ~2 % except F=1); crb_center:67-69, 99-101; SE formatting with 4 significant figures"},
 {"status": "verified", "text": "M3: every \\cite checked against notes A/B and spot-checked in the PDFs (Balzarotti p.6 <0.2 %, p.13 Eq. S13, pp.35/44 beta=1.27/3.8; Masullo SI Table 1 0.94/1.96/3.16; Stefani 2023 scale invariance; Lopez 2023 'same direction of the phase ramp'); no misattribution or invented claim; minor wording: handedness formalization is the project's own, TCP is S24 up to rotation/relabelling, 'Most analyses ... [1]' overgeneralizes (Ref. [1] experiments used a measured PSF), Masullo column's SBR=5/N=500 are inferred"}]
```
