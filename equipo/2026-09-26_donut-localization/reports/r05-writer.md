# Round 5 — Writer (final fix pass on the manuscript text)

Inputs: r05-verifier.md (M1 issues 1–20, M3), r05-code-reviewer.md (F2–F4, F13, F14, F16).
I edited only paper/sections/*.tex and paper/provenance.json. main.pdf was recompiled.
references.bib is unchanged (see below). I did not touch numbers.tex, the scripts, figures.json or the README.
I did not read docs/private/.

## Issue → fix

| # | Issue (source) | Fix (file:line, after edit) |
|---|---|---|
| V1 | "two orders / more than one order" for a filled zero (refuted) | vectorial.tex:36-45. Now "by a factor that depends on L (Table I)", citing the correct, opposite and linear Table I keys at L=50 and L=150. Relative loss is largest at small L. No hand-typed ratios. |
| V2 | "few percent at the largest L" | vectorial.tex:26-29. Cites crb_vec_over_lg300_L150_N100 and adds "difference keeps growing at larger L [Fig. 2(e,f)]". |
| V3 | Curvature justification | vectorial.tex:29-34. States that overall curvature cancels in p_i and that the CRB is set by the departure from quadratic at L/2 (Eq. master). The curvature-matched LG *happens* to reproduce this departure. |
| V4 | crb_center:20 "isotropic" | crb_center.tex:20-21. Now "direction-independent … (the limiting error ellipse itself is anisotropic, see below)". α, β now point to new label eq:alphabeta (crb_center.tex:42; appendix.tex:47). |
| V5 | model:50-52 SBR fixed at every position | model.tex:53-58. States that SBR is held fixed at every position (Eq. S30, constant SBR), implying a position-dependent background. With a fixed background the SBR would vary with position (S29) and at the centre fall as L decreases (S32). Also added to scope item open_points.tex:8-9. NOTE: I did not write the verifier's "fall off-centre". For fixed background the total intensity grows off-centre (Σ\|r−b_i\|² = 4r² + 3R²), so that phrasing looked wrong. |
| V6 | estimators:35-38 (refuted) | estimators.tex:34-45. The expectation of the uncorrected LMS is s times the background-free LMS; the constant cancels because Σb_k=0 (also V20). The plotted LMS/mLMS include 1/s. Beyond x0≈15 nm: bias toward the centre by *linearization compression*, σ<CRB reflects compression. For x0≲12 nm the mLMS σ>CRB, and at 5–15 nm its bias is slightly outward. |
| V7 / F4 | fig5(d) caption, "toward the centre" (refuted) | estimators.tex:72-75. Negative for x0≲8 nm, −0.34 at (2,0) (key), sign change near x0≈9 nm. Also qualified the body at estimators.tex:22-23 ("toward the centre close to it, x0≲8 nm"). |
| V8 | nonidealities:103-106 honest "at the floor" | nonidealities.tex:105-112. "Close to the floor … but slightly above it, by about three SE". At (L/4,0) "further above". "Honest at the floor holds, approximately, only at the centre". |
| V9 | MC vs population slope = weighting + noise floor (refuted) | nonidealities.tex:98-102. Now inverse-variance weighting (favours δ=2 nm, where \|bias\|/δ is smallest) plus 400-pattern sampling. Noise floor removed. |
| V10 | abstract/discussion "about 0.7768" | abstract.tex:25-28: "mean bias of \pnum{…pop_center_d2}–\pnum{…pop_center_d10} times δ (δ=2–10 nm, L=100 nm, SBR=10)". discussion.tex:20-24: same, and "larger at (L/4,0)" with the pop_Lq_d2–d10 keys. The LSQ slope stays in nonidealities.tex with its SE. |
| V11 | abstract (iii)/(iv) conditions; "every number" | abstract.tex:4 "Every computed number". (iii) adds "L=50 nm, N=100". (iv) adds "our protocol (L reduced from 150 to 25 nm, no background)" and "an error … times that of an ideal camera". |
| V12 | discussion:56-59 | discussion.tex:58-64. Now "an error ratio_min–ratio_max times that of the ideal camera". Zero depth, misalignment and background (except one SBR=10 run, iterative_sbr10_sigma_nm ± SE) were absent. |
| V13 | methods:9-11 "<2 %" (refuted) | methods.tex:9-14. "Central results below 2 %. The others are quoted with their SE": misalignment, the MLE bias at (2,0) \src{mle_nobg_bias_x_r2_nm}, and the adaptive fraction \src{adaptive_frac_outside_next_radius}. The discussion now quotes the adaptive fraction with ± SE (discussion.tex:43). |
| V14 | open_points:28-30 population slopes pending | Bullet removed from open_points.tex. nonidealities.tex:88-89 adds "both were reproduced by an independent computation with other seeds" (r5 M2, verified). |
| V15 | open_points:18-21 Caprile | open_points.tex:19-22. "Agree to within a few percent at filling factors 0.5, 2 and uniform, but differ clearly at filling factor 1. We do not claim agreement there." The figures are qualitative because there is no key (see below). |
| V16 | crb_center:67-69 "sizeable fraction" | crb_center.tex:68-71. Gives the leading-order r_c = (√3L/4)e^{−x/2}/√SBR (App. A). It is "not small compared with L" at realistic SBR, and the leading-order estimate there is only indicative. |
| V17 | crb_center:99-101 Masullo | crb_center.tex:100-104. "Match, up to rounding in the last digit at L=150 nm, … Supplementary Table 1". The legend states no SBR/N, so SBR=5, N=500 and this fwhm are our inference. The lit_masullo_table provenance statement was updated to match. |
| V19a | handedness rule | model.tex:23-25. "Circular polarization rotating in the same sense as the phase ramp [7,8]. In our formulation, charge × spin = +1." |
| V19b | TCP vs S24 | model.tex:29-30. "Up to a global rotation and a relabelling of the exposures (centre last)". |
| V19c | "Most analyses … [1]" | introduction.tex:15-17. "Theoretical analyses of MINFLUX, starting with the Supplementary Material of Ref. [1], use …". |
| V20 | LMS shrink wording | See V6. The lms_sbr_shrink provenance statement was updated ("expectation"). |
| F2 | fig4(d) 13-photon triangle | Caption kept (the worker adds the marker). Also added "without background (dashed, nearly coincident with the continuum line)" at crb_center.tex:152-153. |
| F3 | fig8 caption "loss is bias" (refuted) | nonidealities.tex:166-171. "At the centre the naive σ changes little and its loss is bias. At (L/4,0) it grows from \pnum{…sigma_Lq_d0} to \pnum{…sigma_Lq_d10} nm at δ=10." The body heading is now "At the centre the loss is bias" (nonidealities.tex:114), and the body now also cites the d0 value. |
| F13 | fig2(f) grey curve | vectorial.tex:93-97. Vectorial and LG327 are shown relative to LG300; the grey curve is vectorial relative to the curvature-matched LG, *not* to LG300. |
| F14 | fig3(f) dashed curve | crb_center.tex:131-132. "Closed form for fwhm=360 nm dashed". |
| F16 | Bibliography | **No change.** Checked the four PDFs with PyMuPDF, first pages plus a regex over all pages. Caprile2022 is a Word preprint with no journal/volume/DOI. Lopez2023 is the accepted manuscript: DOI only, already in the bib. The Masullo2022 SI prints only "Biophysical Reports, Volume 2", already in the bib. Stefani2023 prints the DOI only (already present), and its metadata gives no volume or pages. I added nothing from memory. |

For figures.json: my caption meanings for fig2(f), fig3(f), fig4(d), fig5(d) and fig8(b) match the code-reviewer and verifier fixes. The worker should mirror them.

## Provenance
- New entry `misalignment_naive_sigma_Lq_d0_nm` (script, compute_paper_numbers.py).
- Statements updated: `lit_masullo_table` (inferred SBR/N) and `lms_sbr_shrink` (expectation).
- File rewritten with the same format (indent 1, sorted keys, LF).

## Missing keys (rephrased qualitatively instead)
- Opposite/correct and linear/correct CRB ratio keys (e.g. `crb_vec_opposite_over_correct_L*`). I cited the Table I values instead.
- CRB difference vectorial vs LG300 at L=250 nm (~11 %). I wrote "keeps growing at larger L".
- Caprile Fig. 10 comparison values (515.5/382.9/377.6/401.4 vs 508/392/380/428 nm). I wrote "a few percent / differ clearly".
- r_c/L at SBR = 5, 10, 20. I gave the formula only.
- misalignment_naive_sigma_Lq_d15_nm (5.02 nm) is not in the registry and is not cited.

## Checks
- Tectonic (from paper/) compiles with no errors: 13 pages, only underfull-hbox warnings. The PDF text has no "??" and no undefined references.
- Pages rendered to work/papercheck/r05w/p01–p13.png. I inspected p1, p4, p8, p11 and p12. Fig. 5(d) and Fig. 8 captions match the panels: the sign change sits between x0=8 and 12 nm, and the naive (L/4,0) σ growth is visible.
- `python scripts/check_provenance.py`: "provenance: 179 tag(s), 179 entry(ies), 0 error(s), 0 warning(s)", then **ERROR: paper/generated/numbers.tex is out of sync with data/paper_numbers.json**. Both files are unchanged since 11:52. The mismatch comes from the worker's uncommitted edit to scripts/compute_paper_numbers.py (SE formatting), which they have not yet regenerated. It clears when the worker reruns compute_paper_numbers.py. After that, main.pdf must be recompiled so it picks up the new number formatting.

## Still open
- The verifier's M2-worker-400 item (needs the worker's seed-42 stream) is not in the manuscript. Nothing to do.
