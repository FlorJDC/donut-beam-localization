# Round 5 — Final check (verifier, short)

I checked the final products themselves, not the writer's or worker's descriptions of their fixes.
Sources: paper/sections/*.tex, the compiled paper/main.pdf, the figure PDFs, structure/figures.json,
README.md, data/paper_numbers.json and paper/generated/numbers.tex. I did not read docs/private/, and
I edited nothing outside work/verify/r05/final/ and this report.

Artefacts in work/verify/r05/final/:
- p01–p13.png: all pages of main.pdf at 110 dpi.
- fig4/5/6/8 *.png: the figures at 150 dpi.
- main.txt: the full text of the PDF.
- chknum.py: my own check of every macro against the registry.

## Sync and checks

- main.pdf (16:19:06) is newer than every section, numbers.tex / paper_numbers.json (16:16:45) and
  every figure (fig4/fig6 at 16:07). The PDF therefore contains the reformatted numbers and the new
  fig4. I confirmed this visually on p. 7: the "13" triangle is present.
- `python scripts/check_provenance.py`: 179 tags, 179 entries, 0 errors, 0 warnings. Claims 39,
  numbers 217. "all checks pass".
- `python -m unittest tests.test_acceptance`: 8 tests OK.

## Numbers: formatting and agreement with the registry

chknum.py parses all 284 `\pnum` / `\pnumse` macros. For each one it checks two things:
- The printed value lies within half a unit of its last printed digit of the value in
  data/paper_numbers.json.
- For every key that has an SE, the value and the SE are printed to the same decimal, and the SE
  has 1–2 significant digits (from 1 to 24 units of the last digit).

Result: 0 real failures. The script flagged one entry, `vectorial_zero_curvature_nm2`, but that is
a bug in my parser (nested braces in `7.042\times10^{-5}`); the value is correct.

Examples as printed:
- 0.991 ± 0.005 (was 0.9912 ± 0.005131);
- 1.943 ± 0.010;
- 0.1606 ± 0.0008;
- −0.344 ± 0.008;
- 0.777 ± 0.005.

Numbers in the abstract, checked against the registry:

| Printed | Registry key | Registry value |
|---|---|---|
| 1.605 / 1.802 nm | | |
| 0.9962 | | |
| 0.8452 | | |
| 0.991 ± 0.005 | | |
| 0.508 nm | | |
| 0.1606 | | |
| 0.747–0.783 | pop_center_d2, pop_center_d10 | 0.74713, 0.78349 |

All are correctly rounded. Table I and Table II print 4 significant figures and keep trailing
zeros (1.600, 2.130, 1.000).

## PDF scan

- The text contains no "??".
- The only doubled words and stray " ," are artefacts of extracting the maths.
- None of the pages I viewed (p5, p7, p10) has overfull text. P5 §IV.B has a few loose, underfull
  justified lines (cosmetic).
- I found no broken sentence. "(Gaussian pedestal, = 100)" in the extracted text is an extraction
  artefact: the tex has `$N=100$`.

## Item by item

1. **vectorial.tex "two orders / more than one order" (V1): fixed.**
   - p. 3 now reads "by a factor that depends on L (Table I)".
   - It cites 1.600→147.9 nm (L = 50) and 5.331→84.20 nm (L = 150), and for linear polarization
     60.77 and 35.72 nm, all matching Table I.
   - "Relative loss largest at small L" is correct: the ratio is 92× vs 16×.
   - V2 (few percent) and V3 (curvature cancels, departure from quadratic at L/2) are also fixed.
2. **fig5(d) bias direction (V7/F4): fixed** in all three places:
   - the caption on p. 8;
   - figures.json fig5;
   - the body text on p. 5 ("toward the centre close to it [for x0 ≲ 8 nm]").

   In the rendered panel, the point at 8 nm is −0.05 and the point at 12 nm is +0.18, so the zero
   crossing is at about 8.8 nm. The text's "near 9 nm" is correct.
3. **Estimators LMS/mLMS below the CRB and shrinkage (V6/V20): fixed in the tex, with a small
   residual.**
   - The 1/s shrink is now stated for the expectation.
   - The text now says the plotted curves include 1/s, and that σ < CRB comes from compression.
   - The mLMS details are correct against data/mc/fig5_sweep.npz:
     - mLMS σ/CRB = 1.36, 1.36, 1.22 and 0.97 at x0 = 0, 5, 10 and 15 nm;
     - mLMS bias +1.65, +2.45 and +1.61 nm at x0 = 5, 10 and 15 nm.

   Residual: "Closer to the centre this does not hold" is not true for the LMS. At x0 = 10 nm the LMS
   already has σ/CRB = 1.828/2.268 = 0.81 and a bias of −1.75 nm. Suggested wording: "for
   x0 ≲ 12 nm the mLMS has σ above the CRB (the LMS drops below it from x0 ≈ 10 nm)". Minor.
4. **methods "<2 %" (V13): fixed.** p. 11 reads: "The central results have a statistical error below
   2%; the others are quoted with their SE". It then lists the misalignment quantities, the bias at
   r = (2,0) and the adaptive fraction. The adaptive fraction is now printed with ± SE on p. 8 and
   p. 10.
5. **MC vs population slope explanation (V9): fixed.** p. 9 now attributes the difference to
   inverse-variance weighting (favouring δ = 2 nm) plus the sampling of 400 patterns. The noise
   floor is no longer mentioned.
6. **Wording issues (V4, V5, V8, V10–V12, V14–V17, V19): all present in the PDF.**
   - "direction-independent … ellipse anisotropic";
   - α and β now point to Eq. (A4);
   - the SBR is held fixed at every position (§II.C and the scope bullet);
   - honest MLE: "slightly above [floor], by about three SE" ((0.1797−0.1643)/0.0051 = 3.0);
   - abstract and discussion give ranges with conditions;
   - "Every computed number";
   - the iterative conditions are stated;
   - "error 0.1577–0.1831 times", plus the SBR = 10 exception;
   - the population-slope open point is removed;
   - Caprile is now "within a few percent … differ clearly at filling factor 1";
   - r_c formula;
   - Masullo "up to rounding … our inference";
   - handedness, S24 rotation and the Ref. [1] wording.

   Residual (minor): figures.json fig8(a) still says "At the centre the honest MLE stays at the
   noise floor". The tex now says "close to … slightly above".
7. **F1 README pip install: fixed.** Quick start installs only `requirements.txt`, which exists. The
   editable install is described separately, with the pip-upgrade caveat.
8. **F2 fig4(d) missing triangle: fixed.** The fig4_scaling.pdf render shows the orange "13"
   triangle next to 10.3 and 15.4. The dashed-with-squares 9×9 no-background line is
   distinguishable. The caption on p. 7 matches ("13.00 with the point value"; "dashed, nearly
   coincident").
9. **F3 fig8 "loss is bias" at (L/4,0): fixed** in the tex caption (p. 11), the body heading (p. 9
   "At the centre the loss is bias") and figures.json. In the render, panel (b) shows the naive
   (L/4,0) σ growing to about 3.7 nm at δ = 10, matching the cited 2.389→3.74.
10. **fig6(a) cosmetics (F15): fixed.** The "10000" tick is inside the panel and the annotation no
    longer collides with the legend.
11. **NEW ERROR, README.md:125 (fig5 row).** It says "the linearized estimators fall below the CRB
    only beyond x0 ~ 15 nm". The LMS is already below it at x0 = 10 nm (σ/CRB = 0.81; 0.97 at
    5 nm). "Only beyond 15 nm" holds for the mLMS alone. Fix: "the mLMS falls below the CRB beyond
    x0 ≈ 15 nm and the LMS from ≈ 10 nm, by compression".

## Verdict

All 5 verifier refutations and all 4 code-reviewer refutations (F1–F4) are fixed in the final
products. Numbers are formatted to SE precision and equal the registry. The PDF is in sync.
Provenance and acceptance pass.

Three minor residuals remain, none affecting a number:
- the README fig5 row (new error);
- one clause about the LMS near the centre in estimators.tex;
- "stays at the noise floor" in figures.json fig8.

```claims
[{"status": "verified", "text": "vectorial.tex 'two orders / more than one order' (V1) fixed: text now gives L-dependent Table I values (1.600->147.9 nm at L=50, 5.331->84.20 at L=150; linear 60.77/35.72), consistent with Table I; V2/V3 wording also fixed"},
 {"status": "verified", "text": "fig5(d) bias direction (V7/F4) fixed in PDF caption, figures.json and body: negative for x0<~8 nm, -0.344+-0.008 at (2,0), sign change near 9 nm (plotted -0.05 at 8, +0.18 at 12)"},
 {"status": "verified", "text": "Estimators LMS/mLMS below-CRB/shrinkage (V6/V20) fixed in tex: 1/s shrink stated for the expectation, plotted curves include 1/s, sigma<CRB attributed to compression; mLMS sigma/CRB 1.36/1.36/1.22/0.97 at x0=0/5/10/15 and outward bias 1.6-2.5 nm at 5-15 nm confirmed from fig5_sweep.npz"},
 {"status": "refuted", "text": "estimators.tex 'Closer to the centre this does not hold' (residual, minor) — false for the LMS: at x0=10 nm LMS sigma/CRB=0.81 and bias -1.75 nm; restrict the clause to the mLMS"},
 {"status": "verified", "text": "methods '<2 %' (V13) fixed: 'central results below 2%; the others quoted with their SE' listing misalignment, bias at (2,0) and adaptive fraction (now printed with +-SE)"},
 {"status": "verified", "text": "MC-vs-population slope explanation (V9) fixed: now inverse-variance weighting favouring delta=2 plus 400-pattern sampling; noise floor removed"},
 {"status": "verified", "text": "Verifier wording items V4,V5,V8,V10,V11,V12,V14,V15,V16,V17,V19 all present in the compiled PDF (honest MLE 'about three SE above' = (0.1797-0.1643)/0.0051=3.0)"},
 {"status": "refuted", "text": "structure/figures.json fig8(a) 'At the centre the honest MLE stays at the noise floor' (residual, minor) — tex now says slightly above (3 SE); figures.json not aligned"},
 {"status": "verified", "text": "F1 README pip install fixed: Quick start installs only requirements.txt (exists); editable install documented separately with pip-upgrade caveat"},
 {"status": "verified", "text": "F2 fig4(d) fixed: S27 '13' triangle present in fig4_scaling.pdf and in main.pdf p.7; caption matches (13.00 with the point value); dashed 9x9 no-bkg line distinguishable (open squares)"},
 {"status": "verified", "text": "F3 fig8 'loss is bias' fixed: restricted to the centre in tex caption, body heading and figures.json; (L/4,0) naive sigma growth 2.389->3.74+-0.16 nm cited and visible in panel (b)"},
 {"status": "verified", "text": "fig6(a) cosmetics fixed: 10000 tick inside panel, annotation clear of legend"},
 {"status": "refuted", "text": "README.md:125 fig5 row 'linearized estimators fall below the CRB only beyond x0 ~ 15 nm' (new error introduced this round) — LMS sigma/CRB is already 0.81 at x0=10 nm (0.97 at 5); only the mLMS crosses near 15 nm"},
 {"status": "verified", "text": "All 284 \\pnum/\\pnumse macros in numbers.tex equal data/paper_numbers.json to printed precision; every value-with-SE is printed to the SE's decimal with a 1-2 significant-digit SE; abstract numbers (1.605/1.802, 0.9962, 0.8452, 0.991+-0.005, 0.508, 0.1606, 0.747-0.783) correctly rounded"},
 {"status": "verified", "text": "main.pdf (16:19) is newer than all sections, numbers.tex, paper_numbers.json and figures; PDF text has no '??', no broken sentence, no visible overfull lines on inspected pages"},
 {"status": "verified", "text": "python scripts/check_provenance.py: 179 tags/179 entries, 0 errors, all checks pass; python -m unittest tests.test_acceptance: 8 OK"}]
```
