# Round 5: worker (final fix pass on code, figures, structure, README, hygiene)

Inputs: r05-code-reviewer.md (F1-F19) and r05-verifier.md (M1-M3). I did not touch paper/main.tex,
paper/sections/*, references.bib, provenance.json, main.pdf or tests/test_acceptance.py, and I did not
read docs/private/. No git commit. Scratch outputs are in `work/r05-worker/` (numbers snapshots, run
log, `numbers_diff.txt`) and `work/figcheck/r05_*.png`.

## Changes per item

1. **README (F1, F9, F10, F17 + verifier numbers).**
   - Quick start now says that no install is needed: tests and scripts put `src/` on `sys.path`.
     `pip install -r requirements.txt` (or `conda env create -f environment.yml`) is enough.
   - New paragraph for importing donutloc from elsewhere: `python -m pip install --upgrade pip`, then
     `pip install -e .`. It explains that pip 20.2 (Python 3.8.6) fails because the project has only
     `pyproject.toml`.
   - The structure block now includes paper/main.pdf ("start here"), OBJECTIVE.md,
     OBJECTIVE.template.md, pyproject/requirements/environment.yml, LICENSE/CITATION.cff, the equipo/
     contents, and job.cmd/job.sh as agent-team launchers (jobs/ is git-ignored).
   - `--quick`: the README now says that its provenance and acceptance steps check the final products.
   - Determinism is documented, and so are the number-formatting rules.
   - Figure table:
     - fig2: wrong hand 16-90x and linear 7-40x (L = 150-50 nm);
     - fig5: linearized estimators fall below the CRB only beyond x0 ~ 15 nm, by compression;
     - fig8: 0.75-0.78 δ at the centre and ~0.80-0.86 δ at (L/4,0) for δ = 2-10 nm (L = 100 nm,
       SBR = 10).
   - The methodology section says that the repo was created from agent-team-template and that the
     template installers were removed.
2. **fig4(d) (F2).**
   - The S27 point-value triangle is added (lms colour), labelled "13": `n5["mf_s27"]` =
     photons_for_5nm_L50_sbrinf = 12.9956.
   - The dashed 9x9 no-background camera line now also has open square markers (markevery=2), so it
     can be told apart from the solid ideal line. It lies 4.09 % above that line: 5.2045 vs 5.000 nm
     at N = 400.
   - The triangles are slightly smaller with a white edge, so the three MINFLUX markers (10.3, 13,
     15.4) stay distinct. Checked on a 600 dpi zoom: `work/figcheck/r05_fig4d_zoom.png`.
3. **structure/figures.json (English).**
   - fig2(f): the grey curve is vectorial/LG 327.14, not relative to LG 300.
   - fig3(f): the dashed fwhm = 360 nm closed form is now mentioned.
   - fig4(d): the S27 triangle and the dashed-with-squares line are described.
   - fig5: the plotted LMS/mLMS include 1/s. The σ/CRB < 1 and the centre-ward bias beyond
     x0 ~ 15 nm come from linearization compression (mLMS σ/CRB ~ 1.3 for x0 ≲ 12 nm). In (d) the
     bias is negative only for x0 ≲ 8 nm, is −0.34 nm at (2,0), changes sign near 9 nm and is
     +0.18 nm at 12 nm.
   - fig8: "loss is bias" now applies only at the centre. At (L/4,0) the naive σ grows from 2.39 to
     3.74 nm at δ = 10. The slope is written as the range 0.75-0.78 / ~0.80-0.86 (LSQ 0.777 / 0.851),
     independently reproduced.
   - **structure/claims.json:**
     - misalignment_naive_bias_population: now `verified`, verified_round 5. The caveat gives the r5
       verifier numbers and the range and significant-figure guidance.
     - misalignment_loss_is_bias: the statement and caveat are restricted to the centre, and the
       (L/4,0) growth is added.
     - misalignment_naive_bias: the "honest at floor" caveat now gives the numbers (9 %/3 SE at the
       centre, 19 %/5 SE at (L/4,0)).
     - mle_superefficient_nobg: caveat on the sign change added.
     - lms_sbr_shrink: the statement now uses the expectation wording, and the caveat says that the
       Fig. 5 compression is not the 1/s shrink.
     - vectorial_crb_opposite_linear: caveat with the L-dependent factors.
     - Result: 38 verified, 1 superseded, 0 pending.
4. **fig6(a) (F15).**
   - The annotation now sits right of the artefact cross, and ylim goes to 1000, so the legend sits
     above it with no overlap.
   - xlim is (180, 1.5e4), so the "10000" tick is inside the panel. Checked:
     `work/figcheck/r05_fig6_iterative.png`.
5. **Number formatting (F12, verifier 18).**
   - `scripts/compute_paper_numbers.py` gains the documented functions `fmt_value`, `se_decimals` and
     `fmt_value_se` (a rule block in the comments). `_fmt` is kept as an alias.
   - A value with SE: the SE is rounded to 2 significant figures if its leading two digits are 10-24,
     otherwise to 1. The value is rounded to the same decimal, and both `\pnum` and `\pnumse` use that
     pair.
   - A value without SE: 4 significant figures with trailing zeros kept (`1.600`, `1.000`).
     Scientific notation below 1e-3 and from 1e5. An exactly integral float with 100 ≤ |v| < 1e5 (a
     set parameter such as fwhm = 300.0) prints as `300`.
   - Tests: `TestNumberFormatting` in tests/test_paper_tooling.py (5 tests).
6. **Determinism (F11).**
   - `_paperstyle.savefig` passes `metadata={"CreationDate": None, "ModDate": None}`.
   - `write_summary`, the paper_numbers.json writer and the iterative_sweep.json writers (in
     compute_paper_numbers and run_iterative_sweep) use `newline="\n"`, `sort_keys=True` and a final
     newline.
   - Test: `TestDeterministicOutputs`. I also checked by hand: regenerating fig4 and fig6 gives
     byte-identical PDFs and summary (sha256 OK), and two full compute_paper_numbers runs give a
     byte-identical data/paper_numbers.json.
7. **Hygiene.**
   - .gitignore: `*.egg-info/`, `build/`, `dist/`, `venv/`, `.venv/`, `.pytest_cache/` and
     `equipo/*/work/**/scratch/` are added (g_break.py regenerates that scratch). I confirmed them
     with `git check-ignore`.
   - Deleted scripts/apply-to-existing.ps1 and .sh.
   - Deleted `work/verify/r03/scratch/` (611 files). The verifier's scripts, json, logs and png/ are
     kept.
   - tests/test_paper_tooling.py: every bare `open()` is replaced by closing helpers. The full suite
     passes under `-W error::ResourceWarning`.
   - reproduce.sh: `$PYTHON` wins. Otherwise it uses `python3` only if `python3 -c "import sys"`
     works, and falls back to `python`. On this machine python3 exits 127, and `sh
     scripts/reproduce.sh --dry-run --no-latex` runs python.exe with all 5 steps OK.

## Regeneration and number changes

- `python scripts/fig_4_scaling.py` and `python scripts/fig_6_iterative.py` used the cached MC and
  data/iterative_sweep.json. Their shared summary keys (14 + 14) equal the registry. The summary JSON
  changed only by the final newline.
- `python scripts/compute_paper_numbers.py` ran twice (400 s):
  - **no value, SE or field changed.** There are 217 keys before and after (`work/r05-worker/numbers_diff.txt`).
  - The paper_numbers.json diff is only the key order inside each entry, from sort_keys.
  - 127 of 284 numbers.tex macros change text, as formatting only. Examples: mle_efficiency_center
    0.9912 ± 0.005131 → 0.991 ± 0.005; iterative_sigma_nm 0.5077 ± 0.002683 → 0.508 ± 0.003;
    misalignment_naive_bias_over_delta_pop_center 0.7768 ± 0.00497 → 0.777 ± 0.005;
    mle_nobg_bias_x_r2_nm −0.3445 ± 0.007702 → −0.344 ± 0.008.
  - For the writer, these values gain trailing zeros: camera_ideal_N400_nm 5 → 5.000,
    lms_nobg_sigma_over_s27_center 1 → 1.000, eps_constant_pedestal_over_s31 1 → 1.000,
    photons_for_5nm_L50_sbrinf 13 → 13.00, crb_center_sbr5_N500_fwhm300_L150_nm 3.37 → 3.370,
    L_opt_over_fwhm_sqrt_eps_eps0p01 0.776 → 0.7760. fwhm_nm and fwhm_masullo_nm stay 300 and 360.

## Test outputs (real, after all changes and the writer's concurrent edits)

- `python -W error::ResourceWarning -m unittest discover -s tests` → `Ran 224 tests in 37.702s  OK (skipped=1)`.
  The skip is the Tectonic compile.
- `python -m unittest tests.test_acceptance -v` → 8/8 ok.
- `python scripts/check_provenance.py` → "provenance: 179 tag(s), 179 entry(ies), 0 error(s),
  0 warning(s) / claims: 39, numbers: 217, provenance entries: 179 / all checks pass".
- sha256 tests/test_acceptance.py = 7d198853e87c7bd94bb8be07a855651146e1d849310281546cd8089dcf262303,
  which matches the guard.

## Not done / open

- The other figure PDFs (fig1-3, 5, 7, 8) were not regenerated. They still carry a CreationDate
  until the next `make_all_figures.py --no-cache`, which will then make them deterministic too.
- F16 (bib volumes/DOIs), F13-F15 in the tex captions and all tex wording belong to the writer. The
  cosmetic items fig8(c) "SBR  =10" and fig3(a) contour labels were outside this task and are
  untouched.
- CITATION.cff is unchanged: the optional `date-released` and `preferred-citation` were not
  requested.
