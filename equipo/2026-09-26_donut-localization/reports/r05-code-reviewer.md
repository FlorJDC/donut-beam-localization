# Round 5: code reviewer (final delivery review)

Scope: D1-D5 of the r5 task, checked against the ENTREGABLE list in OBJECTIVE.md. I edited no project
file. I did not read `docs/private/`. What I produced is in
`work/review/r05/`: `png/` (renders of the 8 figures, a few zooms and the 13 pages of main.pdf),
`quick.log`, `pre_quick.sha`, `pn.txt` (registry dump) and `render.py`. After the run I deleted the
throwaway venv and the fresh clone: both sat untracked inside the parent repo, where a `git add -A`
would have picked them up.

Ledger read: intent.md, OBJECTIVE.md, inbox.jsonl (7 directions), state.json and reports r04-pi.
state.json is still at `round: 3`. It has no r4 history, and `last_acceptance` still says "vacuous".
Its local copy differs from HEAD (pending -> claimed:pending, not committed). See F6.

## D1. Fresh-clone reproducibility

- `git clone https://github.com/FlorJDC/donut-beam-localization`: HEAD is c90a402, the same as the
  local HEAD. The clone has 891 files and 21 MB.
- The clone has no `docs/private/`, no `data/mc/`, no `data/quick/` and no paper PDFs. `git ls-files`
  lists only these `*.pdf`: the 8 figures, `paper/main.pdf`, and 2 of our own draft compilations
  (`equipo/.../work/review/r03/tex/{main,revtex_main}.pdf`). No `*.npz` is tracked.
- `sha256 tests/test_acceptance.py` = `7d198853…2303`. This equals the guard in state.json, so the
  `.gitattributes` eol=lf setting works on a Windows clone.
- README Quick start, followed literally on Windows with Python 3.8.6:
  - `pip install -e .` **FAILS** (F1). The failure happens in a fresh venv, whose pip is the 20.2.1
    bundled with Python 3.8.6. The error is "editable mode currently requires a setup.py based
    build". It works after `python -m pip install --upgrade pip` (pip 25.0.1). Installing is not
    actually needed, because every test and `scripts/_paperconfig.py` puts `src/` on `sys.path`.
    System Python without donutloc installed runs the suite with OK.
  - `python -m unittest discover -s tests`: **218 tests OK (1 skipped: Tectonic not available)**,
    41 s.
  - `python -m unittest tests.test_acceptance -v`: **8/8 OK**.
  - `python scripts/check_provenance.py`: "178 tag(s), 178 entry(ies), 0 error(s)", then "claims: 39,
    numbers: 217", "all checks pass", exit 0.
  - `python scripts/reproduce.py --quick --no-latex`: **all 5 steps OK in 162 s**. I compared the
    sha256 of every tracked file before and after (`sha256sum -c`): identical. The quick products
    went only to `paper/figures/quick/`, `data/quick/` and `paper/generated/quick/`, which are
    git-ignored. The r4 D2 fix is confirmed on a fresh clone.
  - Full regeneration of cheap figures:
    - `fig_7_zero_depth.py --no-cache` took 5 s and `fig_1_schematic.py` took 6 s.
    - The summaries `data/fig1_summary.json` and `data/fig7_summary.json` are identical in content.
      `git diff --ignore-cr-at-eol` is empty, but git reports them as modified because they were
      written with CRLF.
    - The PDFs render pixel-identical to the committed ones. They differ only in the
      `CreationDate` metadata (F11).
- The `pip install -e .` run leaves `src/donutloc.egg-info/` untracked, because it is not in
  .gitignore (F7).

## D2. README accuracy

- Every path in the structure block exists. Every command runs.
- The figure table matches `structure/figures.json`: 8 ids, scripts and PDFs, questions R1-R5. Its
  numbers match `data/paper_numbers.json`:
  - zero depth 0 / 0.845 / 0.372 matches vectorial_zero_depth_*;
  - "L_opt ~ 0.78 fwhm sqrt(eps)" matches L_opt_over_fwhm_sqrt_eps 0.781/0.776.
- The claim "shared keys identical to the registry" holds: 149 shared keys, 0 mismatches.
- "Unit tests (about a minute)" is right (41 s). The "about 25 minutes" for a full run is consistent
  with r4: 17m48s of figures plus 376 s of numbers. I did not rerun the full pipeline.
- Problems found:
  - README.md:19: the `pip install -e .` failure (F1).
  - README.md:22: "fast smoke run". Steps 4-5 of `--quick` check the final products, not the quick
    ones (F17).
  - README.md:40-73: the structure block omits `paper/main.pdf` (the compiled paper, the thing a
    reader wants first), `OBJECTIVE.md`, `OBJECTIVE.template.md`, `environment.yml`,
    `job.cmd`/`job.sh`, `jobs/` and `scripts/apply-to-existing.*` (F9, F10).

## D3. Captions compared with the figures (rendered at 130-450 dpi)

Every figure is `figure*` at `\textwidth`, with PDFs 17.67 cm wide, so no scaling happens. main.pdf
has 13 pages and was built at 11:58, after the final figures (11:29-11:46). There is no `??` in the
text. Labels are readable at the printed size. Problems:

- **F2 (refuted). Fig. 4(d): a triangle is missing.** Both captions say triangles mark the photons
  needed for 5 nm, including "13.0 / 13 with the point value". The two captions are
  `structure/figures.json` fig4 and `paper/sections/crb_center.tex:135ff`, whose rendered text on
  main.pdf p. 7 reads "…MINFLUX 10.31 with the r→0 limit, 13 with the point value and 15.37 with
  SBR = 10". The figure has no S27 triangle: the loop in `scripts/fig_4_scaling.py:174-175` plots
  only `cam_ideal`, `cam_pixbg`, `mf_lim` and `mf_sbr10`. I confirmed this in the zoom
  `png/fig4d_zoom.png`.
  - Fix: add `("mf_s27", S.COLORS["lms"])` to that loop, with a text label. Or remove "13 with the
    point value" from both captions.
  - Also in (d): the dashed "camera 9x9 px, no bkg." line cannot be told apart from the solid ideal
    line (5.20 vs 5.00 nm at N=400). Mention that in the caption, or thin or offset the line.
- **F3 (refuted). Fig. 8 caption over-generalizes.**
  - `structure/figures.json` fig8 says "(b) … the naive sigma changes little, so its loss of accuracy
    is bias". `paper/sections/nonidealities.tex:159` says "the loss of the naive estimator is bias".
  - At (L/4,0) the naive σ grows from 2.389 to 3.74 ± 0.16 nm at δ=10 (+57 %) and to 5.02 nm at
    δ=15. The keys are misalignment_naive_sigma_Lq_d0/d10/d15_nm, and the growth is plainly visible
    in panel (b).
  - The body text (nonidealities.tex:112-120) states this correctly.
  - Fix: "at the centre the naive σ changes little and the loss is bias; at (L/4,0) σ also grows".
- **F4 (refuted). Fig. 5(d) caption: "bias toward the centre" is unqualified.** It appears in
  figures.json fig5 and `estimators.tex:66`.
  - The plotted bias changes sign near x0 ≈ 9 nm. It is +0.18 nm at x0 = 12 nm, pointing away from
    the centre.
  - Fix: "toward the centre for x0 ≲ 8 nm (−0.344 ± 0.008 nm at r=(2,0)), changing sign beyond".
    The same qualifier should go on estimators.tex:22 ("the price is a bias toward the centre").
- F13 (minor). figures.json fig2 (f) says the differences are "with respect to the LG donut with
  fwhm = 300 nm". The grey curve (vectorial / LG 327.1) is not relative to LG 300. The tex caption
  (vectorial.tex:80-84) has it right, so copy it into figures.json.
- F14 (minor). Fig. 3(f) has a dashed orange closed-form curve for fwhm = 360 nm. Neither caption
  mentions it (figures.json fig3, crb_center.tex:127). Add "dashed: fwhm = 360 nm".
- F15 (cosmetic):
  - Fig. 6(a): the grey annotation "no re-centring (search-disk artefact)" sits right under the
    legend and reads as a legend entry. The "10000" x tick label is clipped at the right edge. The
    blue N=1000 point is hidden under the square and diamond markers.
  - Fig. 8(c): the in-panel text reads "SBR  =10" (spacing).
  - Fig. 3(a): the contour labels "2" and "3" are tiny and overlap the zero markers.
- Checked and consistent (numbers in both captions against the registry and the figure):
  - Fig. 1: 327.14, 180.2/192.3, 384.67, <1 %.
  - Fig. 2: 0 / 0.845 / 0.372, <0.4 %.
  - Fig. 3: 1.605 / 1.802 / 1.960, 2/√5.
  - Fig. 4: 1.004, 360 nm, 400, 591.
  - Fig. 5: bump ≈1.5, 0.84, 0.75, 1.123, −0.34.
  - Fig. 6: slope −0.515, ratios 0.183 → 0.158, grey artefact point labelled.
  - Fig. 7: 0.78 only for ε ≲ 0.01, 0.70 at 0.15, 4.9 nm.
  - Fig. 8: 400×200, 0.76δ line, population 0.777 / 0.851.

Typography (F12): `\pnumse` prints 4 significant digits. The abstract of main.pdf p. 1 reads
"σ/σ_CRB = 0.9912 ± 0.005131", and the same happens in the section text. Fix: have
`compute_paper_numbers.py` format the SE to 1-2 significant digits, and the value to match.

## D4. Bibliography (`paper/references.bib`, 9 entries, all cited)

I checked each entry against the first pages of the PDFs (PyMuPDF, read-only):

| Entry | Result |
|---|---|
| Balzarotti2017 | Title and all 8 authors (Gynnå included) match. Science 355, 606-612 (2017) is correct. No DOI is printed in the PDF, and the bib gives none. |
| Gwosch2020 | DOI 10.1038/s41592-019-0688-0 matches the PDF. Nat. Methods 17, 217-224 is correct. |
| Masullo2022 | Authors and title match the SI's first page ("Biophysical Reports, Volume 2"). No article number or DOI. |
| Masullo2022rastmin | Authors (6), title, LSA 11:199 and DOI 10.1038/s41377-022-00896-4 all match. |
| MasulloStefani2022 | LSA 11:70 and DOI 10.1038/s41377-022-00763-2 match. |
| Stefani2023 | Nat. Photonics and DOI 10.1038/s41566-023-01239-4 match. No volume or pages. |
| Lopez2023 | All 7 authors match the accepted manuscript. JOSA B and DOI 10.1364/JOSAB.482413 match. No volume or pages (the PDF is the accepted version). |
| Caprile2022 | Authors and title match. No volume, article number or DOI (none printed). |
| agentteam | Software, URL. |

**No invented or wrong field was found.** Four entries are incomplete (F16):
- Masullo2022, Caprile2022, Stefani2023 and Lopez2023 lack volume and pages.
- Masullo2022 and Caprile2022 also lack a DOI.
- Suggested fix: complete them from the publisher pages. From memory, to be confirmed, not taken from
  the PDFs: Caprile2022 is CPC 275, 108315, and Masullo2022 is Biophys. Rep. 2, 100015.

`papers/README.md` agrees with the bib. For Masullo2021pMINFLUX, the PDF prints DOI
10.1021/acs.nanolett.0c04600; the "21, 840" cannot be checked from this ASAP PDF, but it is not cited
in the paper. The table is in Spanish while the rest of the public docs are in English (cosmetic).
Two hand-typed literature values in the tex are verified against Balzarotti 2017:
- "zero below 0.2 %" (nonidealities.tex:5): "doughnut minimum amounted to <0.2% of the doughnut
  crest";
- mLMS coefficients (1.27, 3.8) (model.tex:81): SI pp. 35/44.

## D5. Repository hygiene

- LICENSE: MIT, 2026, flor-choque. OK.
- CITATION.cff:
  - The YAML is well formed by inspection (PyYAML is not installed). It has cff-version 1.2.0,
    message, title, type, authors, repository-code, url, license, version and keywords.
  - The person is given only by `alias` and `website`, which the 1.2.0 schema allows.
  - Optional additions: `date-released`, `preferred-citation`.
- .gitignore covers `docs/private/`, `data/mc/*.npz`, all three quick directories, LaTeX aux files,
  `__pycache__` and `jobs/*`. It is **missing** `*.egg-info/`, `build/`, `dist/`, `venv/`, `.venv/`
  and `.pytest_cache/` (F7). It also does not cover review scratch such as clones or venvs under
  `equipo/*/work/`.
- Large files: none over 1 MB. The largest is main.pdf at 699 kB. `equipo/` holds 732 tracked files
  (2.2 MB), 611 of them under `work/verify/r03/scratch/` (F8):
  - 16 fixture copies of the whole pipeline: scripts, tex and provenance JSON;
  - about 30 PNGs in `work/figcheck/` and `work/papercheck/`.

  This is legitimate audit trail but heavy. The audit trail is intentionally versioned, so I
  suggest keeping the reports, state and inbox and pruning `scratch/` and the PNGs (or ignoring
  `equipo/*/work/**/scratch/`).
- Stray files (F9): `scripts/apply-to-existing.ps1` and `scripts/apply-to-existing.sh` are
  agent-team template installers. They have nothing to do with the paper pipeline and the README
  does not mention them. Move them under `agent-team/` or delete them.
- Private data:
  - `git grep` for p-minflux, private, unpublished, "no publicad", C_pminflux and "several %" finds
    only the policy statements (CLAUDE.md, OBJECTIVE.md, .gitignore) and earlier reviewers' notes.
  - Nothing from the author's unpublished work shows up in paper/, docs/, src/, scripts/, data/ or
    structure/.
  - `inbox.jsonl:2` still carries the generic phrase "en donas experimentales el mínimo residual
    puede ser de varios % del pico", already flagged in r01 and not specific data.
- structure/claims.json compared with state.json (F5, F6):
  - claims.json has 37 verified, 1 superseded (not cited) and 1 `pending-r4-verification`. The
    pending one is `misalignment_naive_bias_population`.
  - That pending claim's key `misalignment_naive_bias_over_delta_pop_center` (0.7768) **is cited in
    the abstract** and in the fig8 caption. OBJECTIVE ("Todo resultado que entra al manuscrito debe
    haber sido reproducido por el verificador") requires verification before delivery. It depends
    on the r5 verifier. Until then it is unclear.
  - The 4-digit precision is spurious for an estimate with SE ≈ 0.005.
  - check_provenance and acceptance do not look at claim status, so nothing machine-enforces this.
  - state.json must be advanced to r4/r5 (history, statuses of the 12 `claimed:pending` r4 claims,
    last_acceptance) before freezing.
- Other:
  - F18: tests/test_paper_tooling.py:268 opens files without closing them, which floods the test
    output with ResourceWarnings.
  - F19: scripts/reproduce.sh prefers `python3`, which on Windows Git Bash can be the Microsoft
    Store stub.

## r4 tooling claims I could check on the fresh clone

- D2 (the `--quick` isolation): verified. Tracked sha256 values were unchanged after a real
  `reproduce.py --quick`.
- D4: the `numbers.tex` fallback uses `\detokenize{#1}` (numbers.tex:3-4). The Tectonic compile
  test is skipped here because Tectonic is not available.
- D1 (cache hash) and D3 (check_provenance holes): their dedicated tests (TestCacheInvalidation, 28
  check_provenance tests) pass on the clone. I did not re-attack them adversarially this round.
- README.md:88 correctly says code changes are not detected by the cache.

## Summary of requested fixes (by priority)

1. F2: draw the S27 triangle in fig4(d), or drop "13 with the point value" from both captions.
2. F3: fig8 caption, "loss is bias" only at the centre.
3. F4: fig5(d) caption, "toward the centre" only for x0 ≲ 8 nm.
4. F5: have the population misalignment claim verified (r5 verifier) before delivery, or remove it
   from the abstract. Cite it with fewer digits.
5. F1: README Quick start. Add `python -m pip install --upgrade pip`, or say no install is needed, or
   add a `setup.py` shim.
6. F6: update state.json to r4/r5.
7. F7, F9-F12: .gitignore, the stray scripts, the README structure block, deterministic PDF/JSON
   writing (`savefig(..., metadata={"CreationDate": None})`, `open(..., "w", newline="\n")`), and
   SE digits.
8. F8 and F13-F19: cosmetic, or prune when convenient.

```claims
[{"status": "verified", "text": "Fresh clone of https://github.com/FlorJDC/donut-beam-localization (HEAD c90a402 = local HEAD) on Windows/Python 3.8.6: unittest discover 218 OK (1 skip, no Tectonic), tests.test_acceptance 8/8 OK, check_provenance exit 0 (178 tags, 39 claims, 217 numbers), sha256 tests/test_acceptance.py = 7d198853...2303 (guard matches)"},
 {"status": "verified", "text": "reproduce.py --quick --no-latex on the fresh clone: 5/5 steps OK in 162 s; sha256 of every tracked file unchanged; quick outputs only in git-ignored paper/figures/quick, data/quick, paper/generated/quick (r4 D2 fix confirmed)"},
 {"status": "verified", "text": "Full regeneration of fig1 and fig7 on the fresh clone reproduces the committed products: summaries identical in content (CRLF only), PDFs pixel-identical (only CreationDate differs)"},
 {"status": "verified", "text": "Clone hygiene: no docs/private, no data/mc or quick outputs, no *.npz, no third-party paper PDFs (only our 8 figures, main.pdf and 2 r03 draft compilations); no author-private content found by git grep"},
 {"status": "verified", "text": "README paths, commands, figure table and numbers match structure/figures.json and data/paper_numbers.json; 149 shared summary keys identical to the registry"},
 {"status": "verified", "text": "references.bib: 9 entries, all cited; authors/titles/years/journals/volumes/pages/DOIs present agree with the PDFs' first pages; no invented field"},
 {"status": "refuted", "text": "README Quick start 'pip install -e .' (README.md:19) fails with the pip 20.2.1 bundled with Python 3.8.6 ('editable mode currently requires a setup.py based build'); works after 'python -m pip install --upgrade pip'. Install is not needed because tests/scripts add src to sys.path. Fix the README or add a setup.py shim"},
 {"status": "refuted", "text": "Fig. 4(d) captions (structure/figures.json fig4; crb_center.tex caption, main.pdf p.7) say triangles mark MINFLUX 13 photons with the S27 point value, but scripts/fig_4_scaling.py:174-175 never plots mf_s27; no such triangle in the figure"},
 {"status": "refuted", "text": "Fig. 8 captions (figures.json fig8 'naive sigma changes little, so its loss of accuracy is bias'; nonidealities.tex:159 'the loss of the naive estimator is bias') are false at (L/4,0): naive sigma 2.389 -> 3.74+-0.16 nm at delta=10 (+57%), 5.02 nm at delta=15; restrict to the centre as the body text (nonidealities.tex:112-120) does"},
 {"status": "refuted", "text": "Fig. 5(d) captions (figures.json fig5; estimators.tex:66) say the background-free MLE bias points toward the centre, but the plotted bias changes sign near x0~9 nm (+0.18 nm at x0=12 nm); qualify to x0 <~ 8 nm"},
 {"status": "unclear", "text": "misalignment_naive_bias_over_delta_pop_center (0.7768) is cited in the abstract and fig8 caption while its claim in structure/claims.json is 'pending-r4-verification'; it needs r5 verifier confirmation before delivery (OBJECTIVE: only verified results in the manuscript); check_provenance/acceptance do not enforce claim status"},
 {"status": "unclear", "text": "state.json not advanced past round 3 (no r4 history, last_acceptance still 'vacuous', 12 r4 claims 'claimed:pending', local uncommitted edit); must be updated before freezing"},
 {"status": "unclear", "text": "Hygiene items (non-blocking): .gitignore lacks *.egg-info/, venv/, build/, dist/; scripts/apply-to-existing.{sh,ps1} are stray template installers; 611 scratch fixture files under equipo/.../work/verify/r03/scratch; README structure omits paper/main.pdf, OBJECTIVE*.md, job.*, environment.yml; \\pnumse prints 4 sig. digits (abstract '0.9912 +- 0.005131'); figure PDFs/JSON not byte-deterministic (CreationDate, CRLF)"},
 {"status": "unclear", "text": "Minor caption gaps: figures.json fig2(f) says differences are relative to LG300 though the grey curve is vectorial/LG327.1; fig3(f) dashed fwhm=360 curve not mentioned; fig6(a) annotation abuts legend and '10000' tick clipped; fig4(d) dashed 9x9 no-bkg line indistinguishable from the ideal line. Bib incomplete (no volume/pages for Masullo2022, Caprile2022, Stefani2023, Lopez2023; no DOI for Masullo2022, Caprile2022)"}]
```
