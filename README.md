# donut-beam-localization

A reproducible simulation study ("paper companion") of single-emitter localization with
donut-shaped excitation beams (MINFLUX and variants): beam models, Cramer-Rao bounds,
estimators, iterative MINFLUX and experimental non-idealities.

Every number in the manuscript is produced by a script, stored in `data/paper_numbers.json`,
and cited in the LaTeX source through a macro that resolves to that registry, with a provenance
entry saying how to reproduce it. A checker enforces this.

## Quick start

Requirements: Python >= 3.8 with numpy, scipy and matplotlib (numba is optional). Tested on
Windows with Python 3.8, numpy 1.24, scipy 1.10 and matplotlib 3.7.

No installation of the package is needed: the tests and every script put `src/` on `sys.path`
themselves, so numpy, scipy and matplotlib are all you need (`pip install -r requirements.txt`,
or `conda env create -f environment.yml`).

```sh
git clone https://github.com/FlorJDC/donut-beam-localization.git
cd donut-beam-localization
pip install -r requirements.txt           # dependencies only; donutloc itself is not installed

python -m unittest discover -s tests      # unit tests (about a minute)
python scripts/reproduce.py --quick       # fast smoke run: reduced Monte Carlo, separate outputs
python scripts/reproduce.py               # full reproduction (about 25 minutes on a laptop)
sh scripts/reproduce.sh                   # same, from a POSIX shell
```

To import `donutloc` from elsewhere, install it in editable mode. The project is configured by
`pyproject.toml` only (no `setup.py`), which old pip versions cannot install in editable mode:
the pip 20.2 bundled with Python 3.8.6 fails with "editable mode currently requires a setup.py
based build". Upgrade pip first:

```sh
python -m pip install --upgrade pip
pip install -e .
```

`scripts/reproduce.py` runs, in order: the unit tests, `make_all_figures.py --no-cache`,
`compute_paper_numbers.py`, `check_provenance.py` and the acceptance test. If
[Tectonic](https://tectonic-typesetting.github.io/) is on `PATH` (or given by `--tectonic PATH`
or the `TECTONIC` environment variable) it also compiles `paper/main.tex`. `--with-sweep` reruns
the iterative-MINFLUX sweep first; without it the committed `data/iterative_sweep.json` is
used.

`--quick` writes only to `paper/figures/quick/`, `data/quick/` and `paper/generated/quick/`
(git-ignored). The quick numbers file carries `"quick": true`, and `check_provenance.py` rejects
such a file in the final location. A quick run can never overwrite the paper's products. Its last
two steps (provenance check and acceptance test) still check the committed final products, not
the quick ones.

Regeneration is byte-deterministic: figure PDFs carry no creation or modification date, and the
JSON files are written with sorted keys and LF line endings, so rerunning a script on unchanged
inputs leaves `git status` clean.

## Short report (informe)

`informe/informe.pdf` (5 pages, Spanish) summarizes the project, reviews it against SimuFLUX
(Marin & Ries, Nat. Commun. 17:246, 2026) and adds two studies for non-iterative pulsed
interleaved MINFLUX (p-MINFLUX): fluorophore flickering (sequential vs interleaved excitation)
and lifetime cross-talk between TCSPC windows (`src/donutloc/pminflux.py`,
`scripts/fig_9_pminflux_timing.py`), and estimator bias from unmodelled background with a fixed
TCP (`src/donutloc/background.py`, `scripts/fig_10_background_bias.py`). Rebuild the PDF from
`informe/informe.html` with `node informe/build_pdf.js` (Playwright/Chromium), or print the HTML
from any browser. The team record of that work is in `equipo/2026-09-28_sintesis-pdf/`.

## Repository structure

```
src/donutloc/          the Python package (units: nm)
  beams.py               scalar beams: LG01 donut (Balzarotti Eq. S17), Gaussian, quadratic zero
  vectorial.py           Richards-Wolf vortex donut (handedness, linear polarization, zero depth)
  patterns.py            exposure patterns: TCP (3 donuts on a circle of diameter L + centre),
                         polygons, random misalignment
  photons.py             multinomial photon model with background (SBR)
  fisher.py              Fisher information and CRB (with the r -> 0 limit at a perfect zero)
  closed_forms.py        closed-form centre CRBs (Eq. S27, S31, the r -> 0 limit, finite zero depth)
  estimators.py          MLE (vectorized grid + pattern search), LMS, mLMS
  montecarlo.py          Monte Carlo driver, bootstrap standard errors
  camera.py              camera CRB (ideal and pixelated, background conventions)
  experiments.py         drivers: iterative MINFLUX, eps x L sweeps, misalignment studies
scripts/
  _paperconfig.py        every shared parameter (seed 42, fwhm 300 nm, L, N, MC sizes, ...)
  _paperstyle.py         figure style, parameter-hashed Monte Carlo cache, quick-mode routing
  fig_<n>_<name>.py      one script per figure -> paper/figures/fig<n>_<name>.pdf + data/fig<n>_summary.json
  make_all_figures.py    runs every figure script and reports script -> pdf -> status
  run_iterative_sweep.py iterative MINFLUX vs photon budget -> data/iterative_sweep.json
  compute_paper_numbers.py  the ONLY writer of data/paper_numbers.json and paper/generated/numbers.tex
  check_provenance.py    manuscript provenance check (see below)
  verify_crb_closed_forms.py  numerical check of the closed forms
  reproduce.py / .sh     the full pipeline
tests/                 unit tests + tests/test_acceptance.py (hash-pinned definition of "done")
data/                  paper_numbers.json, fig*_summary.json, iterative_sweep.json (data/mc/: cache)
structure/             claims.json (claim -> numbers -> status) and figures.json (figure -> script -> caption)
paper/                 main.pdf (the compiled manuscript: start here), main.tex (revtex4-2),
                       sections/, references.bib, provenance.json, generated/numbers.tex, figures/
docs/derivations/      derivation of the centre CRB of the TCP
docs/literature/       notes extracted from the published literature (equations and page numbers)
papers/README.md       the reference list (the PDFs are not redistributed)
OBJECTIVE.md           the author's statement of the study: questions R1-R5 and the deliverable
OBJECTIVE.template.md  the blank template OBJECTIVE.md was written from
pyproject.toml, requirements.txt, environment.yml   package metadata and dependencies (pip / conda)
LICENSE, CITATION.cff  MIT licence and citation metadata
equipo/2026-09-26_donut-localization/   the audit trail of how the study was produced (see below):
                       intent.md, state.json (ledger), inbox.jsonl, reports/, work/
job.cmd, job.sh        launchers of the agent-team `job` CLI (Windows / POSIX); jobs/ holds its
                       local, git-ignored job data
agent-team/, AGENTS.md, CLAUDE.md, .claude/, .agent-team/   the agent-team methodology and tooling
```

## Pipeline

| Script | Output |
|---|---|
| `scripts/fig_1_schematic.py` ... `scripts/fig_8_misalignment.py` | `paper/figures/fig<n>_<name>.pdf`, `data/fig<n>_summary.json` |
| `scripts/make_all_figures.py [--no-cache] [--quick]` | all of the above, with a status table |
| `scripts/run_iterative_sweep.py` | `data/iterative_sweep.json` |
| `scripts/compute_paper_numbers.py` | `data/paper_numbers.json`, `paper/generated/numbers.tex` |
| `scripts/check_provenance.py` | pass/fail report for the manuscript |
| `scripts/reproduce.py` | the whole sequence above, plus the tests |

Monte Carlo results of the figure scripts are cached in `data/mc/*.npz`. Each cache stores a
hash of the `_paperconfig` parameters; a cache whose hash no longer matches is recomputed.
Code changes are not detected, so the final products are always regenerated with `--no-cache`.

## Figures

| Id | Script | Question | What it shows |
|---|---|---|---|
| fig1 | `fig_1_schematic.py` | R1 | LG vs vectorial donut profiles, the quadratic zero and its curvature, TCP geometry |
| fig2 | `fig_2_vectorial.py` | R1 | Vectorial donut: handedness and linear polarization (zero depth 0 / 0.845 / 0.372), centre CRB vs L against LG beams; the wrong hand raises the centre CRB 16-90x and linear polarization 7-40x (L = 150-50 nm) |
| fig3 | `fig_3_crb_maps.py` | R2 | CRB maps of the TCP; the discontinuity at the centre (r -> 0 limit vs the Eq. S27 point value) and its removal by background |
| fig4 | `fig_4_scaling.py` | R2 | Scaling of the centre CRB with L, N and SBR; MINFLUX vs camera (ideal and pixelated) |
| fig5 | `fig_5_estimators.py` | R3 | MLE, LMS and mLMS bias and sigma/CRB inside and outside the TCP (off-centre the linearized estimators fall below the CRB by compression - the LMS already a few nm from the centre, the mLMS only beyond x0 ~ 15 nm); background-free MLE superefficiency and bias |
| fig6 | `fig_6_iterative.py` | R4 | Iterative MINFLUX vs photon budget against the ideal camera; per-iteration precision |
| fig7 | `fig_7_zero_depth.py` | R5 | Finite zero depth: optimal L, validity of L_opt ~ 0.78 fwhm sqrt(eps), off-centre CRB minimum |
| fig8 | `fig_8_misalignment.py` | R5 | TCP misalignment: bias and precision of the naive vs the honest MLE; mean naive bias 0.75-0.78 δ at the centre and ~0.80-0.86 δ at (L/4, 0) for δ = 2-10 nm (L = 100 nm, SBR = 10) |

The questions (from `OBJECTIVE.md`): R1 beam model (scalar LG vs vectorial donut), R2
Cramer-Rao bound (maps, closed forms, scaling, camera comparison), R3 estimators (MLE vs
LMS/mLMS), R4 iterative MINFLUX vs camera, R5 non-idealities (zero depth, background,
misalignment). The captions are in `structure/figures.json`.

## `data/` and the number registry

`data/paper_numbers.json` maps each key to
`{"value", "unit", "script", "description"[, "se"]}`. `se` is the Monte Carlo standard error,
from a bootstrap for standard deviations. All parameters come from `scripts/_paperconfig.py`
(seed 42), and the file is fully regenerated, deterministically, by
`scripts/compute_paper_numbers.py`. The same script writes `paper/generated/numbers.tex`, which
defines `\pnum{key}` (value) and `\pnumse{key}` (standard error). An unknown key typesets as a
visible `??key??` marker. Formatting (`fmt_value` / `fmt_value_se` in
`compute_paper_numbers.py`, unit-tested): a value without SE gets 4 significant figures with
trailing zeros kept (`1.600`), scientific notation below 1e-3 and from 1e5; a value with SE has
the SE rounded to 1-2 significant figures and the value rounded to the same decimal place
(`0.991 ± 0.005`).

`structure/claims.json` groups the keys into claims. Each claim records its statement, numbers,
script, figure, question, verification status and caveat. Every key belongs to at least one
claim, and every claim has at least one number.

The figure summaries `data/fig*_summary.json` repeat the numbers printed by each figure script.
Their shared keys are identical to the registry.

## Testing, acceptance and provenance

```sh
python -m unittest discover -s tests          # full suite
python -m unittest tests.test_acceptance -v   # the acceptance test ("done")
python scripts/check_provenance.py            # manuscript provenance
```

`tests/test_acceptance.py` was written before the work started and is pinned by its sha256 in
the team ledger. It recomputes the key physics with its own minimal implementation, independent
of `donutloc`.

`scripts/check_provenance.py` does the following:

- It flattens `paper/main.tex`, following `\input{...}`, `\include{...}` and the brace-less
  `\input name`.
- It requires every `\src{key}` to resolve to a reproducible entry of `paper/provenance.json`.
- It requires every `\pnum{k}` or `\pnumse{k}` to be a registry key (with an `se` for
  `\pnumse`).
- It requires every cited number to be backed by a `\src{e}` whose provenance entry lists `k`.
- It checks that `numbers.tex` is in sync with the registry.
- It checks that claims and registry cover each other.

## How this study was made

The study was carried out with the [agent-team](https://github.com/matiaszaldarriaga/agent-team)
methodology (see `AGENTS.md`). A bounded team of AI agents worked in short rounds: a PI
planning each round, workers computing and coding, and a permanent verifier and code reviewer.
The team produced one deliverable and stopped after each round for the author's decision.

- **Independent verification.** A result entered the manuscript only after a verifier with no
  stake in it reproduced it by an independent route. The verifier used its own Fisher and MLE
  code, its own seeds, limits and special cases. Refuted or unclear results stay open until
  resolved; the manuscript lists the remaining ones as open points.
- **Traceability.** The full record is versioned in `equipo/2026-09-26_donut-localization/`:
  - `state.json`: every claim with its status and round, the plan and the check results;
  - `inbox.jsonl`: every human direction;
  - `reports/`: what each role wrote in each round.
- **Provenance rule.** No number in the text is typed by hand. It is `\pnum{key}\src{key}`,
  checked by `scripts/check_provenance.py`, and the acceptance test runs that check.

The repository was created from the agent-team template (`agent-team-template`), which provides
`agent-team/`, `AGENTS.md`, `CLAUDE.md`, `.claude/`, `.agent-team/`, `job.cmd`/`job.sh` and
`OBJECTIVE.template.md`. The template's own installer scripts (`apply-to-existing.*`) are not
part of this study and were removed.

The literature notes in `docs/literature/` cite equations and pages of published papers. The
PDFs themselves are not redistributed (see `papers/README.md`).

## License

MIT. See `LICENSE`. For citation metadata, see `CITATION.cff`.
