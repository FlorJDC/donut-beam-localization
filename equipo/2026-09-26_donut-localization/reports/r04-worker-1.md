# Round 4: Worker 1 (tooling, data, README)

Scope: item "Tarea 1 / Worker 1" of `r04-pi.md`, parts (a) to (k). I touched only my own files: `scripts/**`, `src/**`, `tests/**` except `test_acceptance.py`, `data/**`, `structure/**`, `numbers.tex` (written only through `compute_paper_numbers`), figure PDFs (written only through `make_all_figures`), `README.md`, `CITATION.cff` and `.gitignore`. I did not touch the writer's files or `docs/private/`, and I made no commit.

## Order of work

I did (f) first, then (d) and (c), because the writer depends on them.

The population keys went into `data/paper_numbers.json` early, through a merge that calls `compute_paper_numbers.misalignment_population_section`. The full regeneration in (k) reproduces them bit for bit.

## (a) D1: cache invalidation (`scripts/_paperstyle.py`)

- `cached(name, fn, quick, force, params=None)` stores `_cache_key` in each npz. The key is the first 16 hex characters of the sha256 of three things:
  - every upper-case constant of `_paperconfig`, paths excluded;
  - the caller's `params`;
  - the quick flag.
- A missing or different key triggers a recompute. Legacy npz files without a key are also recomputed.
- Code changes are not detected (documented), so the final run uses `--no-cache`.
- Tests in `tests/test_paper_tooling.py::TestCacheInvalidation`, 6 tests:
  - reuse when nothing changed;
  - recompute when `params` changes;
  - recompute when `_paperconfig.MIS.n_patterns` changes: the exact r3 trap, including switching back;
  - legacy npz without a key;
  - force/quick separation;
  - the hash ignores paths.

## (b) D2: `--quick` never touches the final outputs

- **Figure scripts.** `_paperstyle.parse_args` calls `set_quick`. In quick mode:
  - `savefig` writes to `paper/figures/quick/`;
  - `write_summary` writes to `data/quick/` and adds `"quick": true`.
- **`make_all_figures --quick`** expects its PDFs in `paper/figures/quick/`.
- **`compute_paper_numbers --quick`** writes `data/quick/paper_numbers.json` (with top-level `"quick": true`) and `paper/generated/quick/numbers.tex`.
  - `numbers_tex` skips entries that are not dicts.
  - `check_provenance` fails if the final `paper_numbers.json` contains `"quick"`.
- `.gitignore` ignores the quick directories and `data/mc/*_quick.json`.
- **Tests** (`TestQuickOutputs`, 6 tests):
  - savefig/summary routing;
  - `parse_args`;
  - `make_all_figures` paths;
  - `compute_paper_numbers.main(["--quick"])` with a mocked `compute` writes only the quick files;
  - a subprocess run of `fig_1_schematic.py --quick` leaves the final PDF and summary byte-identical.
- **Real run.** `python scripts/reproduce.py --quick --skip-tests --no-latex` exits 0 in 109 s. The sha256 of every `paper/figures/*.pdf`, every `data/*.json` and `numbers.tex` is unchanged after it.

## (c) D3: `scripts/check_provenance.py`

- **Brace-less input.** The flattener now resolves `\input name`: the file name ends at whitespace, `{`, `}`, `\` or `%`. `\inputencoding` is not matched.
- **Whitespace.** `normalize()` removes whitespace between `\src`, `\pnum` or `\pnumse` and `{`. The agent-team checker therefore sees `\src {k}`.
  - `\pnum{ k}` (space inside the braces) is an error: the csname would include the space.
- **Macro parameters.** Arguments such as `#1` or `##1` inside macro bodies, and `\input{#1}`, are skipped.
- **New global rule.** Every `\pnum{k}` or `\pnumse{k}` needs some `\src{e}` in the flattened tex whose provenance entry has `k` in `number_keys` or `detail == k`. Comma lists like `\src{a,b}` are handled.
- **Claims coverage, both directions (h).** Every claim must have a non-empty `numbers` list, and every registry key must belong to at least one claim.
- **Quick output.** A `--quick` numbers file is rejected.
- **Tests** in `tests/test_check_provenance.py`, 28 tests, 18 of them new. They build copies in a temporary `--root` and cover:
  - brace-less input, as a hidden bad file, a good file and a missing file;
  - `\src {nowhere}`;
  - `\pnum {k}\src {k}` passes;
  - `\pnum{ k}` fails;
  - macro `#1` and `##1`;
  - `\pnum` or `\pnumse` without `\src` fails;
  - a global `\src` elsewhere passes;
  - a `\src` whose entry does not cover the key fails;
  - `detail` alone covers the key;
  - comma lists;
  - quick marker, empty claim, key without a claim;
  - `normalize`.
- **Writer convention.** `\pnum{k}\src{k}` with `detail=k` and `number_keys=[k]` passes. On the real manuscript: 178 tags, 178 entries, all checks pass.

## (d) D4: `numbers.tex` fallback

- The fallback is now `\textbf{??\detokenize{#1}??}` (and `??se:\detokenize{#1}??`).
- `TestNumbersTexFallback::test_unknown_key_compiles_with_tectonic` compiles `\pnum{no_such_key}`, `\pnumse{no_such_key}` and `\pnum` in math mode. It uses Tectonic from `$TECTONIC` or `PATH` and skips when neither is set. Result: exit 0, no "Missing $".
- For contrast, the old macro gives exit 1 in the same document; I ran that by hand.
- A string test checks that the `\detokenize` form is used.

## (e) Constants moved to `_paperconfig`

| Constant | Replaces |
|---|---|
| `N_REP_BIAS0 = 40000` | `rb`, and the bias-scan repetitions in fig 5 |
| `MLE_RADIUS_CENTRE_OVER_L = 1.0`, `MLE_RADIUS_SWEEP_OVER_L = 2.0` | the MLE search radii L / 2L |
| `EPS_TRANSITION = 0.002` | `zero_depth_transition_scale_nm` and `EPS_C` in fig 7 |

- Also new: `MIS_POP_N_PATTERNS = 4000`, `MIS_POP_DELTAS = (2, 5, 10)` and `QUICK_DIRNAME`.
- They are used in `compute_paper_numbers.py`, `fig_5_estimators.py` and `fig_7_zero_depth.py`.
- `TestSharedConstants` checks that the literals are gone from those scripts.
- The regeneration confirms that no value changed.

## (f) Noise-free population slope (inbox r4)

- **New function** `donutloc.experiments.misalignment_population_bias`. It uses the same misalignment model, pattern draws (`default_rng(seed)` re-created for each delta) and naive MLE (disk of radius 0.75 L) as `misalignment_study`.
  - The difference: the MLE is applied to the expected counts `N p_true(r)`, so there is no Poisson noise.
  - SE = std over patterns / √P.
  - `slope` = sum(delta·b) / sum(delta²), unweighted through the origin over delta > 0, with its SE taken from the per-pattern slopes.
  - It takes about 9 s for 4000 patterns × 3 deltas × 2 positions.
- **Keys**, from `compute_paper_numbers.misalignment_population_section`, with description, `script` and SE. The last column is the r3 verifier's independent value (own code, other patterns).

| key | value ± SE | r3 verifier |
|---|---|---|
| misalignment_pop_n_patterns | 4000 | 4000 |
| misalignment_naive_bias_over_delta_pop_center | 0.77683 ± 0.00497 | — |
| misalignment_naive_bias_over_delta_pop_center_d2 | 0.74713 ± 0.00521 | 0.746 |
| misalignment_naive_bias_over_delta_pop_center_d5 | 0.75493 ± 0.00513 | 0.755 |
| misalignment_naive_bias_over_delta_pop_center_d10 | 0.78349 ± 0.00495 | 0.784 |
| misalignment_naive_bias_over_delta_pop_Lq | 0.85070 ± 0.00742 | — |
| misalignment_naive_bias_over_delta_pop_Lq_d2 | 0.80229 ± 0.00595 | 0.794 |
| misalignment_naive_bias_over_delta_pop_Lq_d5 | 0.80734 ± 0.00610 | 0.802 |
| misalignment_naive_bias_over_delta_pop_Lq_d10 | 0.86347 ± 0.00810 | 0.871 |

- **Warning for the PI and writer: the slopes are not 0.75 and 0.83.** The slope defined in the plan (unweighted LSQ through the origin over δ = 2, 5, 10) gives **0.777 at the centre and 0.851 at (L/4, 0)**. The plan's reference values "≈0.75 / ≈0.83" do not come out of this definition. The fit is dominated by δ = 10, where the ratio is largest.
  - The per-δ values agree with the verifier within about 1 SE (independent patterns).
  - For comparison, the plain mean of the three ratios would be 0.762 and 0.824.
  - I implemented the plan's definition. The writer should quote the per-δ ratios next to the slope, as the claim caveat says.
- **Paired check against the MC.** With the same first 400 patterns as the 400×200 MC, the noise-free |bias| is:
  - centre: 1.492 / 3.759 / 7.760 nm, against MC 1.492 / 3.731 / 7.708 (±0.033 / 0.080 / 0.158);
  - (L/4, 0): 1.654 / 4.158 / 9.151, against MC 1.691 / 4.202 / 9.304.
- **Conclusion from the paired check.** Poisson noise barely changes the naive bias. The spread between the MC draws at (L/4, 0) (0.814 / 0.836 / 0.862) is pattern sampling: the population slope over the first 400 patterns is 0.896, over 1000 it is 0.876, and over 4000 it is 0.851. This resolves the r3 unclear about 0.862 ± 0.012: that SE understated the pattern-draw spread.
- The old MC key `misalignment_naive_bias_over_delta_Lq` is kept, and is not cited: its claim is `superseded`.
- **Tests** (`tests/test_experiments.py::TestMisalignmentPopulation`, 5 tests):
  - delta = 0 gives zero bias;
  - convergence as delta → 0: the ratios at 0.25 and 0.5 agree to < 2 % and are within 5 % of the δ = 2 value;
  - slope definition, and pattern identity with `misalignment_study`;
  - the paired comparison with the stored MC centre values, within 3 SE for each δ and for the MC-weighted slope;
  - input validation.

## (g) `structure/figures.json` and figure scripts

- **fig5, rewritten.** Four panels:
  - (a) and (b) are the sweep, and mention the MLE outlier bump (σ/CRB ≈ 1.5 at x0 = 10–20 nm, N = 100), which depends on the 2L search radius and disappears at N = 1000;
  - (c) is the background-free superefficiency;
  - (d) is the bias at r = (2, 0), about −0.34 nm.
- **fig6, rewritten descriptively:**
  - (a) σ(N) against the camera and the CRB, with slope −0.515 (N ≥ 500) and the SBR = 10 and adaptive points; the grey no-re-centring cross is labelled an artefact;
  - (b) σ/σ_cam from 0.183 to 0.158;
  - (c) σ per iteration against L_k.
- **fig7, rewritten descriptively.** 0.78 holds only for ε ≲ 0.01. Panel (c) is in nm and marks the 4.9 nm scale.
  - `fig_7_zero_depth.py` already labels (c) "σ_CRB (nm)". The verifier looked at an older PDF (07:01), and the regenerated PDF has the unit. The only change to that script is the constant.
- **fig8, rewritten:**
  - 400 × 200;
  - the drawn line is 0.76 δ, a weighted fit over 7 deltas including δ = 15;
  - the principal slope is the population value (≈0.78 centre, ≈0.85 L/4, growing with δ);
  - the loss is bias;
  - the honest MLE is at the floor only at the centre;
  - at L/4 the honest σ is 1–2.4 % above the CRB.
- The `fig_8_misalignment.py` docstring now says 400 patterns. I did not add the optional L/4 line.

## (h) `structure/claims.json`: 39 claims

- The 15 claims marked `pending-r3-verification` are now `verified` (`verified_round` 3), each confirmed by the r3 verifier's X1-* claims in `state.json`.
- **Empty `numbers` lists filled:**

| Claim | Numbers added |
|---|---|
| `crb_discontinuity_origin` | the S27 and limit keys |
| `lms_sbr_shrink` | new key `lms_sbr10_shrink_factor` = 0.909091 = 10/11 |
| `eps_constant_pedestal_is_background` | new keys `eps_constant_pedestal_sbr_eff_L100_eps0p05` = 2.9075, `eps_constant_pedestal_over_s31_L100_eps0p05` = 1.0000000000000004, `eps_gaussian_pedestal_over_s31_L100_eps0p05` = 1.00532 |
| `misalignment_sigma_se_fix` | the d10 sigma keys |

- The two claims that gained new keys carry `numbers_added_round: 4`. Their statements were already verified in r1, but these key values are new and should be checked in r5.
- **Coverage:** every one of the 217 registry keys belongs to a claim, checked by `check_provenance`.
- **New claims:**
  - `misalignment_naive_bias_population`, status `pending-r4-verification`;
  - `misalignment_naive_bias_Lq_mc_slope`, status `superseded`.
- All caveats are now in English, taken from the PI's "Caveats" list, including the `\pnumse` notes for keys with SE > 2 %. The 0.79 caveat is marked superseded. No non-ASCII characters remain.

## (i) Reproduction scripts and citation

- `scripts/reproduce.py` runs, in order:
  1. tests;
  2. optionally `--with-sweep`;
  3. `make_all_figures --no-cache`;
  4. `compute_paper_numbers`;
  5. `check_provenance`;
  6. acceptance;
  7. `tectonic paper/main.tex` if found on `PATH`, `$TECTONIC` or `--tectonic`.
- Options: `--quick`, `--skip-tests`, `--no-latex`, `--dry-run`.
- `scripts/reproduce.sh` is a POSIX wrapper.
- `TestReproduce` has 3 tests.
- `CITATION.cff` (CFF 1.2.0): author alias flor-choque, website github.com/FlorJDC, repository URL, MIT. It has no affiliation.

## (j) `README.md`

- The new English README replaces the Spanish agent-team template. I deleted the template: `AGENTS.md` and the upstream template repository already cover it.
- Sections:
  - title and summary, quick start;
  - repository structure;
  - pipeline table and figure table (id, script, R1–R5, content);
  - `data/` and the registry, testing and acceptance and provenance;
  - "How this study was made" (agent-team, independent verification, `equipo/.../state.json`, `inbox.jsonl`, `reports/`, provenance rule);
  - license.
- It has nothing from `docs/private`.

## (k) Regeneration and checks (real outputs)

- **`python scripts/make_all_figures.py --no-cache`:** 8/8 figures OK in 17 min 48 s. fig5 took 328 s and fig8 664 s. Every `data/fig*_summary.json` is **identical** to its pre-regeneration copy: 0 changed keys, 0 removed.
- **`python scripts/compute_paper_numbers.py`:** 217 keys in 376 s, exit 0.
- **Diff of `paper_numbers.json`** against the pre-round-4 copy (204 keys):
  - **0 keys changed** in value, SE or unit; 0 removed;
  - **13 new keys**: the 9 population keys above plus `lms_sbr10_shrink_factor`, `eps_constant_pedestal_sbr_eff_L100_eps0p05`, `eps_constant_pedestal_over_s31_L100_eps0p05` and `eps_gaussian_pedestal_over_s31_L100_eps0p05`.
  - The regenerated values of the new keys equal the early merged values bit for bit.
- **Summaries against the registry:** 149 shared keys, 0 mismatches.
- **`python scripts/check_provenance.py`**, exit 0:
  ```
  provenance: 178 tag(s), 178 entry(ies), 0 error(s), 0 warning(s)
  claims: 39, numbers: 217, provenance entries: 178
  all checks pass
  ```
- **`python -m unittest tests.test_acceptance -v`:** Ran 8 tests, OK. `test_manuscript_provenance_passes` now runs on the writer's real manuscript.
- **`python -m unittest discover -s tests`** (with `TECTONIC` set): **Ran 218 tests in 38.0 s, OK.** Without `TECTONIC` or tectonic on `PATH`, the one compile test skips.
- **sha256 of `tests/test_acceptance.py`:** 7d198853e87c7bd94bb8be07a855651146e1d849310281546cd8089dcf262303, unchanged.

## Not resolved / notes

- The population slope definition gives 0.777 and 0.851, not the "≈0.75 / ≈0.83" written in the plan (see f). The r5 verifier and the PI should confirm that the fit definition is the one intended.
- The four new demonstration keys for r1-verified claims need an independent check in r5.
- `numbers.tex` changed through the fallback macro and the 13 new keys, so `main.pdf` should be recompiled (orchestrator step 5.2).

```claims
[{"status": "pending", "text": "D1: _paperstyle.cached stores a sha256[:16] of all upper-case _paperconfig constants (paths excluded) + caller params + quick flag in each npz and recomputes on mismatch or missing key; tests TestCacheInvalidation (6) incl. changing MIS.n_patterns"},
 {"status": "pending", "text": "D2: --quick writes only paper/figures/quick/, data/quick/ (summaries and paper_numbers.json with top-level quick=true) and paper/generated/quick/numbers.tex; reproduce.py --quick left the sha256 of every paper/figures/*.pdf, data/*.json and numbers.tex unchanged; check_provenance rejects a final paper_numbers.json containing 'quick'"},
 {"status": "pending", "text": "D3: check_provenance flattens brace-less '\\input name', tolerates whitespace in '\\src {k}' / '\\pnum {k}', errors on '\\pnum{ k}', skips macro parameters (#1, ##1), and requires every \\pnum/\\pnumse key to be covered by some \\src{e} whose provenance entry has it in number_keys or as detail; 18 new tests in tests/test_check_provenance.py (28 total, all pass)"},
 {"status": "pending", "text": "D3/h: check_provenance fails on a claim with empty numbers or a paper_numbers key in no claim; on the real project: 39 claims, 217 numbers, 178 tags/entries, 'all checks pass', exit 0"},
 {"status": "pending", "text": "D4: numbers.tex fallback is \\textbf{??\\detokenize{#1}??}; with Tectonic a document using \\pnum{no_such_key}, \\pnumse{no_such_key} and \\pnum in math compiles (exit 0, no 'Missing $'), while the old macro exits 1"},
 {"status": "pending", "text": "Constants moved to _paperconfig (N_REP_BIAS0=40000, MLE_RADIUS_CENTRE_OVER_L=1, MLE_RADIUS_SWEEP_OVER_L=2, EPS_TRANSITION=0.002) and used by compute_paper_numbers.py, fig_5_estimators.py and fig_7_zero_depth.py; the full --no-cache regeneration changed no value"},
 {"status": "pending", "text": "Population (noise-free) naive-MLE bias, 4000 patterns, seed 42, L=100, SBR=10: |bias|/delta = 0.7471/0.7549/0.7835 (centre) and 0.8023/0.8073/0.8635 (L/4,0) for delta=2/5/10 (SE ~0.005-0.008); unweighted LSQ slope through the origin 0.7768+-0.0050 (centre) and 0.8507+-0.0074 (L/4,0); keys misalignment_naive_bias_over_delta_pop_* and misalignment_pop_n_patterns=4000"},
 {"status": "pending", "text": "With the same first 400 patterns as the 400x200 MC, the noise-free naive |bias| (centre 1.492/3.759/7.760, L/4 1.654/4.158/9.151 nm) matches the MC (1.492/3.731/7.708 and 1.691/4.202/9.304) within 1 MC SE; the (L/4,0) MC slope spread (0.814-0.862) is pattern sampling (population slope 0.896/0.876/0.851 over 400/1000/4000 patterns)"},
 {"status": "pending", "text": "New demo keys: lms_sbr10_shrink_factor=0.9090909 (=10/11), eps_constant_pedestal_over_s31_L100_eps0p05=1.0000000000000004, eps_gaussian_pedestal_over_s31_L100_eps0p05=1.00532, eps_constant_pedestal_sbr_eff_L100_eps0p05=2.90751"},
 {"status": "pending", "text": "claims.json: 39 claims (37 verified, 1 pending-r4-verification, 1 superseded), every claim has numbers, all 217 keys covered, 15 pending-r3 set to verified round 3, caveats in English with no non-ASCII"},
 {"status": "pending", "text": "Regeneration: make_all_figures --no-cache 8/8 OK (17m48s), all fig*_summary.json identical to the pre-round versions; compute_paper_numbers 217 keys (376 s): 0 of the 204 existing keys changed value/SE/unit, 13 new; 149 shared summary keys match the registry"},
 {"status": "pending", "text": "Checks: unittest discover 218 tests OK (38.0 s, TECTONIC set); tests.test_acceptance 8/8 OK; check_provenance exit 0 'all checks pass'; sha256 tests/test_acceptance.py = 7d198853...2303 unchanged"}]
```
