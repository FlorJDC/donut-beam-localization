# Ronda 3 — code-reviewer

Qué revisé: el diff `7b7fe98..61efca2`, el árbol de trabajo (la fig 8 regenerada con 400 patrones, a las 07:26) y el ledger (`state.json`, `inbox.jsonl`, `r03-pi.md`). Scratch en `work/review/r03/`: `r2fixes.py`, `prov_scen.py`, `regen.py`, `regen_mis.py` y `tex/`.

No edité código del proyecto ni leí `docs/private/`.

## Corridas reales

- **Hash.** `sha256 tests/test_acceptance.py` = `7d198853…2303`: sin cambios.
- **Suite.** `python -m unittest discover -s tests`: **177 tests OK**, en 38.9 s.
- **Aceptación.** `python -m unittest tests.test_acceptance -v`: **8/8 OK**, sin SKIP.
- **Procedencia.** `python scripts/check_provenance.py` imprime `all checks pass`: 37 claims, 204 números y 0 entradas de procedencia.
- **Regeneración en memoria de todas las secciones de `compute_paper_numbers.py`** (sin escribir archivos), comparada con el JSON commiteado:
  - 137 claves no-desalineación: **0 diferencias**, valores y SE bit a bit, en 70 s;
  - 67 claves de desalineación: **0 diferencias**, en 505 s.
  - Las 204 claves son reproducibles de forma determinista con semilla 42.
- **Fig 8 fresca.** `fig8_n_patterns = 400`. Las 65 claves compartidas con `paper_numbers.json` coinciden exactamente (valor y SE).
- **Cruce del resto de las figuras.** Las claves compartidas de fig1–fig7 coinciden exactamente con `paper_numbers.json`, incluidas `mle_efficiency_center` y `mle_nobg_bias_x_r2_nm` de fig5.

## Las 4 refutadas de r2: escenario original re-ejecutado

1. **`misalignment_study` SE.** Escenario: δ=10, centro, 20×200.
   - Resultado: `sigma_se` honesto 0.0454 frente a `within` 0.0144; ingenuo 0.0538 frente a 0.0156.
   - Ahora el SE es la std entre patrones dividida por √P. Coincide con lo que r2 midió (0.038 y 0.056).
   - Cubierto por `test_misalignment_sigma_se_between_patterns` y `_formula_exact`, que lo recalcula de forma independiente.
   - **Resuelta.**
2. **`crb(..., zero_policy="limit")` con N en array.**
   - `crb(p,[0,0],[100,400],"limit")` da `[1.60510, 0.80255]`, igual a `crb_limit` en cada N; forma (2,).
   - También probé: r (2,2) con N (2,) emparejados, y la grilla con N escalar. Correcto.
   - Tests en `TestLimitPolicyBroadcastN`.
   - **Resuelta.**
3. **Docstring de `l_schedule`.**
   - Ahora dice que σ_k es el CRB central y lo subestima (3.47 frente a 4.26), y da ~97 % de cobertura (2.65 % afuera).
   - `adaptive_frac_outside_next_radius` = 0.0265 ± 0.0016 en `paper_numbers`.
   - Test: `test_adaptive_schedule_coverage_is_about_97_percent`.
   - **Resuelta.**
4. **Beam vectorial lineal.** Ahora es por defecto una tabla de armónicos angulares con spline cúbico en ρ².
   - CRB(7,3): 62.35509 (harmonic) frente a 62.35509 (exact). `crb_limit`: 60.77420 frente a 60.77420. Diferencia relativa < 1e-9.
   - Probé también fuera del test, con `rho_max=300`, contra `mode="exact"`:
     - `pol_angle=0.7`: 9e-11;
     - `charge=2`: 3e-11 (el límite |n|≤4 vale para toda carga, porque los órdenes de campo van de l−2 a l+2);
     - `eps=0.01`: 5e-11;
     - cerca de ρ=0 (extrapolación del nodo 0): exacto.
   - `linear_method="cartesian"` avisa con `UserWarning`, y hay un test del sesgo.
   - **Resuelta.**

Unclear de r2 resuelta: el SE del iterativo ahora es bootstrap (`bootstrap_sigma_se`, 2000 remuestreos). Los tests cubren el caso gaussiano ±10 %, colas pesadas, NaN y reproducibilidad.

## Defectos encontrados (con escenario)

### D1. Cache de figuras sin invalidación: trampa real, que ya ocurrió en esta ronda
`_paperstyle.cached(name, …)` carga `data/mc/<name>.npz` solo por nombre. No compara `_paperconfig` ni los parámetros de `compute()`.

Escenario **observado** a las 07:17, antes de la regeneración:
- `_paperconfig.MIS` ya tenía `n_patterns=400`;
- `data/mc/fig8_misalignment.npz` (06:56) tenía `n_patterns=50`;
- `fig8_summary.json` y el PDF salieron con 50 patrones;
- `misalignment_naive_sigma_Lq_d10_nm` daba 3.128 en la figura frente a 3.739 en `paper_numbers`: **−16 %**. `naive_bias_abs_Lq_d10` daba −11.6 %.

Se arregló solo porque alguien corrió `--no-cache` a mano. Lo mismo le pasa a fig2–fig7 si cambia FWHM, L_REF, N_REP_MLE, EPS_LIST, etc.

Arreglo sugerido: guardar en el npz un hash de los parámetros relevantes (o de `_paperconfig` más el `compute` usado) y recalcular si no coincide. Como mínimo, `make_all_figures` debería pasar `--no-cache` en la corrida final, y R4/verifier deberían regenerar así.

### D2. `--quick` sobrescribe los productos finales commiteados
- **Figuras.** `savefig(fig, "fig8_misalignment")` y `write_summary("fig8")` usan la misma ruta con o sin `--quick`: solo el npz lleva el sufijo `_quick`.
  - Entonces `python scripts/make_all_figures.py --quick` reemplaza los 8 PDF de `paper/figures/` y los 8 `data/fig*_summary.json` versionados por versiones con 1000 reps o 10×100 patrones.
  - `make_all_figures` los reporta como "OK": el chequeo de frescura compara mtime.
- **Números.** `compute_paper_numbers.py --quick` escribe **`data/paper_numbers.json` y `numbers.tex`** con MC reducido: `mle_efficiency_center` con 1000 reps y SE ~1.6 %, desalineación con 10×100.
  - El JSON no lleva ninguna marca `quick`, así que nada lo detecta después (la aceptación seguiría pasando).
  - `run_iterative_sweep.py` sí lo hace bien: `--quick` va a `data/mc/iterative_sweep_quick.json`.
- Arreglo sugerido: en quick, escribir a `*_quick.pdf`, `*_quick_summary.json` y `data/mc/paper_numbers_quick.json`, o como mínimo agregar `"quick": true` y que `check_provenance` falle si está presente.

### D3. `check_provenance.py`: huecos en el aplanado y en la exigencia de `\src`
Probados con el árbol temporal de `tests/test_check_provenance.py`:
- **`\input sections/other` (sin llaves, sintaxis válida de TeX) no se aplana.** El archivo incluido tenía `\src{nowhere}` y `\pnum{zzz}`, y el resultado fue **ok=True**. `_INPUT_RE` solo reconoce `\input{…}`. Tampoco reconoce `\subfile`, `\import`, ni `\input` dentro de macros.
- **`\pnum{k}` sin `\src{…}` pasa.** El plan dice "R4 escribe `\pnum{k}\src{k}`", pero no hay chequeo que lo exija. Un número citado sin entrada en `provenance.json` pasa en verde, contra el principio 3.
  - Sugerencia: exigir que cada `\pnum{k}` tenga una entrada de procedencia cuyo `number_keys` contenga k (o `\src{k}` adyacente).
- **Menor.** `\pnum{ crb_x_nm}` (con espacio) pasa el checker porque hace `.strip()`, pero en LaTeX el csname incluye el espacio y el número no se resuelve.
- Lo que sí funciona, verificado:
  - `\src` sin resolver falla (lo delega a `agent-team/bin/check_provenance.py`, invocado como subproceso con el tex aplanado);
  - `\pnum` desconocido falla;
  - `\pnumse` sin `se` falla;
  - `numbers.tex` desincronizado falla;
  - `\input` inexistente falla;
  - los scripts faltantes en `claims.json` fallan;
  - los comentarios (y `\%`) se tratan bien;
  - recursión y ciclos, bien.
- Nota: las rutas `reproduce` de `provenance.json` se resuelven relativas al cwd del proyecto o a `/tmp` (donde está el aplanado) y sus padres, **no** relativas a `paper/`. Una ruta escrita relativa a `paper/` falla. R4 debe usar rutas relativas a la raíz.

### D4. Fallback de `\pnum` para una clave inexistente: rompe la compilación en vez de marcarla
Compilé con Tectonic (XeTeX), en `work/review/r03/tex/`.
- `main.tex` (article) y una copia del `paper/main.tex` real (revtex4-2) con una sección que usa `\pnum` en texto, modo matemático, `\section` (que pasa por `\MakeUppercase`) y `\caption`, más `\pnumse`: **compilan bien**. Salida: "1.605", "7.042 × 10⁻⁵", "0.9912±0.005131", "CRB OF 1.605 NM".
- El diseño con `\ifcsname`/`\csname` y `\providecommand` funciona en XeTeX (e-TeX). pdfTeX también tiene `\ifcsname`; no pude probarlo porque no hay pdflatex local.
- Pero `\pnum{no_such_key}` da **`! Missing $ inserted`**: el fallback `\textbf{??#1??}` tipea el `_` de la clave en modo texto, y todas las claves tienen `_`. Con el checker verde no debería ocurrir, pero el marcador visible no funciona tal como está diseñado.
- Arreglo: `\textbf{??\detokenize{#1}??}`.

### D5. Menores: sin escenario de número incorrecto hoy
- **Constantes duplicadas fuera de `_paperconfig`.**
  - `rb = 40000` (sesgo del MLE en r=(2,0)) está a mano en `compute_paper_numbers.py` y en `fig_5_estimators.py`;
  - `eps=0.002` está literal en `zero_depth_transition_scale_nm`;
  - la búsqueda del MLE con radio `L`/`2L` está duplicada en ambos scripts.
  - Hoy coinciden (regeneración bit a bit), pero cambiar uno solo desincroniza la figura y el número.
- `_fmt` imprime valores como `-0.49999999999999983` → "-0.5". Bien. Enteros guardados como float (`fwhm_nm` = 300.0) salen "300". Bien.
- `claims.json`:
  - 4 claims sin `numbers` (`crb_discontinuity_origin`, `lms_sbr_shrink`, `eps_constant_pedestal_is_background`, `misalignment_sigma_se_fix`);
  - 34 claves de `paper_numbers` no están en ningún claim;
  - `photons_for_5nm` y otros caveats están en castellano, con caracteres no ASCII (UTF-8 correcto).
  - El esquema es consistente: 30 con `figure`, 7 sin él, todas las `figure` existen en `figures.json` y no hay claves duplicadas. 15 claims llevan `pending-r3-verification`.
- **SE por encima del 2 %.** Los `bias_abs` de desalineación están entre 2.1 y 3.0 %; `misalignment_naive_sigma_Lq_d10` en 4.3 %; `mle_nobg_bias_x_r2_nm` en 2.2 %; `adaptive_frac_outside` en 6 % (binomial).
  - Está cubierto por la directiva del inbox r3: se citan con su SE explícito. R4 debe usar `\pnumse` en esos casos.
  - El `|bias|` promedio honesto (≈0.17–0.22 nm) es el piso MC (≈0.16–0.21), no un sesgo.
- **`make_all_figures.py`.** Detecta fallas (exit≠0 o PDF no renovado) y sale con 1. Bien, salvo D2.
- **Semillas.** Todas las figuras usan `C.SEED`, sin RNG sin semilla. Solo se usa la API pública de `donutloc`.
- `scripts/mc/` es un directorio vacío suelto (no versionado); inofensivo.

## Privacidad
Busqué en todo el diff (`p-minflux`, `private`, `unpublished`, mediciones) y no encontré contenido sospechoso versionado. Las únicas menciones a `docs/private/` son las de `OBJECTIVE.md` y el inbox de r1, que ya existían y solo nombran la ruta.

## Contra el intent: lo que falta para R4 (manuscrito)
Faltan los siguientes entregables de `OBJECTIVE.md`:
- `paper/sections/*.tex`: hoy el `\input` está comentado y `sections/` está vacío;
- `paper/references.bib`;
- `paper/provenance.json`: hoy `{}`, con 0 entradas y 0 tags;
- `scripts/reproduce.sh`;
- `README.md`: **todavía es el README de la plantilla agent-team**; falta el quick start, la estructura y la tabla de figuras.

La aceptación pasa trivialmente en `test_manuscript_provenance_passes`, porque el manuscrito está vacío.

Recomendaciones para R4:
- arreglar D2 y D3 (exigir procedencia por cada `\pnum`) **antes** de escribir, y D4 (una línea);
- regenerar figuras y números con `--no-cache` como paso final (D1);
- los 15 claims `pending-r3-verification` deben pasar por el verificador antes de entrar al texto.

```claims
[{"status": "verified", "text": "r2 fix misalignment_study.sigma_se: re-run original scenario (delta=10, centre, 20x200) gives SE between patterns 0.0454/0.0538 vs within 0.0144/0.0156 (honest/naive); tests test_misalignment_sigma_se_between_patterns and _formula_exact cover it"},
 {"status": "verified", "text": "r2 fix fisher.crb zero_policy='limit' with array N: crb(p,[0,0],[100,400],'limit') = [1.60510,0.80255] = crb_limit at each N, shape (2,); paired (2,2)x(2,) and grid cases correct; TestLimitPolicyBroadcastN"},
 {"status": "verified", "text": "r2 fix l_schedule docstring: now states centre CRB 3.47 underestimates real 4.26 nm and ~97% coverage (2.65% outside); tested and matched by adaptive_frac_outside_next_radius=0.0265+-0.0016"},
 {"status": "verified", "text": "r2 fix linear vectorial beam: harmonic tabulation default gives CRB(7,3)=62.35509 and crb_limit=60.77420 identical to mode='exact' (<1e-9); also exact to ~1e-10 for pol_angle=0.7, charge=2, eps=0.01 at rho_max=300; cartesian path warns"},
 {"status": "verified", "text": "Bootstrap SE (montecarlo.bootstrap_sigma_se, error_stats, run_mc sigma_se_boot, iterative sigma_se) implemented and tested (Gaussian +-10%, heavy tails, NaN rows, reproducibility); resolves r2 unclear on iterative SE"},
 {"status": "verified", "text": "Suite 177 tests OK (38.9 s); acceptance 8/8 OK; test_acceptance.py sha256 7d198853...2303 unchanged; check_provenance prints 'all checks pass'"},
 {"status": "verified", "text": "compute_paper_numbers.py is the sole writer of data/paper_numbers.json and numbers.tex; in-memory regeneration of all 204 keys (137 + 67 misalignment) reproduces the committed values and SEs bit for bit (seed 42, _paperconfig params)"},
 {"status": "verified", "text": "Figure summaries fig1-fig8 (fig8 regenerated, 400 patterns) agree exactly with paper_numbers.json on every shared key (fig8: 65 keys)"},
 {"status": "verified", "text": "numbers.tex macro design (\\providecommand + \\ifcsname/\\csname) compiles under XeTeX (Tectonic) with article and the real revtex4-2 main.tex: \\pnum in text, math, \\section, \\caption and \\pnumse render correctly (1.605, 7.042x10^-5, 0.9912+-0.005131)"},
 {"status": "refuted", "text": "_paperstyle.cached has no cache invalidation (keyed only by name): observed this round, fig8 npz with n_patterns=50 was reused after _paperconfig.MIS changed to 400, so fig8 disagreed with paper_numbers (naive_sigma_Lq_d10 3.128 vs 3.739, -16%) until a manual --no-cache; same trap for fig2-fig7 on any _paperconfig change"},
 {"status": "refuted", "text": "--quick overwrites committed final outputs: fig scripts write the same paper/figures/fig*.pdf and data/fig*_summary.json in quick mode (make_all_figures --quick replaces all 8 and reports OK), and compute_paper_numbers.py --quick overwrites data/paper_numbers.json + numbers.tex with reduced MC and no 'quick' marker, undetectable by acceptance/check_provenance"},
 {"status": "refuted", "text": "check_provenance.py flatten misses brace-less \\input: '\\input sections/other' whose file contains \\src{nowhere} and \\pnum{zzz} passes (ok=True); and \\pnum{k} without any \\src/provenance entry passes, so a cited number can lack provenance"},
 {"status": "refuted", "text": "numbers.tex fallback for an unknown key breaks compilation instead of showing a marker: \\pnum{no_such_key} -> '! Missing $ inserted' (underscore typeset in \\textbf{??#1??}); fix with \\detokenize{#1}"},
 {"status": "unclear", "text": "Constants duplicated outside _paperconfig (rb=40000 bias reps, MLE search radii L/2L, eps=0.002 literal) in both compute_paper_numbers.py and fig_5_estimators.py; consistent today, but a one-sided edit would desynchronise figure and paper number"},
 {"status": "unclear", "text": "structure/claims.json: 4 claims with empty numbers, 34 paper_numbers keys referenced by no claim, 15 claims pending-r3-verification; misalignment |bias| (2.1-3.0%), naive_sigma_Lq_d10 (4.3%), mle_nobg_bias (2.2%) exceed 2% SE and must be cited with \\pnumse per the r3 inbox directive"},
 {"status": "unclear", "text": "Missing vs OBJECTIVE for round 4: paper/sections/*.tex (input commented out), references.bib, provenance.json entries (currently {}), scripts/reproduce.sh, and README.md (still the agent-team template); acceptance provenance test passes only because the manuscript is empty"}]
```
