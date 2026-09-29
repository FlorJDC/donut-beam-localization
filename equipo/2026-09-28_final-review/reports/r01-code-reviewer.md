# r01 — code-reviewer (revisión final, rama final-review-and-explorer)

Revisé el código de los archivos, no las descripciones. Corrí la suite completa: `python3 -m unittest discover -s tests`
dio **254 tests OK (1 skipped) en 22.5 s**. Entorno: numpy 2.4.6, scipy 1.17.1, Python 3.11. No hay Python 3.8
en esta máquina, así que la compatibilidad con 3.8 y numpy 1.24 la revisé con grep y no la ejecuté. Busqué
`np.trapezoid`, `in1d`, `row_stack`, `float_`, `NaN`, `removeprefix`, `list[...]`, `| None`, `:=`, `zip(strict=)`,
`copy=` y no apareció ninguno en `src/`, `scripts/` ni `tests/`.

## [BLOQUEANTE]

1. **`scripts/make_all_figures.py:25-35` + `scripts/fig_9_pminflux_timing.py` / `scripts/fig_10_background_bias.py`:
   `reproduce.py` falla en el paso 2.** `discover()` toma todo `fig_*.py` y espera `paper/figures/fig9_pminflux_timing.pdf`
   y `paper/figures/fig10_background_bias.pdf`. Los scripts nuevos escriben en `informe/figures/fig_*.pdf`.
   El resultado es `FAIL (no new pdf)` en los dos, `make_all_figures` devuelve 1 y
   `python scripts/reproduce.py` falla, que es el comando que el README da como reproducción completa. Además suma el
   Monte Carlo de fig_10 (N_REP=4000 × 14 puntos × 5 MLE) al tiempo de "about 25 minutes". Lo confirmé
   ejecutando `discover()` y `pdf_for()`: los dos PDFs esperados no existen.
   Arreglo mínimo: excluir los scripts del informe en `discover()` (por ejemplo, solo `n <= 8`), o
   renombrarlos a `informe_fig_*.py` / moverlos a `informe/scripts/`, y documentarlos aparte.

2. **`scripts/fig_9_pminflux_timing.py:299-333`: ignora `--quick`.** `make_all_figures --quick` le pasa `--quick` y el
   script corre completo y **sobrescribe `informe/data/pminflux_timing.json` y `informe/figures/fig_pminflux_timing.*`**
   (productos finales). Eso contradice la garantía del README:43-52 ("A quick run can never overwrite ...").
   Arreglo: excluirlo de `make_all_figures` (ver punto 1) o implementar `--quick` con salida en `informe/*/quick/`.

## [IMPORTANTE]

3. **Archivos `--quick` versionados pese a `.gitignore`:** `informe/data/quick/background_bias.json`,
   `informe/figures/quick/fig_background_bias.{pdf,png}` (entraron en 07c84b5, antes de la regla de ignore).
   Si alguien corre `fig_10 --quick`, se ensucia `git status`. Además, un producto que no es final queda publicado.
   Arreglo: `git rm --cached -r informe/data/quick informe/figures/quick`. En `.gitignore` sobran las 3 líneas
   explícitas bajo `informe/data/quick/`, porque ya las cubre el directorio.

4. **README desactualizado (estructura):**
   - README.md:72-83: faltan `pminflux.py` (p-MINFLUX: cross-talk de lifetime, flickering secuencial vs intercalado) y
     `background.py` (fondo como parámetro libre: MLE (x, y, b), CRB 3×3).
   - README.md:71-111: el árbol no incluye `informe/` (informe.html, informe.pdf, build_pdf.js, figures/, data/).
   - README.md:87: "`fig_<n>_<name>.py` one script per figure -> paper/figures/fig<n>_<name>.pdf" es falso para fig_9 y
     fig_10 (escriben en `informe/`). Lo mismo pasa en la tabla Pipeline, README.md:117 (solo fig_1..fig_8, bien), pero
     `make_all_figures` sí los incluye (punto 1).
   - README.md:106: solo lista `equipo/2026-09-26_donut-localization/`. También existen `equipo/2026-09-28_sintesis-pdf/`
     (que se cita en :67) y `equipo/2026-09-28_final-review/`. Propuesta: "equipo/<fecha>_<id>/ the audit trails (one per job)".
   - README.md:13-14: "Tested on Windows with Python 3.8, numpy 1.24 ...". Es correcto, pero conviene agregar que la suite
     también pasa con Python 3.11 / numpy 2.4.6 / scipy 1.17.1 (verificado hoy), porque requirements no tiene cota superior.
   - README.md:25: "(about a minute)". Aquí tardó 22 s: es aceptable, no es un error.
   - Correctos: "5 pages" (informe.pdf tiene 5 páginas), y `docs/literature/` está descrito (:100).

## [MENOR] (validación de entradas, con escenario)

5. `src/donutloc/pminflux.py:95`: `crosstalk_matrix(np.inf)` pasa `tau >= 0` y devuelve una matriz de **NaN** (q=1 → 0/0).
   Arreglo: `if not (np.isfinite(tau) and tau >= 0)`.
6. `src/donutloc/pminflux.py:109-126`: `crosstalk_mc` no valida `T` ni `tau`. Con `T=0`, `t/T=inf`, el cast a int64 da basura y
   devuelve en silencio una matriz con toda la masa en la fila 0 (lo reproduje). Arreglo: la misma validación que `crosstalk_matrix`.
7. `src/donutloc/pminflux.py:178`: `S = 8 + 4 t_total/min(t_on,t_off)` no tiene cota de memoria. Por ejemplo, `t_total=400, t_on=0.01, n=1e4`
   da un array (1e4, 160008) de ~13 GB → MemoryError. Además, el `while` redibuja todo el lote cuando falla una sola
   trayectoria (condiciona a < S cambios). El sesgo es despreciable con el S elegido (P(Poisson(4) ≥ 24) ~ 1e-10), pero no está
   documentado.
8. `src/donutloc/pminflux.py:265-287`: `simulate_flicker_counts` no valida `p.shape == (K,)` ni `n_mean >= 0`
   (`n_mean<0` → ValueError de numpy poco claro). Las filas con w=0 y `fixed_N` dan conteos 0 (documentado, OK).
9. `src/donutloc/background.py:65`: `np.any(b < 0)` deja pasar `b = NaN` y devuelve p = NaN en silencio.
   Arreglo: `if not np.all(b >= 0)`.
10. `src/donutloc/background.py:205-215`: `fisher_matrix_bg` / `crb_free_bg` con b=0 en un cero exacto de una dona
    (p. ej. r=(0,0) con el TCP: p_centro=0) divide por cero (RuntimeWarning) y da F con inf/NaN. No está documentado:
    la CRB en la frontera b=0 no es regular. Arreglo: documentarlo o lanzar ValueError si `p.min() == 0`.
11. `src/donutloc/background.py:170-176`: el refinamiento de `mle_free_bg` crea `(M,125,3)` y `(M,125,K)` **sin** trocear,
    y `mem_budget` solo aplica a la grilla global. Con M=1e6 son ~4 GB. fig_10 usa M=4000, así que no le afecta hoy.
12. `scripts/fig_10_background_bias.py:196`: `open(dpath, "w")` sin `encoding="utf-8", newline="\n"`, a diferencia de
    los otros scripts. `.gitattributes` (`eol=lf`) normaliza al hacer commit y `ensure_ascii=True`, así que no rompe nada.
    Es solo una inconsistencia.

## Cobertura de tests

`tests/test_pminflux.py` (15 tests) y `tests/test_background.py` (15 tests) cubren la forma cerrada frente a MC, los límites,
el orden, la invariancia del fondo uniforme, el conteo exacto de pulsos (fuerza bruta), el MLE libre (lote = individual,
ceros, disco, b≥0, MC cerca de la CRB) y la validación. **No hay tests** para los casos 5-10 ni para los scripts
fig_9/fig_10: ningún test detecta el fallo del punto 1.

## Coherencia del repo

- `requirements.txt` / `pyproject.toml`: `numpy>=1.20` sin cota superior. Con numpy 2.4.6 **todo pasa** (254 OK), así que no hay rotura hoy.
- `git ls-files | grep __pycache__|\.pyc`: vacío. Hay archivos `quick/` versionados (punto 3), además de
  `equipo/.../r05/pre_quick.sha` (es del ledger, OK).
- `git ls-files '*.pdf'`: **no hay PDFs de terceros.** Solo productos propios (paper/, informe/) y dos PDFs de
  compilación del ledger (`equipo/2026-09-26.../work/review/r03/tex/{main,revtex_main}.pdf`), que son propios, no papers.
  Tampoco hay rutas `private` ni `p-minflux-main` versionadas.

## informe/informe.html

Hay tres `src=` (figures/fig_background_bias.png, figures/fig_pminflux_timing.png, figures/repo_fig1_schematic.png) y 10
rutas de repo citadas en el texto. **Todas existen.** "254 tests, 276 números" coincide con el estado real (254 tests;
`data/paper_numbers.json` tiene 276 claves).

## Contra el intent

El intent pide un explorador standalone guiado y una guía HTML. En el árbol no hay ningún `explorador*` ni carpeta nueva
de scripts (el `git status` solo muestra `equipo/2026-09-28_final-review/`). **Ese entregable todavía no existe.**
Esa parte es del worker; la reporto como faltante.
