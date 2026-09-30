# Ronda 6 — Consolidación (worker)

Fecha: 2026-09-30. Rama `consolidate-r06`. Sin commit (lo hace el orquestador).

Contexto: el trabajo cerró en R5 (f411211). Después, el commit 7d41d4d agregó la comparación con
SimuFLUX (secciones `bgphys_*`, `naive_*`, `simuflux_*` de `scripts/compute_paper_numbers.py`, 4
afirmaciones nuevas en `structure/claims.json` con `verified_round` 6) sin dejar rastro en `equipo/`.
Esta pasada deja el rastro y consolida el estado. La verificación independiente de las 4
afirmaciones de r06 la hace otro agente, en `reports/r06-verifier.md`; aquí **no** se re-verifican.

Decisión de la autora (inbox, ronda 6): consolidar. Las tres decisiones que el orquestador tomó en
su nombre quedan como preguntas abiertas.

## 1. Procedencia: 20 entradas huérfanas, ahora 0 warnings

Antes: `provenance: 222 tag(s), 242 entry(ies), 0 error(s), 20 warning(s)`.
Después: `provenance: 242 tag(s), 242 entry(ies), 0 error(s), 0 warning(s)`; `claims: 43, numbers:
276, provenance entries: 242`; `all checks pass`.

No cambió ningún valor numérico. No se tocaron `data/paper_numbers.json`,
`paper/generated/numbers.tex` ni `paper/provenance.json`. Las 20 entradas quedaron citadas
(opción (a)); no hizo falta borrar ninguna.

| Entradas | Decisión | Dónde |
|---|---|---|
| 16 `bgphys_{fixed,phys}_*` (σ/CRB de LMS x25/x50, sesgo MLE x25/x50, σ/CRB MLE x50, sesgo mLMS x50, σ/CRB mLMS x25/x50) | (a) `\src` (invisible) detrás de la frase "the biases by less than 1.3 nm and σ/CRB by less than 8 % at the points compared". Esas cotas se derivan justamente de estas entradas. | `paper/sections/estimators.tex` |
| `bgphys_iter_matchL150_sbr_first` (=10.0) | (a) `\src` junto al "$10$" de "if it is 10 in the first iteration". El 10 es la condición de diseño; la entrada documenta que la corrida lo cumple. | `paper/sections/iterative.tex` |
| `bgphys_iter_matchL25_sbr_last` (=10.0) | (a) igual, en "so that the SBR is 10 in the last iteration". | `paper/sections/iterative.tex` |
| `naive_eps0p002_crb_x10_nm`, `naive_eps0p002_crb_x20_nm` | (a) se citan con `\pnum{k}\src{k}` junto a los CRB de ε=0.01 ("… and 2.401 and 4.461 nm (ε=0.002)"). El párrafo ya compara el sesgo con el CRB para ε=0.002 (N de cruce), así que gana precisión. | `paper/sections/nonidealities.tex` |

Chequeo de las cotas antes de etiquetar (sobre `data/paper_numbers.json`, fixed → phys): el máximo
|Δ sesgo| es 1.222 nm (mLMS en x0=50), menor que 1.3. El máximo cambio relativo de σ/CRB es 7.13 %
(LMS en x0=50), menor que 8 %.

Nota para el verificador de r06: las cotas "1.3 nm" y "8 %" están tipeadas a mano en
`estimators.tex` y en `discussion.tex`. Methods dice "no computed result is typed by hand". Son
cotas derivadas, no parámetros, así que es una desviación menor de la regla. Ahora están
respaldadas por `\src`, pero no son números del registro. Queda en el backlog, porque convertirlas
exige agregar claves al registro.

## 2. Preguntas abiertas en `paper/sections/open_points.tex`

- (i) SBR=10 para `mle_efficiency_center`: ya estaba. Sin cambios.
- (ii) Nuevo ítem "Headline number for the misalignment bias". Dice que la pendiente poblacional
  sin ruido (`\pnum{misalignment_naive_bias_over_delta_pop_center}\src{…}`, 0.777) es el número
  principal y que el MC queda como complemento. Esa elección la hizo el coordinador en nombre de
  la autora y falta su confirmación.
- (iii) Nuevo ítem "Scope of Fig. 8". Dice que la figura se aclaró como error de posicionamiento
  del patrón (posiciones de los ceros, forma de la dona intacta) y no como desalineación de la
  máscara de fase. También es una decisión del coordinador pendiente de confirmación.

No se tipeó ningún número a mano. Las tres preguntas también figuran en el `backlog` de `state.json`.

## 3. PDF (Tectonic)

Compilé `paper/main.tex` con el tectonic.exe del scratchpad, el mismo comando que el paso LaTeX de
`reproduce.py`. No regeneré figuras ni el Monte Carlo.

- Exit 0, 0 Overfull. Hay algunos Underfull (cosméticos) y el aviso habitual de revtex.
- `paper/main.pdf`: **14 páginas** (en R5 eran 13).
- No hay "??" en el texto extraído.
- No hay referencias ni citas indefinidas en `main.log`. BibTeX da solo el warning habitual de
  `jnrlst`.

Rendericé y miré las páginas 6 (estimators: fondo constante por exposición), 8 (Fig. 5 + iterative
SimuFLUX), 10 (estimadores mal especificados + discusión), 11 (Fig. 7/8 + discusión SimuFLUX) y 12
(open points).

- Encontré y corregí un defecto que yo mismo había introducido: "compared(same seed" (faltaba un
  espacio por el `%` de fin de línea). Después de recompilar, el texto dice "points compared (same
  seed".
- Visto en p. 10 (menor, anterior a esta ronda, no corregido): en `nonidealities.tex` los números
  negativos en modo texto salen con guion, no con signo menos ("-1.463", "-0.3411"). Queda en el
  backlog.

## 4. Ledger

- Agregué `resolution` a las 62 afirmaciones con status refuted (27), unclear (21), claimed:pending
  (12) y claimed:unclear (2). No borré ninguna.
- Comprobé las afirmaciones de código ejecutando los escenarios originales contra el código actual:
  - `crb_map` sin salto: [1.6191, 1.6086, 1.6051, 1.6086, 1.6191];
  - `mle` con tol=0, tol<0 o grid_step=0 → ValueError; conteos negativos → ValueError;
  - `crb_tcp_center_limit(power=1.1)` emite un UserWarning;
  - campo lejano → p=NaN y crb=NaN;
  - `crb(..., N=array([100,400]), zero_policy='limit')` = [1.6051, 0.8025].
- Las demás las comprobé con grep sobre el código, los tests, el tex, `figures.json`, el README y
  `references.bib`.
- Resultado: 55 resueltas, 3 parciales y 4 abiertas.
  - Abiertas: (11) la normalización del rmse frente a Masullo, que es un open point del paper; (13)
    y (52) el SBR=10, pregunta a la autora; (143) la coincidencia "mismos 400 patrones", que no se
    reprodujo y no se cita.
  - Parciales:
    - (114) mitad verificada, mitad = (143);
    - (137) PDFs de figuras no deterministas byte a byte (CreationDate);
    - (138) `references.bib` sin pages/doi para Masullo2022 ni volume/pages/doi para Caprile2022.
- También precisé con el commit f411211 las `resolution` que ya existían en 159, 163 y 168, después
  de comprobarlas en el texto actual.
- `round` 5 → 6. `status` sigue "done". Nuevo campo `status_note`: la verificación independiente
  de r06 está en `reports/r06-verifier.md`.
- Agregué dos entradas de history: la ampliación SimuFLUX (commit 7d41d4d) y esta consolidación.
  Actualicé `last_checks` y `last_acceptance`.
- `inbox.jsonl`: agregué la entrada de la ronda 6 con la decisión de la autora de hoy.
- El formato JSON no cambió: indent=1, ensure_ascii=False, sin salto de línea final (un
  round-trip previo dio idéntico).

## 5. README

Dos cambios mínimos:

- La estructura de `docs/literature/` ahora lista A_, B_, C_insilico_minflux.md y
  C_insilico_vs_donutloc.md.
- Agregué un párrafo sobre la comparación con SimuFLUX (fondo constante por exposición,
  estimadores mal especificados, convenciones) y el alcance de la Fig. 8.

No encontré conteos desactualizados en el README.

## 6. Chequeos (resultados reales)

- `python -m unittest discover -s tests` → Ran 224 tests, OK (skipped=1).
- `python -m unittest tests.test_acceptance` → Ran 8 tests, OK.
- `python scripts/check_provenance.py` → 242 tags, 242 entries, 0 errors, 0 warnings; all checks
  pass.

Archivos cambiados: `README.md`, `paper/main.pdf`, `paper/sections/{estimators,iterative,
nonidealities,open_points}.tex`, `equipo/2026-09-26_donut-localization/{state.json,inbox.jsonl}` y
este reporte.

---

# Segunda pasada: correcciones del verificador r06

Fuente: `reports/r06-verifier.md`. Resultado: 4 verificadas y 5 unclear. Ese reporte no trae un
bloque `claims`, así que los registré en `state.json` como 9 afirmaciones de ronda 6 con source
"verifier". Las 5 unclear llevan `resolution`: 4 quedaron resueltas y 1 sigue abierta, porque el
código de SimuFLUX no está en esta máquina.

## Cambios

1. **Precisión.**
   - `scripts/compute_paper_numbers.py` suma `_polish_mle`, un pulido de Newton del
     log-likelihood con gradiente y hessiano por diferencias centrales, h=1e-3. Devuelve el punto
     de partida si sale del disco, si el hessiano no es definido negativo o si no converge.
     `_mle_noisefree` es `estimators.mle` seguido de ese pulido. Lo usan todos los MLE sin ruido de
     `naive_*`.
   - Tests nuevos: `tests/test_paper_tooling.py::TestNoiseFreePolish`, 4 tests.
     - El MLE correcto recupera la verdad a menos de 1e-6 nm.
     - El sesgo ingenuo con eps=0.002 en x0=20 vale 0.4220555 ± 2e-6, el valor independiente del
       verificador.
     - El pulido nunca baja el likelihood.
     - Si sale del disco o no converge, devuelve el punto de partida.
   - La clave del MLE honesto se describe ahora como un cero numérico. El texto dice "to numerical
     precision (below 10^-4 nm)" y ya no imprime el valor.
2. **nonidealities.tex / discussion.tex.** Los valores sin ruido se presentan como el límite de
   N grande del sesgo, que no se anula al crecer N y es mayor a N finito. No se imprimen los MC
   del verificador. También cambió la descripción de `naive_*_mle_N_bias_eq_crb_*`.
3. **Cotas 1.3 nm / 8 %.** Ahora son claves nuevas del registro, calculadas en
   `physical_background_section` sobre MLE, LMS y mLMS en x0 = 25 y 50 nm:
   - `bgphys_max_abs_change_bias_x_nm` = 1.2216 (se imprime 1.222);
   - `bgphys_max_rel_change_sigma_over_crb_pct` = 7.1272 (se imprime 7.127).

   Se citan con `\pnum\src` en `estimators.tex` y `discussion.tex`. La discusión agrega las
   condiciones (x0 = 25 y 50 nm, L=50 nm, N=100) y aclara que solo se comparó la componente x del
   sesgo. Hay 2 entradas nuevas en `paper/provenance.json`, y ambas claves se sumaron a la
   afirmación `bgphys_estimators_offcentre` de `structure/claims.json`.
4. **iterative.tex.**
   - Un SBR fijo es una idealización.
   - Subir la potencia solo ayuda si el fondo no escala con ella.
   - Una rampa moderada "only partially compensates" la caída.
   - La atribución a SimuFLUX se suavizó: "according to our reading of its public code, which we
     have not run".
   - `claims.json` tiene una advertencia nueva en `bgphys_iterative`.
5. **claims.json, `naive_misspecified_estimators`.** La advertencia sobre eps=0.05 ahora se
   restringe a x0=20: el MLE llega al borde del disco, |r|=37.5 nm; en x0=10 no llega. También
   corregí el enunciado ("límite de N grande", "< 1e-4 nm").
6. **Signos menos.** Los 4 `\pnum` negativos que estaban en modo texto
   (`naive_*_lms_extra_bias_x_*`) quedaron entre `$…$`. Un barrido de todos los `\pnum` negativos
   fuera de modo matemático ahora da 0.

## Valores que cambiaron (regeneración completa con `compute_paper_numbers.py`, 454 s)

Hay 2 claves nuevas (arriba) y ninguna eliminada. Cambiaron 13 valores, todos `naive_*`; fuera de
`naive_*` no cambió ningún valor, SE ni unidad. También cambiaron 17 descripciones, todas
`naive_*`, y 17 enunciados de procedencia se sincronizaron con ellas.

| clave | antes | después | impreso antes → después |
|---|---|---|---|
| naive_eps0p002_mle_bias_abs_x10_nm | 0.2687176 | 0.2687116 | 0.2687 → 0.2687 |
| naive_eps0p002_mle_bias_abs_x20_nm | 0.4220012 | 0.4220559 | **0.4220 → 0.4221** |
| naive_eps0p002_mle_N_bias_eq_crb_x10 | 7984.514 | 7984.873 | 7985 → 7985 |
| naive_eps0p002_mle_N_bias_eq_crb_x20 | 11173.39 | 11170.50 | **11173 → 11171** |
| naive_eps0p01_mle_bias_abs_x10_nm | 1.3887139 | 1.3886896 | 1.389 → 1.389 |
| naive_eps0p01_mle_bias_abs_x20_nm | 2.2894023 | 2.2893648 | 2.289 → 2.289 |
| naive_eps0p01_mle_N_bias_eq_crb_x10 | 429.8748 | 429.8898 | 429.9 → 429.9 |
| naive_eps0p01_mle_N_bias_eq_crb_x20 | 515.9977 | 516.0146 | 516.0 → 516.0 |
| naive_sbr20_mle_bias_abs_x10_nm | 0.3579192 | 0.3579223 | 0.3579 → 0.3579 |
| naive_sbr20_mle_bias_abs_x20_nm | 0.7396392 | 0.7396053 | 0.7396 → 0.7396 |
| naive_sbrinf_mle_bias_abs_x10_nm | 0.9157147 | 0.9157233 | 0.9157 → 0.9157 |
| naive_sbrinf_mle_bias_abs_x20_nm | 1.6428069 | 1.6428598 | 1.643 → 1.643 |
| naive_honest_mle_max_abs_bias_nm | 3.05e-5 | 5.70e-8 | ya no se imprime ("below 10^-4 nm") |

Los valores nuevos coinciden con el cálculo independiente del verificador (`deterministic.json`)
a menos de 1e-6 relativo. Por ejemplo, 0.4220559 contra 0.4220555, y N_eq 11170.50 contra
11170.52. En `numbers.tex` cambiaron 3 líneas y se agregaron 2.

## Chequeos después de la segunda pasada

- `python -m unittest discover -s tests` → Ran 228 tests, OK (skipped=1).
- `python -m unittest tests.test_acceptance` → Ran 8, OK. sha256 7d198853…2303 sin cambios.
- `python scripts/check_provenance.py` → 244 tags, 244 entries, 0 errors, 0 warnings; claims 43,
  numbers 278; all checks pass.
- Tectonic: exit 0, 0 Overfull, sin indefinidas. **14 páginas**, sin "??". Revisé visualmente las
  pp. 6, 8, 10 y 11.

---

# Tercera pasada: hallazgos del code-reviewer r06

Fuente: `reports/r06-code-reviewer.md`. Las 8 afirmaciones quedan en `state.json` con ronda 6 y
source "code-reviewer": 4 verified y 4 unclear, las 4 con `resolution`.

1. **`_polish_mle`.**
   - Tiene un parámetro nuevo, `return_reason=False`. Con True devuelve `(punto, motivo)`, donde
     el motivo es converged, not_concave, ill_conditioned, left_disk, max_iter o ll_decreased. Los
     llamadores existentes no cambian.
   - Guardas nuevas:
     - número de condición del hessiano mayor que `POLISH_MAX_COND` = 1e10: devuelve el inicio;
     - log-verosimilitud final menor que la inicial: devuelve el inicio.
   - Tests por rama en `TestNoiseFreePolish`, que ahora tiene 7 tests:
     - left_disk: start (9.95,0), center (9.9,0), radius 0.06. El contraste con radius 0.5
       converge a (10,0) con error menor que 1e-6.
     - not_concave: arranque en (0.5,0).
     - max_iter real: max_iter=2 desde (5,0). Con el valor por defecto converge.
     - ill_conditioned: se fuerza con `POLISH_MAX_COND` ≈ 1.
     - retorno por defecto: solo el punto.
   - Queda sin test dedicado la rama ll_decreased.
   - En los 8 casos reales de eps el motivo es "converged". Los 21 números sin ruido no
     cambiaron.
2. **SE de las claves `bgphys_max_*`.** Se toma el SE del par que da el máximo:
   - Sesgo: SEs en cuadratura, 0.0446 nm (mLMS en x0=50).
   - Cociente: SEs relativos propagados, 0.806 puntos porcentuales (LMS en x0=50).

   Es conservador, porque ambas corridas usan la misma semilla. Ahora se citan con `\pnumse` como
   "1.22 ± 0.04 nm" y "(7.1 ± 0.8) %". El coordinador esperaba ±0.05; el valor calculado es
   0.0446. El texto dice "the largest change … at the points compared" en `estimators.tex` y en
   `discussion.tex`. También actualicé los enunciados de procedencia y la afirmación
   `bgphys_estimators_offcentre` en `claims.json`.
3. **Test nuevo, `TestBgphysMaxChanges`.** Recalcula valor y SE de las dos claves desde las 24
   claves fuente de `data/paper_numbers.json`, compara con el helper `bgphys_max_changes` y
   comprueba que hay 24 fuentes y todas tienen SE.
4. **Regeneración completa (438 s).** Comparado con el registro de la segunda pasada, solo
   cambiaron las 2 claves `bgphys_max_*`: se les agregó el SE y cambió su descripción. Ningún
   valor cambió. En `numbers.tex` cambiaron 3 líneas y se agregaron 4 (7 inserciones y 3
   borrados en el `git diff`).

## Chequeos

- `python -m unittest discover -s tests` → Ran 232, OK (skipped=1).
- `python -m unittest tests.test_acceptance` → Ran 8, OK. sha256 sin cambios.
- `python scripts/check_provenance.py` → 244/244, 0 errors, 0 warnings; claims 43, numbers 278;
  all checks pass.
- Tectonic: exit 0, 0 Overfull, 0 "undefined", **14 páginas**, sin "??". Revisé visualmente las
  pp. 6 y 11, que son las que cambiaron.
