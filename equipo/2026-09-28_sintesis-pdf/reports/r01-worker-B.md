# r01 — worker-B: efectos de temporización de p-MINFLUX

## Qué hice

- `src/donutloc/pminflux.py` (nuevo; no toqué otros módulos ni `__init__.py`):
  `crosstalk_matrix` (forma cerrada), `crosstalk_mc` (Monte Carlo independiente),
  `make_model_crosstalk` (p_obs = M p, contrato `p_fn` compatible con `fisher`/`estimators`),
  `telegraph_on_intervals`, `sequential_on_times`, `interleaved_on_pulses` (conteo exacto de
  pulsos "on"), `exposure_weights`, `simulate_flicker_counts` (Poisson o N fijo).
- `tests/test_pminflux.py`: 15 tests, 0.17 s (M cerrada vs MC con 5 SE para 2 órdenes y 3 τ;
  filas y columnas suman 1; τ=0 → identidad, τ→∞ → 1/4; entradas explícitas; desplazamiento
  cíclico del orden no cambia M; simetría espejo entre órdenes (0,1,2,3) y (0,2,1,3); MLE con M
  recupera la posición sin ruido; bloques secuenciales a mano; conteo de pulsos entrelazados vs
  fuerza bruta; entradas inválidas).
- `scripts/fig_9_pminflux_timing.py` → `informe/figures/fig_pminflux_timing.{png,pdf}` (17 cm,
  3 paneles, dpi 200, español) y `informe/data/pminflux_timing.json` (sort_keys, `params.seed` = 42).
  Tiempo de corrida ≈ 61 s.

## Modelo (supuestos explícitos)

- Cross-talk: ranura s del período (P = 50 ns, T = 12.5 ns) dispara la exposición `order[s]`;
  ventana TCSPC s = [sT, (s+1)T), alineada con el pulso, sin IRF/jitter, decaimiento
  monoexponencial, sin saturación. Con q = exp(−T/τ):
  M[i,j] = (1−q) q^((slot(i)−slot(j)) mod 4) / (1−q^4). Circulante en orden de ranura ⇒
  doblemente estocástica (un fondo uniforme en el tiempo queda uniforme).
- **Orden de pulsos**: solo importa el orden *cíclico* (un corrimiento cíclico da la misma M; test).
  Para el TCP cualquier orden cíclico se lleva a otro por una simetría del patrón (rotación 3-fold /
  espejo), así que el orden cambia la *orientación* del mapa de sesgo, no el CRB en el centro ni
  su promedio en el disco. Probé A = (0,1,2,3) [3 donas en sentido antihorario desde el vértice
  superior, centro último] y B = (0,2,1,3) [espejo x→−x]. CRB idéntico entre A y B (a 1e-12); el
  sesgo sobre +x difiere (B en +x = espejo de A en −x).
- Flickering: telegrafía on/off exponencial, t_on = t_off = 100 µs, arranque estacionario,
  dwell 400 µs, N medio 100, sin fondo, TCP L = 100 nm, fwhm = 300 nm, 6000 localizaciones por
  caso. Secuencial: 4r bloques consecutivos en orden (0,1,2,3). Entrelazado: 8000 pulsos por
  exposición; se cuentan exactamente los que caen en intervalos "on".
  Tres medidas de σ_fl: `pop` (MLE sobre conteos esperados, límite N→∞: σ_fl puro); `fixedN`
  (multinomial N = 100; σ_fl² = STD² − STD_ctrl², control = mismo MLE sin flickering); `poisson`
  (N variable, se descartan N = 0, control con los mismos N_i). STD por eje = sqrt((var_x+var_y)/2).

## Resultados (JSON `informe/data/pminflux_timing.json`)

### Flickering (`flicker.<pos>.<esquema>`)

σ_fl (nm), centro / x = 20 nm; σ_CRB(N=100, límite) = 3.363 / 4.130 nm.

| esquema | pop centro | fixedN centro | pop x20 | fixedN x20 |
|---|---|---|---|---|
| seq r=1  | 16.07 | 16.02 | 25.21 | 24.00 |
| seq r=2  | 12.37 | 12.35 | 19.45 | 19.18 |
| seq r=5  | 6.63 | 6.64 | 10.51 | 10.54 |
| seq r=10 | 3.99 | 3.98 | 6.08 | 6.25 |
| seq r=25 | 2.22 | 2.32 | 3.47 | 3.49 |
| entrelazado | 0.018 | 0 (σ_fl² = −0.13 ± 0.14 nm²) | 0.032 | 0 (−0.13 ± 0.44 nm²) |

- σ_fl secuencial ∝ ~1/sqrt(r) aprox. (16.1 → 2.2 nm de r = 1 a 25 en el centro, razón 7.2 vs
  sqrt(25) = 5); `pop` y `fixedN` coinciden dentro de ~0.3 nm ⇒ STD² = σ_fl² + STD_ctrl² se cumple
  (σ_fl independiente de N, como en SimuFLUX).
- Entrelazado: la dispersión relativa de los pesos entre las 4 exposiciones es 2.9e-4
  (`weight_rel_std`), vs 0.79 (r=1) y 0.073 (r=25); σ_fl(pop) = 0.018 / 0.032 nm. Con MC
  (N = 100) σ_fl² es compatible con 0 en ambas posiciones y ambos modelos de conteo (Poisson x20:
  1.5 ± 2.7 nm²). **Resultado esperado confirmado: el entrelazado hace σ_fl ≈ 0 (< 0.04 nm).**
- Poisson (N variable) da σ_fl algo menor que fixedN a r grande (p.ej. r=25 centro 1.68 vs 2.32),
  con SE mayor (σ_fl² SE ≈ 0.9–3 nm²).
- Advertencia: en el centro, STD_ctrl (MLE, N=100) = 2.72 nm < σ_CRB = 3.36 nm (el MLE está
  sesgado hacia el centro a N bajo). Por eso uso el control empírico, no el CRB, como referencia;
  `fixedN_sigma_fl_vs_crb_nm` (contra CRB) queda en el JSON pero es menos fiable en el centro
  (r=25: 1.10 vs 2.32 nm).

### Cross-talk (`crosstalk.order{A,B}_tau{1..5}`)

q = exp(−12.5/τ): 3.7e-6, 1.9e-3, 0.0155, 0.044, 0.082 para τ = 1..5 ns.

(a) CRB(con M, modelo correcto)/CRB(ideal), idéntico para A y B:

| τ (ns) | centro vs límite r→0 | centro vs valor puntual S27 | disco r ≤ L/2 (razón de medias) | disco (media de razones) |
|---|---|---|---|---|
| 1 | 1.1389 | 1.0000 | 1.0004 | 1.0004 |
| 2 | 1.1411 | 1.0019 | 1.0063 | 1.0060 |
| 3 | 1.1570 | 1.0159 | 1.0398 | 1.0388 |
| 4 | 1.1923 | 1.0469 | 1.1043 | 1.1025 |
| 5 | 1.2446 | 1.0928 | 1.1914 | 1.1886 |

Nota: el centro exacto sin fondo es discontinuo (cualquier piso en la exposición central lleva el
CRB en r=0 al valor puntual S27 = 38.31/√N nm, vs el límite 33.63/√N). La columna "vs S27" es la
comparación limpia del efecto de M; el salto 1.139 a τ = 1 ns es ese artefacto, no cross-talk.
El promedio en el disco (grid 2.5 nm, 1257 puntos) es la cifra robusta: +0.04 % a τ=1 ns,
+19 % a τ=5 ns.

(b) Sesgo del MLE que ignora M (modelo ideal), N = 500, 2000 repeticiones por x, x = 0..50 nm:
- Sesgo poblacional (sin ruido) |b|: máx en x ∈ [0, 50] = 0.05, 1.17, 3.41, 5.95, 8.45 nm (orden A)
  y 0.05, 1.17, 3.41, 6.69, 12.97 nm (orden B) para τ = 1..5. En x = 0 es igual para A y B en
  módulo (A: (−6.03, 5.92) nm a τ = 5; B: (+6.03, 5.92)).
- MC N = 500: el sesgo medio en x = 0 es menor que el poblacional (τ=5, A: (−0.17, 3.78) nm vs
  (−6.03, 5.92)), porque cerca del centro la dirección del desplazamiento la fija el ruido de las
  exposiciones periféricas; lejos del centro MC y poblacional coinciden (fig. panel c).
  Claves: `naive_pop_bias_{x,y}_nm`, `naive_mc_bias_nm`, `naive_mc_rmse_axis_nm`.

(c) MLE con M en el modelo:
- Sin ruido: error máximo < 4e-6 nm en todas las x, τ, órdenes (`withM_pop_max_abs_bias_nm`)
  ⇒ consistente.
- MC N = 500: max |sesgo| por componente 0.46–0.75 nm. **No es cero estadísticamente** (χ² 99–251
  / 22 contra 0), pero la referencia sin cross-talk (τ = 0, modelo ideal) tiene el mismo sesgo de
  N finito: max 0.455 nm, χ² = 253/22 (`crosstalk.ref_tau0`), concentrado en x ≥ 35 nm (borde del
  TCP). Diferencia withM − ref: max |z| entre 1.5 y 5.0 (χ² 13–54 / 22); para A τ=4 y B τ=5 el χ²
  (46.5, 54.2) indica diferencias reales pero sub-nm. Conclusión honesta: el MLE con M elimina el
  sesgo sistemático del cross-talk (de hasta 8–13 nm a ≤ 0.75 nm) y queda con el sesgo de N finito
  propio del MLE, del mismo orden que sin cross-talk; no es "exactamente insesgado".

## Checks corridos

- `python3 -m unittest tests.test_pminflux`: 15 tests OK (0.17 s).
- Suite completa: NO corrida (ver Pendiente).
- M cerrada vs MC (10⁶ fotones por exposición): |M − M̂| máx 1.3e-6 (τ=1), 3.3e-4 (τ=3),
  7.1e-4 (τ=5).

## Pendiente / no resuelto

- **No pude re-correr el script después de los últimos retoques cosméticos de la figura** (títulos
  "(a)/(b)/(c)" con loc="left", posición de la etiqueta σ_CRB, xlabels más cortos): el
  clasificador de permisos de Bash falló repetidamente. El PNG/PDF actual en `informe/figures/`
  es de la corrida anterior: mismos datos, pero las etiquetas de panel "(a)(b)(c)" se superponen
  con los títulos y "σ_CRB(N=100)" se superpone con la curva en el panel (a). Hay que correr
  `python3 scripts/fig_9_pminflux_timing.py` (~61 s) una vez más; el JSON sale idéntico (semilla 42,
  el cómputo no cambió).
- **No corrí la suite completa** (`python3 -m unittest discover -s tests`) por el mismo motivo;
  solo `tests.test_pminflux` (15 OK). No toqué ningún otro módulo ni test, así que no espero
  regresiones, pero no está comprobado.
- No modelé IRF, jitter, ni ventanas desalineadas, ni tiempos muertos del detector; tampoco
  fondo con flickering. Sesgo con fondo: es tarea de worker-C.
- El régimen de flickering más rápido (t_on ~ µs, comparable con 50 ns × pocos pulsos) no se
  exploró; allí el entrelazado sigue siendo exacto por construcción del conteo de pulsos.
