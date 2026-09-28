# r01 — worker-C: sesgo por fondo en TCP fijo (análogo TCP de SimuFLUX Fig. 2d)

## Qué había antes (para no duplicar)
- `photons.probabilities` ya soporta fondo con SBR fija (Eq. S30) o `bg_per_exposure` (Eq. S28);
  `estimators.mle` acepta cualquier `p_fn`. El paper (`estimators.tex`) solo estudia MLE con el
  **modelo correcto** (SBR=10 conocido) y MLE sin fondo **sin fondo en los datos**. No había:
  error de modelo de fondo (datos con fondo, estimador sin fondo o con b erróneo), ni fondo como
  parámetro libre, ni CRB con b como parámetro molesto. Eso es lo nuevo.

## Qué hice
- `src/donutloc/background.py` (nuevo): `bg_from_sbr`, `probabilities_bg` (θ=(x,y,b)),
  `mle_free_bg` (MLE vectorizada en (x,y,b), b=u²≥0, grilla global + pattern search 3-D),
  `fisher_matrix_bg` (Fisher 3×3 por diferencias finitas) y `crb_free_bg` (bloque xy de F⁻¹;
  también devuelve el CRB con b conocido). Solo agregué una línea al docstring de `__init__.py`.
- `tests/test_background.py`: 15 tests (1.7 s). Cubren: ida y vuelta SBR↔b contra
  `photons.sbr_at`; `probabilities_bg` == `photons.probabilities(bg_per_exposure=b)`; bloque xy de
  la Fisher 3×3 == `fisher.fisher_matrix` del modelo con b conocido; CRB libre ≥ CRB conocido;
  desacople exacto F_xb=F_yb=0 en el centro del TCP; escala 1/√N; b=0 (diferencia hacia
  adelante); recuperación sin ruido de (x,y,b); con b̂>0 el ajuste reproduce exactamente n/N
  (4 exposiciones, 3 parámetros); log-verosimilitud del ajuste libre ≥ la del MLE con b conocido;
  lote == individual; cuentas nulas → (center, 0); estimaciones dentro del disco y b̂≥0; MC
  (N=500, SBR₀=5, x=20) sin sesgo (<3.5 SE) y σ dentro de 10 % del CRB libre; validaciones.
- `scripts/fig_10_background_bias.py` → `informe/figures/fig_background_bias.{png,pdf}` (17 cm,
  2 paneles, español; png dpi 200) y `informe/data/background_bias.json` (sort_keys).
  Corre en ~68 s (4000 rep/punto). `--quick` escribe en subcarpetas `quick/`.

## Convención
TCP fijo L=100 nm (`patterns.tcp_centers`, vértice arriba, centro último), dona LG fwhm=300 nm,
N=500 fotones **totales** (señal+fondo), multinomial condicionada a N. Fondo b constante por
exposición (Eq. S28), igual en las 4 ventanas; b fijado para SBR₀ = Σλ_j(0)/(4b) = 5 o 20 en
el centro: b = 0.029075 (SBR₀=5) y 0.0072688 (SBR₀=20) en unidades del pico del anillo. La SBR
real crece con x: en x=50 nm es 10.27 (SBR₀=5) y 41.06 (SBR₀=20) (`sbr_at_x`). Posiciones
(x,0), x = 0, 5, 10, 20, 30, 40, 50 nm. Mismas cuentas para los 5 estimadores
(`default_rng([42, i_sbr, i_x])`). Disco de búsqueda radio 150 nm. RMSE = sqrt(mean|e|²/2),
σ = sqrt((var_x+var_y)/2) (convención del CRB del repo). `bias_pop_*` = MLE sobre N·p (sin ruido).

## Resultados (claves en `informe/data/background_bias.json`)
MLE **sin fondo** (sesgo hacia afuera: el fondo eleva la fracción central y el modelo sin fondo
la interpreta como mayor r²):
- SBR₀=5: bias_x = +5.75±0.05 nm en x=5, +4.76 en x=10, +1.95 en x=30, +5.44±0.12 en x=50;
  bias_y ≈ +2.0…+3.6 nm (el eje x no es eje de simetría del TCP). RMSE 7.93 nm en x=0 (CRB 2.01;
  en x=0 el sesgo medio es 0 por simetría pero las estimaciones caen en un anillo) y 3.42–7.40
  en x≥5 frente a CRB 2.02–5.16. `sbr5.sin_fondo.*`.
- SBR₀=20: bias_x = +2.20±0.03 nm en x=5, +1.35 en x=50; RMSE 4.22 nm en x=0 (CRB 1.79).
- Sesgo poblacional concuerda con el MC (p. ej. SBR₀=5, x=10: pop 4.92 vs MC 4.76).
MLE **fondo exacto**: sesgo poblacional 0 (≤3e-5 nm, resolución del optimizador); σ/CRB = 1.004
(x=0) a 1.037 (x=50) con SBR₀=5. Sesgo MC de muestra finita pequeño pero significativo en x=50:
bias_y = −0.63±0.05 (SBR₀=5), −0.58±0.05 (SBR₀=20).
MLE **fondo libre (x,y,b)**: sesgo poblacional 0; σ = 2.02/2.09/2.40/5.74 nm en x=0/10/20/50
(SBR₀=5) vs CRB libre 2.01/2.08/2.36/6.17. El CRB libre supera al conocido en solo 0–4 % hasta
x=20 y 20 % (SBR₀=5) / 26 % (SBR₀=20) en x=50; en x=0 son iguales (desacople por simetría).
Con SBR₀=20 b̂ cae en la frontera b=0 en 7 % (x=10) a 30 % (x=50) de los casos
(`fondo_libre_frac_b_zero`), lo que da b̂ medio sesgado (0.0113 vs 0.0073 en x=50) y sesgo de
posición −0.60±0.09 (x) y −1.16±0.06 nm (y) en x=50.
MLE **b ±30 %**: sesgo |bias| ≤ 1.6 nm en todo el rango (SBR₀=5: −30 % da bias_x hasta +1.54,
+30 % da bias_y hasta −1.63 en x=50); RMSE dentro de ~10 % del CRB. Con +30 % el RMSE en x=0
queda por debajo del CRB (1.94 vs 2.01; 1.63 vs 1.79): estimador sesgado (encoge), no contradice
el CRB.

## Límites / dudas (no sobre-afirmar)
- El CRB libre supone b sin restricción; con b≥0 activo (SBR alta, x grande) el σ MC puede quedar
  por debajo del CRB libre (SBR₀=5, x=50: 5.74 < 6.17). No es un bound válido ahí.
- En x=0 el MLE sin fondo tiene máximos degenerados (3 simétricos); `bias_pop` ahí (−10.4 nm en y
  con SBR₀=5) es el desempate de la grilla, no un sesgo físico. No citarlo.
- Análogo a SimuFLUX Fig. 2d solo cualitativamente: ellos usan estimadores LSQ y 200+200 fotones;
  aquí es MLE, N=500, SBR₀ 5 y 20. No verifiqué su figura numéricamente.
- Suite completa: `python3 -m unittest discover -s tests` → 254 tests OK (skipped=1) en 25.7 s
  (incluye tests de otro worker, p. ej. test_pminflux). No identifiqué cuál test se saltea (el
  Bash dejó de responder al final); no es de `test_background.py` (15/15 OK, ninguno se saltea).
