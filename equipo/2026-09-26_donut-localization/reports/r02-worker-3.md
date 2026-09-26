# Ronda 2 — Worker 3: drivers de simulación (R4, R5)

Archivos: `src/donutloc/experiments.py` (nuevo), `tests/test_experiments.py` (nuevo). Solo se usa la
API pública (`patterns`, `beams`, `photons.make_model`/`sbr_at`, `fisher.crb`/`crb_limit`,
`estimators.mle`/`lms`, `montecarlo.run_mc`). No se importan `camera` ni `vectorial`. No se usa
el default de `crb_map`. No se hizo commit. Nada de `docs/private/` se usó.

## API implementada
- `error_stats(est, r_true)` -> bias, bias_abs, bias_se, sigma = sqrt((var_x+var_y)/2), sigma_se = sigma/(2 sqrt R), rmse = sqrt(mean|e|²/2).
- `l_schedule(n_iter=4, L0=150, L_min=None(=25), rule="fixed"|"adaptive", N_k, beam, sbr, bg_per_exposure, kappa=6)`.
  - `fixed`: geométrica de L0 a L_min.
  - `adaptive`: L_{k+1} = min(L_k, max(L_min, κ σ_k)), con σ_k = `crb_limit` en el centro de la TCP k. κ = 6 hace que L/2 = 3σ, el cuantil radial 99 % de un error gaussiano 2D. Es una elección propia: las fuentes no traen receta (nota A §5).
- `iterative_minflux(N_total, L_schedule=None, n_iter=4, L0=150, L_min=None, photon_split="equal"|weights, beam=None, sbr=None, r0_spread=None(=L0/4), n_rep=500, seed=42, recenter=True, rule="fixed", kappa=6, bg_per_exposure=None, sigma_psf=100, search_margin=0.25, mle_kw=None)`.
  - Devuelve: sigma, sigma_se, bias, bias_abs, rmse, L, N_k, sigma_iter, rmse_iter, bias_iter, crb_center_iter, sbr_center, camera_sigma = σ_PSF/√N_total, ratio_to_camera, crb_all_photons_Lmin, r_true, estimates, params (incluye la semilla).
  - Los fotones sobrantes van a la última iteración. La estimación final usa solo la última iteración, sin fusionar las anteriores.
  - El MLE busca en un disco de radio 0.75 L_k alrededor del centro de la TCP.
  - Todo va en lote gracias a la invariancia por traslación del haz (`photons.intensities` usa beam(x−cx)): se simula en coordenadas relativas a la TCP.
- `iterative_vs_photons(N_list, sigma_psf=100, **kw)` -> N, sigma, sigma_se, rmse, camera_sigma, ratio_to_camera, slope y slope_rmse (ajuste log-log), runs, params.
- `crb_center(L, N=100, fwhm=300, eps=0, sbr=None, zero_model="gaussian", bg_per_exposure=None, return_method=False)`. Usa el punto `fisher.crb(p,0,N)` si eps > 0 o si hay fondo. Con eps = 0 y sin fondo usa `fisher.crb_limit`.
- `eps_L_sweep(eps_list, L_list, N=100, fwhm=300, sbr=None, zero_model="gaussian", fov_radius=None, bg_per_exposure=None, fov=False)` -> crb_center (n_eps × n_L), method, L_best y, opcionalmente, crb_fov_mean: el promedio pesado por área en un disco de radio L/4 (regla del punto medio polar, sin el origen).
- `optimal_L(eps, N=100, fwhm=300, sbr=None, zero_model, bg_per_exposure, L_bounds=(1, 1.5 fwhm), n_grid=60)` -> L_opt, crb_opt, interior, params. Hace un barrido log de 60 puntos y después `minimize_scalar` acotado en log L.
- `misalignment_study(displacement_list, L=100, N=500, fwhm=300, sbr=10, n_patterns=20, n_rep=200, seed=42, estimator="mle"|"lms", positions=None(=[(0,0),(L/4,0)]), beam=None, search_margin=0.25)`.
  - Números aleatorios comunes: las mismas direcciones de `perturb_centers` para todos los δ (`default_rng(seed)` recreado en cada δ) y las mismas semillas de `run_mc` por (patrón, posición). Los estimadores honesto e ingenuo ven exactamente los mismos conteos.
  - Devuelve, para "honest" y "naive": bias_abs_mean, bias_abs_se (std entre patrones/√P), bias_chi2 = media de |b|²/(se_x²+se_y²) (≈1 si no hay sesgo), sigma_mean, sigma_se, rmse. Además crb_honest (`crb_limit` en r_true), bias_noise_floor = σ√(π/(2R)) e identical.

## Tests
`python -m unittest discover -s tests -p "test_experiments.py"`:
```
.................
----------------------------------------------------------------------
Ran 17 tests in 2.608s

OK
```
Suite completa (con los archivos de W1/W2 tal como estaban en disco en ese momento): `Ran 128 tests in 4.605s — OK (skipped=1)`.

## Resultados (configuraciones de referencia; semilla 42; LG de Eq. S17 con fwhm = 300 nm; TCP con el centro al final; multinomial; SE = sigma/(2√R))

### (a) MINFLUX iterativo contra cámara
Configuración de referencia:
- N_total = 1000, 4 iteraciones con reparto igual (250 fotones por iteración).
- `rule="fixed"`: L = 150, 82.55, 45.43 y 25 nm.
- r0_spread = 37.5 nm (uniforme en disco), recentrado, MLE con search_margin = 0.25, n_rep = 500.
- σ_PSF = 100 nm (Balzarotti p. 1, p. 32).

| esquema | SBR | iterative_sigma_nm | SE | rmse | camera_sigma_nm | cociente |
|---|---|---|---|---|---|---|
| fixed | ∞ | 0.4966 | 0.0111 | 0.4962 | 3.1623 | 0.157 |
| fixed | 10 | 0.6132 | 0.0137 | 0.6127 | 3.1623 | 0.194 |
| adaptive (κ=6) | ∞ | 0.4828 | 0.0108 | 0.4823 | 3.1623 | 0.153 |
| adaptive (κ=6) | 10 | 0.5765 | 0.0129 | 0.5762 | 3.1623 | 0.182 |

- **Candidato a `paper_numbers`:** `iterative_sigma_nm` = 0.497 ± 0.011 frente a `camera_sigma_nm` = 3.162 (fixed, sin fondo). Con SBR = 10: 0.613 ± 0.014.
- σ por iteración (fixed, sin fondo) = 4.344, 1.798, 0.976 y 0.497 nm. El CRB límite en el centro por iteración da 3.470, 1.721, 0.920 y 0.502. La primera iteración queda arriba del CRB central porque el emisor está fuera del centro (hasta L/4).
- Si se gastaran los 1000 fotones en L = 25 (requiere conocer la posición de antemano), el CRB sería 0.251 nm (0.305 con SBR = 10).
- Adaptive sin fondo: L = 150, 25, 25, 25.
- Sin recentrado: σ por iteración = 4.34, 3.20, 3.23 y 8.06 nm. El zoom sin recentrar falla.
- Con fondo fijo por exposición λ_b = 0.0297 (SBR = 10 en L = 150): la SBR en el centro cae a 10.0, 3.42, 1.07 y 0.33 (Eq. S32), y σ por iteración = 4.83, 2.47, 2.07 y 2.26 nm. Achicar L hasta 25 nm empeora el resultado.

Curva en función de N (fixed, n_rep = 500):
- Sin fondo, N = 250, 500, 1000, 2000, 4000 y 8000: σ = 1.094, 0.725, 0.497, 0.362, 0.233 y 0.180 nm. Pendiente log-log −0.526. El cociente con la cámara queda en 0.15–0.17.
- SBR = 10: σ = 1.612, 0.910, 0.613, 0.414, 0.292 y 0.220 nm. Pendiente −0.568. El cociente queda en 0.18–0.25.
- La cámara tiene pendiente −0.5 por construcción. Como el esquema tiene L_min fijo, la pendiente asintótica es −0.5. La parte más empinada en N chico viene de que la primera iteración no alcanza.

### (b) Profundidad finita del cero
N = 100, fwhm = 300, `zero_model="gaussian"`.

CRB en el centro (nm); eps = 0 es el límite, eps > 0 el valor puntual:

| eps \ L | 10 | 25 | 50 | 75 | 100 | 150 | 200 | 300 |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.316 | 0.794 | 1.605 | 2.454 | 3.363 | 5.487 | 8.382 | 25.23 |
| 0.002 | 0.747 | 1.047 | 1.884 | 2.829 | 3.877 | 6.453 | 10.26 | 34.70 |
| 0.01 | 2.317 | 1.683 | 2.213 | 3.060 | 4.063 | 6.607 | 10.42 | 35.23 |
| 0.05 | 10.31 | 4.908 | 3.879 | 4.230 | 5.007 | 7.394 | 11.27 | 38.07 |
| 0.15 | 31.37 | 13.40 | 8.262 | 7.311 | 7.500 | 9.487 | 13.54 | 46.60 |

L óptimo (argmin continuo en [1, 450] nm):

| eps | L_opt (sin fondo) | CRB_opt | L_opt (SBR = 10) | CRB_opt |
|---|---|---|---|---|
| 0 | en la cota (1 nm; monótono) | 0.032 | en la cota (1 nm) | 0.038 |
| 0.002 | 10.48 | 0.746 | 10.55 | 0.815 |
| 0.01 | 23.28 | 1.679 | 23.41 | 1.835 |
| 0.05 | 50.34 | 3.878 | 50.60 | 4.242 |
| 0.15 | 80.81 | 7.287 | 81.17 | 7.974 |

- Escala empírica: L_opt/(fwhm √eps) = 0.78, 0.78, 0.75 y 0.70. Es aproximadamente L_opt ≈ 0.78 fwhm √eps para eps chico. Es una observación numérica, no está derivada.
- Una SBR = 10 fija (Eq. S30) casi no mueve L_opt.
- Con eps = 0 y SBR fija el CRB es monótono en L. Para eps = 0, el L óptimo por fondo requiere fondo fijo por exposición (Eq. S32), y eso está disponible con `bg_per_exposure`.
- Chequeos: con eps = 0.05 y L = 50, el punto (3.878508) coincide con `crb_limit` (3.878508). Con eps = 0, el CRB es ≈ √0.1·L/√N a 1 % para L ≤ 40, y en L = 50 da 1.6051.

### (c) Desalineación: estimador honesto contra ingenuo
Configuración:
- L = 100, N = 500, SBR = 10, fwhm = 300.
- 20 patrones × 200 reps, MLE.
- Cada uno de los 4 ceros se desplaza δ en una dirección aleatoria.

Valores en el centro (0,0) y en (25,0):

| δ (nm) | pos | honesto: mean\|b\| ± se | naive: mean\|b\| ± se | honesto: rmse | naive: rmse | CRB honesto | naive σ |
|---|---|---|---|---|---|---|---|
| 0 | c | 0.187 ± 0.024 | idéntico | 1.850 | idéntico | 1.863 | 1.847 |
| 2 | c | 0.189 ± 0.025 | 1.575 ± 0.137 | 1.847 | 2.198 | 1.855 | 1.851 |
| 5 | c | 0.169 ± 0.029 | 4.008 ± 0.352 | 1.846 | 3.570 | 1.831 | 1.884 |
| 10 | c | 0.176 ± 0.025 | 8.306 ± 0.684 | 1.838 | 6.549 | 1.822 | 1.981 |
| 0 | L/4 | 0.217 ± 0.027 | idéntico | 2.387 | idéntico | 2.363 | 2.384 |
| 2 | L/4 | 0.235 ± 0.025 | 1.507 ± 0.141 | 2.429 | 2.643 | 2.381 | 2.381 |
| 5 | L/4 | 0.231 ± 0.026 | 3.589 ± 0.376 | 2.509 | 3.681 | 2.438 | 2.387 |
| 10 | L/4 | 0.275 ± 0.047 | 7.067 ± 0.757 | 2.687 | 6.760 | 2.634 | 2.952 |

- Piso de ruido MC de mean|b| para un estimador insesgado: σ√(π/2R) ≈ 0.16 (centro) y 0.21–0.23 (L/4).
- Con δ = 0, honesto e ingenuo dan estimaciones idénticas bit a bit.
- El sesgo ingenuo es ≈ 0.8 δ en el centro y ≈ 0.7 δ en L/4. La dispersión cambia poco: la pérdida de precisión del ingenuo es casi toda sesgo.
- bias_chi2 honesto = 1.06–1.52 (esperado 1 ± ~0.22 por celda). El ingenuo da 46–2222.

### Tiempos (esta máquina, Python 3.8)
- `iterative_minflux` de referencia (n_rep = 500, 4 iteraciones): 0.37–0.40 s.
- `iterative_vs_photons` con 6 valores de N: 1.9–2.4 s.
- Tabla de `eps_L_sweep` (5 × 8, con fov): 0.06 s. `optimal_L`: 0.02–0.04 s cada uno (0.4 s en la primera llamada).
- `misalignment_study` de referencia (4 δ × 20 patrones × 2 posiciones × 200 reps × 2 estimadores): 11–12 s.
- `test_experiments.py` completo: 2.6 s.

## Qué no quedó resuelto
- **bias_chi2 del honesto.** Tiende a quedar algo por encima de 1: la media de las 8 celdas es ≈ 1.24, y las celdas están correlacionadas porque se usan números aleatorios comunes. Hay unos 40 pares (patrón, posición) independientes por δ, así que la desviación es de ~1.5 SE. Puede ser un sesgo O(1/N) real del MLE, de ~0.1 nm, o una fluctuación. Para decidirlo haría falta más n_rep.
- **Pendiente frente a N.** −0.53 y −0.57 no son una ley de escala derivada; dependen del esquema (L_min fijo) y del rango de N.
- **Iterativo.** No se fusionan las estimaciones de iteraciones anteriores ni se usa un prior, y la regla adaptativa se basa en el CRB central, no en el error real. Son elecciones propias y están documentadas en los docstrings.
