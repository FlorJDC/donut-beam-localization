# Ronda 3 — Worker 1: fixes, W11, pipeline de números, claims y procedencia

No hice commit. No leí `docs/private/`. El sha256 de `tests/test_acceptance.py` sigue siendo `7d198853e87c7bd94bb8be07a855651146e1d849310281546cd8089dcf262303`.

## Archivos
- **Nuevos:**
  - `scripts/_paperconfig.py` (literal del plan; después, por orden del orquestador, `MIS.n_patterns` pasó de 50 a 400);
  - `scripts/run_iterative_sweep.py`, `scripts/compute_paper_numbers.py`, `scripts/check_provenance.py`;
  - `data/iterative_sweep.json`, `data/paper_numbers.json` (204 claves), `structure/claims.json` (37 claims);
  - `paper/main.tex` (stub), `paper/provenance.json` (`{}`), `paper/generated/numbers.tex`, `paper/sections/.gitkeep`;
  - `tests/test_check_provenance.py`.
- **Modificados:**
  - `src/donutloc/montecarlo.py`, `experiments.py`, `fisher.py`, `vectorial.py`;
  - `tests/test_fisher.py`, `test_montecarlo.py`, `test_experiments.py`, `test_vectorial.py`.
- `data/mc/iterative_sweep_quick.json` (la salida de `--quick`) lo borré porque `.gitignore` no lo cubre.

## (a) Las 4 refutadas de r2

**1. `misalignment_study` sigma_se.**
- `sigma_se` ahora es std(ddof=1) entre patrones de las σ por patrón, dividida por √P. Con P=1 da NaN.
- La fórmula vieja queda como `sigma_se_within`.
- El docstring está actualizado.
- Tests:
  - escenario δ=10 en el centro, 20×200: `sigma_se > 1.5·sigma_within` para honesto e ingenuo;
  - un test recalcula la fórmula a mano (P=3) con la siembra documentada y coincide a 1e-10.

**2. `fisher.crb(..., zero_policy="limit")` con N en array.**
- Ahora hay broadcast conjunto de r y N (`np.broadcast_shapes`).
- `crb(p,[0,0],array([100,400]),"limit")` devuelve shape (2,) igual a `crb_limit` en cada N (1.605096 y 0.802548).
- También funciona el broadcast externo, puntos (2,) × N (3,1) → (3,2).

**3. Docstring de `l_schedule` adaptivo.**
- Ahora explica que σ_k es el CRB central y subestima el error real.
- Medido con los defaults (10⁴ emisores, semilla 42):

  | cantidad | valor |
  |---|---|
  | CRB central de la iteración 0 | 3.47 nm |
  | error real de la iteración 0 | 4.26 nm |
  | fracción fuera de L_1/2 = 12.5 nm (semilla 42) | 2.65 % |
  | fracción fuera de L_1/2 = 12.5 nm (semilla 7) | 2.72 % |
  | fracción fuera de 3σ_k solo | ~6 % |

- La cobertura es ~97 %, no 99 %.
- Test: con 3000 repeticiones, la fracción fuera cae en (1.5 %, 4.5 %).

**4. Beam vectorial lineal: se implementó la opción preferida.**
- Se tabulan los armónicos azimutales de I:
  - el campo tiene órdenes −1..3, así que I tiene |n| ≤ 4;
  - numéricamente solo n = 0 y ±2 son no nulos.
- Representación: `I = g0(ρ²) + 2 Re Σ g_n(ρ²)(x+iy)^n`, con g_n suaves.
- La interpolación es con **spline cúbico en ρ²**. Probé interpolación lineal en ρ² y no alcanzaba: en L=100, punto (7,3), daba −0.25 % con d_rho=1 y +0.11 % con 0.5, porque los gradientes quedan constantes a trozos.
- Resultado (n=1.518, N=100) contra `mode="exact"`:

  | L (nm) | CRB en (7,3) (nm) | crb_limit (nm) | diferencia relativa |
  |---|---|---|---|
  | 50 | 62.35509 | 60.77420 | ≤ 4e-9 |
  | 100 | 39.47879 | 38.06131 | ≤ 4e-9 |
  | 150 | 37.47657 | 35.72298 | ≤ 4e-9 |

  La intensidad coincide a ~1e-10.
- El interp cartesiano sigue disponible solo con `linear_method="cartesian"` y emite `UserWarning`.
- Costo: construir el beam con rho_max=1500 tarda ~21 s.
- Tests: CRB y crb_limit con n_theta=401 a <0.5 % de exact, más 62.36 y 60.77 ±0.1; el cartesiano avisa y sesga >2 %.

## (b) Bootstrap
- `montecarlo.bootstrap_sigma_se(errors, n_boot=1000, seed=42, return_samples=False)` con σ=√((var_x+var_y)/2), ddof=1:
  - descarta filas NaN;
  - trabaja en bloques;
  - lanza ValueError con forma incorrecta, R<2 o n_boot<2.
- También se agrega `montecarlo.sigma_of_errors`.
- `run_mc` devuelve `sigma_se_boot`; con `n_boot=0` da NaN.
- `experiments.error_stats`: `sigma_se` pasa a ser el bootstrap y el viejo queda como `sigma_se_gauss`. `misalignment_study` usa `n_boot=0` internamente.
- `iterative_minflux` devuelve además `sigma_se_gauss` y `sigma_se_iter`.
- Tests:
  - con datos gaussianos (R=500 y 4000) coincide con σ/(2√R) a ±10 %;
  - con t de Student (ν=4) da > 1.3× la fórmula gaussiana;
  - reproducible, invariante a traslaciones e independiente del tamaño de bloque.

## (c) W11: `data/iterative_sweep.json`
- **Protocolo:** `ITER`, sin fondo, semilla 42 en cada N, 10⁴ repeticiones, SE por bootstrap con N_BOOT=2000. El IC de la pendiente sale de remuestrear las repeticiones de cada N. Tardó 128 s.

| N | σ (nm) | SE boot | SE gauss | σ/cámara |
|---|---|---|---|---|
| 250 | 1.1580 | 0.0137 | 0.0058 | 0.1831 |
| 500 | 0.7359 | 0.0042 | 0.0037 | 0.1646 |
| 1000 | 0.5077 | 0.0027 | 0.0025 | 0.1606 |
| 2000 | 0.3544 | 0.0019 | 0.0018 | 0.1585 |
| 4000 | 0.2495 | 0.0013 | 0.0012 | 0.1578 |
| 8000 | 0.1763 | 0.0010 | 0.0009 | 0.1577 |

- **Pendiente log-log:**

  | rango | pendiente | IC 95 % |
  |---|---|---|
  | todo el rango | −0.5365 ± 0.0029 | [−0.5424, −0.5312] |
  | N ≥ 500 | −0.5148 ± 0.0025 | [−0.5196, −0.5096] |

- **Lectura honesta:** la pendiente **no es compatible con −0.5 al 95 %** en ninguno de los dos rangos. El iterativo mejora algo más rápido que 1/√N porque el cociente con la cámara cae de 0.183 a 0.158 y se aplana desde N≈2000. Con N ≥ 1000, un cálculo a mano da ≈ −0.51. El claim queda como "−0.515…−0.537 según el rango; cociente 0.158–0.183".
- **Coherencia con el verificador de r2:** mis valores caen dentro de su rango (1.180±0.020 en 250; 0.2514±0.0024 en 4000; pendientes −0.545 y −0.523). El 1.094 que dio el worker en N=250 y el 0.233 en N=4000 no se reproducen.
- **Extras en N=1000:**
  - SBR=10: 0.6165±0.0033;
  - adaptivo (L=150,25,25,25): 0.4795±0.0027;
  - `recenter=False`: 8.14±0.05, marcado `artefact: true` con nota.
- **Formato del JSON:** `N, sigma, sigma_se, sigma_se_gauss, rmse, camera_sigma, ratio_to_camera, per_N[...], slope_all, slope_N_ge_500, ratio_range, extras_N_ref{sbr10, adaptive, no_recenter}, params`. Cada entrada de `per_N` trae `sigma_iter, sigma_se_iter, L, N_k, crb_center_iter, frac_err_gt_5nm, max_err_nm, ratio_se`.

## (d) `compute_paper_numbers.py`
- Es el único escritor de `data/paper_numbers.json` y `paper/generated/numbers.tex`. Esta última define `\pnum{k}`, `\pnumse{k}` y `\csname pnum@k\endcsname`, sin marca de tiempo, así que la salida es determinista.
- Usa solo la API pública y `_paperconfig`.
- La corrida final, sin `--quick`, tardó 375 s; la desalineación con 400 patrones se lleva 338 s.
- **Puntos compartidos con las figuras:** mismos parámetros y semilla, y mismos nombres de clave que fig 1/2/3/5/6/7/8.
  - `mle_efficiency_center` = 0.99123 ± 0.00513 (σ 1.94287, CRB 1.96006) es **idéntico bit a bit** a `data/fig5_summary.json`.
  - También coinciden exactamente `mle_nobg_*`, `lms_nobg_sigma_over_s27` (N=100), `mle_nobg_bias_x_r2_nm` y `fig5_mle_sbr10_eff_max_dev_inside` (0.51163), este último igual a `mle_sbr10_sigma_over_crb_x15_N100` − 1.
  - Los números del iterativo salen de `iterative_sweep.json`, el mismo archivo que usa fig 6.
- **Desalineación:** con `MIS` = 400 patrones × 200. La fig 8 tiene que volver a correrse con `--no-cache` para coincidir. Mi pendiente del sesgo usa solo δ ∈ {2,5,10} y lleva una clave distinta, `misalignment_naive_bias_over_delta_*`. La fig 8 ajusta una grilla de δ más densa bajo `misalignment_naive_bias_slope_*`.

### Tabla clave → valor → SE (principales; la lista completa de 204 claves está en data/paper_numbers.json)

**Aceptación**

| clave | valor | SE |
|---|---|---|
| fwhm_nm | 300 | |
| crb_center_lg_L50_N100_nm | 1.605096 | |
| crb_exponent_L (L_FIT 5–40) | 1.00435 | |
| crb_exponent_N | −0.500000 | |
| mle_efficiency_center | 0.99123 | 0.00513 |
| vectorial_zero_depth_correct | 0.0 (cero de máquina) | |
| vectorial_zero_depth_wrong_handedness | 0.845241 | |
| iterative_sigma_nm | 0.50774 | 0.00268 |
| camera_sigma_nm (ideal, declarado) | 3.16228 | |

**CRB en el centro y escalamiento**

| clave | valor |
|---|---|
| crb_center_point_S27_L50_N100_nm | 1.802472 |
| crb_center_lg_L50_N100_closed_nm | 1.605096 |
| limit_to_point_ratio_L5 / L50 / L100 / L150 | 0.89439 / 0.89050 / 0.87805 / 0.85527 |
| limit_to_point_ratio_quadratic | 0.894427 |
| limit_to_point_ratio_quadratic_numeric | 0.894427 |
| crb_limit_axis_par / perp, L100 N100 cuadrático | 2.73861 / 3.53553 nm |
| crb_limit_axis_par / perp, L50 N100 | 1.37977 / 1.80247 nm |
| crb_center_sbr10_L50_N100_nm | 1.960059 |
| photons_for_5nm_L50_sbr inf / 50 / 20 / 10 / 5 | 13.00 / 13.45 / 14.16 / 15.37 / 17.93 |
| photons_for_5nm_L50_limit_nobg | 10.31 |
| crb_center_sbr5_N500_fwhm360_L50 / 100 / 150 | 0.94129 / 1.96237 / 3.16727 nm |
| crb_center_sbr5_N500_fwhm300_L50 / 100 / 150 | 0.94694 / 2.01241 / 3.37012 nm |
| crb_1d_quadratic_L50_N100_nm | 1.25 |
| multifotón c = 1 / 2 / 3 | 3.5355 / 1.7678 / 1.1785 nm |

**Cámara**

| clave | valor |
|---|---|
| camera_pixelated_9x9_N400_nm | 5.20455 |
| camera_pixelated_9x9_N1000_nm | 3.29165 |
| camera_perpixel_sbr500_N600_nm | 4.96418 |
| camera_photons_for_5nm_perpixel_sbr500 | 591.43 |
| camera_total_sbr500_N600_nm | 4.26772 (alternativa, SBR total) |
| camera_photons_for_5nm_total_sbr500 | 437.12 |
| camera_sbr_perpixel_to_total | 6.1728 |

**Dona vectorial**

| clave | valor |
|---|---|
| vectorial_zero_depth_linear | 0.371648 |
| vectorial_peak_to_peak_diameter_nm | 384.666 |
| vectorial_zero_curvature_nm2 | 7.04222e-5 |
| vectorial_lg_fwhm_curvature_nm | 327.141 |
| vectorial_lg_fwhm_diameter_nm | 320.256 |

| CRB en el centro (N=100) | L=50 | L=100 | L=150 |
|---|---|---|---|
| vectorial correcta (nm) | 1.60041 | 3.32324 | 5.33104 |
| vectorial opuesta (nm) | 147.886 | 89.588 | 84.197 |
| vectorial lineal (nm) | 60.774 | 38.061 | 35.723 |
| LG 300 (nm) | 1.60510 | 3.36341 | 5.48647 |
| LG igualada por curvatura (nm) | 1.60124 | 3.32964 | 5.35144 |
| vectorial / LG curvatura | 0.99949 | 0.99808 | 0.99619 |
| vectorial / LG 300 | 0.99708 | 0.98806 | 0.97167 |

**Estimadores**

| clave | valor | SE |
|---|---|---|
| mle_sigma_center_sbr10_nm | 1.94287 | 0.0101 |
| mle_nobg_sigma_over_crb_limit_center | 0.83970 | 0.00415 |
| mle_nobg_sigma_over_s27_center | 0.74775 | 0.00369 |
| mle_nobg_bias_x_r2_nm (40000 repeticiones) | −0.34445 | 0.0077 |
| mle_sbr10_sigma_over_crb_x15_N100 | 1.5116 | 0.029 |
| mle_sbr10_sigma_over_crb_x15_N1000 | 1.0042 | 0.0057 |
| lms_nobg_sigma_over_s27_center | 0.99998 | 0.0049 |
| lms_nobg_sigma_over_crb_limit_center | 1.12294 | 0.0055 |

**Profundidad finita del cero**

- V9, CRB(eps)/S27 (constante | gaussiano):

  | eps | L=50 | L=100 |
  |---|---|---|
  | 0.002 | 1.0454 \| 1.0455 | 1.0120 \| 1.0121 |
  | 0.01 | 1.2268 \| 1.2277 | 1.0602 \| 1.0606 |
  | 0.05 | 2.1300 \| 2.1518 | 1.3002 \| 1.3072 |
  | 0.15 | 4.3817 \| 4.5835 | 1.8985 \| 1.9580 |

- W12, pedestal gaussiano, N=100:

  | eps | L_opt (nm) | CRB(L_opt) (nm) | L_opt/(fwhm√eps) |
  |---|---|---|---|
  | 0.002 | 10.484 | 0.7458 | 0.7815 |
  | 0.01 | 23.280 | 1.6787 | 0.7760 |
  | 0.05 | 50.337 | 3.8784 | 0.7504 |
  | 0.15 | 80.810 | 7.2867 | 0.6955 |

- zero_depth_transition_scale_nm = 4.887; eps0_crb_divergence_L_nm = 360.34.

**Iterativo**

| clave | valor | SE |
|---|---|---|
| iterative_ratio_to_camera_N1000 | 0.16056 | 0.00085 |
| iterative_slope_all | −0.5365 | 0.0029 (IC [−0.5424, −0.5312]) |
| iterative_slope_N_ge_500 | −0.5148 | 0.0025 (IC [−0.5196, −0.5096]) |
| iterative_ratio_min / max | 0.1577 / 0.1831 | |
| iterative_sbr10_sigma_nm | 0.6165 | 0.0033 |
| iterative_adaptive_sigma_nm | 0.4795 | 0.0027 |
| iterative_no_recenter_sigma_nm_ARTEFACT | 8.145 | 0.045 |
| iterative_crb_all_photons_L25_N1000_nm | 0.2509 | |
| adaptive_frac_outside_next_radius | 0.0265 | 0.0016 |
| adaptive_iter0_sigma_nm | 4.2635 | 0.027 |
| adaptive_iter0_crb_center_nm | 3.4699 | |

**Desalineación (400 × 200), en el centro**

| δ (nm) | sesgo ingenuo (nm) | sesgo honesto (nm) | piso MC (nm) | σ ingenuo (nm) | σ honesto (nm) | rmse ingenuo (nm) | CRB honesto (nm) |
|---|---|---|---|---|---|---|---|
| 0 | 0.173±0.005 | 0.173±0.005 | 0.164 | 1.855±0.003 | 1.855±0.003 | 1.857 | 1.863 |
| 2 | 1.492±0.033 | 0.169±0.004 | 0.164 | 1.855±0.003 | 1.851±0.004 | 2.180 | 1.861 |
| 5 | 3.731±0.080 | 0.177±0.005 | 0.164 | 1.884±0.004 | 1.845±0.006 | 3.433 | 1.845 |
| 10 | 7.708±0.158 | 0.180±0.005 | 0.164 | 2.015±0.011 | 1.854±0.009 | 6.228 | 1.842 |

- Pendiente del sesgo del ingenuo por el origen (δ ≥ 2): 0.755 ± 0.009 en el centro y 0.862 ± 0.012 en (L/4,0).
- En (L/4,0), sesgo del ingenuo: 1.691, 4.202 y 9.304±0.258 para δ = 2, 5 y 10.

## (e) `structure/claims.json`: 37 claims
- Cubren R1–R5.
- `status` = `"verified"` para los resultados de r1/r2 cuyos números son deterministas y ya estaban verificados. Los números MC recalculados en r3 y los fixes llevan `"pending-r3-verification"`.
- Los caveats del plan van literales, con una excepción que explico abajo. También agregué el caveat "MLE not efficient across the whole TCP at N=100" (claim `mle_offcenter_sbr10_tail`) y el de SBR=10 pendiente de la autora.
- **Conflicto con la indicación del orquestador.** Pidió que el caveat pase a "≈0.8 δ (slope 0.79±0.02)". Con los 400 patrones que él mismo fijó, mi ajuste da 0.755±0.009 en el centro y 0.862±0.012 en L/4 (δ ∈ {2,5,10}). El 0.79 viene de la fig 8 con 50 patrones. No escribí un número que los datos contradicen: el caveat dice "≈0.75–0.86 δ según la posición", con los valores, y deja el 0.79 para re-evaluar cuando la fig 8 se corra con 400 patrones. **El orquestador tiene que decidir la redacción final.**

## (f) `check_provenance.py`
- Aplana `main.tex` (`\input`/`\include` recursivos, sin comentarios), corre `agent-team/bin/check_provenance.py` en un temporal con cwd=ROOT y además verifica:
  - `\pnum` y `\pnumse` (este último exige `se`);
  - `number_keys` de `provenance.json`;
  - `numbers` y `script` de `claims.json`;
  - la sincronía de `numbers.tex`.
- Imprime `all checks pass` y sale con 0 solo si todo pasa.
- Estado actual: `claims: 37, numbers: 204, provenance entries: 0 — all checks pass`.
- `tests/test_check_provenance.py` tiene 11 tests:
  - pasa el caso OK;
  - fallan un `\src` sin resolver, un `\pnum` desconocido, `\pnumse` sin se, `number_keys` faltante, un número o script de claims faltante, `numbers.tex` desincronizado y un `\input` inexistente;
  - salida de la CLI 0 en el proyecto real y 1 con un `\src` roto;
  - formato de `numbers.tex`.

## (g) Cierre: salidas reales
- `python -m unittest discover -s tests` → `Ran 177 tests in 31.254s — OK`.
- `python -m unittest tests.test_acceptance -v` → los **8 tests OK**, incluido el de figuras (existen las 8 PDFs de W2/W3 y `structure/figures.json`) y el de procedencia.

## Sin resolver / avisos
- **SE por encima de 2 % en algunos números de desalineación del ingenuo:**

  | número | SE relativo |
  |---|---|
  | sesgo del ingenuo, centro, δ=10 | 2.05 % |
  | sesgo del ingenuo, L/4, δ=10 | 2.8 % |
  | σ del ingenuo, L/4, δ=10 | 4.3 % |
  | sesgo del ingenuo, centro, δ=2 | 2.2 % |

  La dispersión la domina la variación entre patrones. Para bajar a < 2 % harían falta ~500–800 patrones, o citar esos números con su SE.
- La pendiente del iterativo no es −0.5 (ver (c)).
- El caveat del sesgo del ingenuo espera la decisión del orquestador (ver (e)).
- `numbers.tex` no se compiló porque no hay pdflatex local. Las macros usan `\ifcsname` (e-TeX).
