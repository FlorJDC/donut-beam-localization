# r01 — worker-A: auditoría del repositorio (para el PDF de síntesis ≤ 5 pp.)

Solo lectura; no se editó nada fuera de este informe. Base: commit `f411211`, árbol limpio salvo
`equipo/2026-09-28_sintesis-pdf/` (sin trackear).

## 1. Suite de tests (resultado real)

Entorno de esta sesión: **Linux, Python 3.11.15, numpy 2.4.6, scipy 1.17.1, matplotlib 3.11.2**
(no es el entorno declarado en CLAUDE.md/README: Windows, Py 3.8, numpy 1.24).

`python3 -m unittest discover -s tests` → **Ran 224 tests, FAILED (errors=1, skipped=1)**.
- ERROR: `test_vectorial.TestEnergy.test_flux_same_for_both_hands` (`tests/test_vectorial.py:126`):
  `np.trapz` no existe en numpy ≥ 2.0 (`AttributeError`). Es un defecto de compatibilidad del
  **test**, no del paquete: `grep trapz src scripts` no encuentra nada más. Arreglo trivial:
  `getattr(np, "trapezoid", None) or np.trapz`. (No lo edité: tarea de solo lectura.)
- SKIP: `test_paper_tooling...test_unknown_key_compiles_with_tectonic` (no hay `tectonic`).
- `python3 -m unittest tests.test_acceptance` → **8 tests OK**.
- El cierre r05 del trabajo anterior declaraba 224 tests pasando (en Py 3.8/numpy 1.24); el
  recuento coincide (224), la diferencia es solo el error de numpy 2.

## 2. Resultados verificados principales (todos `status: verified` en `structure/claims.json`)

Parámetros por defecto: LG fwhm = 300 nm (`fwhm_nm` = 300.0), TCP = 3 donas en círculo de
diámetro L + dona central, fotones multinomiales con N fijo, 4 exposiciones de igual peso.

R-1. **CRB en el centro: límite r→0 < valor puntual de Balzarotti Eq. S27 (sin fondo).**
  L=50, N=100: `crb_center_lg_L50_N100_nm` = 1.6050955712820285 nm (forma cerrada
  `crb_center_lg_L50_N100_closed_nm` = 1.6050955676488428) vs `crb_center_point_S27_L50_N100_nm`
  = 1.8024719062965908 nm. Razón `limit_to_point_ratio_L50` = 0.8904968571447623; cero cuadrático
  `limit_to_point_ratio_quadratic` = 0.8944271909999159 = 2/√5 (numérico fwhm=1e9:
  `limit_to_point_ratio_quadratic_numeric` = 0.8944271930298965). Con SBR=10 la discontinuidad
  desaparece: `crb_center_sbr10_L50_N100_nm` = 1.9600587058960406 vs límite
  `crb_center_sbr10_limit_L50_N100_nm` = 1.9600587033709. Fig. 3 (`scripts/fig_3_crb_maps.py`).
  (Chequeo propio: 1.60509557/1.80247191 = 0.89050, coincide con `limit_to_point_ratio_L50`.)

R-2. **Escalamiento con L, N, SBR y fotones para 5 nm.** Pendientes log-log: `crb_exponent_L` =
  1.0043477933205587 (solo L ≪ fwhm, ajuste L = 5–40 nm; caveat: 1.07 en 25–150 nm),
  `crb_exponent_N` = -0.49999999999999983. Fotones para 5 nm, L=50: `photons_for_5nm_L50_sbrinf`
  = 12.995619891953863, `photons_for_5nm_L50_sbr10` = 15.367320522235445,
  `photons_for_5nm_L50_sbr5` = 17.933955450896335. Con ε=0 el CRB diverge en
  `eps0_crb_divergence_L_nm` = 360.33672263593496 nm (diámetro del anillo). Fig. 4
  (`scripts/fig_4_scaling.py`).

R-3. **Cámara vs MINFLUX (L fijo).** Cámara ideal σ_PSF/√N: `camera_ideal_N400_nm` = 5.0 nm
  (400 fotones para 5 nm); cámara pixelada 9×9: `camera_pixelated_9x9_N400_nm` =
  5.204549139040088 nm; con SBR_c=500 por píxel, `camera_photons_for_5nm_perpixel_sbr500` =
  591.4337185203674 fotones. MINFLUX L=50 SBR=10 necesita 15.37 (`photons_for_5nm_L50_sbr10`):
  ≈ 26× menos fotones que la cámara ideal (cociente 400/15.37 = 26.0, **derivado**, no es clave).
  Fig. 4(d) (`scripts/fig_4_scaling.py`). La comparación con 1000 fotones del paper
  (`iterative_ratio_to_camera_N1000` = 0.160561253867469) es del protocolo **iterativo** (ver R-8).

R-4. **Dona vectorial vs LG.** Con la handedness correcta el cero es exacto
  (`vectorial_zero_depth_correct` = 0.0) y la LG igualada en curvatura (fwhm
  `vectorial_lg_fwhm_curvature_nm` = 327.1411191076221) reproduce el CRB: `crb_vec_over_lgcurv_L50_N100`
  = 0.9994871462333158, `_L100_` = 0.9980764218786861, `_L150_` = 0.9961877338914108. Respecto de
  la LG fwhm=300: `crb_vec_over_lg300_L150_N100` = 0.9716701164811385. Handedness opuesta:
  `vectorial_zero_depth_wrong_handedness` = 0.8452409316543332, CRB L=50
  `crb_center_vec_opposite_L50_N100_nm` = 147.88600722506212 nm; lineal:
  `vectorial_zero_depth_linear` = 0.37164824493965715, `crb_center_vec_linear_L50_N100_nm` =
  60.77420112870313 nm (vs `crb_center_vec_correct_L50_N100_nm` = 1.600414754550864). Fig. 2
  (`scripts/fig_2_vectorial.py`).

R-5. **Estimadores: MLE vs LMS/mLMS.** Con SBR=10 el MLE es eficiente en el centro:
  `mle_efficiency_center` = 0.9912292806217669 (σ = `mle_sigma_center_sbr10_nm` =
  1.9428675785187703 nm; L=50, N=100). Sin fondo es superficiente y sesgado:
  `mle_nobg_sigma_over_crb_limit_center` = 0.8397049526504217, `mle_nobg_bias_x_r2_nm` =
  -0.34445140228271487 nm. LMS sin fondo en el centro = valor S27: `lms_nobg_sigma_over_s27_center`
  = 0.9999785960805649 (1.123× el límite, `lms_nobg_sigma_over_crb_limit_center` =
  1.1229445514533927). Cola de outliers del MLE fuera del centro con N bajo:
  `mle_sbr10_sigma_over_crb_x15_N100` = 1.5116262832636163, desaparece a N=1000
  (`..._N1000` = 1.0042031326322736). LMS/mLMS con fondo, fuera del centro, fuertemente sesgados
  hacia el centro — números SOLO en `data/fig5_summary.json` (no en `paper_numbers.json`):
  sesgo x a x0=25 nm: LMS `fig5_lms_sbr10_bias_x_x25_nm` = -14.35 nm, mLMS
  `fig5_mlms_sbr10_bias_x_x25_nm` = -5.19 nm, MLE `fig5_mle_sbr10_bias_x_x25_nm` = 0.68 nm; mLMS en
  el centro `fig5_mlms_sbr10_sigma_over_crb_x0` = 1.357. Fig. 5 (`scripts/fig_5_estimators.py`).
  Caveat abierto: SBR=10 como condición de `mle_efficiency_center` pendiente de confirmar por la autora.

R-6. **Cero imperfecto ε ⇒ L óptimo.** Pedestal gaussiano, N=100, fwhm=300: `L_opt_eps0p002_nm` =
  10.484412357992904, `L_opt_eps0p01_nm` = 23.280315521660803, `L_opt_eps0p05_nm` =
  50.336650731671014, `L_opt_eps0p15_nm` = 80.80970073749592 nm; L_opt/(fwhm√ε) =
  0.7814619578870331 (ε=0.002), 0.7760105173886934 (ε=0.01) → "0.78·fwhm·√ε" vale solo para
  ε ≲ 0.01 (0.695 a ε=0.15). CRB en L_opt: `crb_opt_eps0p01_nm` = 1.6786740447081847 nm.
  Degradación a L=50, ε=0.01: `crb_eps_degradation_L50_eps0p01_gaussian` = 1.2277033606831182.
  Un pedestal constante equivale exactamente a fondo (`eps_constant_pedestal_over_s31_L100_eps0p05`
  = 1.0000000000000004). Fig. 7 (`scripts/fig_7_zero_depth.py`).

R-7. **Desalineación del TCP (δ por cero, dirección aleatoria; L=100, N=500, SBR=10).** El MLE
  ingenuo (asume TCP ideal) tiene sesgo ∝ δ: población sin ruido en el centro
  `misalignment_naive_bias_over_delta_pop_center_d2` = 0.7471281876667419,
  `..._d10` = 0.7834885724752224 (pendiente `misalignment_naive_bias_over_delta_pop_center` =
  0.7768263559835796); en (L/4,0) `misalignment_naive_bias_over_delta_pop_Lq` = 0.8506956555160232.
  MC: `misalignment_naive_bias_abs_center_d10_nm` = 7.708458376248982 nm vs honesto (patrón
  verdadero) `misalignment_honest_bias_abs_center_d10_nm` = 0.17970451698485712 nm (piso MC
  `misalignment_mc_floor_center_d10_nm` = 0.16433435969870003). rmse centro δ=10: ingenuo
  `misalignment_naive_rmse_center_d10_nm` = 6.228370069068393 vs honesto 1.8646607028352271 nm.
  Fig. 8 (`scripts/fig_8_misalignment.py`). Nota: la clave `misalignment_naive_bias_over_delta_Lq`
  está `superseded` en claims.json (no usarla).

R-8. **MINFLUX iterativo (R4).** 4 iteraciones L=150→25 nm, N=1000, sin fondo:
  `iterative_sigma_nm` = 0.507739266193721 nm vs cámara ideal `camera_sigma_nm` =
  3.1622776601683795 nm (razón 0.160561253867469); pendiente N≥500 `iterative_slope_N_ge_500` =
  -0.5148302737955811 (IC95 [-0.5196, -0.5096]). Fig. 6 (`scripts/fig_6_iterative.py`).

## 3. Relevancia para el caso objetivo: p-MINFLUX NO iterativo (4 exposiciones de igual duración, 20 MHz, L fijo, sin control del pulso central)

Nota de modelo: el repo ya modela exactamente "4 exposiciones de igual peso + N fijo multinomial"
(`src/donutloc/photons.py`, Eq. S4/S30), sin temporización de pulsos, vida media, tiempo muerto ni
fotones por ventana. Eso es lo que cubren worker-B/C.

| Resultado | Relevancia | Por qué |
|---|---|---|
| R-1 CRB centro (límite vs S27, 2/√5) | **Alta** | Es el CRB del TCP a L fijo; con fondo real (SBR finito) vale el valor Eq. S31, continuo. |
| R-2 Escalamiento L/N/SBR | **Alta** | Dimensiona el experimento a L fijo (σ ∝ L/√N; costo del fondo). |
| R-6 L óptimo con ε | **Alta** | Con L fijo, elegir L es LA decisión de diseño; depende de ε medido. |
| R-7 Desalineación | **Alta** | En p-MINFLUX los 4 haces vienen de caminos/retardos distintos; un error de posición de ceros produce sesgo ~0.75–0.85·δ si no se calibra el patrón (el MLE "honesto" lo elimina). |
| R-5 Estimadores | **Alta/media** | Con L fijo el emisor puede estar lejos del centro: LMS/mLMS muy sesgados ahí; MLE con modelo correcto es lo recomendable; la cola de outliers a N=100 importa en p-MINFLUX con pocos fotones. |
| R-3 Cámara vs MINFLUX | Media | Contexto (ganancia en fotones a L fijo, ~26× para 5 nm con SBR=10). |
| R-4 Vectorial vs LG | Media | Justifica usar LG igualada por curvatura; la handedness/polarización es chequeo experimental. |
| R-8 Iterativo | **Baja** | Protocolo con L decreciente y re-centrado: fuera del caso objetivo; solo sirve como referencia de techo. |

Casos límite: el repo supone SBR fijo en toda posición (open points del paper); en un TCP fijo el
fondo por exposición constante hace el SBR dependiente de la posición (Eq. S28/S29, soportado por
`bg_per_exposure` en `photons.probabilities`, pero no explotado en las figuras) → hueco para
worker-C.

## 4. Figuras a reusar en el PDF breve (máx. 3)

1. `paper/figures/fig3_crb_maps.pdf` — mapas de CRB a L fijo (50/100 nm, sin fondo y SBR=10) +
   discontinuidad y razón 2/√5: el resultado central y es exactamente la geometría p-MINFLUX a L fijo.
2. `paper/figures/fig7_zero_depth.pdf` — L óptimo vs ε: la decisión de diseño con L fijo.
3. `paper/figures/fig8_misalignment.pdf` — sesgo por desalineación: efecto sistemático realista
   clave para 4 haces entrelazados.
Alternativa si hace falta un panel de escalamiento/cámara: `fig4_scaling.pdf` (panel d). No
reusar `fig6_iterative.pdf` (baja relevancia). `fig5` es densa (4 paneles) para un PDF de 5 pp.

## 5. Huecos y defectos concretos

1. **Test roto con numpy ≥ 2**: `tests/test_vectorial.py:126` usa `np.trapz`;
   `requirements.txt`/`pyproject.toml` piden `numpy>=1.20` sin cota superior, así que una
   instalación nueva rompe la suite (224 tests, 1 error).
2. **README desactualizado respecto del entorno**: dice "Tested on Windows with Python 3.8, numpy
   1.24…"; no menciona numpy 2 ni Python 3.11; no menciona `informe/` (el entregable nuevo;
   `informe/figures/` existe y está vacío).
3. No hay `tectonic`/pdflatex en este entorno: `paper/main.pdf` no se puede recompilar aquí (1 skip).
4. PDFs de builds intermedios versionados: `equipo/2026-09-26_donut-localization/work/review/r03/tex/main.pdf`
   y `revtex_main.pdf` (propios, no papers de terceros; ruido en el repo, no violan privacidad).
5. `docs/private/` no existe en este checkout (normal: gitignored).
6. Números de LMS/mLMS citados en el texto de Fig. 5 viven solo en `data/fig5_summary.json`, no en
   `data/paper_numbers.json` (si el informe los cita, necesita esa procedencia o nuevas claves vía
   `compute_paper_numbers.py`).
7. Backlog del trabajo anterior (`equipo/2026-09-26_donut-localization/state.json`):
   (a) confirmar SBR=10 para `mle_efficiency_center`; (b) open points: normalización rmse vs
   Masullo Eq. 4.2 (texto principal no disponible), diámetro del anillo vs Caprile a F=1 (428 vs
   401 nm, sin acuerdo), cola pesada del iterativo a N=250; (c) extensiones: dona 3D/top-hat,
   aberraciones de Zernike como origen de ε, dipolo de absorción, cámara realista, regla óptima de
   reducción de L, MLE con prior en el iterativo.
8. Open points del paper (`paper/sections/open_points.tex`): todo 2D, N fijo, fondo igual en todas
   las exposiciones con SBR fijo en toda posición, sin orientación del dipolo, sin aberraciones más
   allá de ε, **sin deriva**, cámara ideal, sin 3D ni desalineación axial.
9. Ausentes para p-MINFLUX (no están en el repo): temporización de pulsos a 20 MHz y ventanas
   (TCSPC/gating), vida media y crosstalk entre ventanas, tiempo muerto del detector, blinking /
   fotoblanqueo durante las 4 exposiciones, N Poisson (no fijo), fondo por exposición no uniforme.
   Ninguna comparación con SimuFLUX (no hay nota de literatura sobre Marin & Ries 2026 en
   `docs/literature/`).

## Qué no hice / límites
- No regeneré figuras ni `paper_numbers.json` (tiempo); los números se copiaron del JSON versionado.
- No verifiqué los resultados del paper; solo su estado `verified` en `structure/claims.json` y un
  cociente de control (R-1).
