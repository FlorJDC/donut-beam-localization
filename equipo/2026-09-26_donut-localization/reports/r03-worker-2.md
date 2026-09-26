# Ronda 3: Worker 2 (estilo, figuras 1 a 4, make_all_figures, figures.json)

## Qué hice
- **`scripts/_paperstyle.py`**: sigue la API del plan al pie de la letra.
  - `COLORS` (Okabe-Ito, claves exactas), `CYCLE`, `CMAP_SEQ="viridis"`, `CMAP_DIV="RdBu_r"`, `COL1=3.4`, `COL2=7.0`.
  - `apply()` corre al importar: sans 8 pt, mathtext dejavusans, ejes de 0.6, ticks hacia adentro, pdf.fonttype 42, leyenda sin marco, bbox tight, 300 dpi.
  - Funciones: `figure()` (constrained_layout), `panel_label()`, `savefig()`, `report()`, `write_summary()`, `parse_args()` (`--quick`, `--no-cache`) y `cached()` (npz en `data/mc/`, `*_quick.npz` aparte).
  - Importa `_paperconfig`, así que `src/` queda en `sys.path`. Worker 3 ya lo usa.
- **Figuras** `scripts/fig_{1_schematic,2_vectorial,3_crb_maps,4_scaling}.py` → `paper/figures/fig{1..4}_*.pdf`, más `data/fig{1..4}_summary.json` y los caches `data/mc/fig{1..4}_*.npz`.
  - Usan solo la API pública de donutloc y `_paperconfig`.
  - Todas son de ancho doble (7.0 in), con letras de panel.
  - Revisé cada una renderizada a PNG en `work/figcheck/` y corregí superposiciones de leyendas y textos, etiquetas cortadas y escalas saturadas.
- **`scripts/make_all_figures.py`**:
  - Descubre `scripts/fig_*.py` y los ordena por número.
  - Corre cada uno en un subproceso (`--quick`, `--no-cache`, `--only`).
  - Exige exit 0 y un PDF nuevo en la ruta derivada del nombre, e imprime la tabla script → pdf → tiempo → estado.
  - Sale con 1 si falla alguno.
- **`structure/figures.json`**: las 8 entradas `{id, pdf, script, question, caption}` con los nombres exactos del plan.
  - Las captions 5 a 8 salen de la descripción del plan. Conviene que W3 o R4 las ajusten al contenido final.
  - Las captions respetan los caveats de texto del state:
    - el límite r→0 frente a S27;
    - el CRB monótono solo para L<360 nm;
    - el SBR de la cámara por píxel = 500/81;
    - los 10.3 fotones en el límite;
    - 0.78 solo para eps≲0.01;
    - 4.9 nm es una escala;
    - `recenter=False` es un artefacto;
    - el sesgo ingenuo ≈0.75 δ;
    - el sesgo sin fondo ≈−0.34.
  - No se afirma concordancia con Caprile.

## Contenido de cada figura
- **Fig1.**
  - (a) Perfiles radiales LG 300, LG 327.1 (curvatura igualada) y vectorial correcta.
  - (b) I/ρ² en eje log de ρ, que muestra la curvatura del cero: LG 300 = 8.37e-5 frente a la vectorial 7.042e-5 nm⁻².
  - (c) Geometría del TCP (L=100 nm para que se vea) sobre el mapa de la exposición 0.
- **Fig2.**
  - (a–c) Mapas 121² en ±500 nm: circular correcta, opuesta y lineal (x), todos con `make_vectorial_beam(mode="exact")`, sin grilla cartesiana.
  - (d) Perfiles con las profundidades del cero.
  - (e) CRB centro (límite) frente a L para los 5 haces.
  - (f) Diferencia en % respecto de LG 300.
- **Fig3.**
  - (a–d) `crb_map(zero_policy="limit")` en ±L, grilla de 161², para L=50 y 100, sin fondo y con SBR=10, en escala log con contornos.
  - (e) Cortes en x e y: límite 1.605 frente al punto S27 1.802 en r=0; SBR=10 continuo en 1.960.
  - (f) `limit_to_point_ratio` frente a L (fwhm 300 y 360), con puntos numéricos y la asíntota 2/√5.
- **Fig4.**
  - (a) CRB frente a L: límite numérico, forma cerrada y S27; ajuste en `L_FIT` sombreado; divergencia en 360.3 nm y la rama L>360 punteada.
  - (b) CRB frente a N.
  - (c) CRB frente a SBR (S31 más numérico en r=0 y r→0; tiende a S27, no al límite sin fondo).
  - (d) Cámara ideal, pixelada 9×9 sin fondo y con SBR_c=500 por píxel, frente a MINFLUX; fotones para 5 nm.

## Tiempos (Windows, Python 3.8)
| script | frío (sin cache) | `--quick` | con cache |
|---|---|---|---|
| fig_1 | 9.4 s | 6.9 s | 3–6 s |
| fig_2 | 56.6 s (50 s de cálculo vectorial exacto) | 25.7 s | 4–7 s |
| fig_3 | 6.1 s | 5.5 s | 4–5 s |
| fig_4 | 8.9 s | 6.8 s | 5–7 s |

`make_all_figures.py` completo, con las 8 figuras y cache: 32 s, 8/8 OK, exit 0.

## Números impresos (NUMBER), y su cruce con `data/paper_numbers.json`
- **Cruce.** Todas las claves que comparto con `paper_numbers.json` coinciden a ≤1e-12 relativo:
  - fig1: 6 claves;
  - fig2: 25;
  - fig3: 9;
  - fig4: 14.
- **Claves propias.** Las que no están en `paper_numbers` llevan el prefijo `figN_`.
- **Fig1.**
  - fwhm_nm = 300
  - lg_ring_diameter_nm = 360.337
  - vectorial_peak_to_peak_diameter_nm = 384.666
  - vectorial_zero_curvature_nm2 = 7.04222e-05
  - vectorial_lg_fwhm_curvature_nm = 327.141
  - vectorial_lg_fwhm_diameter_nm = 320.256
  - fig1_vec_vs_lgcurv_max_rel_dev_rho_le_100nm = 0.00811
  - fig1_tcp_L_nm = 100
- **Fig2.**
  - vectorial_zero_depth_correct = 0
  - vectorial_zero_depth_wrong_handedness = 0.845241
  - vectorial_zero_depth_linear = 0.371648
  - CRB centro (N=100, límite) en L = 50 / 100 / 150:

    | haz | L=50 | L=100 | L=150 |
    |---|---|---|---|
    | vec correcta | 1.60041 | 3.32324 | 5.33104 |
    | opuesta | 147.886 | 89.5884 | 84.1967 |
    | lineal | 60.7742 | 38.0613 | 35.723 |
    | LG 300 | 1.60510 | 3.36341 | 5.48647 |
    | LG 327 | 1.60124 | 3.32964 | 5.35144 |

  - Cocientes en L = 50 / 100 / 150:
    - vec/LG327 = 0.999487 / 0.998076 / 0.996188;
    - vec/LG300 = 0.997084 / 0.988056 / 0.97167.
- **Fig3.**
  - crb_center_lg_L50_N100_nm = 1.60510
  - crb_center_point_S27_L50_N100_nm = 1.80247
  - crb_center_sbr10_L50_N100_nm = 1.96006
  - crb_center_sbr10_limit_L50_N100_nm = 1.96006
  - Valores en r=0 numéricos contra S27 y S31: desviación relativa 2.1e-9.
  - Mapas (centro / mínimo / máximo dentro del círculo del TCP):

    | L | fondo | centro | mínimo | máximo dentro del TCP |
    |---|---|---|---|---|
    | 50 | sin fondo | 1.6051 | 1.5919 | 5.2497 |
    | 50 | SBR10 | 1.9601 | 1.9403 | 6.6863 |
    | 100 | sin fondo | 3.3634 | 3.3175 | 9.6642 |
    | 100 | SBR10 | 4.1654 | 4.0664 | 13.683 |

  - limit_to_point_ratio en L = 5 / 50 / 100 / 150: 0.894388 / 0.890497 / 0.878049 / 0.855267.
  - limit_to_point_ratio_quadratic = 0.894427.
  - Numérico frente a forma cerrada, desviación máxima: 1.0e-8.
- **Fig4.**
  - crb_exponent_L = 1.00435
  - crb_exponent_N = −0.5000
  - crb_center_lg / S27 / sbr10 = 1.60510 / 1.80247 / 1.96006
  - Límite en L = 355 / 365 / 700 nm: 301.9 / 350.3 / 6.94 nm.
  - camera_sigma_nm = 3.16228
  - camera_pixelated_9x9_N400_nm = 5.20455
  - camera_perpixel_sbr500_N600_nm = 4.96418
  - camera_sbr_perpixel_to_total = 6.17284
  - Fotones para 5 nm:
    - cámara ideal 400;
    - pixelada 433.40;
    - pixelada con SBR_c=500 por píxel: 591.43;
    - MINFLUX L=50, S27: 12.996;
    - MINFLUX, límite: 10.305;
    - MINFLUX, SBR=10: 15.367.

## Estado de la aceptación (al cierre)
- `test_every_declared_figure_exists_with_its_script`: OK.
- `python -m unittest tests.test_acceptance`: 7/8 OK. Falla solo `test_manuscript_provenance_passes`, porque todavía no existe `structure/claims.json`, que es de W1.
- sha256 de `tests/test_acceptance.py` = 7d198853…2303, intacto.

## Límites / pendientes
- **Tests.** No agregué tests unitarios del tooling de figuras (`_paperstyle.cached`, `make_all_figures.pdf_for`), porque `tests/test_*.py` es de W1 según el plan. Queda sugerido para R4.
- **Figuras de W3.** `make_all_figures` las corrió desde su cache en el estado de ese momento. Si W3 regenera fig8 con `--no-cache` después del fix de W1, basta volver a correr.
- **`--quick`.** Sobrescribe los mismos PDF, porque solo el cache es aparte. La versión final tiene que salir de una corrida sin `--quick`, y los PDF actuales salieron así.
- **Fig1 (c).** Usa L=100 nm por visibilidad, no L_REF=50. Está declarado en el panel y en la caption.
