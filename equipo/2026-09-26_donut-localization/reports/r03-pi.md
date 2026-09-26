# Ronda 3 — PI: plan

## Estado leído
- **Ledger.** 48 verificadas. Hay 5 afirmaciones vivas:
  - 4 refutadas de r2: SE de la desalineación, `crb(limit)` con N en array, docstring de `l_schedule` y beam vectorial lineal interpolado.
  - 1 unclear: W11, la pendiente del iterativo frente a N.
- **Otras unclear, todas no bloqueantes:**
  - SE del iterativo: se resuelve con bootstrap en esta ronda.
  - `recenter=False`: es un caveat de texto.
  - SBR=10 para `mle_efficiency_center`: decisión del orquestador (inbox r1), pendiente de confirmación de la autora.
  - `rmse`, Masullo Eq. 4.2: queda documentado que no se cita.
- **Checks:** 142 tests OK.
- **Aceptación:** SKIP, porque `data/` y `structure/` están vacíos. `paper/figures` y `paper/sections` también están vacíos.
- **Inbox.** Las 4 directivas siguen vigentes:
  - límite r→0 frente a S27, reportar los dos;
  - priorizar eps y desalineación en R5;
  - fwhm = 300 como Balzarotti S17, y la tabla de Masullo se reproduce con ~360;
  - SBR=10 para la eficiencia, y la superficiencia sin fondo es un RESULTADO.
- **Caveats de texto que el paper debe respetar** (van a `structure/claims.json` como `caveat`):
  - Caprile F=1: 428 frente a 401 nm, NO concuerda.
  - SBR de la cámara por píxel: 500/81.
  - W9 en el límite r→0 da 10.3 fotones.
  - Con eps=0, el CRB es monótono solo para L<360 nm.
  - L_opt ≈ 0.78 fwhm √eps vale solo para eps ≲ 0.01.
  - Sesgo ingenuo ≈ 0.75 δ.
  - Los 8 nm de `recenter=False` son un artefacto del disco de búsqueda.
  - 4.9 nm (V10) es una escala, no el argmin.
  - El sesgo del MLE sin fondo es ≈ −0.34 nm.

## Objetivo de R3
Dejar todo listo para que R4 solo escriba el manuscrito:
- fixes y W11 resueltos;
- pipeline de números, registro de claims y chequeo de procedencia funcionando con un `paper/` stub;
- 8 figuras PDF;
- la aceptación en verde. Debería poder pasar ya en R3 con el stub.

## Decisiones (defaults declarados)
- **Parámetros canónicos.** Un único módulo, `scripts/_paperconfig.py` (dueño: W1). Su contenido va abajo y es literal, para que W2 y W3 lo usen desde el primer minuto.
- **Estilo.** `scripts/_paperstyle.py` (dueño: W2). La API está fijada abajo. W3 la importa. Si todavía no existe cuando la necesita, arranca por los cálculos y deja el ploteo para después.
- **Nombres de archivo.** Fijos: `scripts/fig_<n>_<name>.py` → `paper/figures/fig<n>_<name>.pdf`.
- **Números en la salida de las figuras.** Cada figura imprime `NUMBER <clave> = <valor> <unidad>`. Si el número es una clave de `data/paper_numbers.json`, usa esa misma clave. Además escribe `data/fig<n>_summary.json`, que es chico y se commitea.
- **MC pesado.** Se cachea en `data/mc/*.npz` (gitignored, ya en `.gitignore`). Cada script acepta `--quick` (cache aparte `*_quick.npz`) y `--no-cache`.
- **Consistencia.** Los números MC que cita el paper salen de `compute_paper_numbers.py` o del sweep de W1, con semilla 42 y los parámetros de `_paperconfig`. Una figura que muestra el mismo punto usa los mismos parámetros y semilla, así que da el mismo valor. El verificador lo cruza.
- **Macros de números para R4.** `compute_paper_numbers.py` también escribe `paper/generated/numbers.tex` con `\pnum{clave}` (valor formateado). Es el único escritor de ambos archivos. R4 escribe `\pnum{k}\src{k}`.

### `scripts/_paperconfig.py` (literal, W1 lo crea primero)
```python
# -*- coding: utf-8 -*-
"""Canonical parameters shared by compute_paper_numbers.py and every fig_*.py (nm units)."""
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
SEED = 42
FWHM = 300.0                     # LG donut size parameter (Balzarotti Eq. S17)
FWHM_MASULLO = 360.0             # reproduces Masullo table
WAVELENGTH, NA, N_MEDIUM, FILLING = 640.0, 1.4, 1.518, 5.0 / 3.0
L_REF, N_REF = 50.0, 100
SBR_MLE = 10.0                   # mle_efficiency_center convention
L_FIT = (5.0, 10.0, 20.0, 40.0)  # crb_exponent_L: L << fwhm
N_FIT = (100, 300, 1000, 3000, 10000)
L_LIST = (50.0, 100.0, 150.0)
N_REP_MLE = 10000
SIGMA_PSF = 100.0                # camera: ideal sigma_PSF/sqrt(N)
CAM_PIXEL, CAM_NPIX = 100.0, 9
ITER = dict(n_iter=4, L0=150.0, L_min=25.0, photon_split="equal", recenter=True, rule="fixed")
ITER_N_TOTAL = 1000              # iterative_sigma_nm / camera_sigma_nm at the same photons
ITER_N_LIST = (250, 500, 1000, 2000, 4000, 8000)
ITER_N_REP = 10000
N_BOOT = 2000
EPS_LIST = (0.002, 0.01, 0.05, 0.15)
MIS = dict(L=100.0, N=500, sbr=10, n_patterns=50, n_rep=200)
MIS_DELTAS = (0.0, 2.0, 5.0, 10.0)
DATA = os.path.join(ROOT, "data"); MC = os.path.join(DATA, "mc"); FIGDIR = os.path.join(ROOT, "paper", "figures")
```

### API de `scripts/_paperstyle.py` (W2 la implementa tal cual; W3 la usa)
- `import _paperconfig` al comienzo, así `src/` queda en `sys.path`.
- **Colores.** `COLORS`: dict Okabe-Ito, apto para daltónicos. Claves y valores:

  | clave | color |
  |---|---|
  | `lg` | `#0072B2` |
  | `vec` | `#D55E00` |
  | `opp` | `#009E73` |
  | `lin` | `#CC79A7` |
  | `cam` | `#000000` |
  | `crb` | `#000000` |
  | `mle` | `#0072B2` |
  | `lms` | `#E69F00` |
  | `mlms` | `#009E73` |
  | `honest` | `#0072B2` |
  | `naive` | `#D55E00` |
  | `sbr` | `#56B4E9` |
  | `extra` | `#F0E442` |

  - `CYCLE`: lista Okabe-Ito sin el amarillo al comienzo.
  - `CMAP_SEQ = "viridis"`, `CMAP_DIV = "RdBu_r"`.
- **Tamaños.** `COL1 = 3.4`, `COL2 = 7.0` (pulgadas).
- **`apply()`.** Fija los rcParams:
  - fuente sans 8 pt, mathtext `"dejavusans"`;
  - `axes.linewidth` 0.6, ticks hacia adentro;
  - `pdf.fonttype` 42;
  - leyenda sin marco;
  - `savefig.bbox="tight"`, dpi 300.

  Se llama al importar.
- **`figure(width="single"|"double", aspect=0.75, nrows=1, ncols=1, **subplots_kw) -> (fig, axes)`.**
- **`panel_label(ax, letter)`.** Pone "(a)" en negrita, arriba a la izquierda, en coordenadas de ejes.
- **`savefig(fig, name) -> path`.** Escribe `paper/figures/<name>.pdf` (crea el directorio), cierra la figura e imprime `wrote <path>`.
- **`report(key, value, unit="", se=None)`.** Imprime `NUMBER key = value [± se] unit` y acumula en un dict interno.
- **`write_summary(fig_id)`.** Vuelca lo acumulado a `data/<fig_id>_summary.json`.
- **`parse_args(description) -> args`.** Usa argparse, con `--quick` y `--no-cache`.
- **`cached(name, fn, quick=False, force=False) -> dict`.** Carga `data/mc/<name>[_quick].npz` si existe y no es `force`. Si no, llama a `fn()`, que devuelve un dict de arrays o escalares, guarda con `np.savez` e informa qué hizo.

## Tareas (propiedad de archivos disjunta)

1. **Worker 1: fixes, W11, pipeline de números, claims y procedencia.**
   - **Es dueño de:**
     - `src/donutloc/*`;
     - `tests/test_*.py` salvo `test_acceptance.py`, que NO se toca (hash 7d1988…2303);
     - `scripts/_paperconfig.py`: créalo PRIMERO, literal, desde este plan;
     - `scripts/compute_paper_numbers.py`, `scripts/run_iterative_sweep.py`, `scripts/check_provenance.py`;
     - `data/paper_numbers.json`, `data/iterative_sweep.json`;
     - `structure/claims.json`;
     - `paper/main.tex` (stub), `paper/provenance.json` (stub), `paper/generated/numbers.tex`, `paper/sections/.gitkeep`.

   **(a) Las 4 refutadas de r2.** Cada una lleva test que reproduce el escenario original:
   - **`misalignment_study`.**
     - `sigma_se` pasa a ser la std entre patrones de las σ por patrón dividida por √P.
     - La fórmula vieja queda como `sigma_se_within`.
     - El escenario de prueba es δ=10 en el centro con 20×200.
   - **`fisher.crb(..., zero_policy="limit")`.**
     - Debe transmitir r y N juntos (`np.broadcast_shapes`).
     - `crb(p, [0,0], np.array([100,400]), zero_policy="limit")` tiene que coincidir con `crb_limit` en cada N, y la forma resultante debe ser (2,).
   - **Docstring de `l_schedule` adaptivo.**
     - Hay que decir la verdad: σ_k es el CRB central y subestima el error real de la iteración (emisores hasta L0/4 fuera del centro).
     - Hay que dar la cobertura medida (~97 % con los defaults), no 99 %.
     - Un test chico mide la fracción fuera de L_1/2.
   - **Beam vectorial `polarization="linear"`.**
     - **Preferido:** tabular en armónicos angulares, con tablas radiales en ρ² de los coeficientes cos/sin(mφ) de I, que tiene un número finito de armónicos. Eso da un error ≲1e-4 como en el caso circular. Test: CRB en (7,3) y `crb_limit` a <0.5 % de `mode="exact"` (62.36 y 60.77 nm con L=50, N=100).
     - **Fallback mínimo:** `warnings.warn` en interp lineal cartesiano más un test que documenta el sesgo.

   **(b) Bootstrap del SE.**
   - Añadí `montecarlo.bootstrap_sigma_se(errors, n_boot=1000, seed=42) -> float`, con `errors` de forma (R,2) y σ = √((var_x+var_y)/2) con ddof=1.
   - `experiments.error_stats` devuelve `sigma_se` por bootstrap. La clave no cambia y el viejo queda como `sigma_se_gauss`.
   - `run_mc` agrega `sigma_se_boot`.
   - Test: con datos gaussianos coincide con σ/(2√R) a ±10 %.

   **(c) W11.** `scripts/run_iterative_sweep.py` escribe `data/iterative_sweep.json`.
   - **Protocolo:** el mismo que W10/W11 (`experiments.iterative_minflux`, `_paperconfig.ITER`, sin fondo, semilla 42), con `ITER_N_LIST` e `ITER_N_REP=10000`.
   - **Salidas por N:**
     - σ, SE por bootstrap (`N_BOOT`), rmse;
     - σ por iteración, L_k, cámara ideal, cociente.
   - **Salidas globales:**
     - pendiente log-log en todo el rango y para N≥500, con IC 95 % por bootstrap (remuestreo de las repeticiones por N);
     - rango del cociente con la cámara;
     - también en N=1000: SBR=10, `rule="adaptive"` y `recenter=False`, este último marcado como artefacto.
   - `--quick` usa 1000 repeticiones.
   - En el reporte, di honestamente si la pendiente es compatible con −0.5 y cuál es el rango. Si no converge, el claim queda como "−0.52…−0.55 según el rango".

   **(d) `scripts/compute_paper_numbers.py`.**
   - Es el único escritor de `data/paper_numbers.json`. Formato: `{key: {"value", "unit", "script", "description"}}`, más `"se"` en los números MC.
   - Solo API pública más `_paperconfig`.
   - `--quick` con opción de MC reducido; la corrida final es sin `--quick`.
   - Claves mínimas de la aceptación:
     - `fwhm_nm` = 300;
     - `crb_center_lg_L50_N100_nm`: límite r→0, ≈1.6051;
     - `crb_exponent_L`: ajuste log-log con `L_FIT`, N=100, límite;
     - `crb_exponent_N`: con `N_FIT`, L=50;
     - `mle_efficiency_center`: σ_MLE/CRB en el centro, SBR=10, L=50, N=100, `N_REP_MLE`, SE < 2 %;
     - `vectorial_zero_depth_correct`, `vectorial_zero_depth_wrong_handedness`;
     - `iterative_sigma_nm`: leído de `data/iterative_sweep.json` en N=1000; si falta, se llama al sweep;
     - `camera_sigma_nm` = SIGMA_PSF/√1000 = 3.162, ideal y declarado.
   - Además, cada número que el paper va a citar:
     - S27 puntual 1.8025 y el cociente límite/puntual en L=5/50/100/150;
     - 2/√5, con el caso cuadrático;
     - ejes de la elipse límite (V4);
     - S31 con SBR=10 = 1.960;
     - fotones para 5 nm (W9, incluido el límite de 10.3);
     - V12 con fwhm=360;
     - cámara pixelada 9×9 (5.2045) y 591 fotones con la convención por píxel;
     - vectorial: lineal 0.3716, D_pp 384.67, curvatura, fwhm equivalente 327.14/320.26 y CRB vectorial/LG en `L_LIST` (W6);
     - superficiencia del MLE sin fondo (≈0.84) y sesgo en r=(2,0) (≈−0.34, con repeticiones suficientes para SE ≤ 0.01);
     - LMS centro = S27;
     - V9 (tabla de degradación por eps) y W12 (L_opt y CRB por eps, cociente L_opt/(fwhm√eps));
     - desalineación (`MIS`, `MIS_DELTAS`): sesgo del ingenuo y del honesto por δ, pendiente sesgo/δ, σ, piso MC, con el SE ya corregido;
     - iterativo: SE, pendientes con IC, cociente, SBR=10 y adaptivo.
   - Después de escribir, también genera `paper/generated/numbers.tex`, con `\newcommand` por clave vía `\csname` y `\pnum{clave}`.

   **(e) `structure/claims.json`.**
   - Es una lista de `{key, statement, numbers:[claves], script, figure (si hay), question: "R1".."R5", verified_round, status, caveat}`.
   - Hay que cubrir todas las verificadas relevantes de r1 y r2 más las nuevas de r3. Estas últimas llevan `status` `"pending-r3-verification"`.
   - Los caveats de texto de este plan van literales en `caveat`.

   **(f) `scripts/check_provenance.py`.**
   - Aplana `paper/main.tex`, resolviendo recursivamente `\input` y `\include` de `sections/` y `generated/`, a un archivo temporal.
   - Corre `agent-team/bin/check_provenance.py <flat> paper/provenance.json` como subproceso.
   - Además verifica:
     - cada `\pnum{k}` del tex existe en `paper_numbers.json`;
     - las claves `number_keys` de cada entrada de `provenance.json` existen;
     - cada `numbers` de `claims.json` existe, y cada `script` existe;
     - `numbers.tex` está sincronizado con `paper_numbers.json`.
   - Imprime `all checks pass` y sale con 0 solo si todo pasa.
   - Crea stubs mínimos válidos:
     - `main.tex` con `\input{generated/numbers}` y un `\input{sections/...}` vacío o comentado;
     - `provenance.json` como `{}`.
   - Así funciona ya y R4 solo llena `paper/`.
   - Agrega un test unitario del wrapper con un tex temporal que tenga un `\src` sin resolver, que debe fallar, y otro OK.

   **(g) Cierre.**
   - `python -m unittest discover -s tests` en verde.
   - `python -m unittest tests.test_acceptance -v`: todo pasa salvo, si hace falta, el test de figuras mientras W2 y W3 terminan. Reporta el estado real.
   - En el reporte, una tabla clave → valor → SE.

2. **Worker 2: estilo, figuras 1 a 4, `make_all_figures` y `figures.json`.**
   - **Es dueño de:**
     - `scripts/_paperstyle.py`: primero, con la API exacta del plan;
     - `scripts/fig_1_schematic.py`, `fig_2_vectorial.py`, `fig_3_crb_maps.py`, `fig_4_scaling.py`;
     - `paper/figures/fig1_schematic.pdf` … `fig4_scaling.pdf`;
     - `data/fig1..4_summary.json`, `data/mc/fig1..4*`;
     - `scripts/make_all_figures.py`, `structure/figures.json`.
   - **Reglas:**
     - Solo la API pública de `donutloc` más `_paperconfig`.
     - Cada script imprime sus números con `report()`.
     - Estilo publicación: colores de `COLORS`, ejes con unidades en nm, panel labels (a)(b)…
   - **Fig1: esquema.**
     - (a) Perfil radial normalizado LG (fwhm 300) frente a vectorial circular correcto (`vectorial.radial_profile`), con zoom logarítmico del cero y curvatura igualada (LG 327.14).
     - (b) Geometría del TCP: 3 donas en un círculo de diámetro L más el centro, con mapa de la suma o de una dona y los ceros marcados.
   - **Fig2: dona vectorial.**
     - (a–c) Mapas de intensidad para circular correcta, opuesta y lineal (grilla modesta ~121², cacheada).
     - (d) Perfiles y profundidades del cero: 0 (reportar `vectorial_zero_depth_correct`), 0.845 y 0.372.
     - (e) CRB centro frente a L: vectorial correcta, LG 300 y LG 327.
     - (f) Opcional: CRB lineal y opuesta.
     - Para la lineal usa `mode="exact"` o el beam arreglado por W1. NO uses el interp cartesiano de grilla 5 nm.
   - **Fig3: mapas de CRB y discontinuidad.**
     - (a,b) `crb_map` con `zero_policy="limit"` en el FOV ±L, para L=50 y 100, N=100, sin fondo y SBR=10.
     - (c) Corte radial que muestra el valor puntual S27 en r=0 frente al límite.
     - (d) Cociente límite/puntual frente a L (`closed_forms.limit_to_point_ratio`), con la asíntota 2/√5.
   - **Fig4: escalamientos y cámara.**
     - CRB frente a L (limit y S27, con el rango del ajuste `L_FIT` marcado y la divergencia en 360 nm con eps=0).
     - CRB frente a N (pendiente −1/2).
     - CRB frente a SBR.
     - Comparación con la cámara ideal σ_PSF/√N y pixelada 9×9 con la convención por píxel (500/81) frente a N, y fotones para 5 nm.
   - **`make_all_figures.py`.**
     - Descubre `scripts/fig_*.py` por glob.
     - Los corre en subprocesos, pasando `--quick` si se da.
     - Falla si alguno falla e imprime una tabla script → pdf → OK.
   - **`structure/figures.json`.** Lista de las 8 figuras `{id, pdf, script, question, caption}` con estos nombres exactos:
     - fig1_schematic, fig2_vectorial, fig3_crb_maps, fig4_scaling;
     - fig5_estimators, fig6_iterative, fig7_zero_depth, fig8_misalignment.
     - Los pdf van en `paper/figures/<id>.pdf` y los scripts en `scripts/fig_<n>_<name>.py`.
     - Las captions de 5 a 8 se escriben desde la descripción de este plan.
   - **Cierre.** Corre `make_all_figures.py` completo al final, cuando existan los scripts de W3. Si no existen, reporta cuáles faltan. Reporta los NUMBER impresos.

3. **Worker 3: figuras 5 a 8 (MC).**
   - **Es dueño de:**
     - `scripts/fig_5_estimators.py`, `fig_6_iterative.py`, `fig_7_zero_depth.py`, `fig_8_misalignment.py`;
     - `paper/figures/fig5_estimators.pdf` … `fig8_misalignment.pdf`;
     - `data/fig5..8_summary.json`, `data/mc/fig5..8*`.
   - **Reglas:**
     - Usa `_paperconfig` y `_paperstyle` (API en este plan; si aún no existe, arranca por los cálculos cacheados).
     - Solo API pública.
     - `--quick`, cache en npz, `report()` de cada número.
     - SE por bootstrap: `montecarlo.bootstrap_sigma_se`, que W1 agrega. Hasta que exista, deja el hook y vuelve a correr al final.
   - **Fig5: estimadores.** L=50, N=100, semilla 42.
     - (a) Sesgo y (b) σ/CRB frente a la posición sobre x, de 0 a L (dentro y fuera del TCP), para MLE, LMS (`lms_tcp`) y mLMS (`mlms_tcp`) con SBR=10.
     - (c) MLE sin fondo en el centro: σ/CRB_límite frente a N (superficiencia ≈0.84, frente a S27 ≈0.75), y sesgo hacia el centro en r=(2,0) (≈−0.34 nm).
     - El punto centro con SBR=10 usa `N_REP_MLE` y debe coincidir con `mle_efficiency_center`.
   - **Fig6: iterativo frente a cámara.**
     - Lee `data/iterative_sweep.json` de W1. Si no existe, calcula una versión quick con el mismo protocolo y avisa.
     - σ final ± SE frente a N_total, con la cámara ideal σ_PSF/√N y la pendiente ajustada.
     - Inset o panel con σ por iteración frente a L_k.
     - Puntos de SBR=10 y adaptivo.
     - `recenter=False` no se grafica como física: a lo sumo, en gris y rotulado como artefacto del disco de búsqueda.
   - **Fig7: profundidad finita del cero.**
     - CRB centro frente a L para `EPS_LIST` más eps=0 (`experiments.eps_L_sweep`, gaussiano, N=100), con L_opt marcado (`optimal_L`).
     - L_opt/(fwhm√eps) frente a eps, mostrando que 0.78 vale solo para eps ≲ 0.01.
     - Opcional: el mínimo fuera del centro con eps=0.002 (V10, 4.9 nm como escala).
   - **Fig8: desalineación.**
     - `misalignment_study` con `MIS` y `MIS_DELTAS`: |sesgo| del ingenuo y del honesto frente a δ, con el piso de ruido MC y la recta ≈0.75 δ; σ y rmse frente a δ; CRB honesto.
     - Las barras de error usan `bias_abs_se` y `sigma_se`, que ahora es entre patrones (fix de W1).
     - Corre la versión final con `--no-cache` DESPUÉS de que el fix de W1 esté en `src`. Revisá el docstring de `sigma_se`.
   - **Cierre.** Reporta los tiempos de cómputo, los NUMBER y los SE de cada número citado; el objetivo es SE < 2 %.

## Para el verificador de R3
- Reproducir independientemente:
  - los 4 fixes con sus escenarios;
  - el bootstrap;
  - W11, con otra semilla;
  - las claves de aceptación, con la aceptación corriendo en verde;
  - un muestreo de los NUMBER de las figuras contra `paper_numbers.json`, en los puntos compartidos.
- Code-reviewer: leer el diff; ver que `check_provenance` falle cuando corresponde; ver que las figuras usen solo la API pública.

## LO QUE MÁS SE PODRÍA HACER
1. **Confirmación explícita de la autora sobre SBR=10 para `mle_efficiency_center`.** Cierra la unclear de r1. Costo ~0.
2. **Mapas de MC del MLE en todo el FOV** (eficiencia por celda con SE < 2 %, acelerado con numba). Completa R3 más allá del corte en x. Costo: ½ ronda.
3. **Dona 3D top-hat y TCP 3D** (nota B §3). Extiende R1 a z. Costo: 1 ronda.
4. **Aberraciones de Zernike → profundidad del cero** (nota B §4.2). Da una eps física en lugar de un pedestal. Costo: 1 ronda.
5. **Cámara realista** (lectura y exceso EM ×2). Hace la comparación con la cámara más justa. Costo: ½ ronda.
6. **Regla óptima de reducción de L derivada analíticamente.** Es un aporte original sobre el iterativo. Costo: 1 ronda con verificación.
7. **MLE con prior o fusión de iteraciones en el iterativo** (Balzarotti p. 6). Costo: ½ ronda.
8. **Contraste residual de polarización → eps con la dona vectorial elíptica** (López). Costo: ½ ronda.
