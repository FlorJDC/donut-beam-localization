# C. SimuFLUX: simulación realista de MINFLUX (Marin & Ries): nota para simulación

**Cita.** Z. Marin y J. Ries, "Evaluating MINFLUX experimental performance in silico",
*Nature Communications* **17**, 246 (2026). DOI: [10.1038/s41467-025-66952-w](https://doi.org/10.1038/s41467-025-66952-w).
Recibido el 22/07/2025, aceptado el 19/11/2025, © 2025, acceso abierto CC BY 4.0. El código
(MATLAB y Python, GPLv3) está en `github.com/ries-lab/SimuFLUX`; la versión de registro, en
Zenodo 10.5281/zenodo.17523756.

**Fuentes leídas** (carpeta `Papers\Evaluating MINFLUX experimental performance in silico`,
extraídas con PyMuPDF; las figuras se miraron renderizadas):

| Clave | Archivo | Páginas |
|---|---|---|
| [MR-main] | `s41467-025-66952-w.pdf` | 8 (artículo con Métodos) |
| [MR-SI] | `41467_2025_66952_MOESM1_ESM.pdf` | 10 (Supplementary Information: SI Figs. 1–10 y Tablas 1–2) |
| [MR-PR] | `41467_2025_66952_MOESM2_ESM.pdf` | 5 (Peer Review File) |
| [MR-nb] | `41467_2025_66952_MOESM3_ESM.zip` → `SimuFLUX_estimators.nb` / `.pdf` | 11 (notebook de Mathematica con la derivación de los estimadores) |
| [MR-code] | repositorio público `ries-lab/SimuFLUX` (rama principal, consultado el 28/09/2026) | solo se leyó; no se copió nada |

**Convención.** "p. N" es la página del PDF (1-based). Las ecuaciones del artículo se citan como
"Eq. n". [MR-code] no forma parte del artículo revisado: lo que sale de ahí se marca
**(código)** y podría cambiar entre versiones. El SI no trae texto de métodos aparte de las
leyendas; casi toda la información del modelo está en los Métodos del artículo (pp. 4–7).

---

## 1. Qué es y qué modela (resumen)

SimuFLUX simula el microscopio completo: barre una PSF en un patrón o recorre un conjunto de PSFs
en posición fija, cuenta fotones en cada exposición, estima la posición y re-centra el patrón
con un escáner rápido (EOD) y uno lento (galvo o piezo) [MR-main, p. 1]. La Tabla SI 1
[MR-SI, p. 10] lo compara con otros simuladores publicados (Masullo 2021 y 2022, He 2022,
Slenders 2023, Srambickal 2024, Rosati 2024, Liu 2025). Es el único de la lista que incluye a la
vez PSF vectorial, 3D, pinhole, búsqueda ("scouting"), MINFLUX iterativo con re-centrado,
tiempos reales del hardware, parpadeo, blanqueo, movimiento, varios fluoróforos, CRB,
RMSE/sesgo y archivos de secuencia de Abberior.

Mensaje central: el CRB **no** capta la degradación que producen las imperfecciones reales
(desalineación de la máscara de fase, fondo, sesgo del estimador, parpadeo, tiempos muertos).
Hace falta simular el experimento completo y comparar RMSE, STD y sesgo con el CRB
[MR-main, p. 2].

## 2. Modelo óptico (PSF)

- **PSFs analíticas** `PsfGauss2D` y `PsfDonut2D`, y una **vectorial** `PsfVectorial` basada en
  Leutenegger et al. 2006 (integral de difracción de Debye vectorial calculada con FFT)
  [MR-main, p. 5]. La máscara de fase del plano focal trasero (vórtice, top-hat para la dona 3D,
  half-moon/PhaseFLUX) se pasa como `phasepattern`. Admite aberraciones (Zernike) y PSFs
  experimentales calibradas [MR-main, p. 5].
- **Dona analítica que se usa para derivar los estimadores** [MR-main, p. 6, Eq. 7], con la
  integral normalizada a A:
  $f_d(\mathbf x)=A\,\frac{x^2}{4\pi\sigma^4}e^{-x^2/2\sigma^2}$.
  Su aproximación cuadrática es $f_q=A\,x^2/(4\pi\sigma^4)$ (Eq. 9), válida para $|x|<L/2$
  (SI Fig. 6g). Equivalencia con nuestra convención: $\sigma=\mathrm{fwhm}/\sqrt{8\ln2}$, y
  entonces $L_\sigma^2\equiv L^2/8\sigma^2=L^2\ln2/\mathrm{fwhm}^2$.
- **(código)** `PsfDonut2D`: $I=4\ln2\,(r^2/\mathrm{fwhm}^2)\,e^{-4\ln2\,r^2/\mathrm{fwhm}^2}+\texttt{zerooffset}$,
  con fwhm = 310 nm por defecto. **Sin el factor e**: el pico del anillo vale $1/e\approx0.368$,
  es decir, la dona está normalizada a la **misma potencia** que una gaussiana de pico 1 (las dos
  integrales valen $\pi\,\mathrm{fwhm}^2/4\ln2$). La PSF vectorial se normaliza igual: se divide
  por el máximo del foco gaussiano (sin máscara) de la misma potencia. En consecuencia,
  `zerooffset` se mide **relativo al pico gaussiano**, no al pico del anillo.
- **Cero imperfecto**: un pedestal constante sumado a la PSF de excitación (`psf.zerooffset`). En
  SI Fig. 5 se usa 0.01 [MR-SI, p. 4]; según la normalización anterior equivale a ≈2.7 % del pico
  del anillo (inferencia nuestra a partir del código). **Tamaño de la esfera**: la PSF se convoluciona
  con una esfera (`beadsize`), y una esfera grande "llena" el cero, lo que actúa como un fondo
  aparente (SI Fig. 7c) [MR-main, p. 5].
- **Detección**: la señal colectada es la excitación × la eficiencia de detección a través de un
  pinhole **centrado en el centro del patrón**. El EOD no está des-escaneado en detección, así que
  el pinhole no sigue a cada exposición [MR-main, p. 4]. El pinhole puede desalinearse (SI Fig. 2).
- **(código)** Parámetros vectoriales por defecto (`default_microscope.yaml`): λ = 635 nm,
  NA = 1.35, inmersión y muestra con n = 1.406, lente con n = 1.518, polarización circular,
  cintura del haz de 7 mm sobre una apertura de 6.5 mm, píxel de 10 nm, emisión a 700 nm. La
  desalineación de la máscara se da con `maskshift` (relativa a la apertura, 0–1).

## 3. Geometría del patrón y esquema iterativo

- Patrones en órbita: K_o = 3, 4 o 6 puntos sobre un círculo de **diámetro L**, con o sin
  exposición central [MR-main, p. 6; MR-nb, p. 5]. Las secuencias de Abberior usan hexágono
  (6 + centro opcional), cuadrado o triángulo [MR-SI, p. 9, Fig. 10i].
- **Secuencia** = lista de iteraciones. Cada una tiene patrón, L, tiempo de permanencia
  (`patDwellTime`), repeticiones (`patRepeat`), potencia (`pwrFactor`), límite de fotones
  (`phtLimit`), `ccrLimit` (CFR) y umbral de fondo `bgcThreshold`. Los parámetros globales son
  `ctrDwellFactor`, `damping`, `headstart`, `stickiness`, `loclimit` y `fieldGeoFactor`
  [MR-main, p. 5].
- **Parámetros de Abberior citados en el artículo**:
  - Pista experimental (Atto647N, secuencia de tracking 2D), última iteración: L = 75 nm,
    patRepeat = 3, patDwellTime = 0.3 ms, phtLimit = 500, stickiness = 4, ccrLimit = 0.9,
    potencia láser 280 µW [MR-main, p. 7].
  - Tiempos muertos medidos en un Abberior MINFLUX (iMSPECTOR 16.3.21317): movimiento del EOD
    0.011 ms, cálculo del estimador 0.015 ms, movimiento del galvo 0.04 ms. Con firmware anterior
    (Vogler et al. 2025): EOD 0.005 ms y galvo 0.02 ms [MR-main, pp. 4–5].
  - Última iteración de tracking optimizada (SI Tabla 2) [MR-SI, p. 10]: patGeoFactor 0.28 → 0.4,
    phtLimit 20 → 10, patDwellTime 100 → 50 µs, patrón hexágono → cuadrado. D_max pasa de 2.5 a
    4.2 µm²/s (+70 %) [MR-main, p. 4; Fig. 2l].
- **(código)** `L = patGeoFactor × 360 nm` (`sim_sequencefile.py`). Esto sale del código y no del
  artículo; el YAML anotado dice `patGeoFactor × 360 × (λ/642)`. Las secuencias de ejemplo del
  repositorio dan:
  - `Tracking_2D`: L ≈ 284, 302, 151, 101 nm; phtLimit 40, 20, 20, 20; dwell 100 µs;
    pwrFactor 1, 1, 2, 3; ctrDwellFactor 0.16; CFR solo en la iteración 3 (ccrLimit 0.8).
  - `Imaging_2D`: L ≈ 288, 288, 151, 76, 40 nm; phtLimit 160, 150, 100, 100, 150; patRepeat
    1, 5, 5, 5, 5; pwrFactor 1, 1, 2, 4, 6; dwell 1 ms.

  **La potencia sube a medida que L baja**, para compensar la menor señal cerca del cero.
- **Tiempo de permanencia del centro (código)**: la exposición central dura
  `t_punto × K_o × ctrDwellFactor`, de modo que las exposiciones **no tienen el mismo peso**. El
  CFR se calcula como $\mathrm{CFR}=(n_c/\texttt{ctrDwellFactor})/\sum_{\rm órbita}n_i$.
- **Lazo** [MR-main, p. 4]: búsqueda en una grilla hexagonal con una PSF **gaussiana** (aproxima
  la "órbita del pinhole" de Abberior; subestima el blanqueo durante la búsqueda). Cada candidato
  se mide hasta alcanzar `phtLimit`, hasta blanquearse o perderse (fotones < `bgcThreshold`) o
  hasta que falle el chequeo de CFR (CFR > `ccrLimit`). Cada chequeo puede fallar `stickiness`
  veces. Un filtro adicional típico en los análisis es cfr < 0.5. Los reposicionamientos llevan
  amortiguamiento (código: factor $2^{-\texttt{damping}}$).

## 4. Estadística de fotones y fondo

- Fotones esperados por exposición [MR-main, p. 5, Eq. 2]: $I_i=I_0\,(f(\mathbf x-\mathbf b_i)+bg)$.
  Probabilidades: $p_i=I_i/\sum_jI_j$ (Eq. 3). La normalización de la PSF y el factor $I_0$ se
  cancelan.
- **(código)** Los conteos por exposición son **Poisson** con media = brillo (kHz) × PSF ×
  potencia × dwell, más el fondo de autofluorescencia `background` (kHz) × dwell × potencia. **El
  fondo escala con la potencia láser y con el dwell, pero no con la posición del emisor ni con el
  patrón**: es un fondo fijo por exposición (tipo Balzarotti Eq. S28). El número total de fotones
  **no** está fijo: se repiten barridos hasta llegar a `phtLimit`.
- **Fuentes de fondo estudiadas** (SI Fig. 5) [MR-SI, p. 4]: cero imperfecto (zerooffset 0.01),
  autofluorescencia (30 kHz), un segundo fluoróforo de igual brillo en [50, 50, 400] nm o barrido
  lateral y axialmente; L = 75 nm y 100 fotones/localización. En DNA-PAINT, las hebras difusas
  (~1/µm³) generan un fondo que varía en la escala espacio-temporal del barrido (Fig. 2f,g).
- **Blanqueo** [MR-main, p. 5]: presupuesto total de fotones con distribución exponencial, del que
  se restan los fotones emitidos (antes del pinhole). **Parpadeo**: tiempos on/off exponenciales
  (t_on, t_off), con la fracción de tiempo "on" en cada ventana. Reactivaciones: PAINT ∞,
  (d)STORM ≈ 2 en promedio, PALM 0. **Movimiento**: difusión con saltos
  $W_k=\sqrt{2D\,dt}\,\mathrm{randn}$ (Eq. 1) o pasos con tiempos exponenciales. Las vibraciones
  se modelan como fluoróforos que oscilan. Cada medición puntual se hace con el fluoróforo en
  posición fija, así que dt debe ser menor que el tiempo de medición puntual [MR-main, p. 5].

## 5. Estimadores [MR-main, pp. 5–6; MR-nb]

Todos se pueden implementar en una FPGA. El fondo se resta de los conteos (Eq. 5), así que las
fórmulas se derivan sin fondo.

- **LSQ linealizado** en x = 0 (Balzarotti Eqs. S43–S48): $\hat{\mathbf x}=(J^TJ)^{-1}J^T(\hat{\mathbf p}-\mathbf p(0))$
  (Eq. 6).
  - Dona: $\hat{\mathbf x}=\frac{1}{L_\sigma^2-1}\sum_i\hat p_i\mathbf b_i$ (Eq. 8). Es idéntico a
    Balzarotti S50 y a nuestro `lms_tcp`; se verificó numéricamente (§ comparación).
  - Aproximación cuadrática: $\hat{\mathbf x}=-\sum_i\hat p_i\mathbf b_i$ (Eq. 10).
  - Gaussiana: $\frac{1}{L_\sigma^2}\sum\hat p_i\mathbf b_i$ sin centro (Eq. 12) y
    $(K_o+e^{L_\sigma^2})/(K_oL_\sigma^2)\sum\hat p_i\mathbf b_i$ con centro (Eq. 13).
- **LSQ iterativo**: se re-linealiza en $\mathbf x_{0,i}=\hat{\mathbf x}_{i-1}$. Converge en pocas
  iteraciones y extiende el rango sin sesgo. Los jacobianos analíticos para K_o = 3, 4, 6, con y
  sin centro, están en [MR-nb, pp. 6–8]. (Código: hasta 15 iteraciones, tolerancia 0.1 nm, NaN si
  |x| > 1.2 L.) Con aproximación cuadrática queda sesgo más allá de ~30 nm (SI Fig. 6e).
- **1D**:
  - MLE cuadrático con dos puntos en ±L/2 (Balzarotti):
    $\hat x=\frac L2\frac{\sqrt{n_1}\mp\sqrt{n_2}}{\sqrt{n_1}\pm\sqrt{n_2}}$ (Eq. 14). La medición
    central sirve para elegir la raíz.
  - La solución directa con K = 2 coincide con el MLE (Eq. 15) [MR-nb, p. 4].
  - Con K = 3 (0, ±L/2) y **el fondo como parámetro libre**:
    $\hat x=\frac L4\frac{\hat n_2-\hat n_1}{2\hat n_0-\hat n_1-\hat n_2}$ (Eq. 16, ref. Rosati 2024).
- Resultados (SI Fig. 6, 500 fotones; 6f con 15 fotones de señal):
  - El estimador simple tiene un sesgo fuerte fuera del centro, y su STD es "baja pero sin
    sentido".
  - El iterativo y el LSQ gaussiano son los mejores con poca luz.
  - **El MLE tiene el sesgo más fuerte con pocos fotones** [MR-SI, p. 5].
- Fondo (Fig. 2d, 200 fotones de señal y 200 de fondo en x = 0): si el fondo no se considera, el
  sesgo es fuerte. Desaparece si el fondo se estima bien y se resta, o si se ajusta como
  parámetro libre [MR-main, p. 2; SI Fig. 5d–f]. La STD sigue siendo mayor.
- "Creeping" (SI Fig. 7): durante las primeras iteraciones la posición estimada converge
  lentamente, con un error de decenas de nm. Aparece cuando la búsqueda es imprecisa, cuando se
  usan esferas grandes (100 nm) o con fondo alto (30 kHz) mal estimado; se observó
  experimentalmente en 50 moléculas.

## 6. Figuras de mérito [MR-main, p. 7]

- **CRB** numérico según Masullo et al. 2022: la información de Fisher se arma con las derivadas
  por diferencias finitas de las probabilidades normalizadas en $\mathbf x\pm\varepsilon/2$, e
  incluye el fondo, otros fluoróforos y las aberraciones.
  - Se reporta **por eje** ("x: σ_CRB").
  - **(código)** ε = 1 nm; el denominador lleva un regularizador, $p_i+10^{-4}$;
    $\sigma=\sqrt{\mathrm{diag}(F^{-1})}$, y después se divide por $\sqrt{\bar N}$ con $\bar N$ el
    número medio de fotones (fondo incluido) de las localizaciones.
  - Referencia rápida (código): $L/\sqrt{8N}$.
- **RMSE** (Eq. 17), **STD** (Eq. 18) y **sesgo** (Eq. 19), por coordenada, sobre Q
  fluoróforos. **Errata**: la Eq. 17 del PDF no tiene raíz cuadrada (así escrita es el MSE),
  mientras que la Eq. 18 sí la tiene. Hay que leerla como RMSE = $\sqrt{\cdot}$.
- Localizaciones **válidas/inválidas** (pasan los chequeos de fotones y CFR), filtro cfr < 0.5,
  fracción de pistas exitosas (> 90 localizaciones consecutivas con error < 100 nm, SI Fig. 10) y
  D_max (el D con el que la mitad de las pistas son exitosas).

## 7. Resultados cuantitativos útiles

- **Desalineación de la máscara de vórtice** (SI Fig. 1) [MR-SI, p. 1]: L = 75 nm, 6 puntos +
  centro, 0.1 ms por punto, 45 fotones (caso alineado), estimador iterativo. Un corrimiento de
  0.1–0.5 mm produce sesgo dependiente de z (hasta ~−50 nm fuera de foco) y RMSE de hasta
  ~20× σ_CRB alineado; la STD crece hasta ~4–5×. **El CRB no refleja esa degradación**
  [MR-main, p. 2].
- **Pinhole desalineado** hasta 400 nm: efecto despreciable salvo la pérdida de eficiencia de
  detección (STD/σ_CRB entre 0.9 y 1.1, sesgo nulo) [MR-SI, p. 2].
- **Aberraciones** de 0.15 rad (esférica, coma, astigmatismo): afectan mucho más a la dona 3D
  (top-hat) que a la 2D o a PhaseFLUX; con fondo pasa lo mismo [MR-main, p. 2; SI Fig. 3;
  Fig. 2b].
- **Vibraciones**: error extra del orden de la amplitud y casi independiente de la frecuencia
  (0.1–10 kHz, amplitud de 5 nm → STD/σ_CRB ≈ 2.5) [MR-SI, p. 3].
- **Parpadeo** (Fig. 2e, SI Fig. 8; t_on = t_off = 100 µs; dwell del patrón 400 µs; 100 fotones):
  $\mathrm{STD}(r,N)^2=\sigma_{fl}(r)^2+\sigma_{\rm CRB}(N)^2$, donde $\sigma_{fl}$ no depende de N
  sino de la cinética y del número r de repeticiones del patrón. Error extra de hasta ~10 nm; puede
  hacer falta r > 20. En la figura, σ_CRB ≈ 2.8 nm con N = 100 y ≈ 0.9 nm con N = 1000 (valores
  leídos del gráfico) [MR-main, p. 2].
- **Densidad**: fluoróforos activos cercanos "tapan" el cero y producen errores grandes
  (SI Fig. 5g,h: sesgo de hasta ~16 nm con un vecino a ~50 nm). Hay que optimizar t_off
  (SI Fig. 9, dSTORM con 5000 fotones y 2 reactivaciones).
- **Tracking**: los tiempos muertos son el límite principal. La optimización lleva D_max de 2.5 a
  4.2 µm²/s. El D se subestima (el tracking MINFLUX da 1.5 frente a 2.9 µm²/s reales en Fig. 2k).

## 8. Revisión por pares [MR-PR]

- El revisor #3 objeta la falta de validación experimental de la PSF y de la cinética. Los
  autores responden que la PSF sigue a Leutenegger 2006 y que la cinética on/off es práctica
  estándar.
- Un estimador LSQ 1D en FPGA tarda 300 ns, y el iterativo (10 iteraciones) unos 3 µs [MR-PR,
  p. 5].
- La curva roja de Fig. 2d corresponde al estimador iterativo **con fondo no corregido**
  [MR-PR, pp. 3–4].

## 9. Dudas y lo que no se pudo leer

- La Tabla SI 1 y las figuras son imágenes: se leyeron renderizadas y los valores numéricos son
  aproximados (leídos del gráfico).
- La leyenda de SI Fig. 6e remite a "Methods, Eq. 16" para la aproximación cuadrática, pero la
  Eq. 16 es la solución directa con fondo libre. Parece una referencia mal numerada.
- Qué patrón y qué L se usaron en Fig. 2e y SI Fig. 8 no se dice explícitamente (el ejemplo 4 del
  repositorio no se inspeccionó).
- Hay una discrepancia entre el artículo y el código: el artículo menciona L = 75 nm para las
  simulaciones de SI, mientras que `Tracking_2D` da una última L ≈ 101 nm. No se resolvió.
