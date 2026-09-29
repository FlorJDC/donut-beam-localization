# C'. SimuFLUX (Marin & Ries 2026) frente a `donutloc`: comparación, recomendaciones y checklist

Esta nota compara el modelo de SimuFLUX (resumido en `C_insilico_minflux.md`, claves [MR-*])
con nuestra implementación (`src/donutloc`, `scripts/fig_*.py`, `paper/sections/*.tex`) y termina
con un checklist de errores frecuentes en simulaciones de MINFLUX. Todo sale de la literatura
publicada, del código público de SimuFLUX (marcado **(código)**) o de cálculos propios, que
están marcados **[chequeo propio]**. Esos chequeos **no** están verificados por un verificador
independiente: son orientativos hasta que se reproduzcan (ver § 4).

---

## 1. Chequeos numéricos cruzados [chequeo propio]

Scripts en el scratchpad de la sesión (`xcheck.py`, `xcheck2.py`), que solo usan la API pública de
`donutloc`. No se modificó ningún archivo del repositorio, y
`python -m unittest discover -s tests` pasa (OK, 1 skipped).

**1a. Un número suyo reproducido con nuestro paquete.** Fig. 2e y SI Fig. 8 de [MR] dan
σ_CRB ≈ 2.8 nm con N = 100 y ≈ 0.9 nm con N = 1000 (valores leídos del gráfico), con L = 75 nm y
dona de fwhm = 310 nm (default del código). Con `donutloc` (dona LG, fwhm = 310, L = 75, N = 100,
sin fondo, en el centro):

| patrón | valor puntual (S27) | límite r→0 | estilo SimuFLUX* | L/√(8N) |
|---|---|---|---|---|
| TCP (3 + centro) | 2.764 | 2.449 | 2.765 | 2.652 |
| hexágono sin centro | 2.764 | 2.764 | 2.765 | 2.652 |
| hexágono + centro (dwell igual) | 2.764 | 2.568 | 2.765 | 2.652 |
| cuadrado sin centro | 2.764 | 2.764 | 2.765 | 2.652 |

\*Reimplementación de la lógica de `calculateCRBdirect` (código): diferencias finitas con
ε = 1 nm, denominador $p_i+10^{-4}$, σ por eje $\sqrt{(F^{-1})_{xx}}$. Da 2.765 nm para cualquier
patrón, lo que coincide con su ≈ 2.8 nm; con N = 1000 da 0.874 nm, frente a ≈ 0.9.
**Conclusión: su σ_CRB en el centro es el *valor puntual*, no el límite.** Con la exposición
central, el valor puntual queda 13 % por encima de nuestro "centre CRB" (límite) para el TCP
(2.764 frente a 2.449).

**1b. LMS.** Nuestro `estimators.lms_tcp` y su `est_donutLSQ1_2D` / Eq. 8 coinciden a precisión
de máquina (diferencia máxima 1.8e-15 nm).

**1c. CRB por eje fuera del centro** (x = 20 nm, L = 75, N = 100). La anisotropía importa si se
comparan barridos en x:

| patrón | $\sigma_x$ | $\sigma_y$ | nuestro escalar $\sqrt{(\Sigma_{xx}+\Sigma_{yy})/2}$ |
|---|---|---|---|
| TCP | 3.65 | 3.42 | 3.54 |
| hexágono sin centro | 3.87 | 3.19 | 3.55 |

**1d. Peso de la exposición central** (hexágono + centro, `ctrDwellFactor` de Abberior): el CRB
límite vale 2.57 nm con 0.16 (tracking) y 2.23 nm con 1.0 (imaging). Cambia un 15 %.

**1e. Fondo físico frente a SBR fija** (TCP, L = 50, fwhm = 300, N = 100, fondo constante por
exposición elegido para que SBR(0) = 10):

| x (nm) | 0 | 10 | 15 | 25 | 37.5 | 50 |
|---|---|---|---|---|---|---|
| SBR(x) con fondo constante | 10.0 | 12.1 | 14.6 | 22.6 | 37.6 | 57.5 |
| CRB con SBR fija = 10 | 1.96 | 2.27 | 2.96 | 5.99 | 8.71 | 12.20 |
| CRB con fondo constante | 1.96 | 2.24 | 2.85 | 5.77 | 8.34 | 11.39 |

El CRB cambia ≤ 7 % en x ≤ L, pero la SBR real varía en un factor ~6.

**1f. Pedestal no modelado → sesgo** (L = 75, fwhm = 310). El pedestal es constante, ε = 0.0272
del pico del anillo, que equivale a su `zerooffset` = 0.01 porque su pico vale 1/e. El sesgo es
poblacional: se estima sobre los conteos esperados, sin ruido. MLE en un disco de ≈ 0.95 L.

| x (nm) | 0 | 10 | 20 | 30 | 40 |
|---|---|---|---|---|---|
| MLE ingenuo (modelo con ε = 0) | 0 | +4.2 | +1.8 | +3.8 | +9.3 |
| MLE honesto (conoce ε) | 0 | 0 | 0 | 0 | 0 |
| LMS ingenuo (S50 sin 1/s) | 0 | −3.0 | −8.4 | −16.5 | −26.4 |
| CRB con ε (N = 100) | 3.59 | 3.84 | 4.84 | 7.28 | 11.05 |

Es cualitativamente lo mismo que SI Fig. 5c: un cero imperfecto no considerado sesga el
estimador, y el sesgo es del orden del CRB o mayor.

---

## 2. Efecto por efecto: (a) ya lo modelamos · (b) no lo modelamos y convendría · (c) posible conflicto

| Efecto | SimuFLUX [MR] | Nosotros | Clase | Comentario |
|---|---|---|---|---|
| Forma de la dona | LG analítica (fwhm 310, **pico 1/e**, normalización a igual potencia que la gaussiana) o vectorial Debye/FFT (Leutenegger) | LG con pico del anillo = 1 (Balzarotti S17), `beams.py`; vectorial Richards–Wolf, `vectorial.py` (λ 640, NA 1.4, n 1.518) | (a) | Para $p_i$ da lo mismo. Sus defaults vectoriales difieren de los nuestros (λ 635, NA 1.35, n 1.406): no comparar diámetros sin ajustar parámetros. |
| Profundidad finita del cero | `zerooffset` constante, **relativo al pico gaussiano** (0.01 ≈ 2.7 % del anillo) | ε relativo al pico del anillo; modelos "constant" y "gaussian" (Fig. 7 usa el gaussiano) | (a) / (c) | Nuestro modelo constante es el suyo. **Conflicto de convención: ε_nuestro = e·zerooffset.** Declararlo si se citan sus números. |
| Esfera (bead) que llena el cero | convolución con una esfera | — | (b) | Relevante si se calibra con esferas; esfuerzo bajo (convolución 2D de la LG). |
| Polarización / handedness | vectorial circular | vectorial circular/lineal, handedness (Fig. 2) | (a) | Aquí vamos más allá: estudiamos la mano incorrecta y la polarización lineal. |
| Aberraciones (Zernike) | sí; la dona 3D es la más sensible | no | (b) | Fuera de alcance salvo para la discusión; esfuerzo medio (fase en la pupila del módulo vectorial). |
| Desalineación de la **máscara de fase** | corrimiento en el BFP → dona deformada, cero lleno, sesgo dependiente de z; RMSE hasta ~20× σ_CRB | **no**: nuestra Fig. 8 desplaza las *posiciones* de los ceros | (b) / (c) | **Conflicto terminológico**: nuestra "misalignment" es un error de posicionamiento del patrón (calibración EOD/galvo), no desalineación óptica. Conviene renombrarla o aclararlo en `nonidealities.tex`. |
| Pinhole / PSF de detección | excitación × detección con pinhole centrado en el patrón; su desalineación es inocua | no hay detección | (a) por cancelación | Con un pinhole común a todas las exposiciones de un patrón, el factor de detección se cancela en $p_i$. Solo afecta a la SBR y al presupuesto de fotones. Coincide con su SI Fig. 2. Declararlo como hipótesis. |
| Geometría | órbitas K_o = 3, 4, 6 ± centro; L = diámetro; Abberior usa hexágono | TCP (3 + centro), L = diámetro; `polygon_centers` permite M arbitrario | (a) | Misma convención de L. |
| Peso/dwell de la exposición central | `ctrDwellFactor` (0.16 en tracking): centro con otro dwell, y CFR normalizado | todas las exposiciones con igual peso | (b) | Cambia el CRB un 15 % (§ 1d). Esfuerzo bajo: vector de pesos $w_i$ en `photons.make_model`, $p_i\propto w_iI_i$. |
| Estadística | Poisson por exposición (brillo × dwell × potencia); N aleatorio hasta `phtLimit` | multinomial con N fijo (también Poisson en `sample_counts`) | (a) | Equivalentes condicionando en N. La regla de parada por umbral no está modelada. |
| **Fondo** | autofluorescencia constante por exposición (kHz × dwell × potencia), **independiente de la posición del emisor**; más cero imperfecto, vecinos y PAINT | **SBR fija en cada posición** (Eq. S30) en Figs. 5, 6 y 8; `bg_per_exposure` existe pero casi no se usa | (c) | Nuestra elección está declarada (`model.tex`), pero no es física: implica un fondo que crece con la señal. El CRB cambia poco (§ 1e); lo que se afecta son los sesgos y la corrección 1/s del LMS fuera del centro, y la SBR a L chica en el iterativo. |
| Fondo no considerado / mal estimado en el estimador | sesgo fuerte; se corrige restando o ajustando el fondo | nuestros estimadores siempre conocen la SBR exacta | (b) | Falta el caso mal especificado; § 1f muestra el efecto. Esfuerzo bajo. |
| Parpadeo / flickering | on/off exponencial; $\mathrm{STD}^2=\sigma_{fl}^2+\sigma_{\rm CRB}^2$ | no | (b) | Fuera de R1–R5; mencionarlo en la discusión. Esfuerzo medio si se quiere una figura. |
| Blanqueo, reactivaciones, densidad, vecinos | sí | no | (b) | Solo para la discusión. |
| Movimiento / vibración / deriva | difusión, pasos, oscilaciones | no (sin deriva, `open_points.tex`) | (b) | Una vibración de amplitud a suma ~a al error. Discusión. |
| Tiempos muertos EOD/galvo/FPGA | 0.011 / 0.04 / 0.015 ms | no | (b) | Solo importa para tracking y tiempos; fuera de alcance. |
| Estimadores | LSQ simple (= Balzarotti S50), LSQ iterativo, 1D MLE/directo, fondo libre (Eq. 16), MLE | MLE (grilla + refinamiento), LMS (S48, S50), mLMS (S51) | (a) / (b) | LMS idéntico (§ 1b). No tenemos el LSQ iterativo, que es simple y extiende el rango sin sesgo: esfuerzo bajo con el `lms` general re-linealizado. |
| Sesgo del MLE con pocos fotones | "el MLE tiene el mayor sesgo con 15 fotones" | MLE sesgado y superficiente sin fondo (Fig. 5c,d) | (a) | Consistente. |
| CRB | numérico (Masullo 2022), por eje; en el cero exacto da el **valor puntual** (regularizador $10^{-4}$) | $\sqrt{(\Sigma_{xx}+\Sigma_{yy})/2}$; "centre CRB" = **límite** r→0 | (c) de convención | Nuestro "centre CRB" es un 13 % más bajo que el que reportarían SimuFLUX o Abberior para el mismo TCP sin fondo (§ 1a). No es un error (el paper define ambos), pero hay que decirlo al comparar con la literatura. |
| STD / RMSE / sesgo | por eje x (Eqs. 17–19; la Eq. 17 sin raíz es una errata) | promedio por eje de x e y | (a) / (c) | Fuera del centro nuestro σ escalar mezcla ejes anisótropos (§ 1c). Para comparar con sus barridos en x usar `fisher.crb_axes` y la varianza en x. |
| MINFLUX iterativo | secuencias Abberior: L ≈ 288 → 40 nm (imaging) o 284 → 101 nm (tracking), `phtLimit` por iteración, **la potencia sube al bajar L**, CFR, stickiness, damping, búsqueda gaussiana | 4 iteraciones, L = 150 → 25 nm geométrico, reparto igual de N, MLE, re-centrado, SBR fija (Fig. 6) | (a) parcial / (c) | Con SBR fija ignoramos que, a potencia y fondo fijos, la SBR en el centro cae ∝ L² (Eq. S32). Abberior compensa con la potencia. Nuestra L final de 25 nm está por debajo de lo usado en la práctica (40–100 nm). El "creeping" por una primera iteración imprecisa (SI Fig. 7) es análogo a nuestra advertencia sobre la regla adaptativa. |
| Validez (CFR, umbral de fotones) | localizaciones válidas/inválidas, filtro cfr < 0.5 | no | (b) | Filtrar cambia la estadística (selección). Esfuerzo medio. |

---

## 3. Posibles errores o inconsistencias en nuestro trabajo (prioridad para revisar)

1. **Fondo "SBR fija en cada posición"** (`photons.py`, `model.tex`, Figs. 5, 6 y 8). No es un
   bug y está declarado, pero no es el fondo físico que usa SimuFLUX (ni Balzarotti S28/S32). Hay
   tres afirmaciones sensibles:
   - Fig. 5a,b: sesgo y σ/CRB del LMS/mLMS fuera del centro, con la corrección 1/s calculada para
     SBR = 10 constante.
   - Fig. 6c y `iterative_sbr10_sigma_nm`: con fondo fijo, la SBR en L = 25 nm sería mucho menor
     que en L = 150 nm.
   - Fig. 8: sesgo del estimador ingenuo con SBR = 10.

   → Repetir con `bg_per_exposure` (ya implementado) y reportar las dos variantes o justificar la
   elección.
2. **Comparabilidad del "centre CRB"**. Nuestro número por defecto es el límite; SimuFLUX (y en
   la práctica cualquier CRB numérico regularizado) da el valor puntual, un 13 % mayor para el
   TCP sin fondo. → En `crb_center.tex` o `discussion.tex`, una frase que diga qué convención
   usan otros trabajos y cuál es el factor.
3. **"Misalignment" (Fig. 8)**. Lo que simulamos es el error de posición de los ceros; en la
   literatura, "misalignment" es sobre todo desalineación de la máscara de fase, que deforma la
   dona y llena el cero, un efecto cualitativamente distinto y mucho mayor. → Renombrar a
   "pattern-positioning error" o aclararlo en `nonidealities.tex` y en la leyenda de Fig. 8.
4. **Normalización de ε**. Nuestro ε es relativo al pico del anillo; su `zerooffset`, al pico
   gaussiano de igual potencia (factor e ≈ 2.72). Si en la discusión se compara la Tabla de ε con
   valores de SimuFLUX o Abberior, hay que convertir.
5. **Anisotropía fuera del centro**. El σ y el rmse que reportamos promedian x e y. Es correcto
   internamente, pero no es comparable con curvas "x: STD" o "x: RMSE" de otros trabajos (§ 1c).
   `open_points.tex` ya menciona el tema de la normalización del rmse; conviene agregar SimuFLUX
   (por eje x).

No se encontraron errores de fórmula: el LMS, la geometría y la definición de $p_i$ coinciden
exactamente.

---

## 4. Recomendaciones priorizadas

**Estado (28/09/2026):** los ítems 1–5 y 10 están implementados en el manuscrito (secciones `bgphys_*`, `naive_*` y `simuflux_*` de `scripts/compute_paper_numbers.py`). El ítem 6 fue descartado. Los ítems 7–9 quedan para más adelante.

| # | Qué | Archivo o figura afectada | Esfuerzo | Por qué |
|---|---|---|---|---|
| 1 | Variante con **fondo físico** (`bg_per_exposure`) para Fig. 5 (off-centre) y Fig. 6 (iterativo con SBR(L) según S32), o al menos los números clave | `scripts/fig_5_estimators.py`, `fig_6_iterative.py`, `compute_paper_numbers.py`, `model.tex`, `iterative.tex` | bajo (1–2 h de cómputo, API existente) | Elimina la objeción principal de realismo. |
| 2 | Panel o números de **estimador ingenuo frente a honesto con ε no modelado** y con **SBR mal estimada** (p. ej. se asume SBR = ∞ o SBR = 20 cuando la real es 10) | `fig_7_zero_depth.py` (nuevo panel d) o `fig_8_misalignment.py`; `nonidealities.tex` | bajo (reusa `misalignment_population_bias`: MLE sobre conteos esperados) | § 1f: sesgos de +4 a +9 nm (MLE) y de −3 a −26 nm (LMS) con ε ≈ 2.7 %; es el resultado central de [MR] SI Fig. 5. |
| 3 | Frase y número de **conversión de convenciones**: valor puntual frente a límite (13 %), ε frente a zerooffset (×e), CRB por eje frente a escalar | `crb_center.tex`, `discussion.tex`, `open_points.tex` | muy bajo | Evita comparaciones erróneas con la literatura. |
| 4 | Renombrar o aclarar la **desalineación** de Fig. 8 y citar [MR] para la desalineación de la máscara | `nonidealities.tex`, leyenda de Fig. 8, `references.bib` | muy bajo | Precisión terminológica. |
| 5 | Citar SimuFLUX en la introducción y la discusión (Tabla SI 1 como panorama de simuladores) | `introduction.tex`, `discussion.tex`, `references.bib` | muy bajo | Contexto; hoy no se cita. |
| 6 | ~~**Pesos de dwell**~~ **Descartado por la autora (28/09/2026):** en su experimento todas las exposiciones tienen el mismo tiempo; se declara como hipótesis en `open_points.tex`. Texto original: por exposición ($p_i\propto w_iI_i$) y un número con `ctrDwellFactor` | `photons.py` (+ tests), eventualmente `crb_center.tex` | bajo | Acerca el TCP a la práctica comercial (§ 1d: 15 %). |
| 7 | **LSQ iterativo** (re-linealización de `lms`) como cuarto estimador en Fig. 5 | `estimators.py` (+ tests), `fig_5_estimators.py` | bajo–medio | Estimador simple de FPGA con rango sin sesgo; comparación directa con [MR] Fig. 2c. |
| 8 | Esquema iterativo "tipo Abberior": L = 288 → 40 nm, umbral de fotones por iteración, rampa de potencia con fondo fijo | `experiments.py`, `fig_6_iterative.py` | medio | Hace la Fig. 6 comparable con instrumentos reales. |
| 9 | Desalineación de la **máscara de vórtice** en el módulo vectorial (vórtice descentrado en la pupila) → CRB y sesgo del estimador ingenuo | `vectorial.py` (requiere el campo 2D sin simetría), nueva figura | medio–alto | Es el efecto de alineación que más importa según [MR] SI Fig. 1. |
| 10 | Parpadeo ($\sigma_{fl}$), vecinos, vibraciones y esferas: solo mencionarlos en la discusión con cita a [MR] | `discussion.tex`, `open_points.tex` | muy bajo | Delimita el alcance con honestidad. |

Antes de pasar al manuscrito, los números de § 1 deberían reproducirse por un verificador
independiente (metodología del proyecto).

---

## 5. Checklist: errores frecuentes en simulaciones de MINFLUX (que surgen de [MR])

Pensado para auditar simulaciones propias. Cada ítem es una pregunta de sí o no.

**Modelo óptico y de fondo**
1. ¿Está explícita la **normalización de la dona**: pico del anillo = 1 (Balzarotti) o igual
   potencia que la gaussiana (SimuFLUX, pico 1/e)? ¿El "zero offset" o ε se refiere al mismo
   pico? (Diferencia de un factor e.)
2. ¿El **cero imperfecto** está en el modelo *y* en el estimador? Si solo está en los datos, hay
   sesgo (SI Fig. 5c; § 1f).
3. ¿El **fondo** es físico, constante por exposición y proporcional a la potencia y al dwell? ¿O
   es una SBR fija, que implica un fondo que depende de la posición? ¿La SBR baja al achicar L
   (S32)?
4. ¿El estimador **conoce, resta o ajusta el fondo**? Si no, espere un sesgo fuerte hacia el
   centro (LSQ) que crece con |x| (Fig. 2d).
5. ¿El **pinhole** es común a todas las exposiciones (EOD no des-escaneado)? Si se modela la
   detección por exposición, ¿está centrada en el patrón y no en cada haz?
6. ¿Se distingue la desalineación de la **máscara de fase** (deforma la dona) del error de
   **posición del patrón**? ¿Se consideró la dependencia con z?
7. Si se calibra con **esferas**, ¿se tuvo en cuenta que su tamaño llena el cero?

**Estadística y patrón**
8. ¿Los **dwell times y pesos** por exposición son los reales (centro con `ctrDwellFactor`)? ¿El
   CFR está normalizado por el dwell?
9. ¿N fijo (multinomial) o Poisson con umbral de fotones? ¿Se reportan los fotones de señal y de
   fondo por separado?
10. ¿Se aplican los filtros reales (CFR, umbral de fondo, stickiness)? ¿Se reporta la fracción de
    localizaciones válidas y el sesgo de selección que introducen?

**Estimadores y figuras de mérito**
11. ¿El estimador linealizado se usa **fuera de su rango** (|x| ≳ L/4)? Entonces una STD baja
    "no significa nada" si el sesgo es grande (SI Fig. 6).
12. ¿Se reportan **sesgo, STD y RMSE** con el estimador realmente usado, y no solo el CRB? El CRB
    no ve el sesgo por mala especificación.
13. ¿El MLE se evaluó con **pocos fotones**? Puede tener el mayor sesgo (SI Fig. 6f). ¿Se
    reportan el disco de búsqueda y la elección de raíz (1D: dos soluciones)?
14. ¿El **fwhm/σ del estimador** coincide con la PSF real? (Los estimadores por defecto de
    SimuFLUX fijan fwhm = 310 nm.)
15. ¿El RMSE lleva la **raíz**? ¿El error es por eje o promediado sobre ejes? ¿Coincide con la
    convención del CRB?
16. En el cero exacto sin fondo, ¿el CRB es el **valor puntual** o el **límite**? ¿Un
    regularizador ($p_i+10^{-4}$) o el paso de diferencias finitas cambian el resultado?

**Esquema iterativo y dinámica**
17. ¿La primera iteración o la búsqueda es realista? Un arranque impreciso produce "creeping"
    (SI Fig. 7), y el error de la iteración k es mayor que el CRB centrado.
18. ¿La potencia y la SBR cambian entre iteraciones como en el instrumento?
19. ¿Hay **parpadeo** en la escala del barrido del patrón? Suma $\sigma_{fl}$, que es independiente
    de N ($\mathrm{STD}^2=\sigma_{fl}^2+\sigma_{\rm CRB}^2$).
20. ¿Se consideró el **blanqueo**, restando los fotones emitidos antes del pinhole, y la
    **densidad** de emisores activos, dado que un vecino a ~50 nm sesga hasta ~16 nm?
21. ¿El emisor se mueve durante la medición (difusión o vibración, con dt < dwell)? ¿Se
    incluyeron los **tiempos muertos** del hardware si se simula tracking?
