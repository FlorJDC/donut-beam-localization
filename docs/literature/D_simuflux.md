# D. SimuFLUX: simulación realista de experimentos MINFLUX — nota para este proyecto

**Fuente leída**: Z. Marin y J. Ries, "Evaluating MINFLUX experimental performance in silico",
*Nat. Commun.* **17**, 246 (2026), doi:10.1038/s41467-025-66952-w, con su Información
Suplementaria (SI) y el Peer Review File. Clave: [Marin2026]. El PDF no se sube al repo.
"Fig. 2x" es la figura del artículo principal; "SI Fig. n", la de la Información Suplementaria.

## Qué es
SimuFLUX es un simulador de código abierto (MATLAB y Python; numpy+scipy; GPL-3.0) que
reproduce un experimento MINFLUX completo: PSF (analítica o vectorial tipo Leutenegger, con
pinhole, aberraciones, desalineación, cero imperfecto, tamaño de bead), fluoróforos (switching,
flickering, bleaching, movimiento, varios fluoróforos, fondo), microscopio (escáner lento y EOD,
tiempos muertos) y estimadores (LSQ linealizado, LSQ iterativo, solución directa, MLE, con y sin
fondo). Puede leer archivos de secuencia de Abberior [Marin2026, pp. 1–2, 4–5].

## Resultados que importan aquí
- **El CRB no captura las imperfecciones.** Desalinear la máscara de fase unos 100 µm produce
  sesgo dependiente de z y más STD, sin que el CRB lo refleje (Fig. 2a; SI Fig. 1). El pinhole
  desalineado casi no afecta, salvo la eficiencia de detección (SI Fig. 2).
- **Aberraciones y fondo** pesan mucho más en la dona 3D que en la 2D o en la half-moon
  (Fig. 2b; SI Fig. 3). Las vibraciones agregan un error del orden de su amplitud (SI Fig. 4).
- **Fondo → sesgo.** Fondo de cero imperfecto, autofluorescencia o un fluoróforo vecino no solo
  empeora la precisión: sesga fuertemente el estimador. Estimar y restar el fondo, o ajustarlo como
  parámetro libre, elimina el sesgo (Fig. 2d, 200 fotones de señal + 200 de fondo; SI Fig. 5).
- **Estimadores.** Proponen un LSQ iterativo y una solución directa en 1D que amplían el campo
  sin sesgo respecto del estimador lineal estándar. Un sesgo fuerte puede dar una STD baja pero sin
  sentido. Con pocos fotones (15) el MLE es el más sesgado (Fig. 2c; SI Fig. 6).
- **Flickering.** En MINFLUX secuencial, el parpadeo agrega un error σ_fl que no depende de N:
  STD(r, N)² = σ_fl(r)² + σ_CRB(N)², con r el número de repeticiones del patrón. Puede llegar a
  ~10 nm y exige muchas repeticiones (> 20) para mitigarlo (Fig. 2e; SI Fig. 8; t_on = t_off =
  100 µs).
- **Experimentos completos.** DNA-PAINT con hebras que difunden genera un fondo variable y
  localizaciones inválidas, y las hebras "quenched" lo mitigan. Los PAFP se blanquean durante el
  scouting. En tracking, los tiempos muertos del hardware limitan el D máximo (optimizado de 2.5 a
  4.2 µm²/s) (Fig. 2f–l; SI Figs. 9–10).

## Relación con este repo y con p-MINFLUX
- La SI Tabla 1 compara simuladores y cita p-MINFLUX (Masullo et al., *Nano Lett.* **21**, 840,
  2021) y el marco común (Masullo et al. 2022).
- SimuFLUX está centrado en el MINFLUX **secuencial e iterativo** de Abberior. No modela la
  excitación **pulsada entrelazada** de p-MINFLUX: la asignación de fotones por ventanas temporales
  (TCSPC), el cruce entre ventanas por el tiempo de vida del fluoróforo, ni el hecho de que las 4
  exposiciones se muestrean de forma casi simultánea. Esos efectos se estudian en este repo
  (`src/donutloc/pminflux.py`, `scripts/fig_9_pminflux_timing.py`).
- El sesgo por fondo no modelado se estudia aquí para el TCP fijo con MLE
  (`src/donutloc/background.py`, `scripts/fig_10_background_bias.py`). Es el análogo cualitativo
  de su Fig. 2d, no una reproducción numérica (ellos usan LSQ con 200 + 200 fotones).

## Ver también
`C_insilico_minflux.md` (resumen detallado de SimuFLUX con claves [MR-*]) y `C_insilico_vs_donutloc.md`
(comparación numérica y checklist) cubren el mismo paper con más detalle. Esta nota D se centra en su
relación con p-MINFLUX.
