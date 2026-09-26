# Ronda 4 — PI

## Estado leído
- Ledger r3: núcleo, figuras 1-8, 204 claves y check_provenance están verificados. Suite 177 OK, aceptación 8/8. La prueba de procedencia pasa sin comprobar nada porque el manuscrito está vacío.
- Siguen vivas (refutadas en r3; se resuelven en esta ronda, W1):
  - D1: caché sin invalidación.
  - D2: `--quick` pisa las salidas finales.
  - D3: huecos de check_provenance (`\input` sin llaves, `\src {k}` con espacio, `\pnum` dentro de macros, `\pnum` sin `\src`).
  - D4: el fallback de numbers.tex rompe la compilación.
  - Leyendas de fig8 y fig5. La caveat 0.79 de claims.json.
- Unclear de r3 que se resuelven ahora:
  - Leyendas de fig6 y fig7 escritas como plan (W1).
  - Constantes duplicadas (W1).
  - claims.json incompleto (W1).
  - `misalignment_naive_bias_over_delta_Lq = 0.862±0.012`. Por la decisión del inbox r4 se reemplaza por la estimación poblacional sin ruido (W1) y no se cita.
- Unclear de r1 que se mantienen, con su razón. Van al párrafo "Open points" (writer):
  - Normalización de `montecarlo.rmse` frente a Masullo Ec. 4.2: el texto no está en docs/literature. Se declara la normalización por eje que usamos.
  - La elección SBR=10 para `mle_efficiency_center` la hizo el orquestador en nombre de la autora. Se declara como elección y la autora la confirma al cierre.
- Decisiones del orquestador (inbox r4 y tarea):
  - El manuscrito y el README van en inglés.
  - Se usa revtex4-2 y se compila con Tectonic.
  - `\pnum{k}\src{k}` en toda afirmación cuantitativa.
  - Solo entran afirmaciones VERIFICADAS.
  - El autor es "flor-choque".
  - La metodología agent-team se menciona en Acknowledgments/Methods.

## Decisiones de coordinación (declaradas por el PI)
1. **Nombres de las claves nuevas (poblacionales, sin ruido)**. Las fija el PI para que el writer las use por nombre sin esperar a W1:
   - `misalignment_naive_bias_over_delta_pop_center`: pendiente de ajuste por mínimos cuadrados por el origen sobre δ=2, 5, 10.
   - `misalignment_naive_bias_over_delta_pop_Lq`: ídem en (L/4, 0).
   - `misalignment_naive_bias_over_delta_pop_center_d{2,5,10}` y `misalignment_naive_bias_over_delta_pop_Lq_d{2,5,10}`: la razón |sesgo|/δ para cada δ, que muestra el crecimiento no lineal.
   - `misalignment_pop_n_patterns`: número de patrones, ≥4000.
   - `misalignment_naive_bias_over_delta_Lq` (MC 400×200) se conserva en el JSON como complemento, pero **no se cita**.
   - Valores de referencia del verificador r3: centro ≈0.746/0.755/0.784; L/4 ≈0.794/0.802/0.871 (≈0.83).
2. **Regla nueva de check_provenance**. Cada `\pnum{k}` o `\pnumse{k}` del tex aplanado exige que en el mismo tex exista algún `\src{e}` cuya entrada de provenance.json tenga `k` en `number_keys` o `detail == k`. La regla es global, no por proximidad. Convención del writer: `\pnum{k}\src{k}`, con la entrada `k` = {statement, type:"script", reproduce:"scripts/compute_paper_numbers.py", detail:k, number_keys:[k]}. Las citas `\cite` no llevan `\src`.
3. **Propiedad de archivos**:
   - W1 es dueño de `scripts/**`, `src/**`, `tests/**` (salvo `test_acceptance.py`, intocable), `data/**`, `structure/**`, `paper/generated/numbers.tex` (solo vía compute_paper_numbers), `paper/figures/*.pdf` (solo vía make_all_figures), `README.md`, `CITATION.cff` y `.gitignore`.
   - El writer es dueño de `paper/main.tex`, `paper/sections/**`, `paper/references.bib`, `paper/provenance.json` y `paper/main.pdf`.
   - Ninguno toca los archivos del otro.
4. **Leyendas**. Las dos partes trabajan con la misma especificación de contenido (abajo). El texto de figures.json (W1) y el `\caption` (writer) deben decir lo mismo. La leyenda del writer puede usar `\pnum`. La coherencia la revisa el verificador en r5.
5. **Orden al cerrar la ronda** (lo hace el orquestador tras ambos workers):
   1. `python scripts/compute_paper_numbers.py` (si W1 no lo dejó en sincronía).
   2. Recompilar main.pdf con Tectonic.
   3. `python scripts/check_provenance.py`.
   4. Suite completa y aceptación.
   - Si falta una clave que el writer necesitó, W1 o r5 la agrega.

## Plan de la ronda 4
1. **Worker 1 — herramientas, datos, README.**
   - Cerrar las refutadas D1-D4 con tests.
   - Mover las constantes duplicadas a `_paperconfig`.
   - Agregar las claves poblacionales de desalineación.
   - Corregir las leyendas y completar claims.json.
   - Escribir reproduce.sh/.py, CITATION.cff y el README en inglés.
   - Regenerar todo con `--no-cache` y reportar los números que cambien.
2. **Writer — manuscrito.** main.tex revtex4-2 con secciones, references.bib, provenance.json y main.pdf compilado con Tectonic. Pasa check_provenance. Solo lo verificado; el resto va a "Open points / limitations".

## Especificación de leyendas (común a W1 y al writer)
- **fig5**: 4 paneles. Hay que mencionar el bulto de outliers del MLE en (b): σ/CRB ≈1.5 para x0 = 10-20 nm, N=100, que depende del radio del disco de búsqueda (2L), y hay que declararlo. El sesgo en r=(2,0) sin fondo (≈−0.34 nm) está en el panel **(d)**, no en (c).
- **fig6**:
  - (a) σ en función de N frente a la cámara y a la CRB.
  - (b) σ/σ_cámara, de 0.183 a 0.158.
  - (c) σ por iteración en función de L_k.
  - El punto sin re-centrado es un artefacto del disco de búsqueda; va en gris y rotulado.
  - Redacción descriptiva, nada de "optional" ni "where shown".
- **fig7**: descripción de lo que muestra. El eje y de (c) va en nm. Escala de transición 4.9 nm. La razón 0.78 vale solo para ε ≲ 0.01.
- **fig8**:
  - 400 patrones × 200 repeticiones.
  - La recta dibujada es 0.76·δ, ajustada sobre 7 valores de δ que incluyen δ=15.
  - La pendiente principal es la poblacional sin ruido: ≈0.75 en el centro y ≈0.83 en L/4, creciente con δ.
  - La pérdida del estimador ingenuo es sesgo.
  - El honesto está en el piso MC solo en el centro. En L/4, σ_honesto queda 1-2.4 % sobre la CRB.

## Caveats que el texto debe llevar (ledger r3)
- El exponente en L (1.004) vale solo para L ≪ fwhm: 1.07 sobre 25-150 nm.
- La cola de outliers en x0=15 depende del radio de búsqueda; declarar 2L.
- "Honesto en el piso MC" vale solo en el centro.
- La SE bootstrap a N=250 es inestable (cola pesada).
- La regla adaptativa cubre ≈97 %, no 99 %. La CRB del centro (3.47 nm) subestima el error real (4.26 nm).
- La pendiente iterativa excluye −0.5: −0.515 para N≥500.
- Caprile F=1: 428 frente a 401 nm, no concuerda. Se declara como punto abierto.
- Las cantidades con SE > 2 % se citan con `\pnumse`:
  - `misalignment_*_bias_abs_*`
  - `misalignment_naive_sigma_Lq_d10_nm`
  - `mle_nobg_bias_x_r2_nm`
  - `adaptive_frac_outside_next_radius`

---

## Tareas (autocontenidas)

1. **Worker 1 (tooling / infraestructura / README).**

   Raíz: `C:\Users\BANGHO\Documents\GithubPRO\donut-beam-localization`. Leer antes:
   - `equipo/2026-09-26_donut-localization/reports/r03-code-reviewer.md` (D1-D5)
   - `reports/r03-verifier.md` (X2, X3)
   - este reporte (claves nuevas, regla de procedencia, especificación de leyendas)

   Propiedad exclusiva: `scripts/**`, `src/**`, `tests/**` (NO `tests/test_acceptance.py`), `data/**`, `structure/**`, `paper/generated/numbers.tex`, `paper/figures/*.pdf`, `README.md`, `CITATION.cff`, `.gitignore`. No tocar `paper/main.tex`, `paper/sections/`, `paper/references.bib` ni `paper/provenance.json`.

   (a) **D1**. `_paperstyle.cached` indexa la caché por nombre + hash (sha256 corto) de los parámetros de `_paperconfig` que usa cada cálculo, guardado en el npz y comparado al leer. Test: al cambiar un parámetro se recalcula.

   (b) **D2**. `--quick` escribe a rutas separadas:
   - `paper/figures/quick/…`
   - `data/quick/fig*_summary.json`
   - `data/quick/paper_numbers.json` y `paper/generated/quick/numbers.tex`, o equivalente, con un campo `"quick": true`
   - también `make_all_figures --quick`

   Test: `--quick` no modifica los archivos finales.

   (c) **D3**, en `scripts/check_provenance.py`:
   - aplanar `\input nombre` sin llaves
   - regex tolerante a espacios en `\src {k}` y `\pnum {k}`
   - ignorar `#1` dentro de definiciones de macros
   - regla nueva: cada `\pnum{k}` o `\pnumse{k}` exige un `\src{e}` en el tex cuya entrada de provenance.json tenga `k` en `number_keys` o `detail == k`

   Un test para cada hueco, con copias en `--root` temporal. La regla nueva no debe romper el manuscrito del writer, que usa `\pnum{k}\src{k}` con `detail=k` y `number_keys=[k]`.

   (d) **D4**. El fallback de `numbers_tex` pasa a `\textbf{??\detokenize{#1}??}`. Test de compilación con Tectonic si está disponible (ruta en la tarea del orquestador); si no, un test de string.

   (e) **Constantes a `_paperconfig`**: rb=40000, radios de búsqueda del MLE L/2L, eps=0.002 de `zero_depth_transition_scale_nm`. Las usan compute_paper_numbers.py y fig_5_estimators.py.

   (f) **Pendiente poblacional sin ruido** (decisión del inbox r4):
   - MLE sobre conteos esperados, sin ruido de Poisson.
   - Promedio sobre ≥4000 patrones de desalineación al azar, semilla 42.
   - En el centro y en (L/4,0), para δ=2, 5, 10.
   - Claves exactas:
     - `misalignment_naive_bias_over_delta_pop_center`
     - `misalignment_naive_bias_over_delta_pop_Lq` (pendiente LSQ por el origen sobre δ=2, 5, 10)
     - `misalignment_naive_bias_over_delta_pop_{center,Lq}_d{2,5,10}`
     - `misalignment_pop_n_patterns`
   - Con descripción y `script`. Referencia r3: centro ≈0.75, L/4 ≈0.83.
   - Mantener la clave MC antigua sin citarla.
   - Test unitario: convergencia a δ→0 y comparación con la MC del centro dentro de 3 SE.
   - Si el cálculo es caro, usar caché con el hash de (a).

   (g) **structure/figures.json, en inglés**:
   - Corregir las leyendas de fig5 y fig8 y reescribir las de fig6 y fig7 como descripciones, según la especificación de leyendas de este reporte.
   - Unidad (nm) en el eje y de fig7(c): editar `fig_7_zero_depth.py`.
   - La recta de L/4 en fig8 es opcional; si se agrega, que use las claves nuevas.

   (h) **structure/claims.json**:
   - Cada claim tiene `numbers` no vacío.
   - Cada clave de paper_numbers pertenece a algún claim. La clave antigua `misalignment_naive_bias_over_delta_Lq` va con status `superseded`.
   - Statuses y `verified_round` actualizados desde `state.json`: las pending-r3 que el verificador r3 confirmó pasan a `verified` (round 3).
   - Caveats en inglés con los textos de "Caveats" de este reporte. La caveat 0.79 queda superada.
   - Las claves poblacionales nuevas: `pending-r4-verification`.
   - Agregar a check_provenance (o a un test) la comprobación de cobertura en ambas direcciones.

   (i) **Scripts de reproducción**:
   - `scripts/reproduce.sh` y `scripts/reproduce.py`, multiplataforma, con la secuencia tests → make_all_figures `--no-cache` → compute_paper_numbers → check_provenance → aceptación.
   - Opción `--quick`, y compilación opcional de LaTeX si hay `tectonic` en PATH.
   - Además, `CITATION.cff`: autor flor-choque, GitHub FlorJDC; sin afiliación inventada.

   (j) **README.md en inglés**, que reemplaza la plantilla de agent-team. Mover la plantilla a `docs/` solo si tiene valor; si no, eliminarla. Estilo del repo pta-gwb-anisotropy. Contenido:
   - título y resumen de una línea
   - quick start (install, tests, reproduce)
   - estructura del repo
   - tabla del pipeline (script → salida)
   - tabla de figuras (id, script, pregunta R1-R5, qué muestra)
   - `data/` y paper_numbers
   - testing y aceptación
   - "How this study was made": metodología agent-team, verificación independiente, trazabilidad en `equipo/2026-09-26_donut-localization/` (state.json, reports) y regla de procedencia
   - licencia (LICENSE existente)
   - sin nada de docs/private

   (k) **Regenerar** todo con `make_all_figures.py --no-cache` y después `compute_paper_numbers.py`. Correrlo en segundo plano si tarda.
   - Diff de paper_numbers.json antes y después: listar cada clave cuyo valor cambie (con los valores viejo y nuevo) y las claves nuevas.
   - Verificar que los summaries coinciden con paper_numbers.
   - Correr `python -m unittest discover -s tests`, `python -m unittest tests.test_acceptance -v` y `python scripts/check_provenance.py`, y reportar la salida real.
   - Confirmar el sha256 de test_acceptance.py = 7d198853…2303.

   Reporte en `reports/r04-worker-1.md`, con un bloque `claims` de afirmaciones verificables.

2. **Writer (manuscrito).**

   Raíz: la misma. Leer antes:
   - `state.json`: solo las claims `verified`. Las refuted/unclear no entran salvo como punto abierto.
   - `inbox.jsonl` (todas las decisiones)
   - `structure/claims.json` y `structure/figures.json`
   - claves de `data/paper_numbers.json`
   - `docs/derivations/crb_tcp_center.md`
   - `docs/literature/A_minflux_theory.md` y `B_donut_optics.md`
   - `papers/README.md`
   - este reporte (claves nuevas, caveats, especificación de leyendas, regla de procedencia)
   - NADA de `docs/private/`

   Propiedad exclusiva: `paper/main.tex`, `paper/sections/*.tex`, `paper/references.bib`, `paper/provenance.json`, `paper/main.pdf`. No tocar scripts, data, structure, numbers.tex ni las figuras.

   Contenido:
   - `main.tex` en revtex4-2, que conserva `\input{generated/numbers}` y `\newcommand{\src}[1]{}`.
   - Autor "flor-choque", sin afiliación inventada.
   - `\input` con llaves de cada sección:
     - `abstract`
     - `introduction`
     - `model`: haces LG, TCP, estadística multinomial de fotones con fondo y definición de SBR.
     - `crb_center`: valor puntual S27 frente al límite r→0, forma cerrada, razón 2/√5, elipse límite, escalamientos en L y N con su caveat, SBR.
     - `vectorial`: handedness, profundidad del cero 0/0.845/0.372, LG con curvatura igualada, CRB vectorial frente a LG.
     - `estimators`:
       - MLE eficiente con SBR=10
       - superficiencia sin fondo (0.84 respecto del límite, sesgo hacia el centro, modelo no regular)
       - LMS igual a S27 en el centro, sesgo de LMS/mLMS
       - cola de outliers fuera del centro a N=100, con el radio de búsqueda 2L declarado
     - `camera`: ideal, pixelada, fondo por píxel frente al total.
     - `iterative`: σ(N), pendiente −0.515 que excluye −0.5, razón con la cámara 0.183→0.158, regla adaptativa ≈97 %, artefacto sin re-centrado.
     - `nonidealities`:
       - profundidad finita del cero y L óptimo ∝ fwhm√ε (0.78 solo para ε ≲ 0.01)
       - divergencia a 360 nm
       - desalineación honesta frente a ingenua: pendiente POBLACIONAL por nombre de clave, MC con `\pnumse`, la pérdida es sesgo, honesto en el piso solo en el centro
     - `discussion`: qué importa experimentalmente, ordenado por impacto numérico.
     - `open_points`: limitaciones y todo lo no verificado:
       - normalización de rmse frente a Masullo Ec. 4.2
       - Caprile F=1
       - SE bootstrap a N=250
       - SBR=10 como elección declarada
       - claves `pending-r4-verification` si el verificador no las confirma
     - `methods`: reproducibilidad, semilla 42, `reproduce` y check de procedencia. Una o dos frases de que el estudio se hizo con un equipo de agentes (agent-team) con verificación independiente y trazabilidad en `equipo/`.
     - `appendix`: derivación de `docs/derivations/crb_tcp_center.md`.
     - Acknowledgments.
   - Toda afirmación cuantitativa lleva `\pnum{k}\src{k}`, o `\pnumse{k}` para las claves con SE > 2 % (lista en "Caveats"). Usar las claves por NOMBRE. Nunca un número escrito a mano, salvo constantes matemáticas exactas de una derivación (2/√5, 1/8, 1/10) en el apéndice.
   - `paper/provenance.json`: una entrada por clave usada = {statement (en inglés), type:"script", reproduce:"scripts/compute_paper_numbers.py", detail:k, number_keys:[k]}. Pueden sumarse entradas type "derivation" (reproduce: `docs/derivations/crb_tcp_center.md`) para el apéndice.
   - Figuras: `\includegraphics` de `figures/fig1_schematic.pdf` … `fig8_misalignment.pdf`, con leyendas en inglés que cumplan la especificación de leyendas de este reporte y coincidan con figures.json.
   - `references.bib`: solo las obras citadas, con título, autores y año tomados de `papers/README.md` y de las notas de literatura (no inventar DOIs; si no están, se omiten).
   - Compilar con `C:\Users\BANGHO\AppData\Local\Temp\claude\C--Users-BANGHO-Documents-GithubPRO\36219331-e4c3-4b34-80d5-6a75c445dbc9\scratchpad\tectonic\tectonic.exe paper/main.tex`: sin `??`, sin referencias indefinidas, y dejar `paper/main.pdf`.
   - Correr `python scripts/check_provenance.py` hasta "all checks pass". Las claves poblacionales nuevas pueden faltar mientras W1 no termine: listarlas como esperadas, sin reemplazarlas por números.
   - Si hace falta una clave que no existe en paper_numbers.json, no escribir el número: listarla en el reporte con su definición propuesta.

   Reporte en `reports/r04-writer.md`:
   - secciones escritas
   - claves usadas
   - claves faltantes
   - salida de Tectonic y de check_provenance
   - bloque `claims`

## LO QUE MÁS SE PODRÍA HACER
1. **R5**, revisión independiente del manuscrito completo (verifier y code-reviewer): cada `\pnum` contra el ledger, las leyendas contra los PDFs y las claves poblacionales nuevas contra el valor r3. Imprescindible antes de `[[DONE]]`. Costo: 1 ronda.
2. **Confirmación explícita de la autora** de SBR=10 como definición de `mle_efficiency_center`: es una decisión de juicio, no técnica. Costo: minutos.
3. **Normalización de rmse frente a Masullo Ec. 4.2**: extraer la ecuación del PDF a docs/literature y cerrar la unclear de r1. Costo: bajo.
4. **Caprile F=1 (428 frente a 401 nm)**: comparar con PyFocus como verificación externa y cerrar la discrepancia del diámetro del anillo. Costo: medio.
5. **Cola pesada a N=250 en el iterativo**: más repeticiones o un estimador robusto de σ, para dar una SE estable. Costo: medio (cómputo).
6. **Dona 3D y desalineación axial**: fuera del alcance de R1-R5. Una extensión natural del paper. Costo: alto.
