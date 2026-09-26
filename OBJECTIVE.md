OBJETIVO
========
Construir un estudio de simulación reproducible, al estilo de un "paper companion", sobre la
**localización de emisores fluorescentes individuales con haces estructurados tipo dona**
(MINFLUX y variantes). Al terminar debe existir un repositorio con un paquete Python
(`src/donutloc`), scripts que generan cada figura, tests, un registro de números
(`data/paper_numbers.json`) y un manuscrito LaTeX (`paper/`) en el que cada afirmación
cuantitativa está respaldada por un script o un check ejecutable.

El manuscrito responde, con simulaciones verificadas, cinco preguntas:

R1. **Modelo del haz.** ¿Cuánto difiere la dona escalar paraxial (Laguerre-Gauss LG01) de la
    dona vectorial de alta apertura numérica (Richards-Wolf, máscara de vórtice + polarización
    circular) cerca del cero, y cuánto cambia eso la precisión de localización (CRB)? ¿Qué
    pasa con la handedness incorrecta o con polarización lineal (profundidad del cero)?
R2. **Límite de Cramér-Rao.** Mapas de CRB del TCP (3 donas en un círculo de diámetro L +
    centro) en el campo de visión; forma cerrada en el centro frente al cálculo numérico;
    escalamiento con L, N y SBR; comparación con la localización por cámara (σ_PSF/√N).
R3. **Estimadores.** MLE frente a LMS/mLMS: sesgo y precisión por Monte Carlo, comparados
    con el CRB, dentro y fuera del TCP.
R4. **MINFLUX iterativo.** Precisión en función del presupuesto total de fotones al reducir
    L por iteraciones, frente a la localización por cámara.
R5. **No idealidades.** Profundidad finita del cero, fondo (SBR) y errores de posición de las
    donas (desalineación del TCP): sesgo y pérdida de precisión.


FUENTES
=======
Los PDFs están en `C:\Users\BANGHO\Documents\Doctorado\Papers` (NO se suben al repo por
copyright; `papers/README.md` lista las referencias). Notas extraídas, con ecuaciones y
páginas, en `docs/literature/`:
- `A_minflux_theory.md`: modelo probabilístico, estimadores, Fisher/CRB, iterativo
  (Balzarotti 2017; Masullo 2022 marco común; Masullo & Stefani 2022; Stefani 2023).
- `B_donut_optics.md`: dona LG y vectorial, dona 3D, imperfecciones (López 2023; Caprile
  2022 PyFocus; Tarkowski; Gwosch 2020).
- `docs/private/C_pminflux_practice.md` (LOCAL, no versionado): p-MINFLUX, parámetros
  experimentales y el trabajo previo no publicado de la autora. Puede orientar prioridades, pero
  nada de su contenido específico (mediciones, resultados, planes) se copia al repo público.


CERCO -- no leer nada más
==========================
No hay cerco estricto. No copiar código de terceros (PyFocus, repos de MINFLUX): la
implementación es propia, y los paquetes externos solo se usan, si acaso, como verificación
independiente.


LO QUE LAS FUENTES NO DICEN
===========================
Decisiones que el equipo toma y declara en el manuscrito:
- Convención de la dona LG: I(r) = 4 e ln2 (r²/fwhm²) exp(-4 ln2 r²/fwhm²) (pico = 1),
  con fwhm como parámetro de tamaño (fwhm_nm registrado en `data/paper_numbers.json`,
  entre 200 y 500 nm; default sugerido 300 nm).
- TCP de referencia: K = 4 exposiciones, 3 donas equiespaciadas en un círculo de diámetro L
  más una en el centro. Probabilidades p_i = I_i / Σ_j I_j sin fondo; con fondo, mezcla con
  un término uniforme según la SBR (definir exactamente qué SBR se usa).
- CRB reportado como σ_CRB = sqrt((Σ⁻¹_xx + Σ⁻¹_yy)/2), con Σ = F la información de Fisher
  para N fotones multinomiales: F = N Σ_i ∇p_i ∇p_iᵀ / p_i.
- En el centro exacto del TCP la exposición central tiene p = 0 (cero perfecto, sin fondo):
  el CRB "en el centro" se define como el límite r → 0 (es isótropo para el TCP simétrico).
- La eficiencia del MLE que exige la definición de terminado (`mle_efficiency_center`) se mide
  con SBR = 10: sin fondo el modelo no es regular en el centro (p_centro ∝ r²) y el MLE resulta
  sesgado y "superficiente"; ese comportamiento se reporta aparte como resultado.
- Parámetros ópticos por defecto: λ = 640 nm, NA = 1.4, n = 1.518 (ajustables).
- Monte Carlo: semilla fija (42) registrada en cada salida; número de repeticiones suficiente
  para que el error estadístico de cada número citado sea < 2 %.


DISCIPLINA
==========
- Entorno: Python ≥ 3.8 con numpy, scipy, matplotlib (lo único instalado en la máquina de
  referencia; numba está disponible pero es opcional). Nada de dependencias nuevas
  obligatorias. Tests con `unittest` (compatibles con pytest).
- Cada figura tiene su script `scripts/fig_*.py`; `scripts/make_all_figures.py` las genera
  todas; `scripts/compute_paper_numbers.py` es el ÚNICO que escribe `data/paper_numbers.json`.
- Unidades: nm en todo el código público.
- Todo resultado que entra al manuscrito debe haber sido reproducido por el verificador.
- No modificar `tests/test_acceptance.py` (protegido por hash).


DEFINICIÓN DE TERMINADO
=======================
`python -m unittest tests.test_acceptance` pasa. Ese test recalcula de forma independiente
el CRB en el centro del TCP y compara con `data/paper_numbers.json`, verifica los
escalamientos, la eficiencia del MLE, la profundidad del cero vectorial, que cada figura
declarada en `structure/figures.json` exista junto con su script, y que el chequeo de
procedencia del manuscrito pase. Además, la suite de tests unitarios
(`python -m unittest discover -s tests`) pasa completa.


ENTREGABLE
==========
- `src/donutloc/`: beams (LG y vectorial), patrones (TCP), modelo de fotones y fondo,
  Fisher/CRB, estimadores (MLE, LMS), Monte Carlo, MINFLUX iterativo.
- `scripts/`: `fig_*.py`, `make_all_figures.py`, `compute_paper_numbers.py`,
  `check_provenance.py`, `reproduce.sh`.
- `tests/`: unitarios + `test_acceptance.py`.
- `data/paper_numbers.json`: cada número citado en el paper, con el script que lo produce.
- `structure/figures.json` y `structure/claims.json`: registros de figuras y afirmaciones.
- `paper/main.tex` + `paper/sections/*.tex` + `paper/figures/*.pdf` + `references.bib`;
  cada afirmación cuantitativa con `\src{clave}` resuelta en `paper/provenance.json`
  (el chequeo de `agent-team/bin/check_provenance.py` debe pasar).
- `README.md` con quick start, estructura y tabla de figuras (como pta-gwb-anisotropy).
