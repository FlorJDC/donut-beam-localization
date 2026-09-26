# Ronda 1 — PI (plan)

## Estado leído
- `state.json`: ronda 0, sin afirmaciones, sin plan ni backlog. No hay afirmaciones `refuted` ni `unclear` pendientes.
- Proyecto: todavía no existen `src/`, `data/`, `structure/` ni `paper/`. El test de aceptación da SKIP porque falta `data/paper_numbers.json`.
- `inbox.jsonl` tiene 3 indicaciones humanas, todas incorporadas al plan:
  1. La "CRB en el centro" es el **límite r→0**, y Eq. S27 es el valor puntual en r=0. Hay que verificar el cociente 2/√5 y reportar las dos cosas: esto es la tarea 3.
  2. La profundidad finita del cero y la desalineación del TCP con estimador ingenuo son prioridad en R5. Esta ronda ya deja el parámetro `eps` en el haz, `perturb_centers` en los patrones y la forma cerrada con `eps`.
  3. `fwhm` es el de Eq. S17, con default 300 nm. La tabla de Masullo se reproduce con ~360 nm: se anota en la tarea 3 y en el backlog para el paper.
- **Regla de privacidad (nota del orquestador; el repo es público).** La nota C está ahora en `docs/private/` (gitignored). Nada específico de ella (mediciones, resultados, planes o código de la autora) se copia a archivos versionados: código, docs, paper ni reportes. Solo sirve para fijar prioridades. Toda justificación pública sale de la literatura publicada (Balzarotti 2017, Masullo 2021 y su tesis, Stallinga 2026, nota A) o de nuestras simulaciones. Los workers son agentes genéricos que siguen `.claude/agents/worker.md` y deben respetar esta regla.

## Chequeo rápido (lo hice a mano; no es una afirmación verificada)
En el caso cuadrático sin fondo:
- La información de Fisher puntual en el centro es F = (8N/L²)·I.
- En el límite r→0 el término central agrega (16N/(3L²))·r̂r̂ᵀ.
- Acercándose por x: var_x = 3L²/(40N) y var_y = L²/(8N). El promedio es 0.1·L²/N, que es isótropo.
- Resulta σ_lim/σ_S27 = √0.8 = 2/√5.

Esto coincide con el ejemplo de la nota A (σ_x = 2.739 nm con L=100 y N=100). Falta el caso dona (fwhm finito). No sé si el factor (1−L²ln2/fwhm²)⁻¹ se aplica igual al término central; eso lo tiene que derivar la tarea 3.

## Decisiones declaradas (defaults; las tres tareas deben respetarlas)
- **Unidades:** nm. Las posiciones se pasan como arrays de forma `(..., 2)`.
- **TCP:** 3 donas en un círculo de **diámetro L**, en los ángulos `rotation + 2πk/3` con `rotation = π/2` (vértice arriba), y la dona central **al final** (índice 3). Es exactamente la convención de `tests/test_acceptance.py::_tcp`. Balzarotti numera 0 = centro; donde haga falta se documenta la equivalencia.
- **Dona LG:** `I = 4e ln2 r²/fwhm² · exp(-4 ln2 r²/fwhm²)`, con pico 1 en `r_ring = fwhm/(2√ln2)`.
- **Cero residual `eps`:** es la intensidad relativa al pico del anillo en r=0. Hay dos modelos:
  - `zero_model="gaussian"` (default): `I_LG + eps·exp(-4ln2 r²/fwhm²)`;
  - `zero_model="constant"`: `I_LG + eps`.
  - No se renormaliza, porque la escala global se cancela en p_i. Así, "x % del pico en el cero" corresponde a `eps` = x/100.
- **Cuadrático:** `quadratic = 4e ln2 r²/fwhm²`, la curvatura de la dona en el cero (Eq. S20). Con `fwhm` → curvatura se recupera el límite pequeño.
- **Fondo:** Eq. S30 con SBR fija, que es condicionar en N (default). La alternativa es `bg_per_exposure` fijo, con SBR dependiente de la posición: `SBR = Σλ/(K·λ_b)` (Eq. S29).
- **CRB:** `σ = sqrt((Σxx+Σyy)/2)` con `F = N Σ ∇p ∇pᵀ/p`.
  - Valor **puntual**: los términos con `p_i < p_min` (default 1e-12) se excluyen y queda documentado. Es la convención de S27.
  - **Límite central**: promedio sobre 12 direcciones a r0 = 1e-3 nm, igual que en el test de aceptación.
- **Pegamento entre módulos:** un modelo es un callable `p_fn(r)`. Recibe `r` de forma `(...,2)` en nm, devuelve `(...,K)` con suma 1 en el último eje y **debe estar vectorizado**. `fisher`, `estimators` y `montecarlo` solo conocen `p_fn`. Por eso las tres tareas pueden correr en paralelo sin importarse entre sí.
- **Tests:** unittest. Cada archivo `tests/test_*.py` agrega `src` a `sys.path` con este encabezado de 3 líneas, sin helper compartido para evitar choques de propiedad:
  ```python
  import os, sys
  sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
  ```
  Tienen que pasar con `python -m unittest discover -s tests` desde la raíz.
- **Entorno:** Python 3.8. Eso excluye `X | Y` en anotaciones y `match`. Solo numpy y scipy. Los archivos con caracteres no ASCII llevan `# -*- coding: utf-8 -*-`. Semilla 42.
- **Prohibido:** tocar `tests/test_acceptance.py` y copiar código de terceros.

## Tareas (una por worker; propiedad de archivos exclusiva)

1. **Worker A — núcleo físico: haces, patrones, fotones.**
   Archivos propios: `src/donutloc/__init__.py`, `src/donutloc/beams.py`, `src/donutloc/patterns.py`, `src/donutloc/photons.py`, `pyproject.toml` (mínimo, `src`-layout, sin dependencias nuevas), `tests/test_beams.py`, `tests/test_patterns.py`, `tests/test_photons.py`.
   - `__init__.py` solo contiene el docstring y `__version__ = "0.1.0"`, sin imports de submódulos: los otros workers escriben módulos hermanos en paralelo.
   - `beams.py`:
     - `lg_donut(x, y, fwhm=300.0, eps=0.0, zero_model="gaussian")`, `gaussian(x, y, fwhm=300.0)` y `quadratic(x, y, fwhm=300.0)`, con las definiciones de "Decisiones".
     - `ring_radius(fwhm)`.
     - `make_beam(kind="donut", fwhm=300.0, eps=0.0, zero_model="gaussian", power=1)`, que devuelve `f(x, y)` y admite `power` como exponente multifotón (`I**power`).
   - `patterns.py`:
     - `tcp_centers(L, rotation=np.pi/2, center=True)`, que devuelve `(4,2)` con el centro al final.
     - `polygon_centers(L, M, rotation=np.pi/2, center=True)`.
     - `perturb_centers(centers, displacement, rng=None, mode="random_direction")`, que desplaza cada dona `displacement` nm en una dirección al azar (modelo de Masullo para la desalineación). También acepta un array `(K,2)` de desplazamientos explícitos.
   - `photons.py`:
     - `intensities(r, centers, beam)`, que devuelve `(...,K)`.
     - `probabilities(r, centers, beam, sbr=None, bg_per_exposure=None)`. `sbr=None` o `inf` significa sin fondo. Se usa Eq. S30 o S28/S29, y pasar los dos argumentos da `ValueError`.
     - `sbr_at(r, centers, beam, bg_per_exposure)`.
     - `make_model(centers, beam, sbr=None, bg_per_exposure=None)`, que devuelve `p_fn` según el contrato de "Decisiones".
     - `sample_counts(p, N, size=None, rng=None, mode="multinomial")`. Con `"poisson"`, `N` es la media total. Devuelve enteros `(size,K)`.
     - Todo debe estar vectorizado sobre `r` de forma `(...,2)`.
   - Tests mínimos:
     - pico 1 en `ring_radius` y diámetro pico a pico ≈ 1.2·fwhm;
     - `lg_donut(0,0,eps=e) == e`;
     - `quadratic` coincide con `lg_donut` a r ≪ fwhm;
     - `tcp_centers(L)` coincide exactamente con `_tcp(L)` del test de aceptación (copiar la fórmula al test, sin importar el test de aceptación);
     - Σp = 1, y p_centro = 0 en r=0 sin fondo;
     - con `sbr`, p_i(0) cumple Eq. S30;
     - pedestal constante ≡ fondo efectivo (derivación propia a verificar: `1/SBR_eff = 1/SBR + K·eps/ΣI`, comparar Eq. S28–S30);
     - la media de `sample_counts` converge a N·p (con semilla 42).
   - Reportar: la API final (firmas) y las afirmaciones verificables en un bloque `claims`.

2. **Worker B — inferencia: Fisher/CRB numérico, estimadores, Monte Carlo.**
   Archivos propios: `src/donutloc/fisher.py`, `src/donutloc/estimators.py`, `src/donutloc/montecarlo.py`, `tests/test_fisher.py`, `tests/test_estimators.py`, `tests/test_montecarlo.py`.
   **No importar** `beams`, `patterns` ni `photons`, que los escribe otro worker en paralelo. Todo trabaja sobre un callable `p_fn(r)`: recibe `(...,2)` nm y devuelve `(...,K)` con suma 1, vectorizado. En los tests hay que definir un modelo LG+TCP propio mínimo de ~10 líneas:
   - LG `4e ln2 r²/fwhm² exp(-4ln2 r²/fwhm²)`;
   - centros en ángulos `π/2+2πk/3` sobre un radio L/2, más el centro al final;
   - SBR opcional según Eq. S30.
   - `fisher.py`:
     - `fisher_matrix(p_fn, r, N, h=1e-3, p_min=1e-12)`, que devuelve `(...,2,2)`. Usa diferencias centradas. Los términos con `p_i<p_min` se excluyen: es el valor puntual, convención S27, y debe quedar documentado.
     - `crb(p_fn, r, N, **kw)`, que devuelve σ = sqrt(tr(F⁻¹)/2), de forma `(...)`.
     - `crb_axes(p_fn, r, N, **kw)`, que devuelve `(σx, σy, isotropía)`.
     - `crb_limit(p_fn, N, r_center=(0,0), r0=1e-3, n_dir=12)`: el promedio de `crb` sobre `n_dir` direcciones a distancia r0 (`h = r0·1e-2`). Es la definición del test de aceptación.
     - `crb_map(p_fn, xs, ys, N)`.
   - `estimators.py`:
     - `neg_loglike(r, counts, p_fn)`, es decir −Σ nᵢ ln pᵢ con p recortado a 1e-300.
     - `mle(counts, p_fn, search_radius, center=(0,0), grid_step=None, refine=True)`. Acepta `counts` de forma `(K,)` o `(M,K)` y devuelve `(2,)` o `(M,2)`. Primero una grilla vectorizada en un disco (paso por defecto `search_radius/25`) y después un refinamiento continuo con `scipy.optimize.minimize` (Nelder-Mead o L-BFGS-B) acotado al disco. Tiene que ser robusto con n_centro = 0 y con N chico, y rápido: unas 2000 repeticiones en menos de ~30 s. Hay que documentar la rama elegida cuando hay máximos múltiples.
     - `lms(counts, p_fn, r_lin=(0,0), h=1e-3)`: la forma general de Eq. S48 con el Jacobiano numérico completo de `p_fn`. Si `p_fn` incluye fondo, la corrección 1/s (nota A §3.2) sale sola.
     - `lms_tcp(counts, L, fwhm, sbr=None, rotation=np.pi/2)`: Eq. S49–S50 en nuestra convención (centro al final), con el factor 1/s si hay `sbr`.
     - `mlms_tcp(counts, L, fwhm, beta=(1.27, 3.8), sbr=None, rotation=np.pi/2)`, que implementa Eq. S51.
   - `montecarlo.py`:
     - `run_mc(estimator, p_fn, r_true, N, n_rep, seed=42, mode="multinomial")`. Muestrea con `rng.multinomial` directamente. `estimator` es un callable `counts(M,K) -> (M,2)`.
     - Devuelve un dict con `bias` (2,), `std` (2,), `sigma = sqrt((var_x+var_y)/2)`, `rmse` (Eq. 4.2 de Masullo, respecto de la posición verdadera), `sigma_err` (error estándar ≈ sigma/sqrt(2·n_rep)), `n_rep`, `seed` y `estimates`.
   - Tests:
     - `crb` puntual en el centro = Eq. S27, `L/(2√(2N))·(1−L²ln2/fwhm²)⁻¹`, con rtol 1e-4. Se prueba con L=50, N=100 y fwhm=300.
     - Con SBR=10, el centro da Eq. S31.
     - `crb_limit` sin fondo coincide con `_crb_center` del test de aceptación (reimplementado) a 1e-3.
     - Escala exacta ∝ N^-1/2 y ≈ ∝ L para L ≪ fwhm.
     - `mle` recupera la posición verdadera en datos sin ruido (conteos esperados enteros grandes).
     - `lms_tcp` = `lms` general.
     - Monte Carlo en el centro, L=50, N=100, SBR=inf: std(MLE)/`crb_limit` entre 0.9 y 1.2 con n_rep ≥ 2000 (en el test chico basta n_rep=500 y tolerancia amplia; reportar el número grande en el reporte).
   - Reportar: la API final, la eficiencia MC medida con su error, los tiempos de ejecución y los claims.

3. **Worker C — derivación analítica del CRB en el centro del TCP + formas cerradas + script de verificación.**
   Archivos propios: `docs/derivations/crb_tcp_center.md`, `src/donutloc/closed_forms.py`, `tests/test_closed_forms.py`, `scripts/verify_crb_closed_forms.py`.
   Hay que derivar a mano (sympy solo si está instalado, no se agrega como dependencia) para el TCP con L = diámetro, centro al final y dona LG `4e ln2 r²/fwhm² exp(-4ln2 r²/fwhm²)`:
   - **(a)** El valor puntual en r=0 exacto, excluyendo el término con p_centro = 0. Tiene que dar Eq. S27, `L/(2√(2N))·(1−L²ln2/fwhm²)⁻¹`.
   - **(b)** El límite r→0 sin fondo, que es la definición del test de aceptación: el término central `(∇p₀)²/p₀` depende de la dirección, y luego se toma el promedio isótropo de `sqrt(tr/2)`. Hay que verificar la hipótesis humana: en el límite cuadrático, σ²_lim = 0.1·L²/N frente a L²/(8N), con cociente **2/√5**. Después se obtiene la forma cerrada para fwhm finito: ¿se conserva el cociente o se corrige con un factor que dependa de L²ln2/fwhm²? Hay que explicar la discontinuidad en r=0.
   - **(c)** El fondo con SBR fija (Eq. S31). Es continuo y no hay discontinuidad. Hay que discutir por qué SBR→∞ en S31 da el valor **puntual** (a) y no el límite (b): los límites no conmutan.
   - **(d)** La profundidad finita del cero `eps` en los dos modelos:
     - `zero_model="constant"`: `I_LG + eps`;
     - `zero_model="gaussian"`: `I_LG + eps·exp(-4ln2 r²/fwhm²)`. Aquí `eps` es la intensidad relativa al pico del anillo en r=0.
     - Forma cerrada (exacta para "constant"; para "gaussian", exacta o con el orden de error declarado) usando la equivalencia pedestal ≡ fondo efectivo: `1/SBR_eff = 1/SBR + K·eps/ΣI_j(0)`.
     - Hay que dar una tabla del factor de degradación CRB(eps)/CRB(0) para eps ∈ {0.002 (Balzarotti p. 6), 0.01, 0.03, 0.05, 0.1, 0.15}, L ∈ {50, 100, 150} y fwhm ∈ {300, 360}. Solo con valores genéricos: no citar datos no publicados.
   - **(e)** Opcional si queda tiempo: el multifotón con exponente c, `L/(2c√(2N))` (nota A §6).
   - `closed_forms.py` (solo numpy):
     - `crb_tcp_center_point(L, N, fwhm=np.inf, sbr=np.inf)`, para S27 y S31;
     - `crb_tcp_center_limit(L, N, fwhm=np.inf)`, para (b);
     - `crb_tcp_center_eps(L, N, fwhm, eps, sbr=np.inf, zero_model="constant")`;
     - `sbr_center_vs_L(L, L0, sbr0, fwhm)`, para Eq. S32;
     - `crb_1d_center(L, N, fwhm=np.inf, kind="donut"|"quadratic"|"gaussian")`, para S22c, S22f y S23c.
   - Tests (`tests/test_closed_forms.py`): hay que comparar contra una Fisher numérica **propia inline**, sin importar `donutloc.fisher`, que se escribe en paralelo:
     - valor puntual;
     - límite (con el promedio de 12 direcciones a r0 = 1e-3, como el test de aceptación);
     - SBR ∈ {5, 10};
     - eps ∈ {0.01, 0.1};
     - L ∈ {50, 100}, fwhm ∈ {300, 360}.
     - Además, los checks numéricos de la nota A §8.2: 3.536 nm; 1.792 nm; 1.817 nm; 0.941/1.962/3.167 nm con fwhm=360 (tabla de Masullo) frente a 0.947/2.012/3.370 nm con 300; 10.82 nm.
   - `scripts/verify_crb_closed_forms.py` (con `# -*- coding: utf-8 -*-`):
     - imprime una tabla "closed form | numerical | rel. error" para todos los casos, con su propia Fisher numérica;
     - **si** `donutloc.fisher`, `donutloc.photons` y `donutloc.patterns` son importables, agrega una columna con la ruta del paquete (con `try/except ImportError`);
     - termina con exit code 1 si algún error relativo supera 1e-3 (límite) o 1e-6 (puntual).
     - Imprime también el número que el test de aceptación va a exigir: `crb_center_lg_L50_N100_nm` con fwhm=300 (límite r→0).
   - El `.md` tiene que estar autocontenido: pasos algebraicos, supuestos, resultados en caja y una tabla numérica.
   - Reportar las fórmulas finales y los números en un bloque `claims`.

## Backlog (no entra en esta ronda; nada se descarta)
- **R2**
  - Integración: un test que construya `p_fn` con `photons.make_model` y compare `fisher.crb_limit` con `closed_forms.crb_tcp_center_limit` y con el `_crb_center` de aceptación.
  - Si W3 lo deja como función, pasar el límite a `fisher`.
  - CRB de cámara ideal (Eq. S59–S63, cuidando SBR_c frente a SBR total).
  - Mapas de CRB en el FOV.
  - Escalamientos L, N, SBR y L óptimo con λ_b fijo.
  - Monte Carlo MLE frente a LMS/mLMS dentro y fuera del TCP.
  - `mle_efficiency_center` con error < 2 %.
- **R3**
  - Dona vectorial Richards-Wolf (máscara de vórtice + circular, handedness correcta e incorrecta, lineal; λ=640, NA=1.4, n=1.518) en `beams_vectorial.py`. Depende de `B_donut_optics.md`.
  - `vectorial_zero_depth_*`.
  - MINFLUX iterativo (diseño propio `L_{k+1} ∝ σ_k`) frente a cámara: `iterative_sigma_nm` y `camera_sigma_nm`.
- **R3/R4 (R5 del paper, prioridad humana)**
  - Barrido de eps (0.1–15 %) × L, y L_opt(eps, SBR, N).
  - Desalineación con estimador honesto frente a ingenuo; validar contra Masullo (Nano Lett. 2021 y tesis publicada): <0.5–1 nm de degradación con r<15 nm, N=1000, SBR=20, L=100.
  - Pedestal no uniforme, que produce sesgo.
- **R4**
  - `scripts/fig_*.py` (≥6: R1–R5 + esquema), `make_all_figures.py`, `compute_paper_numbers.py` (único escritor de `data/paper_numbers.json`, con `fwhm_nm=300`) y `structure/figures.json`, `structure/claims.json`.
- **R5**
  - `paper/` (main.tex, secciones, bib), `paper/provenance.json`, `scripts/check_provenance.py` (debe imprimir "all checks pass"), `reproduce.sh`, README con la tabla de figuras.
  - En el paper: la discontinuidad S27 frente al límite y la nota fwhm 300 frente a 360 (tabla de Masullo).
- **Opcionales:** numLMS (Eq. S54–S58); estimador fasor y MLE de 4 parámetros (Stallinga); sesgo por dipolo; crosstalk temporal.

## LO QUE MÁS SE PODRÍA HACER
1. **Test de integración entre los tres módulos de R1** (fisher, photons/patterns, closed_forms) más una revisión de API unificada. Importa porque las piezas se escribieron en paralelo contra un contrato y así se detectan diferencias de convención. Costo: bajo, ~1 tarea chica.
2. **CRB de cámara + escalamientos + L óptimo con fondo fijo.** Da los números de R2 y los exponentes que pide la aceptación. Costo: medio, 1 worker.
3. **Monte Carlo MLE/LMS/mLMS con error < 2 %.** Da `mle_efficiency_center` y la figura R3. Costo: medio (cómputo), 1 worker.
4. **Dona vectorial Richards-Wolf** (propia, validada contra PyFocus solo como chequeo). Es la parte R1 de la aceptación y lo más caro técnicamente. Costo: alto, 1 worker durante 1–2 rondas.
5. **No idealidades: eps × L, L_opt, desalineación ingenua.** Es prioridad humana y aporte original. Costo: medio a alto.
6. **MINFLUX iterativo frente a cámara.** Hay que diseñar la regla de zoom, que es aporte propio. Costo: medio.
7. **Figuras, `paper_numbers`, procedencia, manuscrito y README.** Es el cierre de la aceptación. Costo: medio a alto, 1–2 rondas con writer.
8. **Extras:** numLMS, fasor y MLE de 4 parámetros, dipolo fijo, crosstalk. Aportan valor científico pero no los pide la aceptación. Costo: medio cada uno.
