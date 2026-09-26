# Ronda 2 — PI (plan)

## Estado leído
- `state.json`: en R1 quedaron 24 afirmaciones verificadas, 4 refutadas y 5 unclear. Los checks están en verde: 84 tests OK y 1 skip (aceptación, porque falta `data/paper_numbers.json`). El hash de aceptación está intacto.
- `inbox.jsonl`: las 4 indicaciones siguen vigentes.
  - "Centro" = límite r→0; se reportan el valor puntual S27 y el límite.
  - Profundidad finita del cero y desalineación con estimador ingenuo: **prioridad** (va en esta ronda, Worker 3).
  - fwhm según Eq. S17, default 300 nm; la tabla de Masullo se reproduce con ~360 nm.
  - `mle_efficiency_center` con SBR=10; la superficiencia sin fondo es un resultado del paper.
  - **Privacidad:** nada de `docs/private/` va a archivos versionados, y ningún worker necesita leerlo.
- Refutadas vivas (las 4 se resuelven en la tarea 1):
  - pico S27 en `crb_map` en el origen;
  - `mle` con tol ≤ 0 se cuelga;
  - memoria de `mle` (chunk fijo);
  - `crb_tcp_center_limit` con 1 < c < 2.
- Unclear vivas:
  1. Falta el test de integración en la suite → tarea 1.
  2. Normalización de `rmse` / Eq. 4.2 de Masullo. **Decisión PI:** se mantiene la normalización por eje sqrt(mean|Δr|²/2). Es consistente con la definición del CRB del proyecto; la Eq. 4.2 no está en las notas, así que no se cita. Se documenta en la tarea 1.
  3. NaN en campo lejano y conteos negativos → tarea 1.
  4. SBR=10 fijado por el orquestador. Se mantiene: es una indicación del inbox y no bloquea. Sigue pendiente la confirmación explícita de la autora, que se pide en el resumen de la ronda.
  5. Faltan piezas del entregable → R2–R5 (ver backlog).
- Nota B: para la dona vectorial se usan las integrales I_A, I_B, I_C de §2.2 y la regla de handedness ℓσ. **No** se usa la E_y de la Eq. 15 impresa de PyFocus. Números de chequeo en §6.3 y §8.
- Nota A §4.4 y §9: la cámara de Balzarotti está en las Eq. S59–S63. Cuidado con las convenciones de SBR: SBR_c = K_pix · SBR_total.

## Decisiones declaradas (defaults)
- **Parámetros ópticos por defecto:**
  - λ = 640 nm, NA = 1.4, n = 1.518 (OBJECTIVE).
  - Llenado de pupila F = w0/h con default **5/3**: es la referencia de la nota B, que da w0 = 5 mm y h = 3 mm.
  - Los tests que reproducen los números de la nota B usan **n = 1.5**, como la nota.
- **Todo modelo es un `p_fn(r)`**, de forma (...,2) → (...,K), igual que en R1. Esto incluye la cámara: K = n_pix², con un píxel por "exposición". Así `fisher.crb` sirve para la cámara sin código nuevo.
- **Un "beam" es un `f(x, y)` vectorizado con pico 1**, compatible con `photons.make_model`. La dona vectorial se entrega como beam interpolado radialmente (casos con simetría de revolución) o como beam por grilla 2D (polarización lineal).
- **Referencia de cámara en los drivers del Worker 3:** la ideal σ_PSF/√N con σ_PSF = 100 nm (Balzarotti p. 32), calculada inline. El W3 **no importa** `camera.py`, que se escribe en paralelo. En R3, `compute_paper_numbers` elige qué número de cámara se publica.
- **Compatibilidad hacia atrás:** el Worker 1 puede agregar kwargs a `fisher`, `estimators`, `montecarlo`, `closed_forms` y `photons`, pero **no** cambiar ni quitar firmas ni defaults existentes. Hay una excepción: `crb_map` pasa a usar por defecto la política "limit" en los ceros (ver tarea 1). `crb(p_fn, 0, N)` sigue dando S27 por defecto.
- **Entorno:** Python 3.8, solo numpy y scipy, unittest, semilla 42, unidades nm. Encabezado de `sys.path` de 3 líneas en cada test, como en R1. Prohibido tocar `tests/test_acceptance.py` y copiar código de terceros (incluido PyFocus).
- **Tests rápidos:** cada archivo de test nuevo debe correr en menos de ~20 s. Los números grandes van en el reporte, no en los tests.

## Tareas (propiedad de archivos exclusiva)

1. **Worker 1: corregir los defectos de R1, agregar el test de integración y el CRB de cámara.**
   **Archivos propios:**
   - `src/donutloc/fisher.py`, `estimators.py`, `montecarlo.py`, `closed_forms.py`, `photons.py`;
   - `src/donutloc/camera.py` (nuevo);
   - `tests/test_fisher.py`, `test_estimators.py`, `test_montecarlo.py`, `test_closed_forms.py`, `test_photons.py`;
   - `tests/test_integration.py` (nuevo) y `tests/test_camera.py` (nuevo).

   No toca `beams.py`, `patterns.py` ni `__init__.py`. Los Workers 2 y 3 usan sus APIs públicas en paralelo, así que **no se rompen firmas ni defaults**.

   **(a) Pico S27 en los mapas.** Agregar a `fisher.crb` y a `crb_map` el kwarg `zero_policy in {"point","limit"}`.
   - Con `"limit"`, en los puntos donde algún p_i ≤ p_min se reemplaza el valor por `crb_limit` alrededor de ese punto (r0 = 1e-3, n_dir = 12).
   - Defaults: `crb` → `"point"`, que conserva S27 y los tests actuales; `crb_map` → `"limit"`.
   - Test: `crb_map(p, linspace(-2,2,5), [0], 100)` sin salto en el centro (≈1.6051 en L=50, fwhm=300). Con `"point"` se reproduce 1.8025.

   **(b) `mle` con tol ≤ 0.** Tiene que dar `ValueError` (también con `grid_step` ≤ 0 y `search_radius` ≤ 0). Además se valida:
   - que no haya conteos negativos ni no finitos;
   - que K coincida con `p_fn`, con un mensaje claro.

   **(c) Memoria de `mle`.** El chunk se adapta al tamaño de la grilla, `chunk = max(1, floor(mem_budget/(8·G)))`, con un kwarg `mem_budget` de ~64 MB por defecto. Se mantiene el kwarg `chunk` como tope. Test: `grid_step` fino (G ~ 2e5) con pocas filas, sin reventar la memoria y con el mismo resultado que el grueso refinado.

   **(d) `closed_forms.crb_tcp_center_limit` con 1 < c < 2.**
   - Docstring: "límite = puntual para todo c > 1 (término central ∝ r^(2c−2)). La convergencia con r0 finito es lenta para c → 1⁺".
   - Emitir `warnings.warn` si 1 < c < 1.5.
   - Test con c = 1.5 y c = 2 contra `fisher.crb_limit` (4e-6 y 1e-9).

   **(e) Pendientes menores.**
   - `montecarlo.sigma_err = sigma/(2·sqrt(n_rep))`. Es el valor que corrigió el verificador (V13); hay que documentar la normalización por eje de `rmse`.
   - NaN en campo lejano: si ΣI = 0 sin fondo, `photons.probabilities` devuelve NaN documentado, y `fisher.crb` devuelve **NaN**, no inf, en esos puntos.
   - Test con un haz gaussiano a r = 6000 nm.

   **(f) `tests/test_integration.py`.** Usa el camino real: `patterns.tcp_centers` + `beams.make_beam` + `photons.make_model`, hacia `fisher.crb_limit`, `closed_forms.crb_tcp_center_limit` y una copia local de `_crb_center` del test de aceptación (copiar las fórmulas, no importar). Se prueban (L, fwhm) ∈ {(50,300), (100,300), (50,360)} con tolerancia 1e-6. Además:
   - `run_mc(mle)` con SBR=10, L=50, N=100 y n_rep = 1000: eficiencia en [0.9, 1.1];
   - `lms` coincide con `lms_tcp` usando `make_model`.

   **(g) `src/donutloc/camera.py` (nuevo), cámara ideal según Balzarotti Eq. S59–S63 (nota A §4.4).**
   - `make_camera_model(sigma_psf=100.0, pixel=100.0, n_pix=9, sbr=None, sbr_convention="total")`:
     - Devuelve un `p_fn` de forma (...,2) → (..., n_pix²) con p por píxel integrado con erf y la PSF gaussiana centrada en r. La grilla está centrada en (0,0).
     - Fondo uniforme por píxel. `sbr_convention="total"` es señal total / fondo total (Masullo); `"per_pixel"` es SBR_c de Balzarotti (señal total / fondo por píxel), con SBR_c = K·SBR_total.
     - Hay que documentar la señal fuera del campo: se renormaliza dentro de la ventana, y eso tiene que quedar declarado.
   - `crb_camera_ideal(sigma_psf, N)`, que da σ_PSF/√N.
   - `crb_camera(sigma_psf, N, pixel, n_pix, sbr=None, r=(0,0), sbr_convention="total")`, vía `fisher.crb` sobre el `p_fn`.
   - Tests:
     - píxel chico y ventana grande sin fondo → σ_PSF/√N a 1 %;
     - la pixelación (a = σ_PSF) empeora el CRB, y el fondo lo empeora de forma monótona;
     - check de Balzarotti p. 1: σ_PSF = 100 nm y N = 400 dan ≈ 5 nm;
     - check de nota A §9: SBR_c = 500 con 81 píxeles ≡ SBR total ≈ 6.2 (el mismo CRB por las dos convenciones).
   - Número para el reporte (no para el test): el CRB de cámara con σ_PSF = 100, a = 100, 9×9, SBR_c = 500 y N = 600 frente al MINFLUX de L = 50 con S31, para comentar el "22×".

   **Afirmaciones a reportar** (en el bloque `claims`):
   - cada defecto corregido, con el escenario de R1 y el resultado nuevo;
   - la salida de la suite completa;
   - los números de cámara, con sus convenciones.

2. **Worker 2: dona vectorial Richards-Wolf para R1 del OBJECTIVE.**
   **Archivos propios:** `src/donutloc/vectorial.py` (nuevo) y `tests/test_vectorial.py` (nuevo).
   - Implementación propia a partir de la nota B §2.2 y §7.2: integrales 1D I_A, I_B, I_C, con g(θ) = a(θ)√cosθ sinθ e^{ikz cosθ}, u = kρ sinθ, k = 2πn/λ y a(θ) = exp(−(f sinθ/w0)²).
     - Se parametriza con F = w0/h, sin h absoluto: f sinθ/h = sinθ/sinα, con sinα = NA/n.
     - Trapecio con ≥ 801 puntos en θ.
   - **Handedness** según la regla ℓσ de §2.2:
     - σ = +1 (misma mano), cero perfecto: E_z ∝ J2 y las transversales ∝ J1, J3.
     - σ = −1 (mano opuesta): E_z ∝ J0 con máximo central y las transversales ∝ J1. **Derivar y documentar** los órdenes de Bessel y los prefactores para σ = −1 en el docstring, y chequearlos contra una integral 2D directa de la Eq. 1 de Caprile (N_θ = N_φ = 200) en un test chico.
   - **Polarización lineal:** sin simetría de revolución. Se calcula con una integral 2D directa en (θ, φ'), o con la identidad de Bessel por componentes armónicas si se deriva, y se chequea contra la integral 2D.
   - **API:**
     - `focal_field(rho, phi, z=0, wavelength=640., NA=1.4, n=1.518, filling=5/3, polarization="circular", handedness=+1, charge=1, n_theta=801)`. `polarization` ∈ {"circular", "linear"}, con ángulo opcional `pol_angle`. Devuelve (Ex, Ey, Ez) sin normalizar.
     - `intensity(x, y, z=0, **opt)`: intensidad total |E|², vectorizada.
     - `radial_profile(rho, **opt)`: solo para casos con simetría de revolución (circular ±). Da error claro para "linear".
     - `zero_depth(**opt)`, que devuelve I(0)/I_max. I_max se busca en una grilla 2D para la lineal.
     - `peak_to_peak_diameter(**opt)` y `zero_curvature(**opt)`: el coeficiente c de I/I_max ≈ c ρ² cerca del cero (en nm⁻²); para lineal, el promedio azimutal.
     - `make_vectorial_beam(rho_max=1500., d_rho=1., eps=0.0, **opt)`: devuelve un beam `f(x, y)` con pico 1 (interpolación radial con `np.interp`, cero fuera de rho_max), compatible con `photons.make_model`.
       - Para "linear", la interpolación es 2D sobre una grilla cartesiana (paso ≤ 5 nm), con `scipy.interpolate.RegularGridInterpolator`.
       - `eps` opcional como en LG: `+ eps·(gaussiano normalizado)`, o bien omitido y documentado.
     - `lg_equivalent_fwhm(beam_opts, match="curvature"|"diameter")`: el fwhm de la LG del proyecto que iguala la curvatura en el cero (4e ln2/fwhm² = c) o el diámetro pico a pico (= 2·fwhm/(2√ln2)).
     - `compare_crb_vectorial_vs_lg(L_list, N=100, **opt)`: usa **solo** la API pública (`patterns.tcp_centers`, `photons.make_model`, `fisher.crb_limit`). Devuelve un dict con el CRB del límite central para la dona vectorial y para la LG equivalente por curvatura y por diámetro.
   - **Tests** (con n = 1.5 donde se reproduzca la nota B):
     - mano correcta: zero_depth < 1e-10;
     - mano opuesta: 0.869 ± 0.01 con F = 5/3;
     - lineal: 0.382 ± 0.01 con F = 5/3;
     - D_pp uniforme (F grande, p. ej. 100): ≈ 380 nm ± 3 %; con F = 2, ≈ 392 nm ± 3 %;
     - curvatura ≈ 7.0e-5 nm⁻² (F = 5/3) ± 5 %;
     - simetría de revolución de la circular: la intensidad 2D a φ = 0, 1 y 2 rad coincide a 1e-10;
     - la versión con I_A, I_B, I_C coincide con la integral 2D directa a 1e-6 en unos pocos puntos;
     - conservación de energía: el flujo transversal ∫ I 2πρ dρ en z = 0 es igual para σ = ±1 dentro de 1e-3, y se documenta si no lo es por la componente z;
     - `make_vectorial_beam` da pico 1 y se enchufa a `make_model` con Σp = 1.
   - **Afirmaciones a reportar**, con los defaults n = 1.518 y F = 5/3 y también con los parámetros de la nota B:
     - zero_depth correcta, opuesta y lineal;
     - D_pp, curvatura y fwhm LG equivalente (por curvatura y por diámetro);
     - la tabla del CRB del límite central, vectorial frente a LG, con L ∈ {50, 100, 150} y N = 100. Esta es la respuesta a R1.
     - Tiempos de cálculo. Los candidatos a `paper_numbers` son `vectorial_zero_depth_correct` y `vectorial_zero_depth_wrong_handedness`.

3. **Worker 3: drivers de simulación para R4 y R5.**
   **Archivos propios:** `src/donutloc/experiments.py` (nuevo) y `tests/test_experiments.py` (nuevo).
   - Solo usa la API pública existente: `patterns.tcp_centers` y `perturb_centers`, `beams.make_beam`, `photons.make_model` y `sample_counts`, `fisher.crb`, `crb_limit` y `crb_map` (pasando kwargs explícitos), `estimators.mle` y `lms`, `montecarlo.run_mc`.
   - No importa `camera` ni `vectorial`, que se escriben en paralelo. Todo beam entra como argumento (default LG fwhm = 300), para que en R3 se pueda pasar la dona vectorial.
   - Son funciones puras que devuelven dicts de arrays y registran los parámetros y la semilla (42).
   - **(a) `iterative_minflux(N_total, L_schedule=None, n_iter=4, L0=150., L_min=None, photon_split="equal", beam=None, sbr=None, r0_spread=None, n_rep=500, seed=42, recenter=True, rule="fixed")`.** MINFLUX iterativo con recentrado.
     - La posición verdadera se sortea uniforme en un disco de radio `r0_spread` (default L0/4) alrededor del origen. La primera TCP se centra en el origen.
     - En la iteración k, la TCP de diámetro L_k se centra en la estimación previa, se muestrean N_k fotones y se estima con MLE en un disco de radio ~L_k/2 + margen alrededor del centro de la TCP.
     - `rule="fixed"`: L_k geométrico de L0 a L_min. `rule="adaptive"`: L_{k+1} = max(L_min, κ·σ_k) con σ_k del CRB. Se documenta κ.
     - El fondo es SBR fija por iteración (`sbr`) y queda declarado. Opcional: fondo por exposición fijo, para mostrar que la SBR cae con L chico.
     - Devuelve:
       - σ final (por eje, respecto de la verdad), sesgo y rmse;
       - σ por iteración;
       - `camera_sigma = sigma_psf/sqrt(N_total)` (con `sigma_psf` = 100 como kwarg);
       - la curva σ_final(N_total) para una lista de N_total (función `iterative_vs_photons(N_list, ...)`);
       - la relación con σ_PSF/√N a igual presupuesto de fotones.
   - **(b) `eps_L_sweep(eps_list, L_list, N=100, fwhm=300., sbr=None, zero_model="gaussian", fov_radius=None)`.** Devuelve la matriz del CRB en el centro. Con eps > 0 el CRB es continuo, así que se usa `crb` en r = 0; con eps = 0 se usa `crb_limit`, y se documenta. Opcionalmente también el promedio en un disco de radio `fov_radius` (default L/4). Se agrega `optimal_L(eps, ...)`, que da el argmin del CRB sobre L, con y sin fondo. El test muestra que con eps > 0 el L óptimo es finito y crece con eps.
   - **(c) `misalignment_study(displacement_list, L=100., N=500, fwhm=300., sbr=10, n_patterns=20, n_rep=200, seed=42, estimator="mle")`.** Para cada desplazamiento δ:
     - se sortean `n_patterns` TCP perturbadas con `perturb_centers` (dirección aleatoria);
     - se simulan conteos con los centros verdaderos perturbados;
     - se estima con el modelo **honesto** (conoce los centros perturbados) y con el **ingenuo** (asume la TCP ideal).
     - Devuelve, para cada estimador, el sesgo medio |b|, σ y rmse, promediados sobre patrones, junto con el CRB honesto.
     - Evaluar en el centro y en uno o dos puntos fuera de él (p. ej. (L/4, 0)).
   - **Tests** con configuraciones chicas (n_rep ≤ 200, menos de ~20 s en total):
     - (a) σ del iterativo es menor que `camera_sigma` para N_total = 1000, L0 = 150 → L_min = 25, SBR = None o 10, y σ decrece con N_total.
     - (b) con eps = 0 el CRB crece con L (L/√N); con eps = 0.05 existe un mínimo interior.
     - (c) con δ = 0, honesto = ingenuo exactamente; con δ > 0 el sesgo ingenuo es mayor que el honesto (el honesto tiene sesgo compatible con 0 dentro de 3 SE).
     - Reproducibilidad con la semilla 42.
   - **Afirmaciones a reportar**, con errores estándar calculados como sigma/(2√R):
     - el número candidato `iterative_sigma_nm` frente a `camera_sigma_nm` en una configuración de referencia declarada (N_total = 1000, σ_PSF = 100, esquema elegido);
     - la curva frente a N con su pendiente log-log;
     - la tabla de L óptimo frente a eps (0.002, 0.01, 0.05, 0.15), con y sin SBR = 10;
     - el sesgo y la pérdida de precisión del ingenuo frente al honesto para δ ∈ {0, 2, 5, 10} nm;
     - los tiempos.

## Backlog (nada se descarta)
- **R3:**
  - `scripts/compute_paper_numbers.py` (el único que escribe `data/paper_numbers.json`). Lleva `fwhm_nm` y `crb_center_lg_L50_N100_nm` (límite); `crb_exponent_L` y `crb_exponent_N` (ajuste log-log; L ≪ fwhm para el de L); `mle_efficiency_center` (SBR=10, con error < 2 %); `vectorial_zero_depth_*`; `iterative_sigma_nm` y `camera_sigma_nm`.
  - `structure/figures.json` y `structure/claims.json`.
  - `scripts/check_provenance.py` como envoltorio de `agent-team/bin/check_provenance.py`, que debe imprimir "all checks pass".
  - Mapas de CRB (con `zero_policy="limit"`) y MC de MLE frente a LMS/mLMS en el FOV (R3 del OBJECTIVE).
- **R3–R4:** al menos 6 figuras, `scripts/fig_*.py` y `make_all_figures.py`:
  - esquema;
  - R1: perfil vectorial frente a LG y profundidad del cero por handedness;
  - R2: mapa de CRB y escalamientos, junto con la cámara;
  - R3: estimadores;
  - R4: iterativo frente a cámara;
  - R5: eps × L y desalineación.
- **R4–R5:** `paper/main.tex` + `sections/` + `references.bib` + `provenance.json` con `\src{}`, `README.md` y `reproduce.sh`.
- **Para el paper:**
  - sesgo del MLE sin fondo ≈ −0.34 ± 0.01 nm;
  - V10: 4.9 nm es una escala, no el argmin;
  - la discontinuidad S27/límite;
  - la superficiencia;
  - la advertencia sobre la E_y de la Eq. 15 de PyFocus.
- Pendiente con la autora: confirmar SBR=10 para `mle_efficiency_center` (no bloquea).

## LO QUE MÁS SE PODRÍA HACER
1. Dona top-hat / 3D (K0–K2 de la nota B §3) con TCP 3D y CRB en z. Extiende R1 a 3D. Costo aproximado: 1 ronda de un worker.
2. Aberraciones de Zernike (nota B §4.2) con integral 2D: profundidad del cero frente a la amplitud de la aberración, más realista que eps constante. Costo aproximado: 1 ronda de un worker.
3. Contraste residual de polarización del 1–5 % (Lopez), mapeado a eps con la dona vectorial elíptica. Conecta eps con un parámetro medible. Costo aproximado: media ronda.
4. Cámara realista con ruido de lectura y exceso EM (factor 2–3 de Masullo) en lugar de la ideal. Da una comparación más justa. Costo aproximado: media ronda.
5. Regla óptima de reducción de L (L_{k+1} ∝ σ_k) derivada analíticamente, con la ley de escala del presupuesto total. Es un aporte original porque no está en las fuentes. Costo aproximado: 1 ronda de derive con verificación.
6. Estimador bayesiano o MLE con prior en el iterativo, frente a recentrado puro, siguiendo el comentario de Balzarotti p. 6. Costo aproximado: media ronda.
7. Acelerar con numba el MC de mapas en el FOV. Hace falta para cumplir el error < 2 % en todas las celdas. Costo aproximado: media ronda.
