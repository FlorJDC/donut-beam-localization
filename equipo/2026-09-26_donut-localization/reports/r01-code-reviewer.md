# Ronda 1 — code-reviewer

## Qué hice
- Leí el ledger: `intent.md`, `state.json` (ronda 0, sin claims), `inbox.jsonl` (4 indicaciones), `reports/r01-pi.md`, los reportes de los workers A/B/C y `work/r01-claims-to-verify.md`.
- Leí completo el código nuevo (sin `git diff`, porque todo es untracked): `src/donutloc/{__init__,beams,patterns,photons,fisher,estimators,montecarlo,closed_forms}.py`, los 7 tests nuevos y `scripts/verify_crb_closed_forms.py`. Revisé la privacidad de `docs/derivations/crb_tcp_center.md`.
- **Hash de aceptación:** `sha256(tests/test_acceptance.py) = 7d198853e87c7bd94bb8be07a855651146e1d849310281546cd8089dcf262303`. Es correcto y el archivo no cambió.
- Escribí y corrí dos scripts de integración propios, que no son código del proyecto:
  - `work/review/r01_integration.py`: `make_model` → `fisher` / `estimators` / `montecarlo` / `closed_forms` / `test_acceptance._crb_center`.
  - `work/review/r01_edge.py`: casos borde.

## Salidas reales
- `python -m unittest discover -s tests` (Python 3.8.6) → `Ran 84 tests in 1.651s — OK (skipped=1)`. El test omitido es el de aceptación: `data/paper_numbers.json` todavía no existe.
- `python scripts/verify_crb_closed_forms.py` → `all 89 exact checks pass (max rel. error point 2.13e-07, limit 2.23e-07)`, exit 0, 0.7 s. Imprime `crb_center_lg_L50_N100_nm = 1.605096` (límite r→0, fwhm=300).
- Integración (`r01_integration.py`):
  - **A. Límite r→0 por las tres rutas.** `fisher.crb_limit(make_model(tcp_centers(L), make_beam(fwhm)))` coincide con `test_acceptance._crb_center` con diferencia relativa **0.0 exacta** y con `cf.crb_tcp_center_limit` a ≤2.3e-9. Probé (L,fwhm) = (50,300), (100,300), (50,360) y (150,250).
    - L=50, fwhm=300: límite 1.6050956 nm, puntual 1.8024719 nm (= S27).
  - **B. `eps` con fondo.** Probé `eps` ∈ {0.002, 0.05}, los dos modelos de cero y SBR ∈ {None, 5}, pasando por `photons` (convención "beam"). `fisher.crb(0)` coincide con `cf.crb_tcp_center_eps(..., "beam")` a ≤5e-10, y límite = puntual a ≤1e-8, así que es continuo. **Las convenciones de SBR y `eps` coinciden entre `photons` y `closed_forms`.**
  - **C. LMS.** Con conteos esperados en r=(0.5,−0.3), `lms` general (con `make_model`) coincide con `lms_tcp` a 1e-8, con y sin SBR=10.
  - **D. MLE en el centro** (Monte Carlo con `make_model`, L=50, N=100, 5000 repeticiones, semilla 42):

    | caso | eficiencia σ/crb_limit | sesgo |
    |---|---|---|
    | SBR=10 | **0.996 ± 0.010** | compatible con 0 |
    | sin fondo | **0.837 ± 0.008** | — |

    Los dos casos tardan unos 0.9 s. Coinciden con lo que reportó B: 0.988 ± 0.005 y 0.838 ± 0.004 con 20000 repeticiones.
    - Error estándar de σ por bootstrap: 0.0151 frente a `sigma_err` = 0.0195. La fórmula σ/√(2n) resulta conservadora en este caso.

Conclusión: **el pegamento `p_fn` funciona sin fricción.** Los tres módulos escritos en paralelo se enchufan tal cual, con las mismas formas, convenciones y números.

## Hallazgos (defectos con escenario concreto)
Ninguno es bloqueante para la ronda 1. Los ordeno por relevancia para las rondas siguientes.

1. **Pico falso en el centro de los mapas de CRB (`fisher.crb` / `crb_map` con el `p_min=1e-12` por defecto).**
   - `crb(p, [r,0], 100)` da 1.605096 para r ≥ 5e-5 nm y **1.802472 (S27) para r ≤ 1e-5 nm**.
   - Concretamente, `crb_map(p, np.linspace(-2,2,5), [0.0], 100)` devuelve `[1.619, 1.609, **1.802**, 1.609, 1.619]`. Cualquier grilla simétrica que incluya el origen (por ejemplo `linspace(-50,50,101)`) mete un píxel 12 % más alto en el centro. Además, "el CRB en el centro leído del mapa" sale igual a S27, cuando el número de aceptación es el límite.
   - Está documentado como convención, pero es una trampa para las figuras de R2. El script de mapas debería evitar el origen o usar `crb_limit` en ese píxel, y tendría que tener un test.
2. **`estimators.mle(..., tol<=0)` se cuelga.** `_pattern_search` divide `s` hasta que llega a 0.0 y la condición `s < tol` nunca se cumple. Reproducido: `mle([10,10,10,1], p, 50., tol=0.0)` no termina y lo maté a los 20 s. Falta validar `tol > 0`.
3. **Riesgo de memoria en `mle` con `grid_step` fino.** El chunk de la grilla global es fijo (4096 filas) y no depende del tamaño G de la grilla: la matriz `C[chunk] @ logp.T` ocupa 4096×G×8 B.

   | search_radius | grid_step | G (puntos) | memoria |
   |---|---|---|---|
   | 50 | 2 (default) | 1961 | 0.06 GB |
   | 50 | 0.5 | 31417 | **1.0 GB** |
   | 50 | 0.2 | 196321 | **6.4 GB** |

   Con el default no hay problema; un barrido con grilla fina en R3 revienta la RAM. El chunk debería escalar como ~cte/G.
4. **`crb_tcp_center_limit` con exponente no entero 1 < c < 2.** El código devuelve el valor puntual para todo `power > 1`, aunque el docstring dice c ≥ 2. Matemáticamente es correcto (el término central ∝ r^(2c−2) → 0), pero la convergencia es tan lenta que la definición numérica del proyecto (r0 = 1e-3) no la alcanza:

   | c | `fisher.crb_limit` (numérico) | forma cerrada |
   |---|---|---|
   | 1.1 | 3.4038 | 3.4823 (−2.3 %) |
   | 1.5 | coinciden a 4e-6 | |
   | 2 | coinciden a 1e-9 | |

   `make_beam` y `crb_tcp_center_point` aceptan c no entero, así que falta declararlo o restringirlo. Es menor.
5. **Probabilidades NaN en campo lejano.** Con un haz gaussiano a r = 6000 nm todas las intensidades hacen underflow a 0: `make_model` devuelve `p = [nan]*4` y `fisher.crb` devuelve `inf` en silencio, porque det = NaN cae en la rama "singular". Con LG pasa recién a más de ~5 µm. Es menor: solo importa si algún mapa o FOV llega tan lejos.
6. **Validaciones que faltan en `mle`.** Estos casos no rompen nada con entradas sanas:
   - conteos negativos se aceptan en silencio (`mle([-5,10,10,0])` → `[0, 25]`);
   - con K distinto del de `p_fn` el error es críptico, de `matmul`.

## Cobertura de tests
- Todo el código nuevo tiene tests y pasan.
- **Falta un test de integración en la suite** entre `photons.make_model` → `fisher` → `closed_forms` → `_crb_center`. Hoy solo existe como columna del script de verificación y en mi script de revisión. Es el ítem 1 del backlog de R2 y es barato. Todos los tests de fisher, estimators y montecarlo usan un modelo inline propio (a propósito), así que ningún unittest ejercita `make_model` con `fisher`, `mle` ni `run_mc`.
- `montecarlo.rmse` usa una normalización por eje, sqrt(mean|Δr|²/2). El plan pedía "Eq. 4.2 de Masullo", que nadie verificó porque el texto no está en `docs/literature/`. Queda abierto.
- `test_montecarlo` cubre el hallazgo de superficiencia con 500 repeticiones y una banda amplia [0.75, 0.92]. Está bien como regresión, pero no sirve como número del paper (hace falta un error < 2 %).

## Contra el intent (lo que falta del entregable)
Todo esto está planificado para R2–R5, pero hoy no existe:
- `data/paper_numbers.json` y `scripts/compute_paper_numbers.py`, y por eso la aceptación da SKIP;
- `structure/figures.json` y `structure/claims.json`;
- `scripts/fig_*.py`, `make_all_figures.py`, `scripts/check_provenance.py` y `reproduce.sh`;
- la dona vectorial Richards-Wolf (`vectorial_zero_depth_*`);
- el MINFLUX iterativo (`iterative_sigma_nm`);
- el CRB de cámara (`camera_sigma_nm`);
- los escalamientos (`crb_exponent_L/N`);
- `paper/` y la tabla de figuras del README.

La ronda 1 cumplió su parte del plan: el núcleo del paquete y la teoría del centro.

## Proceso y privacidad
- En los archivos versionables (`src/`, `tests/`, `scripts/`, `docs/derivations/`, reportes) no encontré contenido del trabajo no publicado de la autora. La única cita numérica "de práctica" es eps = 0.002 de Balzarotti p. 6. No leí `docs/private/`.
- `inbox.jsonl` (versionado) nombra la ruta `docs/private/C_pminflux_practice.md` y dice "en donas experimentales el mínimo residual puede ser de varios % del pico". Es una afirmación genérica y no revela datos, pero conviene que la autora confirme que está bien que sea pública.
- La cuarta línea del inbox y el cambio en `OBJECTIVE.md` fijan `mle_efficiency_center` con SBR=10. Esa decisión la tomó el orquestador "en nombre de la autora" **después** de que el caso sin fondo cayó fuera de la banda de aceptación. El test de aceptación no cambió (el hash está bien) y el comportamiento sin fondo se reporta, así que no es trampa. Aun así, **la autora debería confirmarla explícitamente**, porque cambia qué mide un número de la definición de terminado.
- `state.json` sigue en ronda 0 y sin claims. Su actualización le corresponde al orquestador.

```claims
[{"status": "verified", "text": "Suite completa: python -m unittest discover -s tests -> Ran 84 tests, OK (skipped=1, aceptación por falta de data/paper_numbers.json), Python 3.8.6; hash de tests/test_acceptance.py = 7d198853...2303 intacto"},
 {"status": "verified", "text": "Integración make_model -> fisher: crb_limit(make_model(tcp_centers(L), make_beam(fwhm))) == test_acceptance._crb_center (rel. 0.0) y == closed_forms.crb_tcp_center_limit (<=2.3e-9) para (L,fwhm) en {(50,300),(100,300),(50,360),(150,250)}; L=50,fwhm=300: 1.6050956 nm (límite), 1.8024719 nm (puntual = S27)"},
 {"status": "verified", "text": "Convención de SBR/eps coincide entre photons (sbr sobre la señal total del haz) y closed_forms.crb_tcp_center_eps(sbr_ref='beam'): <=5e-10 para eps en {0.002,0.05}, modelos constant/gaussian, SBR en {inf,5}; con eps>0 límite = puntual (continuo)"},
 {"status": "verified", "text": "estimators.lms (general, Jacobiano numérico, con make_model) == lms_tcp a 1e-8 con y sin SBR=10; lms_tcp tiene el signo y la escala correctos (recupera r=(0.5,-0.3) de conteos esperados)"},
 {"status": "verified", "text": "MLE vía make_model en el centro, L=50, N=100, 5000 reps, semilla 42: SBR=10 eficiencia 0.996+-0.010; sin fondo 0.837+-0.008 (superficiente), consistente con los números de Worker B; ~0.9 s por corrida; sigma_err es conservador (bootstrap 0.015 vs 0.0195)"},
 {"status": "verified", "text": "scripts/verify_crb_closed_forms.py: exit 0, 89 checks exactos pasan (max rel. 2.1e-7 punto, 2.2e-7 límite), columna de paquete activa"},
 {"status": "refuted", "text": "fisher.crb/crb_map con p_min=1e-12 por defecto: crb_map(p, linspace(-2,2,5), [0.0], 100) da [1.619,1.609,1.802,1.609,1.619]; el píxel del origen salta al valor puntual S27 (+12%) para r<=1e-5 nm; trampa para los mapas de R2 (documentado pero sin test ni salvaguarda)"},
 {"status": "refuted", "text": "estimators.mle con tol<=0 se cuelga: _pattern_search divide s hasta 0.0 y s<tol nunca se cumple (mle([10,10,10,1], p, 50., tol=0.0) matado a los 20 s); falta validar tol>0"},
 {"status": "refuted", "text": "estimators.mle, memoria: el chunk fijo de 4096 no escala con el tamaño G de la grilla; search_radius=50, grid_step=0.2 -> G=196321 -> 6.4 GB por chunk (grid_step=0.5 -> 1.0 GB)"},
 {"status": "refuted", "text": "closed_forms.crb_tcp_center_limit con power no entero 1<c<2 devuelve el valor puntual (el docstring dice c>=2); con c=1.1 la definición numérica del proyecto (r0=1e-3) da 3.4038 frente a 3.4823 de la forma cerrada (-2.3%); menor: declarar o restringir"},
 {"status": "unclear", "text": "Falta en la suite un test de integración photons.make_model -> fisher/estimators/montecarlo -> closed_forms/_crb_center; hoy solo lo ejercitan el script de verificación y el script de revisión (ítem 1 del backlog de R2)"},
 {"status": "unclear", "text": "montecarlo.rmse usa la normalización por eje sqrt(mean|dr|^2/2); el plan pedía la Eq. 4.2 de Masullo, que nadie verificó (el texto no está en docs/literature)"},
 {"status": "unclear", "text": "Probabilidades NaN en campo lejano (haz gaussiano a r=6000 nm: p=[nan]*4) y fisher.crb devuelve inf en silencio; conteos negativos aceptados en silencio por mle; menores"},
 {"status": "unclear", "text": "El orquestador fijó 'mle_efficiency_center' con SBR=10 en nombre de la autora después de que el caso sin fondo cayó fuera de la banda; el test de aceptación no cambió, pero la autora debería confirmarlo explícitamente"},
 {"status": "unclear", "text": "Contra el intent aún faltan (previsto para R2-R5): data/paper_numbers.json, compute_paper_numbers.py, structure/*.json, fig_*.py, check_provenance.py, dona vectorial, MINFLUX iterativo, CRB de cámara, escalamientos, paper/ y README"},
 {"status": "verified", "text": "Privacidad: ningún archivo versionable nuevo (src, tests, scripts, docs/derivations, reportes) contiene material del trabajo no publicado de la autora; inbox.jsonl solo nombra la ruta privada y hace una afirmación genérica"}]
```
