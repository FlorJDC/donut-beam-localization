# Ronda 2 — code-reviewer

## Qué revisé
- **El diff real.** Mientras revisaba, el orquestador hizo commit del change set como `df7a97c` ("Round 2 (execution, verification pending)"), así que el árbol quedó limpio. Revisé `git diff HEAD~1 HEAD`:
  - los cambios en `fisher`, `estimators`, `montecarlo`, `closed_forms` y `photons`, con sus tests;
  - los archivos nuevos `camera.py`, `vectorial.py` y `experiments.py`;
  - `test_camera`, `test_vectorial`, `test_experiments` y `test_integration`.
- **Contexto que leí:** OBJECTIVE, CLAUDE/AGENTS, intent, state.json, inbox, r02-pi, el reporte de W3, `work/r02-claims-to-verify.md` y `tests/test_acceptance.py`.
- **Scripts propios** en `work/review/r02/`:
  - `r1_fixes.py`: escenarios de R1;
  - `vec_checks.py`;
  - `cam_checks.py`;
  - `exp_checks.py`;
  - `mis_se.py`: SE entre patrones;
  - `edge.py`.
- No leí `docs/private/`.

## Suite y hash
- `python -m unittest discover -s tests` → `Ran 142 tests in 17.744s — OK (skipped=1)`. El skip es la aceptación, porque falta `data/paper_numbers.json`.
- sha256 de `tests/test_acceptance.py` = `7d198853e87c7bd94bb8be07a855651146e1d849310281546cd8089dcf262303`. Intacto.

## Las 4 refutadas de R1 (re-corridas con su escenario original)
1. **Pico S27 en `crb_map`.** Corregido.
   - `crb_map(p, linspace(-2,2,5), [0], 100)` da `[1.619, 1.609, 1.60510, 1.609, 1.619]`. Con `zero_policy="point"` se reproduce `1.80247`.
   - El barrido cerca del origen (|x| de 1e-7 a 1e-3) es continuo: 1.6050956 ± 3e-8.
   - En un mapa 2D con el origen y el cero periférico (0,25) exactos, todo es finito.
   - `crb(p,0,N)` sigue dando S27 por defecto.
2. **`mle` con tol ≤ 0.** Corregido.
   - tol = 0, tol < 0 y grid_step = 0 dan `ValueError`.
   - tol = 1e-300 termina en 0.14 s.
3. **Memoria de `mle`.** Corregido para el bloque de verosimilitud.
   - Con M = 2000 y grid_step = 0.2 (G ≈ 1.96e5), el pico es de 135 MiB y tarda 2.4 s.
   - Residual menor: la evaluación del modelo sobre la grilla es O(G·K) y no está acotada por `mem_budget`. Con radius = 200 y step = 0.2 (G ≈ 3.1e6) el pico llega a 623 MiB. El docstring es correcto: solo promete acotar el bloque.
4. **`crb_tcp_center_limit` con 1 < c < 2.** El docstring quedó corregido y hay warning para 1 < c < 1.5. Diferencia relativa de la forma cerrada contra `fisher.crb_limit`:

   | c | 1.1 | 1.5 | 2 |
   |---|---|---|---|
   | diferencia relativa | −2.16e-2 | −7.1e-6 | 5.2e-9 |

Unclear de R1, ya resueltas:
- El test de integración existe y pasa (tres caminos a 1e-6).
- NaN en campo lejano: `crb` da NaN y no inf, también en `crb_map` con la política limit.
- Los conteos negativos o no finitos y un K distinto ahora dan `ValueError`.
- `sigma_err = σ/(2√n)`.
- Normalización de `rmse`: documentada.

## Compatibilidad hacia atrás
- `experiments` solo usa `crb` con la política "point" por defecto: en puntos sin cero (FOV sin origen), o en r = 0 únicamente cuando eps > 0 o hay fondo. En el caso perfecto usa `crb_limit`. Es correcto.
- Las firmas y los defaults de R1 se mantienen, con la excepción declarada de `crb_map`.
- **Defecto nuevo en la política `"limit"`** (menor, pero real):
  - `fisher.crb(p, [0,0], np.array([100,400]), zero_policy="limit")` → `ValueError: cannot broadcast a non-scalar to a scalar array`.
  - Con `"point"` la misma llamada devuelve `[1.802, 0.901]`.
  - El caso de uso natural es barrer N en un punto para `crb_exponent_N`. La salida es pasar r repetido o usar `crb_limit`.

## vectorial.py
- **Algebra.** Rederivé la reducción de Bessel (E0 = a e^{imφ'}(isφ̂'+θ̂)/√2, órdenes m±1 y m). Coincide con la nota B §2.2: P = −I_B, Q = I_A, R = I_C para m = 2, y los signos de Ex, Ey y Ez son los de la nota.
- **2D directa contra Bessel.** El error relativo por componente es ≤ 1.1e-15 en 5 casos: correcta, opuesta, lineal con pol_angle = 0.3, charge = 2 y z = 200. Eso valida también la superposición lineal (se suman campos, no intensidades).
- **Convergencia en n_theta.** Valores a 201 → 3201: opuesta 0.8452440 → 0.8452407; lineal 0.3716495 → 0.3716482; curvatura estable a 1e-7; D_pp 384.6663 → 384.6662. Con 801 alcanza de sobra.
- **Handedness.** La condición lσ = +1 da cero perfecto para charge = ±1; `charge=-1, s=+1` da 0.845. Con charge = 0 la profundidad es 1 (máximo en el centro). La profundidad lineal no depende de pol_angle.
- **Números de la nota B (n = 1.5):** 0.86896 y 0.38231.
- **Beam interpolado circular** (lineal en ρ²):
  - interp/exacto = 1 − 2.5e-5 cerca del cero, y se conserva I ∝ ρ² (I/(cρ²) = 0.99998);
  - `crb_limit`: diferencia relativa de 2.4e-5;
  - CRB fuera del centro (7,3): 3e-4 a 6e-4.
- **Beam interpolado lineal: defecto como trampa.** La interpolación bilineal con grid_step = 5 da gradientes constantes por celda. En L = 50, N = 100:
  - CRB en (7,3): 66.99 interpolado frente a 62.36 exacto (+7.4 %);
  - `crb_limit`: 58.77 frente a 60.77 (−3.3 %).

  No hay warning ni test. `compare_crb_vectorial_vs_lg` usa `mode="exact"` y no está afectado.
- **Rendimiento.** Con `mode="exact"` el modelo tarda 10.2 s para 2000 puntos. Es inviable para MLE o mapas en R3, que tienen que usar el modo interp.
- **Borde no relevante con NA = 1.4:** `_max_search` busca solo hasta ρ = 1000 nm. Con NA ≲ 0.2 el anillo queda fuera y fallaría en silencio.

## camera.py
- **Probabilidades por píxel.** Coinciden con `norm.cdf` independiente a 1e-16. El orden es row-major, con y exterior (verificado con un emisor en +x).
- **Fisher separable independiente.** Da 5.2045491 nm, igual que el paquete (9×9, a = 100, N = 400).
- **Límite ideal.** Con píxel de 5 nm y ventana 241×241 da 1.0001 × σ/√N.
- **Convenciones de SBR.** Por píxel 500 equivale a total 500/81: 4.96418 por los dos caminos.
- **Renormalización en la ventana.** Está declarada, pero da un CRB finito aunque el emisor esté fuera: 10913 nm en r = 2000. Más lejos da NaN (r = 1e5).
- **Casos degenerados.** `sbr=0` y `n_pix=1` dan inf. Es correcto.

## experiments.py
- **Honesto contra ingenuo.** El honesto sí usa los centros perturbados (`p_true`). Los 4 ceros, incluido el central, se mueven exactamente |δ|.
- **Números aleatorios comunes.** Con δ = 0 los dos estimadores son idénticos bit a bit. Las direcciones son las mismas para todos los δ.
- **Iterativo, corrida independiente** (semilla 7, n_rep = 2000):
  - fijo: σ = 0.508 ± 0.0066 (SE por bootstrap), consistente con el 0.497 ± 0.011 de W3;
  - SBR = 10: 0.631;
  - adaptivo: 0.475;
  - cámara: 3.162.
  - El reparto de fotones es correcto: 1003/4 → [250, 250, 250, 253].
- **Radio de búsqueda.** 0.75 L_k es razonable. En la iteración 0 (L = 150, N_k = 250) el error |e| tiene mediana 4.7 nm, cuantil 99 % de 15.5 nm y máximo de 25.7 nm.
- **Regla adaptiva: documentación incorrecta.** El docstring dice que κ = 6 hace que el patrón siguiente contenga la verdad con probabilidad ~99 %. Pero σ_k es el CRB central (3.47), mientras que el error real de la iteración 0 es 4.29, porque el emisor está hasta L0/4 fuera del centro. En consecuencia, el 2.65 % de los emisores queda fuera de L_1/2 = 12.5 nm, no el 1 %.
- **SE de la desalineación: defecto estadístico.**
  - `sigma_se = sigma_mean/(2√(R·P))` ignora la varianza entre patrones. Con δ = 10 en el centro y 20 patrones de 200 repeticiones:
    - honesto: SE de la fórmula 0.0144, SE real entre patrones 0.0380 (2.6×);
    - ingenuo: 0.0157 frente a 0.0561 (3.6×).
  - Las σ por patrón tienen una dispersión de 0.17 a 0.25 nm.
  - Por eso la significancia de "σ del ingenuo 1.85 → 1.98" está sobreestimada. En cambio, `bias_abs_se` sí usa la dispersión entre patrones.
- **SE del iterativo.** La fórmula gaussiana σ/(2√R) subestima el bootstrap en ~16 % (0.0057 frente a 0.0066). La causa es que los errores son una mezcla heterocedástica, con curtosis de 3.3. Para el requisito de "< 2 %" conviene el bootstrap.
- **`recenter=False`.** Esa comparación está dominada por el truncamiento del disco de búsqueda: el radio es 0.75 L_k alrededor del origen, y r0_spread = 37.5 nm es mayor que 18.75 nm cuando L = 25. El 8.06 nm que da no es un efecto físico del zoom.
- **Rendimiento para R3.** El iterativo con n_rep = 2000 tarda 2 s. La desalineación con 1 δ × 5 patrones × 625 repeticiones tarda 2.7 s. No hay problemas de escala salvo el beam vectorial en modo exact.

## Contra el intent (lo que falta para R3)
Todavía no existen:
- `data/paper_numbers.json`;
- `scripts/compute_paper_numbers.py`;
- `structure/figures.json` y `claims.json`;
- `scripts/check_provenance.py`, que debe envolver `agent-team/bin/check_provenance.py` (existe) e imprimir "all checks pass";
- los `fig_*.py`, `make_all_figures.py`, `reproduce.sh`, `paper/*.tex` y el README.

Las claves de la aceptación están bien definidas y el código para producirlas ya existe. Hay una trampa:
- `crb_exponent_L` depende del rango de L en el ajuste log-log, y la tolerancia es ±0.05:

  | rango de L (nm) | exponente |
  |---|---|
  | [10..50] | 1.008 |
  | [10, 100] | 1.027 |
  | [25..150] | 1.074 (fuera de la tolerancia) |

  Hay que fijar L ≪ fwhm.

Las demás claves:
- `camera_sigma_nm` tiene que declarar si es la ideal σ_PSF/√N (3.162) o la pixelada 9×9 (≈3.29 para N = 1000). Las dos quedan por encima de `iterative_sigma_nm`.
- `mle_efficiency_center` (SBR = 10) y `vectorial_zero_depth_*` están listos.

## Privacidad
- Ningún archivo nuevo versionado tiene contenido de `docs/private/`, que está en `.gitignore`.
- Los reportes solo mencionan que no lo leyeron.

```claims
[{"status":"verified","text":"Suite completa: python -m unittest discover -s tests -> Ran 142 tests in 17.744s, OK (skipped=1, aceptación sin paper_numbers); sha256 de tests/test_acceptance.py = 7d198853...2303 intacto (change set revisado como HEAD~1..HEAD = df7a97c)"},
 {"status":"verified","text":"Fix R1 crb_map: escenario original crb_map(p, linspace(-2,2,5), [0], 100) -> [1.619,1.609,1.60510,1.609,1.619] (sin salto); zero_policy='point' reproduce 1.80247; continuo para |x| 1e-7..1e-3 (1.6050956±3e-8); crb(p,0,N) sigue dando S27 por defecto"},
 {"status":"verified","text":"Fix R1 mle tol<=0: tol=0, tol<0, grid_step=0 -> ValueError; tol=1e-300 termina en 0.14 s; conteos negativos/no finitos y K distinto -> ValueError"},
 {"status":"verified","text":"Fix R1 memoria de mle: M=2000, grid_step=0.2 (G~1.96e5) -> pico de 135 MiB, 2.4 s; el bloque (chunk, G) queda acotado por mem_budget (residual menor: la evaluación del modelo O(G*K) no está acotada; radio 200, step 0.2 -> 623 MiB)"},
 {"status":"verified","text":"Fix R1 closed_forms 1<c<2: docstring corregido y warning para 1<c<1.5; diferencia relativa contra fisher.crb_limit: c=1.1 -2.16e-2, c=1.5 -7.1e-6, c=2 5.2e-9"},
 {"status":"verified","text":"tests/test_integration.py existe y pasa (patterns+beams+make_model -> crb_limit = copia de _crb_center = closed form a 1e-6; MLE SBR=10 eficiente; lms = lms_tcp); NaN en campo lejano -> crb NaN (también en crb_map con la política limit)"},
 {"status":"verified","text":"vectorial.py: reducción de Bessel rederivada (coincide con nota B §2.2); 2D directa contra Bessel <=1.1e-15 en correcta, opuesta, lineal pol_angle=0.3, charge=2 y z=200 (valida la superposición lineal de campos); n_theta 801 contra 3201: opuesta 0.8452409/0.8452407, lineal 0.3716482, curvatura 7.042221e-5, D_pp 384.666; ls=+1 da cero exacto para charge ±1; nota B n=1.5: 0.86896/0.38231"},
 {"status":"verified","text":"make_vectorial_beam circular (interp lineal en rho^2): interp/exacto = 1-2.5e-5 cerca del cero, conserva I∝rho^2; crb_limit diferencia relativa 2.4e-5 (L=50,100), CRB en (7,3) 3e-4..6e-4"},
 {"status":"verified","text":"camera.py: p por píxel = norm.cdf independiente a 1e-16, orden row-major (y exterior); Fisher separable independiente 5.2045491 nm (9x9, a=100, N=400) = paquete; SBR_c=500 por píxel equivale a total 500/81 (4.96418); límite ideal 1.0001 sigma/sqrt(N); renormalización en la ventana declarada"},
 {"status":"verified","text":"experiments.misalignment_study: el honesto usa los centros perturbados (p_true); los 4 ceros (incluido el central) se mueven exactamente |delta|; con delta=0 honesto e ingenuo son idénticos bit a bit (números aleatorios comunes)"},
 {"status":"verified","text":"experiments.iterative_minflux, corrida independiente (semilla 7, n_rep=2000): fijo sigma=0.508±0.0066 (bootstrap), consistente con 0.497±0.011 de W3; SBR=10 0.631; adaptivo 0.475; cámara 3.162; reparto de fotones correcto (1003/4 -> [250,250,250,253])"},
 {"status":"verified","text":"Privacidad: ningún archivo versionado nuevo contiene material de docs/private/ (que está en .gitignore)"},
 {"status":"refuted","text":"experiments.misalignment_study sigma_se = sigma_mean/(2 sqrt(R P)) ignora la varianza entre patrones: con delta=10 en el centro, 20x200, SE de la fórmula 0.0144/0.0157 contra SE real entre patrones 0.0380/0.0561 (honesto/ingenuo, 2.6x/3.6x); la significancia de 'sigma ingenuo 1.85->1.98' está sobreestimada"},
 {"status":"refuted","text":"fisher.crb con zero_policy='limit' falla cuando N se transmite más allá de r: crb(p,[0,0],np.array([100,400]),zero_policy='limit') -> ValueError 'cannot broadcast a non-scalar to a scalar array', mientras que con 'point' devuelve [1.802,0.901] (barrido de N en un punto, relevante para crb_exponent_N)"},
 {"status":"refuted","text":"Docstring de experiments.l_schedule (adaptivo): 'kappa=6 contiene la verdad con ~99 %' es falso porque sigma_k es el CRB central (3.47) y el error real de la iteración 0 es 4.29 (emisor hasta L0/4 fuera del centro); el 2.65 % de los emisores queda fuera de L_1/2=12.5 nm (L0=150, N_k=250)"},
 {"status":"refuted","text":"make_vectorial_beam polarization='linear' (bilineal, grid_step=5) no sirve para el CRB y no avisa: con L=50, N=100 el CRB en (7,3) da 66.99 interpolado contra 62.36 exacto (+7.4 %) y crb_limit 58.77 contra 60.77 (-3.3 %); sin test ni warning (compare_crb usa mode='exact' y no está afectado)"},
 {"status":"unclear","text":"SE del iterativo sigma/(2 sqrt R): subestima el bootstrap en ~16 % (0.0057 contra 0.0066, errores heterocedásticos, curtosis 3.3); para el requisito '<2 %' de los números del paper conviene usar el bootstrap"},
 {"status":"unclear","text":"iterative_minflux recenter=False: el 8.06 nm de la última iteración está dominado por el truncamiento del disco de búsqueda (radio 0.75 L_k alrededor del origen < r0_spread=37.5 nm para L=25), no por la física; hay que interpretarlo con cuidado en el paper"},
 {"status":"unclear","text":"Pendiente para R3: faltan data/paper_numbers.json, compute_paper_numbers.py, structure/figures.json y claims.json, scripts/check_provenance.py (envoltorio de agent-team/bin/check_provenance.py), fig_*.py, paper y README; crb_exponent_L depende del rango (1.008 en [10..50], 1.074 en [25..150], fuera de ±0.05): fijar L<<fwhm; declarar si camera_sigma_nm es la ideal (3.162) o la pixelada (~3.29); la dona vectorial en mode='exact' es inviable para MLE o mapas (10 s por 2000 puntos): usar interp"},
 {"status":"unclear","text":"SBR=10 para mle_efficiency_center sigue pendiente de confirmación explícita de la autora (no bloquea)"}]
```
