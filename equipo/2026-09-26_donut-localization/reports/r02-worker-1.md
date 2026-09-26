# Ronda 2 — Worker 1 (fixes R1 + integración + cámara)

Lo que ejecuté es la tarea 1 de `reports/r02-pi.md`. Solo edité los archivos que tenía asignados:
- `src/donutloc/{fisher,estimators,montecarlo,closed_forms,photons}.py`;
- `src/donutloc/camera.py` (nuevo);
- `tests/test_{fisher,estimators,montecarlo,closed_forms,photons}.py`;
- `tests/test_integration.py` y `tests/test_camera.py` (nuevos).

No toqué `beams.py`, `patterns.py`, `__init__.py` ni `tests/test_acceptance.py`. El hash sigue siendo `7d198853…2303`. No leí `docs/private/`. No hice commit.

## Cambios de API (exactos)
Todos son aditivos. La única excepción es el default de `crb_map`, como estaba planificado.

**`fisher.crb`**
- Firma nueva: `crb(p_fn, r, N, zero_policy="point", limit_r0=1e-3, limit_n_dir=12, **kw)`.
- Con `"limit"`, en los puntos donde algún `p_i <= p_min` (el `p_min` de `kw`, default 1e-12) el valor se reemplaza por `crb_limit(p_fn, N_punto, r_center=punto, r0=limit_r0, n_dir=limit_n_dir)`. Ese cálculo usa h = r0·1e-2 y p_min = 0.
- Si `zero_policy` no es válido, da `ValueError`.

**`fisher.crb_map`**
- Firma nueva: `crb_map(p_fn, xs, ys, N, zero_policy="limit", **kw)`. **El default cambió** de "point" a "limit".

**`fisher._cov_from_F`** (afecta a `crb` y `crb_axes`)
- Si det = NaN, devuelve NaN. `inf` queda solo para una Fisher singular finita.

**`estimators.mle`**
- Firma nueva: `mle(..., tol=1e-4, chunk=4096, mem_budget=64*2**20)`.
- El chunk efectivo es `min(chunk, max(1, floor(mem_budget/(8G))))`.
- Da `ValueError` en estos casos:
  - `search_radius`, `tol` o `grid_step` que no sean finitos y > 0 (incluye NaN e inf);
  - `chunk < 1` o `mem_budget <= 0`;
  - conteos negativos o no finitos;
  - `counts.ndim` distinto de 1 o 2;
  - `center` con forma distinta de (2,);
  - K de los conteos distinto de `len(p_fn(center))`, con el mensaje "counts have K = 3 exposures but p_fn returns K = 4 probabilities".

**`montecarlo.run_mc`**
- `sigma_err` pasa a ser `sigma/(2·sqrt(n_valid))`; antes era `sigma/sqrt(2 n_valid)`.
- El docstring documenta la normalización por eje de `rmse`, que es la decisión del PI. No se cita la Eq. 4.2 de Masullo.

**`closed_forms.crb_tcp_center_limit`**
- La firma no cambia.
- El docstring ahora dice "límite = puntual para todo c > 1 (término central ∝ r^(2c−2)); convergencia lenta con r0 finito para c → 1⁺".
- Emite `UserWarning` si algún valor de `power` cumple 1 < c < 1.5.

**`photons.probabilities`**
- Solo cambia el docstring: documenta que da NaN en campo lejano (también con `sbr`) y `1/K` con `bg_per_exposure > 0`.

**Nuevo `donutloc.camera`**
- `pixel_edges(pixel=100., n_pix=9)`.
- `make_camera_model(sigma_psf=100., pixel=100., n_pix=9, sbr=None, sbr_convention="total")`:
  - p_fn de forma (...,2) → (..., n_pix²), con los píxeles en orden row-major (iy afuera, ix adentro);
  - erf estable (usa erfc en las colas);
  - renormaliza dentro de la ventana (declarado);
  - da NaN lejos de la ventana;
  - `"per_pixel"` corresponde a SBR_c = K·SBR_total.
- `crb_camera_ideal(sigma_psf, N)`, que da σ/√N.
- `crb_camera(sigma_psf, N, pixel=100., n_pix=9, sbr=None, r=(0,0), sbr_convention="total")`, que usa `fisher.crb`.

## Defectos de R1 → estado nuevo
1. **Pico S27 en los mapas.** `crb_map(make_model(tcp_centers(50), LG300), linspace(-2,2,5), [0], 100)` ahora da `[1.61906, 1.60858, 1.60510, 1.60858, 1.61906]`. Con `zero_policy="point"` reproduce el 1.80247 del centro. `crb(p, 0, 100)` sigue dando 1.8024719, que es S27.
2. **`mle` con tol ≤ 0.** Ahora da `ValueError` para `tol` ∈ {0, −1e-3, NaN}, `grid_step` ∈ {0, −1}, `search_radius` ∈ {−5, inf}, `chunk=0` y `mem_budget=0`. Ya no se cuelga.
3. **Memoria de `mle`.** Con `search_radius=50`, `grid_step=0.2` (G = 196321) y 200 filas, el pico de tracemalloc es de 134.8 MiB. Forzando el chunk viejo de 4096 (mem_budget=1e12) sube a 309.6 MiB. El resultado difiere del cálculo grueso refinado en 8.5e-5 nm como máximo. El test exige < 220 MiB y coincidencia a 2e-3 nm.
4. **`crb_tcp_center_limit` con 1 < c < 2.** Comparado con `fisher.crb_limit` (L=100, N=100, fwhm=300):
   - c = 1.5: −4.39e-6;
   - c = 2: 1.20e-9;
   - c = 1.1: −2.25 %, con warning.

## Pendientes menores
- **`sigma_err`.** Hice 60 semillas con n = 400 (LMS). La dispersión de σ entre semillas coincide con el `sigma_err` nuevo dentro de ±30 %, y es < 0.85 × la fórmula vieja. Está como test.
- **Campo lejano.** Con un haz gaussiano y r = (6000, 0):
  - p = NaN;
  - `crb` da NaN con las dos políticas;
  - `crb_axes` da NaN;
  - con `bg_per_exposure` el resultado es uniforme, 0.25.

## Integración (`tests/test_integration.py`)
- Comparé `fisher.crb_limit(make_model(tcp_centers(L), make_beam(fwhm)))` con la copia local de `_crb_center` y con `crb_tcp_center_limit`:

  | (L, fwhm) | valor | vs `_crb_center` | vs forma cerrada |
  |---|---|---|---|
  | (50, 300) | 1.6050956 | diferencia 0.0 | 2.3e-9 |
  | (100, 300) | 3.3634111 | diferencia 0.0 | 5.3e-10 |
  | (50, 360) | 1.5976973 | diferencia 0.0 | 2.3e-9 |

- `crb_map` en el origen da el mismo número.
- `run_mc(mle)` con make_model, SBR=10, L=50, N=100, 1000 reps y semilla 42 da σ = 1.8973 ± 0.0300. La eficiencia σ/crb_limit es 0.968; el CRB es S31 = 1.96006 a 1e-6.
- `lms(make_model)` coincide con `lms_tcp` a 1e-8 con conteos esperados y a 1e-7 con 50 muestras multinomiales, con y sin SBR=10.

## Cámara (`tests/test_camera.py`, más los números del reporte)
**Checks a N = 400, sin fondo:**
- Píxel de 10 nm en una ventana de 121×121: CRB/(σ/√N) = 1.00042. Con píxel de 5 nm y ventana de 241×241: 1.00010.
- Balzarotti p. 1: `crb_camera_ideal(100, 400)` = 5.000. Con píxel de 100 nm y 9×9 sin fondo, 5.2045 nm (el efecto de la pixelación a = σ).

**Fondo** (a = 100, 9×9, N = 400, convención total). El CRB crece de forma monótona:

| SBR | CRB (nm) |
|---|---|
| inf | 5.2045 |
| 100 | 5.2898 |
| 10 | 5.7867 |
| 5 | 6.2514 |
| 1 | 9.4173 |

**Convenciones (nota A §9).** SBR_c = 500 con K = 81 equivale a SBR_total = 500/81 = 6.173. Las dos convenciones dan el mismo CRB (diferencia relativa 0.0). Las probabilidades coinciden a 1e-13.

**Número para el paper ("22×")**
- Cámara con σ_PSF = 100, a = 100, 9×9, SBR_c = 500 y N = 600: **4.964 nm**. La ideal σ/√N a N = 600 da 4.082.
- Para 5 nm la cámara necesita 591 fotones. Esto reproduce el "~600 fotones para 5 nm" de Balzarotti (Fig. 3).
- MINFLUX con L = 50, fwhm = 300 y la Eq. S31 en el centro (valor puntual):

  | SBR | fotones para 5 nm | cámara/MINFLUX |
  |---|---|---|
  | inf | 13.0 | 45.5× |
  | 50 | 13.5 | 44.0× |
  | 20 | 14.2 | 41.8× |
  | 10 | 15.4 | 38.5× |
  | 5 | 17.9 | 33.0× |

  - Con fwhm = 360 los valores son casi iguales: 12.8 fotones y 46×.
  - Con el límite sin fondo hacen falta 10.3 fotones.
  - σ de MINFLUX a N = 600: 0.736 nm (sin fondo) y 0.800 nm (SBR = 10).
- Conclusión: el "22×" del paper (600/27) usa ~27 fotones experimentales. La teoría S31 da entre 33× y 46× según la SBR, así que el 22× es conservador respecto del CRB ideal.

## Suite
`python -m unittest discover -s tests` (Python 3.8, desde la raíz) → **`Ran 128 tests in 4.540s — OK (skipped=1)`**.
- El test omitido es el de aceptación, porque falta `data/paper_numbers.json`.
- La corrida incluyó `test_experiments.py`, que ya existía (es del Worker 3). `test_vectorial.py` todavía no estaba.

Mis archivos solos (test_fisher, test_estimators, test_montecarlo, test_closed_forms, test_photons, test_integration, test_camera) → `Ran 94 tests in 2.801s — OK`.

`scripts/verify_crb_closed_forms.py` sigue pasando los 89 checks.

## Lo que queda abierto
- `crb_axes` no tiene `zero_policy`: en un cero perfecto sigue devolviendo el valor puntual. Es aditivo si hace falta en R3.
- La eficiencia de 0.968 en la integración viene de 1000 reps (SE ≈ 1.6 %). En R1, con 5000 reps, el valor fue 0.996, así que no hay contradicción.

```claims
[{"status":"unclear","text":"fisher.crb_map ahora usa zero_policy='limit' por defecto: crb_map(make_model(tcp_centers(50),LG fwhm300), linspace(-2,2,5),[0],100) = [1.61906,1.60858,1.60510,1.60858,1.61906] (sin salto central); con zero_policy='point' el centro da 1.80247 (S27); crb(p,0,100) por defecto sigue = 1.8024719"},
 {"status":"unclear","text":"fisher.crb(..., zero_policy='limit') en un punto con p_i<=p_min es igual a crb_limit centrado en ese punto (rel 1e-12, también en el cero periférico (0,25) y con N vectorial); con SBR=10 ambas políticas coinciden a 1e-12"},
 {"status":"unclear","text":"fisher.crb devuelve NaN (no inf) donde p_fn da NaN: haz gaussiano fwhm300, TCP L=50, r=(6000,0) -> p=[nan]*4, crb=NaN con 'point' y 'limit', crb_axes=NaN; Fisher singular finita sigue dando inf"},
 {"status":"unclear","text":"estimators.mle levanta ValueError (sin colgarse) para tol en {0,-1e-3,nan}, grid_step en {0,-1}, search_radius en {-5,inf}, chunk=0, mem_budget=0, conteos negativos o no finitos, y K de conteos != K de p_fn (mensaje explícito)"},
 {"status":"unclear","text":"estimators.mle memoria: search_radius=50, grid_step=0.2 (G=196321), 200 filas -> pico tracemalloc 134.8 MiB con mem_budget=64 MiB por defecto vs 309.6 MiB con el chunk fijo de R1; resultado igual al grueso refinado a 8.5e-5 nm; mem_budget=1 (chunk=1) y chunk=5 dan resultados idénticos bit a bit"},
 {"status":"unclear","text":"closed_forms.crb_tcp_center_limit vs fisher.crb_limit (make_beam power=c, L=100, N=100, fwhm=300): c=1.5 -> -4.39e-6, c=2 -> 1.20e-9, c=1.1 -> -2.25% con UserWarning (emitido para 1<c<1.5, no para c=1, 1.5, 2)"},
 {"status":"unclear","text":"montecarlo.run_mc sigma_err = sigma/(2 sqrt(n_valid)); consistente con la dispersión entre 60 semillas (LMS, n=400) dentro de ±30%; rmse documentado con normalización por eje sqrt(mean|dr|^2/2)"},
 {"status":"unclear","text":"Integración make_model -> fisher.crb_limit == copia local de test_acceptance._crb_center (dif. 0.0) y == crb_tcp_center_limit (<=2.3e-9) para (L,fwhm)=(50,300):1.6050956, (100,300):3.3634111, (50,360):1.5976973"},
 {"status":"unclear","text":"run_mc(mle) con make_model, SBR=10, L=50, N=100, 1000 reps, semilla 42: sigma=1.8973±0.0300, eficiencia sigma/crb_limit = 0.968 (CRB = S31 = 1.96006); lms(make_model) == lms_tcp a 1e-8 con y sin SBR=10"},
 {"status":"unclear","text":"camera: pixel 10 nm, ventana 121x121, sin fondo, N=400 -> CRB/(sigma/sqrtN) = 1.00042; a=100 nm 9x9 sin fondo -> 5.2045 nm (ideal 5.000, Balzarotti p.1); fondo total SBR inf/100/10/5/1 -> 5.2045/5.2898/5.7867/6.2514/9.4173 nm (monótono)"},
 {"status":"unclear","text":"camera: SBR_c=500 por píxel con K=81 == SBR_total=500/81=6.17, CRB idéntico (rel 0.0); sigma_PSF=100, a=100, 9x9, SBR_c=500, N=600 -> 4.964 nm (591 fotones para 5 nm, reproduce ~600 de Balzarotti Fig.3)"},
 {"status":"unclear","text":"MINFLUX L=50 fwhm300 Eq. S31 en el centro necesita 13.0/13.5/14.2/15.4/17.9 fotones para 5 nm (SBR inf/50/20/10/5) -> razón cámara/MINFLUX 45.5/44.0/41.8/38.5/33.0x; el '22x' del paper (600/27) usa ~27 fotones experimentales"},
 {"status":"unclear","text":"Suite completa python -m unittest discover -s tests: Ran 128 tests, OK (skipped=1, aceptación); incluye test_experiments.py de W3, sin test_vectorial.py; mis 7 archivos: 94 tests OK en 2.8 s; hash de test_acceptance.py intacto (7d198853...2303)"}]
```
