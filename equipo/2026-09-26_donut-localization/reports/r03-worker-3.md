# Ronda 3 — Worker 3: figuras 5 a 8 (MC)

## Qué hice
- Escribí cuatro scripts. Usan solo la API pública de `donutloc`, más `_paperconfig` y `_paperstyle`. Todos aceptan `--quick` y `--no-cache`, cachean en `data/mc/fig{5..8}_*.npz` (gitignored), imprimen `NUMBER` con `report()` y escriben `data/fig{5..8}_summary.json`:
  - `scripts/fig_5_estimators.py` → `paper/figures/fig5_estimators.pdf`
  - `scripts/fig_6_iterative.py` → `paper/figures/fig6_iterative.pdf`
  - `scripts/fig_7_zero_depth.py` → `paper/figures/fig7_zero_depth.pdf`
  - `scripts/fig_8_misalignment.py` → `paper/figures/fig8_misalignment.pdf`
- Rendericé cada PDF a PNG con PyMuPDF (`work/figcheck/render.py`, PNG en `work/figcheck/`) y los revisé a ojo. Corregí solapamientos de leyendas y etiquetas, ticks log sin rótulo y etiquetas cortadas.
  - Todos los PDF son vectoriales, de 7.0 in de ancho (doble columna), Okabe-Ito de `COLORS`, con unidades y panel labels.
- Los logs de las corridas finales están en `work/figlogs/fig{5,7,8}.log`.
- SE de σ: bootstrap con `montecarlo.bootstrap_sigma_se` (N_BOOT=2000, semilla 42), que W1 ya había agregado antes de la corrida final; el fallback gaussiano no se usó.
- En Fig 8, `sigma_se` es el SE entre patrones: verifiqué el fix de W1 en `src/donutloc/experiments.py`, con `sigma_se_within` para la fórmula vieja, y recién después corrí con `--no-cache`.
- No escribí tests: `tests/` es de W1 y los scripts de figuras no son código de biblioteca. Los cuatro scripts corren en `--quick` y en modo completo, y la ruta de fallback de Fig 6 (sin json) también: la probé con 1000 reps y reproduce el protocolo, σ(N=1000)=0.506.

## Tiempos de cómputo (corrida final, esta máquina)
| figura | cómputo | notas |
|---|---|---|
| Fig 5 | 355 s | sweep 137 s + centro 45 s + sesgo cerca del centro 172 s |
| Fig 6 | 3 s | lee `data/iterative_sweep.json` de W1: 10000 reps, no quick; su sweep tardó 128 s |
| Fig 7 | 7 s | determinista |
| Fig 8 | 101 s | corrida con `--no-cache` tras el fix de W1 |

Con el caché, rehacer las 4 figuras tarda menos de 10 s.

## Números (valor ± SE)
### Fig 5: estimadores (L=50, N=100, fwhm=300, semilla 42)
**Centro, SBR=10** (MLE con `search_radius=L`, 10000 reps):
- `mle_efficiency_center` = **0.9912 ± 0.0051**. σ = 1.9429 ± 0.0101 nm y CRB = 1.9601 nm.
  - Es consistente con r1: 0.988 ± 0.005 con 20000 reps.
  - `compute_paper_numbers.py` todavía no existe. Para que coincida, W1 debe usar `search_radius=L`, `run_mc(seed=42)` y n_rep=10000.

**Centro, sin fondo, MLE** (10000 reps por N):

| N | σ/CRB_lim | σ/S27 |
|---|---|---|
| 100 | 0.8397 ± 0.0041 | 0.7478 ± 0.0037 |
| 300 | 0.8258 ± 0.0041 | 0.7354 ± 0.0036 |
| 1000 | 0.8314 ± 0.0042 | 0.7404 ± 0.0038 |
| 3000 | 0.8316 ± 0.0041 | 0.7405 ± 0.0037 |
| 10000 | 0.8286 ± 0.0040 | 0.7379 ± 0.0036 |

- En todo el rango, el MLE queda en σ/CRB_lim ≈ 0.83 y σ/S27 ≈ 0.74.
- El LMS sin fondo da σ/S27 = 1.000/0.992/1.001/1.002/0.998 (±0.005): coincide con S27, como verificó V16.

**Sesgo del MLE sin fondo en r=(2,0), N=100:** `mle_nobg_bias_x_r2_nm` = **−0.3445 ± 0.0077 nm** (40000 reps). Coincide con V14 (−0.338 ± 0.010).
- El sesgo en x, en función de x0 = 0.5/1/2/3/4/6/8/12 nm, vale −0.13/−0.24/−0.34/−0.33/−0.26/−0.14/−0.05/+0.18 nm: hacia el centro para x0 ≲ 8 nm y hacia afuera en 12 nm.

**Corte en x con SBR=10** (10000 reps por punto; MLE con radio 2L porque x0 llega a L):

| estimador | x0=0: sesgo_x (nm) | x0=0: σ/CRB | x0=25: sesgo_x (nm) | x0=25: σ/CRB | x0=50: sesgo_x (nm) | x0=50: σ/CRB |
|---|---|---|---|---|---|---|
| MLE | 0.04 ± 0.02 | 0.991 | 0.68 ± 0.07 | 1.035 | −2.34 ± 0.13 | 1.052 |
| LMS | — | 1.000 | −14.35 | 0.278 | −42.57 | 0.136 |
| mLMS | — | 1.357 | −5.19 | 0.514 | −34.38 | 0.284 |

- LMS y mLMS dan σ/CRB < 1 porque están muy sesgados (encogen hacia el centro). Eso no es superioridad: el caption tiene que decirlo.
- **Hallazgo:** con N=100 el MLE tiene σ/CRB de hasta 1.51 en x0=15 nm, dentro del TCP.
  - La causa es una cola no asintótica: ~4 % de estimaciones caen a más de 4 CRB, por máximos secundarios de la verosimilitud.
  - Con N=1000 en el mismo punto, σ/CRB = 0.999 (control aparte, 2000 reps, no está en la figura).
  - El paper no debe decir "MLE eficiente en todo el TCP" con N=100.

### Fig 6: iterativo frente a cámara (de `data/iterative_sweep.json`, 10000 reps, SE bootstrap)
**σ final por N_total:**

| N_total | σ (nm) | σ/cámara |
|---|---|---|
| 250 | 1.158 ± 0.014 | 0.1831 ± 0.0022 |
| 500 | 0.7359 ± 0.0042 | 0.1646 |
| 1000 | 0.5077 ± 0.0027 | 0.1606 ± 0.0008 |
| 2000 | 0.3544 ± 0.0019 | 0.1585 |
| 4000 | 0.2495 ± 0.0013 | 0.1578 |
| 8000 | 0.1763 ± 0.0010 | 0.1577 |

- `iterative_sigma_nm` = 0.5077 ± 0.0027 nm; `camera_sigma_nm` = 3.1623 nm (ideal).

**Pendiente log-log:**
- En todo el rango: −0.5365, IC95 [−0.5424, −0.5312].
- Con N ≥ 500: −0.5148, IC95 [−0.5196, −0.5096].
- Ninguna de las dos es compatible con −0.5: el cociente con la cámara todavía baja hasta N≈4000.
- Rango del cociente: 0.158–0.183.

**N_total=1000:**
- SBR=10: 0.6165 ± 0.0033 nm.
- Adaptivo: 0.4795 ± 0.0027 nm.
- CRB con los 1000 fotones en L=25: 0.2509 nm. El iterativo queda a 2.02× de ese valor.
- `recenter=False`: 8.14 ± 0.05 nm. Se dibuja solo como una × gris rotulada "no re-centring (search-disk artefact)", sin línea ni leyenda de física.

### Fig 7: profundidad finita del cero (gaussiano, N=100, determinista)
| eps | L_opt (nm) | CRB_opt (nm) | L_opt/(fwhm√eps) |
|---|---|---|---|
| 0.002 | 10.484 | 0.7458 | 0.7815 |
| 0.01 | 23.280 | 1.6787 | 0.7760 |
| 0.05 | 50.337 | 3.8784 | 0.7504 |
| 0.15 | 80.810 | 7.2867 | 0.6955 |

- Idéntico a W12.
- Con eps=1e-4 el cociente vale 0.7828 (límite de eps chico).
- Con eps=0 la divergencia está en fwhm/√ln2 = 360.34 nm. El centro con eps=0 y L=50 da 1.6051 nm, igual a `crb_center_lg_L50_N100_nm`.

**Panel (c): eps=0.002, L=100:**
- Centro: 3.8768 nm.
- Mínimo según la dirección:
  - 0°: 3.6087 nm en r=6.7 nm;
  - 45°: 3.5331 nm en 8.0 nm;
  - 90°: 3.6937 nm en 5.1 nm.
- La escala √(eps/(e a)) = 4.887 nm se marca como escala, no como argmin (caveat de V10).

### Fig 8: desalineación (MIS: L=100, N=500, SBR=10, 50 patrones × 200 reps, `--no-cache` post-fix)
| δ (nm) | ingenuo \|sesgo\| centro | honesto \|sesgo\| | piso MC | σ ingenuo | σ honesto | CRB honesto |
|---|---|---|---|---|---|---|
| 0 | 0.187 ± 0.013 | 0.187 | 0.164 | 1.846 ± 0.008 | 1.846 | 1.863 |
| 2 | 1.551 ± 0.082 | 0.167 ± 0.014 | 0.164 | 1.853 ± 0.009 | 1.849 ± 0.011 | 1.859 |
| 5 | 3.919 ± 0.211 | 0.162 ± 0.014 | 0.163 | 1.876 ± 0.012 | 1.843 ± 0.016 | 1.842 |
| 10 | 8.069 ± 0.431 | 0.167 ± 0.013 | 0.164 | 1.987 ± 0.028 | 1.845 ± 0.026 | 1.840 |

**En (L/4,0) con δ=10:**
- Ingenuo: |sesgo| = 8.22 ± 0.68 nm y σ = 3.13 ± 0.40 nm.
- Honesto: |sesgo| = 0.29 ± 0.03 nm y σ = 2.64 ± 0.06 nm.
- CRB honesto: 2.59 nm.

**Pendiente del sesgo ingenuo** (WLS por el origen, δ ≥ 2):
- Solo con `MIS_DELTAS`: centro 0.788 ± 0.024, (L/4,0) 0.799 ± 0.032.
- Con los 8 δ, incluidos 1/3.5/7.5/15: 0.797 ± 0.017 y 0.801 ± 0.023.

**Otros resultados:**
- Con δ=0, honesto e ingenuo coinciden bit a bit.
- El honesto en el centro tiene σ/CRB = 0.991–1.005 para todos los δ, y su sesgo queda en el piso MC.
- El `sigma_se` entre patrones crece con δ: en el centro, 0.0079 con δ=0 y 0.0276 con δ=10. `sigma_se_within` queda en ~0.0093–0.0099. Es el efecto esperado del fix.

## Lo que no quedó resuelto
1. **SE < 2 % en Fig 8.** La meta no se cumple para el sesgo ingenuo por δ.
   - SE relativo: 5.3 % con δ=10 en el centro, 5.4 % con δ=2 y 8 % en (L/4,0).
   - No es ruido MC por repetición: es la dispersión intrínseca entre patrones aleatorios de desalineación, con n_patterns=50 canónico.
   - La pendiente sí queda a 2.2–3 %.
   - Para llegar al 2 % en cada δ habría que subir `MIS.n_patterns` a ~350, con un costo de ~6 min. Es decisión del PI y W1 (`_paperconfig`), y no la toqué.
   - El resto de los números MC quedan todos < 1.2 %. El peor es σ(N=250) del iterativo, con 1.2 %.
2. **Caveat "sesgo ingenuo ≈ 0.75 δ".** Con este sorteo de 50 patrones la pendiente es 0.79–0.80 ± 0.02–0.03.
   - Es compatible con 0.75 dentro de la dispersión entre sorteos de patrones que midió el verificador en r2 (±10 % con 20 patrones).
   - Propongo reformular el claim como "≈ 0.8 δ (0.79 ± 0.02 para el sorteo canónico)" o "0.75–0.80 δ".
3. **Consistencia de Fig 5 con `compute_paper_numbers.py`.** El script todavía no existía al cierre.
   - El punto centro con SBR=10 usa `search_radius=L`, `run_mc(seed=42)` y n_rep=`N_REP_MLE`.
   - Si W1 usa otro radio, el valor de `mle_efficiency_center` puede diferir por ruido MC en el tercer decimal.
   - El verificador debe cruzar `data/fig5_summary.json` con `data/paper_numbers.json`.
4. **Fig 6 depende de `data/iterative_sweep.json`.** Si W1 regenera el sweep, basta con volver a correr `fig_6_iterative.py` (3 s).

## Afirmaciones verificables (una por línea)
```claims
[{"status":"proposed","text":"fig_5_estimators.py: MLE at the TCP centre, L=50, N=100, SBR=10, search_radius=L, run_mc seed 42, 10000 reps: sigma/crb_limit = 0.9912 ± 0.0051 (bootstrap), sigma = 1.9429 nm, CRB = 1.9601 nm"},
 {"status":"proposed","text":"fig_5: background-free centre MLE (search_radius=L, 10000 reps, seed 42) sigma/crb_lim = 0.840/0.826/0.831/0.832/0.829 (±0.004) and sigma/S27 = 0.748/0.735/0.740/0.741/0.738 for N=100/300/1000/3000/10000; LMS sigma/S27 = 1.000/0.992/1.001/1.002/0.998 (±0.005)"},
 {"status":"proposed","text":"fig_5: background-free MLE bias along x at r=(2,0), N=100, 40000 reps, seed 42: -0.3445 ± 0.0077 nm; the bias is toward the centre for x0<=8 nm (-0.13 at 0.5, -0.24 at 1, -0.33 at 3, -0.26 at 4, -0.14 at 6, -0.05 at 8) and +0.18 nm at x0=12"},
 {"status":"proposed","text":"fig_5: with SBR=10, N=100, L=50 the MLE (search radius 2L) has sigma/CRB up to 1.51 at x0=15 nm (inside the TCP) because of a ~4% tail of estimates beyond 4 CRB; at N=1000 the same point gives sigma/CRB=0.999 (2000 reps)"},
 {"status":"proposed","text":"fig_5: SBR=10, L=50, N=100: LMS bias_x = -14.35 nm at x0=25 and -42.57 nm at x0=50; mLMS -5.19 and -34.38 nm; their sigma/CRB < 1 (0.28/0.14 LMS, 0.51/0.28 mLMS) is due to the shrinkage, and mLMS has sigma/CRB = 1.357 at the centre"},
 {"status":"proposed","text":"fig_6 (from data/iterative_sweep.json, 10000 reps): iterative sigma at N_total=1000 = 0.5077 ± 0.0027 nm, ratio to the ideal camera 0.1606; slope -0.5365 [95% CI -0.5424,-0.5312] over 250..8000 and -0.5148 [-0.5196,-0.5096] for N>=500 (neither contains -0.5); ratio range 0.158-0.183; SBR=10 0.6165 ± 0.0033; adaptive 0.4795 ± 0.0027; no-recentre 8.14 nm drawn only as a grey artefact marker"},
 {"status":"proposed","text":"fig_7 (eps_L_sweep/optimal_L, gaussian pedestal, N=100, fwhm=300): L_opt = 10.484/23.280/50.337/80.810 nm, CRB_opt = 0.7458/1.6787/3.8784/7.2867 nm, L_opt/(fwhm sqrt eps) = 0.7815/0.7760/0.7504/0.6955 for eps = 0.002/0.01/0.05/0.15, and 0.7828 at eps=1e-4"},
 {"status":"proposed","text":"fig_7c: eps=0.002, L=100, N=100: centre CRB 3.8768 nm; the minimum along 0/45/90 deg is 3.6087 nm at 6.7 nm, 3.5331 nm at 8.0 nm and 3.6937 nm at 5.1 nm; sqrt(eps/(e a)) = 4.887 nm"},
 {"status":"proposed","text":"fig_8 (misalignment_study, MIS = L100/N500/SBR10/50 patterns x 200 reps, seed 42, post-fix, --no-cache): naive centre |bias| = 1.551 ± 0.082 / 3.919 ± 0.211 / 8.069 ± 0.431 nm for delta = 2/5/10; honest 0.167/0.162/0.167 nm against an MC floor of 0.164 nm; naive sigma 1.846 -> 1.987 ± 0.028 nm; honest sigma/CRB 0.991-1.005"},
 {"status":"proposed","text":"fig_8: WLS slope through the origin of the naive |bias| vs delta (delta>=2, MIS_DELTAS) = 0.788 ± 0.024 (centre) and 0.799 ± 0.032 at (L/4,0); over 8 deltas up to 15 nm 0.797 ± 0.017 / 0.801 ± 0.023 — the text caveat '≈0.75 delta' should read ≈0.8 delta (0.75-0.80)"},
 {"status":"proposed","text":"fig_8: the naive per-delta |bias| has a relative SE of 5-8% (pattern-to-pattern spread with n_patterns=50), above the 2% target; ~350 patterns would be needed"}]
```
