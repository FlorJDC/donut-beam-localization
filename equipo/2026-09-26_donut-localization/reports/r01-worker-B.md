# Ronda 1 — Worker B (inferencia: Fisher/CRB, estimadores, Monte Carlo)

## Qué hice
Escribí `src/donutloc/fisher.py`, `src/donutloc/estimators.py` y `src/donutloc/montecarlo.py`, con sus tests `tests/test_fisher.py`, `tests/test_estimators.py` y `tests/test_montecarlo.py`.
- No importé `beams`, `patterns` ni `photons`. Cada test define su propio modelo LG+TCP de unas 10 líneas (centro al final, SBR según Eq. S30).
- No creé `__init__.py`: ya existía (lo escribió A).
- No toqué archivos de otros workers ni `tests/test_acceptance.py`.

## API tal como quedó implementada
```
fisher.fisher_matrix(p_fn, r, N, h=1e-3, p_min=1e-12) -> (...,2,2)   # Eq. S11; se excluyen los términos con p_i <= p_min
fisher.crb(p_fn, r, N, **kw) -> (...)                                  # sqrt(tr F^-1 / 2), Eq. S13; F singular -> inf
fisher.crb_axes(p_fn, r, N, **kw) -> (sigma_x, sigma_y, isotropia)     # isotropía = sqrt(lmin/lmax) de Sigma, Eq. S14
fisher.crb_limit(p_fn, N, r_center=(0,0), r0=1e-3, n_dir=12, p_min=0.0) -> float   # h = r0*1e-2; definición del test de aceptación
fisher.crb_map(p_fn, xs, ys, N, **kw) -> (len(ys), len(xs))

estimators.neg_loglike(r, counts, p_fn) -> (...)                       # -sum n ln p, con p recortado a 1e-300
estimators.mle(counts, p_fn, search_radius, center=(0,0), grid_step=None, refine=True, tol=1e-4, chunk=4096) -> (2,) | (M,2)
estimators.lms(counts, p_fn, r_lin=(0,0), h=1e-3)                      # Eq. S48 con el Jacobiano numérico completo
estimators.lms_tcp(counts, L, fwhm, sbr=None, rotation=pi/2)           # Eq. S49-S50, centro al final, factor 1/s si hay sbr
estimators.mlms_tcp(counts, L, fwhm, beta=(1.27,3.8), sbr=None, rotation=pi/2)   # Eq. S51

montecarlo.run_mc(estimator, p_fn, r_true, N, n_rep, seed=42, mode="multinomial") -> dict
    claves: bias, std, sigma, rmse, sigma_err, n_rep, n_valid, seed, estimates, counts
montecarlo.sample_counts_from_p(p, N, n_rep, rng, mode)                # helper: rng.multinomial / rng.poisson
```

### Diferencias con el plan (a propósito)
1. **`crb_limit` usa `p_min=0.0` por defecto** (el plan no lo fijaba). Así se conserva siempre la p_centro chica pero finita a distancia r0, y el límite no cae en silencio al valor puntual si alguien pasa un r0 chico. Solo se excluyen los ceros exactos. `fisher_matrix` mantiene el default 1e-12.
2. **El refinamiento de `mle` no usa por defecto `scipy.optimize.minimize` por muestra.** Hace una búsqueda de patrón vectorizada de grueso a fino: una grilla local de 5×5 de semiancho 2·paso, luego el paso se divide por 2, hasta paso < `tol` = 1e-4 nm. Los candidatos fuera del disco se descartan. La razón es el rendimiento:
   - 2000 repeticiones tardan 0.47 s con la búsqueda de patrón y 7.8 s con Nelder-Mead por muestra, que para 500 repeticiones ya llega a 7.8 s.
   - La opción `refine="scipy"` sigue disponible: aplica Nelder-Mead sobre el óptimo de la grilla.
   - Las dos coinciden a menos de 2e-3 nm (test).
3. **Rama con máximos múltiples:** se refina el máximo **global** de la grilla. Los empates exactos se resuelven a favor del punto más cercano a `center`, porque la grilla está ordenada por radio. Con todos los conteos en cero se devuelve `center`. El refinamiento es local.
4. **`rmse` = sqrt(mean|r̂−r|²/2)**, normalizado por eje para que sea comparable con σ y con el CRB. No tuve acceso a la Eq. 4.2 de Masullo (no está en `docs/literature/`), así que la normalización está declarada en el docstring. `sigma_err` = σ/sqrt(2·n_valid), como pide el plan; es conservador.
5. **mLMS con `sbr`:** el factor 1/s se aplica a la parte LMS y el polinomio usa p̂_0 crudo.

## Hallazgo principal: el MLE es superEFICIENTE en el centro con cero perfecto (sin fondo)
El test del plan pedía std(MLE)/`crb_limit` ∈ [0.9, 1.2] en el centro con L=50, N=100 y SBR=inf. **Da 0.838 y no pasa.** No es un error del optimizador:
- La búsqueda de patrón, Nelder-Mead y una grilla bruta de 0.02 nm independiente dan el mismo σ (1.2671 frente a 1.2673 nm en 400 muestras; diferencia máxima 0.01 nm).
- Fuera del centro el MLE es eficiente, así que el problema es propio del centro.

Mi explicación, que el verificador debería revisar:
- Sin fondo, el modelo **no es regular** en r=0: p_centro ∝ r² y la información de Fisher es discontinua.
- Con n_centro≈0, el factor exp(−N·c·r²) de la verosimilitud actúa como un "prior" gaussiano que tira hacia el centro. El MLE queda sesgado hacia el centro (a r=(2,0), N=100 el sesgo es −0.29 nm) y su varianza en r=0 cae por debajo de la cota para estimadores insesgados.
- Como N·r² ~ L² en la escala σ ~ L/√N, la razón no depende de N.

Con fondo (SBR=10) el modelo es regular y el MLE es eficiente.

Consecuencia para `mle_efficiency_center` (aceptación, [0.9, 1.2]): **tiene que definirse con fondo**, por ejemplo SBR=10, donde da 0.988 ± 0.005. Otra opción es declarar explícitamente el caso regular. Sin fondo, frente al límite da 0.838 y frente al valor puntual S27 da 0.747, y ninguno de los dos entra en [0.9, 1.2]. **Lo decide el PI.**

Cambié el test del plan por tres tests:
- (a) SBR=10 en el centro: eficiencia ∈ [0.9, 1.2];
- (b) sin fondo en (10,0) con N=1000: ∈ [0.9, 1.15];
- (c) una regresión documentada del hallazgo: sin fondo en el centro ∈ [0.75, 0.92], más el sesgo hacia el centro en (2,0).

## Resultados de los tests (salida real)
- `python -m unittest discover -s tests -p "test_fisher.py"` → `Ran 14 tests in 0.036s  OK`
- `python -m unittest discover -s tests -p "test_estimators.py"` → `Ran 14 tests in 0.381s  OK`
- `python -m unittest discover -s tests -p "test_montecarlo.py"` → `Ran 7 tests in 1.241s  OK`
- La suite completa (`discover -s tests`, que incluye los tests de A) da `Ran 84 tests in 1.635s  OK (skipped=1)`. El test omitido es el de aceptación, porque todavía falta `data/paper_numbers.json`.

## Números (scripts en el scratchpad; semilla 42; L=50, N=100, fwhm=300; MLE con search_radius=L)
| caso | σ_MC (nm) | n_rep | crb_limit | CRB puntual | σ/crb_limit | tiempo |
|---|---|---|---|---|---|---|
| SBR=inf, centro | 1.3456 ± 0.0067 | 20000 | 1.60510 | 1.80247 (S27) | 0.8383 ± 0.0042 | 4.26 s |
| SBR=10, centro | 1.9362 ± 0.0097 | 20000 | 1.96006 | 1.96006 (S31) | 0.9878 ± 0.0049 | 4.31 s |

Barrido en N con 5000 repeticiones (σ/crb_limit):

| | N=100 | N=1000 | N=10000 |
|---|---|---|---|
| SBR=inf | 0.837 | 0.831 | 0.830 |
| SBR=100 | 0.862 | 0.963 | 0.999 |
| SBR=10 | 0.996 | 0.989 | 1.006 |

Fuera del centro, sin fondo (5000 repeticiones, σ/CRB puntual):

| r (nm) | N=100 | N=1000 |
|---|---|---|
| (2,0) | 0.966 | 1.011 |
| (5,0) | 1.044 | 1.010 |
| (10,0) | 1.101 | 0.999 |
| (0,10) | 1.045 | 1.010 |

Con L=50, N=100, sin fondo, el LMS (`lms_tcp`, 20000 repeticiones) da σ/crb_limit = 1.121 y σ/S27 = 0.998.

## Afirmaciones verificables
```claims
[{"status":"proposed","text":"fisher.crb(p_fn,[0,0],N=100) para LG+TCP con L=50, fwhm=300 y sin fondo da 1.8024719 nm, igual a Eq. S27 (1.8024719) con error relativo <1e-8; también coincide con L=100 a <1e-4 (el término con p=0 exacto se excluye)."},
 {"status":"proposed","text":"Con SBR=10 fijo (Eq. S30), fisher.crb en el centro coincide con Eq. S31 a <1e-4 (1.96006 nm con L=50, N=100, fwhm=300), y crb_limit da el mismo valor (con fondo el CRB es continuo)."},
 {"status":"proposed","text":"fisher.crb_limit(p_fn,100) con L=50 y fwhm=300 da 1.6050956 nm, idéntico (diferencia relativa 0.0) a tests/test_acceptance._crb_center(50,100,300); también coincide a <1e-3 con una reimplementación inline para (L,fwhm) = (100,300) y (50,360)."},
 {"status":"proposed","text":"crb_limit / CRB puntual = 0.8905 con fwhm=300 y L=50, y 2/sqrt(5)=0.8944 a <1e-3 en el régimen cuadrático (fwhm=1e6)."},
 {"status":"proposed","text":"El CRB escala exactamente como N^-1/2 (rtol 1e-12) y como L para L<<fwhm (crb(L=10)/crb(L=5)=2 a menos de 1%, tanto el valor puntual como el límite)."},
 {"status":"proposed","text":"estimators.mle recupera la posición verdadera a <5e-3 nm con conteos esperados (1e7 fotones) en 5 posiciones dentro del TCP con L=50, y con SBR=5, fwhm=360, L=100 y rotación 0.3."},
 {"status":"proposed","text":"lms_tcp (Eq. S49-S50, centro al final, 1/s con fondo) coincide con el lms general (Eq. S48, Jacobiano numérico) a rtol 1e-5 en 4 configuraciones (con y sin SBR, fwhm finito e infinito, rotaciones pi/2, 0 y 1)."},
 {"status":"proposed","text":"Sin el factor 1/s, el LMS con SBR=10 contrae la estimación en s=10/11 (test con conteos esperados en r=(0.4,-0.3))."},
 {"status":"proposed","text":"MLE en el centro con L=50, N=100, fwhm=300 y SBR=10: sigma_MC=1.9362±0.0097 nm, eficiencia sigma/crb_limit=0.988±0.005 (n_rep=20000, semilla 42)."},
 {"status":"proposed","text":"MLE en el centro con L=50, N=100, fwhm=300 y SIN fondo: sigma_MC=1.3456±0.0067 nm, sigma/crb_limit=0.838±0.004 y sigma/S27=0.747±0.004 (n_rep=20000, semilla 42). Es superEFICIENTE, fuera de [0.9,1.2]; la razón es ~0.83 para N=100, 1e3 y 1e4 (5000 repeticiones cada uno)."},
 {"status":"proposed","text":"La superEFICIENCIA sin fondo en el centro viene de que el modelo no es regular y de un sesgo del MLE hacia el centro: en r=(2,0) con N=100 el sesgo en x es -0.29 nm (5000 repeticiones). En r=(10,0) con N=1000, sigma/CRB=0.999. Una grilla bruta independiente de 0.02 nm da el mismo sigma (1.2671 frente a 1.2673 nm, 400 muestras)."},
 {"status":"proposed","text":"Rendimiento: run_mc con mle (grilla + búsqueda de patrón vectorizada, tol=1e-4 nm) tarda 0.47 s para 2000 repeticiones y 4.1-4.3 s para 20000 (K=4, search_radius=50). refine='scipy' tarda 7.8 s para 500 repeticiones."}]
```

## Lo que queda abierto
- Hay que decidir la definición de `mle_efficiency_center` para `paper_numbers`: con fondo (recomiendo SBR=10), porque sin fondo el MLE queda por debajo del CRB. Esto va al PI.
- No verifiqué la normalización de la Eq. 4.2 de Masullo para el RMSE (no tengo el texto).
- `sigma_err` es la fórmula conservadora del plan. Para bajar el error de la eficiencia por debajo de 2 % (backlog R2), 20000 repeticiones ya alcanzan: el error es 0.5 %.
