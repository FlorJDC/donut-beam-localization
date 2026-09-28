# r01 — verifier (worker B: temporización p-MINFLUX)

No leí `reports/r01-worker-B.md`. Todo lo reproduje con código propio que no importa `donutloc`
(beam LG, TCP, Fisher, MLE grilla+Nelder-Mead, telegrafía y conteo de pulsos escritos de cero):
`work/verify/verify_pminflux_B.py` (+ `.out`), `work/verify/verify_extra.py` (+ `.out`).
`work/verify/diag_crb.py` sí importa `donutloc.fisher`, pero solo para diagnosticar una discrepancia de convención.
Cómputo total: unos 3 min.

## 1. Matriz de cross-talk
Forma cerrada: el fotón del pulso en la ranura s llega en s·T + t, con t ~ Exp(τ); su ventana es
floor(((sT+t) mod P)/T). P(m) = Σ_k ∫_{(m+4k)T}^{(m+4k+1)T} e^{-t/τ}/τ = (1−q)q^m/(1−q⁴).
τ = 5 ns: [0.9180, 0.0754, 0.0062, 0.0005]. El MC de 2·10⁶ fotones da 0.918/0.0754/0.0061/0.0005.
La fuga va a la ventana **siguiente**: M[1,0] = 0.0754 y M[3,0] = 0.0005, y el JSON usa la misma
convención (fila = exposición observada). Es doblemente estocástica. VERIFICADO.

## 2. CRB con p_obs = M p
- En el centro, relativo al valor puntual S27 (38.305 nm; el límite r→0 da 33.634): obtengo 1.0000, 1.0019,
  1.0159, 1.0469, 1.0928. Coincide exactamente. VERIFICADO.
- Disco r ≤ 50 (razón de medias): con una malla fina desplazada (paso 0.5 nm, 31428 puntos, que no cae
  sobre ningún cero) obtengo 1.00002, 1.0060, 1.0399, 1.1044, 1.1915. Para τ = 2..5 coincide con 1.0063/1.0398/1.1043/1.1914
  (diferencias ≤ 3e-4). **Para τ = 1, el 1.0004 es un artefacto de la malla**: la malla de 2.5 nm centrada en el
  origen contiene el cero exacto de la dona 0 en (0, 50). Ahí el ideal con `zero_policy="limit"` vale 95.8 nm,
  y con M la fuga (3.7e-6) deja el valor puntual (≈124 nm). Solo ese punto aporta ≈ +3.6e-4. En el continuo
  la razón es 1.0000. Con la misma malla y la política "point" en todos los puntos obtengo 1.0000, 1.0059,
  1.0394, 1.1039, 1.1909. Es decir, la cuarta cifra depende de la convención en los ceros.
- Independencia del orden cíclico: los órdenes (0,1,2,3), (0,2,1,3) y (3,0,1,2) dan resultados idénticos. (1,0,3,2) difiere
  en ≤1e-4, por la asimetría de la malla cuadrada. VERIFICADO.

## 3. MLE que ignora M (conteos esperados, y = 0)
En la malla de 11 puntos, máximo |sesgo|: τ = 3 → 3.410 nm (en x = 0, ambos órdenes). τ = 5, orden A → 8.448 nm (x = 0).
τ = 5, orden B → 12.971 nm (x = 50). En una malla de 1 nm el orden B alcanza 12.993 nm en x = 49. El MLE que incluye M
recupera la posición exacta (error < 1e-5 nm). VERIFICADO.

## 4. Flickering (σ_fl poblacional en el centro)
Resultados con mi telegrafía (estado estacionario; realizaciones sin señal excluidas; MLE restringido a r ≤ 75 nm):

| Esquema | Realizaciones | σ_fl propio (nm) | Afirmación / JSON (nm) |
|---|---|---|---|
| Secuencial r = 1 | 2000 | 16.29 | ≈ 16.1 (JSON 16.07) |
| Secuencial r = 5 | 2000 / 1500 | 6.48 / 6.91 | JSON 6.63 |
| Secuencial r = 25 | 2000 | 2.15 | ≈ 2.2 |
| Entrelazado | 2000 / 6000 | 0.014 / 0.011 | ≈ 0.02 (JSON 0.018) |

El secuencial cae dentro del ruido de muestreo (±5 %). En el entrelazado, el orden de magnitud es correcto
(10⁻² nm, unas 10³ veces menor que el secuencial con r = 1). Pero el valor no está convergido: lo dominan las
pocas realizaciones con muy poco tiempo "on", y varía entre 0.011 y 0.018 según la muestra. Conviene
escribir "≲ 0.02 nm" y no "≈ 0.02".

## 5. SimuFLUX (STD² = σ_fl² + σ²_ctrl)
Con conteos multinomiales y N = 100:

| Esquema | STD (nm) | ctrl (nm) | √(STD² − ctrl²) (nm) | σ_fl poblacional (nm) |
|---|---|---|---|---|
| Secuencial r = 1 | 16.44 | 2.74 | 16.21 | 16.29 |
| Secuencial r = 5 | — | — | 6.46 | 6.48 |

Independencia de N (r = 5, mismas realizaciones): con N = 100 obtengo 6.917 nm y con N = 1000, 6.918 nm, frente a un σ_fl
poblacional de 6.913 nm. En el entrelazado, la resta da 0.75 ± ~0.4 nm, compatible con 0 (el JSON trae un
σ_fl² negativo). VERIFICADO.

Observación aparte (no es una afirmación del worker): en el centro, la STD de control con N = 100 es 2.7 nm. Queda por debajo
de CRB/√N = 3.36 nm (límite) o 3.83 nm (puntual). El MLE en el cero exacto de la dona central no es eficiente ni
insesgado (casi siempre n_c = 0 y eso lo atrae al origen). Si el informe compara σ con el CRB en el centro,
debe aclararlo.

## Tests
`python3 -m unittest tests.test_pminflux`: 15 tests, OK (0.17 s).

## Supuestos puestos a prueba
- Dirección de la fuga: correcta (hacia la ranura siguiente).
- Ventana: [sT, (s+1)T), alineada con el pulso. Si hubiera IRF o las ventanas estuvieran desplazadas, todos los números de cross-talk cambiarían.
  Hoy son válidos solo para ese modelo, que además supone excitación no saturada.
- σ por eje = sqrt((σx²+σy²)/2): usé la misma convención y los números coinciden. La convención radial daría √2 veces más.

## No resuelto
- Si la razón en el disco para τ = 1 se reporta como 1.0004, está inflada por la convención en los ceros. En el continuo vale 1.0000.
- El σ_fl del entrelazado no está convergido en su tercera cifra decimal.
