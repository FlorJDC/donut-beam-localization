# r06 — verificador (comparación con SimuFLUX, commit 7d41d4d)

Alcance: las 4 afirmaciones con `verified_round: 6` en `structure/claims.json`
(`bgphys_estimators_offcentre`, `bgphys_iterative`, `naive_misspecified_estimators`,
`simuflux_point_value_convention`), sus números en `data/paper_numbers.json` y el texto nuevo en
`paper/sections/{estimators,iterative,nonidealities,discussion,model,open_points}.tex`.
No leí `docs/private/` ni el repo de la autora. No toqué código ni paper.

## Camino propio
`work/verify/r06/vcore.py`: implementación mínima sin `donutloc`: dona LG, TCP (centro al final),
p con SBR fija (S30) o fondo constante por exposición (S28), pedestal constante `I_LG+eps`, Fisher con
derivadas por paso complejo (exactas a precisión de máquina), MLE por lotes con grilla global en el
disco y 7-10 niveles de grillas locales (hasta ~1e-4 nm), LMS S49-S50 con 1/s, mLMS (1.27, 3.8).
Scripts: `v_deterministic.py` (CRB, SBR, SimuFLUX, MLE sin ruido con pulido Nelder-Mead),
`v_polish.py` (Newton con gradiente complejo), `v_mc.py` (MC, semilla 20260930, 20000 reps en los
puntos de la fig 5, 10000 en el iterativo), `v_finiteN.py` (sesgo del MLE ingenuo a N finito).
Resultados: `deterministic.json`, `mc_fig5.json`, `mc_iter.json`, `finiteN.json`.

## Resultados
**Fondo constante por exposición, puntos de la fig 5 (L=50, N=100).** SBR 22.5814 (x0=25) y 57.5202
(x0=50): coinciden exactamente. CRB 5.98543/5.77286 (x25) y 12.19679/11.38952 (x50): coinciden a 1e-9.
MC propio frente al publicado, 24 cantidades (sesgo x y σ/CRB, 3 estimadores, 2 posiciones, 2 modelos):
todas dentro de |z| ≤ 1.37 del SE combinado. Cambios phys − fixed en mi MC: sesgo x ≤ 1.26 nm
(mLMS x50), σ/CRB ≤ 6.3 % (LMS x50). Se cumple "< 1.3 nm" y "< 8 %". Los valores "fixed" son
idénticos, bit a bit, a `data/mc/fig5_sweep.npz` (misma semilla y código): "reproduce la fig 5
exactamente" es cierto.

**Iterativo (N_tot=1000, 4 iteraciones L=150→25 geométricas, 250 fotones cada una, MLE en 0.75L).**
SBR del centro: con ajuste en L=150 → [10, 3.417, 1.074, 0.3287]; con ajuste en L=25 → primera 304.18.
Ambos coinciden con la forma cerrada de S32, L² e^{-ln2 L²/fwhm²}. σ en mi MC: 2.228±0.015
(publicado 2.232±0.014), 0.609±0.003 (0.611±0.003), SBR fija 10: 0.613±0.003 (0.616±0.003), sin
fondo: 0.507±0.003 (0.508±0.003). Todo dentro de 0.7 SE.

**Estimadores mal especificados (sin ruido).** |sesgo| del MLE ingenuo (pulido Newton, gradiente
~1e-16): eps=0.002 → 0.268712 / 0.422056 nm; eps=0.01 → 1.38869 / 2.28936 nm; SBR supuesta 20 →
0.35792 / 0.73961; sin fondo → 0.91572 / 1.64286. CRB 2.40116/4.46073/2.87928/5.20051. N_eq
429.9 / 516.0 / 7985 / 11170.5. Desplazamiento extra del LMS: −0.3411/−0.3005/−1.4634/−1.3517
(exacto). MLE honesto: 5e-7 nm. Con eps=0.05 el MLE ingenuo llega al borde del disco en x0=20
(|r̂|=37.5), pero no en x0=10 (r̂=(14.4, 5.0)); la advertencia es cierta solo en parte.

**SimuFLUX (fwhm=310, L=75, N=100, sin fondo).** Valor puntual 2.763782 (Fisher sin la exposición
central y forma cerrada propia), límite 2.448757 (r0=1e-4…1e-2 nm, en 3 direcciones). Cociente 1.129.
El valor puntual es el mismo para un hexágono sin centro, así que no depende del patrón de anillo.
Diferencia central al estilo del Methods de [MR] (x±ε/2, p+1e-4, eje x): 2.7648 (N=100) y 0.8743
(N=1000). Esto concuerda con ~2.8 / ~0.9 leídos de la Fig. 2e / SI Fig. 8 **si L=75**; con L=70/80 el
valor puntual sería 2.57/2.97. Factor e: la dona 4ln2 (r/fwhm)² e^{-4ln2 r²/fwhm²} tiene máximo 1/e
y la misma potencia que la gaussiana de pico unidad (cociente 1.0000000001), y la nuestra es
exactamente e veces esa: ε = e·zerooffset es correcto **dada** esa forma.

## Problemas encontrados (ninguno invalida una afirmación)
1. **Precisión impresa de más.** `naive_eps0p002_mle_bias_abs_x20_nm` vale 0.42200 en el JSON y
   0.422056 con precisión de máquina: el texto dice 0.4220 y debería decir 0.4221. En consecuencia
   `naive_eps0p002_mle_N_bias_eq_crb_x20` es 11170, no "11173". Su tolerancia del MLE (~3e-5 nm)
   no alcanza para el 4.º dígito. Conviene imprimir 3 cifras (0.422; 1.12×10⁴). El
   "3.052×10⁻⁵ nm" del MLE honesto es la tolerancia del optimizador, no un resultado: mejor
   "< 10⁻⁴ nm".
2. **"The bias does not decrease with N"** vale para el sesgo asintótico (pseudo-verdadero), que
   por construcción no depende de N. En MC (eps=0.01, x0=20, 8000 reps) el sesgo medio del MLE
   ingenuo es 3.84 (N=100), 3.61 (500), 2.67 (2000) y 2.32 nm (20000). Sí baja con N y tiende a
   2.29. El honesto da 1.72 / 0.53 / 0.10 / 0.01 nm. La conclusión (el pedestal ignorado domina por
   encima de unos cientos de fotones) se sostiene, y a N finito el efecto es incluso mayor. Aun así,
   el texto debería decir que 2.29 nm es el límite de N grande.
3. **Discusión:** "a background constant per exposure changes the off-centre biases … by less than
   1.3 nm and σ/CRB by less than 8 %" omite "en los dos puntos comparados (x0=25, 50; L=50, N=100)",
   que sí dice `estimators.tex`. Además solo se reporta el sesgo en x: el sesgo en y del MLE es del
   mismo orden (−1.3 a −2.7 nm), aunque su cambio también es < 0.4 nm en mi MC.
4. **Iterativo, frase de interpretación.** Subir la potencia al achicar L mantiene la SBR solo si el
   fondo no escala con la potencia (cuentas oscuras, luz ambiente), no si es autofluorescencia o
   fondo por excitación. Además, rampas de 1→3 o 1→6 no compensan una caída de ~30× (0.33 frente a
   10). El texto sugiere que las secuencias de SimuFLUX "mantienen" la SBR, y no pude comprobar los
   `pwrFactor` porque el código de SimuFLUX no está local. Sugiero "partially compensate" y declarar
   el supuesto sobre el fondo.
5. `calculateCRBdirect` y `PsfDonut2D`: el código de SimuFLUX no está en la máquina, así que la
   regularización p+1e-4 y la forma exacta de `zerooffset` no se pueden comprobar aquí. El Methods
   publicado (derivadas centrales en x±ε/2) sí implica que la exposición central aporta 0 en el
   centro exacto, es decir, el valor puntual.

## Afirmaciones vivas previas
Las `refuted`/`unclear` antiguas del ledger (r01-r02, código) están fuera del alcance de esta
verificación. No las re-examiné.
