# Ronda 2 — Verificador (independiente)

Todo el código es mío y está en `work/verify/r02/`. No importé `src/donutloc`. Reutilicé solo mi `work/verify/vfisher.py` de la ronda 1 (Fisher multinomial con gradientes analíticos, LG, TCP). Leí los reportes de los workers **después** de tener mis números, y solo para diagnosticar W11.

## Scripts
- `rw.py`: campo focal de Richards-Wolf por **cuadratura 2D directa** de la pupila, como suma de ondas planas. No hay reducción a Bessel.
  - Cuadratura: θ por Gauss-Legendre en [0, θmax] con sinθmax = NA/n; φ uniforme, que es espectralmente exacto.
  - Pupila: amplitud exp(−(sinθ/(F sinθmax))²) (gaussiana de amplitud con ρ' = f sinθ, h = f sinθmax), por √cosθ, por e^{iφ}.
  - Polarización: E_inc se descompone en ρ̂ → θ̂ y φ̂ → φ̂.
  - Gradiente exacto (factor ik), así que la curvatura y dI/dρ salen sin diferencias finitas.
- `w1_5.py` y `w2b.py`: W1–W4. `w5.py`: W5. `w6.py`: W6. `cam.py`: W7–W9. `w12.py`: W12. `w10.py` y `w11.py`: W10–W11. `w13.py` y `w13b.py`: W13.

Precisión: los resultados vectoriales son idénticos en todos los dígitos mostrados con grillas (150–400) × (64–160). La simetría de rotación da una dispersión de 1e-15.

## Resultados

**W1 — verificada.** Con ℓ = 1 y σ = +1, I(0)/Imax vale 1e-32 en z = 0, 300 y −500 nm (cero de máquina). Tiene que ser así: con σ = +1 todas las componentes llevan e^{iφ}, e^{2iφ} o e^{3iφ}, que se integran a 0 en φ para todo z.

**W2 — verificada.** Profundidad I(0)/Imax, con Imax buscado en 2D:

| caso | n = 1.518 | n = 1.5 |
|---|---|---|
| mano opuesta | 0.845241 | 0.868963 |
| lineal x | 0.371648 | 0.382305 |

Con lineal y obtuve el mismo valor que con lineal x.

**W3 — verificada.** D_pp = 384.666 nm (n = 1.518, F = 5/3). Con n = 1.5: 515.52, 401.35, 382.87 y 377.62 nm para F = 0.5, 1, 2 y 100. F = 1e4 da 377.61.
- Observación: [Caprile2022, Fig. 10], según la nota de literatura, da 508, 428, 392 y 380 nm. En F = 1 la diferencia es de 27 nm. Mis números coinciden con los del worker por un camino independiente, así que es probable que la diferencia venga de cómo se leyó la figura o de una convención de F distinta en PyFocus. Conviene no afirmar concordancia con Caprile en F = 1.

**W4 — verificada.**
- c = |∂E/∂x|²/Imax = 7.04223e-5 nm⁻². Igual en y. Coincide con I/(Imax ρ²) en ρ = 0.1 nm.
- fwhm por curvatura = 327.141 nm; fwhm por diámetro = 320.256 nm.

**W5 — verificada.**
- Con σ = −1, E_x(0) y E_y(0) son ~1e-17 y E_z(0) ≠ 0. E_z(0, z) coincide con −(2π/√2)∫ g sinθ e^{ikz cosθ} dθ (g = a√cosθ sinθ) a 1e-15 en z = 0, 200 y 500 nm. Fuera del eje ese núcleo es J0(kρ sinθ).
- Fracción de energía en |E_z|² en el plano focal con σ = +1: 21.76 % por Parseval, que es exacto. La integral directa en el plano hasta 3 µm da 21.36 %, porque faltan las colas. Otros llenados: F = 1 da 17.1 % y F = 100 da 24.4 %.

**W6 — verificada.** Fisher multinomial directo con el patrón vectorial real (gradiente exacto), evaluado en |r| = 1e-3 nm en 3 direcciones. También lo comprobé con la fórmula general del límite F/N = (f'L)²/(2f²)·Id + (4c/3f) r̂r̂, con f' = df/d(d²).

| haz | L = 50 | L = 100 | L = 150 |
|---|---|---|---|
| vectorial correcta | 1.60041 | 3.32324 | 5.33104 |
| vectorial opuesta | 147.886 | 89.588 | 84.197 |
| lineal (patrón 2D real, polarización x fija) | 60.774 | 38.061 | 35.723 |
| LG fwhm = 327.14 | 1.60124 | 3.32964 | 5.35145 |
| LG fwhm = 300 | 1.60510 | 3.36341 | 5.48647 |

Diferencias de la LG igualada por curvatura respecto de la vectorial: 0.05, 0.19 y 0.38 %.

**W7 — verificada.**
- 9×9, píxel 100 nm, N = 400: 5.20455 nm. Sin normalizar a la ventana da 5.20458; con el emisor en la esquina del píxel, 5.2057.
- Píxel 10 nm, ventana 121×121: CRB/(σ/√N) = 1.000417. Con píxel 1 nm da 1.000004.

**W8 — verificada con la convención explícita.** Uso SBR_c = señal total / fondo por píxel, es decir SBR_total = SBR_c/81 = 6.17.
- N = 600 da 4.9642 nm; para 5 nm hacen falta 591.4 fotones.
- Si SBR = 500 se interpreta como total, el CRB es 4.268 nm y basta con 437 fotones. El texto tiene que fijar la convención.

**W9 — verificada.** Eq. S31 en el valor puntual r = 0 da 13.00, 13.45, 14.16, 15.37 y 17.93 fotones para SBR ∞, 50, 20, 10 y 5.
- Matiz: sin fondo el número corresponde al valor puntual (Eq. S27). El límite r → 0 da 10.3 fotones.

**W10 — verificada.** Implementé mi propio lazo: L geométrico 150 → 25, recentrado, MLE por grilla con refinamiento dentro del disco 0.75 L_k, 500 repeticiones, semillas mías.

| caso | σ por iteración (nm) | σ final (nm) |
|---|---|---|
| sin fondo, 2 semillas | 4.19/1.86/0.93/0.51 | 0.506 ± 0.012 y 0.521 ± 0.012 |
| SBR = 10 | — | 0.606 ± 0.014 y 0.607 ± 0.013 |
| sin recentrado, semilla A | 4.34/3.16/3.17/7.94 | — |
| sin recentrado, semilla B | 4.11/3.44/3.20/8.53 | — |

Todo es compatible con lo afirmado. La cámara da 3.162 nm (trivial).

**W11 — unclear.** Lo cualitativo se reproduce: pendiente cercana a −1/2 y un cociente con la cámara de ~0.16. Los números precisos no. Con estadística alta (3000–6000 repeticiones; 500–1000 en N = 1000) obtengo:

| N | σ (nm) | cociente con la cámara |
|---|---|---|
| 250 | 1.180 ± 0.020 | 0.187 |
| 500 | 0.748 | 0.167 |
| 1000 | 0.509 | 0.161 |
| 2000 | 0.355 | 0.159 |
| 4000 | 0.2514 ± 0.0024 | 0.159 |
| 8000 | 0.1736 | 0.155 |

- La pendiente log-log es −0.545 en todo el rango y −0.523 si se excluye N = 250.
- El worker obtuvo 1.094 en N = 250 (500 repeticiones). Ninguno de mis 12 bloques independientes de 500 repeticiones baja de 1.103: el rango es 1.103–1.256.
- En N = 4000 el worker da 0.233 frente a mi 0.251 ± 0.002.
- La distribución en N = 250 tiene colas pesadas: el 0.4 % de los errores supera 5 nm, y en una semilla con 1500 repeticiones un solo valor atípico de 72 nm llevó σ a 1.75. Con 500 repeticiones, la pendiente y el borde superior del cociente no son robustos.
- Para decidir: repetir con muchas más repeticiones y el mismo protocolo, o informar "pendiente ≈ −0.52 a −0.55 según el rango; cociente 0.155–0.19".

**W12 — verificada, con matices.** El pedestal es gaussiano, I = LG + eps·e^{−a r²}, con N = 100.

| eps | L_opt (nm) | CRB(L_opt) (nm) |
|---|---|---|
| 0.002 | 10.484 | 0.7458 |
| 0.01 | 23.280 | 1.6787 |
| 0.05 | 50.337 | 3.8784 |
| 0.15 | 80.810 | 7.2867 |

- Con eps = 0, el CRB (tanto el valor puntual como el límite) crece monótonamente en L **solo para L menor que el diámetro del anillo**, fwhm/√ln2 = 360.4 nm. Ahí diverge (f' = 0) y más allá decrece: en L = 2000 nm vale 2.37 nm. Es una rama sin sentido físico sin fondo. La afirmación debe acotarse a L < D_pp.
- La ley empírica L_opt/(fwhm√eps) da 0.781, 0.776, 0.750 y 0.695. No depende de fwhm (lo comprobé con 200 y 400 nm). El valor 0.78 vale para eps ≲ 0.01; para eps = 0.15 es 0.70.
- Con pedestal constante los números cambian poco: L_opt = 10.48, 23.24, 49.94 y 79.38 nm.

**W13 — verificada, con matiz estadístico.** Implementé mi propio MLE con SBR = 10, L = 100, N = 500, 20 patrones × 200 repeticiones, desplazando los 4 ceros, semilla mía.

| δ (nm) | sesgo ingenuo (nm) | sesgo honesto (nm) | σ ingenuo (nm) |
|---|---|---|---|
| 0 | 0.19 | 0.19 | 1.84 |
| 2 | 1.45 | 0.16 | 1.84 |
| 5 | 3.81 | 0.20 | 1.87 |
| 10 | 7.76 | 0.16 | 2.02 |

- Calculé el sesgo asintótico con conteos esperados: la media sobre 20 patrones tiene una dispersión entre sorteos de patrones de ±10 %. Para δ = 2, 5 y 10 da 1.49 ± 0.15, 3.78 ± 0.35 y 7.84 ± 0.72.
- Los valores afirmados (1.58, 4.01, 8.31) caen dentro de esa dispersión. Son un sorteo particular, no la media poblacional (≈ 0.75 δ).
- El piso del honesto (≈ 0.17) y "la pérdida es sesgo" se confirman.

## Sin resolver
- W11: diferencias estadísticas en N = 250 y N = 4000 (ver arriba).
- Discrepancia de W3 con la figura de Caprile en F = 1 (428 frente a 401 nm): la nuestra es reproducible; no es un error de W3.

```claims
[{"status": "verified", "text": "W1 — Vectorial RW donut, l=1 with sigma=+1 circular: on-axis intensity exactly 0 (independent direct 2D pupil quadrature gives I(0)/Imax ~1e-32 at z=0, 300, -500 nm; all components carry e^{i phi}, e^{2i phi}, e^{3i phi})."},
 {"status": "verified", "text": "W2 — Zero depth I(0)/Imax (Imax searched in 2D): opposite handedness 0.845241 (n=1.518) / 0.868963 (n=1.5); linear 0.371648 (n=1.518) / 0.382305 (n=1.5); x and y linear identical."},
 {"status": "verified", "text": "W3 — Peak-to-peak diameter 384.666 nm (n=1.518, F=5/3); with n=1.5: 515.52/401.35/382.87/377.62 nm for F=0.5/1/2/100 (independent 2D quadrature). Caveat: differs from Caprile2022 Fig.10 as transcribed (508/428/392/380), notably at F=1; do not claim agreement with Caprile at F=1."},
 {"status": "verified", "text": "W4 — Curvature at the zero c = 7.04223e-5 nm^-2 (n=1.518; exact via |dE/dx|^2/Imax, isotropic); curvature-matched LG fwhm = 327.141 nm; diameter-matched fwhm = 320.256 nm."},
 {"status": "verified", "text": "W5 — sigma=-1: on-axis Ex,Ey = 0 and Ez(0,z) = -(2pi/sqrt2) int g sin(theta) e^{ikz cos theta} dtheta != 0 (J0 kernel; matched 1e-15 at z=0,200,500 nm); sigma=+1: |Ez|^2 carries 21.76% of focal-plane energy (Parseval, exact; direct plane integral to 3 um 21.36%)."},
 {"status": "verified", "text": "W6 — r->0 centre CRB, N=100, no background, n=1.518, L=50/100/150: vectorial correct 1.60041/3.32324/5.33104; LG fwhm=327.14 1.60124/3.32964/5.35145; LG fwhm=300 1.60510/3.36341/5.48647; opposite 147.886/89.588/84.197; linear (true 2D pattern, fixed x pol) 60.774/38.061/35.723; curvature-matched LG within 0.05/0.19/0.38% (direct Fisher with exact field gradient at |r|=1e-3 in 3 directions + general limit formula F/N=(f'L)^2/(2f^2) Id + 4c/(3f) rr)."},
 {"status": "verified", "text": "W7 — Pixelated Gaussian PSF Fisher: sigma_PSF=100, pixel 100, 9x9, no bkg, N=400: CRB=5.2045 nm; pixel 10 nm, 121x121: CRB/(sigma/sqrtN)=1.000417."},
 {"status": "verified", "text": "W8 — With SBR_c=500 defined as total signal / background per pixel (SBR_total = 500/81 = 6.17), 9x9, N=600: CRB=4.964 nm and 591 photons for 5 nm. Caveat: if SBR=500 were total, CRB=4.268 nm and 437 photons; the convention must be stated."},
 {"status": "verified", "text": "W9 — MINFLUX S31 at the centre (point value), L=50, fwhm=300: 13.00/13.45/14.16/15.37/17.93 photons for 5 nm at SBR inf/50/20/10/5 (at SBR=inf this is the Eq. S27 point value; the r->0 limit would give 10.3)."},
 {"status": "verified", "text": "W10 — Iterative MINFLUX (N_total=1000, 4x250, L=150->82.55->45.43->25, recentring, MLE in disk 0.75 L_k, 500 reps, own seeds): final sigma 0.506±0.012 / 0.521±0.012 nm without bkg, 0.606 / 0.607±0.014 with SBR=10, consistent with 0.497±0.011 / 0.613±0.014; without recentring 4.34/3.16/3.17/7.94 and 4.11/3.44/3.20/8.53 nm, consistent with 4.34/3.20/3.23/8.06; camera 3.162 nm."},
 {"status": "unclear", "text": "W11 — sigma_iter vs N_total slope -0.526 and camera ratio 0.15–0.17 — qualitatively reproduced, not at the stated precision: with 3000–6000 reps I get sigma = 1.180±0.020 (N=250, ratio 0.187), 0.748, 0.509, 0.355, 0.2514±0.0024 (N=4000; worker 0.233), 0.1736 nm, slope -0.545 over the full range (-0.523 for N>=500), ratio 0.155–0.187. N=250 is heavy-tailed (none of 12 blocks of 500 reps reaches 1.094). Needs more reps under the same protocol, or restating as slope ≈ -0.52…-0.55 depending on range, ratio 0.155–0.19."},
 {"status": "verified", "text": "W12 — Gaussian pedestal eps, N=100, fwhm=300, centre CRB: L_opt = 10.484/23.280/50.337/80.810 nm, CRB = 0.7458/1.6787/3.8784/7.2867 nm for eps=0.002/0.01/0.05/0.15. Caveats: with eps=0 the CRB is monotonic increasing only for L < ring diameter fwhm/sqrt(ln2)=360 nm (diverges there, then decreases); L_opt/(fwhm sqrt eps) = 0.781/0.776/0.750/0.695 (fwhm-independent), so 0.78 holds only for eps <~ 0.01."},
 {"status": "verified", "text": "W13 — Misalignment (L=100, N=500, SBR=10, all 4 zeros displaced by delta in random directions, 20 patterns x 200 reps): naive MLE centre bias 1.45/3.81/7.76 nm (own draw) vs honest 0.16–0.20 (MC floor ≈0.17); naive sigma 1.84->2.02: the loss is bias. The claimed 1.58/4.01/8.31 are inside the ±10% pattern-sampling spread (population means from noise-free MLE 1.49±0.15 / 3.78±0.35 / 7.84±0.72 over 20-pattern draws, ≈0.75 delta)."}]
```
