# Ronda 1: reporte del verificador

Todo se hizo de forma independiente. No importé `src/donutloc`, y no leí el código ni la nota de
derivación de los workers antes de tener mis propios resultados. Mis scripts están en
`work/verify/`:
- `vfisher.py`: Fisher con gradientes analíticos. Cubre la dona LG S17 y el límite cuadrático,
  SBR según S30, pedestal constante o gaussiano y el exponente c. Los gradientes coinciden con
  diferencias finitas centrales a 3e-10 (relativo).
- `check_analytic.py`: V1–V12.
- `vmc.py`, `vmc_bias.py`: Monte Carlo (V13–V17). Uso multinomial con semillas propias (20260926,
  1000+N, 777, 4242, 11–13, 99). El MLE es una grilla gruesa (paso σ_S27/8, ±10σ_S27) seguida de 4
  refinamientos ×10. Lo contrasté con Nelder-Mead desde 5 arranques: la grilla nunca quedó peor, y
  la diferencia máxima fue de 2e-5 nm.

## Derivaciones a mano (resumen)
Uso a = 4ln2/fwhm² y x = a L²/4.

**Valor puntual en r = 0.** En el origen los tres haces periféricos valen I = e·x·e^{-x}. Sus
gradientes son ∇I_i = −2 e a e^{-x} g r_i y suman 0. El haz central tiene gradiente 0. Resulta
∇p_i = −8 g r_i/(3L²). Como Σ r_i r_iᵀ = (3L²/8)·1, la Fisher es F = (8Ng²/L²)·1, que es la Eq. S27.

**Límite r → 0.** Con p₃ ≈ e a r²/S, el término central vale 4 e a/S · r̂r̂ᵀ = (16 eˣ/(3L²)) r̂r̂ᵀ,
que es β/N. Queda F_lim = α·1 + β r̂r̂ᵀ, y de ahí se obtienen V2 y V4 en forma cerrada. Además
√p₃ ∝ |r|, que no es diferenciable en 0 (V5).

**Fondo.** Con S30, p_out = (3+s)/12. Resulta σ² = S27²·(1+1/SBR)(1+3/(4SBR)), que es la Eq. S31.

**Pedestal constante.** En r = 0 se cumple ∇ΣI = 0, así que el peso de mezcla dependiente de r no
contribuye. El resultado es S31 con SBR_eps = ΣI/(4eps). Las dos reglas de combinación de V7 son
álgebra exacta.

**Pedestal gaussiano.** dI/dd² = a e^{-x}(e·g − eps). Con eso se obtiene exactamente el cociente
de V8.

**Multifotón (exponente c).** ∇p_i = c∇I_i/(3I), así que F escala como c². El término central va
como r^{2c−2}, que tiende a 0 para todo c > 1.

**LMS.** Σ_out p_i r_i = −g r a primer orden, de donde r̂ = −(1/g)Σ p_i r_i. En r = 0 su
covarianza es exactamente L²/(8Ng²)·1, así que σ_LMS(0) = S27 de forma exacta. El término de fondo
(1−s)/4·Σ r_i se anula, por lo que el LMS sin el factor 1/s es exactamente s veces el LMS
correcto. Esto vale para toda r, no solo en el régimen lineal.

## Resultados numéricos (los míos)
- **V1.** Obtengo 1.802472, igual a la fórmula hasta la precisión de máquina. Probé también L=100
  y (150, fwhm=360, N=1000). Con fwhm → ∞ da L/(2√(2N)).
- **V2.** El límite numérico en |r| = 1e-6, en 4 direcciones, coincide con la fórmula: 1.605096.
  El resultado es idéntico en todas las direcciones.
- **V3.** Obtengo ρ = 0.89439 / 0.89050 / 0.87805 / 0.85527, exactamente los valores de la
  afirmación. Fuera del caso cuadrático ρ ≠ 2/√5; por ejemplo, con L=150 y fwhm=200 da 0.798. En
  el caso cuadrático, ρ = 2/√5 y σ²·10N/L² = 1.0000000.
- **V4.** Obtengo σ_∥ = 2.73861 y σ_⊥ = 3.53553. El autovector del autovalor mayor es r̂.
- **V5.** El término central vale 0.2174807 a |r| = 1e-2, 1e-4 y 1e-6, frente a β/N = 0.2174807.
- **V6.**
  - Con SBR = 5 y 10 y L = 50 y 100, el valor puntual, el límite y S31 coinciden en 1e-15.
  - En r = r_c, la parte de señal de p₃ es igual a la parte de fondo, y el término central vale
    β/2. Lo comprobé con SBR = 1e4 y 1e6.
- **V7.** El pedestal da el valor S31 con SBR_eps hasta 1e-14 (8 casos). Las dos reglas de
  combinación de SBR coinciden en 1e-15.
- **V8.** El cociente numérico coincide con la fórmula en 1e-15; el límite r → 0 es isótropo y da
  el mismo valor. No equivale a un fondo: con la SBR que reproduce el mismo vector p en el centro
  (L=50, eps=0.05), el CRB es 3.8778, frente a 3.8785 con el pedestal gaussiano.
- **V9.** La tabla coincide dígito a dígito en los 16 valores.
- **V10.** El CRB vale 3.8767 en r = 0 y 3.6275 en r = (5, 0). La transición √(eps/(e a)) =
  4.887 nm no depende de L por construcción. Hay dos matices:
  - El CRB en r = 5 depende de la dirección: 3.628 sobre x, 3.694 sobre y y 3.583 a 45°.
  - La posición del mínimo sí depende de L. Sobre x queda en 3.5 / 6.7 / 9.0 nm para L = 50 / 100 /
    150, y el mínimo sobre x con L = 100 es 3.609 en r = 6.7. O sea, 4.9 nm es una escala, no el
    argmin.
- **V11.** Obtengo 3.5355 / 1.7678 / 1.1785. Para c = 2 y 3 el límite es igual al valor puntual,
  y también para c = 1.5.
- **V12.**
  - Con fwhm = 360 obtengo 0.941 / 1.962 / 3.167, y con 300, 0.947 / 2.012 / 3.370. El 1D da 1.25.
  - Los valores publicados (0.94 / 1.96 / 3.16) los tomé de `docs/literature/A_minflux_theory.md`;
    no abrí el PDF.
- **V13.**
  - Con mi semilla (8000 reps): σ = 1.9231 ± 0.0110, eficiencia 0.981 ± 0.006. Es compatible con
    0.988 ± 0.005.
  - Con la semilla 42 (20000 reps) y mi MLE obtengo 1.93614, lo que reproduce su 1.9362.
  - Nota: su error estándar σ/√(2R) sobreestima el error en √2. Con x e y agrupados, el correcto
    es ≈ σ/(2√R), es decir 0.0070 en lugar de 0.0097. Es conservador, no es un problema.
- **V14.**
  - Con mi semilla, σ/σ_lim vale 0.840 / 0.822 / 0.828 (±0.005) para N = 1e2 / 1e3 / 1e4, y
    σ/S27 = 0.748 con N = 100. Con la semilla 42 obtengo 1.34557, lo que reproduce su 1.3456.
  - El sesgo en (2, 0) tiene el signo correcto (hacia el centro). Su magnitud, con 24000 reps de
    semillas propias, es −0.338 ± 0.010 nm. Con la semilla 42 y 5000 reps, mi MLE da exactamente
    −0.290 ± 0.022. Es decir, −0.29 es una fluctuación de esa muestra (a 2σ) y el valor mejor
    estimado es ≈ −0.34 nm.
  - El n₃ fue siempre 0, como se espera.
- **V15.** Obtengo una eficiencia de 1.000 ± 0.007 (σ = 0.6264 frente a un CRB de 0.6267).
- **V16.** En MC, σ/S27 = 0.997 ± 0.004 y σ/σ_lim = 1.120. El valor analítico exacto es 1.000 y
  1/ρ = 1.1230, así que 0.998 y 1.121 son consistentes.
- **V17.** El cociente E[LMS sin 1/s] / E[LMS sin fondo] vale 0.90909 = 10/11 en (2,0), (0,3) y
  (4,−2). La contracción es exacta, no solo lineal. Por separado, el LMS sin fondo también contrae
  por no linealidad (0.96–0.99 a r = 2–4.5 nm), cosa que la afirmación no menciona.

## Pendiente
No queda nada refutado. Hay tres matices para el paper:
1. Reportar el sesgo como ≈ −0.34 ± 0.01 nm.
2. El 4.9 nm de V10 es una escala de transición, no la posición del mínimo, que depende de L y de
   la dirección.
3. Los errores estándar del MC están inflados por un factor √2.

## Bloque de afirmaciones (copiado por el orquestador desde la respuesta del verificador)

```claims
[{"status": "verified", "text": "V1 — Point value at r=0 without background (central p=0 term excluded) equals Eq. S27, sigma = L/(2 sqrt(2N)) / g; 1.802472 nm for L=50, N=100, fwhm=300 (hand derivation + analytic-gradient Fisher; fwhm->inf gives L/(2sqrt(2N)))."},
 {"status": "verified", "text": "V2 — r->0 limit without background: sigma_lim^2 = L^2/(8N g^2)(3g^2+e^x)/(3g^2+2e^x), sqrt(tr/2) independent of approach direction; 1.605096 nm for L=50, N=100, fwhm=300 (derived from F_lim = alpha I + beta r r^T; numeric at |r|=1e-6 in 4 directions agrees)."},
 {"status": "verified", "text": "V3 — rho = sigma_lim/sigma_pt = 2/sqrt5 exactly only in the quadratic case (sigma^2 = L^2/(10N)); at fwhm=300: 0.89439 (L=5), 0.89050 (L=50), 0.87805 (L=100), 0.85527 (L=150)."},
 {"status": "verified", "text": "V4 — The limiting error ellipse is not round: sigma_par = (alpha+beta)^-1/2 along r-hat, sigma_perp = sigma_S27, alpha = 8Ng^2/L^2, beta = 16N e^x/(3L^2); L=100, N=100, quadratic: 2.7386 and 3.5355 nm."},
 {"status": "verified", "text": "V5 — The discontinuity comes from the central Fisher term tending to (beta/N) r r^T, which has no matrix limit (checked: 0.2174807 at |r|=1e-2..1e-6 vs beta/N=0.2174807 for L=50); equivalently sqrt(p3) ∝ |r| is not differentiable at 0."},
 {"status": "verified", "text": "V6 — With finite SBR (5, 10) point value = limit = Eq. S31 (agreement 1e-15); the limits SBR->inf and r->0 do not commute; at r_c = (sqrt3 L/4) e^{-x/2}/sqrt(SBR) the signal and background parts of p3 are equal and the central term is beta/2 (checked at SBR=1e4, 1e6)."},
 {"status": "verified", "text": "V7 — Constant pedestal eps is exactly Eq. S31 at the centre with SBR_eps = 3e x e^{-x}/(4 eps) (because grad of sum I = 0 at r=0); combination rules 1/SBR_eff = 1/SBR + 1/SBR_eps (SBR against pure LG) and 1+1/SBR_eff = (1+1/SBR)(1+1/SBR_eps) (SBR against total beam signal) are exact (numerically 1e-15)."},
 {"status": "verified", "text": "V8 — Gaussian pedestal: F_eps/F_S27 = 3x^2(e g - eps)^2/[g^2(e x+eps)(3(e x+eps)+eps e^x)] (hand derivation, numerically 1e-15, point = limit); not equivalent to a background (at the SBR reproducing the same centre p-vector the CRB differs, e.g. 3.8778 vs 3.8785 at L=50, eps=0.05)."},
 {"status": "verified", "text": "V9 — Degradation CRB(eps)/CRB_S27, fwhm=300, constant|gaussian: L=50: 1.045|1.045, 1.227|1.228, 2.130|2.152, 4.382|4.584; L=100: 1.012|1.012, 1.060|1.061, 1.300|1.307, 1.898|1.958 for eps = 0.002/0.01/0.05/0.15 (all 16 values reproduced)."},
 {"status": "verified", "text": "V10 — With eps=0.002, L=100, fwhm=300, N=100 the CRB minimum is not at the centre: 3.877 nm at r=0 vs 3.628 nm at r=(5,0); transition scale sqrt(eps/(e a)) = 4.887 nm, L-independent by construction. Caveat: the CRB at r=5 depends on direction (3.628 along x, 3.694 along y, 3.583 at 45 deg), and the argmin position depends on L (3.5/6.7/9.0 nm along x for L=50/100/150), so 4.9 nm is a scale, not the minimum location."},
 {"status": "verified", "text": "V11 — Multiphoton exponent c: sigma_pt = L/(2c sqrt(2N)) / g; for c>=2 (in fact any c>1) limit = point value (central term ~ r^{2c-2}); L=100, N=100, quadratic: 3.536/1.768/1.179 nm for c=1/2/3."},
 {"status": "verified", "text": "V12 — Eq. S31 with N=500, SBR=5 at the centre gives 0.941/1.962/3.167 nm with fwhm=360 and 0.947/2.012/3.370 with fwhm=300 (published 0.94/1.96/3.16, taken from docs/literature note, PDF not opened); 1D L/(4sqrtN) = 1.25 nm for L=50, N=100."},
 {"status": "verified", "text": "V13 — Centre MLE with SBR=10 (L=50, N=100, fwhm=300) is efficient: own seed, 8000 reps: sigma = 1.9231 ± 0.0110, efficiency 0.981 ± 0.006, consistent with the claimed 1.9362 ± 0.0097 / 0.988 ± 0.005 (seed 42 reproduced exactly: 1.93614); CRB 1.96006 = Eq. S31. Note: the claimed SE sigma/sqrt(2R) is conservative by sqrt2."},
 {"status": "verified", "text": "V14 — Centre MLE without background is superefficient: sigma/sigma_lim = 0.840/0.822/0.828 (±0.005) for N=1e2/1e3/1e4 (own seeds), sigma/S27 = 0.748 at N=100; seed 42 reproduces sigma = 1.3456. Biased toward the centre: at r=(2,0), N=100 the bias in x is -0.338 ± 0.010 nm (24000 reps); the quoted -0.29 is a 2-sigma fluctuation of the 5000-rep seed-42 sample (reproduced exactly as -0.290 ± 0.022), so the better value to report is ≈ -0.34 nm."},
 {"status": "verified", "text": "V15 — MLE without background at r=(10,0), N=1000: efficiency 1.000 ± 0.007 (sigma 0.6264 vs CRB 0.6267), consistent with 0.999."},
 {"status": "verified", "text": "V16 — LMS without background at the centre: sigma_LMS = S27 exactly (analytic covariance L^2/(8N g^2)), so sigma/S27 = 1 and sigma/sigma_lim = 1/rho = 1.123; MC 0.997 ± 0.004 and 1.120; claimed 0.998/1.121 consistent."},
 {"status": "verified", "text": "V17 — Without the 1/s factor the LMS with SBR=10 shrinks the estimate by exactly s = 10/11 relative to the background-free LMS, at any r, not only linearly (the constant term (1-s)/4 * sum r_i vanishes); checked 0.90909 at (2,0), (0,3), (4,-2)."}]
```
