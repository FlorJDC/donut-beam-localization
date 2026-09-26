# Ronda 1 — Worker C: CRB en el centro del TCP (derivación, formas cerradas, verificación)

## Qué hice
- `docs/derivations/crb_tcp_center.md`: derivación completa en español.
  - Fórmula maestra (M) para cualquier haz radial.
  - Secciones (a) punto S27, (b) límite r→0, (c) fondo S31 y no conmutatividad, (d) eps en los
    modelos constant y gaussian, (e) multifotón.
  - 1D y S32, perfil radial y tabla de degradación.
- `src/donutloc/closed_forms.py` (solo numpy, vectorizado):
  - `crb_tcp_center_point(L,N,fwhm=inf,sbr=inf,power=1)`;
  - `crb_tcp_center_limit(L,N,fwhm=inf,power=1)`;
  - `crb_tcp_center_limit_axes`;
  - `limit_to_point_ratio`;
  - `crb_tcp_center_eps(L,N,fwhm,eps,sbr=inf,zero_model="constant",sbr_ref="beam"|"total"|"lg")`;
  - `sbr_eff_constant_pedestal`;
  - `crossover_radius`;
  - `sbr_center_vs_L`;
  - `crb_1d_center(kind=donut|quadratic|gaussian)`;
  - `x_param`.
- `tests/test_closed_forms.py`: 21 tests con Fisher numérica inline propia. Cubren punto, límite
  (12 direcciones, r0 = 1e-3), SBR 5/10, eps 0.01/0.1 en los dos modelos y las dos convenciones,
  L 50/100, fwhm 300/360/∞, los checks de la nota A §8.2, los límites no conmutativos, el
  multifotón, errores y vectorización.
- `scripts/verify_crb_closed_forms.py`: 89 casos exactos más 6 filas "approx" informativas. Tiene
  una columna del paquete, y `donutloc.fisher/photons/patterns/beams` importaron: los 89 casos
  coinciden con el paquete. Termina con exit 1 si el error supera 1e-6 (punto) o 1e-3 (límite).
- No toqué archivos de otros workers. `src/donutloc/__init__.py` ya existía (Worker A).

## Salidas reales
- `python -m unittest discover -s tests -p "test_closed_forms.py"` → `Ran 21 tests ... OK`.
- `python scripts/verify_crb_closed_forms.py` → `all 89 exact checks pass (max rel. error point
  2.13e-07, limit 2.23e-07)`, EXIT=0. El máximo sale de L=5; con L ≥ 50 el error es ≤ 2.3e-9.
- Suite completa `python -m unittest discover -s tests` → 84 tests OK (1 skip, el de aceptación).

## Resultados (fórmulas)
Notación: x = L² ln2/fwhm², g = 1−x, R = L/2.
- (M) En el centro del TCP, para cualquier haz radial h(u) con u = r²:
  F = 6N R² h'(R²)² / [h(R²)(3h(R²)+h(0))] · 1.
  - Fondo b por exposición: h(R²) → h(R²)+b y 3h(R²)+h(0) → 3h(R²)+h(0)+4b.
- (a) σ_pt = L/(2√(2N))·g⁻¹. Es Eq. S27.
- (b) σ_lim² = L²/(8N g²) · (3g²+eˣ)/(3g²+2eˣ).
  - ρ = σ_lim/σ_pt = √((3g²+eˣ)/(3g²+2eˣ)).
  - Cuadrático: ρ = 2/√5 y σ² = L²/(10N).
  - Para x ≪ 1: ρ ≈ (2/√5)(1−9x/40).
  - El valor de sqrt(tr/2) es el mismo para toda dirección de aproximación. La elipse no es
    isótropa: σ_∥ = (α+β)^−1/2, σ_⊥ = σ_pt, con α = 8Ng²/L² y β = 16N eˣ/(3L²).
- (c) Con fondo: S31. Es continuo porque el término central es O(r²).
  - Radio de transición: r_c = (√3 L/4) e^{−x/2}/√SBR (orden dominante).
  - Límites: lim_{SBR→∞} lim_{r→0} = S27, y lim_{r→0} lim_{SBR→∞} = σ_lim.
- (d) eps, modelo constante (exacto): S31 con SBR_eps = 3e x e^{−x}/(4 eps).
  - "lg": 1/SBR_eff = 1/SBR + 1/SBR_eps.
  - "beam"/"total" (convención de photons.probabilities): 1+1/SBR_eff = (1+1/SBR)(1+1/SBR_eps).
    Coincide con la fórmula que relayó Worker A.
- (d) eps, modelo gaussiano (exacto, NO es un fondo efectivo):
  F_eps/F_pt = 3x²(e g − eps)² / [g²(e x + eps)(3(e x+eps) + eps eˣ)].
  - A primer orden en eps coincide con el modelo constante, con error relativo 3x²/7 en el
    coeficiente.
- (e) σ^(c)_pt = L/(2c√(2N))·g⁻¹. Para c ≥ 2 el límite coincide con el punto: no hay
  discontinuidad.

## Número de aceptación
crb_center_lg_L50_N100_nm (fwhm = 300, límite r→0) = **1.605096 nm**. El valor puntual S27 es
1.802472 nm.

## Hallazgos no previstos
- Con un cero casi perfecto (eps = 0.002 o SBR alto), **el mínimo del CRB no está en el centro**:
  está en un anillo de radio ~ r_c.
  - Perfil numérico con L = 100, eps = 0.002: 3.877 nm en r = 0 y 3.628 nm en r = 5.
  - El radio de transición por eps es r_c = √(eps/(e a)) = 4.9 nm (fwhm = 300, eps = 0.002) y no
    depende de L.
- La fórmula de transición es solo de orden dominante:
  - error 1e-4 con r_c = 0.4 nm;
  - 1.6 % con r_c = 4.9 nm;
  - falla (12 % en r = r_c) con SBR = 10 (r_c = 13 nm, L = 100).

  Ahí se reporta como "approx" y no como afirmación exacta.
- Con eps > 0, cuando eps → 0 el CRB en el centro tiende a **S27**, no al límite (b). Es la misma
  no conmutatividad. La tabla de degradación da los dos cocientes.

## No resuelto / pendiente
- La transición cerca del centro para r_c ~ L/10 no tiene forma cerrada: haría falta llevar los
  términos O(r/L) de las exposiciones periféricas.
- La sección 1D está verificada numéricamente (S21); allí solo cito S22c/S22f/S23c y no las
  re-derivo en detalle.

```claims
[
{"text": "Fórmula maestra: en el centro del TCP (L=diámetro), para cualquier haz radial h(u=r^2), F = 6 N R^2 h'(R^2)^2 / [h(R^2)(3h(R^2)+h(0))] * Id (término central nulo o excluido); fondo b por exposición: h(R^2)->h(R^2)+b, denominador +4b."},
{"text": "(a) Valor puntual en r=0 sin fondo (término central p=0 excluido) = L/(2 sqrt(2N)) (1 - L^2 ln2/fwhm^2)^-1 = Eq. S27; verificado contra Fisher numérica a <=2.1e-9 (L>=50) para fwhm 300/360/inf."},
{"text": "(b) Límite r->0 sin fondo: sigma_lim^2 = L^2/(8 N g^2) (3g^2+e^x)/(3g^2+2e^x), x=L^2 ln2/fwhm^2, g=1-x; independiente de la dirección de aproximación; verificado a <=2.3e-9 (promedio 12 direcciones, r0=1e-3)."},
{"text": "Hipótesis 2/sqrt(5) CONFIRMADA solo en el límite cuadrático (sigma_lim^2 = L^2/(10N)); con fwhm finito rho = sqrt((3g^2+e^x)/(3g^2+2e^x)) ≈ (2/sqrt5)(1-9x/40): 0.89439 (L=5), 0.89050 (L=50), 0.87805 (L=100), 0.85527 (L=150) con fwhm=300."},
{"text": "Covarianza límite anisótropa: sigma_par = (alpha+beta)^-1/2, sigma_perp = sigma_S27, alpha=8Ng^2/L^2, beta=16N e^x/(3L^2); L=100, N=100 cuadrático: 2.7386 / 3.5355 nm."},
{"text": "crb_center_lg_L50_N100_nm (fwhm=300, límite r->0) = 1.605096 nm; el valor puntual S27 es 1.802472 nm."},
{"text": "(c) Con SBR fija (Eq. S30) el CRB en el centro es Eq. S31 y es continuo (punto = límite a <=1.3e-9 para SBR 5 y 10); lim_{SBR->inf} de S31 da S27, no el límite (b): los límites no conmutan; el término central se enciende en r_c = (sqrt3 L/4) e^{-x/2}/sqrt(SBR) (orden dominante)."},
{"text": "(d) Pedestal constante eps: exacto = S31 con SBR_eps = 3 e x e^{-x}/(4 eps); con fondo adicional: 'lg' 1/SBR_eff = 1/SBR + 1/SBR_eps, 'beam'/'total' (convención de photons) 1+1/SBR_eff = (1+1/SBR)(1+1/SBR_eps); verificado a <=2e-9 (numérico propio y paquete)."},
{"text": "(d) Pedestal gaussiano eps (no equivale a un fondo): F_eps/F_S27 = 3x^2 (e g - eps)^2 / [g^2 (e x+eps)(3(e x+eps)+eps e^x)], exacto; verificado a <=2e-9; a primer orden en eps difiere del constante en un factor 1+3x^2/7 del coeficiente."},
{"text": "Degradación CRB(eps)/CRB_S27 (sin fondo, constante|gaussiano), fwhm=300: L=50: eps=0.002 1.045|1.045, 0.01 1.227|1.228, 0.05 2.130|2.152, 0.15 4.382|4.584; L=100: 1.012, 1.060|1.061, 1.300|1.307, 1.898|1.958; tabla completa (fwhm 300/360, L 50/100/150) en docs/derivations/crb_tcp_center.md."},
{"text": "(e) Multifotón exponente c: sigma_pt = L/(2c sqrt(2N)) g^-1, y para c>=2 límite = punto (sin discontinuidad); L=100,N=100: 3.536/1.768/1.179 nm (cuadrático), 3.831/1.915/1.277 nm (fwhm=300)."},
{"text": "Checks nota A §8.2 reproducidos: 3.536, 3.735, 1.792, 1.817, Masullo 0.941/1.962/3.167 (fwhm 360) vs 0.947/2.012/3.370 (fwhm 300), 10.82, 1.25, S32 2.602/0.657, 2.739."},
{"text": "Con cero casi perfecto (eps=0.002, L=100, fwhm=300) el CRB mínimo no está en el centro: 3.877 nm en r=0 vs 3.628 nm en r=5 nm (numérico, 12 direcciones)."}
]
```
