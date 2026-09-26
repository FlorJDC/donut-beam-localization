# Ronda 1 — afirmaciones a verificar (solo las afirmaciones; sin el razonamiento de los workers)

Convenciones: dona LG Eq. S17 con parámetro fwhm; TCP = 3 donas a ángulos π/2+2πk/3 en un círculo
de diámetro L + centro (índice 3); p_i = I_i/ΣI (sin fondo), Eq. S30 con SBR; Fisher
F = N Σ ∇p_i∇p_iᵀ/p_i; σ_CRB = sqrt(tr F⁻¹ / 2). Notación x = L² ln2/fwhm², g = 1−x.

## Física / matemática
V1. Valor puntual en r=0 sin fondo (excluyendo el término central con p=0) = Eq. S27:
    σ_pt = L/(2√(2N))·g⁻¹. Para L=50, N=100, fwhm=300: 1.802472 nm.
V2. Límite r→0 sin fondo: σ_lim² = L²/(8N g²)·(3g²+eˣ)/(3g²+2eˣ); isótropo en sqrt(tr/2) respecto
    de la dirección de aproximación. Para L=50, N=100, fwhm=300: 1.605096 nm.
V3. Razón ρ = σ_lim/σ_pt = 2/√5 exactamente solo para cero cuadrático (σ² = L²/(10N)); a fwhm=300:
    0.89439 (L=5), 0.89050 (L=50), 0.87805 (L=100), 0.85527 (L=150).
V4. La elipse de error en el límite no es redonda: σ_∥ = (α+β)^−½, σ_⊥ = σ_S27, con α = 8Ng²/L²,
    β = 16N eˣ/(3L²). L=100, N=100, cuadrático: σ_∥ = 2.7386, σ_⊥ = 3.5355 nm.
V5. La discontinuidad en r=0 se debe a que el término central tiende a (β/N) r̂r̂ᵀ, que no tiene
    límite matricial; equivalentemente √p₃ ∝ |r| no es diferenciable en 0 (falla la regularidad).
V6. Con SBR finita (5, 10) valor puntual = límite = Eq. S31 (sin discontinuidad). Los límites
    SBR→∞ y r→0 no conmutan; radio de transición (orden dominante) r_c = (√3 L/4)e^{−x/2}/√SBR.
V7. Pedestal constante eps (I = I_LG + eps): exactamente Eq. S31 con SBR_eps = 3e·x·e^{−x}/(4 eps).
    Con SBR medida contra la señal LG pura: 1/SBR_eff = 1/SBR + 1/SBR_eps. Con SBR medida contra
    la señal total del haz: 1+1/SBR_eff = (1+1/SBR)(1+1/SBR_eps).
V8. Pedestal gaussiano eps·exp(−4ln2 r²/fwhm²): F_eps/F_S27 = 3x²(e·g−eps)² / [g²(e·x+eps)(3(e·x+eps)+eps·eˣ)];
    no equivale a un fondo.
V9. Degradación CRB(eps)/CRB_S27, fwhm=300, constante|gaussiano: L=50: eps 0.002→1.045|1.045,
    0.01→1.227|1.228, 0.05→2.130|2.152, 0.15→4.382|4.584; L=100: 1.012|1.012, 1.060|1.061,
    1.300|1.307, 1.898|1.958.
V10. Con eps=0.002, L=100, fwhm=300 (N=100) el CRB mínimo NO está en el centro: 3.877 nm en r=0 vs
    3.628 nm en r=5 nm; radio de transición √(eps/(e·a)) ≈ 4.9 nm, independiente de L.
V11. Multifotón (exponente c): σ_pt = L/(2c√(2N)) g⁻¹; para c ≥ 2 límite = valor puntual.
    L=100, N=100, cuadrático: 3.536/1.768/1.179 nm para c=1/2/3.
V12. Reproducción de la literatura: tabla de Masullo (N=500, SBR=5, centro) 0.941/1.962/3.167 nm con
    fwhm=360 (0.947/2.012/3.370 con 300) vs publicados 0.94/1.96/3.16; 1D L/(4√N) = 1.25 nm (L=50,
    N=100).

## Estadística (Monte Carlo, semilla 42, L=50, N=100, fwhm=300)
V13. MLE en el centro con SBR=10 (20000 reps): σ = 1.9362 ± 0.0097 nm, eficiencia σ/CRB = 0.988 ± 0.005
    (CRB = 1.96006 nm = Eq. S31).
V14. MLE en el centro SIN fondo es superficiente: σ = 1.3456 ± 0.0067 nm, σ/σ_lim = 0.838 ± 0.004,
    σ/S27 = 0.747; independiente de N (0.837/0.831/0.830 para N=1e2/1e3/1e4); sesgado hacia el
    centro (sesgo −0.29 nm en x para r=(2,0), N=100).
V15. MLE sin fondo en r=(10,0), N=1000: eficiencia 0.999.
V16. LMS sin fondo en el centro: σ/σ_lim = 1.121, σ/S27 = 0.998.
V17. Sin el factor 1/s el LMS con SBR=10 contrae la estimación por s = 10/11.
