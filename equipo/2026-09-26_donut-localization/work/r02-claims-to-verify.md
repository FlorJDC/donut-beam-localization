# Ronda 2 — afirmaciones a verificar (solo afirmaciones; sin el razonamiento de los workers)

Convenciones de la Ronda 1 (ya verificadas): LG Eq. S17 con fwhm; TCP 3+centro; Fisher multinomial;
σ = sqrt(tr F⁻¹/2); límite r→0 en ceros perfectos. Ópticos: λ = 640 nm, NA = 1.4, n = 1.518
(salvo que se diga n = 1.5), pupila gaussiana con filling F = w0/h = 5/3, apodización aplanática.

## Dona vectorial (Richards-Wolf, vórtice de carga 1)
W1. Con polarización circular de handedness correcta (lσ = +1) la intensidad en el eje es exactamente
    0 (en z = 0 y z = 300 nm).
W2. Profundidad del cero I(0)/Imax: handedness opuesta 0.84524 (n = 1.518) / 0.86896 (n = 1.5);
    polarización lineal 0.37165 (n = 1.518) / 0.38231 (n = 1.5).
W3. Diámetro pico a pico (handedness correcta): 384.67 nm (n = 1.518). Con n = 1.5: 515.5 (F = 0.5),
    401.3 (F = 1), 382.9 (F = 2), 377.6 nm (F = 100).
W4. Curvatura en el cero (límite ρ→0 de I/Imax = cρ²): c = 7.0422e-5 nm⁻² (n = 1.518); fwhm de la LG
    que iguala la curvatura (4e ln2/fwhm² = c) = 327.14 nm; la que iguala el diámetro (fwhm/√ln2 = D_pp)
    = 320.26 nm.
W5. Para σ = −1 (l = 1): E_z en el eje ∝ ∫g sinθ J0 ≠ 0 (de ahí el cero lleno); para σ = +1 la
    componente z lleva ~22 % del flujo total en el plano focal.
W6. CRB en el centro del TCP (límite r→0, N = 100, sin fondo, n = 1.518):
    L = 50/100/150: vectorial correcta 1.60041/3.32324/5.33104; LG fwhm=327.14 1.60124/3.32964/5.35144;
    LG fwhm=300 1.60510/3.36341/5.48647; vectorial opuesta 147.886/89.588/84.197; lineal 60.774/38.061/35.723 nm.
    ⇒ la LG igualada por curvatura reproduce el CRB vectorial a 0.05/0.19/0.38 %.

## Cámara (Balzarotti Eq. S59–S63)
W7. PSF gaussiana σ_PSF = 100 nm, píxel 100 nm, ventana 9×9, sin fondo, N = 400: CRB = 5.2045 nm
    (ideal σ/√N = 5.000). Con píxel 10 nm y ventana 121×121: CRB/(σ/√N) = 1.00042.
W8. Con SBR_c = 500 por píxel (9×9, N = 600) CRB = 4.964 nm; para 5 nm hacen falta 591 fotones
    (Balzarotti Fig. 3: ~600).
W9. MINFLUX (L = 50, fwhm = 300, S31 en el centro) necesita 13.0/13.5/14.2/15.4/17.9 fotones para 5 nm
    con SBR inf/50/20/10/5.

## Experimentos (semilla 42; LG fwhm = 300; multinomial)
W10. MINFLUX iterativo, N_total = 1000, 4 iteraciones de 250 fotones, L = 150→82.55→45.43→25 nm
    (geométrico), posición verdadera uniforme en disco de radio 37.5 nm, recentrado del TCP en la
    estimación previa, MLE en disco de radio 0.75 L_k, n_rep = 500: σ final = 0.497 ± 0.011 nm sin
    fondo (0.613 ± 0.014 con SBR = 10) vs cámara σ_PSF/√N_total = 3.162 nm. Sin recentrado σ por
    iteración = 4.34/3.20/3.23/8.06 nm.
W11. σ iterativo vs N_total (250…8000, sin fondo): pendiente log-log −0.526; cociente con la cámara
    0.15–0.17.
W12. L óptimo con cero imperfecto (pedestal gaussiano eps, N = 100, sin fondo, CRB en el centro):
    eps = 0.002/0.01/0.05/0.15 → L_opt = 10.48/23.28/50.34/80.81 nm, CRB(L_opt) = 0.746/1.679/3.878/7.287 nm;
    con eps = 0 el CRB es monótono creciente en L (sin mínimo interior). Empírico: L_opt ≈ 0.78 fwhm √eps.
W13. Desalineación (L = 100, N = 500, SBR = 10, cada cero desplazado δ en dirección aleatoria, 20
    patrones × 200 reps, MLE): en el centro, sesgo medio del estimador ingenuo (supone TCP ideal)
    1.58/4.01/8.31 nm para δ = 2/5/10, frente al honesto (conoce los centros) 0.19/0.17/0.18 nm (≈ piso
    de ruido MC 0.16); σ del ingenuo casi no cambia (1.85→1.98): la pérdida es sesgo.
