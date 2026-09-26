# Ronda 2 — Worker 2: dona vectorial Richards-Wolf (R1)

Archivos (propios, nuevos): `src/donutloc/vectorial.py`, `tests/test_vectorial.py`. No se tocó
ningún otro archivo; no hay commit. Implementación propia desde la Eq. 1 de Caprile2022 tal como
la transcribe la nota B §2.1 y las integrales de §2.2 (no se usó la E_y de la Eq. 15 impresa ni
código de PyFocus). Nada de `docs/private/`.

## API implementada
- `focal_field(rho, phi, z=0, wavelength=640., NA=1.4, n=1.518, filling=5/3, polarization="circular", handedness=+1, charge=1, n_theta=801, pol_angle=0., method="bessel", n_phi=200)` → (Ex, Ey, Ez) complejos sin normalizar (prefactor −ikf e^{−ikf} omitido). `method="2d"`: cuadratura directa 2D (θ, φ') de la Eq. 1 (control independiente).
- `intensity(x, y, z=0, **opt)`: |E|² vectorizado (broadcast); valida las claves de opción (TypeError si hay una desconocida).
- `radial_profile(rho, **opt)`: solo circular (±); ValueError para "linear".
- `zero_depth(**opt)` = I(0)/I_max (búsqueda radial para circular: grilla de 4 nm + refinamiento acotado; lineal: grilla polar 2D + Nelder-Mead). I_max memoizado según las opciones.
- `peak_to_peak_diameter(**opt)` = 2 ρ del máximo de I(ρ) (lineal: del promedio azimutal).
- `zero_curvature(**opt)`: c de I/I_max ≈ c0 + c ρ² + d ρ⁴, ajuste LSQ en ρ = 1..10 nm (lineal: promedio azimutal de 36 direcciones).
- `lg_equivalent_fwhm(beam_opts, match="curvature"|"diameter")`: curvatura → fwhm = sqrt(4e ln2 / c); diámetro → fwhm = D_pp·sqrt(ln2), es decir D_pp = 2·beams.ring_radius(fwhm).
- `make_vectorial_beam(rho_max=1500., d_rho=1., eps=0., mode="interp"|"exact", grid_step=5., fwhm_eps=None, **opt)`: beam f(x, y) con pico 1.
  - Circular: `np.interp` **en s = ρ²** (así I ≈ cρ² se conserva exacto debajo del primer nodo, que importa para crb_limit con r0 = 1e-3 nm); vale 0 para ρ > rho_max.
  - Lineal: `RegularGridInterpolator` lineal sobre una grilla cartesiana de paso grid_step, con 0 fuera de ella.
  - `mode="exact"` evalúa las integrales en cada llamada.
  - `eps` suma eps·exp(−4 ln2 ρ²/fwhm_eps²) después de normalizar; fwhm_eps es por defecto el fwhm LG equivalente por diámetro.
- `compare_crb_vectorial_vs_lg(L_list=(50,100,150), N=100, sbr=None, **opt)`:
  - Usa solo `patterns.tcp_centers`, `photons.make_model`, `fisher.crb_limit` (defaults) y `beams.make_beam`.
  - La dona vectorial entra con `mode="exact"`.
  - Devuelve un dict con: L, crb_vectorial, crb_lg_curvature, crb_lg_diameter, crb_lg_300, fwhm_curvature, fwhm_diameter, zero_depth, N, sbr, opts.

## Derivación (en el docstring del módulo)
Parto de e = (x̂ + i s ŷ)/√2 y un vórtice e^{ilφ'}.
- Proyecciones sobre la base local: e·ρ̂' = e^{isφ'}/√2 y e·φ̂' = i s e^{isφ'}/√2. Por eso E_0 = a e^{imφ'}(i s φ̂' + θ̂)/√2, con m = l + s.
- Aplico la identidad de la Eq. 12 de Caprile y defino tres integrales:
  - P = ∫g (cosθ − s) J_{m+1};
  - Q = ∫g (cosθ + s) J_{m−1};
  - R = ∫g sinθ J_m.
- Componentes del campo:
  - Ex = [i^{m+1}e^{i(m+1)φ}P + i^{m−1}e^{i(m−1)φ}Q]/(2√2);
  - Ey = (i/2√2)[−i^{m+1}e^{i(m+1)φ}P + i^{m−1}e^{i(m−1)φ}Q];
  - Ez = −i^m e^{imφ}R/√2.
- Intensidad: I = (|P|² + |Q|²)/4 + |R|²/2, sin dependencia en φ.

Casos:
- **l = 1, s = +1 (m = 2).** Da exactamente las fórmulas de la nota B §2.2, con P = −I_B, Q = I_A, R = I_C. Órdenes: J1 y J3 en las transversales, J2 en la axial.
- **l = 1, s = −1 (m = 0), mano opuesta.**
  - Ex = (i/2√2)(I_A e^{iφ} − I_B' e^{−iφ}).
  - Ey = −(1/2√2)(I_A e^{iφ} + I_B' e^{−iφ}).
  - Ez = −I_C'/√2.
  - Integrales: I_A = ∫g(1+cosθ)J1, I_B' = ∫g(1−cosθ)J1, I_C' = ∫g sinθ J0.
  - Las transversales son J1 puras y la axial es J0, así que I(0) = |I_C'(0)|²/2.
- **Lineal.** Uso (cos γ, sin γ) = e^{−iγ}(x̂+iŷ)/2 + e^{iγ}(x̂−iŷ)/2, así que E_lin = [e^{−iγ}E(s=+1) + e^{iγ}E(s=−1)]/√2.
  - Es una suma de campos (no de intensidades) con armónicos en φ explícitos, y es exacta.
  - **Desvío declarado:** la tarea pedía una integral 2D (ρ', φ') para la lineal. Esa integral está implementada (`method="2d"`) y es el control, pero por defecto uso esta reducción armónica exacta. La coincidencia entre las dos es de ~1e-15, y la armónica es mucho más rápida (la tabla 2D del beam sería inviable con la cuadratura 2D).

## Resultados (NA 1.4, λ = 640 nm, F = 5/3, n_theta = 801)

| | n = 1.518 (default) | n = 1.5 (nota B) |
|---|---|---|
| zero_depth, mano correcta | 0.0 exacto (I(0) = 0 en punto flotante) | 0.0 exacto |
| zero_depth, mano opuesta | 0.84524 | 0.86896 (nota: 0.869) |
| zero_depth, lineal (γ = 0) | 0.37165 | 0.38231 (nota: 0.382) |
| D_pp, correcta | 384.67 nm | 385.32 nm (nota: 385) |
| curvatura c, correcta | 7.0422e-5 nm⁻² | 7.0181e-5 nm⁻² (nota: 6.96e-5) |
| fwhm LG equivalente, por curvatura | 327.14 nm | 327.70 nm |
| fwhm LG equivalente, por diámetro | 320.26 nm | 320.80 nm |

D_pp en función de F (n = 1.5):

| F | D_pp | Nota B (propia) | Caprile p. 20 |
|---|---|---|---|
| 0.5 | 515.5 nm | 516 nm | 508 nm |
| 1 | 401.3 nm | 401 nm | 428 nm |
| 2 | 382.9 nm | 383 nm | 392 nm |
| 100 (≈ uniforme) | 377.6 nm | 378 nm | 380 nm |

- **Curvatura.** Mi c es el límite analítico exacto ρ → 0, (k²/16)|∫g(1+cosθ)sinθ|²/I_max. El ajuste coincide con ese límite a 2e-7 relativo.
  - El 6.96e-5 de la nota corresponde a I/(I_max ρ²) evaluado en ρ ≈ 20 nm: 6.949e-5 allí, y 6.863e-5 en 30 nm. La diferencia de 0.8 % se explica por el rango de ajuste, no es un error.
- **Mano opuesta.** I(ρ) no es una dona sino un hundimiento leve: 0.871 en 0, con máximo en ρ = 129 nm.
  - El "D_pp" que da la función (258.2 nm con n = 1.5; 266.8 nm con n = 1.518) es 2× ese radio.
  - Sus fwhm "equivalentes" (650 y 215 nm con n = 1.5) no tienen sentido físico: la LG no tiene fondo en el cero.
- **Lineal** (n = 1.518): D_pp del promedio azimutal = 344.4 nm; c = 3.77e-5; fwhm equivalentes 447.1 y 286.7 nm. Con n = 1.5: 343.1 nm y c = 3.65e-5.

### CRB en el límite central (fisher.crb_limit), TCP, N = 100, sin fondo, nm

"LG equivalente": LG del proyecto (Balzarotti S17) cuyo fwhm iguala **(a)** la curvatura en el cero (4e ln2/fwhm² = c) o **(b)** el diámetro pico a pico (fwhm/√ln2 = D_pp). Se agrega la LG de fwhm = 300 como referencia.

Defaults (n = 1.518):

| L | vectorial correcta | LG por curvatura (327.14) | LG por diámetro (320.26) | LG 300 | vectorial opuesta | vectorial lineal |
|---|---|---|---|---|---|---|
| 50 | 1.60041 | 1.60124 | 1.60212 | 1.60510 | 147.886 | 60.774 |
| 100 | 3.32324 | 3.32964 | 3.33734 | 3.36341 | 89.588 | 38.061 |
| 150 | 5.33104 | 5.35144 | 5.38180 | 5.48647 | 84.197 | 35.723 |

Nota B (n = 1.5):

| L | vectorial correcta | LG por curvatura (327.70) | LG por diámetro (320.80) | LG 300 | opuesta | lineal |
|---|---|---|---|---|---|---|
| 50 | 1.60035 | 1.60117 | 1.60205 | 1.60510 | 169.805 | 70.540 |
| 100 | 3.32265 | 3.32904 | 3.33671 | 3.36341 | 103.296 | 43.824 |
| 150 | 5.32872 | 5.34907 | 5.37929 | 5.48647 | 98.498 | 40.806 |

Respuesta a R1:
- Con la mano correcta, la dona vectorial y la LG equivalente por curvatura dan CRB centrales que difieren 0.05 % en L = 50, 0.19 % en L = 100 y 0.38 % en L = 150. Por diámetro: 0.11 %, 0.42 % y 0.95 %. Contra la LG de fwhm = 300: 0.3 %, 1.2 % y 2.9 %.
- Con la mano opuesta o con polarización lineal, el fondo en el cero (84.5 % y 37.2 % del pico) destruye la ventaja de MINFLUX: el CRB central es 92× y 38× peor que el de la LG en L = 50.
- En estos casos el CRB **baja** al crecer L, porque el fondo domina y un L mayor da más contraste.

### Otros chequeos
- **Bessel frente a integral 2D directa** (n_theta = 201, n_phi = 200, mismos nodos en θ).
  - Casos: circular ±1 y lineal (γ = 0.4); ρ ∈ {0, 40, 150, 330} nm; z ∈ {0, 250} nm.
  - Diferencia máxima por componente / max|E|: ≤ 1e-15 (en exploración con 5 puntos, 6.1e-16 a 9.3e-16).
- **Interpolación del beam** (mano correcta, crb_limit N = 100). La interpolación con d_rho = 1 frente a "exact" difiere 2.4e-5 relativo en L = 50 (1.600453 frente a 1.600415) y 2.8e-5 en L = 150.
- **Energía.** ∫|E|² 2πρ dρ en z = 0 coincide para s = ±1: diferencia relativa 3e-5 a R = 2 µm y 3e-7 a 20 µm.
  - Converge al valor de Parseval (2π/k²)∫a² sinθ dθ (0.969 a 2 µm, 0.997 a 20 µm, cola ~1/R). Parseval no depende de s porque |E_0|² = a² para las dos manos.
  - También la parte transversal sola coincide (0.7807 frente a 0.7807 de Parseval a 20 µm). La componente z lleva ~22 % del flujo |E|² en ambos casos.
- **Regla general ls = +1:** charge = −1 con handedness = −1 da cero perfecto (en el test).

## Tests
`python -m unittest discover -s tests -p "test_vectorial.py"` (salida real):
```
..............
----------------------------------------------------------------------
Ran 14 tests in 13.190s

OK
```
Cubren:
- las profundidades del cero (< 1e-10; 0.869 ± 0.01; 0.382 ± 0.01) y el cero fuera de foco (z = 300);
- los órdenes de Bessel de la mano opuesta, contra sumas independientes con J0 y J1;
- D_pp (F = 100: 380 ± 3 %; F = 2: 392 ± 3 %);
- la curvatura (7.0e-5 ± 5 % y contra el límite analítico);
- la simetría de revolución a φ = 0, 1, 2 rad (1e-10);
- la asimetría de la lineal y el error de `radial_profile` con "linear";
- los argumentos inválidos;
- Bessel = 2D (< 1e-6) en tres polarizaciones y dos z;
- el flujo igual para ±1 (1e-3) y cercano a Parseval (5 %);
- el beam: pico 1, 0 en el centro y fuera de rho_max, I/ρ² → c, Σp = 1, eps;
- el beam lineal en grilla;
- las fórmulas de fwhm equivalente;
- `compare_crb` en L = 50 (vectorial / LG por curvatura dentro de 1 %).

## Tiempos (esta máquina)
- zero_depth circular: ~0.5 s en la primera llamada (memoizado después). Lineal: ~3 s.
- compare_crb_vectorial_vs_lg para L ∈ {50, 100, 150}: ~1 s para circular y ~16 s para lineal.
- Suite test_vectorial: 13.2 s.

## Candidatos a paper_numbers
- `vectorial_zero_depth_correct` = 0.0 (I(0) exactamente 0; se reporta < 1e-10).
- `vectorial_zero_depth_wrong_handedness` = 0.845 (n = 1.518) / 0.869 (n = 1.5).
- Hay que decidir qué n se publica: el default del proyecto es 1.518.

## Abierto / límites
- D_pp de Caprile en F = 1 (428 nm) no se reproduce: da 401 nm, como ya observó la nota B. El test de D_pp usa solo F = 2 y F → ∞.
- `peak_to_peak_diameter` y `lg_equivalent_fwhm` para la mano opuesta o la lineal devuelven números definidos pero sin significado de "dona".
- El beam interpolado de la lineal (grilla de 5 nm, lineal) tiene un error de ~1e-3 relativo cerca del pico. Para el CRB conviene `mode="exact"`.
