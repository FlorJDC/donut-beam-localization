# B. Óptica de haces dona (2D) y dona 3D / top-hat para simulación tipo MINFLUX

Nota de literatura. Fuentes (número de página = página del PDF, empezando en 1):

| Clave | Archivo |
|---|---|
| [Lopez2023] | `Lopez et al. - 2023 - Optimization and characterization of toroidal foci ... tut.pdf` (preprint aceptado en JOSA B; 16 pp.) |
| [Caprile2022] | `Caprile et al. - 2022 - PyFocus – A Python package for vectorial calculations ...pdf` (42 pp.) |
| [Tarkowski] | `Tarkowski  - MINFLUX excitation patterns.pdf` (18 pp.) |
| [Gwosch2020] | `MINFLUX_3D_Gwosh.pdf` (Nature Methods 17, 217–224; 14 pp.) |
| [Pape2020] | `MINFLUX_3D_pnas.2009364117.pdf` (PNAS 117, 20607; 8 pp.) |

Convenciones de esta nota: $\lambda$ es la longitud de onda en vacío, $n$ el índice del medio de enfoque, $k = 2\pi n/\lambda$ el número de onda en el medio, $\alpha$ el semiángulo máximo de apertura y $\mathrm{NA} = n\sin\alpha$. Lo marcado **[propio]** es cálculo o derivación mía, no del paper. Lo marcado **[estándar, no en los papers asignados]** es teoría de libro de texto que ningún PDF asignado trae explícitamente; un verificador no lo va a encontrar en esas fuentes.

> Advertencia general: [Lopez2023] no tiene ecuaciones; es un tutorial experimental. [Pape2020] casi no tiene óptica (los detalles quedan en un SI Appendix que no está en la carpeta). Ningún paper asignado da la fórmula del modo Laguerre-Gauss ni el radio del disco de la máscara top-hat. Las ecuaciones cuantitativas vienen de [Caprile2022] (vectorial), [Tarkowski] (modelo analítico de la dona) y [Gwosch2020] (modelos cuadráticos y parámetros).

---

## 1. Aproximación escalar paraxial: modo LG$_0^1$ (dona)

### 1.1 Modelo LG$_0^1$ [estándar, no en los papers asignados]

Campo del modo Laguerre-Gauss con $p=0$ y $\ell=1$ (carga topológica 1), en el foco ($z=0$):
$$u(r,\varphi) \propto \frac{\sqrt2\, r}{w}\, e^{-r^2/w^2}\, e^{i\varphi}.$$

Intensidad normalizada a potencia $P$:
$$I(r) = \frac{4P}{\pi w^4}\, r^2\, e^{-2r^2/w^2}, \qquad \int_0^\infty I\,2\pi r\,dr = P.$$

- Radio del máximo: $r_{pk} = w/\sqrt2$. Distancia pico-pico: $D_{pp} = \sqrt2\,w$.
- Intensidad de pico: $I_{pk} = 2P/(\pi e\, w^2)$.
- Forma normalizada: $I/I_{pk} = (2r^2/w^2)\,e^{\,1-2r^2/w^2}$.
- Cerca del cero: $I \simeq \dfrac{4P}{\pi w^4}\, r^2$, o sea $I/I_{pk} \simeq 2e\, r^2/w^2$.
- Fuera del foco: $w(z) = w_0\sqrt{1+(z/z_R)^2}$ con $z_R = \pi w_0^2 n/\lambda$; la fase de Gouy vale $2\arctan(z/z_R)$ (orden $|\ell|+2p+1 = 2$). El cero sobre el eje se mantiene para todo $z$: la dona 2D **no** da información axial.

### 1.2 Modelo analítico usado en [Tarkowski]

[Tarkowski, p. 4, Eq. 4] usa (en 3D)
$$I(\mathbf r) = I_0\left(\frac{x^2}{d_x^2}+\frac{y^2}{d_y^2}+\frac{z^2}{d_z^2}\right)\exp\!\left[-4\ln(2)^{(3?)}\left(\frac{x^2}{d_x^2}+\frac{y^2}{d_y^2}+\frac{z^2}{d_z^2}\right)\right],$$
donde $d_{x,y,z}$ son las "FWHM de las crestas de la dona" en cada dirección [Tarkowski, p. 5]. Para 2D usa $\mathbf r = (x,y,0)$ con $d_x = d_y = 360$ nm [Tarkowski, p. 5]. Para 3D isotrópico usa $d_x = d_y = d_z = 360$ nm [Tarkowski, p. 6], y para el caso anisotrópico $d_z = \alpha\cdot 360$ nm con $\alpha = 1.5$ [Tarkowski, p. 13]. El exponente impreso es ambiguo; ver §7.

**Equivalencia con LG$_0^1$ [propio].** Si el exponente es $-4\ln 2\,(r^2/d^2)$, la Eq. 4 es exactamente un LG$_0^1$ con $w = d/\sqrt{2\ln 2}$, es decir, $d$ es la FWHM de un gaussiano de cintura $w$. Consecuencias:
- $r_{pk} = d/(2\sqrt{\ln 2}) \approx 0.60\,d$.
- $D_{pp} \approx 1.20\,d$, lo que da 432 nm para $d = 360$ nm.
- Cerca del cero: $I/I_{pk} \simeq 4e\ln 2\; r^2/d^2 \approx 7.54\, r^2/d^2$, que da $5.8\times10^{-5}\,\mathrm{nm}^{-2}$ para $d = 360$ nm.

### 1.3 Aproximación cuadrática cerca del cero

- [Gwosch2020, p. 1] parte de un "cero de intensidad flanqueado por un perfil cuadrático" para derivar la cota $\sigma \ge L/(4\sqrt N)$.
- Para la dona 3D asumen la PSF $I(x,y,z) = a\,(x^2+y^2+z^2)$, con $a$ constante [Gwosch2020, p. 10]. La misma forma aparece en la tabla de funciones de haz, fila "3D donut": $f = (x-x_b)^2+(y-y_b)^2+(z-z_b)^2$ [Gwosch2020, p. 9].
- Para la dona 2D la tabla usa la PSF 2D experimental por un gaussiano en $z$: $f = f_{PSF2d}(x-x_b, y-y_b)\,e^{-(z-z_b)^2/2\sigma_Z^2}/(\sqrt{2\pi}\sigma_Z)$, con $\sigma_Z = 1000\,\mathrm{nm}/(2\sqrt{2\ln 2})$ [Gwosch2020, p. 9].

---

## 2. Enfoque vectorial (Richards-Wolf / Debye)

### 2.1 Integral general [Caprile2022]

Campo cerca del foco de una lente aplanática sin aberraciones [Caprile2022, p. 4, Eq. 1]:
$$\mathbf E_f(\rho,\phi,z) = \frac{-ikf\,e^{-ikf}}{2\pi}\int_0^{2\pi}\!\!\int_0^{\alpha}\mathbf E_0\, e^{i\mathbf k\cdot\mathbf r}\sqrt{\cos\theta}\,\sin\theta\,d\theta\,d\phi'.$$

Ingredientes:
- Producto escalar: $\mathbf k\cdot\mathbf r = k\,[z\cos\theta + \rho\sin\theta\cos(\phi'-\phi)]$ [Caprile2022, p. 4].
- Condición del seno: el rayo que entra a radio $\rho'$ sale con ángulo $\theta$ tal que $\rho' = f\sin\theta$ [Caprile2022, p. 4].
- Amplitud en el frente focal: $\mathbf E_0 = (\mathbf E_i\cdot\hat{\boldsymbol\phi}')\,\hat{\boldsymbol\phi}' + (\mathbf E_i\cdot\hat{\boldsymbol\rho}')\,\hat{\boldsymbol\theta}$. La componente $s$ conserva su dirección y la componente $p$ pasa a apuntar según $\hat{\boldsymbol\theta}$ [Caprile2022, p. 4].
- Vectores unitarios: $\hat{\boldsymbol\rho}' = (\cos\phi', \sin\phi', 0)$, $\hat{\boldsymbol\phi}' = (-\sin\phi', \cos\phi', 0)$, $\hat{\boldsymbol\theta} = (\cos\theta\cos\phi', \cos\theta\sin\phi', -\sin\theta)$ [Caprile2022, p. 5].
- Intensidad: $I = |E_{fx}|^2+|E_{fy}|^2+|E_{fz}|^2$ [Caprile2022, p. 5]. El paper escribe $E^2$; ver §7.
- **Apodización aplanática:** el factor $\sqrt{\cos\theta}$ de la Eq. 1 [Caprile2022, p. 4].

Parámetros geométricos:
- $\mathrm{NA} = n\sin\alpha$ con $\alpha = \arcsin(h/f)$, donde $h$ es el radio de apertura [Caprile2022, p. 11].
- $f = h\,n/\mathrm{NA}$ [Caprile2022, p. 14].

Polarizaciones base:
- Polarización elíptica general $\mathbf E_i = E_i(c_x\hat x + c_y\hat y)$, y $\mathbf E_f = c_x\mathbf E_f^{(x)} + c_y\mathbf E_f^{(y)}$ [Caprile2022, p. 5]. Las integrales 2D para $x$ e $y$ son las [Caprile2022, p. 5, Eqs. 2–3].
- Haz gaussiano: $\mathbf E_i = \mathbf A\, e^{-(\rho'/w_0)^2}$, con $w_0$ el radio donde la **amplitud** cae a $e^{-1}$. Con la condición del seno queda $\mathbf A\,e^{-(f\sin\theta/w_0)^2}$ [Caprile2022, p. 6, Eq. 4].
- Parametrización de la polarización de entrada: $\mathbf E_e = E_e(\cos\gamma\,\hat x + \sin\gamma\, e^{i\beta}\hat y)$. $\gamma = 45^\circ$, $\beta = 90^\circ$ es "circular derecha (vista desde la fuente)" [Caprile2022, p. 12].

Máscara de fase:
- Sin propagación entre máscara y objetivo: $\mathbf E_i = \mathbf E_e\, e^{i\psi(\rho',\phi')}$ [Caprile2022, p. 7, Eq. 6].
- Con propagación: integral de Fraunhofer [Caprile2022, p. 6, Eq. 5].
- Vórtice 0–2π: $e^{i\psi} = e^{\pm i\phi''}$, donde el signo es la handedness; PyFocus usa $+$ [Caprile2022, p. 7].
- Sin propagación: $\mathbf E_i = \mathbf A\, e^{-(\rho'/w_0)^2}e^{i\phi'}$ [Caprile2022, p. 8, Eq. 8].
- La distancia de propagación $L$ (0.2, 2, 10 y 50 m) cambia mucho el perfil incidente pero **no** cambia apreciablemente el foco toroidal, porque el factor $e^{i\phi}$ sobrevive a la propagación. Por eso `propagation=False` es el caso por defecto [Caprile2022, pp. 18–19, Fig. 9].

Identidad que reduce la integral en $\phi'$ a funciones de Bessel [Caprile2022, p. 31, Eq. 12]:
$$\int_0^{2\pi} e^{im\phi'}e^{i\epsilon\cos(\phi'-\phi)}\,d\phi' = 2\pi\, i^m e^{im\phi} J_m(\epsilon).$$

- Gaussiano sin máscara: [Caprile2022, p. 31, Eq. 13] e integrales $I_0, I_1, I_2$ con $J_0$, $J_1$ y $J_2$ [Caprile2022, pp. 31–32, Eq. 14].
- Gaussiano con VP: [Caprile2022, p. 32, Eq. 15] con integrales $I'_1 \ldots I'_5$ [Caprile2022, p. 32, Eq. 16; p. 33, Eq. 17 sin propagación]. **Atención:** encontré un error en la componente $E_{fy}$ de la Eq. 15 tal como está impresa; ver §2.2 y §7.

### 2.2 Fórmulas implementables: vórtice $\ell=+1$ y polarización circular "correcta" [propio, verificado numéricamente]

Derivé las fórmulas de la Eq. 1 de [Caprile2022] con la identidad de la Eq. 12, para:
- $\mathbf E_i = A\,a(\theta)\,e^{i\phi'}\,(\hat x + i\hat y)/\sqrt2$, es decir $\gamma = 45^\circ$, $\beta = +90^\circ$ en la convención PyFocus;
- una amplitud de pupila $a(\theta)$ arbitraria con simetría de revolución; para un gaussiano, $a(\theta) = e^{-(f\sin\theta/w_0)^2}$.

Defino
$$g(\theta) = a(\theta)\sqrt{\cos\theta}\,\sin\theta\; e^{ikz\cos\theta}, \qquad u = k\rho\sin\theta,$$
$$I_A = \int_0^\alpha g\,(1+\cos\theta)\,J_1(u)\,d\theta,\qquad I_B = \int_0^\alpha g\,(1-\cos\theta)\,J_3(u)\,d\theta,\qquad I_C = \int_0^\alpha g\,\sin\theta\,J_2(u)\,d\theta.$$

Entonces, con $C_0 = -ikf e^{-ikf}A$:
$$E_x = C_0\,\frac{i}{2\sqrt2}\left(I_A e^{i\phi} + I_B e^{3i\phi}\right),\qquad E_y = -C_0\,\frac{1}{2\sqrt2}\left(I_A e^{i\phi} - I_B e^{3i\phi}\right),\qquad E_z = C_0\,\frac{1}{\sqrt2}\, I_C\, e^{2i\phi},$$
$$\boxed{\,I(\rho,z) = |C_0|^2\left[\tfrac14\left(|I_A|^2+|I_B|^2\right) + \tfrac12|I_C|^2\right]\,}$$

Propiedades:
- La intensidad tiene simetría de revolución exacta: no depende de $\phi$.
- $I(0,z) = 0$ para **todo** $z$, porque $J_1(0) = J_2(0) = J_3(0) = 0$.
- Para comparar: $E_x$ coincide con la Eq. 15 de [Caprile2022] ($I'_1 \equiv I_A$, $I'_3 \equiv I_B$, $I'_5 \equiv I_C$ salvo constantes), y $E_z$ también, pero $E_y$ no. Ver §8.

Regla general de handedness [propio, consistente con Lopez y Caprile]:
- Con vórtice $e^{i\ell\phi'}$ y polarización $(\hat x + i\sigma\hat y)/\sqrt2$, con $\sigma = \pm1$, la componente $p$ lleva $e^{i(\ell+\sigma)\phi'}$. Por eso $E_z \propto J_{\ell+\sigma}$, y las transversales $\propto J_{\ell+\sigma\pm1}$.
- Para $\ell = 1$ y $\sigma = +1$: $E_z \propto J_2$, que se anula en el centro. Cero perfecto.
- Para $\ell = 1$ y $\sigma = -1$: $E_z \propto J_0$, máximo en el centro. Las transversales (órdenes $\ell+\sigma\pm1 = \pm1$, es decir $\propto J_1$) siguen siendo nulas en el centro.

### 2.3 Handedness y polarización: qué dicen los papers

- Hace falta polarización circular "en la dirección de la fase helicoidal" / "con la misma mano que la VPP" [Lopez2023, pp. 3–4; p. 7; p. 9; p. 10].
- Las componentes laterales tienen mínimo en el centro para **cualquier** polarización. La axial solo tiene mínimo para circular de la misma mano; con la mano opuesta tiene un **máximo central**. Con elíptica, el mínimo se mantiene pero su valor queda por encima del fondo. A NA baja el efecto es menor porque $E_z$ es relativamente débil [Lopez2023, p. 11].
- Experimentalmente la mano de la VPP no se ve a simple vista. Se ajusta la circularidad y, si no aparece el cero, se da vuelta la VPP [Lopez2023, p. 11].
- En PyFocus, con la VP de mano positiva hay que usar circular derecha (`gamma=45, beta=90`). Con circular izquierda, $E_z$ tiene un máximo central que degrada mucho el contraste [Caprile2022, p. 17; p. 35, App. A.4, Fig. 18].
- Figuras de referencia: [Lopez2023, p. 12, Fig. 5] (simulación y experimento para circular derecha, circular izquierda, elíptica y lineal).
- Contraste de polarización medido con analizador rotante: $(V_{max}-V_{min})/(V_{max}+V_{min})$. El mínimo alcanzable en la práctica es del 1 al 5 % [Lopez2023, p. 8].
- [Gwosch2020, p. 9]: una lámina $\lambda/4$ acromática para polarización circular "ayudó a minimizar la intensidad en el mínimo de la dona".

### 2.4 Dependencia con NA, $n$, $\lambda$ y llenado de pupila

- Factor de llenado: $F = w_0/h$. $F$ grande equivale a iluminación uniforme (se usa toda la NA); bajar $F$ equivale a bajar la NA efectiva [Caprile2022, p. 20].
- Diámetro pico-pico con NA 1.4, $n = 1.5$, $\lambda = 640$ nm, $h = 3$ mm [Caprile2022, p. 20, Fig. 10]:

  | $F$ | Diámetro pico-pico |
  |---|---|
  | 0.5 | 508 nm |
  | 1.0 | 428 nm |
  | 2.0 | 392 nm |
  | uniforme ($F\to\infty$) | 380 nm |

- Cerrar el iris reduce la NA efectiva, agranda la dona y suele mejorar la simetría y la calidad del cero, porque los bordes de un objetivo de NA alta concentran las aberraciones. Ejemplo: NA 1.4 contra NA efectiva 1.2 [Lopez2023, p. 5; pp. 12–13, Fig. 6].
- Escalas: todas las longitudes transversales van como $\lambda/\mathrm{NA}$ vía $u = k\rho\sin\theta = 2\pi n\rho\sin\theta/\lambda$; la forma depende de $\alpha$ y de $F$. [propio]

---

## 3. Dona 3D / top-hat / bottle beam (localización en z)

- **Máscaras en [Gwosch2020, p. 2, Fig. 1a, recuadro]:**
  - "Flat": foco regular.
  - "Vortex": rampa 0–2π, da la dona 2D.
  - "Tophat": disco central con fase $\pi$ y anillo exterior con fase 0; da la dona 3D.
  - Se generan con un SLM, y se conmuta entre foco regular, dona 2D y dona 3D con EOMs [Gwosch2020, p. 9].
  - El dibujo de la PSF top-hat en el plano $xz$ muestra un mínimo central encerrado por lóbulos arriba y abajo en $z$ [Gwosch2020, p. 2, Fig. 1a].
- **Radio del disco $\pi$:** ningún paper asignado lo da.
  - [estándar] En aproximación escalar con pupila uniforme se elige $\rho_0 = R/\sqrt2$, para que las dos zonas lleven igual potencia y el campo en el foco se cancele.
  - [propio] En cálculo vectorial (NA 1.4, $n = 1.5$, circular), la condición de cancelación en el foco da $\rho_0 = 0.7145\,h$ con pupila uniforme y $\rho_0 = 0.682\,h$ con $F = 5/3$ ($w_0 = 5$ mm, $h = 3$ mm).
- **Fórmulas vectoriales top-hat con polarización circular y sin vórtice [propio, verificado contra la integral 2D]:**
  - Amplitud: $a(\theta) = a_G(\theta)\cdot\mathrm{sgn}(f\sin\theta - \rho_0)$, es decir fase $\pi$ para $\rho' < \rho_0$.
  - Integrales: $K_0 = \int g\,(1+\cos\theta)J_0(u)$, $K_2 = \int g\,(1-\cos\theta)J_2(u)$, $K_1 = \int g\sin\theta\,J_1(u)$.
  - Intensidad: $I \propto \tfrac14(|K_0|^2+|K_2|^2)+\tfrac12|K_1|^2$. En el eje solo sobrevive $K_0$.
- **Perfil axial y cuadrático [propio]** (NA 1.4, $\lambda = 640$ nm, $n = 1.5$, $F = 5/3$):
  - Máximos axiales en $z \approx \pm 510$ nm, simétricos.
  - En el plano focal hay un anillo en $\rho \approx 261$ nm con solo **~25 %** de la intensidad de los lóbulos axiales.
  - Cerca del cero, $I/I_{max} \simeq c_z z^2 + c_\rho\rho^2$ con $c_z \approx 1.0\times10^{-5}$ nm$^{-2}$ y $c_\rho \approx 5.1\times10^{-6}$ nm$^{-2}$, o sea $c_z/c_\rho \approx 2$.
  - Conclusión: la top-hat pura **no** es isotrópica. El confinamiento lateral es débil.
- **Cómo se obtiene el mínimo 3D isotrópico:** en la práctica se ajustan las contribuciones de vórtice y top-hat; eso también permite hacer el mínimo anisotrópico a propósito [Tarkowski, p. 13 (cita Gwosch y Wildanger 2009)]. [Gwosch2020] modela su dona 3D como cuadrática isotrópica $a(x^2+y^2+z^2)$ [Gwosch2020, p. 10; p. 9, tabla].
- **Uso operativo** [Gwosch2020, p. 4]:
  - Primero se localiza en $xy$ con el foco regular.
  - Después se ubica la dona 3D arriba y abajo de la posición estimada en $z$.
  - Por último se exploran pares en $x$, $y$, $z$ más el centro.
  - El posicionamiento axial usa una lente varifocal: rango ~400 nm en 50 µs.

---

## 4. Imperfecciones

### 4.1 Profundidad finita del cero
- Criterio experimental: el mínimo central debe igualar el fondo de la imagen y ser angosto y rotacionalmente simétrico [Lopez2023, p. 9]. Con circular derecha lo consiguen: el cero coincide con el fondo experimental [Lopez2023, p. 11].
- Los papers no dan un porcentaje de profundidad de cero. Ver §8 para cuantificaciones **[propias]** en función de la ellipticidad, la descentración y las aberraciones.
- [Gwosch2020, p. 9] midió el frente de onda aberrado con un esquema de segmentación de pupila "para asegurar un mínimo profundo".

### 4.2 Aberraciones
- PyFocus modela aberraciones como máscara $e^{i(\phi' + A_\Phi\Phi(\rho',\phi'))}$ con polinomios de Zernike normalizados, $r = \rho'/h$ [Caprile2022, p. 25]:
  - astigmatismo $Z_2^{-2} = \sqrt6\,r^2\sin 2\phi'$;
  - "esférica" $Z_2^0 = \sqrt3\,2r^2$ (en realidad es desenfoque; ver §7);
  - coma $Z_3^{-1} = \sqrt8(3r^2-2)r\sin\phi'$;
  - trébol $Z_3^{-3} = \sqrt8\, r^3\sin3\phi'$.
- Muestran focos con $A_\Phi = 1.5$ rad, cada uno fuertemente distorsionado [Caprile2022, p. 25, Fig. 14].
- Recomendaciones de [Lopez2023, p. 4]:
  - Dicroicos gruesos (3–5 mm, planitud $\le 0.25$ ondas/pulgada P-V), pegados y no apretados con tornillos, para evitar astigmatismo.
  - Espejos de plata a 45°; los dieléctricos alteran los estados mixtos de polarización.

### 4.3 Desalineación
- Haz y VPP descentrados juntos una distancia $d$ (con $D = d/h$): se reemplazan las coordenadas por $\tilde\rho = \sqrt{\rho'^2+d^2-2\rho' d\cos\phi'}$ y $\tilde\phi = \arctan[\rho'\sin\phi'/(\rho'\cos\phi'-d)]$ [Caprile2022, p. 21, Eq. 10].
- Campo incidente: $\mathbf E_i = \mathbf E_G(\tilde\rho)\,e^{i\tilde\phi}$ [Caprile2022, p. 21]. El foco se distorsiona al crecer $D$ [Caprile2022, pp. 22–23, Fig. 12]; con $w_0\to\infty$ solo se desplaza la fase vortical [Caprile2022, p. 23].
- Haz inclinado un ángulo $\alpha$: término de fase $e^{ik\rho'\sin\alpha\cos\phi'}$, y el foco se desplaza $d = f\tan\alpha$ [Caprile2022, p. 23].
  - Ángulos de 0.0009°, 0.0018° y 0.0027° dan desplazamientos de 50, 100 y 150 nm sin deformar la dona [Caprile2022, pp. 23–24, Fig. 13].
  - Con ángulos grandes el foco además se inclina [Caprile2022, p. 24].
- Protocolo de alineación de la VPP [Lopez2023, pp. 9–10, Fig. 4]:
  - Primero, alineación gruesa $(x, y)$ mirando el haz con cámara.
  - Después, ajuste fino de $x, y, \theta, \varphi$ mirando la PSF medida con perlas.
  - La inclinación no se nota en la alineación gruesa.
- VPP usada: fused silica, 64 escalones 0–2π, diseñada para 633 nm y usada a 642 nm [Lopez2023, p. 5; p. 4]. Eficiencia: >95 % para la VPP contra 70–90 % para un SLM [Lopez2023, p. 3].

### 4.4 Cubreobjetos / interfaz / multicapa
- Campo con una multicapa de espesor despreciable ubicada en $z_{int}$ [Caprile2022, p. 8, Eq. 9]:
  - $\mathbf E_f = \mathbf E_{f0} + \mathbf E_r$ para $z < z_{int}$;
  - $\mathbf E_f = \mathbf E_t$ para $z > z_{int}$.
- Los coeficientes $r_{s,p}$ y $t_{s,p}$ salen de una matriz de transferencia; para una sola interfaz son los de Fresnel [Caprile2022, p. 9].
- Integrales explícitas para $\mathbf E_r$ y $\mathbf E_t$ [Caprile2022, pp. 34–35, Eqs. 19–21], con $\sin\theta_t = (n_1/n_2)\sin\theta$ y $\cos\theta_t = \sqrt{1-(n_1/n_2)^2\sin^2\theta}$ [Caprile2022, p. 35].
- Ejemplo: vidrio $n_1 = 1.5$ / agua $n_2 = 1.33$ con la interfaz en el plano focal [Caprile2022, pp. 26–27, Fig. 15].
- Dona en TIR: se obstruyen los ángulos $\rho' < f\sin\theta_c$ [Caprile2022, p. 27, Eq. 11]. La dona toroidal se mantiene en el plano focal, con polarización circular hasta ~150 nm del centro [Caprile2022, p. 29, Fig. 17].
- Desajuste de índice en 3D MINFLUX [Gwosch2020]:
  - El desajuste vidrio-agua comprime $L_z$ del TCP; eso hace que $\sigma_z$ salga algo mejor que $\sigma_{xy}$ [Gwosch2020, p. 4].
  - Aplicaron un **factor de escala de 0.70** a todas las $z$ estimadas, confirmado con simulaciones [Gwosch2020, p. 9].
- [Lopez2023, p. 5] mide la PSF con perlas en medio acuoso para reproducir las condiciones reales y reducir el desajuste de índices.

---

## 5. Valores típicos de parámetros

| Parámetro | Valor | Fuente |
|---|---|---|
| $\lambda$ excitación | 642 nm (láser); 640 nm (simulación) | [Lopez2023, p. 4; p. 12]; [Caprile2022, p. 15] |
| $\lambda$ excitación MINFLUX | 642 nm (o 560 nm) | [Gwosch2020, p. 9]; [Pape2020, p. 2] |
| NA / $n$ | 1.4, inmersión en aceite; $n = 1.5$ en simulación | [Lopez2023, p. 5; p. 12]; [Caprile2022, p. 15]; [Gwosch2020, p. 9] |
| Radio de apertura $h$ | 3 mm | [Caprile2022, p. 15]; [Lopez2023, p. 12] |
| $w_0$ de entrada | 5 mm ($F \approx 1.67$) | [Caprile2022, p. 15] |
| Diámetro pico-pico de la dona | 380–508 nm según $F$ | [Caprile2022, p. 20] |
| NA efectiva reducida | 1.2 | [Lopez2023, p. 13] |
| FWHM del foco regular | 360 nm (estimador mLMSE); 300 nm (estimador gaussiano); 330 nm lateral y 800 nm axial (tabla de simulación) | [Gwosch2020, p. 9] |
| "FWHM de cresta" en el modelo analítico | $d = 360$ nm | [Tarkowski, p. 5] |
| Potencia | 3–5 µW en el BFP para medir la PSF con perlas | [Lopez2023, p. 7] |
| Potencia MINFLUX | 20–60 µW, es decir 10–50 kW cm$^{-2}$ en el pico de la dona | [Gwosch2020, p. 6; p. 10] |
| Área de la dona | ~3 veces mayor que la del foco confocal | [Gwosch2020, p. 6] |
| Contraste residual de polarización | 1–5 % | [Lopez2023, p. 8] |
| Sonda para medir la PSF | perlas de 40 nm; sonda ≥3 veces menor que los rasgos da ~5 % de error; píxel 5 nm, FOV 1.2 µm | [Lopez2023, p. 5; p. 10] |
| Perlas en MINFLUX 3D | FluoSpheres de 20 nm | [Gwosch2020, p. 10] |
| Rango axial de la dona 3D | ~400 nm (lente varifocal) | [Gwosch2020, p. 4; p. 9] |
| Tamaños de TCP (2D) | $L$ = 300 (regular) → 150 → 90 → 40 nm | [Gwosch2020, p. 2, Fig. 1] |
| Tamaños de TCP (3D) | $L$ = 300 (regular) → 400 → 150 → 90 → 40 nm | [Gwosch2020, p. 5, Fig. 3] |
| Precisión 3D obtenida | $\sigma_{xy} = 2.4$ nm, $\sigma_z = 1.7$ nm | [Gwosch2020, p. 4] |
| Precisión 3D obtenida | 3–6 nm (dos colores); ~3 nm (un color) | [Pape2020, p. 2; p. 5] |

---

## 6. Para la simulación

### 6.1 Modelo analítico (rápido)

- **2D:** LG$_0^1$ con $I/I_{pk} = (2\rho^2/w^2)e^{1-2\rho^2/w^2}$. Equivale a [Tarkowski, p. 4, Eq. 4] con $d = w\sqrt{2\ln2}$.
  - Elegir $w$ de modo que $D_{pp} = \sqrt2 w$ reproduzca el valor vectorial: para NA 1.4 y $\lambda = 640$ nm, $D_{pp} \approx 380$–$390$ nm, lo que da $w \approx 270$ nm.
  - Opción equivalente: calibrar por la curvatura en el cero (§6.3).
  - Para una dona 2D en 3D, multiplicar por un gaussiano axial con $\sigma_Z = 1000/(2\sqrt{2\ln2})$ nm [Gwosch2020, p. 9], o usar $w(z)$ de LG.
- **3D:** $I = I_0\,q\,e^{-4\ln 2\, q}$ con $q = x^2/d_x^2+y^2/d_y^2+z^2/d_z^2$ [Tarkowski, p. 4, Eq. 4], o el cuadrático puro $a(x^2+y^2+z^2)$ [Gwosch2020, p. 10] para $|r| \ll d$.
- **Profundidad del cero:** sumar un término $I_{bg} = \epsilon\, I_{pk}$ (ver §8 para valores de $\epsilon$ razonables).

### 6.2 Modelo vectorial numérico

- **2D, vórtice + circular:** usar las integrales 1D $I_A, I_B, I_C$ de §2.2. No usar la Eq. 15 de [Caprile2022] tal cual para $E_y$.
- **Top-hat:** usar $K_0, K_1, K_2$ de §3.
- **Casos no simétricos** (elíptica, descentración, Zernike, inclinación): integrar la Eq. 1 de [Caprile2022] en 2D, con $\mathbf E_0$ construido como en [Caprile2022, p. 4].
- **Parámetros de referencia:** NA = 1.4, $n = 1.5$, $\lambda = 640$ nm, $h = 3$ mm, $f = hn/\mathrm{NA}$, $w_0 = 5$ mm, $k = 2\pi n/\lambda$.
- **Discretización:**
  - En $\theta$: trapecio con ≥800 puntos. [propio] Con 801 y 4001 puntos $I_{max}$ coincide a $<10^{-4}$ relativo. PyFocus usa `scipy.integrate.quad` con error relativo $< 1.5\times10^{-8}$ [Caprile2022, p. 14].
  - Integral 2D: $N_\theta = N_\phi = 200$ alcanza; 10 no alcanza [Caprile2022, pp. 36–37, Fig. 19].
  - Grilla focal: píxel ≤5 nm, como la resolución de simulación de [Lopez2023, p. 12].
  - Para MINFLUX conviene tabular $I(\rho)$ en 1D (por simetría) e interpolar. [propio]
- **Interfaz vidrio/agua:** Eqs. 19–21 de [Caprile2022, pp. 34–35], o el factor empírico 0.70 en $z$ [Gwosch2020, p. 9].

### 6.3 Números de chequeo

1. **Cero exacto.** Vórtice $+1$ con $(\hat x+i\hat y)/\sqrt2$: $I(0,z) = 0$ para todo $z$.
   - Papers: [Caprile2022, p. 18, Fig. 8c; p. 17]; [Lopez2023, p. 11].
   - [propio] $I(0)/I_{max} \sim 10^{-33}$ numérico.
2. **Mano opuesta.** Circular izquierda da máximo central de $E_z$ [Caprile2022, p. 35].
   - [propio] Con $F = 5/3$: $I(0)/I_{max} = 0.87$, polarización lineal da $0.38$, elíptica ($\gamma = 45^\circ$, $\beta = 60^\circ$) da $0.063$.
3. **Diámetro pico-pico** (NA 1.4, $n = 1.5$, $\lambda = 640$ nm, $h = 3$ mm) [Caprile2022, p. 20]:
   - $F = 0.5$: 508 nm;
   - $F = 1$: 428 nm;
   - $F = 2$: 392 nm;
   - uniforme: 380 nm.
   - Mi cálculo difiere en algunos casos; ver §8.
4. **Curvatura en el cero** [propio]: $I/I_{max} \approx 7.0\times10^{-5}\,\rho^2$ ($\rho$ en nm) para $F = 5/3$, y $7.2\times10^{-5}$ para pupila uniforme.
   - Compatible con el modelo de [Tarkowski, Eq. 4] con $d = 360$ nm: $5.8\times10^{-5}$.
   - Compatible con un LG de $w = 272$ nm: $7.35\times10^{-5}$.
5. **FWHM del hueco central** [propio]: ≈189 nm para $F = 5/3$, ≈185 nm para pupila uniforme.
6. **Inclinación.** 0.0009° produce un desplazamiento de 50 nm sin deformación [Caprile2022, p. 23].
7. **Escala axial en 3D con desajuste de índice:** 0.70 [Gwosch2020, p. 9].
8. **Contraste residual de polarización** 1–5 % [Lopez2023, p. 8]. [propio] Da $I(0)/I_{max} \approx 2\times10^{-5}$ a $6\times10^{-4}$ (§8).
9. **Top-hat** [propio]: lóbulos axiales en $\pm 510$ nm; anillo lateral al 25 %; $c_z/c_\rho \approx 2$.

---

## 7. Ambigüedades de notación y posibles erratas

1. **Exponente de [Tarkowski, p. 4, Eq. 4].** Está impreso como $-4\,\ln(2)^3(\ldots)$.
   - Si se lee $-4\ln2$, la Eq. 4 es un LG$_0^1$ con $d$ = FWHM del gaussiano equivalente, y $D_{pp} = 432$ nm para $d = 360$ nm, en línea con los 380–430 nm vectoriales.
   - Si se lee $-4(\ln 2)^3$, daría $D_{pp} \approx 624$ nm, poco realista.
   - Además, $d$ se llama "FWHM de las crestas" [Tarkowski, p. 5], lo que no coincide con ninguna de las dos lecturas. Recomiendo $-4\ln2$ y verificar con los autores o el código (github.com/stefani-lab).
2. **Error en $E_{fy}$ de [Caprile2022, p. 32, Eq. 15]** [propio]. Para circular $(c_x, c_y) = (1, i)/\sqrt2$:
   - La Eq. 15 da $E_y \propto -\tfrac{i}{\sqrt2}[I'_1e^{i\phi} - I'_3 e^{3i\phi} + I'_2 e^{-i\phi}]$.
   - La integración directa de la Eq. 1 da $E_y \propto \tfrac{i}{\sqrt2}[I'_1 e^{i\phi} - I'_3 e^{3i\phi}]$, sin término $I'_2$.
   - Con la Eq. 15, $|E_y|^2$ no es rotacionalmente simétrico: error de hasta 40 % en $|E_y|^2$ a $\rho = 100$ nm.
   - $E_x$ y $E_z$ sí coinciden con la integral directa. El cero central no se ve afectado.
   - Puede ser una errata tipográfica del paper y no del código PyFocus; no revisé el código.
3. En [Caprile2022, p. 32] se dice "using equation 10", pero la identidad usada es la Eq. 12.
4. En la Eq. 17 [Caprile2022, p. 33] (sin propagación) queda el prefactor $k/L$, que no tiene sentido sin propagación. Es irrelevante para la intensidad normalizada.
5. $I = E_{fx}^2+\ldots$ [Caprile2022, p. 5] debe leerse como $|E|^2$.
6. **Significado de $k$.** En la Eq. 1 es el número de onda en el medio ($2\pi n/\lambda$). En las máscaras custom, PyFocus define $k = 2\pi/\lambda$ [Caprile2022, p. 21].
7. **Tipo de radio de $w_0$.** En PyFocus $w_0$ es el radio de **amplitud** $1/e$ [Caprile2022, p. 6], que coincide con el $1/e^2$ de intensidad. En LG, $w$ es el radio $1/e^2$ de intensidad del gaussiano asociado.
8. **Handedness.** "Derecha" en PyFocus es $(\hat x + i\hat y)/\sqrt2$ "desde la fuente" [Caprile2022, p. 12]. [Lopez2023, p. 11] dice "right-handed" para el caso bueno. Eso depende de la convención temporal ($e^{-i\omega t}$ o $e^{+i\omega t}$) y del sentido de observación. Lo robusto es la condición $\ell\sigma = +1$ en la notación de §2.2, que aquí se lee como "misma mano".
9. **"Spherical" de PyFocus.** $Z_2^0 = \sqrt3\,2r^2$ [Caprile2022, p. 25] es desenfoque sin pistón, no aberración esférica primaria ($Z_4^0 = \sqrt5(6r^4-6r^2+1)$).
10. **FWHM del foco regular en [Gwosch2020, p. 9].** Aparece como 360 nm, 300 nm y 330 nm (lateral) en distintos lugares.
11. **Nombre $\alpha$.** Designa a la vez el semiángulo de apertura [Caprile2022, p. 11], el ángulo de inclinación [Caprile2022, p. 23] y el factor de anisotropía $d_z/d_{xy}$ [Tarkowski, p. 13].
12. **Radio del disco top-hat:** no especificado en ninguna fuente asignada.

---

## 8. Verificación rápida (CÁLCULO PROPIO)

Scripts en el scratchpad de la sesión: `verify.py`, `dscan.py`, `tophat.py`, `thchk.py`, `imperf.py`, `signs.py`. Parámetros: NA 1.4, $n = 1.5$, $\lambda = 640$ nm, $h = 3$ mm y, salvo indicación, $w_0 = 5$ mm. La referencia es la integración 2D directa de la Eq. 1 de [Caprile2022] con $\mathbf E_0$ físico (256 puntos en $\phi'$ y 4001 en $\theta$).

**Fórmulas cerradas.**
- Las de §2.2 ($I_A, I_B, I_C$) reproducen la integral 2D componente por componente hasta el redondeo.
- Las de top-hat de §3 también.
- La Eq. 15 impresa en [Caprile2022] reproduce $E_x$ y $E_z$ pero no $E_y$. Ejemplo a $\rho = 100$ nm, $\phi = 0$: $|E_y|^2 = 0.0146$ con la Eq. 15 contra $0.0102$ con la integral.

**Cero.** Con la mano correcta, $I(0)/I_{max} \approx 10^{-33}$, y es nulo sobre el eje a $z = 0$, 200 y 500 nm.

**Otras polarizaciones** ($F = 5/3$), valor de $I(0)/I_{max}$:

| Polarización | $I(0)/I_{max}$ | Composición |
|---|---|---|
| Circular opuesta | 0.869 | todo $|E_z|^2$ |
| Lineal | 0.382 | |
| Elíptica, $\beta = 60^\circ$ | 0.063 | |

**Diámetro pico-pico contra $F$:**

| $F$ | Mío | Caprile p. 20 |
|---|---|---|
| 0.33 | 750 nm | — |
| 0.5 | 516 nm | 508 nm |
| 1.0 | 401 nm | 428 nm |
| 5/3 | 385 nm | — |
| 2.0 | 383 nm | 392 nm |
| uniforme | 378 nm | 380 nm |

Coincide bien en los extremos, pero a $F = 1$ la diferencia es de 27 nm. No encontré un reescaleo único de $w_0$ que reconcilie los tres valores, así que puede ser lectura de píxel o figura en el paper. Conviene tomar ~380–400 nm como chequeo robusto para $F \gtrsim 1$.

**Curvatura y hueco central.**

| $F$ | $I/I_{max}$ cerca del cero | FWHM del hueco central |
|---|---|---|
| 5/3 | $6.96\times10^{-5}\rho^2$ | 189 nm |
| 1 | $6.45\times10^{-5}\rho^2$ | 196 nm |
| uniforme | $7.23\times10^{-5}\rho^2$ | 185 nm |

($\rho$ en nm.)

**Profundidad del cero por polarización imperfecta**, con $\beta = 90^\circ$ y contraste del analizador $C = |\cos 2\gamma|$ [def. de contraste de Lopez2023, p. 8]. Se ajusta a $I_{min}/I_{max} \approx 0.23\,C^2$, con el mínimo siempre en el centro.

| $C$ | $I_{min}/I_{max}$ |
|---|---|
| 0.01 | $2.3\times10^{-5}$ |
| 0.02 | $9.4\times10^{-5}$ |
| 0.05 | $5.9\times10^{-4}$ |
| 0.10 | $2.3\times10^{-3}$ |
| 0.20 | $9.5\times10^{-3}$ |

**Descentración VPP + haz** (Eq. 10 de Caprile), con $D = d/h$. El cero se conserva casi perfecto ($I_{min}/I_{max} \le 10^{-5}$) pero **se desplaza** en el plano focal.

| $D$ | Desplazamiento del cero | $I(0)/I_{max}$ en el centro original |
|---|---|---|
| 0.01 | 1.8 nm | $2.4\times10^{-4}$ |
| 0.03 | 5.5 nm | $2.1\times10^{-3}$ |
| 0.10 | 18 nm | $2.4\times10^{-2}$ |

Si se descentra solo la VPP, el desplazamiento es algo mayor: 2.4, 7.1 y 24 nm. Implicancia: la descentración es sobre todo un **sesgo de posición** del cero, más que una pérdida de profundidad.

**Aberraciones sobre el vórtice** (amplitud $A$ en rad, Zernike normalizados), valor de $I_{min}/I_{max}$:

| Aberración | $A = 0.1$ | $A = 0.3$ | Comentario |
|---|---|---|---|
| Astigmatismo | $5.2\times10^{-3}$ | $4.4\times10^{-2}$ | cero relleno, centrado |
| Coma | $1.2\times10^{-5}$ | $9\times10^{-4}$ | cero desplazado 5 y 16 nm |
| Trébol | $6.7\times10^{-4}$ | $5.5\times10^{-3}$ | |
| Esférica $Z_4^0$ | 0 | 0 | la simetría de revolución preserva el cero en el plano focal; afecta la estructura axial, no evaluada aquí |

En resumen: el astigmatismo es la aberración que más degrada la profundidad; la coma sobre todo desplaza el cero.

**Top-hat vectorial** ($F = 5/3$):
- $\rho_0 = 0.682\,h$ para cero exacto; con pupila uniforme, $0.7145\,h$ contra el $1/\sqrt2 = 0.7071$ escalar.
- Máximos axiales en $\pm 510$ nm.
- Anillo focal en $\rho = 261$ nm al 25 % del máximo axial.
- Curvaturas: $c_z = 1.02\times10^{-5}$ y $c_\rho = 5.1\times10^{-6}$ nm$^{-2}$; cociente ≈2.

**Consistencia de modelos analíticos:**
- LG$_0^1$ con $w = 272$ nm: $D_{pp} = 385$ nm y curvatura $7.35\times10^{-5}$ nm$^{-2}$, cercana a la vectorial ($6.96\times10^{-5}$).
- [Tarkowski, Eq. 4] con $d = 360$ nm y exponente $-4\ln2$: $D_{pp} = 432$ nm y curvatura $5.8\times10^{-5}$ nm$^{-2}$. Es una dona algo más ancha que la de NA 1.4 con pupila llena.
