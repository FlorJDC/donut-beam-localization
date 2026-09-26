# A. Teoría de MINFLUX y localización con mínimos de luz: nota para simulación

**Fuentes leídas** (PDFs en `C:\Users\BANGHO\Documents\Doctorado\Papers`, extraídas con PyMuPDF):

| Clave | Archivo | Páginas del PDF |
|---|---|---|
| [Balzarotti2017] | `Balzarotti et al. - 2017 - Science - MINFLUX.pdf` | 79 (pp. 1–7 artículo; pp. 8–75 Material Suplementario (SM); pp. 76–79 referencias del SM) |
| [MasulloSI2022] | `Masullo et al. - 2022 - A common framework ... - SI.pdf` | 10 (**solo la Información Suplementaria**; el texto principal no está en la carpeta) |
| [RASTMIN2022] | `paper_rastmin.pdf` (Masullo et al., Light Sci. Appl. 11:199, 2022) | 9 |
| [MS2022] | `Masullo y Stefani - 2022 - Multiphoton ... light minima.pdf` (News & Views) | 4 |
| [Stefani2023] | `Stefani - 2023 - Tracking nanoscopic motion with minima of light.pdf` (News & Views) | 2 |

**Convención de citas.** "p. N" es siempre la página **del PDF** (1-based), no el número impreso. En el SM de Balzarotti, página del PDF = número impreso del SM + 7 (p. ej., la página impresa "4" del SM es la p. 11 del PDF). Las ecuaciones del SM se citan como "Eq. S*n*".
Los resultados marcados **[derivación propia]** no aparecen en los papers: los derivé yo y los verifiqué numéricamente (script en el scratchpad, basado en Eq. S11/S26 con diferencias finitas). Hay que tratarlos como afirmaciones que se deben volver a verificar.

---

## 1. Modelo probabilístico

### 1.1 Exposiciones, intensidad y tasa de emisión
- El emisor, ubicado en $\bar r_m\in\mathbb R^d$, se expone secuencialmente a $K$ intensidades $\{I_0(\bar r),\dots,I_{K-1}(\bar r)\}$ y se registran los conteos $\bar n=\{n_0,\dots,n_{K-1}\}$ [Balzarotti2017, p. 11].
- Sin saturación, la media de Poisson de cada exposición es lineal en la intensidad local:
  $\lambda_i = c_e\, q_e\, \sigma_a\, I_i(\bar r_m)$, donde $c_e$ es la eficiencia de colección, $q_e$ el rendimiento cuántico y $\sigma_a$ la sección eficaz de absorción [Balzarotti2017, p. 11, Eq. S1].
- Todas las exposiciones usan el mismo haz desplazado: $I_i(\bar r)=I(\bar r-\bar r_{b_i})$ [Balzarotti2017, p. 14, Eq. S15]; es la misma forma que usa [MasulloSI2022, p. 3, Eq. S1]: $I(\mathbf r_E-\mathbf r_i)$.
- En el caso multifotónico (o de otra no linealidad), $\lambda_i \propto I^c(\mathbf r-\mathbf r_i)$, con $c=2$ para dos fotones y $c=3$ para tres; también se admiten exponentes fraccionarios [MS2022, p. 2].
- En 1D, el artículo principal escribe la fluorescencia como $f(x)=C\,I(x_m-x)$. $C$ agrupa brillo, detección y orientación molecular, y se cancela en la solución, de modo que la orientación del dipolo no sesga la estimación [Balzarotti2017, p. 2].

### 1.2 Estadística: Poisson → multinomial
- Cada $n_i\sim\mathrm{Poisson}(\lambda_i)$. Condicionando en el total $N=n_0+\dots+n_{K-1}$ [Balzarotti2017, p. 11, Eq. S2] se obtiene una multinomial:
  $$P(\bar n|N)=\frac{N!}{n_0!\cdots n_{K-1}!}\prod_{i=0}^{K-1}p_i^{\,n_i}$$ [Balzarotti2017, p. 11, Eq. S3]
- Sin fondo:
  $$p_i^{(0)}(\bar r_m)=\frac{\lambda_i}{\sum_j\lambda_j}\approx\frac{I_i(\bar r_m)}{\sum_j I_j(\bar r_m)}$$ [Balzarotti2017, p. 11, Eq. S4]. El brillo molecular se cancela. La expresión equivalente en [MasulloSI2022] está en la p. 3, Eq. S1.
- Solo $K-1$ de los $p_i$ son independientes (espacio $\bar p$ "reducido"), porque $p_{K-1}=1-\sum_{j\le K-2}p_j$ [Balzarotti2017, p. 12, Eq. S5]. El SM advierte que ignorar esta restricción da una matriz de Fisher incorrecta [Balzarotti2017, p. 12].
- En 1D con $K=2$ la distribución es binomial: $p_0=f_0/(f_0+f_1)=I_0/(I_0+I_1)$ [Balzarotti2017, p. 2].

### 1.3 Fondo y SBR (definición exacta)
- El fondo se modela como Poisson con media $\lambda_{b_i}$ por exposición:
  $$p_i=\frac{\lambda_i+\lambda_{b_i}}{\sum_j(\lambda_j+\lambda_{b_j})}$$ [Balzarotti2017, p. 19, Eq. S28]; [MasulloSI2022, p. 3, Eq. S2].
- **SBR** = señal total sobre fondo total, sumados sobre las $K$ exposiciones:
  $$\mathrm{SBR}(\bar r_m)=\frac{\sum_{j=0}^{K-1}\lambda_j}{\sum_{j=0}^{K-1}\lambda_{b_j}}$$ [Balzarotti2017, p. 19, Eq. S29].
  Con fondo igual en todas las exposiciones ($\lambda_{b_i}=\lambda_b$) queda $\mathrm{SBR}=\sum_j\lambda_j/(K\lambda_b)$ [MasulloSI2022, p. 3, Eq. S3; p. 4]. La SBR **depende de la posición del emisor y del patrón**, no solo del sistema [Balzarotti2017, p. 19].
- Forma que conviene implementar:
  $$p_i(\bar r_m)=\frac{\mathrm{SBR}}{\mathrm{SBR}+1}\,p_i^{(0)}(\bar r_m)+\frac{1}{\mathrm{SBR}+1}\,\frac1K$$ [Balzarotti2017, p. 19, Eq. S30]; [MasulloSI2022, p. 4, Eq. S6].
  Se supone que el fondo depende solo de la potencia de excitación y es igual en todas las exposiciones [Balzarotti2017, p. 19].
- Si $\lambda_b$ es fijo y se achica $L$, el emisor ubicado en el centro recibe menos intensidad total y la SBR cae. Para la dona de Eq. S17 [Balzarotti2017, p. 20, Eq. S32]:
  $$\mathrm{SBR}(\bar0,L)=\frac{L^2}{L_0^2}\exp\!\Big(\frac{\ln 2}{\mathrm{fwhm}^2}(L_0^2-L^2)\Big)\,\mathrm{SBR}(\bar0,L_0)$$
- Estimación experimental: $\mathrm{SBR}=N/(\lambda_{b}Q)-1$, donde $Q$ es el número de ciclos de multiplexado y $\lambda_b$ el fondo medio por ciclo, con $K=4$ [Balzarotti2017, p. 42, Eq. S76].

---

## 2. Patrones de excitación

### 2.1 Perfiles de haz [Balzarotti2017, p. 14]
- Cuadrático (ideal, "físicamente imposible"): $I_{quad}=A_{quad}\,r^2$ (Eq. S16).
- **Dona** (sección radial de un LG$_{01}$, normalizada a pico $A_0$):
  $$I_{donut}(r)=A_0\,\frac{4e\ln2\,r^2}{\mathrm{fwhm}^2}\,e^{-4\ln2\,r^2/\mathrm{fwhm}^2}$$ (Eq. S17).
  "fwhm" es un **parámetro de tamaño**, no el FWHM del anillo. El diámetro entre los picos del anillo es $\mathrm{fwhm}/\sqrt{\ln 2}\approx1.2\,\mathrm{fwhm}$ [p. 14]. Verificación propia: el máximo está en $r=\mathrm{fwhm}/(2\sqrt{\ln2})$ y vale exactamente $A_0$.
- Onda estacionaria: $I_{sw}=A_i\sin^2(\bar k\cdot\bar r)$ (Eq. S18). Gaussiano: $I_{gauss}=A_0e^{-4\ln2\,r^2/\mathrm{fwhm}^2}$ (Eq. S19).
- Aproximación cuadrática de la dona: $A_{quad}=4e\ln2\,A_0/\mathrm{fwhm}^2$. Es válida si $4\ln2\,r^2\ll\mathrm{fwhm}^2$ (Eq. S20) y si $L\ll\mathrm{fwhm}/\sqrt{\ln2}$ [Balzarotti2017, pp. 14–16].
- Profundidad del cero en el experimento: el mínimo de la dona era $<0.2\%$ del máximo del anillo [Balzarotti2017, p. 6] (fig. S7). En los papers leídos **no** hay un modelo analítico de cero imperfecto. Si se quiere simular, una opción es sumar un término $\epsilon A_0$ constante o gaussiano; eso es una propuesta, no algo tomado de estas fuentes.
- En el experimento se ajustó la zona central de la PSF medida (región de 200×200 nm², o 300×300 nm² si $L>100$ nm) con un polinomio 2D de orden 4 (Eq. S72), y los $p_i$ se construyeron con esa PSF medida, restándole el offset [Balzarotti2017, p. 38, Eq. S72–S73].

### 2.2 1D con dos exposiciones
Los ceros están en $x_{b0}=-L/2$ y $x_{b1}=+L/2$. **$L$ es la separación total** (el rango de sondeo es $-L/2<x<L/2$) [Balzarotti2017, p. 2; p. 15].

### 2.3 TCP 2D con K = 4 (patrón de Balzarotti)
[Balzarotti2017, p. 17, Eq. S24]:
$$\bar r_{b0}=[0,0]^T,\qquad \bar r_{b_i}=\frac L2[\cos\alpha_i,\sin\alpha_i]^T,\quad \alpha_i=i\,\frac{2\pi}{3},\ i=1,2,3$$
- **$L$ es el DIÁMETRO** del círculo que contiene los tres ceros periféricos. El radio es $L/2$ y el lado del triángulo es $L\sqrt3/2$ (lo último por geometría; no está en el texto).
- Con esa indexación, $\alpha_3=2\pi$, así que la exposición 3 queda sobre el eje $+x$, la 1 en 120° y la 2 en 240°. Esto importa para la fórmula LMS (Eq. S49).
- Hacen falta al menos 3 exposiciones en 2D por la restricción de Eq. S2 [Balzarotti2017, p. 17]. La cuarta exposición central elimina máximos múltiples de la verosimilitud porque aporta información radial. En la fig. S2, $\bar n=[6,8,12]$ con 3 donas (fwhm = 200 nm) da una verosimilitud mal comportada; con $\bar n=[0,6,8,12]$ aparece un máximo único. Un conteo de **cero** fotones en la exposición central es informativo [Balzarotti2017, pp. 17–18; p. 54].
- El artículo principal llama "campo de visión" (FOV) a la región de diámetro ~$L$ [Balzarotti2017, p. 2].

### 2.4 Otros patrones (marco común SML-SSI)
- [MasulloSI2022, p. 9] compara MINFLUX, OTMIN, RASTMIN, RASTMAX, OT/MINSTED y SMCT (tabla en §8). El texto principal, que define con precisión las geometrías de OTMIN y OT, **no está disponible**.
- RASTMIN: barrido raster cuadrado de lado $L$ con $K=6\times6$ píxeles. Los autores recomiendan $K=16$–$100$ y reportan que la precisión mejora poco más allá de $K\approx36$ [RASTMIN2022, pp. 2–3].
- Según [MS2022, p. 2] y [Stefani2023, p. 2], distintos arreglos de mínimos (polígonos, franjas, raster) dan una eficiencia fotónica prácticamente igual.

---

## 3. Estimadores

### 3.1 MLE
- $\hat{\bar r}_m^{MLE}=\arg\max\mathcal L(\bar r|\bar n)$ con $\mathcal L=P(\bar n|N,\bar r)$ [Balzarotti2017, p. 21, Eq. S33]. Basta maximizar $\ell(\bar r|\bar n)=\sum_i n_i\ln p_i$ [p. 22, Eq. S37].
- 1D con dos parábolas:
  $$p_0=\frac{(1+2x/L)^2}{2(1+4x^2/L^2)},\quad p_1=\frac{(1-2x/L)^2}{2(1+4x^2/L^2)}$$ [p. 21, Eq. S34].
  Solución cerrada, eligiendo la rama $-L/2<x<L/2$:
  $$\hat x_m^{MLE}=-\frac L2+\frac{L}{1+\sqrt{n_1/n_0}}$$ [p. 22, Eq. S36]. El artículo principal da la misma fórmula [p. 2].
  La otra raíz, $-L/2+L/(1-\sqrt{n_1/n_0})$, queda fuera del rango de sondeo.
- 2D con 4 donas: no hay forma cerrada. Se usan búsquedas en grilla sucesivas [p. 22]. En el experimento hubo 4 grillas con pasos de 5, 1, 0.1 y 0.01 nm; la primera cubría un diámetro de 240 nm [p. 42].
- **Convergencia al CRB**: en el origen, con SBR = 10, el MLE alcanza el CRB desde $N\approx100$ para todos los $L$ mostrados. Lejos del centro necesita más fotones: en $x=50$ nm con $L=75$ nm, recién desde $N\approx500$ [Balzarotti2017, pp. 22–23; p. 29: "$N\gtrsim100$–$500$"]. El artículo principal dice $N=\sum n_i\gtrsim100$ [p. 3].
- 1D con dos gaussianas: $\hat x^{MLE}=\frac{\mathrm{fwhm}^2}{8\ln2\,L}\ln(n_0/n_1)$, con la convención $\hat x\equiv0$ si algún $n_i=0$ [p. 23, Eq. S40]. Fuera del centro, el MLE sesgado puede quedar **por debajo** del CRB [p. 24].

### 3.2 LMS linealizado [Balzarotti2017, pp. 25–26]
- Se linealiza en el origen: $p_i(\bar r)\approx p_i(\bar 0)+\sum_j r_j\,\partial p_i/\partial r_j$ (Eq. S43). Con $\hat p_i=n_i/N$ (MLE de la multinomial, Eq. S45) se resuelve por mínimos cuadrados:
  $$\hat{\bar r}_{LMS}=(\mathcal J^T\mathcal J)^{-1}\mathcal J^T(\hat{\bar p}-\bar p(\bar0))$$ (Eq. S48), donde $\mathcal J$ es el Jacobiano completo ($K\times d$), no el reducido.
- Para el TCP con dona (Eq. S49–S50):
  $$\hat{\bar r}_{LMS}=-\frac{1}{1-L^2\ln2/\mathrm{fwhm}^2}\sum_{i=1}^{3}\hat p_i\,\bar r_{b_i}
  =\frac L2\,\frac{1}{1-L^2\ln 2/\mathrm{fwhm}^2}\begin{bmatrix}-\hat p_3+\tfrac12(\hat p_1+\hat p_2)\\ \tfrac{\sqrt3}{2}(\hat p_2-\hat p_1)\end{bmatrix}$$
  Nota de chequeo: con $\alpha_3=0$, la componente $x$ de $-\sum\hat p_i\bar r_{b_i}$ es $\frac L2[\frac12(\hat p_1+\hat p_2)-\hat p_3]$, lo que coincide con Eq. S49. La componente $y$ también coincide. Para la cuadrática ($\mathrm{fwhm}\to\infty$) verifiqué que $(\mathcal J^T\mathcal J)^{-1}\mathcal J^T$ da exactamente $-\sum\hat p_i\bar r_{b_i}$ **[derivación propia]**.
- El LMS no usa $n_0$: $p_0$ no tiene término lineal en el origen [p. 26].
- **Fondo [derivación propia]**: la Eq. S43 se construye con el $\bar p$ de Eq. S4, o sea sin fondo. Si se usa el $p_i$ de Eq. S30, el Jacobiano se multiplica por $s=\mathrm{SBR}/(\mathrm{SBR}+1)$ y el término constante se cancela en $\mathcal J^T(\cdot)$, porque $\sum_i\nabla p_i^{(0)}=0$. El LMS con fondo es entonces el de Eq. S50 dividido por $s$. Si no se hace esta corrección, la estimación queda contraída hacia el centro en un factor $s$.

### 3.3 mLMS (LMS modificado) [Balzarotti2017, p. 26, Eq. S51]
$$\hat{\bar r}^{(k)}_{mLMS}(\hat{\bar p},\bar\beta)=-\frac{1}{1-L^2\ln2/\mathrm{fwhm}^2}\Big(\sum_{j=0}^{k}\beta_j\hat p_0^{\,j}\Big)\sum_{i=1}^3\hat p_i\,\bar r_{b_i}$$
- Incorpora información radial a través de $\hat p_0$. En el tracking en vivo (FPGA) se usó $k=1$ con **$\beta_0=1.27$, $\beta_1=3.8$** [Balzarotti2017, p. 35; p. 44; p. 73].
- **Es sesgado**. Su media se puede escribir en forma cerrada (Eq. S53). Para $k=2$ [p. 27, Eq. S54]:
  $$\bar R^{(2)}=\frac{-1}{1-L^2\ln2/\mathrm{fwhm}^2}\sum_{i=1}^3\Big(\beta_0p_i+\tfrac{N-1}{N}\big[\beta_1+\tfrac{\beta_2}{N}\big]p_0p_i+\beta_2\tfrac{(N-1)(N-2)}{N^2}p_0^2p_i\Big)\bar r_{b_i}$$
  Esta expresión usa $E(\hat p_0^j\hat p_i)=E(n_0^jn_i)/N^{j+1}$ y los momentos factoriales de la multinomial [p. 27].
- **numLMS** (el que se usó en el post-procesamiento): se optimiza $\bar\beta$ para minimizar el sesgo medio dentro de un radio $\mathcal R_1$ (con $\beta_1=0$ por defecto), y después se invierte numéricamente $\bar R^{(2)}$ con un interpolante $\mathcal F$ [pp. 28–29, Eq. S55–S58]. Parámetros usados: $\mathcal R_2=L/2$, grilla de 1 nm, $M=10^4$. El estimador funciona en $|\bar r|\lesssim L/2$ y se probó con SBR = 2–60 [p. 29]. Hace falta conocer $N$ y la SBR [p. 29].
- Para $N<100$ (tracking) conviene el numLMS en lugar del MLE [Balzarotti2017, p. 3; p. 4, leyenda de la Fig. 3E].

---

## 4. Información de Fisher y CRB

### 4.1 Fórmula general (patrones arbitrarios)
- $F_{\bar p}$ sobre el espacio reducido: $\{F_{\bar p}\}_{ij}=N\big(1/p_{K-1}+\delta_{ij}/p_i\big)$ con $i,j\in[0,K-2]$ [Balzarotti2017, p. 12, Eq. S8]. Se reparametriza como $F_{\bar r}=\mathcal J^{*T}F_{\bar p}\mathcal J^*$ [p. 12, Eq. S9; p. 13, Eq. S10].
- Forma simplificada, válida para **cualquier** patrón $\{I_i\}$, y que es la que conviene implementar:
  $$F_{\bar r_m}=N\sum_{i=0}^{K-1}\frac{1}{p_i}\,\nabla p_i\,\nabla p_i^{T}$$ [Balzarotti2017, p. 13, Eq. S11].
- Cota: $\Sigma(\bar r_m)\ge\Sigma_{CRB}=F^{-1}_{\bar r_m}$ [p. 13, Eq. S12]. Métricas:
  $$\tilde\sigma_{CRB}=\sqrt{\tfrac1d\,\mathrm{tr}(\Sigma_{CRB})}$$ [p. 13, Eq. S13; p. 74, Tabla S1] e isotropía $\mathbb I=\min\sigma_i/\max\sigma_i$ [p. 13, Eq. S14].
  El texto llama a $\tilde\sigma_{CRB}$ "media aritmética de los autovalores". En rigor es la raíz de la media de los autovalores $\sigma_i^2$, es decir, un σ RMS por eje; la convención del repo coincide.
- Forma explícita en 2D: [p. 18, Eq. S25–S26].
- **El fondo entra solo a través de $p_i$** (Eq. S30) dentro de Eq. S11. El marco de Masullo usa exactamente la misma construcción con $p_i$ de Eq. S6 [MasulloSI2022, p. 4]; el CRB completo está en el texto principal, que no está disponible.

### 4.2 1D, dos exposiciones [Balzarotti2017, pp. 14–16]
- $\tilde\sigma_{CRB}=\frac{1}{\sqrt N}\frac{\sqrt{p_0(1-p_0)}}{|dp_0/dx|}$ (Eq. S21).
- Cuadrática: $\tilde\sigma_{CRB}(x)=\frac{L}{4\sqrt N}\big[1+(2x/L)^2\big]$ (Eq. S22b), con **$\tilde\sigma_{CRB}(0)=L/(4\sqrt N)$** (Eq. S22c; texto principal p. 2).
- Dona: $\tilde\sigma_{CRB}(0)=\frac{L}{4\sqrt N}\big(1-\ln2\,L^2/\mathrm{fwhm}^2\big)^{-1}$ (Eq. S22f). La expresión general en $x$ es Eq. S22e.
- Onda estacionaria: $\tilde\sigma(0)=\frac{1}{2k\sqrt N}\tan\frac{kL}{2}$ (Eq. S22i). Gaussiana: $\tilde\sigma(0)=\frac{\mathrm{fwhm}^2}{4\ln2\,L\sqrt N}$, que escala como $L^{-1}$ (Eq. S23c).

### 4.3 TCP 2D, K = 4, forma cerrada en el centro
- **Sin fondo** [Balzarotti2017, p. 18, Eq. S27]:
  $$\tilde\sigma_{CRB}(\bar0)=\frac{L}{2\sqrt{2N}}\Big(1-\frac{L^2\ln2}{\mathrm{fwhm}^2}\Big)^{-1}$$
- **Con fondo** [Balzarotti2017, p. 20, Eq. S31]:
  $$\tilde\sigma_{CRB}(\bar0)=\frac{L}{2\sqrt{2N}}\Big(1-\frac{L^2\ln2}{\mathrm{fwhm}^2}\Big)^{-1}\sqrt{\Big(1+\frac{1}{\mathrm{SBR}(\bar0)}\Big)\Big(1+\frac{3}{4\,\mathrm{SBR}(\bar0)}\Big)}$$
  El factor $3/4$ viene de $K=4$ con 3 exposiciones periféricas. Verifiqué que Eq. S31 coincide con el cálculo numérico de Eq. S11 + Eq. S30 hasta $10^{-9}$ **[derivación propia]**.
- En el centro el CRB es isótropo: $\sigma_x=\sigma_y=\tilde\sigma_{CRB}$, verificado numéricamente. Fuera del centro la elipse deja de ser isótropa [p. 19; fig. S3B, p. 55].
- **Dependencias**: lineal en $L$ para $L\ll\mathrm{fwhm}/\sqrt{\ln2}$ [p. 18]; proporcional a $N^{-1/2}$; y con fondo fijo existe un $L$ óptimo, porque la SBR cae con $L$ según Eq. S32 y por debajo de ese límite achicar $L$ empeora la precisión [p. 20; fig. S3D, pp. 55–56]. [MS2022, p. 3] y [RASTMIN2022, p. 2] repiten que la mejora con $L$ está limitada por la SBR.
- Relación 1D/2D (cuadrática, sin fondo): $L/(4\sqrt N)$ en 1D frente a $L/(2\sqrt2\sqrt N)$ en 2D. Por eje, el 2D es $\sqrt2$ peor para el mismo $L$ y el mismo $N$ total (comparación entre Eq. S22c y Eq. S27).
- **Sutileza numérica [derivación propia]**: con fondo nulo, $p_0(\bar0)=0$ y $\nabla p_0(\bar0)=0$, así que el término $(\nabla p_0)^2/p_0$ de Eq. S11 es $0/0$ en el origen. En el límite $r\to0$ ese término vale $\frac{16N}{3L^2}\hat r\hat r^T$ (cuadrática), que no es cero y depende de la dirección de aproximación. Eq. S27 corresponde a tomarlo como 0 exactamente en $\bar r=\bar0$. Ejemplo con $L=100$, $N=100$: $\tilde\sigma(\bar0)=3.536$ nm, pero en $x=10^{-3}$ nm se obtiene $\sigma_x=2.739$ nm. Con SBR finita el problema desaparece, porque $p_0>0$ y el término vale 0 de forma continua. **Recomendación: evaluar el CRB en el origen con SBR finita o excluir los $p_i=0$ de forma explícita y documentada.**

### 4.4 Comparación con localización por cámara
- Regla general: $\sigma\approx\sigma_{PSF}/\sqrt N$ con $\sigma_{PSF}\approx100$ nm, así que $N=400$ da $\sigma\approx5$ nm [Balzarotti2017, p. 1]. [Stefani2023, p. 2] escribe $\sigma\propto\lambda/\sqrt N$: para $\sigma=1$–$2$ nm hacen falta $N=2500$–$40000$ según la SBR. Con mínimos, $\sigma\propto L/\sqrt N$ con $L$ sub-difracción [Stefani2023, p. 2].
- Modelo de cámara ideal (sin ruido de lectura ni exceso EM): PSF gaussiana con $\sigma_{PSF}$ [Balzarotti2017, p. 30, Eq. S59]; $p_i$ por píxel integrado con funciones erf [p. 30, Eq. S60]; fondo mediante $\mathrm{SBR}_c=\lambda_{signal}/(\lambda_{bg}/K)$ [p. 31, Eq. S61] y
  $$p_i=\frac{1}{K+\mathrm{SBR}_c}+\frac{\mathrm{SBR}_c}{K+\mathrm{SBR}_c}p_i^{(0)}$$ [p. 31, Eq. S62–S63]. El CRB sale de Eq. S11 y Eq. S13.
- Parámetros usados: $\sigma_{PSF}=100$ nm (Figs. 3–4) u 87 nm (Fig. 5); píxel $a=100$ nm; ROI de 9×9 ($K=81$); $\mathrm{SBR}_c$ = 50 o 500 (Fig. 3), 500 (Fig. 4), 20 o 40 (Fig. 5) [Balzarotti2017, p. 32].
- Masullo usa la misma construcción con su SBR "total/total" y $\mathrm{FWHM}\approx2.35\,\sigma_{PSF}$ [MasulloSI2022, pp. 5–6, Eq. S7–S10]. Muestra que RASTMAX (raster con un máximo gaussiano) y la cámara ideal tienen la misma forma matemática y precisión prácticamente idéntica ($N=500$, SBR = 5, FWHM = 300 nm) [MasulloSI2022, p. 6; p. 8]. Las cámaras reales son 2–3 veces peores [p. 6].

---

## 5. MINFLUX iterativo y presupuesto de fotones
- [Balzarotti2017] **no** implementa ni cuantifica el esquema iterativo. Lo propone cualitativamente: empezar al límite de difracción y reducir $L$ ("zoom") a medida que se gana información, lo que también reduce las anisotropías que aparecen con $L$ grande. Basta corregir aberraciones en el último paso, el de $L$ más chico [Balzarotti2017, pp. 6–7]. Achicar $L$ **no** equivale a usar un prior bayesiano [p. 6].
- La base del zoom es la invariancia de escala: con mínimos parabólicos, $p(x)$ es invariante ante $L$ (solo depende de $x/L$), mientras que con máximos se vuelve cada vez menos sensible [Stefani2023, p. 1, Fig. 1; p. 2]. Eq. S34 lo muestra explícitamente: depende solo de $2x/L$.
- El compromiso es entre precisión en el centro y FOV: el CRB crece como $[1+(2x/L)^2]$ fuera del centro [Balzarotti2017, p. 15, Eq. S22b; p. 16].
- **Leyes de escala con el presupuesto total de fotones en esquemas iterativos: no aparecen en ninguno de los 5 PDFs.** Lo único disponible es $\sigma\propto L/\sqrt N$ por paso, más la cota que impone la SBR (Eq. S31–S32). Cualquier ley iterativa, por ejemplo elegir $L_{k+1}\propto\sigma_k$, habrá que derivarla y simularla en el proyecto.
- Tracking experimental (referencia): frecuencia de localización de 8 kHz; ~9 fotones por localización en promedio y precisión media <48 nm; $L=130$ nm; ~5800 fotones y 742 localizaciones por traza [Balzarotti2017, pp. 4–5]. El texto extraído dice "$\Delta t$ de 125 ms", pero 8 kHz corresponde a 125 µs, así que probablemente se perdió el símbolo µ al extraer el texto. Hay que verificarlo en el PDF.
- Kinesina (News & Views): 2–3 nm con resolución sub-ms usando unas pocas decenas de fotones [Stefani2023, p. 2]. No da parámetros de simulación.

---

## 6. Mínimos de orden superior / multifotón
- Modelo: $\lambda_i\propto I^c(\mathbf r_E-\mathbf r_i)$ [MS2022, p. 2]. Cerca de un cero cuadrático, $\lambda_i\propto|\mathbf r-\mathbf r_i|^{2c}$.
- Resultado citado: 2p-MINFLUX mejora ~2 veces la precisión en el centro, equivalente a necesitar 1/4 de los fotones, dado que $\sigma\propto1/\sqrt N$ [MS2022, p. 2] (el trabajo original es Zhao et al., eLight 2022, que **no** está en la carpeta). Con multifotón la precisión es más heterogénea en el FOV de MINFLUX; RASTMIN lo mitiga [MS2022, p. 2; Fig. 1 en p. 3: $L=100$ nm, $N=500$, SBR = 4, $\lambda$ = 647/800/1300 nm para 1p/2p/3p].
- **[derivación propia, verificada numéricamente]**: con un cero de orden $2c$ (ya sea $c$-fotón sobre un cero cuadrático, o un cero intrínseco $\propto r^{2c}$) y sin fondo:
  - 1D, 2 exposiciones: $\sigma(0)=\dfrac{L}{4c\sqrt N}$; con $c=1$ se recupera Eq. S22c.
  - 2D TCP, K = 4: $\tilde\sigma(\bar0)=\dfrac{L}{2c\sqrt{2N}}$; con $c=1$ se recupera Eq. S27 en el límite cuadrático. Con fondo aparece el mismo factor $\sqrt{(1+1/\mathrm{SBR})(1+3/(4\mathrm{SBR}))}$, pero la SBR en el centro cambia con $c$.
  - Numéricamente ($L=100$, $N=100$): 3.536, 1.768 y 1.179 nm para $c$ = 1, 2 y 3, lo que coincide con el "~2×" de [MS2022].
  - Para $c>1$ el término central $(\nabla p_0)^2/p_0\propto r^{2c-2}\to0$, así que la sutileza de §4.3 solo existe para $c=1$.
- [MS2022, p. 3] señala otros dos puntos: la no linealidad podría compensar la caída de SBR al reducir $L$, y la fotofísica (fotoconmutación cis/trans) permitiría una no linealidad efectiva con dosis bajas.

---

## 7. Valores típicos de parámetros

| Parámetro | Valor | Fuente |
|---|---|---|
| $L$ (nanoscopía, origami) | 70 nm (origami de 11 nm), 50 nm (origami de 6 nm) | [Balzarotti2017, p. 3] |
| $L$ (mapa de precisión) | 100 nm; grilla de 35×35 píxeles cada 3 nm; ≲2 cuentas/píxel | [Balzarotti2017, p. 3; p. 38] |
| $L$ (tracking) | 130 nm | [Balzarotti2017, p. 6 (leyenda Fig. 5H); p. 44] |
| $L$ en figuras teóricas | 25/50/150 nm (1D, fig. S1) | [Balzarotti2017, p. 53] |
| $N$ | ≥500 y ≥1000 (umbrales de origami); 500 para 2 nm con $L=100$; ~9 en tracking | [Balzarotti2017, pp. 3, 5] |
| SBR | 13.6 en el píxel central (Fig. 3); 10 (simulaciones de fig. S9); 3.7 (simulación de tracking); ⟨SBR⟩ ≈ 2.7 en E. coli; 17.75 en un evento de origami | [Balzarotti2017, p. 4; p. 63; p. 35; p. 46; p. 67] |
| SBR (Masullo/RASTMIN) | 5 (tablas); 4 (RASTMIN, multifotón) | [MasulloSI2022, p. 9]; [RASTMIN2022, pp. 2, 4]; [MS2022, p. 3] |
| Parámetro "fwhm" de la dona (Eq. S17) | 300 nm (fig. S1), 360 nm (fig. S3, 2D), 200 nm (fig. S2), 450 nm (simulación de tracking) | [Balzarotti2017, pp. 53, 55, 54, 35] |
| FWHM de la dona (RASTMIN) | 300 nm (**definición de FWHM no explicitada en lo disponible**) | [RASTMIN2022, p. 2] |
| $\lambda$ de excitación | 642 nm (Alexa 647 / ATTO 647N) y 560 nm; activación a 405 nm | [Balzarotti2017, p. 3; p. 47] |
| $\lambda$ (RASTMIN) | 640 nm pulsado; objetivo 100×/1.4 NA de aceite | [RASTMIN2022, p. 4; p. 7] |
| NA (Balzarotti) | "objetivo de inmersión en aceite"; **NA no indicada en el texto** | [Balzarotti2017, p. 48] |
| Profundidad del cero | <0.2% del máximo del anillo | [Balzarotti2017, p. 6] |
| Pinhole | 420 nm (imaging, APD1) y 2.5 µm (tracking, APD2) | [Balzarotti2017, p. 48] |
| Cámara ideal | $\sigma_{PSF}$ = 100 o 87 nm; $a$ = 100 nm; 9×9 píxeles | [Balzarotti2017, p. 32] |
| Deriva / estabilización (RASTMIN) | $\sigma_{DC}$ = 0.8–1.3 nm lateral; $\sigma_{TOT}=\sqrt{\sigma_{DC}^2+\sigma_{CRB}^2}$ | [RASTMIN2022, p. 4] |

---

## 8. Para la simulación

### 8.1 Qué implementar
1. **Haces**: `donut(r; fwhm)` según Eq. S17 [p. 14] y `quadratic(r)` según Eq. S16. Opcionalmente, gaussiana (Eq. S19) y onda estacionaria (Eq. S18). Agregar el exponente $c$ para multifotón ($\lambda\propto I^c$) [MS2022, p. 2].
2. **Patrón**: TCP según Eq. S24 [p. 17], con $L$ = diámetro e índice 0 en el centro. Conviene que sea genérico (lista de $\bar r_{b_i}$) para cubrir raster/RASTMIN y órbitas.
3. **Probabilidades**: $p_i^{(0)}$ de Eq. S4 y $p_i$ con fondo de Eq. S30 [p. 19]. La SBR debe poder ser (a) fija o (b) derivada de un $\lambda_b$ fijo, en cuyo caso es dependiente de la posición y de $L$ según Eq. S29/S32.
4. **Generador de datos**: (a) multinomial($N$, $\bar p$) a $N$ fijo (Eq. S3), o (b) Poisson independientes con $\lambda_i+\lambda_b$ ($N$ aleatorio), que es más fiel al experimento.
5. **Fisher/CRB**: Eq. S11 con gradientes analíticos o por diferencias finitas centradas, Eq. S12 y $\tilde\sigma=\sqrt{\mathrm{tr}\Sigma/d}$ (Eq. S13). Manejar el caso $p_i=0$ (§4.3).
6. **Fórmulas cerradas para tests**: Eq. S22c, S22f, S23c, S27, S31 y S32.
7. **Estimadores**: MLE por grilla gruesa-fina (Eq. S37; pasos 5→1→0.1→0.01 nm, diámetro inicial de 240 nm [p. 42]); MLE 1D cerrado (Eq. S36); LMS (Eq. S50, con el factor $1/s$ si hay fondo, §3.2); mLMS $k=1$ con $\beta=(1.27, 3.8)$ (Eq. S51); y opcionalmente numLMS (Eq. S54–S58).
8. **Métricas Monte Carlo**: sesgo (Eq. S41), desviación estándar (Eq. S42), RMSE; comparar con el CRB.
9. **Cámara ideal** (referencia): Eq. S59–S63, con cuidado con la definición de $\mathrm{SBR}_c$ (§9).
10. **Iterativo**: no hay receta en las fuentes (§5); hay que diseñarlo y documentarlo como aporte propio.

### 8.2 Checks numéricos de validación
(Los números con [cálculo propio] los calculé con las fórmulas citadas; los demás están impresos en los papers.)

1. **1D cuadrática, $L=50$ nm, $N=100$**: $\sigma_{CRB}(0)=1.25$ nm; la Fig. 1D lo redondea a "1.3 nm" [Balzarotti2017, p. 2, Eq. S22c en p. 15]. En los bordes $x=\pm L/2$: 2.5 nm, que coincide con el "≤2.5 nm" del texto principal [p. 2].
2. **2D TCP cuadrática sin fondo, $L=100$, $N=100$**: $\tilde\sigma(\bar0)=L/(2\sqrt{2N})=3.536$ nm [Eq. S27, p. 18; cálculo propio]. Con dona de fwhm = 360 nm: 3.735 nm. Con $L=50$ y fwhm = 360: 1.792 nm [cálculo propio con Eq. S27].
3. **Con fondo (Eq. S31)**, $L=100$, $N=500$, SBR = 10, fwhm = 360: 1.817 nm [cálculo propio; Eq. S31 en p. 20].
4. **Tabla de Masullo** (MINFLUX, $N=500$, SBR = 5, centro): **0.94 / 1.96 / 3.16 nm** para $L$ = 50/100/150 nm [MasulloSI2022, p. 9]. Eq. S31 con fwhm = 360 nm da 0.941 / 1.962 / 3.167 nm y con fwhm = 300 nm da 0.947 / 2.012 / 3.370 nm [cálculo propio]. Es decir, la tabla se reproduce con el parámetro de Eq. S17 ≈ 360 nm; ver §9.
5. **Otros valores de la tabla** (centro, $N=500$, SBR = 5): RASTMIN ($K=6\times6$) 0.74 / 1.52 / 2.35 nm; OTMIN 0.96 / 2.00 / 3.23 nm; OT/MINSTED ($L$=50/100) 1.37 / 2.73 nm; RASTMAX ($L=300$) 11.56 nm [MasulloSI2022, p. 9]. OTMIN con $L=100$: $\sigma_x=\sigma_y=2.01$ nm [MasulloSI2022, p. 7]. Para reproducirlos hace falta la geometría del texto principal, que no está disponible.
6. **Gaussianas 1D**: $\sigma(0)=\mathrm{fwhm}^2/(4\ln2\,L\sqrt N)$ [Eq. S23c, p. 16]. Con fwhm = 300, $L=300$, $N=100$: 10.82 nm [cálculo propio].
7. **Escalamiento de la SBR con $L$** (Eq. S32, fwhm = 360, SBR($L_0$=100) = 10): SBR = 2.60 para $L=50$ y 0.66 para $L=25$ [cálculo propio; p. 20].
8. **Cámara**: $\sigma_{PSF}=100$ nm y $N=400$ dan ≈5 nm [Balzarotti2017, p. 1]. Para 5 nm, la cámara ideal ($\mathrm{SBR}_c=500$) necesita ~600 fotones y MINFLUX con $L=50$ nm ~27 fotones, lo que corresponde al "22×" del resumen [p. 4, leyenda de la Fig. 3; p. 1]. Ojo: con Eq. S27 sin fondo, $L=50$ y $N=27$ dan 3.40–3.45 nm; los ~27 fotones son un valor experimental con fondo real.
9. **Convergencia del MLE**: en el origen, con SBR = 10, el MLE alcanza el CRB desde $N\approx100$; en $(50,0)$ nm con $L=75$ nm, desde $N\approx500$ [Balzarotti2017, pp. 22–23].
10. **Experimento**: $L=100$ nm, $N=500$, SBR = 13.6 dan ≈2 nm [Balzarotti2017, pp. 3–4]. Eq. S31 con fwhm = 360 da 1.78 nm [cálculo propio], del mismo orden (la PSF real era la medida). Origami: MINFLUX 2.1 nm frente a PALM ideal 5.4 nm ($N=500$, $L=70$); 1.2 nm frente a 3.8 nm ($N=1000$, $L=50$) [p. 4].
11. **Multifotón [derivación propia]**: $\tilde\sigma(\bar0)\propto1/c$; con $L=100$, $N=100$, sin fondo, cuadrática: 3.536 / 1.768 / 1.179 nm para $c$ = 1/2/3. Coincide con el ~2× de [MS2022, p. 2].
12. **Sutileza en el origen [derivación propia]**: cuadrática sin fondo, $L=100$, $N=100$: $\sigma_x(10^{-3}\,\text{nm},0)=2.739$ nm frente a 3.536 nm en $(0,0)$ exacto.

---

## 9. Discrepancias y ambigüedades de notación
1. **$L$**: en 1D es la separación entre los dos ceros [Balzarotti2017, p. 15]; en el TCP 2D es el **diámetro** del círculo de los ceros periféricos [p. 17, Eq. S24]. En RASTMIN, $L$ es el **lado** del raster cuadrado [RASTMIN2022, p. 2]. En Masullo SI, el FOV promedio es un círculo de diámetro $0.75L$ [MasulloSI2022, p. 9]. Cuando se comparen esquemas "con el mismo $L$" hay que decir qué $L$ se está usando.
2. **Factor 1D vs 2D**: $L/(4\sqrt N)$ en 1D frente a $L/(2\sqrt{2N})$ en 2D TCP. No hay que mezclarlos.
3. **"fwhm" de la dona**: en Eq. S17 es un parámetro de tamaño (anillo con diámetro pico a pico ≈ 1.2·fwhm), no el FWHM de ninguna curva de la dona [Balzarotti2017, p. 14]. RASTMIN habla de "FWHM de 300 nm" [RASTMIN2022, p. 2] sin definirlo en el texto disponible. La tabla de Masullo se reproduce con parámetro ≈ 360 nm en Eq. S17, no con 300 nm (§8.2, check 4); puede que usen otra parametrización o 360 nm, que es el valor de la fig. S3 de Balzarotti. **Esto no está resuelto.**
4. **SBR de MINFLUX vs cámara**: para MINFLUX, SBR = señal total / fondo total (Eq. S29; Masullo Eq. S3). Para la cámara de Balzarotti, $\mathrm{SBR}_c$ = señal total / fondo **por píxel** (Eq. S61), de modo que $\mathrm{SBR}_c=K\cdot\mathrm{SBR}_{total}$: $\mathrm{SBR}_c=500$ con $K=81$ equivale a SBR total ≈ 6.2 (cálculo propio). Masullo usa para la cámara la forma total/total (Eq. S7). Las SBR de cámara de ambos papers no son directamente comparables.
5. **Índices**: Balzarotti usa $i=0..K-1$ con 0 en el centro, y en 1D $p_0$ corresponde al cero en $-L/2$ ($f_0=CI(x_m+L/2)$, que crece con $x$) [p. 2]. Masullo usa $i=1..K$ [MasulloSI2022, p. 3]. Stefani 2023 define $p(x)=I_2/(I_1+I_2)$ con $I_1(x)=I(x+L/2)$, es decir, $p$ del cero en $+L/2$, que decrece con $x$ [Stefani2023, p. 1]. Los signos cambian según la convención.
6. **"Media aritmética de autovalores"** (Eq. S13): la fórmula es $\sqrt{\mathrm{tr}\Sigma/d}$, la raíz de la media de las varianzas, no la media de los σ.
7. **Eq. S43** usa $\bar p$ de Eq. S4 (sin fondo), mientras que el MLE usa Eq. S28/S30. El LMS/mLMS publicado no incluye el factor $1/s$ del fondo (§3.2, derivación propia). El sesgo resultante se absorbe en el numLMS, que sí depende de la SBR [p. 29].
8. **$p_0=0$ en el origen sin fondo**: Eq. S27 lo trata como término nulo, lo que genera una discontinuidad del CRB en $\bar r=\bar0$ (§4.3).
9. **"125 ms" vs 125 µs** en la resolución temporal del tracking (§5): probable artefacto de extracción; hay que confirmarlo en el PDF.
10. **Tabla S1 de MasulloSI2022**: la leyenda describe 4 columnas numéricas, y las columnas 3 y 5 son idénticas (p. ej., 1.29 y 1.29), lo que parece una duplicación. Además, la leyenda asigna SBR = 5 a la columna 3 y SBR = 10000 a la 4, mientras que el encabezado pone "SBR=10000" y "SBR=5" en otras posiciones. La columna del centro (col. 2) es consistente con SBR = 5 y $N=500$ según el check 4.
11. **Texto principal de Masullo et al. 2022 (Biophys. Rep.) ausente**: la "fórmula general del CRB para patrones arbitrarios" de ese paper no pudo verificarse directamente. Aquí se usa la de Balzarotti (Eq. S11), que es general. Por lo que muestra el SI (Eq. S1–S6), la construcción es la misma.
