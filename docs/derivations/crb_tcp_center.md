# CRB en el centro del TCP: valor puntual, límite $r\to0$, fondo, cero finito y multifotón

Derivación analítica autocontenida. Cada resultado en caja está implementado en
`src/donutloc/closed_forms.py`, se contrasta con una matriz de Fisher numérica propia en
`tests/test_closed_forms.py` y se reproduce con `python scripts/verify_crb_closed_forms.py`. Ese
script compara 89 casos exactos y el error relativo máximo es $2\times10^{-7}$.

Referencias (página del PDF, convención de `docs/literature/A_minflux_theory.md`):
Balzarotti et al. 2017, Material Suplementario:
- Eq. S11 y S13 (p. 13): Fisher y métrica;
- Eq. S17 (p. 14): dona;
- Eq. S21–S23 (pp. 15–16): 1D;
- Eq. S24 (p. 17): TCP;
- Eq. S27 (p. 18): centro sin fondo;
- Eq. S28–S30 (p. 19): fondo;
- Eq. S31–S32 (p. 20): centro con fondo y SBR frente a $L$.

---

## 0. Modelo, notación y supuestos

- **Patrón (Eq. S24).** Hay tres ceros periféricos en
  $\mathbf b_k=R(\cos\varphi_k,\sin\varphi_k)$ con $R=L/2$ ($L$ es el **diámetro**) y
  $\varphi_k=\pi/2+2\pi k/3$, $k=0,1,2$. El cero central es $\mathbf b_3=\mathbf 0$ y va al final,
  mientras que Balzarotti lo numera 0. La rotación global no cambia nada en el centro, porque el
  resultado será isótropo.
- **Haz.** Es radial, así que conviene escribirlo como función de $u=|\mathbf r-\mathbf b|^2$:
  $\lambda_i(\mathbf r)=h(|\mathbf r-\mathbf b_i|^2)$. Para la dona LG (Eq. S17, pico 1):
  $$h_{\rm LG}(u)=e\,a\,u\,e^{-au},\qquad a=\frac{4\ln2}{\mathrm{fwhm}^2}.$$
  El número adimensional que controla todo es
  $$x\equiv aR^2=\frac{L^2\ln2}{\mathrm{fwhm}^2},\qquad g\equiv1-x .$$
  El caso cuadrático es $x\to0$, porque la escala global se cancela en $p_i$.
- **Probabilidades (Eq. S4).** $p_i=\lambda_i/S$ con $S=\sum_{j=0}^{3}\lambda_j$.
- **Fisher (Eq. S11).** $F=N\sum_i \nabla p_i\nabla p_i^{T}/p_i$. **Métrica (Eq. S13):**
  $\sigma=\sqrt{\mathrm{tr}(F^{-1})/2}$.
- **Supuestos.** Emisor 2D. $N$ fijo (multinomial). Fondo igual en las 4 exposiciones. $x<1$, es
  decir que los ceros periféricos quedan dentro del máximo del anillo, porque de otro modo la
  pendiente cambia de signo.

**Identidades del triángulo.** Se usan todo el tiempo. Como $\sum_k e^{i\varphi_k}=0$ y
$\sum_k e^{2i\varphi_k}=0$:
$$\sum_k\mathbf b_k=\mathbf 0,\qquad \sum_k\mathbf b_k\mathbf b_k^{T}=\tfrac{3R^2}{2}\,\mathbb 1 .\tag{0.1}$$
Para la segunda: $\sum\cos^2\varphi_k=\tfrac32+\tfrac12\sum\cos2\varphi_k=\tfrac32$, y
$\sum\cos\varphi_k\sin\varphi_k=\tfrac12\sum\sin2\varphi_k=0$.

## 1. Desarrollo cerca de $\mathbf r=\mathbf 0$ y fórmula maestra

Sean $h_R\equiv h(R^2)$, $h'_R\equiv h'(R^2)$ y $h_0\equiv h(0)$. Como
$|\mathbf r-\mathbf b_k|^2=R^2-2\mathbf b_k\!\cdot\!\mathbf r+r^2$,
$$\lambda_k=h_R-2h'_R\,\mathbf b_k\!\cdot\!\mathbf r+O(r^2),\qquad
\lambda_3=h(r^2)=h_0+h'(0)\,r^2+O(r^4).\tag{1.1}$$
Por (0.1), el término lineal de $S$ se cancela: $S=S_0+O(r^2)$ con $S_0=3h_R+h_0$, y
$\nabla S(\mathbf 0)=\mathbf 0$. Entonces, en $\mathbf r=\mathbf 0$:
$$p_k=\frac{h_R}{S_0},\qquad \nabla p_k=\frac{\nabla\lambda_k}{S_0}-\frac{\lambda_k\nabla S}{S_0^2}=-\frac{2h'_R}{S_0}\mathbf b_k,
\qquad \nabla p_3(\mathbf 0)=\mathbf 0\ \text{(exacto, }\lambda_3\text{ es par)}.\tag{1.2}$$
Los tres términos periféricos de Eq. S11, usando (0.1), dan
$$F_{\rm per}=N\sum_k\frac{S_0}{h_R}\frac{4h_R'^2}{S_0^2}\mathbf b_k\mathbf b_k^T
=\boxed{\;F=\frac{6NR^2\,h_R'^2}{h_R\,(3h_R+h_0)}\;\mathbb 1\;}\tag{M}$$
Esta es la **fórmula maestra**. Vale para cualquier haz radial, siempre que el término central
sea nulo en el origen: $\nabla p_3=0$ con $p_3>0$, o bien $p_3=0$ excluido por la convención del
valor puntual. $F$ es isótropa, así que $\sigma_x=\sigma_y$.

**Fondo como constante aditiva.** Sumar un fondo $b$ por exposición equivale a
$h_R\to h_R+b$ y $h_0\to h_0+b$, con $h'$ sin cambios. En (M) basta reemplazar
$h_R\to h_R+b$ y $3h_R+h_0\to S_0+4b$.

## (a) Valor puntual en $r=0$: Eq. S27

Para LG: $h_R=e\,x\,e^{-x}$, $R^2h'_R=e\,x(1-x)e^{-x}$ y $h_0=0$. En (M):
$$F=\frac{6N\,e^2x^2g^2e^{-2x}}{R^2\,e x e^{-x}\,3exe^{-x}}=\frac{2Ng^2}{R^2}\mathbb 1=\frac{8Ng^2}{L^2}\mathbb 1
\;\Rightarrow\;
\boxed{\sigma_{\rm pt}=\frac{L}{2\sqrt{2N}}\Big(1-\frac{L^2\ln2}{\mathrm{fwhm}^2}\Big)^{-1}}\quad\text{(Eq. S27, p. 18)}.$$
El término central se excluyó porque en el origen es $0/0$ ($p_3=0$, $\nabla p_3=0$). Con fwhm
finito aparece el factor $g^{-1}$ porque la pendiente de la dona a distancia $R$ es
$\propto(1-x)$.

## (b) Límite $r\to0$ sin fondo

Para $\mathbf r\neq\mathbf 0$ se tiene $p_3>0$, y el término central no es cero. Con $h_0=0$ y
$h'(0)=ea$, (1.1) da
$$p_3=\frac{ea\,r^2}{S_0}\big(1+O(r^2)\big),\qquad
\nabla p_3=\frac{2ea\,\mathbf r}{S_0}\big(1+O(r^2)\big)$$
(la corrección $-\lambda_3\nabla S/S^2$ es $O(r^3)$). Entonces
$$T(\mathbf r)\equiv\frac{\nabla p_3\nabla p_3^T}{p_3}=\frac{4ea}{S_0}\,\hat{\mathbf r}\hat{\mathbf r}^T+O(r^2)
=\frac{16\,e^{x}}{3L^2}\,\hat{\mathbf r}\hat{\mathbf r}^T+O(r^2),$$
porque $4ea/(3exe^{-x})=4e^x/(3R^2)$. En el caso cuadrático esto da $16/(3L^2)$, lo que coincide
con la nota A §4.3. Los términos periféricos son continuos y solo cambian en $O(r)$; esos
cambios son impares en $\hat{\mathbf r}$ y se cancelan al promediar direcciones. Así:
$$F_{\lim}(\hat{\mathbf r})=\alpha\,\mathbb 1+\beta\,\hat{\mathbf r}\hat{\mathbf r}^T,\qquad
\alpha=\frac{8Ng^2}{L^2},\quad \beta=\frac{16N e^{x}}{3L^2},\quad \frac\beta\alpha=\frac{2e^x}{3g^2}.$$
Los autovalores son $\alpha+\beta$ (a lo largo de $\hat{\mathbf r}$) y $\alpha$ (en la dirección
perpendicular). Entonces $\mathrm{tr}F^{-1}=1/\alpha+1/(\alpha+\beta)$ **no depende de la dirección
de aproximación**, y el promedio sobre 12 direcciones del test de aceptación es redundante: la
dispersión numérica es $<10^{-4}$. Por lo tanto
$$\boxed{\sigma_{\lim}^2=\frac{L^2}{8Ng^2}\,\frac{3g^2+e^{x}}{3g^2+2e^{x}},\qquad
\rho\equiv\frac{\sigma_{\lim}}{\sigma_{\rm pt}}=\sqrt{\frac{3g^2+e^{x}}{3g^2+2e^{x}}}}$$

- **Hipótesis 2/√5: CONFIRMADA en el límite cuadrático.** Con $x=0$ queda $\rho^2=4/5$, es decir
  $\rho=2/\sqrt5=0.894427$, y $\sigma^2_{\lim}=L^2/(10N)$.
- **Con fwhm finito el cociente NO se conserva.** Para $x\ll1$,
  $\rho\simeq\frac{2}{\sqrt5}\big(1-\frac{9x}{40}\big)$. La razón es que el término periférico
  escala como $g^2=(1-x)^2$, por la pendiente reducida (el origen del factor S27), mientras que el
  término central escala como $e^{x}$: la curvatura del cero central es fija ($ea$) y el
  denominador $S_0\propto xe^{-x}$ cae. Los valores con fwhm = 300 nm son $\rho$ = 0.89439
  ($L$=5), 0.89050 ($L$=50), 0.87805 ($L$=100) y 0.85527 ($L$=150). Coinciden con los
  observados, 0.8944 / 0.8905 / 0.8780.
- **Ejes.** $\sigma_\parallel=1/\sqrt{\alpha+\beta}$ y $\sigma_\perp=1/\sqrt\alpha=\sigma_{\rm pt}$.
  La isotropía es $\sqrt{3g^2/(3g^2+2e^x)}$, que vale $\sqrt{3/5}$ en el caso cuadrático. Con
  $L=100$, $N=100$ (cuadrático) queda $\sigma_\parallel=2.7386$ nm, igual a la nota A (2.739),
  y $\sigma_\perp=3.5355$ nm.
- **Número de aceptación.** $L=50$, $N=100$, fwhm = 300 da $\sigma_{\lim}=1.605096$ nm, frente al
  valor puntual S27 de 1.802472 nm.

### Por qué hay una discontinuidad en $r=0$
1. $T(\mathbf r)$ tiene norma fija ($\beta/N$) pero su dirección es $\hat{\mathbf r}$: como matriz
   **no tiene límite** en $\mathbf 0$. Su traza sí lo tiene. Asignarle el valor 0 en el origen
   (convención S27) no coincide con ningún límite direccional.
2. Una forma equivalente es $F=4N\sum_i\nabla\sqrt{p_i}\,\nabla\sqrt{p_i}^T$, con
   $\sqrt{p_3}\simeq\sqrt{ea/S_0}\,|\mathbf r|$. Eso es un **cono**, no diferenciable en el
   origen, de modo que la hipótesis de regularidad de la cota de Cramér-Rao (diferenciabilidad de
   $\sqrt p$) falla exactamente en $r=0$. Allí ningún valor de $F$ es "el correcto". S27 es una
   convención, y el límite describe a un emisor a cualquier distancia pequeña pero no nula.
3. Físicamente, $p_3\propto r^2$ hace que la variación relativa $dp_3/p_3=2\,dr/r$ diverja justo
   cuando $p_3\to0$, y el producto queda finito. La exposición central aporta información solo
   **radial** ($\hat{\mathbf r}\hat{\mathbf r}^T$), la misma a cualquier distancia (invariancia de
   escala de un cero cuadrático), y un conteo nulo es informativo.

## (c) Fondo con SBR fija (Eq. S30) → Eq. S31, y límites que no conmutan

Eq. S30, $p_i=s\,p_i^{(0)}+(1-s)/4$ con $s=\mathrm{SBR}/(\mathrm{SBR}+1)$, es idéntica a sumar
$b=S/(4\,\mathrm{SBR})$ a cada $\lambda_i$, porque $(\lambda_i+b)/(S+4b)=s\lambda_i/S+(1-s)/4$.
Con $h_0=0$ y $S_0=3h_R$: $h_R+b=h_R\big(1+\tfrac{3}{4\mathrm{SBR}}\big)$ y
$S_0+4b=S_0\big(1+\tfrac1{\mathrm{SBR}}\big)$. En (M):
$$\boxed{\sigma(\mathbf 0)=\frac{L}{2\sqrt{2N}}\,g^{-1}\sqrt{\Big(1+\frac1{\rm SBR}\Big)\Big(1+\frac{3}{4\,\rm SBR}\Big)}}\quad\text{(Eq. S31, p. 20)}.$$
Ahora $p_3\ge(1-s)/4>0$ y $\nabla p_3=s\nabla p_3^{(0)}=O(r)$, así que $T=O(r^2)\to0$: el CRB es
**continuo** y el valor puntual coincide con el límite (verificado a $10^{-9}$ con SBR = 5 y 10).

**Transición.** Cerca del origen, $\lambda_3+b\simeq ea\,r^2+b$, y entonces
$$T(\mathbf r)\simeq\frac{4ea}{S_0+4b}\,\frac{r^2}{r^2+r_c^2}\,\hat{\mathbf r}\hat{\mathbf r}^T,\qquad
r_c^2=\frac{b+\epsilon}{ea}\ \ \Rightarrow\ \ r_c=\frac{\sqrt3\,L}{4}\,\frac{e^{-x/2}}{\sqrt{\rm SBR}}\ \ (\epsilon=0)$$
($\epsilon$ es el pedestal constante de §d).

**Los límites no conmutan:**
$$\lim_{\rm SBR\to\infty}\lim_{r\to0}\sigma=\sigma_{\rm pt}\ (\text{S27}),\qquad
\lim_{r\to0}\lim_{\rm SBR\to\infty}\sigma=\sigma_{\lim}=\rho\,\sigma_{\rm pt}.$$
El término central se enciende en una escala $r_c\propto\mathrm{SBR}^{-1/2}$ que se reduce a un
punto. Con SBR fija, en $r=0$ vale exactamente 0 para toda SBR, y por eso S31 → S27. Con $r$ fijo,
SBR → ∞ lo lleva a su valor completo $\beta$.

Chequeo numérico ($L$=50, fwhm=300):
- SBR = $10^4$ y $r_0=10^{-3}\ll r_c$ da el valor puntual a $10^{-4}$;
- SBR = $10^{14}$ y $r_0\gg r_c$ da $\sigma_{\lim}$ a $10^{-3}$.

La fórmula de transición es solo de **orden dominante** ($r_c\ll L$):
- error $10^{-4}$ con $r_c=0.42$ nm;
- 1.6 % con $r_c=4.9$ nm ($L$=100);
- falla con $r_c=13$ nm (SBR = 10), donde dominan los términos $O(r/L)$.

Con SBR realistas (5–20), $r_c\approx9$–$19$ nm para $L$=100: la discontinuidad desaparece en la
práctica y S31 es el valor correcto en el centro.

**Perfil radial numérico** ($L$=100, $N$=100, fwhm=300; $\sigma$ en nm, promedio de 12
direcciones). Con un cero casi perfecto el **mínimo del CRB no está en el centro**, sino en un
anillo de radio $\sim r_c$:

| $r$ (nm) | 0 | 1 | 2 | 5 | 10 | 15 | 20 | 30 |
|---|---|---|---|---|---|---|---|---|
| $\epsilon=0$, SBR=∞ | 3.831 (pt) | 3.365 | 3.371 | 3.410 | 3.553 | 3.800 | 4.165 | 5.341 |
| $\epsilon=0.002$ (const.) | 3.877 | 3.848 | 3.781 | 3.628 | 3.665 | 3.885 | 4.243 | 5.425 |
| $\epsilon=0.01$ (const.) | 4.061 | 4.057 | 4.046 | 3.998 | 4.008 | 4.193 | 4.539 | 5.752 |
| $\epsilon=0$, SBR=10 | 4.165 | 4.163 | 4.156 | 4.129 | 4.162 | 4.367 | 4.749 | 6.105 |

## (d) Profundidad finita del cero $\epsilon$

$\epsilon$ es la intensidad en $r=0$ relativa al pico del anillo.

### Modelo "constant": $h=h_{\rm LG}+\epsilon$ (exacto)
Sumar $\epsilon$ a **todas** las exposiciones es exactamente un fondo $b=\epsilon$ (Eq. S28 con
$\lambda_b=\epsilon$), válido para todo $\mathbf r$. En el origen $\nabla\sum I_{\rm LG}=0$, así que
fondo fijo y SBR fija dan la misma $F$. Con (M) y el reemplazo $b=\epsilon$:
$$\boxed{\sigma_\epsilon(\mathbf 0)=\sigma_{\rm pt}\sqrt{\Big(1+\frac1{\mathrm{SBR}_\epsilon}\Big)\Big(1+\frac{3}{4\,\mathrm{SBR}_\epsilon}\Big)},\qquad
\mathrm{SBR}_\epsilon=\frac{\sum_jI_{{\rm LG},j}(\mathbf 0)}{4\epsilon}=\frac{3e\,x\,e^{-x}}{4\epsilon}}$$
Ejemplo: $L$=100, fwhm=300 y $\epsilon$=0.01 dan $\mathrm{SBR}_\epsilon=14.54$. $p_3(\mathbf 0)>0$,
así que no hay discontinuidad: el valor en el origen es continuo y, cuando $\epsilon\to0$, tiende a
**S27**, no al límite (b). Es la misma no conmutatividad de (c).

**Fondo adicional: dos convenciones.** Las dos son exactas y están en el código como `sbr_ref`.
- `"lg"`: la SBR se mide contra la señal LG pura, con $\lambda_b=\sum I_{\rm LG}(\mathbf 0)/(4\,{\rm SBR})$.
  Queda $b=\epsilon+\lambda_b$, es decir
  $1/\mathrm{SBR}_{\rm eff}=1/\mathrm{SBR}+1/\mathrm{SBR}_\epsilon$ (la fórmula del plan).
- `"beam"`/`"total"` (default): Eq. S30 aplicada sobre el haz **con** pedestal. Es la convención de
  `donutloc.photons.probabilities(sbr=...)`, con
  $\lambda_b=(\sum I_{\rm LG}+4\epsilon)/(4\,\rm SBR)$. Queda
  $1+1/\mathrm{SBR}_{\rm eff}=(1+1/\mathrm{SBR})(1+1/\mathrm{SBR}_\epsilon)$, es decir
  $1/\mathrm{SBR}_{\rm eff}=1/\mathrm{SBR}+\frac{4\epsilon}{\sum I_{\rm LG}}(1+1/\mathrm{SBR})$.

En los dos casos, $\sigma=$ Eq. S31 evaluada en $\mathrm{SBR}_{\rm eff}$.

### Modelo "gaussian": $h=(ea\,u+\epsilon)\,e^{-au}$ (exacto vía (M); NO es un fondo)
Se tiene $h_R=(ex+\epsilon)e^{-x}$, $R^2h'_R=(ex\,g-\epsilon x)e^{-x}$ y $h_0=\epsilon$. En (M),
dividiendo por $F_{\rm pt}=2Ng^2/R^2$ y multiplicando numerador y denominador por $e^{2x}$:
$$\boxed{\frac{F_\epsilon}{F_{\rm pt}}=\frac{3x^2\,\big(e\,g-\epsilon\big)^2}{g^2\,(ex+\epsilon)\,\big(3(ex+\epsilon)+\epsilon e^{x}\big)},\qquad
\sigma_\epsilon=\sigma_{\rm pt}\sqrt{F_{\rm pt}/F_\epsilon}}$$
(con fondo: $h_R\to h_R+b$ y $S_0\to S_0+4b$). **No equivale exactamente a un fondo efectivo**,
por dos razones:
1. el pedestal depende de la exposición: vale $\epsilon$ en la central y $\epsilon e^{-x}$ en las
   periféricas;
2. reduce la pendiente en el factor $(1-\epsilon/(eg))$.

A primer orden en $\epsilon$:
$$\ln\frac{F^{\rm const}_\epsilon}{F_{\rm pt}}=-\frac{7e^{x}}{3ex}\epsilon,\qquad
\ln\frac{F^{\rm gauss}_\epsilon}{F_{\rm pt}}=-\Big[\frac{2}{eg}+\frac{6+e^x}{3ex}\Big]\epsilon .$$
El cociente de los coeficientes es $1+\tfrac{3}{7}x^2+O(x^3)$ (numérico: 1.00016, 1.0027 y
1.0147 para $L$ = 50/100/150 con fwhm = 300). El fondo efectivo del modelo constante aproxima
entonces al gaussiano con error relativo $O(x^2)$ en el término lineal. Las diferencias que
quedan son $O(\epsilon^2)$: 1.6 % en $\sigma$ con $\epsilon=0.1$ y $L=100$.

### Tabla genérica de degradación (sin fondo, $r=0$)
Cada celda es CRB($\epsilon$)/CRB$_{\rm S27}$ en la forma **constante | gaussiano**. Entre
corchetes, el cociente contra el límite $r\to0$ con $\epsilon=0$ (modelo constante), que es la
referencia del test de aceptación. $\epsilon=0.002$ es la profundidad reportada por Balzarotti
(p. 6). Los demás valores son genéricos.

| fwhm | $L$ | 0.002 | 0.01 | 0.03 | 0.05 | 0.10 | 0.15 |
|---|---|---|---|---|---|---|---|
| 300 | 50 | 1.045\|1.045 [1.174] | 1.227\|1.228 [1.378] | 1.679\|1.687 [1.885] | 2.130\|2.152 [2.392] | 3.256\|3.344 [3.657] | 4.382\|4.584 [4.920] |
| 300 | 100 | 1.012\|1.012 [1.153] | 1.060\|1.061 [1.207] | 1.180\|1.183 [1.344] | 1.300\|1.307 [1.481] | 1.600\|1.626 [1.822] | 1.898\|1.958 [2.162] |
| 300 | 150 | 1.006\|1.006 [1.176] | 1.029\|1.030 [1.204] | 1.088\|1.091 [1.272] | 1.147\|1.153 [1.341] | 1.294\|1.312 [1.513] | 1.440\|1.479 [1.684] |
| 360 | 50 | 1.065\|1.065 [1.194] | 1.324\|1.326 [1.485] | 1.971\|1.982 [2.210] | 2.616\|2.647 [2.934] | 4.228\|4.353 [4.741] | 5.838\|6.125 [6.547] |
| 360 | 100 | 1.017\|1.017 [1.151] | 1.085\|1.085 [1.228] | 1.253\|1.257 [1.419] | 1.422\|1.431 [1.610] | 1.843\|1.878 [2.086] | 2.263\|2.342 [2.562] |
| 360 | 150 | 1.008\|1.008 [1.161] | 1.040\|1.041 [1.198] | 1.121\|1.123 [1.291] | 1.201\|1.206 [1.383] | 1.401\|1.421 [1.614] | 1.601\|1.646 [1.844] |

La tabla cuantifica la advertencia de la persona: con $L$ chico, la pérdida por $\epsilon$ es
grande, porque $\mathrm{SBR}_\epsilon\propto x=L^2\ln2/\mathrm{fwhm}^2$. Con $L=50$, fwhm = 300 y
$\epsilon=0.05$, el CRB se duplica. Con $\epsilon\approx$ 0.2 % la pérdida es de 0.6–6.5 %, pero el
radio de transición $r_c=\sqrt{\epsilon/(ea)}$ vale 4.9 nm (fwhm=300) y no depende de $L$.

## (e) Multifotón con exponente $c$

Con $\lambda\propto I^c$ se tiene $h(u)=(ue^{-au})^c$ (la constante se cancela),
$h'_R/h_R=c\,g/R^2$ y $h_0=0$. En (M), $F=2NR^2(h'_R/h_R)^2=2Nc^2g^2/R^2$:
$$\boxed{\sigma^{(c)}_{\rm pt}=\frac{L}{2c\sqrt{2N}}\,g^{-1}}$$
Con $g=1$ se recupera la nota A §6. En el límite, $p_3\propto r^{2c}$ y
$T\propto r^{2c-2}\to0$ para $c\ge2$: **no hay discontinuidad** y el límite coincide con el valor
puntual. Con fondo de SBR fija aparece el mismo factor de S31, con la misma álgebra (el factor de
SBR es independiente de $c$, pero la SBR física en el centro sí depende de $c$). Valores
numéricos con $L$ = 100 y $N$ = 100:
- cuadrático: 3.536 / 1.768 / 1.179 nm para $c$ = 1 / 2 / 3;
- fwhm = 300: 3.831 / 1.915 / 1.277 nm.

## Complementos: 1D y Eq. S32

- **1D** (Eq. S21: $\sigma=\sqrt{p_0(1-p_0)}/(|p_0'|\sqrt N)$, con $p_0=1/2$ en el origen). Los
  casos son:
  - cuadrático $L/(4\sqrt N)$ (S22c);
  - dona $L/(4\sqrt N)\,g^{-1}$ (S22f; el mismo $x$, porque los ceros están a $L/2$);
  - gaussiana $\mathrm{fwhm}^2/(4\ln2\,L\sqrt N)$ (S23c).

  Los tres están verificados contra S21 numérico a $3\times10^{-10}$.
- **Eq. S32**: con $\lambda_b$ fijo, $\mathrm{SBR}(\mathbf 0)\propto\sum I(\mathbf 0)=3exe^{-x}\propto L^2e^{-L^2\ln2/\mathrm{fwhm}^2}$.
  Eso es S32. Con fwhm = 360 y SBR(100) = 10 se obtiene SBR = 2.602 ($L$=50) y 0.657 ($L$=25).

## Tabla numérica de chequeo

Salida de `scripts/verify_crb_closed_forms.py`, $N=100$ salvo donde se indica.

| Caso | Forma cerrada (nm) | Numérico propio | err. rel. |
|---|---|---|---|
| S27, $L$=50, fwhm=300 | 1.802472 | 1.802472 | 2.1e-9 |
| límite, $L$=50, fwhm=300 (**aceptación**) | 1.605096 | 1.605096 | 2.3e-9 |
| S27, $L$=100, fwhm=∞ | 3.535534 | 3.535534 | 5.4e-10 |
| límite, $L$=100, fwhm=∞ ($=L/\sqrt{10N}$) | 3.162278 | 3.162278 | 5.3e-10 |
| S27, $L$=100, fwhm=360 | 3.735312 | 3.735312 | 5.2e-10 |
| S31, $L$=100, fwhm=300, SBR=10 (punto y límite) | 4.165447 | 4.165447 | ≤5.8e-10 |
| $\epsilon$=0.01 const., $L$=100 | 4.060972 | 4.060972 | 4.9e-10 |
| $\epsilon$=0.1 gauss., $L$=100, SBR=5 (beam) | 7.372553 | 7.372553 | 3.3e-10 |
| $c$=2, $L$=100, fwhm=300 | 1.915274 | 1.915274 | 1.1e-9 |

Checks de la nota A §8.2 (tests):
- 3.536;
- 3.735 y 1.792 (fwhm=360);
- 1.817 ($L$=100, $N$=500, SBR=10, fwhm=360);
- tabla de Masullo ($N$=500, SBR=5): 0.941 / 1.962 / 3.167 con fwhm=360, y 0.947 / 2.012 / 3.370
  con fwhm=300;
- 10.82 (1D gaussiana);
- 1.25 (1D cuadrática);
- 2.60 / 0.66 (S32);
- 2.739 ($\sigma_\parallel$).
