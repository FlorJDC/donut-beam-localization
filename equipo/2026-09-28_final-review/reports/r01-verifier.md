# r01 — verifier (final-review, explorador)

Scripts propios: `equipo/2026-09-28_final-review/work/verify/indep.py` (CRB/L_opt/M/p_i desde cero, Fisher multinomial por diferencias finitas) y `work/verify/bg.py` (sesgo MLE sin fondo por minimización de KL sobre cuentas esperadas + CRB con b conocido/libre, 3 parámetros).

## 1. Guía "Qué esperar" frente a la salida de `--all --yes --no-show --save /tmp/expl_verify` (3.7 s en total, 8 PNG)
Todos los números de la guía coinciden con la salida (con el redondeo indicado): modo 1 (180.2 nm, 0.6006, p centro/L/4), modo 2 (1.605/1.802, pendiente 1.01 ≈ 1.0), modo 3 (sesgo (-0.06,-0.01), 3.33/3.34, 1.00), modo 4 (10.5/23.3/50.3; regla 10.5/23.4/52.3), modo 5 (4.17 ≈ 4.2, 0.83·δ), modo 6 (6.84/4.92/2.53/1.57/2.69/5.25; CRB 2.02→5.16 y 2.03→6.17), modo 7 (4.63/30.02/4.55; σfl 29.67/0.00), modo 8 (M diag 0.9561, siguiente 0.042; 3.84 nm; 1.074). **Sin discrepancias.**

## 2. Contraste con valores del proyecto (camino propio)
- CRB centro L=50, N=100, fwhm=300, sin fondo: límite r→0 = 1.6051, puntual S27 (excluyendo la exposición con p=0) = 1.8025. Coincide.
- Guía modo 2 "SBR=10 → 1.96 nm, sin discontinuidad": propio 1.9601 (r=0) y 1.9601 (r→0). Coincide.
- L_opt con pedestal gaussiano I_LG + eps·exp(-4ln2 r²/fwhm²): 10.484/23.280/50.337 nm (paper_numbers.json: 10.484/23.280/50.337). Coincide. (Con pedestal constante sale 10.49/23.35/51.1: el valor de eps=0.05 depende del modelo del cero.)
- Fondo, SBR0=5, N=500, L=100: sesgo x sin fondo 4.919 (x=10), 5.247 (x=50); CRB b conocido/libre en x=50: 5.161/6.173. Coincide. Observación: el sesgo tiene también componente y (≈2.8 nm en x=10, ≈3.6 en x=50) que el explorador no muestra; la tabla dice "sesgo x", así que no es falso.
- M con τ=4 por forma cerrada (1−q)q^m/(1−q⁴), q=e^{−12.5/τ}: 0.95607, 0.04201, 0.00185, 0.00008; orientación columna=pulso, fila=ventana correcta (el fotón retrasado cae en la ventana siguiente). Coincide.
- Modo 1 p_i en x=L/4 con TCP rotado π/2: [0.3197, 0.5063, 0.1050, 0.0691]. Coincide.
No reproduje de forma independiente los Monte Carlo de los modos 3, 5 y 7 (solo confirmé que la salida coincide con la guía).

## 3. Lecciones
- **Modo 2 (refuted, menor):** "S27 mayor que el límite (razón ~2/sqrt(5))". 2/√5 = 0.894 < 1: contradice "mayor". La razón real es 1.8025/1.6051 = 1.123 ≈ √5/2 (1.118). Hay que invertir la fracción (y "~" cubre la diferencia de 0.4 %).
- **Modo 8 (refuted como está escrito):** "CRB +4–19 % para τ = 3–5 ns". Ese rango es la media sobre el disco de `informe/data/pminflux_timing.json` (1.039 / 1.189). La métrica que imprime el explorador (media sobre x=0..50, SBR=20) da 1.026 (τ=3), 1.074 (τ=4), 1.144 (τ=5), así que con τ=5 el usuario ve +14 %, no +19 %, y con τ=4 ve 7.4 %, cuando el rango haría esperar ~10 %. Arreglo: decir "+3–14 % en esta métrica (+4–19 % promediado sobre el disco, informe)", o calcular sobre el disco.
- Modo 6: "≈ +5 nm con SBR = 5" es flojo: la tabla va de 1.6 a 6.8 nm según x. No es falso en x=10 y x=50; sugiero "de 1.6 a 6.8 nm según x".
- Modo 4 (0.78 solo para eps ≲ 0.01; por debajo para eps=0.05: 50.3 < 52.3), modo 5 (0.83 vs 0.78), modo 7 (período de 50 ns = 4×12.5), modo 1 y modo 3: consistentes con la salida.

## 4. Modo interactivo
- `printf '3\n\n\n\n\n100\n0\n' | python3 explore/donut_explorer.py --no-show`: imprime el menú y "elegí un modo [0]", toma el valor por defecto 0 sin leer stdin y **sale con código 0 sin correr nada**. Tal como lo pedía la tarea, usa los valores por defecto sin tty, pero el valor por defecto del menú es "salir", así que la entrada enviada por pipe se ignora en silencio. Sugiero avisar ("stdin no es interactivo: usá --mode/--all --yes").
- `--mode 9`: error de argparse "--mode debe estar entre 1 y 8", código 2. Correcto.
