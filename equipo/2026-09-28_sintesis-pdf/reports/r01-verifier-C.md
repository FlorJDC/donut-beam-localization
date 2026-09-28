# r01 — verificador (worker C: sesgo por fondo en TCP fijo)

No leí `r01-worker-C.md`. Verifiqué contra `informe/data/background_bias.json` y los artefactos.

## Método propio
`work/verify/verify_background.py` (salida en `work/verify/verify_background.out`, unos 4 min). No importa
`donutloc`. Reimplementé la dona LG (`e·a·r²·exp(-a r²)`, con `a=4ln2/fwhm²`) y el TCP (radio L/2, ángulos
π/2+2πk/3, más el centro). Para la multinomial uso p_i=(λ_i+b)/Σ(λ_j+b). El Fisher es N·Σ ∂p∂pᵀ/p con
**derivadas numéricas centradas**, y el worker usa gradiente analítico. El MLE combina una grilla vectorizada
(paso de 2 nm, radio de 150 nm) con un refinamiento Nelder-Mead. Con b libre agrego una grilla en b
(0 y 25 valores log hasta 0.3) y comparo contra el ajuste con b=0 exacto en la frontera. El sesgo
poblacional sale de un MLE sobre cuentas esperadas con 10 arranques distintos. Hice 1500 repeticiones
por punto (600 en el caso de b libre) con semilla `[42,7,x]`, que no es la del worker.

## Resultados
1. b = 0.0290751 / 0.00726877. SBR(50) = 10.2652 / 41.0609. Lo verifiqué además a mano: cada dona del
   anillo aporta λ(0)=0.193834 y Σ=0.581503.
2. El CRB con b conocido da x=0: 2.0124 / 1.7880; x=50: 5.1605 / 4.5135. Coincide con el worker en 4 cifras.
3. El CRB con b libre coincide con el de b conocido en x=0 (F_xb, F_yb ~ 1e-9, o sea cero numérico). En x=50
   da 6.1734 / 5.6905. Los valores intermedios (5, 10, 20 nm) también coinciden.
4. Sesgo poblacional del MLE sin fondo: +4.919 nm (SBR₀=5, x=10), +5.247 (x=50) y +2.688 (SBR₀=20, x=5).
   Coincide con el JSON. En el MC propio obtuve +4.81±0.05 (x=10), +5.43±0.20 (x=50) y +2.22±0.04
   (SBR₀=20, x=5). El RMSE por eje en x=0 es 7.92 nm frente a 7.93 del worker. El signo es hacia afuera (+x).
   - **Matiz 1.** El sesgo **no es solo en x**: en y es de +2.5 a +3.6 nm con SBR₀=5 (MC: 2.60 en x=10,
     2.86 en x=50). La afirmación menciona solo el eje x.
   - **Matiz 2.** En x=0 el ajuste poblacional sin fondo es **degenerado**: hay 3 mínimos equivalentes a
     10.45 nm del centro (SBR₀=5) o a 5.22 nm (SBR₀=20), uno en cada dirección de la simetría C3. El valor
     `bias_pop_y[0] = -10.45` del JSON es uno de ellos, elegido arbitrariamente. No es un sesgo en y. Ese
     número no debe citarse como sesgo; el RMSE de 7.9 nm en x=0 refleja justamente ese salto entre los 3 mínimos.
5. MLE con fondo exacto: el sesgo poblacional da 0 (el JSON trae ~3e-5, que es la tolerancia).
   Mi σ/CRB da 0.970 (x=0), 1.011, 1.021 y 1.012 con SBR₀=5, y 0.963 y 0.993 con SBR₀=20.
   El rango "1.00–1.04" **no es exacto** ni siquiera en el propio JSON: da 0.975 (SBR₀=20, x=0) y 1.057
   (SBR₀=5, x=40). El rango correcto es ≈0.97–1.06, compatible con ruido MC más efectos de N finito.
   Con N finito, σ<CRB es posible porque el MLE queda algo sesgado: por ejemplo, en x=50 el sesgo en y es
   −0.56±0.08 nm.
6. MLE con b libre (600 repeticiones), SBR₀=5: 1.963 vs 2.012, 2.114 vs 2.084 y 2.402 vs 2.364. Cada
   diferencia es ≤ 1.3 veces el error estándar de σ (~2.9 %). La fracción de b̂=0 queda por debajo de 0.3 %
   hasta x=20, así que el CRB libre (interior) es aplicable ahí. **Fuera de ese rango no lo es**: con SBR₀=20
   el JSON reporta b̂=0 en el 16–30 % de los casos desde x≥20, y ahí el CRB libre no es una cota válida.
   La afirmación se limita a x≤20 con SBR₀=5, así que sigue en pie.
7. Con ±30 % de error en el fondo, mi sesgo poblacional máximo por eje es 1.19 nm (SBR₀=5) y 0.47 nm
   (SBR₀=20). En el MC del worker el máximo llega a 1.63 nm (y, +30 %, x=50), pero incluye el sesgo de N
   finito, que ya existe con fondo exacto (−0.63). "≤ ~1.6 nm" se sostiene.

## Ejecución
- `python3 -m unittest tests.test_background`: 15 tests OK (1.6 s).
- `python3 scripts/fig_10_background_bias.py --quick`: corre en 11 s. Escribe en `informe/data/quick/` e
  `informe/figures/quick/` y **no** pisa el JSON final (lo comprobé con `cmp`). Esos subdirectorios
  `quick/` los creó mi corrida (u otra anterior), así que conviene no versionarlos.
- La figura muestra coherentemente el sesgo en x y el RMSE por eje, con σ_CRB conocido y libre.

## No resuelto
- No reproduje por MC propio el caso de b libre en x>20 ni el de ±30 % en MC; en ambos solo verifiqué el
  sesgo poblacional.
