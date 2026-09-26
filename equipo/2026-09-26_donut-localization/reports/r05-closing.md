# Cierre del trabajo — 2026-09-26_donut-localization

Estado: **done**. Presupuesto: 5 rondas (usadas 5), con una pasada de correcciones y un chequeo
final dentro de la ronda 5.

[[DONE]] — condición cumplida: el entregable está completo, `tests/test_acceptance.py` pasa (8/8,
con el manuscrito completo: 179 etiquetas de procedencia, 39 afirmaciones, 217 números), la suite
unitaria pasa (224 tests) y toda afirmación que entra al paper fue reproducida por un verificador
independiente. Las afirmaciones refutadas en la última verificación se corrigieron y se volvieron a
comprobar; los tres retoques residuales del chequeo final (redacción, sin cambios de números) los
aplicó el orquestador.

## Qué se entregó
- `src/donutloc/`: haces (LG, vectorial Richards-Wolf), patrones, modelo de fotones, Fisher/CRB,
  formas cerradas, estimadores, Monte Carlo, cámara, experimentos.
- `scripts/`: 8 figuras, `compute_paper_numbers.py`, `check_provenance.py`, `reproduce.py/.sh`.
- `paper/main.pdf` (13 pp., revtex4-2) con procedencia de cada número.
- `data/paper_numbers.json`, `structure/{claims,figures}.json`, README, CITATION.

## Resultados principales (verificados)
1. El CRB en el centro del TCP es discontinuo sin fondo: el límite r→0 es menor que el valor puntual
   de Balzarotti Eq. S27 (razón 2/√5 para cero cuadrático; forma cerrada exacta para fwhm finito).
2. La dona vectorial con la handedness correcta equivale (CRB <0.4 %) a una LG igualada por
   curvatura; la handedness opuesta o la polarización lineal llenan el cero y degradan el CRB 7–90×.
3. El MLE es eficiente con fondo; sin fondo es superficiente y sesgado hacia el centro.
4. MINFLUX iterativo con 1000 fotones: ~0.51 nm frente a 3.16 nm de una cámara ideal; la precisión
   mejora algo más rápido que 1/√N.
5. Con cero imperfecto hay un L óptimo (≈0.78·fwhm·√ε para ε ≲ 0.01); la desalineación con estimador
   ingenuo produce un sesgo de ~0.75–0.78·δ en el centro, que desaparece con el patrón verdadero.

## Pendiente para la autora
- Confirmar SBR = 10 como condición de `mle_efficiency_center` (declarado en el paper).
- Ver `backlog` en `state.json` y la sección "Open points" del paper.
