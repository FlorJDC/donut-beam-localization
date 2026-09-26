# Intent — donut-beam-localization

El contrato completo es `OBJECTIVE.md` en la raíz del proyecto (leerlo entero).

Tipo: híbrido `feature` + `derive`: el trabajo edita el PROYECTO (workdir = raíz del repo); el
entregable es el repositorio completo tipo "paper companion", con `paper/main.tex` como documento
legible y `paper/provenance.json` como registro de procedencia.

Definición de terminado: `python -m unittest tests.test_acceptance` pasa (hoy: SKIP/rojo) y
`python -m unittest discover -s tests` pasa completo.

Equipo: pi · workers (1-3 según la ronda) · verifier (física: re-deriva y recalcula de forma
independiente) · code-reviewer (código: diff real, tests) · writer (al final: paper + procedencia).

Rondas previstas (presupuesto 5): 1) teoría + núcleo del paquete; 2) CRB, estimadores y Monte
Carlo; 3) dona vectorial, iterativo y no idealidades; 4) figuras + números; 5) manuscrito,
procedencia, README.

Anti-objetivos: no copiar código de terceros; no subir PDFs; no agregar dependencias
obligatorias; no editar `tests/test_acceptance.py`.
