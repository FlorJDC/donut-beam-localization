# Intent — síntesis del proyecto en un PDF breve (≤ 5 páginas)

Pedido de la autora (2026-09-28): el repositorio debe quedar completo hoy con todo el material
producido y **un PDF breve (5 páginas máx.) con texto y figuras** que explique el trabajo realizado.

Tipo: `draft` + un `feature` chico. Entregable: `informe/informe.pdf` (fuente `informe/informe.html`,
figuras en `informe/figures/`) + `out/provenance.json` de este trabajo.

Contenido:
1. Qué hace el repo y qué resultados verificados tiene (todo número sale de `data/paper_numbers.json`
   o de un script nuevo verificado).
2. Revisión contra SimuFLUX (Marin & Ries, Nat. Commun. 17:246, 2026): qué efectos realistas cubre
   SimuFLUX que este repo no, y cuáles importan para el caso objetivo.
3. Caso objetivo: p-MINFLUX (pulsos entrelazados, láser de 20 MHz, 4 exposiciones de igual duración,
   **no iterativo**, L fijo). Las variantes iterativas o de dwell central distinto quedan fuera.
4. Qué falta y qué se puede mejorar (priorizado).

Rondas: 1 (plazo: hoy). Equipo: 3 workers (A auditoría, B simulación de temporización p-MINFLUX,
C sesgo por fondo en TCP fijo) · verifier independiente · writer (orquestador).
Anti-objetivos: no tocar `tests/test_acceptance.py`; no escribir `data/paper_numbers.json` salvo vía
`compute_paper_numbers.py`; nada de trabajo no publicado de la autora; no subir PDFs de papers.
