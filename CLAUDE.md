# CLAUDE.md

<!-- agent-team:inicio (agent-team-template; no editar dentro de este bloque) -->
## Metodología agent-team

@AGENTS.md

## Específico de Claude Code

### Subagentes disponibles (`.claude/agents/`)
`pi`, `worker`, `verifier`, `code-reviewer`, `test-writer`, `writer` — el elenco de agent-team
adaptado a la sesión interactiva. El verificador nunca es la misma instancia que hizo el trabajo:
lánzalo como subagente aparte, con el trabajo a revisar pero **sin** el razonamiento del worker
(para que la comprobación sea independiente).

### Comandos
- `/equipo-nuevo <intent>` — crea un trabajo interactivo en `equipo/<id>/` (sin gasto extra).
- `/equipo-ronda <id> [indicación]` — corre UNA ronda: pi → workers → verifier → writer → checks.
- `/equipo-estado [<id>]` — resume el estado de los trabajos (interactivos y CLI).
- `/job-preparar <tipo> <descripción>` — ayuda a preparar un trabajo del CLI `job`: test de
  aceptación primero, intent corto, y el comando listo para que la persona lo confirme.
- El skill `agent-team` (`.claude/skills/`) sabe manejar el CLI `job` de este proyecto (`.\job.cmd` / `./job.sh`).

### Reglas de la sesión
- Una ronda por pedido. Al terminar una ronda, mostrar el resumen y **esperar** la decisión de la
  persona (continuar con indicación / congelar / abandonar). Nunca encadenar rondas solo.
- Si la tarea es chica y bien entendida, decirlo y ofrecer hacerla directo, sin equipo.
- El estado de un trabajo es de las máquinas (`state.json`); la persona lee el entregable en `out/`.
<!-- agent-team:fin -->

## Sobre este proyecto
<!-- Completa: qué es el proyecto, cómo se corre, cómo se testea, convenciones congeladas,
     límites de datos, qué NO tocar. Esto es lo que ningún test puede expresar. -->
- Descripción: simulación de localización de emisores fluorescentes con haces tipo dona
  (MINFLUX y variantes). Objetivo completo en `OBJECTIVE.md`; forma final tipo "paper
  companion" (paquete `src/donutloc`, `scripts/fig_*.py`, `tests/`, `data/paper_numbers.json`,
  `paper/` con procedencia).
- Entorno: Windows, Python 3.8 (numpy 1.24, scipy 1.10, matplotlib 3.7, numba opcional). Los
  scripts con caracteres no ASCII llevan `# -*- coding: utf-8 -*-`. No hay pdflatex local.
- Cómo correr los tests: `python -m unittest discover -s tests` (desde la raíz; el paquete se
  importa con `src/` en `sys.path` o `pip install -e .`).
- Convenciones: unidades en nm; dona LG I = 4e ln2 r²/fwhm² exp(-4 ln2 r²/fwhm²); TCP = 3 donas
  en círculo de diámetro L + centro; CRB = sqrt((Σxx+Σyy)/2); semilla 42.
- Fuentes: PDFs en `C:\Users\BANGHO\Documents\Doctorado\Papers` (NO subirlos); notas en
  `docs/literature/`.
- No tocar: `tests/test_acceptance.py` (protegido por hash en `equipo/*/state.json`), los PDFs
  de la carpeta de papers, el repo `p-minflux-main` de la autora (solo lectura).
