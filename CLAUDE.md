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
- Descripción:
- Cómo correr los tests:
- Convenciones:
- No tocar:
