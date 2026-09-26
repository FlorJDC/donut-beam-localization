# agent-team-template

Plantilla de arranque para proyectos que trabajan con la metodología de
[`agent-team`](https://github.com/matiaszaldarriaga/agent-team): un equipo acotado de agentes
(Claude Code o Codex) con un **verificador independiente permanente**, que produce **un**
entregable con **procedencia reproducible** y se detiene solo. Tú decides si continuarlo,
congelarlo o abandonarlo.

Trae el toolkit ya integrado (funciona en Windows, Linux y macOS, con Python ≥ 3.8), el elenco
de agentes como subagentes de Claude Code, comandos para trabajar en modo interactivo, y la
estructura mínima para que un proyecto nuevo empiece a andar en minutos.

## Estructura

```
OBJECTIVE.template.md   plantilla del objetivo: cópiala a OBJECTIVE.md y complétala
papers/                 material fuente que el equipo puede leer (papers, specs, datos)
AGENTS.md               la metodología (la leen Claude y Codex)
CLAUDE.md               reglas para Claude Code (importa AGENTS.md) + "Sobre este proyecto"
job.cmd / job.sh        lanzadores del CLI `job` para este proyecto (Windows / bash)
agent-team/             el toolkit agent-team (motor, roles, recetas, tests); ver agent-team/README.md
.agent-team/            configuración propia del proyecto:
  policy.json             límites a lo que el PI puede elegir
  roles/  recipes/        roles y tipos de trabajo propios (tienen prioridad sobre los del toolkit)
  plantillas/             intent de `feature`, test de aceptación, checks.py, checks.sh
  intents/                los intents de cada trabajo
.claude/agents/         el elenco como subagentes: pi, worker, verifier, code-reviewer, test-writer, writer
.claude/commands/       /equipo-nuevo  /equipo-ronda  /equipo-estado  /job-preparar
.claude/skills/         skill agent-team (Claude sabe manejar `job` en este proyecto)
jobs/                   corridas del CLI `job` (gitignored: es estado local)
equipo/                 trabajos del modo interactivo (gitignored)
scripts/                apply-to-existing.ps1 / .sh: llevar esta plantilla a un proyecto existente
```

## Puesta en marcha de un proyecto nuevo

1. Crea el repo del proyecto desde esta plantilla: botón **"Use this template"** en GitHub,
   o clónala:

   ```sh
   git clone https://github.com/FlorJDC/agent-team-template mi-proyecto
   cd mi-proyecto && git remote remove origin      # que el proyecto tenga su propio remoto
   ```

2. Requisitos (una vez por máquina): Python ≥ 3.8, Git (en Windows, Git for Windows: su `bash`
   corre los checks), y los CLIs `claude` y/o `codex` para los backends reales. **No hay nada
   que instalar**: `job.cmd` / `job.sh` usan el toolkit de `agent-team/`.

   ```powershell
   .\job.cmd recipes                  # Windows (cmd / PowerShell)
   ./job.sh recipes                   # bash (Linux, macOS, Git Bash)
   ```

3. Completa la sección **"Sobre este proyecto"** de `CLAUDE.md` (cómo correr tests,
   convenciones, qué no tocar) y el objetivo:

   ```sh
   cp OBJECTIVE.template.md OBJECTIVE.md
   ```

4. Pon el material fuente en `papers/` (o lo que corresponda).

5. Trabaja, en cualquiera de los dos modos.

## Dos modos, la misma metodología

### Modo interactivo (dentro de Claude Code)

Para seguir y entender cada paso. Claude hace de motor: lanza los subagentes, cosecha las
afirmaciones del verificador en `equipo/<id>/state.json`, corre checks y aceptación, y se
detiene después de cada ronda.

```
/equipo-nuevo derive @OBJECTIVE.md          crea equipo/<id>/ (sin correr nada)
/equipo-ronda <id>                          UNA ronda: pi → workers → verifier → writer → checks
/equipo-ronda <id> "revisa el caso límite"  otra ronda, con tu indicación
/equipo-estado                              estado de todos los trabajos
```

### Modo CLI (`job`: el motor completo)

Para trabajos más largos o desatendidos: rondas acotadas, tripwires de progreso, checkpoint
a las 2 rondas, `view.html` para monitorear e inyectar indicaciones.

```powershell
.\job.cmd new derive "@OBJECTIVE.md" --pi --run      # el PI dota el equipo y corre (gasto real)
start jobs\<id>\view.html                            # monitor + caja para inyectar directivas
.\job.cmd resume <id> --say "concéntrate en el canal resonante"
.\job.cmd freeze <id>                                # o: abandon <id>
```

`derive` es una de varias recetas: también hay `feature`, `draft`, `wiki` (y las propias en
`.agent-team/recipes/`). `job roles` y `job recipes` listan lo disponible. En Claude Code,
`/job-preparar` te ayuda a escribir primero el test de aceptación y un intent corto, y arma el
comando. Probar sin gastar: `.\job.cmd new derive "prueba" --backend mock --run --rounds 2`.

## La metodología en 8 reglas

1. **Empieza simple**: si es chico y bien entendido, hazlo directo. El equipo es para cuando la
   verificación independiente *es* el producto.
2. **Primero el test de aceptación, en rojo.** Es la especificación; el equipo puede hacerlo
   pasar, nunca editarlo (queda protegido por hash).
3. **Objetivo corto**: solo lo que un test no puede expresar.
4. **El verificador es permanente** y trabaja sin el razonamiento del worker. Solo lo que él
   reproduce pasa a `verified`, y lo verificado no se re-deriva.
5. **Nada se inventa**: cada afirmación lleva `\src{clave}` / `[src:clave]` respaldada en
   `out/provenance.json`; el check de procedencia bloquea "terminado" si falta algo.
6. **Pocas rondas** (3), leer, y continuar con una indicación.
7. **Nunca un scheduler**: todo acotado por rondas y presupuesto, con kill-switch.
8. **Tú decides** al final de cada corrida. Integrar el resultado al proyecto es un paso aparte.

Detalle completo en [`AGENTS.md`](AGENTS.md), [`agent-team/README.md`](agent-team/README.md) y
el porqué en [`agent-team/docs/DESIGN.md`](agent-team/docs/DESIGN.md).

## Usarla en un proyecto que ya existe

Desde un clon de esta plantilla:

```powershell
.\scripts\apply-to-existing.ps1 -Target C:\ruta\al\proyecto -Git
```
```sh
./scripts/apply-to-existing.sh /ruta/al/proyecto --git
```

No pisa nada: si el proyecto ya tiene `CLAUDE.md`, `AGENTS.md` o `.gitignore`, solo agrega (o
actualiza al re-ejecutar) una sección marcada `<!-- agent-team:inicio -->…<!-- agent-team:fin -->`.

## Actualizar el toolkit

`agent-team/` es una copia del toolkit upstream (commit `87e1368`) con parches para Windows /
Python 3.8 y soporte de configuración por proyecto (`.agent-team/`). Por eso va copiado y no
como submódulo: el upstream todavía no corre en Windows con Python 3.8. Qué cambió y cómo
traer cambios del upstream: [`agent-team/PARCHES-PLANTILLA.md`](agent-team/PARCHES-PLANTILLA.md).
Tests: `cd agent-team && python -m unittest discover -s tests`.

## Qué NO incluye esta plantilla a propósito

Cosas específicas de un proyecto (setup del entorno, el paper puntual, los `OBJECTIVE.md` ya
completados, las corridas en `jobs/` y `equipo/`) no van en la plantilla: nacen en cada
proyecto concreto. Esta plantilla es el esqueleto reutilizable, no un proyecto en sí.

## Licencia

MIT. El toolkit `agent-team/` es © sus autores ([`agent-team/LICENSE`](agent-team/LICENSE)).
