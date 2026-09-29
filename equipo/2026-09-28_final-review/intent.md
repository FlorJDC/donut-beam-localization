# Intent — revisión final + explorador guiado

Pedido de la autora (2026-09-28, 22:55): revisión final de donut-beam-localization con la metodología
agent-team, mejoras en una rama nueva (`final-review-and-explorer`), y un set de scripts ejecutables:
un explorador standalone y guiado (importa `donutloc`) para aprender de las simulaciones, más una guía
HTML en la misma carpeta que explique cada modo y qué esperar. Merge a main solo si se termina antes de
las 23:50 (hora de Argentina).

Equipo: code-reviewer (revisión final del repo, diff real) · worker = orquestador (explorador + guía +
arreglos) · verifier (corre el explorador en todos los modos y contrasta las salidas con los números
verificados). Una ronda. Anti-objetivos: no tocar tests/test_acceptance.py; nada de trabajo no publicado.
