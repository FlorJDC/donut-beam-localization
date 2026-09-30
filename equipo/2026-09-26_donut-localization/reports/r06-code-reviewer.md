# r06 — code-reviewer (consolidación, rama consolidate-r06)

Alcance: `git diff main -- scripts/ tests/ data/paper_numbers.json paper/ structure/`. No edité
código ni paper. Scripts auxiliares propios en `%TEMP%` (fuera del repo).

## Qué comprobé y números

**(c) Registro sin cambios fuera de lo pedido.** Comparé por programa `data/paper_numbers.json`
con `git show main:data/paper_numbers.json`. Ninguna clave desapareció. Las únicas claves nuevas son
`bgphys_max_abs_change_bias_x_nm` y `bgphys_max_rel_change_sigma_over_crb_pct`. Cambiaron de valor
exactamente 13 claves, todas `naive_*`: las 8 `mle_bias_abs` y `N_bias_eq_crb` de eps, las 4 de
SBR y `naive_honest_mle_max_abs_bias_nm`. Otras 4 claves (`naive_eps*_lms_extra_bias_x_*`)
conservan el valor pero cambian la `description`, porque heredan el texto `base`, que ahora
menciona el pulido Newton. Ninguna otra clave cambió en ningún campo.

**(a) `_polish_mle`.**
- Volví a correr `misspecified_estimator_section` y las 21 claves coinciden bit a bit con el JSON.
- Instrumenté los 12 casos: 8 ingenuos y honestos de eps, más 4 de SBR. En ninguno se usó el
  fallback. El desplazamiento desde el óptimo de la búsqueda por patrones es de 6e-6 a 7e-5 nm. La
  log-verosimilitud sube siempre (Δll de 5e-12 a 1e-10). |∇ll| baja de ~1e-6 a ~1e-9.
- El Newton independiente del verificador (`work/verify/r06/v_polish.py`, gradiente por paso
  complejo) da 0.2687115839 y 0.4220558338 nm. El código da 0.2687115970 y 0.4220558629 nm, o sea
  diferencias de 1e-8 a 3e-8 nm. El Nelder-Mead de `deterministic.json` queda a unos 4e-7 nm y es
  la referencia menos precisa.
- El MLE honesto baja a 5.7e-8 nm. En el paper se cita solo como "< 10⁻⁴ nm", que es correcto.
- Ramas: `eigvalsh(H) < 0` estricto evita resolver con un Hessiano singular. Un paso fuera del
  disco devuelve `start`, y agotar `max_iter` también devuelve `start`. En los casos reales la
  búsqueda arranca a unos 1e-5 nm del máximo y no hay riesgo de saltar a otro máximo local. La
  función no verifica que ll(final) ≥ ll(start); es aceptable con estos datos, pero no está
  garantizado en general.

**(b) Claves nuevas.** Las recalculé desde las 24 claves `bgphys_{fixed,phys}_{mle,lms,mlms}_{bias_x,sigma_over_crb}_x{25,50}`:
- 1.2216 nm, que corresponde a mLMS en x0=50.
- 7.127 %, que corresponde a LMS en x0=50.

Las dos coinciden con el JSON. Solo comparan el sesgo x y usan x0 = 25 y 50, como dice la
descripción.

**Tests.** `TestNoiseFreePolish` pasa 4/4. La suite completa (`python -m unittest discover -s tests`)
da OK con 1 skipped. Todos los `\src{}` de `paper/sections/*.tex` resuelven en
`paper/provenance.json`.

**(d) Texto.**
- "independent of N" y "does not decrease with N" pasaron a "large-N limit"; "at finite N its mean
  bias is larger" coincide con el MC del verificador.
- La discusión ya declara los puntos comparados y que solo se compara la componente x.
- `iterative.tex` ahora dice que la rampa de potencia "partially compensates" y declara el supuesto
  sobre el fondo. Esto resuelve los puntos 1 a 4 del reporte r06 del verificador.
- En `open_points.tex` hay 2 ítems nuevos coherentes con la decisión de la autora en el inbox.
- No quedaron números viejos (0.4220, 11173, 3.05e-5, "1.3 nm", "8 %") en paper, README ni claims.

## Hallazgos

1. **Cobertura de tests engañosa en las ramas de fallback.**
   `test_step_leaving_disk_returns_start` arranca en (0.5, 0) con el emisor en (10, 0). Ahí el
   Hessiano por diferencias finitas tiene autovalores (−0.1436, **+0.00062**): la función vuelve
   por la rama "Hessiano no negativo definido", no por la de salida del disco. La rama del disco
   **no está ejercitada**. Existe una entrada que sí la ejercita: `start=(9.95,0)`,
   `center=(9.9,0)`, `radius=0.06` devuelve `start`, y con `radius=0.5` converge a (10,0).
   `test_no_convergence_returns_start` usa `max_iter=0`: el bucle nunca corre, así que no prueba
   una no-convergencia real. La rama del Hessiano sí queda cubierta, pero por accidente.
   Sin efecto en los números; es una laguna de tests.
2. **Precisión impresa de más en las claves nuevas (menor).** 1.222 nm y 7.127 % son máximos de
   diferencias de MC y no tienen SE registrado. El SE de las diferencias que los determinan es de
   ~0.045 nm y ~0.75 puntos porcentuales si se suman como independientes (misma semilla, así que
   la cota es conservadora). El MC propio del verificador dio 1.26 nm y 6.3 %. El texto dice "at
   most 1.222 nm / 7.127 %", con 4 cifras, como si fuera una cota exacta. Conviene 1.2 nm y 7 %,
   como `structure/claims.json`, que ya dice 1.22 y 7.1.
3. **Sin test para las claves `bgphys_max_*`.** Las comprobé a mano. Nada impide que la cuenta se
   rompa en silencio, por ejemplo con un cambio de claves en `BGPHYS_X`.
4. Menor: `radius=np.inf` por defecto en `_polish_mle` deja sin red de seguridad un Hessiano casi
   singular. Solo se llama con radio finito, así que no hay escenario de falla actual.

## Afirmaciones vivas previas
Las refuted/unclear de código de r01-r02 no entran en el alcance de esta ronda y no las
re-examiné.
