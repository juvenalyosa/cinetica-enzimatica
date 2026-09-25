"""Sección 17. resumen."""
from ..celdas import code, md

# ============================================================================ 17. resumen
md(r"""
## 17. Resumen, glosario y ejercicios

> 🎯 **En esta sección** juntamos todas las piezas: el mapa del curso, las ecuaciones en una sola
> página, un glosario y ejercicios para comprobar lo aprendido.

---

### 🗺️ El mapa del curso

[[fig:mapa_curso | Mapa conceptual: estructura, dinámica molecular, QM/MM, de ΔE‡ a ΔG‡, Eyring, Michaelis–Menten, Hill e inhibidores, con temperatura y pH como condiciones]]

**Lo que hemos visto, en una frase cada cosa:**

1. Una enzima acelera una reacción **bajando la barrera** del estado de transición, como un guía
   que conoce el paso más bajo.
2. La estructura cristalina (3FGU) muestra a los reactivos ya alineados; la MD muestra que la
   enzima **los mantiene alineados**.
3. QM/MM permite **ver la química** dentro de la enzima; hay que cuidar la partición, los átomos
   de enlace, la electrostática y la repulsión QM–MM.
4. El estado de transición es un **punto de silla**: se localiza (escaneo → NEB → dímero) y se
   **valida** con una única frecuencia imaginaria y con el camino que baja hacia R y P.
5. De ΔE‡ a ΔG‡: vibraciones, promedio sobre conformaciones y sensibilidad al método. Un
   resultado computacional lleva **incertidumbre**, y nuestra barrera coincide con otro cálculo
   publicado y sobrestima el experimento en ~3–4 kcal/mol.
6. Eyring convierte la barrera en *k*<sub>cat</sub>; 1.4 kcal/mol son un factor 10.
7. Michaelis–Menten **emerge** del mecanismo E + S ⇌ ES → E + P; K<sub>M</sub>, V<sub>max</sub>,
   k<sub>cat</sub> y k<sub>cat</sub>/K<sub>M</sub> se derivan de los datos por ajuste no lineal.
8. La glucoquinasa es **sigmoide** (n ≈ 1.7) sin tener varias subunidades: cooperatividad
   cinética; sus mutantes desplazan el interruptor de la insulina.
9. Cada tipo de inhibidor deja una huella distinta; K<sub>i</sub> se obtiene por ajuste global,
   gráfico secundario o Dixon; los activadores hacen lo contrario.
10. Temperatura y pH modulan k a través de ΔH‡, ΔS‡ y los pK<sub>a</sub> del sitio activo.

---

### 🧮 Las ecuaciones del curso, en una página

**Velocidad inicial** (sección 2): la velocidad es la pendiente inicial de la curva de producto.

$$v_0 = \left(\frac{\Delta[\mathrm{P}]}{\Delta t}\right)_{t\to 0}$$

**Eyring** (secciones 6 y 12): velocidad = intentos por segundo × probabilidad de llegar a la cima.

$$k = \kappa\,\frac{k_B T}{h}\,\exp\!\left(-\frac{\Delta G^{\ddagger}}{RT}\right)$$

**Michaelis–Menten** (sección 13): la caja se satura; con mucho sustrato la enzima trabaja a tope.

$$v_0 = \frac{V_{\max}[\mathrm{S}]}{K_M + [\mathrm{S}]}$$

**Hill** (sección 14): un interruptor; poca actividad hasta un umbral, luego mucha.

$$v_0 = \frac{V_{\max}[\mathrm{S}]^{n}}{S_{0.5}^{\,n} + [\mathrm{S}]^{n}}$$

**Inhibición competitiva** (sección 15): el intruso exige más sustrato para la misma velocidad.

$$K_M^{\mathrm{ap}} = K_M\left(1 + \frac{[\mathrm{I}]}{K_i}\right)$$

**Arrhenius** (sección 16): calentar acelera más cuanto más alta es la barrera.

$$\ln k = \ln A - \frac{E_a}{RT}$$

**pH** (sección 16): la enzima solo funciona si cada grupo del sitio activo tiene la carga correcta.

$$v = \frac{v_{\max}}{1 + 10^{pK_1 - \mathrm{pH}} + 10^{\mathrm{pH} - pK_2}}$$

---

### 📖 Glosario rápido

| Término | Significado |
|---|---|
| *Sustrato* | molécula transformada |
| *Sitio activo* | donde ocurre la química |
| *ES* | complejo enzima‑sustrato |
| *Estado de transición* | cima de la barrera |
| *ΔG‡* | altura de la barrera |
| *k*<sub>cat</sub> | reacciones por segundo por enzima |
| *K*<sub>M</sub> | [S] a media velocidad |
| *k*<sub>cat</sub>/*K*<sub>M</sub> | eficiencia a baja [S] |
| *S*<sub>0.5</sub>, *n* | parámetros de Hill |
| *K*<sub>i</sub> | constante de disociación del inhibidor |
| *QM/MM* | mecánica cuántica para la región reactiva y campo de fuerza para el resto |
| *NEB, dímero* | métodos para hallar el TS |
| *Frecuencia imaginaria* | la firma de un punto de silla |

---

### ✍️ Ejercicios

**Para empezar**

1. Mueve el deslizador `umbral` de la sección 5 a 3.0 y a 4.0 Å y vuelve a ejecutar la celda. ¿Cómo cambia la fracción de conformaciones de
   ataque cercano? ¿Qué pasaría con *k*<sub>cat</sub> si la enzima no cerrara sus dominios?
2. Con el explorador de Eyring, averigua cuánto tendría que bajar la barrera para multiplicar
   *k*<sub>cat</sub> por 1000. Compáralo con la diferencia entre el sitio activo en agua y en la enzima.
3. En la sección 13, mueve el deslizador `k2` a 6 s⁻¹ y a 600 s⁻¹ y vuelve a ejecutar la celda. ¿Cómo cambian K<sub>M</sub> y V<sub>max</sub>?
   ¿Cuándo K<sub>M</sub> ≈ K<sub>d</sub> = k<sub>−1</sub>/k<sub>1</sub>?

**Para profundizar**

4. Simula datos con n = 1.0, 1.4 y 2.0 en la sección 14 y ajústalos con Michaelis–Menten.
   ¿Qué error cometes si ignoras la cooperatividad? ¿Qué mutante de la tabla dejaría de liberar
   insulina a 5 mM?
5. Diseña un experimento (concentraciones de S e I) que distinga un inhibidor competitivo de uno
   mixto con K<sub>i</sub>′ = 3K<sub>i</sub>. Para comprobarlo, añade una celda nueva (botón «+ Código») y usa `cin.fit_inhibition` y el gráfico de Dixon.
6. Con la termoquímica de la sección 10, ¿cuánto cambia ΔG‡ entre 25 y 37 °C? ¿Qué factor en
   *k* supone? Compáralo con lo que predice el explorador de temperatura.

**Reto (avanzado, modo `completo`)**

7. En `enzimas/glucoquinasa.py` saca Lys169 de `QM_RESIDUES` (quedará como cargas MM) y
   recalcula el escaneo. ¿Sube o baja la barrera? Es un "mutante computacional"; Zhang et al.
   (2009) obtuvieron 32 kcal/mol para K169A frente a 18 para la silvestre.

---

### 📚 Referencias
""")

code(r'''
# @title 📚 Referencias
from IPython.display import HTML
import html as _html
refs = ref.get("references", [])
items = "".join(f'<li style="margin:0 0 8px 0;padding-left:4px">{_html.escape(r)}</li>' for r in refs)
viz.mostrar(HTML(
    f'<div style="font-family:Figtree,\'Avenir Next\',\'Segoe UI\',Roboto,Arial,sans-serif;background:{viz.SURFACE};'
    f'border:1px solid {viz.GRID};border-radius:18px;padding:18px 22px;max-width:960px;color:{viz.INK_SECONDARY};font-size:14px;line-height:1.45">'
    f'<div style="font-size:16px;font-weight:650;color:{viz.INK};margin-bottom:10px">📚 Fuentes de los datos del curso</div>'
    f'<ol style="margin:0;padding-left:22px">{items}</ol>'
    f'<div style="margin-top:12px;font-size:13px;color:{viz.INK_MUTED}">Herramientas QM/MM adaptadas de la suite Leonardo (Juvenal Yosa, MIT): '
    f'<a href="https://github.com/juvenalyosa/Leonardo" style="color:{viz.PALETTE[0]}">github.com/juvenalyosa/Leonardo</a></div></div>'))
''')
