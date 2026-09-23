"""Sección 9. TS."""
from .celdas import code, md

# ============================================================================ 9. TS

md(r"""
## 9. El estado de transición

> 🎯 **En esta sección** localizarás el estado de transición de verdad (un punto de silla), lo
> **validarás** con sus frecuencias de vibración y comprobarás que conecta el reactivo con el
> producto.

---

### 💡 La analogía: el paso de montaña, otra vez

El escaneo nos dio una colina aproximada, pero forzando una sola distancia. El verdadero estado
de transición es un **punto de silla**: un máximo a lo largo del camino y un mínimo en todas las
demás direcciones, exactamente como el paso entre dos montañas (si te sales del camino, subes).

[[fig:punto_silla | Mapa de curvas de nivel con dos valles y un collado; las cuentas del NEB siguen el camino y la central trepa al punto de silla]]

**Cómo se localiza (tres pasos):**

1. **NEB con imagen trepadora** (*nudged elastic band*): una cadena de geometrías entre reactivo
   y producto, unidas por muelles, que se relaja hasta dibujar el camino de mínima energía. La
   imagen más alta "trepa" hasta la cima. Es el análogo moderno del método **QST2** (dos
   extremos → TS) de Gaussian y MOPAC.
2. **Método del dímero**: refina la cima usando solo gradientes, hasta que la fuerza es cero.
3. **Frecuencias**: en un punto de silla verdadero hay **exactamente una** frecuencia imaginaria,
   y su vector describe el movimiento de la reacción. Es la prueba definitiva.

---

### 🧮 La prueba en palabras: ¿por qué "imaginaria"?

**¿Para qué queremos una ecuación?** En la cima la fuerza es cero… pero también lo es en el fondo
de un valle. Para distinguirlos miramos **la forma** del terreno alrededor, y eso lo dan las
vibraciones. Cada vibración es como un **resorte**: cuanto más empinadas las paredes, más rápido
oscila.

$$\text{frecuencia} \;\propto\; \sqrt{\text{rigidez del resorte}}$$

En un valle todas las paredes suben: todos los resortes son normales. En la silla, a lo largo del
camino el terreno **baja** a los dos lados: es un "resorte al revés" (rigidez negativa), y la raíz
cuadrada de un número negativo es un número **imaginario**. Los programas lo escriben como una
frecuencia negativa.

[[fig:frecuencia_imaginaria | En un valle todas las direcciones suben; en la silla una dirección baja y da una frecuencia imaginaria]]

---

### 📐 La ecuación completa, término a término

$$\nu \;=\; \frac{1}{2\pi}\sqrt{\frac{{\color{#eb6834}{k}}}{\mu}}$$

| Término | Qué es | En la analogía | En el estado de transición |
|---|---|---|---|
| $\nu$ | frecuencia de la vibración (en cm⁻¹ en química) | cuántas veces por segundo oscila la bola | −149 cm⁻¹ en el modo de la reacción |
| ${\color{#eb6834}{k}}$ | curvatura de la energía en esa dirección | la rigidez del resorte: cuán empinadas son las paredes | **negativa** a lo largo del camino |
| $\mu$ | masa efectiva que se mueve | cuán pesada es la bola | sobre todo el fósforo y sus oxígenos |

---

### 🎛️ Qué pasa si…

| Si al calcular las frecuencias sale… | entonces la geometría es… |
|---|---|
| **ninguna** frecuencia imaginaria | un **mínimo** (reactivo, producto o un intermedio), no un TS |
| **exactamente una** | un **punto de silla de primer orden**: un estado de transición ✔ |
| **dos o más** | una cima en varias direcciones (silla de orden superior): hay que seguir buscando |

---

### 🔬 Los datos reales
""")

code(r'''
# @title ⛰️ Del escaneo al camino de mínima energía (NEB)
neb = datos.csv("qmmm/neb.csv")
perfiles = [(xi_esc, e_esc, "escaneo restringido"),
            (np.r_[R["xi"], neb["xi"].values[1:]], np.r_[0.0, neb["energia_rel_kcal"].values[1:]], "NEB (imagen trepadora)")]
fig = viz.plot_energy_profiles(perfiles, xlabel="ξ (Å)", relative=False, colors=[viz.INK_MUTED, viz.COLORS["reactivo"]],
                               title="El NEB encuentra la misma cima por un camino más natural",
                               subtitle="Gris: escaneo que fuerza una sola distancia. Azul: NEB, que relaja todas las coordenadas a la vez")
fig;
''')

code(r'''
# @title ✅ El estado de transición y su frecuencia imaginaria
ts = res["etapas"].get("dimero")
if ts is None or not ts.get("convergido", True):
    display(viz.mensaje("El dímero no convergió: se usa la imagen trepadora del NEB como TS.", "ojo"))
    ts = ts or {**res["etapas"]["neb"], "xi": neb.loc[res["etapas"]["neb"]["imagen_ts"], "xi"], "d_PG_O6": float("nan"), "d_PG_O3B": float("nan")}
fr = res["etapas"].get("frecuencias", {})
val = fr["validacion"]
freqs = np.asarray(fr["frecuencias_mas_bajas"])
tarjetas_ts = viz.tarjetas(
    [("Barrera ΔE‡", f"{ts['barrera_kcal']:.1f}", "kcal/mol", "desde el reactivo relajado (método del dímero)", "naranja"),
     ("ξ en la cima", f"{ts['xi']:+.2f}", "Å", "casi cero: a medio camino", "naranja"),
     ("Pγ ··· O6 / Pγ ··· O3β", f"{ts['d_PG_O6']:.2f} / {ts['d_PG_O3B']:.2f}", "Å", "el fósforo, equidistante de los dos oxígenos", "azul"),
     ("Frecuencias imaginarias", f"{val['imaginary_mode_count']}", "", "punto de silla de primer orden ✔" if val["ok"] else "revisar ✘",
      "verde" if val["ok"] else "rojo")],
    titulo="El estado de transición, localizado y validado")
fig = viz.plot_frequencies([("TS en la enzima", freqs)],
                           title=f"Una sola frecuencia imaginaria ({freqs.min():.0f} cm⁻¹): es un punto de silla".replace("-", "−"),
                           subtitle="Las seis vibraciones más lentas del estado de transición; las demás son reales (positivas)")
viz.mostrar(tarjetas_ts, fig)
''')

md(r"""
**Qué nos dice la geometría del estado de transición.** En el TS el fósforo está a la misma
distancia de los dos oxígenos (≈ 2.1–2.2 Å de O6 y de O3β): es un TS **concertado y "en
línea"**; el fosfato pasa de un oxígeno al otro invirtiendo sus tres oxígenos como un paraguas.
El protón del O6, en cambio, **todavía no se ha movido** (sigue a ~1.07 Å de O6). Solo después
de cruzar la cima, cuando el enlace P–O6 ya está formado, el protón salta a Asp205: el camino
baja por una "meseta" hasta ξ ≈ 1 Å y allí termina de caer. En este modelo Asp205 actúa de base
**tras** la transferencia del fosforilo. Un estudio QM/MM independiente (Zhang et al., 2009)
describe a Asp205 como base general y a Lys169 como ácido general que protona el fosfato; que
dos modelos discrepen en el *orden* de los pasos, pero coincidan en la barrera (sección 12), es
un buen ejemplo de por qué los mecanismos se estudian con simulaciones y con mutantes.

**El movimiento del estado de transición.** La frecuencia imaginaria corresponde a una vibración
que no oscila sino que "cae" hacia reactivo o hacia producto. Animándola vemos la química: mira
cómo el fósforo va y viene entre el O3β del ATP y el O6 de la glucosa, y cómo el punto de la
gráfica cruza de la zona «hacia el reactivo» a la zona «hacia el producto»:
""")

code(r'''
# @title 🎞️ El movimiento del estado de transición
visor3d.modo_imaginario()
''')

md(r"""
**Bajar de la cima.** Si empujamos el TS un poquito a cada lado a lo largo del modo imaginario
y dejamos que la geometría relaje, descendemos por el camino de mínima energía hasta el
reactivo y hasta el producto. Esto confirma que **este** TS conecta **estos** reactivos con
**estos** productos (la idea del IRC, *intrinsic reaction coordinate*).
""")

code(r'''
# @title ⬇️ Bajar de la cima hacia reactivo y producto
camino = datos.csv("qmmm/camino_descenso.csv")
fig = viz.plot_energy_profile(camino["xi"].values, camino["energia_rel_kcal"].values, relative=False, xlabel="ξ (Å)",
                              ts_index=int(camino["energia_rel_kcal"].idxmax()), smooth=False, sort=False, annotate_states=True,
                              state_names=("reactivo", "TS", "fin del descenso"),
                              title="Desde la cima se baja a los dos valles: este TS conecta R con P",
                              subtitle="Cada punto es una geometría relajada; la meseta tras la cima es el protón que pasa a Asp205")
fig;
''')

code(r'''
# @title 📊 Diagrama de energía de la reacción
niveles = [("E·S (reactivo)", 0.0), ("TS", ts["barrera_kcal"]), ("E·P (producto)", P["dE_reaccion_kcal"])]
fig = viz.plot_energy_levels(niveles, ts_indices=[1], title="La reacción en tres números",
                             subtitle="PM7/Amber, entorno fijo: energías electrónicas relativas (ΔE), todavía no energías libres")
fig;
''')

md(r"""
> ✅ **Para llevar.** El estado de transición es un **punto de silla**: se localiza en pasos
> (escaneo → NEB → dímero) y se **valida** con dos pruebas: exactamente una frecuencia imaginaria
> (aquí −149 cm⁻¹, el fósforo que va y viene) y un camino que, al bajar de la cima, llega a los
> reactivos por un lado y a los productos por el otro.
""")
