"""Sección 15. Inhibidores y activadores."""
from .celdas import code, md

# ============================================================================ 15. inhibición
md(r"""
## 15. Inhibidores y activadores

> 🎯 **En esta sección** distinguirás los tipos de inhibidores por la huella que dejan en
> V<sub>max</sub>, K<sub>M</sub> y el gráfico de Lineweaver–Burk, aprenderás tres maneras de medir
> K<sub>i</sub> y verás los inhibidores, reguladores y activadores reales de la glucoquinasa.

---

### 💡 La analogía: tres maneras de frenar la caja

Volvamos al cajero de la sección 13. Hay tres maneras de frenarlo:

* **Competitivo**: alguien se pone en la cola *sin comprar nada* y ocupa el turno. Con muchos
  clientes reales, el intruso pasa desapercibido: V<sub>max</sub> no cambia, pero hace falta más
  sustrato para llegar a ella (K<sub>M</sub> aparente sube).
* **Acompetitivo**: alguien traba la caja *solo cuando ya hay un cliente*: baja V<sub>max</sub>
  y, curiosamente, también K<sub>M</sub> (atrapa el complejo ES).
* **No competitivo**: alguien apaga la luz de la caja, haya o no cliente: la caja trabaja más
  lento (baja V<sub>max</sub>) pero los clientes entran igual (K<sub>M</sub> no cambia).

[[fig:inhibidores_cajero | Tres tipos de inhibidor representados en la caja del supermercado, con su efecto en Vmax, en KM y su huella en el gráfico de Lineweaver–Burk]]

---

### 🧮 La ecuación en palabras

**¿Para qué queremos una ecuación?** Para saber *cuánto* frena un inhibidor a una concentración
dada, y para medir su potencia con un solo número, K<sub>i</sub>, que no depende del experimento.

Para el competitivo la idea es: el intruso "roba" una parte de las cajas, así que hace falta
**más sustrato** para conseguir la misma ocupación. ¿Cuánto más? Un factor que crece con la
cantidad de intruso:

$$K_M^{\text{aparente}} \;=\; K_M \times \underbrace{\left(1 + \frac{\text{inhibidor}}{\text{potencia del inhibidor}}\right)}_{\text{cuántas veces más sustrato necesitas}}$$

Si hay tanto inhibidor como su K<sub>i</sub>, necesitas el doble de sustrato; si hay 10 veces su
K<sub>i</sub>, 11 veces más.

---

### 📐 La ecuación completa, término a término

**La derivación para el competitivo (las otras son análogas).** El inhibidor I se une a la
enzima libre: $\mathrm{E + I \rightleftharpoons EI}$ con $K_i = [\mathrm{E}][\mathrm{I}]/[\mathrm{EI}]$.
Ahora la enzima total es $[\mathrm{E}]_0 = [\mathrm{E}] + [\mathrm{ES}] + [\mathrm{EI}]$ (*cada caja está
libre, atendiendo a un cliente o bloqueada por el intruso*) y, repitiendo el álgebra de la
sección 13,

$$v_0 = \frac{V_{\max}[\mathrm{S}]}{K_M{\color{#e34948}{\left(1 + \dfrac{[\mathrm{I}]}{K_i}\right)}} + [\mathrm{S}]}
\qquad\Rightarrow\qquad K_M^{\mathrm{ap}} = K_M{\color{#e34948}{\left(1+\frac{[\mathrm{I}]}{K_i}\right)}},\; V_{\max}^{\mathrm{ap}} = V_{\max}.$$

El factor **rojo** (el color del inhibidor en el dibujo) es todo el efecto del intruso.

| Término | Qué es | En la analogía |
|---|---|---|
| [I] | concentración de inhibidor | cuántos intrusos hay en la tienda |
| $K_i$ | constante de disociación del complejo EI | la [I] que bloquea la mitad de las cajas libres: cuanto **menor**, más potente |
| 1 + [I]/K<sub>i</sub> (rojo) | factor de inhibición | cuántas veces más clientes necesitas para la misma ocupación |

| Tipo | Se une a | En la analogía | V<sub>max</sub> aparente | K<sub>M</sub> aparente | Huella en Lineweaver–Burk |
|---|---|---|---|---|---|
| competitivo | E | ocupa el turno | igual | × (1 + [I]/K<sub>i</sub>) | rectas que se cruzan en el eje y |
| acompetitivo | ES | traba la caja con cliente | ÷ (1 + [I]/K<sub>i</sub>′) | ÷ (1 + [I]/K<sub>i</sub>′) | rectas paralelas |
| no competitivo | E y ES por igual | apaga la luz | ÷ (1 + [I]/K<sub>i</sub>) | igual | rectas que se cruzan en el eje x |
| mixto | E y ES, distinto | un poco de cada | ÷ (1 + [I]/K<sub>i</sub>′) | × (1+[I]/K<sub>i</sub>)/(1+[I]/K<sub>i</sub>′) | se cruzan a la izquierda del eje y |

---

### 🎛️ Qué pasa si…

| Si… | entonces… | porque… |
|---|---|---|
| **[I] = K<sub>i</sub>** (competitivo) | la K<sub>M</sub> aparente **se duplica** | 1 + 1 = 2 |
| **[I] = 10 K<sub>i</sub>** (competitivo) | la K<sub>M</sub> aparente se multiplica por **11** | 1 + 10 = 11 |
| **añades mucho sustrato** | un competitivo **se vence**; un no competitivo, **nunca** | con cola larga el intruso no consigue turno; la luz apagada frena igual |

Elige el tipo de inhibidor y mueve [I] y K<sub>i</sub>:
""")

code(r'''
# @title 🎛️ Tipos de inhibidor y su huella
# 🎛️ Elige el tipo de inhibidor y mueve [I] y Ki: mira la curva y su huella en Lineweaver–Burk
interactivo.explorar_inhibicion()
''')

md(r"""
### 🔬 De los datos a K<sub>i</sub>: tres caminos que deben coincidir

1. **Ajuste global**: todos los puntos (varias [S] y varias [I]) contra la ecuación completa
   → V<sub>max</sub>, K<sub>M</sub> y K<sub>i</sub> con sus errores. Es el método recomendado hoy.
2. **Gráfico secundario**: se ajusta cada [I] por separado y se dibuja K<sub>M</sub><sup>ap</sup>
   frente a [I]: es una recta de pendiente K<sub>M</sub>/K<sub>i</sub>.
3. **Gráfico de Dixon**: 1/v₀ frente a [I] para varias [S]; en un competitivo las rectas se
   cruzan en [I] = −K<sub>i</sub>.

**De IC₅₀ a K<sub>i</sub>.** En farmacología se mide a menudo la **IC₅₀**: la concentración que
reduce la actividad a la mitad. En palabras: *cuanto más sustrato compite, más inhibidor hace
falta para frenar a la mitad*, así que la IC₅₀ depende de [S]. La ecuación de **Cheng–Prusoff**
descuenta ese efecto y devuelve K<sub>i</sub> (para un competitivo):

$$K_i = \frac{\mathrm{IC}_{50}}{1 + {\color{#eb6834}{[\mathrm{S}]/K_M}}}$$

El término **naranja** (el color del sustrato) es la "ventaja" que el sustrato le lleva al
inhibidor.
""")

code(r'''
# @title 🔬 Tres caminos para medir K_i
Ki_real = 3.0                                   # mM (inhibidor competitivo hipotético)
S_exp = np.array([1, 2, 4, 8, 16, 32, 64]); I_exp = (0, 2, 5, 10)
filas = []
for I in I_exp:
    v = cin.competitive_inhibition(S_exp, Vmax_real, Km_real, I, Ki_real) * (1 + 0.03 * rng.standard_normal(S_exp.size))
    filas += [(s, vv, I) for s, vv in zip(S_exp, v)]
df_i = pd.DataFrame(filas, columns=["S", "v", "I"])
# 1) ajuste global
aji = cin.fit_inhibition(df_i.S.values, df_i.v.values, df_i.I.values, kind="competitive")
# 2) gráfico secundario: Km aparente frente a [I]
km_ap = np.array([cin.fit_michaelis_menten(df_i[df_i.I == I].S.values, df_i[df_i.I == I].v.values)["km"] for I in I_exp])
sec = cin.secondary_plot_competitive(np.array(I_exp), km_ap, aji["km"])
# 3) Dixon
dx = cin.dixon_plot(df_i.S.values, df_i.v.values, df_i.I.values)
kd = cin.ki_from_dixon(sorted(dx["lines"]), dx, aji["km"], aji["vmax"])
# 4) IC50 -> Ki (Cheng–Prusoff)
S_ensayo = 8.0
ic50 = Ki_real * (1 + S_ensayo / Km_real)

fig, axes = viz.figure(14, 5.6, ncols=2)
colores = viz.sequential_blue(len(dx["lines"]))
for (s_val, (ii, inv_v)), color in zip(sorted(dx["lines"].items()), colores):
    fit = dx["fits"][s_val]
    xx = np.linspace(-1.3 * Ki_real, max(I_exp), 50)
    axes[0].plot(xx, fit["intercept"] + fit["slope"] * xx, color=color, lw=1.4, ls=(0, (4, 3)), zorder=2)
    xs = np.linspace(0, max(I_exp), 20)
    axes[0].plot(xs, fit["intercept"] + fit["slope"] * xs, color=color, lw=2.6, zorder=3)
    viz._markers(axes[0], ii, inv_v, color, label=f"[S] = {s_val:g} mM", size=8)
f_lo = dx["fits"][sorted(dx["lines"])[0]]
axes[0].axvline(0, color=viz.AXIS, lw=0.9, zorder=1); axes[0].axhline(0, color=viz.AXIS, lw=0.9, zorder=1)
viz._keypoint(axes[0], -kd["ki"], f_lo["intercept"] - f_lo["slope"] * kd["ki"], viz.COLORS["ts"], size=10)
axes[0].annotate(f"se cruzan en [I] = −Kᵢ\nKᵢ ≈ {kd['ki']:.2f} mM", (-kd["ki"], f_lo["intercept"] - f_lo["slope"] * kd["ki"]),
                 xytext=(4, -24), textcoords="offset points", ha="left", va="top", fontsize=11, color=viz.INK, fontweight="semibold")
viz._finish(axes[0], "[I] (mM)", "1/v₀ (s/µM)")
axes[0].set_ylim(bottom=-0.18 * axes[0].get_ylim()[1])
viz._legend(axes[0], loc="upper left", fontsize=10)
axes[0].set_title("Gráfico de Dixon", loc="left", fontsize=13, fontweight="semibold", pad=10)
xx = np.linspace(0, max(I_exp), 20)
axes[1].fill_between(xx, aji["km"], aji["km"] * (1 + xx / aji["ki"]), color=viz._tint(viz.PALETTE[7], 0.9), lw=0, zorder=1)
axes[1].plot(xx, aji["km"] * (1 + xx / aji["ki"]), color=viz.INK_SECONDARY, lw=2.4, zorder=2)
viz._markers(axes[1], I_exp, km_ap, viz.COLORS["datos"], size=9)
viz._guide(axes[1], "h", aji["km"])
axes[1].annotate("K_M sin inhibidor", (max(I_exp), aji["km"]), xytext=(-4, -6), textcoords="offset points", ha="right", va="top",
                 fontsize=10.5, color=viz.INK_SECONDARY)
axes[1].annotate(f"pendiente = K_M/Kᵢ\n→ Kᵢ = {sec['ki']:.2f} mM", (0.55 * max(I_exp), aji["km"] * (1 + 0.55 * max(I_exp) / aji["ki"])),
                 xytext=(-14, 10), textcoords="offset points", ha="right", va="bottom", fontsize=11, color=viz.INK, fontweight="semibold")
axes[1].set_ylim(0, None)
viz._finish(axes[1], "[I] (mM)", "K_M aparente (mM)")
axes[1].set_title("Gráfico secundario", loc="left", fontsize=13, fontweight="semibold", pad=10)
viz._fig_title(fig, f"Tres caminos, un mismo Kᵢ ≈ {aji['ki']:.1f} mM (valor usado para simular: {Ki_real:g} mM)",
               "Izquierda: 1/v₀ frente a [I] para varias [S]. Derecha: K_M aparente frente a [I]. Ambos con los mismos datos.")
tabla_dixon = df_i.rename(columns={"S": "[S]", "I": "[I]", "v": "v₀"})[["[I]", "[S]", "v₀"]].copy()
tabla_dixon["1/v₀"] = 1 / tabla_dixon["v₀"]
tabla_dixon["v₀ sin inhibidor (misma [S])"] = [df_i[(df_i.I == 0) & (df_i.S == sv)].v.iloc[0] for sv in tabla_dixon["[S]"]]
tabla_dixon["actividad restante"] = tabla_dixon["v₀"] / tabla_dixon["v₀ sin inhibidor (misma [S])"]
datos_dixon = viz.datos(
    tabla_dixon, "el gráfico de Dixon",
    "Los 28 tubos del experimento: 4 concentraciones de inhibidor × 7 de glucosa (3 % de ruido). En el gráfico de "
    "Dixon cada [S] es una recta: sus puntos son las filas con esa [S], con [I] en el eje x y 1/v₀ en el eje y. "
    "Con mucha glucosa el inhibidor competitivo casi no se nota (actividad restante cerca del 100 %).",
    x="[I]", y="1/v₀",
    calculadas={"1/v₀": "1 ÷ v₀", "actividad restante": "v₀ ÷ v₀ sin inhibidor, con la misma [S]"},
    unidades={"[I]": "mM", "[S]": "mM", "v₀": "µM/s", "1/v₀": "s/µM", "v₀ sin inhibidor (misma [S])": "µM/s"},
    formatos={"[I]": "{:g}", "[S]": "{:g}", "v₀": "{:.2f}", "1/v₀": "{:.3f}", "v₀ sin inhibidor (misma [S])": "{:.2f}",
              "actividad restante": "{:.0%}"},
    resaltar={i: ("[I] = 10 mM", "rojo") for i in tabla_dixon.index[(tabla_dixon["[I]"] == max(I_exp)) & (tabla_dixon["[S]"] <= 2)]},
    barra="1/v₀", max_filas=40, abierta=False)
tabla_sec = pd.DataFrame({"[I]": np.array(I_exp, dtype=float), "K_M aparente": km_ap, "K_M ap / K_M": km_ap / aji["km"],
                          "1 + [I]/Kᵢ (teoría)": 1 + np.array(I_exp) / aji["ki"]})
datos_sec = viz.datos(
    tabla_sec, "el gráfico secundario",
    "Una fila por concentración de inhibidor: a cada grupo de 7 tubos se le ajusta su propia hipérbola, y su K_M "
    "aparente es un punto del gráfico secundario. Si el inhibidor es competitivo, K_M ap/K_M crece como 1 + [I]/Kᵢ.",
    x="[I]", y="K_M aparente",
    calculadas={"K_M aparente": "ajuste de Michaelis–Menten solo a los tubos con esa [I]",
                "K_M ap / K_M": "K_M aparente ÷ K_M sin inhibidor (del ajuste global)",
                "1 + [I]/Kᵢ (teoría)": "lo que predice la ecuación competitiva con el Kᵢ del ajuste global"},
    unidades={"[I]": "mM", "K_M aparente": "mM"},
    formatos={"[I]": "{:g}", "K_M aparente": "{:.2f}", "K_M ap / K_M": "{:.2f}", "1 + [I]/Kᵢ (teoría)": "{:.2f}"},
    barra="K_M aparente")
viz.mostrar(fig, datos_dixon, datos_sec, viz.tarjetas([
    ("1 · ajuste global", f"{aji['ki']:.2f}", f"± {aji['ki_err']:.2f} mM", f"K_M = {aji['km']:.2f} mM, V_max = {aji['vmax']:.2f} µM/s", "azul"),
    ("2 · gráfico secundario", f"{sec['ki']:.2f}", "mM", f"R² = {sec['r2']:.3f}", "agua"),
    ("3 · Dixon", f"{kd['ki']:.2f}", "mM", "por el cruce de las rectas", "naranja"),
    ("4 · Cheng–Prusoff", f"{cin.cheng_prusoff(ic50, S_ensayo, Km_real):.2f}", "mM", f"desde IC₅₀ = {ic50:.1f} mM con [S] = {S_ensayo:g} mM", "violeta"),
], titulo="Kᵢ del mismo inhibidor, medido de cuatro maneras",
   nota="Si los métodos coinciden, el modelo (competitivo) describe bien los datos."))
''')

md(r"""
### 🔬 Los inhibidores y reguladores reales de la glucoquinasa

Tabla con fuentes. Tres cosas que sorprenden:

* La glucosa‑6‑fosfato, que frena a las hexoquinasas I–III, **no** inhibe a la glucoquinasa: por
  eso el hígado sigue fosforilando glucosa aunque el producto se acumule (Viñuela, Salas y Sols, 1963).
* La **manoheptulosa**, que muchos libros llaman "competitiva", resulta de tipo **mixto** cuando se
  mide con cuidado (Scruel 1998), y su K<sub>i</sub> varía 100 veces entre fuentes secundarias: por
  eso no damos un número.
* En el hígado, la **proteína reguladora GKRP** actúa como inhibidor competitivo frente a la glucosa
  (sube S<sub>0.5</sub> sin cambiar V; IC₅₀ ≈ 13 mM de glucosa) y secuestra a la enzima en el núcleo
  cuando la glucosa baja; la fructosa‑6‑fosfato refuerza ese secuestro y la fructosa‑1‑fosfato lo
  deshace.
""")

code(r'''
# @title 📋 Inhibidores y activadores reales de la glucoquinasa
inh = ref.get("inhibitors", [])
tipos = {"competitive": "competitivo", "mixed": "mixto", "none": "ninguno", "slow-binding": "unión lenta",
         "uncompetitive": "acompetitivo", "noncompetitive": "no competitivo"}
act = ref.get("activators", [])
viz.mostrar(
    viz.tabla(pd.DataFrame([(d.get("name"), tipos.get(d.get("kind"), d.get("kind")), d.get("versus", ""), d.get("ki", ""), d.get("note", ""), d.get("source", "")) for d in inh],
                           columns=["inhibidor / regulador", "tipo", "frente a", "Kᵢ / efecto", "nota", "fuente"]),
              titulo="Inhibidores y reguladores de la glucoquinasa"),
    viz.tabla(pd.DataFrame([(d.get("name"), d.get("fold", ""), d.get("ec50", ""), d.get("note", ""), d.get("source", "")) for d in act],
                           columns=["activador", "activación (veces)", "EC₅₀", "nota", "fuente"]),
              titulo="Activadores alostéricos (GKA)"),
)
''')

md(r"""
### 🔬 Activadores alostéricos (GKA)

En la analogía: alguien **aceita la caja** y abre el turno antes. Se unen a un sitio a ~20 Å del
sitio de glucosa (Kamata 2004), bajan S<sub>0.5</sub> y suben V<sub>max</sub>: lo contrario de un
inhibidor. El compuesto de referencia RO‑28‑1675 activa a la enzima silvestre ~16 veces con
EC₅₀ ≈ 7 µM (Sayed 2009) y se ensayó como fármaco para la diabetes tipo 2 (Grimsby 2003). ¿Qué
pasa con la actividad a 5 mM de glucosa cuando bajas S<sub>0.5</sub>? Muévelo:
""")

code(r'''
# @title 🎛️ Un activador alostérico
# 🎛️ Un activador alostérico: baja S0.5 y sube Vmax. ¿Qué pasa con la actividad a 5 mM de glucosa?
interactivo.explorar_activador()
''')

md(r"""
> ✅ **Para llevar.** Cada tipo de inhibidor deja una huella distinta: el competitivo sube
> K<sub>M</sub> (y se vence con sustrato), el acompetitivo baja V<sub>max</sub> y K<sub>M</sub>, el
> no competitivo baja solo V<sub>max</sub>. K<sub>i</sub> se mide por ajuste global, gráfico
> secundario o Dixon, y Cheng–Prusoff convierte una IC₅₀ en K<sub>i</sub>. Los activadores hacen lo
> contrario: bajan S<sub>0.5</sub>.
""")
