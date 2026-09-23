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
Ki_real = 3.0                                   # mM (inhibidor competitivo hipotético)
S_exp = np.array([1, 2, 4, 8, 16, 32, 64]); I_exp = (0, 2, 5, 10)
filas = []
for I in I_exp:
    v = cin.competitive_inhibition(S_exp, Vmax_real, Km_real, I, Ki_real) * (1 + 0.03 * rng.standard_normal(S_exp.size))
    filas += [(s, vv, I) for s, vv in zip(S_exp, v)]
df_i = pd.DataFrame(filas, columns=["S", "v", "I"])
# 1) ajuste global
aji = cin.fit_inhibition(df_i.S.values, df_i.v.values, df_i.I.values, kind="competitive")
print(f"1) Ajuste global:      K_i = {aji['ki']:.2f} ± {aji['ki_err']:.2f} mM   (K_M = {aji['km']:.2f} mM, V_max = {aji['vmax']:.2f} µM/s)")
# 2) gráfico secundario: Km aparente frente a [I]
km_ap = np.array([cin.fit_michaelis_menten(df_i[df_i.I == I].S.values, df_i[df_i.I == I].v.values)["km"] for I in I_exp])
sec = cin.secondary_plot_competitive(np.array(I_exp), km_ap, aji["km"])
print(f"2) Gráfico secundario: K_i = {sec['ki']:.2f} mM (R² = {sec['r2']:.3f})")
# 3) Dixon
dx = cin.dixon_plot(df_i.S.values, df_i.v.values, df_i.I.values)
kd = cin.ki_from_dixon(sorted(dx["lines"]), dx, aji["km"], aji["vmax"])
print(f"3) Dixon:              K_i = {kd['ki']:.2f} mM (por el cruce de las rectas)")
fig, axes = viz.figure(11, 4, ncols=2)
for (s_val, (ii, inv_v)), color in zip(sorted(dx["lines"].items()), viz.sequential_blue(len(dx["lines"]))):
    fit = dx["fits"][s_val]
    xx = np.linspace(-1.3 * Ki_real, max(I_exp), 50)
    axes[0].plot(xx, fit["intercept"] + fit["slope"] * xx, color=color, lw=1.5)
    axes[0].plot(ii, inv_v, "o", color=color, ms=5, label=f"[S] = {s_val:g} mM")
axes[0].axvline(-Ki_real, color=viz.INK_MUTED, ls="--", lw=1); axes[0].set_xlabel("[I] (mM)"); axes[0].set_ylabel("1/v₀ (s/µM)")
axes[0].set_title("Gráfico de Dixon: las rectas se cruzan en [I] = −Kᵢ", loc="left", fontsize=11); axes[0].legend(fontsize=8)
axes[1].plot(I_exp, km_ap, "o", color=viz.COLORS["datos"], ms=6)
xx = np.linspace(0, max(I_exp), 20); axes[1].plot(xx, aji["km"] * (1 + xx / aji["ki"]), color=viz.INK_MUTED, lw=1.5)
axes[1].set_xlabel("[I] (mM)"); axes[1].set_ylabel("K_M aparente (mM)"); axes[1].set_title("Gráfico secundario: pendiente = K_M/Kᵢ", loc="left", fontsize=11)
fig;
# 4) IC50 -> Ki (Cheng–Prusoff)
S_ensayo = 8.0
ic50 = Ki_real * (1 + S_ensayo / Km_real)
print(f"4) Cheng–Prusoff: con [S] = {S_ensayo} mM la IC₅₀ sería {ic50:.1f} mM y devuelve K_i = {cin.cheng_prusoff(ic50, S_ensayo, Km_real):.2f} mM")
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
inh = ref.get("inhibitors", [])
display(pd.DataFrame([(d.get("name"), d.get("kind"), d.get("versus", ""), d.get("ki", ""), d.get("note", ""), d.get("source", "")) for d in inh],
                     columns=["inhibidor / regulador", "tipo", "frente a", "Kᵢ / efecto", "nota", "fuente"]))
act = ref.get("activators", [])
display(pd.DataFrame([(d.get("name"), d.get("fold", ""), d.get("ec50", ""), d.get("note", ""), d.get("source", "")) for d in act],
                     columns=["activador", "activación (veces)", "EC₅₀", "nota", "fuente"]))
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
