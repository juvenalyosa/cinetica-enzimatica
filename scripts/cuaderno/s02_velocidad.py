"""Sección 2. velocidad."""
from .celdas import code, md

# ============================================================================ 2. velocidad
md(r"""
## 2. ¿Qué es la velocidad de una reacción?

> 🎯 **En esta sección** aprenderás qué significa "velocidad" en una reacción química, por qué los
> bioquímicos miden la **velocidad inicial v₀** y cómo se mide en el laboratorio para la
> glucoquinasa.

---

### 💡 La analogía: contar monedas que caen en un frasco

La velocidad de una reacción es **cuánto producto aparece por segundo** (o cuánto sustrato
desaparece). Es como contar las monedas que caen en un frasco: si a los 5 s hay 5 y a los 10 s
hay 10, caen a razón de una por segundo. Si en un tubo con enzima la glucosa‑6‑fosfato pasa de
0 a 10 µM en 10 s, la velocidad es 1 µM/s.

[[fig:velocidad_inicial | Tres frascos que se llenan de producto a ritmo constante, y una curva de progreso cuya pendiente inicial es v0]]

---

### 🧮 La ecuación en palabras

**¿Para qué queremos una ecuación?** Para que "rápido" y "lento" se conviertan en un número con
unidades, que podamos comparar entre experimentos y con la teoría.

$$\text{velocidad} \;=\; \frac{\text{cuánto producto nuevo apareció}}{\text{cuánto tiempo pasó}}$$

### 📐 La ecuación completa, término a término

$$v_0 \;=\; \left(\frac{{\color{#e87ba4}{\Delta[\mathrm{P}]}}}{{\color{#52514e}{\Delta t}}}\right)_{t\,\to\,0}$$

| Término | Qué es | En la analogía | Unidades |
|---|---|---|---|
| $v_0$ | velocidad **inicial** | el ritmo de monedas al principio | µM/s |
| $\Delta[\mathrm{P}]$ | aumento de la concentración de producto | monedas nuevas en el frasco | µM |
| $\Delta t$ | tiempo transcurrido | segundos en el cronómetro | s |
| $t \to 0$ | "medido al comienzo" | contar antes de que se acaben las monedas | — |

**¿Por qué "velocidad inicial"?** Al principio del experimento el sustrato apenas se ha gastado
y el producto apenas se ha acumulado, así que la velocidad es constante y refleja solo la
enzima y la concentración de sustrato que pusimos. Después la curva se dobla (el sustrato se
agota). Por eso los bioquímicos miden la **velocidad inicial, v₀**: la pendiente de la curva
de producto al comienzo.

> ⚠️ **Ojo.** Si mides la pendiente demasiado tarde, cuando la curva ya se dobló, obtendrás una
> velocidad más baja que la real y todos los parámetros que calcules después saldrán mal.

---

### 🔬 Cómo se mide en la glucoquinasa: el ensayo acoplado

No se ve la glucosa‑6‑fosfato directamente: se acopla otra enzima (glucosa‑6‑fosfato
deshidrogenasa, G6PDH) que la convierte y, al hacerlo, produce NADPH, que absorbe luz a 340 nm.
El espectrofotómetro registra la absorbancia frente al tiempo, y la pendiente inicial es v₀ (el
ensayo "acoplado a G6PDH" de casi todos los artículos de la tabla anterior).

[[fig:ensayo_acoplado | La glucoquinasa produce G6P; la G6PDH la convierte y produce NADPH; el espectrofotómetro mide la absorbancia a 340 nm, cuya pendiente es proporcional a v0]]

---

### 🎛️ Qué pasa si… cambio la cantidad de sustrato

Simulemos esa curva para tres concentraciones iniciales de sustrato, [S]₀, y midamos la pendiente
inicial de cada una:
""")

code(r'''
# Curva de progreso simulada: producto frente a tiempo para tres concentraciones de sustrato
E0, k1, k_1, k2 = 0.05, 1.0, 50.0, 60.0          # µM, µM⁻¹s⁻¹, s⁻¹, s⁻¹
fig, ax = viz.figure(7.5, 4.2)
for S0, color in zip((50, 200, 1000), viz.sequential_blue(3)):
    sim = cin.simulate_mechanism(e0=E0, s0=S0, k1=k1, k_minus1=k_1, k2=k2, t_end=40.0)
    ax.plot(sim["t"], sim["P"], color=color, lw=2, label=f"[S]₀ = {S0} µM")
    ventana = (sim["t"] > 0.02) & (sim["t"] < 1.0)
    pend = cin.linear_fit(sim["t"][ventana], sim["P"][ventana])["slope"]
    ax.plot([0, 8], [0, 8 * pend], color=viz.INK_MUTED, ls="--", lw=1)
    ax.annotate(f"v₀ = {pend:.2f} µM/s", (8, 8 * pend), xytext=(6, 0), textcoords="offset points", fontsize=9, color=viz.INK_SECONDARY)
ax.set_xlabel("Tiempo (s)"); ax.set_ylabel("Producto (µM)"); ax.legend()
ax.set_title("Curva de progreso: la pendiente inicial es la velocidad v₀", loc="left")
fig;
print("Con más sustrato la pendiente inicial es mayor, pero no proporcionalmente. Ese 'no proporcionalmente' es toda la cinética enzimática.")
''')

md(r"""
> ✍️ **Para pensar.** Si duplicas [S]₀ de 200 a 400 µM, ¿se duplica v₀? Anota tu predicción; en
> la sección 13 veremos por qué la respuesta es "no del todo".

> ✅ **Para llevar.** La velocidad es la pendiente de la curva de producto frente al tiempo, y se
> mide **al principio** (v₀). En la glucoquinasa se mide con un ensayo acoplado: cada G6P
> formada produce un NADPH que absorbe luz a 340 nm.
""")
