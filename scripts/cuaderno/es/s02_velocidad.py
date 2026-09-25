"""Sección 2. velocidad."""
from ..celdas import code, md

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
# @title 📈 Curvas de progreso y velocidad inicial
# Curva de progreso simulada: producto frente a tiempo para tres concentraciones de sustrato
E0, k1, k_1, k2 = 0.05, 1.0, 50.0, 60.0          # µM, µM⁻¹s⁻¹, s⁻¹, s⁻¹
curvas = []
for S0 in (50, 200, 1000):
    sim = cin.simulate_mechanism(e0=E0, s0=S0, k1=k1, k_minus1=k_1, k2=k2, t_end=150.0, n_points=4000)
    ventana = (sim["t"] > 0.02) & (sim["t"] < 1.0)
    pend = cin.linear_fit(sim["t"][ventana], sim["P"][ventana])["slope"]
    curvas.append((sim["t"], sim["P"], pend, f"[S]₀ = {S0} µM"))
fig = viz.plot_progress_curves(curvas, title="La pendiente inicial es la velocidad v₀",
                               subtitle="Con 50, 200 y 1000 µM de sustrato. Discontinuas: las tangentes al inicio. La curva se dobla cuando el sustrato se agota")
v = [c[2] for c in curvas]
t_tab = np.array([0, 0.5, 1, 2, 5, 10, 20, 40, 80, 150])              # instantes de la tabla de datos (s)
col = {c[3]: f"[P] con {c[3].split('= ')[1]}" for c in curvas}
viz.mostrar(fig, viz.tarjetas(
    [("v₀ con 50 µM", f"{v[0]:.2f}", "µM/s", None, "azul"),
     ("v₀ con 200 µM", f"{v[1]:.2f}", "µM/s", f"4 × más sustrato → × {v[1] / v[0]:.1f} en v₀", "azul"),
     ("v₀ con 1000 µM", f"{v[2]:.2f}", "µM/s", f"20 × más sustrato → × {v[2] / v[0]:.1f} en v₀", "azul")],
    nota="Con más sustrato la pendiente inicial es mayor, pero no proporcionalmente. "
         "Ese «no proporcionalmente» es toda la cinética enzimática."),
    viz.datos(pd.DataFrame({"tiempo": t_tab, **{col[c[3]]: np.interp(t_tab, c[0], c[1]) for c in curvas}}),
              "curvas de progreso",
              "Cada fila es un instante del experimento simulado; cada columna, una de las tres curvas. Al principio "
              "el producto crece a ritmo constante (la pendiente es v₀); después se frena porque se gasta el sustrato.",
              x="tiempo", y=list(col.values()),
              unidades={"tiempo": "s", **{c: "µM" for c in col.values()}},
              formatos={"tiempo": "{:g}", **{c: "{:.1f}" for c in col.values()}},
              resaltar={2: ("v₀", "azul")},
              nota=f"Fila v₀: en el primer segundo, [P] ÷ tiempo ≈ v₀ ({v[0]:.2f}, {v[1]:.2f} y {v[2]:.2f} µM/s)."))
''')

md(r"""
> ✍️ **Para pensar.** Si duplicas [S]₀ de 200 a 400 µM, ¿se duplica v₀? Anota tu predicción; en
> la sección 13 veremos por qué la respuesta es "no del todo".

> ✅ **Para llevar.** La velocidad es la pendiente de la curva de producto frente al tiempo, y se
> mide **al principio** (v₀). En la glucoquinasa se mide con un ensayo acoplado: cada G6P
> formada produce un NADPH que absorbe luz a 340 nm.
""")
