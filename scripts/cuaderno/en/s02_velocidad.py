"""Section 2. Reaction rate."""
from ..celdas import code, md

# ============================================================================ 2. rate
md(r"""
## 2. What is the rate of a reaction?

> 🎯 **In this section** you will learn what "rate" means for a chemical reaction, why
> biochemists measure the **initial rate v₀** and how it is measured in the laboratory for
> glucokinase.

---

### 💡 The analogy: counting coins falling into a jar

The rate of a reaction is **how much product appears per second** (or how much substrate
disappears). It is like counting the coins that fall into a jar: if there are 5 after 5 s and
10 after 10 s, they fall at one per second. If, in a tube with enzyme, glucose‑6‑phosphate goes
from 0 to 10 µM in 10 s, the rate is 1 µM/s.

[[fig:velocidad_inicial | Three jars that fill with product at a constant pace, and a progress curve whose initial slope is v0]]

---

### 🧮 The equation in words

**Why do we want an equation?** So that "fast" and "slow" become a number with units, which we
can compare between experiments and with theory.

$$\text{rate} \;=\; \frac{\text{how much new product appeared}}{\text{how much time went by}}$$

### 📐 The full equation, term by term

$$v_0 \;=\; \left(\frac{{\color{#e87ba4}{\Delta[\mathrm{P}]}}}{{\color{#52514e}{\Delta t}}}\right)_{t\,\to\,0}$$

| Term | What it is | In the analogy | Units |
|---|---|---|---|
| $v_0$ | **initial** rate | the pace of coins at the start | µM/s |
| $\Delta[\mathrm{P}]$ | increase in product concentration | new coins in the jar | µM |
| $\Delta t$ | elapsed time | seconds on the stopwatch | s |
| $t \to 0$ | "measured at the start" | count before the coins run out | — |

**Why "initial rate"?** At the start of the experiment the substrate has barely been used up and
the product has barely accumulated, so the rate is constant and reflects only the enzyme and the
substrate concentration we added. Later the curve bends (the substrate runs out). That is why
biochemists measure the **initial rate, v₀**: the slope of the product curve at the beginning.

> ⚠️ **Watch out.** If you measure the slope too late, when the curve has already bent, you will
> get a rate lower than the real one, and every parameter you compute afterwards will be wrong.

---

### 🔬 How it is measured for glucokinase: the coupled assay

Glucose‑6‑phosphate is not seen directly: a second enzyme is coupled (glucose‑6‑phosphate
dehydrogenase, G6PDH) that converts it and, in doing so, produces NADPH, which absorbs light at
340 nm. The spectrophotometer records the absorbance over time, and the initial slope is v₀ (the
"G6PDH-coupled" assay of almost every paper in the previous table).

[[fig:ensayo_acoplado | Glucokinase produces G6P; G6PDH converts it and produces NADPH; the spectrophotometer measures the absorbance at 340 nm, whose slope is proportional to v0]]

---

### 🎛️ What if… I change the amount of substrate

Let's simulate that curve for three initial substrate concentrations, [S]₀, and measure the
initial slope of each one:
""")

code(r'''
# @title 📈 Progress curves and initial rate
# Simulated progress curve: product versus time for three substrate concentrations
E0, k1, k_1, k2 = 0.05, 1.0, 50.0, 60.0          # µM, µM⁻¹s⁻¹, s⁻¹, s⁻¹
curvas = []
for S0 in (50, 200, 1000):
    sim = cin.simulate_mechanism(e0=E0, s0=S0, k1=k1, k_minus1=k_1, k2=k2, t_end=150.0, n_points=4000)
    ventana = (sim["t"] > 0.02) & (sim["t"] < 1.0)
    pend = cin.linear_fit(sim["t"][ventana], sim["P"][ventana])["slope"]
    curvas.append((sim["t"], sim["P"], pend, f"[S]₀ = {S0} µM"))
fig = viz.plot_progress_curves(curvas, title="The initial slope is the rate v₀",
                               subtitle="With 50, 200 and 1000 µM substrate. Dashed: the tangents at the start. The curve bends when the substrate runs out")
v = [c[2] for c in curvas]
t_tab = np.array([0, 0.5, 1, 2, 5, 10, 20, 40, 80, 150])              # time points for the data table (s)
col = {c[3]: f"[P] with {c[3].split('= ')[1]}" for c in curvas}
viz.mostrar(fig, viz.tarjetas(
    [("v₀ with 50 µM", f"{v[0]:.2f}", "µM/s", None, "azul"),
     ("v₀ with 200 µM", f"{v[1]:.2f}", "µM/s", f"4 × more substrate → × {v[1] / v[0]:.1f} in v₀", "azul"),
     ("v₀ with 1000 µM", f"{v[2]:.2f}", "µM/s", f"20 × more substrate → × {v[2] / v[0]:.1f} in v₀", "azul")],
    nota="With more substrate the initial slope is larger, but not proportionally. "
         "That \"not proportionally\" is what enzyme kinetics is all about."),
    viz.datos(pd.DataFrame({"time": t_tab, **{col[c[3]]: np.interp(t_tab, c[0], c[1]) for c in curvas}}),
              "progress curves",
              "Each row is an instant of the simulated experiment; each column, one of the three curves. At first "
              "the product grows at a constant pace (the slope is v₀); later it slows down because the substrate is used up.",
              x="time", y=list(col.values()),
              unidades={"time": "s", **{c: "µM" for c in col.values()}},
              formatos={"time": "{:g}", **{c: "{:.1f}" for c in col.values()}},
              resaltar={2: ("v₀", "azul")},
              nota=f"Row v₀: in the first second, [P] ÷ time ≈ v₀ ({v[0]:.2f}, {v[1]:.2f} and {v[2]:.2f} µM/s)."))
''')

md(r"""
> ✍️ **Think about it.** If you double [S]₀ from 200 to 400 µM, does v₀ double? Write down your
> prediction; in section 13 we will see why the answer is "not quite".

> ✅ **Takeaway.** The rate is the slope of the product-versus-time curve, and it is measured
> **at the start** (v₀). For glucokinase it is measured with a coupled assay: every G6P formed
> produces one NADPH that absorbs light at 340 nm.
""")
