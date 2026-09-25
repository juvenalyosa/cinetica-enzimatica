"""Section 13. Michaelis–Menten from the mechanism."""
from ..celdas import code, md

# ============================================================================ 13. MM
md(r"""
## 13. Michaelis–Menten from the mechanism

> 🎯 **In this section** you will see why an enzyme's rate **saturates** when there is a lot of
> substrate, where the Michaelis–Menten equation comes from (no need to take it on faith: we will
> derive it and simulate it) and how K<sub>M</sub>, V<sub>max</sub>, k<sub>cat</sub> and
> k<sub>cat</sub>/K<sub>M</sub> are obtained from data.

---

### 💡 The analogy: the supermarket checkout lane

A supermarket checkout lane (the **enzyme**) serves customers (the **substrate**). With few
customers, each one who arrives is served right away: the lane's "rate" grows with the number of
customers. With many customers a queue forms: the cashier works flat out and it makes no difference
if more arrive. The rate **saturates**.

[[fig:cajero_saturacion | Three panels of a checkout lane with few, some and many customers, and the three corresponding stretches of the saturation curve]]

That behavior has the shape of a hyperbola and is described by two numbers:

* **V**<sub>max</sub>: the rate with the lane always busy (all the enzyme in the ES form).
* **K**<sub>M</sub>: the customer concentration at which the lane works at **half** its maximum.
  A small K<sub>M</sub> means the enzyme "fills up" with little substrate.

---

### 🧮 The equation in words

**Why do we want an equation?** To predict how much product comes out per second for any amount of
substrate, and to summarize a whole enzyme in two numbers that can be compared across labs,
mutants and drugs.

The idea is simple: the rate is the **maximum rate** multiplied by the **fraction of time the lane
is busy**:

$$\text{rate} \;=\; \underbrace{\text{rate with the lane always busy}}_{V_{\max}} \;\times\; \underbrace{\text{fraction of time busy}}_{\text{from 0 to 1}}$$

[[fig:mm_anatomia | The Michaelis–Menten equation as the product of the maximum rate and the fraction of enzyme occupied]]

The occupied fraction, $[\mathrm{S}]/(K_M + [\mathrm{S}])$, is 0 without substrate, 1/2 when
$[\mathrm{S}] = K_M$, and approaches 1 (without ever reaching it) when there is a huge amount. That's all.

---

### 📐 The full equation, term by term

**The mechanism behind it (Michaelis and Menten 1913; Briggs and Haldane 1925):**

$$\mathrm{E + S \;\underset{k_{-1}}{\overset{k_1}{\rightleftharpoons}}\; ES \;\overset{k_2}{\longrightarrow}\; E + P}$$

[[fig:mecanismo_mm | Scheme of the mechanism: the enzyme binds the substrate (k1), can release it (k−1) or turn it into product (k2)]]

| Constant | What it describes | In the analogy | Units |
|---|---|---|---|
| $k_1$ | the substrate enters the active site (productive collision) | the customer reaches the lane | M⁻¹ s⁻¹ |
| $k_{-1}$ | the substrate is released without reacting | leaves without buying | s⁻¹ |
| $k_2$ (= k<sub>cat</sub>) | the ES complex crosses the barrier and releases the product | pays and leaves with the purchase | s⁻¹ |

**The derivation, step by step.** (Follow each line: it is just algebra.)

1. Rate of product formation: $v_0 = k_2[\mathrm{ES}]$ — *customers paying per second = cashier's
   speed × customers being served*.
2. **Steady-state assumption**: ES forms as fast as it is consumed, so its concentration barely
   changes: $k_1[\mathrm{E}][\mathrm{S}] = (k_{-1} + k_2)[\mathrm{ES}]$ — *as many customers reach the
   lane as leave it (buying or not)*.
3. The total enzyme is shared out: $[\mathrm{E}]_0 = [\mathrm{E}] + [\mathrm{ES}]$, so $[\mathrm{E}] = [\mathrm{E}]_0 - [\mathrm{ES}]$ —
   *each lane is either free or busy*.
4. Substituting (3) into (2) and solving: $[\mathrm{ES}] = \dfrac{[\mathrm{E}]_0[\mathrm{S}]}{K_M + [\mathrm{S}]}$ with
   $K_M \equiv \dfrac{k_{-1}+k_2}{k_1}$.
5. Plugging (4) into (1):

$$\boxed{v_0 \;=\; {\color{#2a78d6}{V_{\max}}}\;\times\;{\color{#1baf7a}{\frac{[\mathrm{S}]}{K_M + [\mathrm{S}]}}}} \qquad\text{with}\qquad {\color{#2a78d6}{V_{\max} = k_2[\mathrm{E}]_0}}$$

The **blue** factor is the maximum rate; the **green** one, the fraction of enzyme occupied (the
same colors as in the drawing).

| Term | What it is | In the analogy |
|---|---|---|
| $v_0$ | initial rate of product formation | customers leaving with their purchase per second |
| ${\color{#2a78d6}{V_{\max}}}$ | rate with all the enzyme occupied: k<sub>cat</sub>·[E]<sub>0</sub> | every lane serving nonstop |
| [S] | substrate concentration | how many customers are in the store |
| $K_M$ | [S] at which the enzyme is half occupied: (k<sub>−1</sub> + k<sub>2</sub>)/k<sub>1</sub> | how many customers it takes for the lane to be busy half the time |

> ⚠️ **Watch out.** $K_M$ equals the dissociation constant $K_d = k_{-1}/k_1$ only when
> $k_2 \ll k_{-1}$ (the customer leaves without buying far more often than they buy). In general
> K<sub>M</sub> is not a pure "affinity".

---

### 🎛️ What if…

| If… | then… | because… |
|---|---|---|
| [S] ≪ K<sub>M</sub> | v₀ ≈ (V<sub>max</sub>/K<sub>M</sub>)·[S]: **a straight line** | the enzyme is almost empty and every substrate molecule counts; slope = (k<sub>cat</sub>/K<sub>M</sub>)·[E]<sub>0</sub> |
| [S] = K<sub>M</sub> | v₀ = V<sub>max</sub>/2 | the lane is busy half the time |
| [S] ≫ K<sub>M</sub> | v₀ → V<sub>max</sub>: **saturation** | there is always a queue: more customers don't speed up the cashier |
| **K<sub>M</sub> goes up** | the curve stretches to the right | the enzyme "grabs" the substrate less well or releases it faster: more is needed for the same rate |
| **V<sub>max</sub> goes up** | the whole curve scales upward | more enzyme, or a faster enzyme |

---

### 🔬 Let's simulate the mechanism

Instead of taking the formula on faith, **let's solve the mechanism numerically** and watch the
hyperbola appear on its own. First, the movie of the concentrations:
""")

code(r'''
# @title 🎛️ The mechanism in action: E + S ⇌ ES → E + P
# 🎛️ Move k1, k-1, k2, [E]0 and [S]0 and watch the pre-steady state and product formation
interactivo.explorar_mecanismo()
''')

code(r'''
# @title 🔁 From the simulation to the hyperbola
k2 = 60  # @param {type:"slider", min:6, max:600, step:6}
# k2 = k_cat (s⁻¹); try 6 and 600 (exercise 3)
E0, k1, k_1 = 0.05, 1.0, 50.0                   # µM, µM⁻¹s⁻¹, s⁻¹
S = np.array([5, 10, 20, 40, 80, 150, 300, 600, 1200, 2400])       # µM
v0_sim = cin.initial_rates_from_simulation(E0, S, k1, k_1, k2, t_window=(0.01, 0.1))
ajuste = cin.fit_michaelis_menten(S, v0_sim)
teoria = cin.steady_state_parameters(k1, k_1, k2)
fig = viz.plot_initial_rates_from_ode(S, v0_sim, fit=ajuste)
v0_mm = cin.michaelis_menten(S, ajuste["vmax"], ajuste["km"])
tabla_sim = pd.DataFrame({"[S]₀": S, "simulated v₀": v0_sim, "v₀ from the hyperbola": v0_mm,
                          "difference": v0_sim - v0_mm, "occupied fraction": S / (ajuste["km"] + S)})
datos_sim = viz.datos(
    tabla_sim, "from the simulation to the hyperbola",
    "Each row is a complete simulation of the mechanism with a different [S]₀: v₀ is the initial slope of the "
    "product curve. The line in the plot is the hyperbola that best passes through these points.",
    x="[S]₀", y="simulated v₀",
    calculadas={"v₀ from the hyperbola": "V_max·[S]₀ / (K_M + [S]₀), with V_max and K_M from the fit",
                "difference": "simulated v₀ − v₀ from the hyperbola",
                "occupied fraction": "[S]₀ / (K_M + [S]₀)"},
    unidades={"[S]₀": "µM", "simulated v₀": "µM/s", "v₀ from the hyperbola": "µM/s", "difference": "µM/s"},
    formatos={"[S]₀": "{:g}", "simulated v₀": "{:.4f}", "v₀ from the hyperbola": "{:.4f}", "difference": "{:+.1e}",
              "occupied fraction": "{:.0%}"},
    resaltar={int(np.argmin(np.abs(S - ajuste["km"]))): ("≈ K_M", "naranja")}, barra="simulated v₀")
viz.mostrar(fig, datos_sim, viz.tarjetas([
    ("K_M from the fit", f"{ajuste['km']:.1f}", "µM", "fitting the hyperbola to the simulated v₀", "azul"),
    ("K_M from the formula", f"{teoria['km']:.1f}", "µM", "(k₋₁ + k₂)/k₁", "agua"),
    ("V_max from the fit", f"{ajuste['vmax']:.3f}", "µM/s", "plateau of the curve", "azul"),
    ("V_max from the formula", f"{k2 * E0:.3f}", "µM/s", "k₂·[E]₀", "agua"),
], titulo="Do the simulation and the paper derivation agree?",
   nota="If the pairs match, the derivation works: the hyperbola comes out of the mechanism on its own."))
''')

code(r'''
# @title 🎛️ The hyperbola: move V_max and K_M
# 🎛️ The hyperbola: move Vmax and Km
interactivo.explorar_michaelis_menten()
''')

md(r"""
### 🔬 The real data: three numbers that come out of the curve

| Constant | How it is obtained | What it measures | In the analogy | Glucokinase (literature) |
|---|---|---|---|---|
| **K<sub>M</sub>** (or S<sub>0.5</sub>) | fit of v₀ versus [S] | apparent affinity; at what [S] the enzyme is half full | customers needed for half occupancy | ≈ 7.7 mM glucose |
| **k<sub>cat</sub>** = V<sub>max</sub>/[E]<sub>0</sub> | V<sub>max</sub> from the fit ÷ enzyme concentration | turnover: reactions per second per enzyme molecule | customers per second that **one** lane serves flat out | ≈ 62–66 s⁻¹ |
| **k<sub>cat</sub>/K<sub>M</sub>** | ratio of the two above | efficiency at low [S]; physical limit ≈ 10⁸–10⁹ M⁻¹s⁻¹ (diffusion) | how well the lane makes use of an almost empty store | ≈ 8 × 10³ M⁻¹s⁻¹ |

Now an "experiment": data simulated with the literature parameters and 4 % noise, treated the way a
biochemist would: a **nonlinear fit** and, for comparison, the three classic linearizations.

> ⚠️ **Watch out.** **Lineweaver–Burk** (1/v₀ versus 1/[S]) amplifies the error of the
> low-concentration points: today it is used to *visualize*, and the nonlinear fit to *quantify*.
""")

code(r'''
# @title 🔬 A simulated experiment and its fit
rng = np.random.default_rng(7)
S_mM = np.array([0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 40, 60])
E0_uM = 0.2                                     # enzyme in the tube (µM)
Km_real = valor("s_half_mm", 7.5)               # we treat the enzyme as hyperbolic in this exercise
Vmax_real = valor("kcat_s", 60.0) * E0_uM       # µM/s
v_obs = cin.michaelis_menten(S_mM, Vmax_real, Km_real) * (1 + 0.04 * rng.standard_normal(S_mM.size))
aj = cin.fit_michaelis_menten(S_mM, v_obs)
fig = viz.plot_michaelis_menten(S_mM, v_obs, fit=aj, title=f"The enzyme reaches half its maximum at {aj['km']:.1f} mM substrate",
                                subtitle="Data simulated with 4 % noise (points) and nonlinear Michaelis–Menten fit (line)")
kc = cin.kcat_km_from_fit(aj, E0_uM)
v_ajuste = cin.michaelis_menten(S_mM, aj["vmax"], aj["km"])
tabla_exp = pd.DataFrame({"[S]": S_mM, "measured v₀": v_obs, "fitted v₀": v_ajuste, "residual": v_obs - v_ajuste,
                          "occupied fraction": S_mM / (aj["km"] + S_mM)})
datos_exp = viz.datos(
    tabla_exp, "the simulated experiment",
    "Each row is a test tube: a glucose concentration and the measured initial rate (with 4 % noise). "
    "The nonlinear fit looks for the V_max and K_M that make the residuals as small as possible.",
    x="[S]", y="measured v₀",
    calculadas={"fitted v₀": "V_max·[S] / (K_M + [S]), with V_max and K_M from the fit",
                "residual": "measured v₀ − fitted v₀ (what the curve does not explain: the noise)",
                "occupied fraction": "[S] / (K_M + [S]) = fitted v₀ / V_max"},
    unidades={"[S]": "mM", "measured v₀": "µM/s", "fitted v₀": "µM/s", "residual": "µM/s"},
    formatos={"[S]": "{:g}", "measured v₀": "{:.2f}", "fitted v₀": "{:.2f}", "residual": "{:+.2f}", "occupied fraction": "{:.0%}"},
    resaltar={int(np.argmin(np.abs(S_mM - aj["km"]))): ("≈ K_M", "naranja")}, barra="measured v₀")
viz.mostrar(fig, datos_exp, viz.tarjetas([
    ("V_max", f"{aj['vmax']:.2f}", f"± {aj['vmax_err']:.2f} µM/s", "the plateau: all the enzyme occupied", "azul"),
    ("K_M", f"{aj['km']:.2f}", f"± {aj['km_err']:.2f} mM", f"[S] at half rate (R² = {aj['r2']:.4f})", "naranja"),
    ("k_cat = V_max/[E]₀", f"{kc['kcat']:.0f}", "s⁻¹", f"one reaction every {cin.turnover_time(kc['kcat'])*1000:.0f} ms per enzyme", "agua"),
    ("k_cat/K_M", f"{kc['kcat_over_km_molar']:.2e}", "M⁻¹s⁻¹", "efficiency when substrate is scarce", "violeta"),
], titulo="The four numbers that come out of the curve"))
''')

code(r'''
# @title 📐 The three classic linearizations
fig = viz.plot_linearizations(S_mM, v_obs, title="Three ways to turn the hyperbola into a straight line",
                              subtitle="The axis intercepts (gray points) give V_max and K_M; orange circles: the points that distort Lineweaver–Burk the most")
lb = cin.linear_fit(*cin.lineweaver_burk(S_mM, v_obs))
tabla_lin = pd.DataFrame({"[S]": S_mM, "v₀": v_obs, "1/[S]": 1 / S_mM, "1/v₀": 1 / v_obs,
                          "v₀/[S]": v_obs / S_mM, "[S]/v₀": S_mM / v_obs})
datos_lin = viz.datos(
    tabla_lin, "the three linearizations",
    "The same 12 tubes, transformed. Each plot uses a pair of columns: Lineweaver–Burk plots 1/v₀ versus "
    "1/[S]; Eadie–Hofstee, v₀ versus v₀/[S]; Hanes–Woolf, [S]/v₀ versus [S]. Look at the first two rows: "
    "with little glucose, 1/[S] and 1/v₀ shoot up, and those points, the noisiest ones, dominate the Lineweaver–Burk line.",
    calculadas={"1/[S]": "1 ÷ [S]  (x axis of Lineweaver–Burk)", "1/v₀": "1 ÷ v₀  (y axis of Lineweaver–Burk)",
                "v₀/[S]": "v₀ ÷ [S]  (x axis of Eadie–Hofstee; its y axis is v₀)",
                "[S]/v₀": "[S] ÷ v₀  (y axis of Hanes–Woolf; its x axis is [S])"},
    unidades={"[S]": "mM", "v₀": "µM/s", "1/[S]": "mM⁻¹", "1/v₀": "s/µM", "v₀/[S]": "µM·s⁻¹/mM", "[S]/v₀": "mM·s/µM"},
    formatos={"[S]": "{:g}", "v₀": "{:.2f}", "1/[S]": "{:.3f}", "1/v₀": "{:.3f}", "v₀/[S]": "{:.3f}", "[S]/v₀": "{:.3f}"},
    resaltar={0: ("weighs a lot in L‑B", "naranja"), 1: ("weighs a lot in L‑B", "naranja")}, barra="1/v₀")
viz.mostrar(fig, datos_lin, viz.tarjetas([
    ("V_max from Lineweaver–Burk", f"{1/lb['intercept']:.2f}", "µM/s", f"nonlinear fit: {aj['vmax']:.2f}", "gris"),
    ("K_M from Lineweaver–Burk", f"{lb['slope']/lb['intercept']:.2f}", "mM", f"nonlinear fit: {aj['km']:.2f}", "gris"),
], titulo="Lineweaver–Burk versus the nonlinear fit",
   nota="The straight lines help visualize the type of behavior; to quantify, use the nonlinear fit."))
''')

md(r"""
> ✅ **Takeaway.** Rate = maximum rate × fraction of enzyme occupied. The hyperbola **emerges**
> from the mechanism E + S ⇌ ES → E + P; K<sub>M</sub> is the [S] of half occupancy,
> k<sub>cat</sub> the reactions per second of each enzyme, and k<sub>cat</sub>/K<sub>M</sub> its
> efficiency when substrate is scarce. They are measured with a nonlinear fit.
""")
