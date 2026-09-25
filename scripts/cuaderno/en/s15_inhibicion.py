"""Section 15. Inhibitors and activators."""
from ..celdas import code, md

# ============================================================================ 15. inhibition
md(r"""
## 15. Inhibitors and activators

> 🎯 **In this section** you will tell the types of inhibitors apart by the fingerprint they leave
> on V<sub>max</sub>, K<sub>M</sub> and the Lineweaver–Burk plot, learn three ways to measure
> K<sub>i</sub>, and meet the real inhibitors, regulators and activators of glucokinase.

---

### 💡 The analogy: three ways to slow down the checkout lane

Back to the cashier of section 13. There are three ways to slow them down:

* **Competitive**: someone joins the queue *without buying anything* and takes a turn. With many
  real customers, the intruder goes unnoticed: V<sub>max</sub> does not change, but more substrate
  is needed to reach it (the apparent K<sub>M</sub> goes up).
* **Uncompetitive**: someone jams the till *only when there is already a customer*: V<sub>max</sub>
  goes down and, curiously, so does K<sub>M</sub> (it traps the ES complex).
* **Noncompetitive**: someone turns off the lane's light, whether or not there is a customer: the
  lane works more slowly (V<sub>max</sub> goes down) but customers still come in the same way
  (K<sub>M</sub> does not change).

[[fig:inhibidores_cajero | Three types of inhibitor shown at the supermarket checkout lane, with their effect on Vmax, on KM and their fingerprint on the Lineweaver–Burk plot]]

---

### 🧮 The equation in words

**Why do we want an equation?** To know *how much* an inhibitor slows the enzyme at a given
concentration, and to measure its potency with a single number, K<sub>i</sub>, that does not depend
on the experiment.

For the competitive case the idea is: the intruder "steals" part of the lanes, so **more
substrate** is needed to reach the same occupancy. How much more? A factor that grows with the
amount of intruder:

$$K_M^{\text{apparent}} \;=\; K_M \times \underbrace{\left(1 + \frac{\text{inhibitor}}{\text{inhibitor potency}}\right)}_{\text{how many times more substrate you need}}$$

If there is as much inhibitor as its K<sub>i</sub>, you need twice the substrate; if there is 10 times
its K<sub>i</sub>, 11 times more.

---

### 📐 The full equation, term by term

**The derivation for the competitive case (the others are analogous).** The inhibitor I binds to
the free enzyme: $\mathrm{E + I \rightleftharpoons EI}$ with $K_i = [\mathrm{E}][\mathrm{I}]/[\mathrm{EI}]$.
Now the total enzyme is $[\mathrm{E}]_0 = [\mathrm{E}] + [\mathrm{ES}] + [\mathrm{EI}]$ (*each lane is
free, serving a customer, or blocked by the intruder*) and, repeating the algebra of section 13,

$$v_0 = \frac{V_{\max}[\mathrm{S}]}{K_M{\color{#e34948}{\left(1 + \dfrac{[\mathrm{I}]}{K_i}\right)}} + [\mathrm{S}]}
\qquad\Rightarrow\qquad K_M^{\mathrm{app}} = K_M{\color{#e34948}{\left(1+\frac{[\mathrm{I}]}{K_i}\right)}},\; V_{\max}^{\mathrm{app}} = V_{\max}.$$

The **red** factor (the inhibitor's color in the drawing) is the whole effect of the intruder.

| Term | What it is | In the analogy |
|---|---|---|
| [I] | inhibitor concentration | how many intruders are in the store |
| $K_i$ | dissociation constant of the EI complex | the [I] that blocks half of the free lanes: the **smaller**, the more potent |
| 1 + [I]/K<sub>i</sub> (red) | inhibition factor | how many times more customers you need for the same occupancy |

| Type | Binds to | In the analogy | Apparent V<sub>max</sub> | Apparent K<sub>M</sub> | Lineweaver–Burk fingerprint |
|---|---|---|---|---|---|
| competitive | E | takes a turn | unchanged | × (1 + [I]/K<sub>i</sub>) | lines crossing on the y axis |
| uncompetitive | ES | jams the till with a customer | ÷ (1 + [I]/K<sub>i</sub>′) | ÷ (1 + [I]/K<sub>i</sub>′) | parallel lines |
| noncompetitive | E and ES equally | turns off the light | ÷ (1 + [I]/K<sub>i</sub>) | unchanged | lines crossing on the x axis |
| mixed | E and ES, differently | a bit of each | ÷ (1 + [I]/K<sub>i</sub>′) | × (1+[I]/K<sub>i</sub>)/(1+[I]/K<sub>i</sub>′) | cross to the left of the y axis |

---

### 🎛️ What if…

| If… | then… | because… |
|---|---|---|
| **[I] = K<sub>i</sub>** (competitive) | the apparent K<sub>M</sub> **doubles** | 1 + 1 = 2 |
| **[I] = 10 K<sub>i</sub>** (competitive) | the apparent K<sub>M</sub> is multiplied by **11** | 1 + 10 = 11 |
| **you add a lot of substrate** | a competitive inhibitor **is overcome**; a noncompetitive one, **never** | with a long queue the intruder never gets a turn; a light that is off slows the lane all the same |

Choose the type of inhibitor and move [I] and K<sub>i</sub>:
""")

code(r'''
# @title 🎛️ Inhibitor types and their fingerprint
# 🎛️ Choose the inhibitor type and move [I] and Ki: look at the curve and its Lineweaver–Burk fingerprint
interactivo.explorar_inhibicion()
''')

md(r"""
### 🔬 From data to K<sub>i</sub>: three routes that must agree

1. **Global fit**: all the points (several [S] and several [I]) against the full equation
   → V<sub>max</sub>, K<sub>M</sub> and K<sub>i</sub> with their errors. This is the recommended method today.
2. **Secondary plot**: each [I] is fitted separately and K<sub>M</sub><sup>app</sup> is plotted
   versus [I]: it is a straight line with slope K<sub>M</sub>/K<sub>i</sub>.
3. **Dixon plot**: 1/v₀ versus [I] for several [S]; for a competitive inhibitor the lines cross at
   [I] = −K<sub>i</sub>.

**From IC₅₀ to K<sub>i</sub>.** Pharmacology often measures the **IC₅₀**: the concentration that
cuts the activity in half. In words: *the more substrate competes, the more inhibitor it takes to
cut the rate in half*, so the IC₅₀ depends on [S]. The **Cheng–Prusoff** equation discounts that
effect and returns K<sub>i</sub> (for a competitive inhibitor):

$$K_i = \frac{\mathrm{IC}_{50}}{1 + {\color{#eb6834}{[\mathrm{S}]/K_M}}}$$

The **orange** term (the substrate's color) is the "head start" the substrate has over the
inhibitor.
""")

code(r'''
# @title 🔬 Three routes to measure K_i
Ki_real = 3.0                                   # mM (hypothetical competitive inhibitor)
S_exp = np.array([1, 2, 4, 8, 16, 32, 64]); I_exp = (0, 2, 5, 10)
filas = []
for I in I_exp:
    v = cin.competitive_inhibition(S_exp, Vmax_real, Km_real, I, Ki_real) * (1 + 0.03 * rng.standard_normal(S_exp.size))
    filas += [(s, vv, I) for s, vv in zip(S_exp, v)]
df_i = pd.DataFrame(filas, columns=["S", "v", "I"])
# 1) global fit
aji = cin.fit_inhibition(df_i.S.values, df_i.v.values, df_i.I.values, kind="competitive")
# 2) secondary plot: apparent Km versus [I]
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
axes[0].annotate(f"they cross at [I] = −Kᵢ\nKᵢ ≈ {kd['ki']:.2f} mM", (-kd["ki"], f_lo["intercept"] - f_lo["slope"] * kd["ki"]),
                 xytext=(4, -24), textcoords="offset points", ha="left", va="top", fontsize=11, color=viz.INK, fontweight="semibold")
viz._finish(axes[0], "[I] (mM)", "1/v₀ (s/µM)")
axes[0].set_ylim(bottom=-0.18 * axes[0].get_ylim()[1])
viz._legend(axes[0], loc="upper left", fontsize=10)
axes[0].set_title("Dixon plot", loc="left", fontsize=13, fontweight="semibold", pad=10)
xx = np.linspace(0, max(I_exp), 20)
axes[1].fill_between(xx, aji["km"], aji["km"] * (1 + xx / aji["ki"]), color=viz._tint(viz.PALETTE[7], 0.9), lw=0, zorder=1)
axes[1].plot(xx, aji["km"] * (1 + xx / aji["ki"]), color=viz.INK_SECONDARY, lw=2.4, zorder=2)
viz._markers(axes[1], I_exp, km_ap, viz.COLORS["datos"], size=9)
viz._guide(axes[1], "h", aji["km"])
axes[1].annotate("K_M without inhibitor", (max(I_exp), aji["km"]), xytext=(-4, -6), textcoords="offset points", ha="right", va="top",
                 fontsize=10.5, color=viz.INK_SECONDARY)
axes[1].annotate(f"slope = K_M/Kᵢ\n→ Kᵢ = {sec['ki']:.2f} mM", (0.55 * max(I_exp), aji["km"] * (1 + 0.55 * max(I_exp) / aji["ki"])),
                 xytext=(-14, 10), textcoords="offset points", ha="right", va="bottom", fontsize=11, color=viz.INK, fontweight="semibold")
axes[1].set_ylim(0, None)
viz._finish(axes[1], "[I] (mM)", "apparent K_M (mM)")
axes[1].set_title("Secondary plot", loc="left", fontsize=13, fontweight="semibold", pad=10)
viz._fig_title(fig, f"Three routes, the same Kᵢ ≈ {aji['ki']:.1f} mM (value used in the simulation: {Ki_real:g} mM)",
               "Left: 1/v₀ versus [I] for several [S]. Right: apparent K_M versus [I]. Both from the same data.")
tabla_dixon = df_i.rename(columns={"S": "[S]", "I": "[I]", "v": "v₀"})[["[I]", "[S]", "v₀"]].copy()
tabla_dixon["1/v₀"] = 1 / tabla_dixon["v₀"]
tabla_dixon["v₀ without inhibitor (same [S])"] = [df_i[(df_i.I == 0) & (df_i.S == sv)].v.iloc[0] for sv in tabla_dixon["[S]"]]
tabla_dixon["remaining activity"] = tabla_dixon["v₀"] / tabla_dixon["v₀ without inhibitor (same [S])"]
datos_dixon = viz.datos(
    tabla_dixon, "the Dixon plot",
    "The 28 tubes of the experiment: 4 inhibitor concentrations × 7 glucose concentrations (3 % noise). In the "
    "Dixon plot each [S] is a line: its points are the rows with that [S], with [I] on the x axis and 1/v₀ on the y axis. "
    "With plenty of glucose the competitive inhibitor is barely noticeable (remaining activity close to 100 %).",
    x="[I]", y="1/v₀",
    calculadas={"1/v₀": "1 ÷ v₀", "remaining activity": "v₀ ÷ v₀ without inhibitor, at the same [S]"},
    unidades={"[I]": "mM", "[S]": "mM", "v₀": "µM/s", "1/v₀": "s/µM", "v₀ without inhibitor (same [S])": "µM/s"},
    formatos={"[I]": "{:g}", "[S]": "{:g}", "v₀": "{:.2f}", "1/v₀": "{:.3f}", "v₀ without inhibitor (same [S])": "{:.2f}",
              "remaining activity": "{:.0%}"},
    resaltar={i: ("[I] = 10 mM", "rojo") for i in tabla_dixon.index[(tabla_dixon["[I]"] == max(I_exp)) & (tabla_dixon["[S]"] <= 2)]},
    barra="1/v₀", max_filas=40, abierta=False)
tabla_sec = pd.DataFrame({"[I]": np.array(I_exp, dtype=float), "apparent K_M": km_ap, "K_M app / K_M": km_ap / aji["km"],
                          "1 + [I]/Kᵢ (theory)": 1 + np.array(I_exp) / aji["ki"]})
datos_sec = viz.datos(
    tabla_sec, "the secondary plot",
    "One row per inhibitor concentration: each group of 7 tubes gets its own hyperbola fit, and its apparent K_M "
    "is one point of the secondary plot. If the inhibitor is competitive, K_M app/K_M grows as 1 + [I]/Kᵢ.",
    x="[I]", y="apparent K_M",
    calculadas={"apparent K_M": "Michaelis–Menten fit to the tubes with that [I] only",
                "K_M app / K_M": "apparent K_M ÷ K_M without inhibitor (from the global fit)",
                "1 + [I]/Kᵢ (theory)": "what the competitive equation predicts with Kᵢ from the global fit"},
    unidades={"[I]": "mM", "apparent K_M": "mM"},
    formatos={"[I]": "{:g}", "apparent K_M": "{:.2f}", "K_M app / K_M": "{:.2f}", "1 + [I]/Kᵢ (theory)": "{:.2f}"},
    barra="apparent K_M")
viz.mostrar(fig, datos_dixon, datos_sec, viz.tarjetas([
    ("1 · global fit", f"{aji['ki']:.2f}", f"± {aji['ki_err']:.2f} mM", f"K_M = {aji['km']:.2f} mM, V_max = {aji['vmax']:.2f} µM/s", "azul"),
    ("2 · secondary plot", f"{sec['ki']:.2f}", "mM", f"R² = {sec['r2']:.3f}", "agua"),
    ("3 · Dixon", f"{kd['ki']:.2f}", "mM", "from where the lines cross", "naranja"),
    ("4 · Cheng–Prusoff", f"{cin.cheng_prusoff(ic50, S_ensayo, Km_real):.2f}", "mM", f"from IC₅₀ = {ic50:.1f} mM at [S] = {S_ensayo:g} mM", "violeta"),
], titulo="Kᵢ of the same inhibitor, measured four ways",
   nota="If the methods agree, the (competitive) model describes the data well."))
''')

md(r"""
### 🔬 The real inhibitors and regulators of glucokinase

Table with sources. Three surprising things:

* Glucose‑6‑phosphate, which slows down hexokinases I–III, does **not** inhibit glucokinase: that
  is why the liver keeps phosphorylating glucose even as the product accumulates (Viñuela, Salas
  and Sols, 1963).
* **Mannoheptulose**, which many textbooks call "competitive", turns out to be **mixed** when
  measured carefully (Scruel 1998), and its K<sub>i</sub> varies 100-fold between secondary
  sources: that is why we do not give a number.
* In the liver, the **regulatory protein GKRP** acts as a competitive inhibitor with respect to
  glucose (it raises S<sub>0.5</sub> without changing V; IC₅₀ ≈ 13 mM glucose) and sequesters the
  enzyme in the nucleus when glucose falls; fructose‑6‑phosphate reinforces that sequestration and
  fructose‑1‑phosphate undoes it.
""")

code(r'''
# @title 📋 Real inhibitors and activators of glucokinase
inh = ref.get("inhibitors", [])
tipos = {"competitive": "competitive", "mixed": "mixed", "none": "none", "slow-binding": "slow binding",
         "uncompetitive": "uncompetitive", "noncompetitive": "noncompetitive"}
act = ref.get("activators", [])
viz.mostrar(
    viz.tabla(pd.DataFrame([(d.get("name"), tipos.get(d.get("kind"), d.get("kind")), d.get("versus", ""), d.get("ki", ""), d.get("note", ""), d.get("source", "")) for d in inh],
                           columns=["inhibitor / regulator", "type", "versus", "Kᵢ / effect", "note", "source"]),
              titulo="Inhibitors and regulators of glucokinase"),
    viz.tabla(pd.DataFrame([(d.get("name"), d.get("fold", ""), d.get("ec50", ""), d.get("note", ""), d.get("source", "")) for d in act],
                           columns=["activator", "activation (fold)", "EC₅₀", "note", "source"]),
              titulo="Allosteric activators (GKA)"),
)
''')

md(r"""
### 🔬 Allosteric activators (GKA)

In the analogy: someone **oils the till** and opens the lane sooner. They bind to a site ~20 Å from
the glucose site (Kamata 2004), lower S<sub>0.5</sub> and raise V<sub>max</sub>: the opposite of an
inhibitor. The reference compound RO‑28‑1675 activates the wild-type enzyme ~16-fold with
EC₅₀ ≈ 7 µM (Sayed 2009) and was tested as a drug for type 2 diabetes (Grimsby 2003). What happens
to the activity at 5 mM glucose when you lower S<sub>0.5</sub>? Try it:
""")

code(r'''
# @title 🎛️ An allosteric activator
# 🎛️ An allosteric activator: lowers S0.5 and raises Vmax. What happens to the activity at 5 mM glucose?
interactivo.explorar_activador()
''')

md(r"""
> ✅ **Takeaway.** Each type of inhibitor leaves a different fingerprint: the competitive one raises
> K<sub>M</sub> (and is overcome by substrate), the uncompetitive one lowers V<sub>max</sub> and
> K<sub>M</sub>, the noncompetitive one lowers only V<sub>max</sub>. K<sub>i</sub> is measured by a
> global fit, a secondary plot or a Dixon plot, and Cheng–Prusoff turns an IC₅₀ into a
> K<sub>i</sub>. Activators do the opposite: they lower S<sub>0.5</sub>.
""")
