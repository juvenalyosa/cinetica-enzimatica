"""Section 14. Cooperativity (Hill)."""
from ..celdas import code, md

# ============================================================================ 14. Hill
md(r"""
## 14. Cooperativity: glucokinase is a sensor

> 🎯 **In this section** you will see why the glucokinase curve is not a hyperbola but a
> **sigmoid**, what the Hill coefficient *n* measures, and how an enzyme with **a single site** can
> behave like a switch that decides when the pancreas releases insulin.

---

### 💡 The analogy: dimmer versus switch

A **dimmer** (hyperbola) raises the light little by little from the very start. A **threshold
switch** (sigmoid) does almost nothing up to a certain point and then turns on all at once. The
pancreas needs a switch: below ~5 mM glucose it must not release insulin, and above it, it must.
Glucokinase is that switch.

[[fig:regulador_interruptor | A hyperbola (dimmer) compared with a sigmoid (switch): the sigmoid goes from 10 to 90 percent with a much smaller change in substrate]]

Look at the shaded bands: to go from 10 % to 90 % activity, the hyperbola needs glucose to be
multiplied **by 81**; the glucokinase sigmoid, only **by 13**.

---

### 🧮 The equation in words

**Why do we want another equation?** Because Michaelis–Menten can only draw hyperbolas, and the
glucokinase data are not a hyperbola. We need a number that says **how abrupt** the switch is.

The idea is the same as in section 13 (rate = maximum × fraction occupied), but the substrate
concentration enters **raised to a power n**:

$$\text{rate} \;=\; V_{\max} \times \frac{(\text{substrate})^{\,n}}{(\text{threshold})^{\,n} + (\text{substrate})^{\,n}}$$

Raising to *n* > 1 makes concentrations **below** the threshold count even less and those **above**
it count even more: the curve flattens at the start and steepens in the middle. That is the switch.

---

### 📐 The full equation, term by term (Hill)

$$v_0 = \frac{V_{\max}[\mathrm{S}]^{\color{#eb6834}{n}}}{{\color{#2a78d6}{S_{0.5}}}^{\,\color{#eb6834}{n}} + [\mathrm{S}]^{\color{#eb6834}{n}}}$$

| Term | What it is | In the analogy | Glucokinase |
|---|---|---|---|
| $V_{\max}$ | maximum rate | the light at full brightness | k<sub>cat</sub>·[E]<sub>0</sub> |
| ${\color{#2a78d6}{S_{0.5}}}$ | concentration at which half of V<sub>max</sub> is reached (the analog of K<sub>M</sub>) | where the switch's threshold is | ≈ 7.7 mM glucose |
| ${\color{#eb6834}{n}}$ | **Hill coefficient**: how abrupt the change is | how "suddenly" the light turns on | ≈ 1.7 |

$n = 1$ is Michaelis–Menten; $n > 1$ is positive cooperativity.

---

### 🎛️ What if…

| If… | to go from 10 % to 90 % activity, [S] must be multiplied by… | example |
|---|---|---|
| n = 1 | **81** | a Michaelis–Menten enzyme |
| n = 1.7 | **~13** | glucokinase |
| n = 4 | **3** | hemoglobin |

---

### 💡 The surprise: an enzyme with memory

Glucokinase is a **monomer** with a single glucose site: there can be no "communication between
subunits" as in hemoglobin. Its cooperativity is **kinetic**: the enzyme switches between a poorly
active (open) form and an active (closed) form at a rate **comparable to that of catalysis**
(k<sub>ex</sub> ≈ 5–100 s⁻¹ versus k<sub>cat</sub> ≈ 60 s⁻¹, Larion et al. 2012).

[[fig:cooperatividad_cinetica | Glucokinase alternates between a slow open form and a fast closed form; with plenty of glucose it has no time to relax and stays in the fast form]]

With plenty of glucose the enzyme has no time to relax to the slow form between one cycle and the
next, so it works more than it "should". This is the **mnemonic** or **slow-transition** model
(Storer and Cornish‑Bowden 1977; Neet and Ainslie; Cárdenas 1984).

> ✍️ **Think about it.** A competitive inhibitor such as N‑acetylglucosamine *suppresses*
> cooperativity (n → 1), because it keeps the site occupied and breaks the "memory". The MD
> simulations of section 5 (domains that open and close) are the molecular picture of that change.

---

### 🔬 The real data

Move *n* and S<sub>0.5</sub>, and compare with the hyperbola:
""")

code(r'''
# @title 🎛️ Sigmoid versus hyperbola
# 🎛️ Move n and S0.5: compare with the hyperbola and look at the 10 %–90 % window
interactivo.explorar_hill()
''')

code(r'''
# @title 🩸 Glucokinase in the blood glucose range
S_g = np.linspace(0.01, 30, 300)
s_half, n_h = valor("s_half_mm", 7.5), valor("hill_n", 1.7)
fig = viz.plot_hill_vs_mm(S_g, cin.michaelis_menten(S_g, 1.0, s_half), cin.hill(S_g, 1.0, s_half, n_h), n_hill=n_h, s_half=s_half,
                          ylabel="v₀ / V_max", title="In the blood range, glucokinase responds like a switch",
                          subtitle="Between 4 and 10 mM glucose (the physiological range) the sigmoid is much more sensitive than the hyperbola")
fig.set_size_inches(10.5, 5.8)
ax = fig.axes[0]
viz._kband(ax, 4, 7, "fasting blood\nglucose (4–7 mM)", color=viz._tint(viz.COLORS["ts"], 0.9), y_text=0.97)
s10, s90 = cin.substrate_at_fraction(1.0, s_half, n_h, 0.1), cin.substrate_at_fraction(1.0, s_half, n_h, 0.9)
S_tab = np.array(sorted({1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, round(s_half, 2)}), dtype=float)
mm_tab, hill_tab = cin.michaelis_menten(S_tab, 1.0, s_half), cin.hill(S_tab, 1.0, s_half, n_h)
tabla_hill = pd.DataFrame({"[glucose]": S_tab, "[S]/S₀.₅": S_tab / s_half, "hyperbola (n = 1)": mm_tab,
                           f"sigmoid (n = {n_h:g})": hill_tab, "difference": hill_tab - mm_tab})
resalte = {int(np.argmin(np.abs(S_tab - s_half))): ("S₀.₅", "naranja")}
resalte.update({i: ("fasting", "agua") for i, sv in enumerate(S_tab) if 4 <= sv <= 7 and i not in resalte})
datos_hill = viz.datos(
    tabla_hill, "sigmoid versus hyperbola",
    "The two curves evaluated at round glucose concentrations, as a fraction of V_max. With little glucose the "
    "sigmoid runs below (the enzyme barely responds); near S₀.₅ they cross; above it, the sigmoid rises faster.",
    x="[glucose]", y=["hyperbola (n = 1)", f"sigmoid (n = {n_h:g})"],
    calculadas={"[S]/S₀.₅": "[glucose] ÷ S₀.₅",
                "hyperbola (n = 1)": "[S] / (S₀.₅ + [S])  (Michaelis–Menten with K_M = S₀.₅)",
                f"sigmoid (n = {n_h:g})": "[S]ⁿ / (S₀.₅ⁿ + [S]ⁿ)  (Hill)",
                "difference": "sigmoid − hyperbola"},
    unidades={"[glucose]": "mM", "hyperbola (n = 1)": "v₀/V_max", f"sigmoid (n = {n_h:g})": "v₀/V_max", "difference": "v₀/V_max"},
    formatos={"[glucose]": "{:g}", "[S]/S₀.₅": "{:.2f}", "hyperbola (n = 1)": "{:.0%}", f"sigmoid (n = {n_h:g})": "{:.0%}",
              "difference": "{:+.0%}"},
    resaltar=resalte, barra=f"sigmoid (n = {n_h:g})")
viz.mostrar(fig, datos_hill, viz.tarjetas([
    ("S₀.₅", f"{s_half:g}", "mM", "glucose at half activity", "naranja"),
    ("Hill coefficient", f"{n_h:g}", "", "n > 1: sigmoid", "naranja"),
    ("from 10 % to 90 %", f"{s10:.1f} → {s90:.1f}", "mM", f"[S] × {s90/s10:.0f}", "naranja"),
    ("with n = 1", "× 81", "", "the hyperbola would need much more", "azul"),
], titulo="The glucokinase switch in numbers"))
''')

code(r'''
# @title 📐 Hill fit and Hill plot
# Hill fit to simulated data and Hill plot (the slope is n)
S_h = np.array([1, 2, 3, 4, 5, 6, 7.5, 9, 11, 14, 18, 25, 35, 50])
v_h = cin.hill(S_h, 10.0, s_half, n_h) * (1 + 0.03 * rng.standard_normal(S_h.size))
ajh = cin.fit_hill(S_h, v_h)
fig = viz.plot_hill_plot(S_h, v_h, ajh["vmax"], title=f"The data rise more steeply than the reference: n = {ajh['n']:.2f}",
                         subtitle="log[v/(V_max − v)] versus log[S]: the slope is the Hill coefficient; the dashed line is n = 1")
frac = v_h / (ajh["vmax"] - v_h)
tabla_hp = pd.DataFrame({"[S]": S_h, "v₀": v_h, "log₁₀[S]": np.log10(S_h), "v₀/(V_max − v₀)": frac,
                         "log₁₀[v₀/(V_max − v₀)]": np.log10(frac)})
datos_hp = viz.datos(
    tabla_hp, "the Hill plot",
    "The 14 tubes of the simulated experiment. For the Hill plot both columns are transformed: the x axis is "
    "log₁₀[S] and the y axis is log₁₀[v₀/(V_max − v₀)], with V_max from the fit. If the enzyme follows the Hill equation, "
    "the points fall on a straight line of slope n; they cross 0 on the y axis right at [S] = S₀.₅ (half occupied).",
    x="log₁₀[S]", y="log₁₀[v₀/(V_max − v₀)]",
    calculadas={"log₁₀[S]": "base-10 logarithm of [S]",
                "v₀/(V_max − v₀)": "enzyme \"on\" ÷ enzyme \"off\"",
                "log₁₀[v₀/(V_max − v₀)]": "logarithm of the previous column (= n·log[S] − n·log S₀.₅)"},
    unidades={"[S]": "mM", "v₀": "µM/s"},
    formatos={"[S]": "{:g}", "v₀": "{:.2f}", "log₁₀[S]": "{:.3f}", "v₀/(V_max − v₀)": "{:.3f}", "log₁₀[v₀/(V_max − v₀)]": "{:+.3f}"},
    resaltar={int(np.argmin(np.abs(np.log10(frac)))): ("≈ S₀.₅", "naranja")})
viz.mostrar(fig, datos_hp, viz.tarjetas([
    ("V_max", f"{ajh['vmax']:.2f}", "", "plateau of the fit", "azul"),
    ("S₀.₅", f"{ajh['s_half']:.2f}", "mM", "glucose at half activity", "naranja"),
    ("n (Hill)", f"{ajh['n']:.2f}", f"± {ajh['n_err']:.2f}", "n > 1: positive cooperativity", "naranja"),
], titulo="Hill fit to simulated data with 3 % noise"))
''')

md(r"""
### 🔬 Mutations that move the switch: GCK‑MODY and hyperinsulinism

A mutation that raises S<sub>0.5</sub> (or lowers k<sub>cat</sub>) makes the pancreas "see" less
glucose than there is and release insulin late: high blood glucose from birth, but stable
(GCK‑MODY). An activating mutation lowers S<sub>0.5</sub>: too much insulin, congenital
hypoglycemia. In the analogy: the switch shifts to the right or to the left. With the published
values (Valentínová 2012; Sayed 2009) we plot what each mutant does at 5 mM glucose:
""")

code(r'''
# @title 🧬 GCK‑MODY mutants: the switch shifts
mut = ref.get("mutants", [])
if mut:
    S_g = np.linspace(0.01, 30, 300)
    elegidos = [m for m in mut if m["name"] in ("V244G", "G223S", "I110N", "W99L", "M197I")]
    curvas = [("wild type", s_half, n_h, 1.0)] + [(m["name"], m["s_half_mm"], m.get("hill_n", n_h), m.get("kcat_rel", 1.0)) for m in elegidos]
    fig, ax = viz.figure(11.5, 6.4)
    viz._kband(ax, 4, 7, color=viz._tint(viz.COLORS["ts"], 0.9))
    etiquetas = []
    for (nombre, s05, nn, krel), color in zip(curvas, [viz.INK] + list(viz.PALETTE[: len(elegidos)])):
        y = krel * cin.hill(S_g, 1.0, s05, nn)
        ax.plot(S_g, y, color=color, lw=3.2 if nombre == "wild type" else 2.4, zorder=4 if nombre == "wild type" else 3,
                label=f"{nombre}: S₀.₅ = {s05:g} mM, k_cat ×{krel:.2f}")
        y5 = krel * cin.hill(5.0, 1.0, s05, nn)
        ax.plot([5.0], [y5], "o", ms=8, mfc=color, mec="white", mew=1.6, zorder=6)
        etiquetas.append([y[-1], nombre])
    # direct labels at the end of each curve, spaced so they don't overlap
    etiquetas.sort()
    y_max = max(e[0] for e in etiquetas)
    for j in range(1, len(etiquetas)):
        etiquetas[j][0] = max(etiquetas[j][0], etiquetas[j - 1][0] + 0.06 * y_max)
    for y_lab, nombre in etiquetas:
        ax.annotate(nombre, (S_g[-1], y_lab), xytext=(6, 0), textcoords="offset points", ha="left", va="center",
                    fontsize=10.5, color=viz.INK, fontweight="semibold" if nombre == "wild type" else "normal",
                    annotation_clip=False)
    viz._guide(ax, "v", 5.0, color=viz.INK_SECONDARY)
    ax.annotate("5 mM (band: fasting glucose, 4–7 mM)\neach point: the activity with which\nthe pancreas \"sees\" glucose", (5.0, 0.0), xytext=(96, 14),
                textcoords="offset points", ha="left", va="bottom", fontsize=10.5, color=viz.INK_SECONDARY)
    ax.set_xlim(0, 30); ax.set_ylim(bottom=0)
    viz._finish(ax, "[glucose] (mM)", "activity relative to wild-type V_max",
                "Each mutation shifts the insulin switch",
                "Hill curves with the published values: further right or lower = the pancreas \"sees\" less glucose (GCK‑MODY)")
    ax.legend(fontsize=10, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.13), frameon=False)
    tabla = pd.DataFrame([(m["name"], m["kind"], m["s_half_mm"], m.get("hill_n"), m.get("kcat_s"),
                           round(m.get("kcat_rel", 1.0) * cin.hill(5.0, 1.0, m["s_half_mm"], m.get("hill_n", n_h)), 3)) for m in mut],
                         columns=["mutant", "phenotype", "S₀.₅ (mM)", "n", "k_cat (s⁻¹)", "activity at 5 mM (rel.)"])
    tabla.loc[len(tabla)] = ["wild type", "-", s_half, n_h, valor("kcat_s"), round(cin.hill(5.0, 1.0, s_half, n_h), 3)]
    viz.mostrar(fig, viz.tabla(tabla, titulo="The mutants, one by one", nota="Source: " + ref.get("mutants_source", "")))
else:
    viz.mostrar(viz.mensaje("No mutant data in the reference table.", tipo="ojo"))
''')

md(r"""
> ✅ **Takeaway.** The Hill coefficient *n* measures how abrupt the switch is: with n ≈ 1.7
> glucokinase goes from off to on with a change in glucose ~6 times smaller than a
> Michaelis–Menten enzyme. Its cooperativity does not come from several subunits but from a
> **kinetic memory**, and mutations that shift S<sub>0.5</sub> shift the insulin threshold.
""")
