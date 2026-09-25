"""Section 6. Eyring."""
from ..celdas import code, md

# ============================================================================ 6. Eyring
md(r"""
## 6. Energy hills and reaction rate

> 🎯 **In this section** you will see why a lower energy hill means a faster reaction, and how
> the Eyring equation turns the height of that hill into a number of reactions per second.

---

### 💡 The analogy: a ball trying to cross a hill

Picture a ball at the bottom of a valley. Between that valley (the **reactants**: glucose + ATP)
and the neighbouring valley (the **products**: glucose‑6‑phosphate + ADP) there is a hill. The ball
is never still: heat shakes it constantly and throws it uphill again and again. Almost every time
it falls short and rolls back. Once in a long while, an especially strong kick carries it all the
way to the summit, and it drops down the other side: **the reaction has happened**.

[[fig:colina_intentos | A ball in a valley tries to cross a hill; most attempts fail and one reaches the summit]]

In chemistry the hill is called the **activation barrier** and its summit the **transition
state**. An enzyme does not kick the ball harder: **it lowers the hill**.

---

### 🧮 The equation in words

**Why do we want an equation?** To go from "the hill is high" to "the reaction happens 60 times
per second". With a number we can compare the simulation with the experiment.

The idea is the same as in any game of chance: if you roll a die many times, the number of
"sixes" per minute is *(rolls per minute)* × *(probability of rolling a six)*. Here:

$$\text{rate} \;=\; \underbrace{\text{attempts per second}}_{\text{a huge number}} \;\times\; \underbrace{\text{probability of reaching the summit}}_{\text{tiny}}$$

[[fig:eyring_anatomia | The Eyring equation split into two factors: attempts per second and probability of success]]

---

### 📐 The full equation, term by term (Eyring, 1935)

$$k \;=\; {\color{#2a78d6}{\kappa\,\frac{k_B T}{h}}}\;\times\;{\color{#eb6834}{\exp\!\left(-\frac{\Delta G^{\ddagger}}{RT}\right)}}$$

The **blue** factor is the attempts; the **orange** one, the probability that an attempt reaches
the summit. Each symbol:

| Term | What it is | In the analogy | Typical value |
|---|---|---|---|
| $k$ | rate constant | how many balls cross per second | ≈ 60 s⁻¹ for glucokinase |
| $k_B T/h$ | attempt frequency | how many times per second the ball is thrown uphill | 6.2 × 10¹² s⁻¹ at 25 °C |
| $\Delta G^{\ddagger}$ | free energy of activation | the height of the hill | 10–25 kcal/mol |
| $R T$ | available thermal energy | the typical strength of the thermal kicks | 0.59 kcal/mol at 25 °C |
| $\kappa$ | transmission coefficient | fraction of balls that, once at the summit, do not turn back | ≈ 1 |

Why an **exponential**? Because getting high up requires many favourable kicks in a row, and the
probability of them all happening multiplies: every extra step of height divides the probability
by the same factor. That is why what matters is the ratio $\Delta G^{\ddagger}/RT$: how many
"typical kicks" tall the hill is.

---

### 🎛️ What if…

| If… | then… | because… |
|---|---|---|
| **ΔG‡ goes up by 1.36 kcal/mol** | the reaction is **10 times slower** | the exponential drops to one tenth (the golden rule) |
| **T goes from 25 to 37 °C** | for a 15 kcal/mol barrier, the rate **almost doubles** | the thermal kicks (RT) are 4 % stronger and the exponential amplifies it |
| **the enzyme lowers ΔG‡ by 10 kcal/mol** | the reaction runs **10⁷ times** faster | 10 / 1.36 ≈ 7 steps of ×10 |

Try it yourself with the sliders:
""")

code(r'''
# @title 🎛️ Explore the Eyring equation
# 🎛️ Explore the Eyring equation: barrier and temperature
interactivo.explorar_eyring()
''')

code(r'''
# @title 🎛️ Draw your own energy hill
# 🎛️ Energy profile: barrier (uphill) and reaction energy (drop between valleys)
interactivo.explorar_perfil_energia()
''')

code(r'''
# @title 📉 Every 1.36 kcal/mol, a factor of 10
barreras = np.linspace(5, 30, 200)
kcat = valor("kcat_s", 60.0)
fig, ax = viz.figure(10.5, 5.8)
k_curva = cin.eyring_rate(barreras)
ax.fill_between(barreras, k_curva.min() / 1e3, k_curva, color=viz._tint(viz.COLORS["reactivo"], 0.9), zorder=1, lw=0)
ax.semilogy(barreras, k_curva, color=viz.COLORS["reactivo"], lw=2.8, zorder=3)
puntos = [(cin.barrier_from_rate(kcat), f"glucokinase\nk_cat ≈ {kcat:.0f} s⁻¹", viz.COLORS["ts"], (12, 6), "left"),
          (cin.barrier_from_rate(1e-8), "a reaction that\ntakes years", viz.INK_SECONDARY, (-16, -10), "right")]
for dg, texto, color, off, ha in puntos:
    viz._guide(ax, "v", dg, start=k_curva.min() / 1e3, end=cin.eyring_rate(dg))
    viz._keypoint(ax, dg, cin.eyring_rate(dg), color, size=11)
    ax.annotate(f"{texto}\nΔG‡ ≈ {dg:.1f} kcal/mol", (dg, cin.eyring_rate(dg)), xytext=off, textcoords="offset points",
                ha=ha, va="bottom" if off[1] > 0 else "top", fontsize=11, color=viz.INK)
# the staircase: 3 steps of 1.36 kcal/mol up from the glucokinase barrier
dg0 = cin.barrier_from_rate(kcat)
for j in range(3):
    a, b = dg0 + j * 1.364, dg0 + (j + 1) * 1.364
    ax.plot([a, b, b], [cin.eyring_rate(a), cin.eyring_rate(a), cin.eyring_rate(b)], color=viz.INK_MUTED, lw=1.2, zorder=4)
    ax.annotate("÷10", (b, np.sqrt(cin.eyring_rate(a) * cin.eyring_rate(b))), xytext=(4, 0), textcoords="offset points",
                ha="left", va="center", fontsize=10, color=viz.INK_SECONDARY)
ax.set_xlim(5, 30); ax.set_ylim(k_curva.min() / 1e3, k_curva.max() * 10)
viz._finish(ax, "ΔG‡: height of the hill (kcal/mol)", "k (s⁻¹, log scale)",
            "Every 1.36 kcal/mol of barrier changes the rate 10-fold",
            "Each grey step raises the hill by 1.36 kcal/mol and divides the rate by 10.")
dg_tab = np.r_[np.arange(10, 30, 2.0), dg0]
dg_tab.sort()
k_tab = cin.eyring_rate(dg_tab)
viz.mostrar(fig, viz.datos(
    pd.DataFrame({"ΔG‡": dg_tab, "k": k_tab, "time per reaction": [viz.duracion(1 / k) for k in k_tab]}),
    "the Eyring equation",
    "Each row is a point on the line: pick a hill height, compute k with the Eyring equation at 25 °C "
    "and, turning it around, how long one reaction takes on average.",
    x="ΔG‡", y="k", calculadas={"k": "(k_B·T/h) · exp(−ΔG‡/RT), with k_B·T/h = 6.2 × 10¹² s⁻¹ and RT = 0.59 kcal/mol",
                               "time per reaction": "1 / k"},
    unidades={"ΔG‡": "kcal/mol", "k": "s⁻¹"}, formatos={"ΔG‡": "{:.1f}"},
    resaltar={int(np.argmin(np.abs(dg_tab - dg0))): ("k_cat", "naranja")}))
''')

md(r"""
> ✅ **Takeaway.** Rate = attempts × probability of success. The attempts are almost fixed
> (~10¹³ per second); what decides everything is the height of the hill, and every 1.36 kcal/mol
> is a factor of 10. An enzyme speeds up the reaction by **lowering the hill**, not by pushing harder.
""")
