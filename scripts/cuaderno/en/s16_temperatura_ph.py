"""Section 16. Temperature and pH."""
from ..celdas import code, md

# ============================================================================ 16. T and pH
md(r"""
## 16. Temperature and pH

> 🎯 **In this section** you will see why heating speeds up reactions (and by how much), how to
> split the barrier into **enthalpy** and **entropy** with the Arrhenius and Eyring equations, and
> why an enzyme's activity versus pH is **bell-shaped**.

---

### 💡 The analogy: heating = stronger thermal kicks

Go back to the ball of section 6. Temperature is the strength of the thermal kicks: when you heat,
all the molecules move more and the **fraction** that has enough energy to cross the hill **grows**.
Because that fraction is the "tail" of a distribution, a small change in temperature enlarges it
far more than you might expect.

[[fig:temperatura_sacudidas | Distribution of molecular energies at 25 and 37 °C; the tail that exceeds the barrier is larger at 37 °C]]

(Until the enzyme denatures: too much heat unravels the protein.)

---

### 🧮 The equation in words

**Why do we want an equation?** To read the **height of the barrier** from a simple experiment:
measure k at several temperatures. How much the rate rises on heating tells you how tall the hill
is.

$$\text{ln(rate)} \;=\; \text{constant} \;-\; \frac{\text{height of the barrier}}{\text{energy of the thermal kicks}}$$

If you plot ln k versus 1/T you get a **straight line**, and its slope is the barrier: the steeper
it is, the taller the hill and the more sensitive the reaction is to temperature.

---

### 📐 The full equations, term by term

**Arrhenius**: $\ln k = \ln A - {\color{#eb6834}{E_a}}/RT$. In a plot of ln k versus 1/T the slope
is $-{\color{#eb6834}{E_a}}/R$: the **activation energy** is "how much the rate rises on heating".

**Eyring**: $\ln(k/T) = \ln(k_B/h) + {\color{#4a3aa7}{\Delta S^{\ddagger}}}/R - {\color{#eb6834}{\Delta H^{\ddagger}}}/RT$.
It splits the barrier into **enthalpy** (energy that has to be supplied) and **entropy** (order that
has to be imposed: a negative ΔS‡ means the transition state is more ordered than the reactants).

| Term | What it is | In the analogy |
|---|---|---|
| ${\color{#eb6834}{E_a,\ \Delta H^{\ddagger}}}$ | activation energy / activation enthalpy | the height of the hill |
| ${\color{#4a3aa7}{\Delta S^{\ddagger}}}$ | activation entropy | how narrow the pass is: if the summit must be "threaded" precisely, ΔS‡ < 0 |
| $A$ | pre-exponential factor | the attempts per second |
| $RT$ | thermal energy | the typical strength of the thermal kicks |

---

### 🎛️ What if…

| If… | going from 25 to 37 °C multiplies k by… | because… |
|---|---|---|
| E<sub>a</sub> = 12 kcal/mol | **~2.2** | the tail above a tall hill grows a lot with the kicks |
| E<sub>a</sub> = 6 kcal/mol | **~1.5** | over a low hill many already get across: heating adds less |

The higher the barrier, the more sensitive to temperature.
""")

code(r'''
# @title 🎛️ Enthalpy, entropy and temperature
# 🎛️ Activation enthalpy and entropy: look at the Eyring plot and k(T)
interactivo.explorar_temperatura()
''')

code(r'''
# @title 🌡️ Arrhenius and Eyring with data
T_K = np.array([283.15, 288.15, 293.15, 298.15, 303.15, 308.15, 313.15])
dH, dS = 12.0, -8.0                       # kcal/mol, cal/mol/K (illustrative)
k_T = np.array([cin.eyring_rate(dH - T * dS / 1000, T) for T in T_K]) * (1 + 0.03 * rng.standard_normal(T_K.size))
aj_arr, aj_eyr = cin.fit_arrhenius(T_K, k_T), cin.fit_eyring(T_K, k_T)
fig, axes = viz.figure(14, 5.8, ncols=2)
viz.plot_arrhenius(T_K, k_T, fit=aj_arr, ax=axes[0], title="Arrhenius: ln k versus 1000/T")
viz.plot_eyring(T_K, k_T, fit=aj_eyr, ax=axes[1], title="Eyring: ln(k/T) versus 1000/T")
viz._fig_title(fig, "More heat, more speed: the slope measures how much",
               "Seven temperatures from 10 to 40 °C (top axis). A steep line = a reaction very sensitive to temperature.")
tabla_T = pd.DataFrame({"T": T_K, "T (°C)": T_K - 273.15, "measured k": k_T, "1000/T": 1000 / T_K,
                        "ln k": np.log(k_T), "ln(k/T)": np.log(k_T / T_K)})
datos_T = viz.datos(
    tabla_T, "Arrhenius and Eyring",
    "Seven measurements of the rate constant, one per temperature (3 % noise). Both lines use the same "
    "x axis, 1000/T; Arrhenius puts ln k on the y axis and Eyring, ln(k/T). Their slopes give E_a and ΔH‡.",
    x="1000/T", y=["ln k", "ln(k/T)"],
    calculadas={"T (°C)": "T − 273.15", "1000/T": "1000 ÷ T (in K⁻¹ × 1000, so the numbers are convenient)",
                "ln k": "natural logarithm of k  (y axis of Arrhenius; slope = −E_a/R)",
                "ln(k/T)": "natural logarithm of k ÷ T  (y axis of Eyring; slope = −ΔH‡/R)"},
    unidades={"T": "K", "T (°C)": "°C", "measured k": "s⁻¹", "1000/T": "K⁻¹"},
    formatos={"T": "{:.2f}", "T (°C)": "{:.0f}", "measured k": "{:.1f}", "1000/T": "{:.4f}", "ln k": "{:.3f}", "ln(k/T)": "{:.3f}"},
    resaltar={int(np.argmin(np.abs(T_K - 298.15))): ("25 °C", "azul")},
    barra="measured k")
viz.mostrar(fig, datos_T, viz.tarjetas([
    ("E_a (Arrhenius)", f"{aj_arr['ea_kcal']:.1f}", "kcal/mol", "how much k rises on heating", "naranja"),
    ("ΔH‡ (Eyring)", f"{aj_eyr['delta_h_kcal']:.1f}", "kcal/mol", "energy that has to be supplied", "naranja"),
    ("ΔS‡ (Eyring)", f"{aj_eyr['delta_s_cal']:.1f}", "cal/mol/K", "negative: the TS is more ordered", "violeta"),
    ("ΔG‡ (25 °C)", f"{aj_eyr['delta_g_kcal']:.1f}", "kcal/mol", "ΔH‡ − TΔS‡", "azul"),
], titulo="What the two lines say", nota="How they are related: E_a ≈ ΔH‡ + RT (0.6 kcal/mol at 25 °C)."))
''')

md(r"""
---

### 💡 The pH analogy: two groups, two conditions

Every enzyme has an optimal pH: above or below it, the groups that do the chemistry lose the
charge they need. In glucokinase, **Asp205** (the catalytic base) must be **deprotonated** to accept
the proton from glucose, and **Lys169** must be **protonated** to stabilize the phosphate. It is like
a door with two locks: it only opens when both are in the right position.

[[fig:ph_campana | Bell-shaped activity versus pH curve, with the protonation states of Asp205 and Lys169 at low, optimal and high pH]]

### 🧮 The equation in words

$$\text{activity} \;=\; \frac{\text{maximum activity}}{1 + (\text{penalty if the base is missing}) + (\text{penalty if the acid is missing})}$$

Each penalty is almost zero at the right pH and grows tenfold for each pH unit you move away: hence
the bell.

### 📐 The full equation, term by term

$$v = \frac{v_{\max}}{1 + {\color{#eb6834}{10^{\,pK_1 - \mathrm{pH}}}} + {\color{#4a3aa7}{10^{\,\mathrm{pH} - pK_2}}}}$$

| Term | What it is | When it matters |
|---|---|---|
| ${\color{#eb6834}{10^{\,pK_1 - \mathrm{pH}}}}$ | ratio of **protonated** Asp205 (useless as a base) / deprotonated | at low pH: it "switches off" the enzyme |
| ${\color{#4a3aa7}{10^{\,\mathrm{pH} - pK_2}}}$ | ratio of **deprotonated** Lys169 (no positive charge) / protonated | at high pH: it "switches off" the enzyme |
| (pK<sub>1</sub> + pK<sub>2</sub>)/2 | position of the maximum | the optimal pH |

> 🔬 **The real data.** For human glucokinase the measured optimum is pH 8.5–8.7 (Šimčíková and
> Heneberg 2019); for decades it was thought to be lower because ATP acidifies the assay buffers:
> an example of how a technical detail can change a textbook "fact".
""")

code(r'''
# @title 🎛️ The pH bell curve
# 🎛️ pKa of the base and the acid: the pH bell curve
interactivo.explorar_ph()
''')

md(r"""
> ✅ **Takeaway.** Heating enlarges the tail of molecules that clear the barrier: the taller the
> hill, the more the reaction gains with temperature (Arrhenius measures E<sub>a</sub>; Eyring splits
> it into ΔH‡ and ΔS‡). The pH decides whether the catalytic groups carry the right charge: two
> pK<sub>a</sub> values draw a bell with the maximum between them.
""")
