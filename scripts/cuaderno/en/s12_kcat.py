"""Section 12. barrier -> kcat."""
from ..celdas import code, md

# ============================================================================ 12. barrier -> kcat

md(r"""
## 12. From the barrier to *k*<sub>cat</sub>

> 🎯 **In this section** you will turn the calculated barriers into a turnover number
> *k*<sub>cat</sub> and compare it, honestly, with experiment and with another published
> calculation.

---

### 💡 The analogy: a ruler with two scales

A thermometer can have two scales, °C and °F: the same measurement read in two ways. The Eyring
equation is a ruler like that: on top you read the **height of the hill**, underneath the
**reactions per second**. From the height of the hill we estimate how many times per second the
enzyme completes the reaction: that number is **k**<sub>cat</sub>, the "turnover number". We compare
it with the value measured in the lab and with another published QM/MM calculation.

[[fig:escalera_kcat | Ruler with the barrier on top and k_cat underneath; the calculated and experimental barriers are marked]]

---

### 🧮 The equation in words

It is the same equation as in section 6, *attempts × probability of success*, used in both
directions:

* **from simulation to experiment**: given the calculated barrier, how many reactions per second does it predict?
* **from experiment to barrier**: given the measured k<sub>cat</sub>, what hill height does it imply?

---

### 📐 The full equation, term by term

$$k_{\mathrm{cat}} \;=\; {\color{#2a78d6}{\frac{k_BT}{h}}}\;{\color{#eb6834}{e^{-\Delta G^{\ddagger}/RT}}}
\qquad\Longleftrightarrow\qquad
\Delta G^{\ddagger} \;=\; RT\,\ln\!\frac{k_BT}{h\,k_{\mathrm{cat}}}$$

| Term | What it is | In the analogy | Value |
|---|---|---|---|
| $k_{\mathrm{cat}}$ | reactions per second per enzyme molecule | the bottom scale | measured: ≈ 62–66 s⁻¹ |
| ${\color{#2a78d6}{k_BT/h}}$ | attempts per second | — | 6.2 × 10¹² s⁻¹ at 25 °C |
| ${\color{#eb6834}{\exp(-\Delta G^{\ddagger}/RT)}}$ | probability of crossing the summit | — | tiny |
| $\Delta G^{\ddagger}$ | height of the hill | the top scale | experiment: ≈ 15 kcal/mol |

---

### 🎛️ What if… (or why we must be honest)

Our ΔG‡ comes from a semiempirical method (PM7), a frozen environment and a harmonic
approximation. An error of 1.4 kcal/mol is already a factor of 10 in rate. That is why we compare
orders of magnitude and, above all, *trends*, not decimals.

| If the calculated barrier is off by… | k<sub>cat</sub> is off by a factor of… |
|---|---|
| 0.5 kcal/mol | ~2 |
| 1.4 kcal/mol | ~10 |
| 4 kcal/mol | ~850 |

> ⚠️ **Watch out.** A calculated barrier "almost equal" to the experimental one can give a
> k<sub>cat</sub> ten or a hundred times different. Always compare on the energy scale **and** on
> the rate scale.

---

### 🔬 The real data
""")

code(r'''
# @title 🏁 Our barrier versus experiment
T = 298.15
kcat_exp = valor("kcat_s", 60.0)
dG_qmmm = termo.get(T_ref, {}).get("dG_kcal", ts["barrera_kcal"])
filas = [("ΔE‡ PM7, one conformation (this notebook)", ts["barrera_kcal"], viz.COLORS["reactivo"]),
         ("ΔG‡ PM7 + harmonic thermochemistry (this notebook)", dG_qmmm, viz.COLORS["ts"])]
try:
    if len(ok) > 1:
        filas.append(("ΔE‡ averaged over MD snapshots (this notebook)", media, viz.COLORS["reactivo"]))
except NameError:
    pass
filas.append(("published QM/MM ΔE‡, Zhang et al. 2009", valor("qmmm_barrier_literature_kcal", np.nan), viz.PALETTE[6]))
filas.append(("active site in water, COSMO (this notebook)", agua["barrera_kcal"], viz.INK_MUTED))
dG_exp = cin.barrier_from_rate(kcat_exp, T)
etiquetas = {"ΔE‡ PM7, one conformation (this notebook)": "ΔE‡ in the enzyme (PM7)",
             "ΔG‡ PM7 + harmonic thermochemistry (this notebook)": "ΔG‡ with vibrations",
             "ΔE‡ averaged over MD snapshots (this notebook)": "average ΔE‡ (MD)",
             "published QM/MM ΔE‡, Zhang et al. 2009": "Zhang et al. 2009 (QM/MM)",
             "active site in water, COSMO (this notebook)": "active site in water"}
fig = viz.plot_estimates([(etiquetas.get(n, n), b, c) for n, b, c in filas], reference=dG_exp,
                         reference_label=f"experiment (k_cat ≈ {kcat_exp:.0f} s⁻¹)", rate_fn=lambda b: cin.eyring_rate(b, T),
                         title="The QM/MM calculations sit 3–4 kcal/mol above experiment",
                         subtitle="Each point is a barrier estimate; on the right, the rate implied by the Eyring equation")
filas.append(("experiment: ΔG‡ implied by k_cat", dG_exp, None))
tabla = pd.DataFrame([(n, b, viz._sci(cin.eyring_rate(b, T))) for n, b, _ in filas], columns=["estimate", "barrier (kcal/mol)", "k (s⁻¹) from Eyring"])
viz.mostrar(fig, viz.tabla(tabla, formatos={"barrier (kcal/mol)": "{:.1f}"}, titulo="The same estimates, in numbers",
                           nota=f"Experimental k_cat ≈ {kcat_exp:.0f} s⁻¹ ({ref['kcat_s']['source']})."),
            viz.mensaje("Our barrier and that of Zhang et al. (18.3 kcal/mol, a different program and a different QM/MM partition) "
                        "agree within 1 kcal/mol; both overestimate the experimental one (≈15) by 3-4 kcal/mol, that is, a factor "
                        "of ~10²-10³ in k. That is what to expect from a semiempirical method with a fixed environment.", "idea"))
''')

md(r"""
> ✅ **Takeaway.** Eyring turns barriers into rates in both directions. Our barrier and that of
> Zhang et al. agree with each other within 1 kcal/mol and sit 3–4 kcal/mol above the one implied
> by experiment (≈ 15): the simulation gets the mechanism and the order of magnitude of the hill
> right, not the decimals of *k*<sub>cat</sub>.
""")
