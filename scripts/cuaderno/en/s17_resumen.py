"""Section 17. summary."""
from ..celdas import code, md

# ============================================================================ 17. summary
md(r"""
## 17. Summary, glossary and exercises

> 🎯 **In this section** we put all the pieces together: the course map, the equations on a single
> page, a glossary and exercises to check what you have learned.

---

### 🗺️ The course map

[[fig:mapa_curso | Concept map: structure, molecular dynamics, QM/MM, from ΔE‡ to ΔG‡, Eyring, Michaelis–Menten, Hill and inhibitors, with temperature and pH as conditions]]

**What we have seen, one sentence each:**

1. An enzyme speeds up a reaction by **lowering the barrier** of the transition state, like a
   guide who knows the lowest mountain pass.
2. The crystal structure (3FGU) shows the reactants already aligned; MD shows that the enzyme
   **keeps them aligned**.
3. QM/MM lets us **see the chemistry** inside the enzyme; care is needed with the partition, the
   link atoms, the electrostatics and the QM–MM repulsion.
4. The transition state is a **saddle point**: it is located (scan → NEB → dimer) and **validated**
   with a single imaginary frequency and with the path that goes down to R and P.
5. From ΔE‡ to ΔG‡: vibrations, averaging over conformations and sensitivity to the method. A
   computational result carries **uncertainty**, and our barrier agrees with another published
   calculation and overestimates experiment by ~3–4 kcal/mol.
6. Eyring turns the barrier into *k*<sub>cat</sub>; 1.4 kcal/mol is a factor of 10.
7. Michaelis–Menten **emerges** from the mechanism E + S ⇌ ES → E + P; K<sub>M</sub>, V<sub>max</sub>,
   k<sub>cat</sub> and k<sub>cat</sub>/K<sub>M</sub> are obtained from data by a nonlinear fit.
8. Glucokinase is **sigmoidal** (n ≈ 1.7) without having several subunits: kinetic
   cooperativity; its mutants shift the insulin switch.
9. Each type of inhibitor leaves a different fingerprint; K<sub>i</sub> is obtained by a global fit,
   a secondary plot or a Dixon plot; activators do the opposite.
10. Temperature and pH modulate k through ΔH‡, ΔS‡ and the pK<sub>a</sub> values of the active site.

---

### 🧮 The course equations, on one page

**Initial rate** (section 2): the rate is the initial slope of the product curve.

$$v_0 = \left(\frac{\Delta[\mathrm{P}]}{\Delta t}\right)_{t\to 0}$$

**Eyring** (sections 6 and 12): rate = attempts per second × probability of reaching the summit.

$$k = \kappa\,\frac{k_B T}{h}\,\exp\!\left(-\frac{\Delta G^{\ddagger}}{RT}\right)$$

**Michaelis–Menten** (section 13): the checkout lane saturates; with plenty of substrate the enzyme works flat out.

$$v_0 = \frac{V_{\max}[\mathrm{S}]}{K_M + [\mathrm{S}]}$$

**Hill** (section 14): a switch; little activity up to a threshold, then a lot.

$$v_0 = \frac{V_{\max}[\mathrm{S}]^{n}}{S_{0.5}^{\,n} + [\mathrm{S}]^{n}}$$

**Competitive inhibition** (section 15): the intruder demands more substrate for the same rate.

$$K_M^{\mathrm{app}} = K_M\left(1 + \frac{[\mathrm{I}]}{K_i}\right)$$

**Arrhenius** (section 16): heating speeds things up more the higher the barrier.

$$\ln k = \ln A - \frac{E_a}{RT}$$

**pH** (section 16): the enzyme only works if each group of the active site carries the right charge.

$$v = \frac{v_{\max}}{1 + 10^{pK_1 - \mathrm{pH}} + 10^{\mathrm{pH} - pK_2}}$$

---

### 📖 Quick glossary

| Term | Meaning |
|---|---|
| *Substrate* | the molecule that is transformed |
| *Active site* | where the chemistry happens |
| *ES* | enzyme–substrate complex |
| *Transition state* | the summit of the barrier |
| *ΔG‡* | height of the barrier |
| *k*<sub>cat</sub> | reactions per second per enzyme |
| *K*<sub>M</sub> | [S] at half rate |
| *k*<sub>cat</sub>/*K*<sub>M</sub> | efficiency at low [S] |
| *S*<sub>0.5</sub>, *n* | Hill parameters |
| *K*<sub>i</sub> | dissociation constant of the inhibitor |
| *QM/MM* | quantum mechanics for the reactive region and a force field for the rest |
| *NEB, dimer* | methods to find the TS |
| *Imaginary frequency* | the signature of a saddle point |

---

### ✍️ Exercises

**Getting started**

1. Move the `umbral` (threshold) slider in section 5 to 3.0 and to 4.0 Å and run the cell again. How does the fraction of
   near-attack conformations change? What would happen to *k*<sub>cat</sub> if the enzyme did not close its domains?
2. With the Eyring explorer, find out how much the barrier would have to drop to multiply
   *k*<sub>cat</sub> by 1000. Compare it with the difference between the active site in water and in the enzyme.
3. In section 13, move the `k2` slider to 6 s⁻¹ and to 600 s⁻¹ and run the cell again. How do K<sub>M</sub> and V<sub>max</sub> change?
   When is K<sub>M</sub> ≈ K<sub>d</sub> = k<sub>−1</sub>/k<sub>1</sub>?

**Going deeper**

4. Simulate data with n = 1.0, 1.4 and 2.0 in section 14 and fit them with Michaelis–Menten.
   What error do you make if you ignore cooperativity? Which mutant in the table would stop
   releasing insulin at 5 mM?
5. Design an experiment (concentrations of S and I) that tells a competitive inhibitor apart from a
   mixed one with K<sub>i</sub>′ = 3K<sub>i</sub>. To check it, add a new cell ("+ Code" button) and use `cin.fit_inhibition` and the Dixon plot.
6. With the thermochemistry of section 10, how much does ΔG‡ change between 25 and 37 °C? What
   factor in *k* does that mean? Compare it with what the temperature explorer predicts.

**Challenge (advanced, `completo` mode)**

7. In `enzimas/glucoquinasa.py` remove Lys169 from `QM_RESIDUES` (it will become MM charges) and
   recompute the scan. Does the barrier go up or down? It is a "computational mutant"; Zhang et al.
   (2009) obtained 32 kcal/mol for K169A versus 18 for the wild type.

---

### 📚 References
""")

code(r'''
# @title 📚 References
from IPython.display import HTML
import html as _html
refs = ref.get("references", [])
items = "".join(f'<li style="margin:0 0 8px 0;padding-left:4px">{_html.escape(r)}</li>' for r in refs)
viz.mostrar(HTML(
    f'<div style="font-family:Figtree,\'Avenir Next\',\'Segoe UI\',Roboto,Arial,sans-serif;background:{viz.SURFACE};'
    f'border:1px solid {viz.GRID};border-radius:18px;padding:18px 22px;max-width:960px;color:{viz.INK_SECONDARY};font-size:14px;line-height:1.45">'
    f'<div style="font-size:16px;font-weight:650;color:{viz.INK};margin-bottom:10px">📚 Sources of the course data</div>'
    f'<ol style="margin:0;padding-left:22px">{items}</ol>'
    f'<div style="margin-top:12px;font-size:13px;color:{viz.INK_MUTED}">QM/MM tools adapted from the Leonardo suite (Juvenal Yosa, MIT): '
    f'<a href="https://github.com/juvenalyosa/Leonardo" style="color:{viz.PALETTE[0]}">github.com/juvenalyosa/Leonardo</a></div></div>'))
''')
