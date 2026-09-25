"""Section 11. Water."""
from ..celdas import code, md

# ============================================================================ 11. water

md(r"""
## 11. The same reaction **outside** the enzyme

> 🎯 **In this section** you will measure how much of the catalysis is done by *the rest* of the
> protein, by repeating the reaction with the active site immersed in water, and you will meet
> MOPAC's native transition-state finders.

---

### 💡 The analogy: the same song, outside the stadium

How much of a concert's sound is due to the band and how much to the acoustics of the stadium? To
find out, the band plays the same song in a studio. Here the "band" is the active site and the
"stadium" is the rest of the protein.

We take out the active site (the same QM cluster, with its anchors fixed) and put it in **water**
(COSMO implicit solvent), without the electric field of the other 2300 atoms. If the barrier
changes, that change is the work of the protein environment.

[[fig:dentro_fuera | Energy profiles inside the enzyme and of the cluster in water; what the cluster keeps and what it loses]]

---

### 🧮 The equation in words

**Why do we want an equation?** To put a number on the "stadium's contribution":

$$\text{effect of the rest of the protein} \;=\; \text{barrier inside the enzyme} \;-\; \text{barrier of the cluster in water}$$

$$\Delta\Delta E^{\ddagger}_{\text{environment}} \;=\; {\color{#2a78d6}{\Delta E^{\ddagger}_{\text{enzyme}}}} \;-\; {\color{#898781}{\Delta E^{\ddagger}_{\text{water}}}} \;=\; 19.1 - 21.6 \;\approx\; -2.5\ \text{kcal/mol}$$

A **negative** number means the protein environment **lowers** the hill. With the golden rule
(1.36 kcal/mol = ×10), 2.5 kcal/mol means the reaction runs about **70 times** faster inside the
enzyme than in the isolated cluster.

---

### 🧰 The tools: MOPAC's native methods

Outside the enzyme there is no external potential to "freeze", so here we do use MOPAC's native
transition-state finders, with the names under which the Leonardo suite orchestrates them:

| Name in Leonardo | MOPAC keyword | What it does | In the mountain-pass analogy |
|---|---|---|---|
| QST2 | `SADDLE` | Searches for the TS starting from reactant and product (two end points) | walking from both valleys until you meet |
| QST3 / TS | `TS` | Refines the saddle point by following the mode of negative curvature | pinning down the exact position of the pass |
| Validation | `FORCETS` | Frequencies over the free coordinates: there must be one imaginary | checking that only one direction goes down |
| Path | `IRC=1*` | Follows the intrinsic reaction coordinate in both directions | walking down from the pass into both valleys |

A practical detail: `SADDLE` gives an *estimate* of the TS that has to be refined. Here the
refinement was done with the dimer method and then with MOPAC's `TS` keyword, which reached the
same saddle point.

---

### 🔬 The real data
""")

code(r'''
# @title 💧 The active site in water: MOPAC results
agua = res["etapas"]["agua"]
display(viz.tarjetas(
    [("1 · SADDLE (QST2)", f"{agua['barrera_saddle_qst2_kcal']:.1f}", "kcal/mol", f"first TS estimate, ξ = {agua['saddle_xi']:+.2f} Å", "gris"),
     ("2 · Refined TS", f"{agua['barrera_kcal']:.1f}", "kcal/mol", textos.t(agua["metodo_ts"]), "naranja"),
     ("TS geometry", f"{agua['ts_d_PG_O6']:.2f} / {agua['ts_d_PG_O3B']:.2f}", "Å",
      f"Pγ···O6 / Pγ···O3β; O6–H is still {agua['ts_d_O6_H']:.2f} Å", "azul"),
     ("Reaction ΔE", f"{agua['dE_reaccion']:+.1f}", "kcal/mol", "in water the reaction is downhill", "agua")],
    titulo="The same reaction, outside the enzyme (PM7 + COSMO implicit water)"))
conjuntos = [("3 · FORCETS (MOPAC)", agua["frecuencias_mas_bajas"])]
if "frecuencias_ase_mas_bajas" in agua:
    conjuntos.append(("Finer Hessian (ASE)", agua["frecuencias_ase_mas_bajas"]))
fig = viz.plot_frequencies(conjuntos, title="Look at the magnitude, don't just count the imaginary ones",
                           subtitle="Grey: noise from the solvent cavity. Orange: motions that \"fall\". The finer Hessian leaves only one: the reaction")
n_f = max(len(f) for _, f in conjuntos)
tabla_f = pd.DataFrame({"mode": np.arange(1, n_f + 1),
                        **{nombre.split(" (")[0].replace("3 · ", ""): np.r_[np.asarray(f, float), [np.nan] * (n_f - len(f))]
                           for nombre, f in conjuntos}})
viz.mostrar(fig, viz.datos(
    tabla_f, "the lowest frequencies of the TS in water",
    "Each column is one Hessian calculation; each row, one of the slowest modes sorted from lowest to highest. "
    "Values between −40 and +40 cm⁻¹ are numerical noise from the solvent cavity.",
    y=list(tabla_f.columns[1:]), unidades={c: "cm⁻¹" for c in tabla_f.columns[1:]}, formatos={c: "{:.1f}" for c in tabla_f.columns[1:]},
    resaltar={0: ("the reaction", "naranja")}))
''')

md(r"""
**A lesson about validation.** The `FORCETS` Hessian on the COSMO surface shows several small
"imaginary" frequencies (|ν| < 40 cm⁻¹): they are numerical noise from the solvent cavity, not
real motions; that is why it pays to look at the *magnitude* of the frequencies and not just count
them. The Hessian from finer finite differences (ASE) leaves a single clear mode (≈ −160 cm⁻¹, the
jumping phosphorus) and a residual −37 cm⁻¹ attributable to the cavity. MOPAC's `IRC` aborted
because of an internal program error on this solvated surface (that happens in practice too); the
path was obtained with the same descent from the TS as in section 9.
""")

code(r'''
# @title 📈 The reaction path in water
try:
    cam_agua = datos.csv("qmmm/agua_camino_descenso.csv")
    fig = viz.plot_energy_profile(cam_agua["xi"].values, cam_agua["energia_kcal"].values - agua["E_reactivo"], relative=False,
                                  xlabel="ξ (Å)", smooth=False, sort=False, ts_index=int(np.argmax(cam_agua["energia_kcal"].values)),
                                  annotate_states=True, annotate_reaction=True, state_names=("reactant", "TS", "product"),
                                  color=viz.PALETTE[6],
                                  title=f"In water the hill is {agua['barrera_kcal']:.1f} kcal/mol high",
                                  subtitle="Descent path from the refined TS; energies relative to the reactant in water")
    e_rel = cam_agua["energia_kcal"].values - agua["E_reactivo"]
    viz.mostrar(fig, viz.datos(
        pd.DataFrame({"step": cam_agua["cuadro"], "ξ": cam_agua["xi"], "E": e_rel}), "the path in water",
        "Each row is a step of the descent from the TS in water, in the order of the path.",
        x="ξ", y="E", calculadas={"E": "energy of the step − energy of the reactant in water"},
        unidades={"ξ": "Å", "E": "kcal/mol"}, formatos={"ξ": "{:+.2f}", "E": "{:.1f}"},
        resaltar={int(np.argmax(e_rel)): ("TS", "naranja")}, barra="E"))
except FileNotFoundError:
    display(viz.mensaje("Path in water not available.", "ojo"))
''')

md(r"""
The same movie as in section 8, but now **without the rest of the protein**: only the active site
in water. Compare the height of the hill in the plot with the one inside the enzyme.
""")

code(r'''
# @title 🎬 The reaction in water, synchronized with its energy
visor3d.pelicula_reaccion("agua")
''')

code(r'''
# @title ⚖️ Inside the enzyme versus outside it
con = [("E·S", 0.0), ("TS", ts["barrera_kcal"]), ("E·P", P["dE_reaccion_kcal"])]
sin = [("R", 0.0), ("TS", agua["barrera_kcal"]), ("P", agua["dE_reaccion"])]
efecto = ts["barrera_kcal"] - agua["barrera_kcal"]
fig = viz.plot_energy_levels(con, compare=sin, ts_indices=[1], label="inside the enzyme", compare_label="active site in water",
                             title=f"The rest of the protein lowers the hill by {abs(efecto):.1f} kcal/mol",
                             subtitle="Blue/orange: inside the enzyme. Grey: the same active site in water")
viz.mostrar(fig, viz.tarjetas(
    [("Effect of the protein environment", f"{efecto:+.1f}", "kcal/mol", "barrier in the enzyme − barrier in water", "azul"),
     ("In rate", f"× {10 ** (-efecto / 1.364):.0f}", "", "every 1.36 kcal/mol is a factor of 10", "naranja")]),
    viz.datos(pd.DataFrame({"state": ["reactant", "transition state", "product"],
                            "in the enzyme": [e for _, e in con], "in water": [e for _, e in sin],
                            "difference": [a - b for (_, a), (_, b) in zip(con, sin)]}),
              "the two diagrams",
              "The heights of the two diagrams, side by side. The transition-state row is the one that decides the rate.",
              y=["in the enzyme", "in water"], calculadas={"difference": "in the enzyme − in water"},
              unidades={"in the enzyme": "kcal/mol", "in water": "kcal/mol", "difference": "kcal/mol"},
              formatos={"in the enzyme": "{:.1f}", "in water": "{:.1f}", "difference": "{:+.1f}"},
              resaltar={1: ("TS", "naranja")}))
''')

md(r"""
**How to read this comparison.** The cluster "in water" is **not** the uncatalyzed reaction: it
still contains the catalytic base (Asp205), the positive charge of Lys169 and the Mg²⁺, that is,
the essential chemical machinery. That is why the two hills look alike: the difference (a few
kcal/mol) measures what the *rest* of the protein contributes, mainly its electric field and the
fact that it keeps the reactants aligned. The true reaction without enzyme, glucose and ATP alone
in water, is extraordinarily slow: the spontaneous hydrolysis of a phosphate dianion has a
half-life of ~10¹² years (Lad, Williams and Wolfenden, 2003; ΔG‡ ≈ 44 kcal/mol), because no base,
no cation and no electric field are there to help. Between that reaction and k<sub>cat</sub> ≈ 60 s⁻¹
there is a factor of ~10²¹: that is the real size of catalysis.
""")

md(r"""
> ✅ **Takeaway.** The active site alone, in water, already has almost all the chemical machinery
> (Asp205, Lys169, Mg²⁺), which is why its barrier (21.6) resembles the enzyme's (19.1). The rest
> of the protein contributes about 2.5 kcal/mol more (a factor of ~70). Full catalysis, compared
> with the reaction with no help at all, is ~10²¹-fold.
""")
