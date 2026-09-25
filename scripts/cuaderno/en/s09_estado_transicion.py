"""Section 9. TS."""
from ..celdas import code, md

# ============================================================================ 9. TS

md(r"""
## 9. The transition state

> 🎯 **In this section** you will locate the true transition state (a saddle point), **validate**
> it with its vibrational frequencies and check that it connects the reactant with the product.

---

### 💡 The analogy: the mountain pass, again

The scan gave us an approximate hill, but by forcing a single distance. The true transition state
is a **saddle point**: a maximum along the path and a minimum in every other direction, exactly
like the pass between two mountains (if you step off the trail, you go up).

[[fig:punto_silla | Contour map with two valleys and a pass; the NEB beads follow the path and the central one climbs to the saddle point]]

**How it is located (three steps):**

1. **Climbing-image NEB** (*nudged elastic band*): a chain of geometries between reactant and
   product, joined by springs, that relaxes until it traces the minimum-energy path. The highest
   image "climbs" to the summit. It is the modern analogue of the **QST2** method (two end points →
   TS) in Gaussian and MOPAC.
2. **Dimer method**: refines the summit using only gradients, until the force is zero.
3. **Frequencies**: a true saddle point has **exactly one** imaginary frequency, and its vector
   describes the motion of the reaction. This is the definitive test.

---

### 🧮 The test in words: why "imaginary"?

**Why do we want an equation?** At the summit the force is zero… but so it is at the bottom of a
valley. To tell them apart we look at **the shape** of the terrain around the point, and that is
what the vibrations give us. Each vibration is like a **spring**: the steeper the walls, the
faster it oscillates.

$$\text{frequency} \;\propto\; \sqrt{\text{spring stiffness}}$$

In a valley every wall goes up: all the springs are normal. At the saddle, along the path the
terrain goes **down** on both sides: it is an "upside-down spring" (negative stiffness), and the
square root of a negative number is an **imaginary** number. Programs write it as a negative
frequency.

[[fig:frecuencia_imaginaria | In a valley every direction goes up; at the saddle one direction goes down and gives an imaginary frequency]]

---

### 📐 The full equation, term by term

$$\nu \;=\; \frac{1}{2\pi}\sqrt{\frac{{\color{#eb6834}{k}}}{\mu}}$$

| Term | What it is | In the analogy | At the transition state |
|---|---|---|---|
| $\nu$ | vibrational frequency (in cm⁻¹ in chemistry) | how many times per second the ball oscillates | −149 cm⁻¹ in the reaction mode |
| ${\color{#eb6834}{k}}$ | curvature of the energy in that direction | the stiffness of the spring: how steep the walls are | **negative** along the path |
| $\mu$ | effective mass that moves | how heavy the ball is | mostly the phosphorus and its oxygens |

---

### 🎛️ What if…

| If the frequency calculation gives… | then the geometry is… |
|---|---|
| **no** imaginary frequency | a **minimum** (reactant, product or an intermediate), not a TS |
| **exactly one** | a **first-order saddle point**: a transition state ✔ |
| **two or more** | a summit in several directions (higher-order saddle): keep searching |

---

### 🔬 The real data
""")

code(r'''
# @title ⛰️ From the scan to the minimum-energy path (NEB)
neb = datos.csv("qmmm/neb.csv")
perfiles = [(xi_esc, e_esc, "restrained scan"),
            (np.r_[R["xi"], neb["xi"].values[1:]], np.r_[0.0, neb["energia_rel_kcal"].values[1:]], "NEB (climbing image)")]
fig = viz.plot_energy_profiles(perfiles, xlabel="ξ (Å)", relative=False, colors=[viz.INK_MUTED, viz.COLORS["reactivo"]],
                               title="The NEB finds the same summit along a more natural path",
                               subtitle="Grey: scan that forces a single distance. Blue: NEB, which relaxes all coordinates at once")
tabla_neb = pd.DataFrame({"image": neb["imagen"] if "imagen" in neb else np.arange(len(neb)),
                          "ξ": np.r_[R["xi"], neb["xi"].values[1:]], "E": np.r_[0.0, neb["energia_rel_kcal"].values[1:]]})
viz.mostrar(fig, viz.datos(
    tabla_neb, "the NEB (the chain of images)",
    "Each row is a \"bead\" of the elastic chain that joins reactant and product. At the end of the optimization the highest "
    "image has climbed to the summit: that is the estimate of the transition state.",
    x="ξ", y="E", unidades={"ξ": "Å", "E": "kcal/mol"}, formatos={"ξ": "{:+.2f}", "E": "{:.1f}"},
    resaltar={0: ("R", "azul"), int(tabla_neb["E"].idxmax()): ("climbing", "naranja"), len(tabla_neb) - 1: ("P", "agua")},
    barra="E"))
''')

code(r'''
# @title ✅ The transition state and its imaginary frequency
ts = res["etapas"].get("dimero")
if ts is None or not ts.get("convergido", True):
    display(viz.mensaje("The dimer did not converge: the NEB climbing image is used as the TS.", "ojo"))
    ts = ts or {**res["etapas"]["neb"], "xi": neb.loc[res["etapas"]["neb"]["imagen_ts"], "xi"], "d_PG_O6": float("nan"), "d_PG_O3B": float("nan")}
fr = res["etapas"].get("frecuencias", {})
val = fr["validacion"]
freqs = np.asarray(fr["frecuencias_mas_bajas"])
tarjetas_ts = viz.tarjetas(
    [("Barrier ΔE‡", f"{ts['barrera_kcal']:.1f}", "kcal/mol", "from the relaxed reactant (dimer method)", "naranja"),
     ("ξ at the summit", f"{ts['xi']:+.2f}", "Å", "almost zero: halfway", "naranja"),
     ("Pγ ··· O6 / Pγ ··· O3β", f"{ts['d_PG_O6']:.2f} / {ts['d_PG_O3B']:.2f}", "Å", "the phosphorus, equidistant from both oxygens", "azul"),
     ("Imaginary frequencies", f"{val['imaginary_mode_count']}", "", "first-order saddle point ✔" if val["ok"] else "check ✘",
      "verde" if val["ok"] else "rojo")],
    titulo="The transition state, located and validated")
fig = viz.plot_frequencies([("TS in the enzyme", freqs)],
                           title=f"A single imaginary frequency ({freqs.min():.0f} cm⁻¹): it is a saddle point".replace("-", "−"),
                           subtitle="The six slowest vibrations of the transition state; all the others are real (positive)")
viz.mostrar(tarjetas_ts, fig, viz.datos(
    pd.DataFrame({"mode": np.arange(1, len(freqs) + 1), "frequency": freqs,
                  "type": ["imaginary: the geometry \"falls\"" if f < 0 else "real: the geometry vibrates" for f in freqs]}),
    "the lowest frequencies of the TS",
    "The frequencies come from the matrix of second derivatives of the energy (the Hessian). By convention, a negative "
    "curvature is written as a negative (\"imaginary\") frequency.",
    y="frequency", unidades={"frequency": "cm⁻¹"}, formatos={"frequency": "{:.1f}"},
    resaltar={int(np.argmin(freqs)): ("the reaction", "naranja")}))
''')

md(r"""
**What the geometry of the transition state tells us.** At the TS the phosphorus is at the same
distance from both oxygens (≈ 2.1–2.2 Å from O6 and from O3β): it is a **concerted, "in-line"**
TS; the phosphate passes from one oxygen to the other, flipping its three oxygens like an umbrella.
The O6 proton, on the other hand, **has not moved yet** (it is still ~1.07 Å from O6). Only after
crossing the summit, once the P–O6 bond has formed, does the proton jump to Asp205: the path goes
down a "plateau" until ξ ≈ 1 Å and only then finishes falling. In this model Asp205 acts as a base
**after** the phosphoryl transfer. An independent QM/MM study (Zhang et al., 2009) describes Asp205
as a general base and Lys169 as a general acid that protonates the phosphate; the fact that two
models disagree on the *order* of the steps but agree on the barrier (section 12) is a good
example of why mechanisms are studied with both simulations and mutants.

**The motion of the transition state.** The imaginary frequency corresponds to a vibration that
does not oscillate but "falls" towards the reactant or towards the product. Animating it shows us
the chemistry: watch the phosphorus go back and forth between the O3β of ATP and the O6 of glucose,
and the dot on the plot cross from the "towards the reactant" zone to the "towards the product" zone:
""")

code(r'''
# @title 🎞️ The motion of the transition state
visor3d.modo_imaginario()
''')

md(r"""
**Coming down from the summit.** If we push the TS a tiny bit to each side along the imaginary
mode and let the geometry relax, we descend along the minimum-energy path to the reactant and to
the product. This confirms that **this** TS connects **these** reactants with **these** products
(the idea of the IRC, *intrinsic reaction coordinate*).
""")

code(r'''
# @title ⬇️ Down from the summit to reactant and product
camino = datos.csv("qmmm/camino_descenso.csv")
fig = viz.plot_energy_profile(camino["xi"].values, camino["energia_rel_kcal"].values, relative=False, xlabel="ξ (Å)",
                              ts_index=int(camino["energia_rel_kcal"].idxmax()), smooth=False, sort=False, annotate_states=True,
                              state_names=("reactant", "TS", "end of descent"),
                              title="From the summit you come down to both valleys: this TS connects R with P",
                              subtitle="Each point is a relaxed geometry; the plateau after the summit is the proton moving to Asp205")
viz.mostrar(fig, viz.datos(
    camino.rename(columns={"cuadro": "step", "energia_rel_kcal": "E"})[["step", "xi", "E"]].rename(columns={"xi": "ξ"}),
    "the descent from the summit",
    "Each row is a step of the descent: starting from the TS the geometry is nudged a little to each side and allowed to fall. "
    "The steps are in the order of the path, from one valley to the other through the summit.",
    x="ξ", y="E", unidades={"ξ": "Å", "E": "kcal/mol"}, formatos={"ξ": "{:+.2f}", "E": "{:.1f}"},
    resaltar={int(camino["energia_rel_kcal"].idxmax()): ("TS", "naranja")}, barra="E"))
''')

code(r'''
# @title 📊 Energy diagram of the reaction
niveles = [("E·S (reactant)", 0.0), ("TS", ts["barrera_kcal"]), ("E·P (product)", P["dE_reaccion_kcal"])]
fig = viz.plot_energy_levels(niveles, ts_indices=[1], title="The reaction in three numbers",
                             subtitle="PM7/Amber, fixed environment: relative electronic energies (ΔE), not yet free energies")
viz.mostrar(fig, viz.datos(
    pd.DataFrame({"state": [n for n, _ in niveles], "E": [e for _, e in niveles]}), "the three levels",
    "The three heights of the diagram. The barrier is the difference TS − reactant; the reaction energy, product − reactant.",
    y="E", unidades={"E": "kcal/mol"}, formatos={"E": "{:.1f}"},
    resaltar={0: ("R", "azul"), 1: ("TS", "naranja"), 2: ("P", "agua")}))
''')

md(r"""
> ✅ **Takeaway.** The transition state is a **saddle point**: it is located in steps
> (scan → NEB → dimer) and **validated** with two tests: exactly one imaginary frequency (here
> −149 cm⁻¹, the phosphorus going back and forth) and a path that, coming down from the summit,
> reaches the reactants on one side and the products on the other.
""")
