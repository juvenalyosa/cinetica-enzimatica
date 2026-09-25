"""Section 8. R, P and scan."""
from ..celdas import code, md

# ============================================================================ 8. R, P and scan

md(r"""
## 8. Reactant, product and the energy profile

> 🎯 **In this section** you will describe the whole reaction with **a single number**, the
> reaction coordinate ξ, and obtain the first energy profile of the phosphate transfer inside the
> enzyme.

---

### 💡 The analogy: the ribbon in a tug of war

In a tug of war, the ribbon tied to the middle of the rope tells you who is winning: if it is on
one team's side, that team is ahead; if it is in the middle, it is a draw. Here the rope is the
phosphate group (the Pγ), and the two "teams" are the O3β oxygen of ATP (which lets it go) and the
O6 oxygen of glucose (which receives it).

First we relax the QM region at its starting point (the **reactant**: glucose + ATP). Then we
build the **product** (glucose‑6‑phosphate + ADP, with the O6 proton already on Asp205) and relax
it too. The energy difference tells us whether the reaction is "downhill" or "uphill" inside the
enzyme.

[[fig:xi_coordenada | Three panels of the phosphate moving from the O3β of ATP to the O6 of glucose, with the two distances and the value of ξ]]

---

### 🧮 The equation in words

**Why do we want an equation?** A reaction moves dozens of atoms at once; to draw an energy
profile we need a horizontal axis, a number that advances continuously from the reactant to the
product.

$$\xi \;=\; \text{how far the phosphorus has moved away from ADP} \;-\; \text{how far it still is from glucose}$$

---

### 📐 The full equation, term by term

$$\xi \;=\; {\color{#2a78d6}{d(\mathrm{P_\gamma\!-\!O_{3\beta}})}} \;-\; {\color{#1baf7a}{d(\mathrm{P_\gamma\!-\!O_6})}}$$

The **blue** term is the bond that is breaking (reactant side); the **aquamarine** one, the bond
that is forming (product side).

| Situation | d(Pγ–O3β) | d(Pγ–O6) | ξ | In the analogy | In our calculation |
|---|---|---|---|---|---|
| reactant (phosphate on ATP) | short (1.6 Å) | long (3.2 Å) | negative | the ribbon on the ATP side | ξ = −1.56 Å |
| transition state | intermediate | intermediate | ≈ 0 | the ribbon in the middle | ξ = +0.12 Å |
| product (phosphate on glucose) | long | short (1.7 Å) | positive | the ribbon on the glucose side | ξ = +1.75 Å |

---

### 🎛️ What if…

| If… | then… | because… |
|---|---|---|
| **the phosphorus moves 0.1 Å** towards glucose, in line | ξ goes up by about **0.2 Å** | one distance grows while the other shrinks at the same time |
| **ξ = 0** | the phosphorus is equally far from both oxygens, but that **does not guarantee** it is the summit | the summit is the energy maximum, which is searched for separately (section 9) |
| **we fix ξ** at an intermediate value and let everything else relax | we get one point of the **relaxed scan** | the restraint holds the geometry "halfway" |

A **relaxed scan** fixes ξ at intermediate values (with a harmonic restraint) and lets everything
else adjust: the result is an approximate **energy profile**, and its highest point is a first
estimate of the transition state.

---

### 🔬 The real data
""")

code(r'''
# @title ⚗️ Reactant and product inside the enzyme
if MODO == "completo":
    subprocess.run([sys.executable, "scripts/03_qmmm_reaccion.py"], check=True)

res = datos.json_("qmmm/resumen.json")
R, P = res["etapas"]["reactivo"], res["etapas"]["producto"]
viz.tarjetas(
    [("Reactant · ξ", f"{R['xi']:+.2f}".replace("-", "−"), "Å",
      f"Pγ–O3β {R['d_PG_O3B']:.2f} Å (bonded) · Pγ–O6 {R['d_PG_O6']:.2f} Å (free)", "azul"),
     ("Product · ξ", f"{P['xi']:+.2f}", "Å",
      f"Pγ–O6 {P['d_PG_O6']:.2f} Å (bonded) · Pγ–O3β {P['d_PG_O3B']:.2f} Å (free)", "agua"),
     ("The O6 proton", f"{P['d_OD1_H']:.2f}", "Å",
      f"in the product it sits on Asp205 (in the reactant, O6–H = {R['d_O6_H']:.2f} Å)", "violeta"),
     ("Reaction ΔE", f"{P['dE_reaccion_kcal']:+.1f}", "kcal/mol", "uphill inside the enzyme (PM7, fixed environment)", "naranja")],
    titulo="The two ends of the reaction, relaxed in the enzyme")
''')

code(r'''
# @title 📈 The energy profile of the reaction
escaneo = datos.csv("qmmm/escaneo.csv")
xi_esc = np.r_[R["xi"], escaneo["xi"].values]                   # the profile starts at the relaxed reactant (E = 0)
e_esc = np.r_[0.0, escaneo["energia_rel_kcal"].values]
fig = viz.plot_energy_profile(xi_esc, e_esc, relative=False, xlabel="ξ = d(Pγ–O3β) − d(Pγ–O6)  (Å)",
                              ts_index=int(np.argmax(e_esc)), annotate_states=True,
                              state_names=("reactant\nphosphate on ATP", "summit: 1st TS estimate", "product\nphosphate on glucose"),
                              title=f"To move the phosphate you must climb a {e_esc.max():.0f} kcal/mol hill",
                              subtitle="Relaxed scan: ξ is fixed at each point and everything else relaxes (PM7 in the field of the enzyme)")
tabla_esc = pd.DataFrame({
    "point": ["reactant"] + [str(int(p)) for p in escaneo["punto"]],
    "fixed ξ": np.r_[np.nan, escaneo["xi_objetivo"].values],
    "Pγ–O3β": np.r_[R["d_PG_O3B"], escaneo["d_PG_O3B"].values],
    "Pγ–O6": np.r_[R["d_PG_O6"], escaneo["d_PG_O6"].values],
    "ξ": xi_esc,
    "O6–H": np.r_[R.get("d_O6_H", np.nan), escaneo["d_O6_H"].values],
    "E": e_esc})
k_ts = int(np.argmax(e_esc))
viz.mostrar(fig, viz.datos(
    tabla_esc, "the reaction scan",
    "Each row is an optimized geometry. The program fixes ξ at a value (\"fixed ξ\"), relaxes all the other atoms and "
    "measures the energy. Notice how Pγ–O3β (the breaking bond) gets longer while Pγ–O6 (the forming bond) gets shorter.",
    x="ξ", y="E", calculadas={"ξ": "d(Pγ–O3β) − d(Pγ–O6)", "E": "energy of the point − energy of the relaxed reactant"},
    unidades={"fixed ξ": "Å", "Pγ–O3β": "Å", "Pγ–O6": "Å", "ξ": "Å", "O6–H": "Å", "E": "kcal/mol"},
    formatos={c: "{:.2f}" for c in ("fixed ξ", "Pγ–O3β", "Pγ–O6", "ξ", "O6–H")} | {"E": "{:.1f}"},
    resaltar={0: ("R", "azul"), k_ts: ("summit", "naranja"), len(tabla_esc) - 1: ("P", "agua")}, barra="E"))
''')

md(r"""
> **A detail worth noticing.** The profile starts with a jump of ~4 kcal/mol between the reactant
> (ξ = −1.56 Å, E = 0) and the first point of the scan. The initial optimization of the reactant
> stopped on a *shoulder* of the surface; the true minimum appeared later, when descending from the
> transition state (section 9). Practical lesson: a "converged" optimizer does not guarantee the
> lowest minimum; that is why reaction paths are checked in both directions.
""")

md(r"""
### 🎬 Watch it happen

Now the movie: the atoms on the left and, on the right, the **energy profile with a yellow dot
that follows the movie**. Press ▶ or drag the bar. Watch three things:

1. The **pink** bond (Pγ–O3β, on ATP) gets thinner and turns dashed: **it is breaking**.
2. The **aqua** bond (Pγ–O6, towards glucose) appears and thickens: **it is forming**.
3. At the summit of the plot (the transition state) both are half-formed: the phosphorus sits
   **between** the two oxygens.

Only real bonds are drawn: those of the model's topology, plus these four changing ones, with a
thickness proportional to their bond order (solid when formed, dashed when half-formed or
half-broken). The thin green lines are the coordination of Mg²⁺, which is not a covalent bond.
The movie follows the minimum-energy path computed in section 9, from the relaxed reactant to the
product.
""")

code(r'''
# @title 🎬 The reaction in motion, synchronized with its energy
visor3d.pelicula_reaccion("enzima")
''')

md(r"""
> ✅ **Takeaway.** ξ = d(Pγ–O3β) − d(Pγ–O6) sums up the reaction in a single number: negative
> with the phosphate on ATP, close to zero at the summit, positive with the phosphate on glucose.
> A relaxed scan along ξ gives the first energy profile and a first estimate of the barrier,
> which section 9 refines.
""")
