"""Section 7. QM/MM."""
from ..celdas import code, md

# ============================================================================ 7. QM/MM

md(r"""
## 7. Watching the reaction with quantum mechanics: QM/MM

> 🎯 **In this section** you will see why quantum mechanics is needed to watch the phosphate jump
> from ATP to glucose, why it cannot be used for the whole enzyme, and how the **QM/MM** method
> combines the best of both worlds.

---

### 💡 The analogy: photographing a match

To photograph a football match you do not need to focus on all 50 000 spectators: you focus on
the ball and the nearby players, and the stadium stays in the background. In the enzyme, "the
ball" is the jumping phosphate; the stadium is the rest of the protein and the water.

[[fig:qmmm_estadio | Stadium analogy: the ball is in focus (QM region) and the stands stay in the background (MM region)]]

**Why is molecular dynamics not enough?** The classical molecular dynamics of section 5 cannot
break or form bonds: its springs know nothing about electrons. To see the phosphate jump from ATP
to glucose we need **quantum mechanics (QM)**, which does describe the electrons. But it is very
expensive: impossible for 57 000 atoms.

**The solution: QM/MM.** Only the small piece where the chemistry happens (the **QM region**, a
few dozen atoms) is treated with QM, and the rest with the classical force field (**MM region**).
The two regions "see" each other:

* the charges of the MM region create an electric potential that enters the quantum calculation
  (**electrostatic embedding**), and
* the MM region pushes on the QM region with van der Waals forces (**Lennard‑Jones**): without
  this the QM region would "sink" into the protein, because electrostatics alone only attracts.

---

### 🧮 The equation in words

**Why do we want an equation?** To search for the reaction path, the program needs the energy
(and the forces) of every geometry it tries. In QM/MM that energy has two parts:

$$\text{energy} \;=\; \underbrace{\text{quantum energy of the QM region (in the field of the enzyme)}}_{\text{the in-focus photo}} \;+\; \underbrace{\text{the "elbows" between QM and MM}}_{\text{not walking through the crowd}}$$

Back to the stadium: the players we focus on (the QM region) hear the shouting from the stands
(the charges of the enzyme, which attract or repel their electrons), but they cannot push *into*
the crowd: the elbows of the people around them stop them (the short-range repulsion).

[[fig:qmmm_regiones | QM region surrounded by the MM region with point charges; link atom at the boundary; electrostatics and Lennard-Jones]]

---

### 📐 The full equation, term by term

$$E \;=\; {\color{#eb6834}{E_{\mathrm{PM7}}\big[\text{QM in the potential of the enzyme}\big]}} \;+\; {\color{#4a3aa7}{E_{\mathrm{LJ}}(\text{QM–MM})}}$$

| Term | What it is | In the analogy | In our model |
|---|---|---|---|
| ${\color{#eb6834}{E_{\mathrm{PM7}}[\dots]}}$ | energy of the electrons and nuclei of the QM region, computed with the semiempirical PM7 method | the players in the middle, sharply in focus | 78 atoms, charge −2 |
| potential of the enzyme | the electric field of the MM charges, included *inside* the quantum calculation (embedding) | the music and the shouting of the crowd | ~2300 Amber charges within 16 Å, fixed |
| ${\color{#4a3aa7}{E_{\mathrm{LJ}}}}$ | Lennard‑Jones between QM and MM atoms: weak attraction far away, strong repulsion up close | the elbows of the crowd | Amber parameters |
| link atom | an H that "caps" each covalent bond cut by the boundary | — | side chains cut at the boundary |

**Our recipe (tools from the Leonardo suite, MOPAC PM7 + Amber):**

* **QM region** (78 atoms, charge −2): glucose, the C5'–triphosphate fragment of ATP, the side
  chains of Asp205, Lys169 and Thr228, the Mg²⁺ and its coordinated water.
* Wherever a covalent bond crosses the boundary a **link atom** (an H) is placed.
* The **MM environment** (2300 charges within 16 Å) is kept **fixed** (rigid-environment
  approximation).
* Total energy: $E = E_{\mathrm{PM7}}[\text{QM in the potential of the enzyme}] + E_{\mathrm{LJ}}(\text{QM–MM})$.

---

### 🎛️ What if…

| If… | then… |
|---|---|
| **we remove the potential of the enzyme** | the QM region no longer "sees" the protein: it is like taking it out (section 11 does this, in water) |
| **we remove the Lennard‑Jones term** | electrostatics, which only attracts, "sinks" the QM region into the protein and the geometry stops making sense |
| **we enlarge the QM region** | the model is more faithful, but each calculation costs much more |
| **we let the MM environment move** | the calculation would be more realistic and much more expensive; with a fixed environment, the barrier depends on the chosen conformation (section 10) |

---

### 🔬 The real data: our partition
""")

code(r'''
# @title 🔍 The quantum region inside the enzyme
particion = datos.json_("qmmm/particion.json")
print(f"QM region: {particion['n_qm']} atoms + {particion['n_link']} link H; free: {particion['n_free']}")
print(f"MM environment: {particion['n_mm']} point charges; QM charge = {particion['qm_charge']}; MM charge = {particion['mm_charge']:.2f} e")
print("QM atoms per residue:", pd.Series([a["resname"] for a in particion["qm_atoms"]]).value_counts().to_dict())
display(visor3d.region_qm())
''')

md(r"""
**How do we know the program computes the forces correctly?** By comparing the analytic gradient
with the numerical one (move an atom by 0.005 Å and see how much the energy changes). Every
simulation code should pass this check before it is used:
""")

code(r'''
# @title 🧪 Quality control: analytic versus numerical forces
grad = datos.csv("qmmm/comprobacion_gradiente.csv")
grad["relative error"] = (grad["diff"].abs() / grad["analytic"].abs().clip(lower=1e-6)).map(lambda x: f"{x:.1%}")
display(grad.round(3))
''')

md(r"""
> ✅ **Takeaway.** QM/MM focuses quantum mechanics where the chemistry happens (78 atoms) and
> treats the rest of the enzyme as classical charges and "elbows". The energy has two parts: the
> quantum energy of the QM region **inside the electric field of the enzyme** and the
> Lennard‑Jones repulsion that keeps it from sinking into the protein.
""")
