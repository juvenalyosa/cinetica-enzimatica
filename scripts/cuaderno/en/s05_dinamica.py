"""Section 5. MD."""
from ..celdas import code, md

# ============================================================================ 5. MD
md(r"""
## 5. The enzyme moves: molecular dynamics

> 🎯 **In this section** you will see how a protein is "filmed" with Newton's laws and what that
> movie is good for: measuring how much time the reactants spend well positioned to react.

---

### 💡 The analogy: the photo and the movie

The crystal structure is a **photo**: a single pose. But a protein at 37 °C is never still: it
bends, breathes and vibrates. **Molecular dynamics (MD)** is the **movie**: a sequence of frames
2 femtoseconds apart (0.000000000000002 s). One nanosecond of movie is 500 000 frames.

[[fig:pelicula_md | A photo of the crystal versus a strip of molecular dynamics film with frames every 2 fs]]

---

### 🧮 The equation in words

**Why do we want an equation?** To know where each atom will be in the next frame. The recipe is
that of a pool table: if you know how hard a ball is pushed, you know how it accelerates and
where it will be an instant later.

Each atom is a **ball** with mass, joined to its neighbors by **springs** (the bonds), and it also
feels the **charges** and the **collisions** of all the others. With that, each frame is one turn
of the same cycle:

$$\text{forces} \;\rightarrow\; \text{acceleration} \;\rightarrow\; \text{new position} \;\rightarrow\; \text{forces} \;\rightarrow\; \cdots$$

[[fig:newton_ciclo | Ball-and-spring model and Newton's cycle: positions, forces, acceleration and a 2 fs step]]

---

### 📐 The full equation, term by term

Newton's second law for each atom *i*, and one time step:

$${\color{#eb6834}{\mathbf{F}_i}} = m_i\,{\color{#4a3aa7}{\mathbf{a}_i}}
\qquad\Rightarrow\qquad
\mathbf{v}_i \leftarrow \mathbf{v}_i + {\color{#4a3aa7}{\mathbf{a}_i}}\,{\color{#1baf7a}{\Delta t}},
\qquad
{\color{#2a78d6}{\mathbf{r}_i}} \leftarrow {\color{#2a78d6}{\mathbf{r}_i}} + \mathbf{v}_i\,{\color{#1baf7a}{\Delta t}}$$

| Term | What it is | In the analogy | Value here |
|---|---|---|---|
| ${\color{#2a78d6}{\mathbf{r}_i}}$ | position of atom *i* | where the ball is | ~57 000 atoms |
| ${\color{#eb6834}{\mathbf{F}_i}}$ | total force on the atom | pulls from the springs + attraction/repulsion of the charges + collisions | Amber ff14SB force field |
| $m_i$ | mass of the atom | how heavy the ball is | H = 1, C = 12, O = 16 g/mol |
| ${\color{#4a3aa7}{\mathbf{a}_i}}$ | acceleration | how much its velocity changes | $\mathbf{F}_i/m_i$ |
| $\mathbf{v}_i$ | velocity | where it is heading and how fast | ~ the one set by the temperature (300 K) |
| ${\color{#1baf7a}{\Delta t}}$ | time step | time between frames | 2 fs |

The forces come from an **energy** $U$ that adds up springs, charges and collisions: each atom is
pushed "downhill" on that energy, $\mathbf{F}_i = -\partial U/\partial \mathbf{r}_i$. The program
(OpenMM) uses a Langevin integrator, a more careful version of this same cycle that also keeps
the temperature at 300 K.

---

### 🎛️ What if…

| If… | then… | because… |
|---|---|---|
| **Δt were much larger** | the simulation "blows up" | the atoms would jump too far between frames and the springs would stretch out of control; 2 fs is possible because the lengths of bonds to H are held fixed |
| **we simulate for longer** | we see slower motions (domains opening and closing) | 1 ns is 500 000 turns of the cycle; large changes can take much longer |
| **we wanted to break a bond** | classical MD cannot do it | its springs know nothing about electrons: that requires QM/MM (section 7) |

---

### 🔬 The real data

**What was done.** With OpenMM and Amber: minimization, heating up to 300 K, equilibration at
1 bar and 1 ns of production (script `02_dinamica_molecular.py`). There is no chemistry here (bonds
don't break): MD tells us **how much time the active site spends in a reactive geometry**, with
the Pγ of ATP close to the O6 of glucose. Those geometries are called **near-attack
conformations** (NAC).

[[fig:ataque_cercano | Near-attack conformation with O6 less than 3.5 Å from the phosphorus versus a non-reactive conformation]]

**How to read the plot.** Each point is a frame of the movie (one every 10 ps) and the thick line
is the moving average, which removes the jitter and shows the trend.

* **Top, large:** the distance between the γ phosphorus of ATP and the O6 of glucose. The green
  band is the **near-attack zone** (d < 3.5 Å): while the curve is inside it, the reactants "aim at
  each other". On the right, how many frames fall at each distance and the percentage of time
  inside the zone. The black line is the crystal value.
* **Bottom left:** the RMSD, how far the protein has drifted from the initial structure. A flat
  value of ~1 Å is small: the enzyme is stable and does not fall apart.
* **Bottom right:** the distance from O6 to the catalytic base, Asp205. The green band marks a
  hydrogen bond (< 3.2 Å). Notice that in the classical simulation Asp205 stays **outside** the
  band: we will come back to this in section 10.
""")

code(r'''
# @title 🎬 Molecular dynamics: how the enzyme moves
if MODO == "completo":
    subprocess.run([sys.executable, "scripts/02_dinamica_molecular.py"], check=True)

md_df = datos.csv("md/analisis.csv")
md_info = datos.json_("md/md.json")
fig = viz.plot_md_summary(md_df, crystal_value=2.68, crystal_hbond=2.5)
viz.mostrar(fig, viz.tarjetas(
    [("Movie", f"{md_info['produccion_ns']:g}", "ns", f"{len(md_df)} frames, one every {md_info['cuadro_ps']} ps", "gris"),
     ("Simulated atoms", f"{md_info['n_atomos']:,}".replace(",", " "), "", "protein, ligands, water and ions", "gris"),
     ("Mean RMSD", f"{md_df['rmsd_CA_A'].mean():.2f}", "Å", "the protein barely deforms", "azul"),
     ("Mean d(Pγ–O6)", f"{md_df['d_PG_O6_A'].mean():.2f}", "Å", "crystal: 2.68 Å", "agua")]),
    viz.datos(pd.DataFrame({"time": md_df["tiempo_ps"].round(0), "d(Pγ–O6)": md_df["d_PG_O6_A"],
                            "near attack?": np.where(md_df["d_PG_O6_A"] < 3.5, "yes", "no"),
                            "d(O6···Asp205)": md_df["d_O6_OD1asp205_A"], "RMSD Cα": md_df["rmsd_CA_A"]}),
              "molecular dynamics",
              "Each row is a frame of the movie: the program measures three numbers in each one and draws the "
              "three plots with them. The \"near attack?\" column is the one that gives the percentage in the title.",
              x="time", y=["d(Pγ–O6)", "d(O6···Asp205)", "RMSD Cα"],
              calculadas={"near attack?": "\"yes\" if d(Pγ–O6) < 3.5 Å"},
              unidades={"time": "ps", "d(Pγ–O6)": "Å", "d(O6···Asp205)": "Å", "RMSD Cα": "Å"},
              formatos={"time": "{:.0f}", "d(Pγ–O6)": "{:.2f}", "d(O6···Asp205)": "{:.2f}", "RMSD Cα": "{:.2f}"},
              resaltar={int(md_df["d_PG_O6_A"].idxmin()): ("closest", "agua"), int(md_df["d_PG_O6_A"].idxmax()): ("farthest", "rojo")},
              barra="d(Pγ–O6)"))
''')

md(r"""
### 🎬 The movie

Here is the real movie: 100 frames, one every 10 ps, aligned on the protein. The ribbons are the
enzyme's chain (large domain in blue and small domain in pink); in the center, glucose, ATP and
Mg²⁺. On the right, the Pγ–O6 distance frame by frame: the yellow dot tells you whether, at that
instant, the reactants are "ready to react" (green band).
""")

code(r'''
# @title 🎬 The enzyme movie (rotate it and press ▶)
visor3d.pelicula_md()
''')

code(r'''
# @title 🎯 How much time does the enzyme spend ready to react?
umbral = 3.5  # @param {type:"slider", min:2.8, max:4.5, step:0.1}
# umbral (threshold, Å): practical definition of a near-attack conformation
fraccion = (md_df["d_PG_O6_A"] < umbral).mean()
viz.mostrar(
    viz.tarjetas([("Chosen threshold", f"{umbral:.1f}", "Å", "d(Pγ–O6) below it = \"near attack\"", "gris"),
                  ("Time in near attack", f"{fraccion:.0%}", "", f"{int(round(fraccion * len(md_df)))} of {len(md_df)} frames", "agua"),
                  ("Mean d(Pγ–O6)", f"{md_df['d_PG_O6_A'].mean():.2f}", "Å", None, "azul")]),
    viz.mensaje("The enzyme keeps the reactants \"aiming at each other\" most of the time. Without the enzyme, in water, "
                "that encounter would be a stroke of luck. Move the slider and run the cell again.", "idea"))
''')

md(r"""
> ✅ **Takeaway.** MD is a movie made with Newton's laws, one frame every 2 fs. It does not break
> bonds, but it tells us something key: the enzyme **keeps the reactants aiming at each other**
> most of the time, and only from those near-attack conformations can the chemistry happen.
""")
