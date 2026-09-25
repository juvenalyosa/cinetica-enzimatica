"""Section 3. Structure."""
from ..celdas import code, md

# ============================================================================ 3. structure
md(r"""
## 3. The 3D structure of the catalytic complex

> 🎯 **In this section** you will see where every atom of glucokinase is just before the
> reaction, what each piece of the active site does and why the reactants are already "aiming at
> each other".

---

### 💡 The analogy: a clamshell that closes over glucose

To simulate an enzyme we need to know where every atom is. Crystallographers give us that: the
structure **3FGU** in the Protein Data Bank is a "snapshot" of human glucokinase with
**glucose**, an ATP analog and the **Mg²⁺** ion, all inside the active site, just before the
reaction.

Glucokinase has two parts, the **large domain** and the **small domain**, joined by a hinge.
Without glucose the enzyme is **open**, like a half-open clamshell. When glucose enters, the
small domain closes over it: the reactants are enclosed, with no water around them and lined up to
react.

[[fig:almeja_dominios | Glucokinase open and closed: the small domain closes over glucose like a clamshell]]

---

### 🔬 The real data: who is who in the active site

In the closed form, each piece has a role:

| Piece | Where it is | What it does | In the analogy |
|---|---|---|---|
| **O6 of glucose** | ~2.7 Å from the phosphorus Pγ of ATP | will receive the phosphate | the hand that catches the ball |
| **Pγ of ATP** | at the end of the triphosphate | is the phosphorus that jumps | the ball |
| **Asp205** | 2.5 Å from the O6–H hydroxyl | **catalytic base**: it will take the proton | the one who frees the hand to catch |
| **Lys169** and **Mg²⁺** | next to the phosphates | neutralize their negative charges | shock absorbers |

A P–O bond is 1.6 Å long: at 2.7 Å, O6 and Pγ are **very close, but not yet bonded**.

[[fig:sitio_activo | Scheme of the glucokinase active site in the 3FGU crystal with the key distances]]

Explore the real enzyme, atom by atom: **drag** to rotate it, use the **scroll wheel** to zoom in
and the buttons to see the **active site** or the **surface**. The **dashed lines** mark the two
distances that tell us whether everything is ready to react (Pγ ··· O6 and O6 ··· Asp205); their
values are computed in the next cell:
""")

code(r'''
# @title 🧬 Glucokinase in 3D (rotate it with the mouse)
pdb_cristal = pathlib.Path("data/raw/3FGU.pdb").read_text()
visor3d.complejo_cristal(pdb_cristal)
''')

md(r"""
And let's check the two key distances directly in the crystal file:
""")

code(r'''
# @title 📏 Key distances in the crystal
# Key distances in the crystal (Å): is everything ready to react?
def coords_pdb(texto, resn, name, resi=None):
    for l in texto.splitlines():
        if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == resn and l[12:16].strip() == name and l[16] in " A":
            if resi is None or int(l[22:26]) == resi:
                return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])

PG, O6, OD1 = coords_pdb(pdb_cristal, "ANP", "PG"), coords_pdb(pdb_cristal, "BGC", "O6"), coords_pdb(pdb_cristal, "ASP", "OD1", 205)
viz.mostrar(
    viz.tarjetas([("Pγ ··· O6 (glucose)", f"{np.linalg.norm(PG - O6):.2f}", "Å",
                   "a P–O bond is 1.6 Å long: they are very close, but not yet bonded", "agua"),
                  ("O6 ··· OD1 (Asp205)", f"{np.linalg.norm(O6 - OD1):.2f}", "Å",
                   "short hydrogen bond: Asp205 is ready to take the proton", "azul")],
                 titulo="The two key distances, measured in the crystal file"),
    viz.mensaje("These are the same dashed lines as in the 3D viewer above (aqua and blue). Press \"Active site\" "
                "to zoom in and rotate the molecule: you will see that the O6 of glucose points straight at the γ phosphorus of ATP.",
                "idea"))
''')

md(r"""
> ✍️ **Think about it.** The crystal uses AMP‑PNP, in which the oxygen between the β and γ
> phosphorus atoms was replaced by a nitrogen (N3B). That change prevents the reaction and makes it
> possible to capture the complex. To simulate the real reaction, we will put the oxygen back:
> true ATP.

[[fig:ampnp_atp | AMP-PNP with a nitrogen bridge that does not break versus ATP with an oxygen bridge that does break]]

> ✅ **Takeaway.** Glucokinase closes like a clamshell over glucose and leaves O6 2.7 Å from the
> phosphorus it must attack, with Asp205 ready to take its proton and Mg²⁺ and Lys169 holding the
> phosphate charges. The snapshot already shows an enzyme **primed to react**.
""")
