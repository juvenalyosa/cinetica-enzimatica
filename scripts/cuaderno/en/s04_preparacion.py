"""Section 4. Preparation."""
from ..celdas import code, md

# ============================================================================ 4. preparation
md(r"""
## 4. Preparing the enzyme for simulation

> 🎯 **In this section** you will see why a crystal structure cannot be simulated as it is and
> which four fixes turn it into a model ready for molecular dynamics.

---

### 💡 The analogy: from the photo to the scale model

A crystallographic structure is like a photo of a building: it tells you what it looks like, but
to study how it sways in the wind you need to build a **complete scale model**. The photo is
missing things:

* the **hydrogens** (X-rays don't see them),
* **pieces of the chain** that were disordered in the crystal and can't be seen,
* the **rules** for how each atom moves (springs and charges),
* and the **water** with ions that surrounds it in the cell.

[[fig:preparacion_pasos | Four preparation steps: add hydrogens, complete loops, assign parameters and immerse in water with ions]]

---

### 🔬 The real data: what the script `scripts/01_preparar_sistema.py` did

| Step | What was done | Tool |
|---|---|---|
| 1. Complete the chain | Chain A of 3FGU, residues 5–458; three missing internal loops modeled | PDBFixer |
| 2. Restore ATP | AMP‑PNP → ATP (N3B → O3B) | PDB editing |
| 3. Assign parameters | glucose and ATP with GAFF2 and AM1‑BCC charges; protein with Amber ff14SB; Mg²⁺ and K⁺ with Li–Merz parameters | antechamber, tleap |
| 4. Water and ions | crystallographic waters kept; 10 Å TIP3P water box and 21 Na⁺ to neutralize | tleap |

> ⚠️ **Watch out.** In `completo` (full) mode this preparation requires AmberTools (antechamber,
> tleap), which is not available in Colab; that is why the prepared files travel with the
> repository.

Let's look at the result:
""")

code(r'''
# @title 🧰 What the prepared system contains
prep = datos.json_("sistema/preparacion.json")
solv = prep["sistema_solvatado"]
cuentas = solv["residue_counts"]
viz.mostrar(
    viz.tarjetas([
        ("atoms in the box", f"{solv['n_atoms']:,}".replace(",", " "), "", "protein, ligands, ions and water", "azul"),
        ("water molecules", f"{cuentas.get('HOH', cuentas.get('WAT', 0)):,}".replace(",", " "), "", "10 Å TIP3P box around it", "agua"),
        ("total charge", f"{solv['total_charge']:+.3f}", "e", f"neutralized with {cuentas.get('Na+', 0)} Na⁺", "violeta"),
        ("modeled loops", str(len(prep["bucles_modelados"])), "", "pieces the crystal could not see", "naranja"),
    ], titulo="The system, ready to simulate",
       nota=f"Source: {prep['fuente']} · Force field: {prep['campo_de_fuerza']}"),
    viz.tabla(pd.DataFrame(
        [("ligands: charges", ", ".join(f"{k} {v:+.2f}" for k, v in solv["ligand_charges"].items())),
         ("residues", ", ".join(f"{k} × {v}" for k, v in cuentas.items() if k in ("HOH", "WAT", "Na+", "GLC", "ATP", "MG", "K+"))),
         ("loops (index, residues)", "; ".join(str(b) for b in prep["bucles_modelados"]))],
        columns=["what", "detail"]), titulo="In detail"),
)
''')

md(r"""
> ✅ **Takeaway.** A crystal structure is the starting point, not the model: you have to add
> hydrogens, complete the missing pieces, give each atom its springs and charges and immerse it in
> water with ions. Only then can it be simulated.
""")
