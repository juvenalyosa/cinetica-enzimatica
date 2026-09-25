"""Section 0. Cover."""
from ..celdas import code, md

# ============================================================================ 0. cover
md(r"""
# Enzyme kinetics, step by step, with QM/MM simulations

**Enzyme of study: human glucokinase** (hexokinase IV, EC 2.7.1.2)

$$\text{glucose} + \text{ATP} \;\xrightarrow{\;\text{glucokinase}\;}\; \text{glucose‑6‑phosphate} + \text{ADP}$$

[[fig:portada_banner | The glucokinase reaction and the path of the course: structure, molecular dynamics, quantum chemistry and kinetics]]

This notebook explains enzyme kinetics **from scratch**, for undergraduate students. You don't
need to know how to program: just run the cells in order (▶ or `Shift + Enter`) and read.
We use a real enzyme and everything that can be used with it: a crystal structure, molecular
dynamics and quantum-mechanics (QM/MM) simulations, and the numbers measured in the
laboratory for this enzyme, with their sources.

---

### 🧭 How each section is organized

Every section follows the same order, always marked with the same symbols:

| Symbol | Block | What you will find |
|:---:|---|---|
| 🎯 | **In this section** | the question we are going to answer |
| 💡 | **The analogy** | an everyday picture, with a drawing |
| 🧮 | **The equation in words** | why we need an equation and what it says, without symbols |
| 📐 | **The full equation** | the formal version, term by term, with the same colors as the drawing |
| 🎛️ | **What if…** | plots with sliders: move a term and see what changes |
| 🔬 | **The real data** | simulations and values measured in the laboratory |
| ✅ | **Takeaway** | the key idea in two lines |

In addition, the **⚠️ Watch out** boxes warn about a common trap and the **✍️ Think about it** boxes leave you a question.

---

### 🗺️ Contents

| # | Section | What you learn | Tool |
|---|---|---|---|
| | **🧱 Block 1 · The basics** | | |
| 1 | What is an enzyme, and why glucokinase? | basic vocabulary, experimental data | reading |
| 2 | What is the rate of a reaction? | v₀, units, how it is measured | simulation |
| | **🔭 Block 2 · Seeing the enzyme** | | |
| 3 | The 3D structure of the catalytic complex | active site | PDB 3FGU, 3D viewer |
| 4 | Preparing the enzyme for simulation | from the snapshot to the model | PDBFixer, Amber |
| 5 | The enzyme moves: molecular dynamics | reactive conformations | OpenMM |
| | **⚛️ Block 3 · The chemistry, with quantum mechanics** | | |
| 6 | Energy hills and reaction rate | Eyring equation | 🎛️ interactive |
| 7 | Looking at the reaction with quantum mechanics: QM/MM | QM and MM regions | MOPAC PM7 + Amber |
| 8 | Reactant, product and energy profile | reaction coordinate | scan |
| 9 | The transition state | NEB, dimer, imaginary frequency | ASE |
| 10 | From ΔE‡ to ΔG‡: improving the calculation | thermochemistry, averages, method | ASE, MOPAC |
| 11 | The reaction outside the enzyme | SADDLE (QST2), TS, FORCETS | native MOPAC |
| 12 | From the barrier to *k*<sub>cat</sub> | comparison with experiment | Eyring |
| | **📈 Block 4 · Laboratory kinetics** | | |
| 13 | Michaelis–Menten from the mechanism | K<sub>M</sub>, V<sub>max</sub>, k<sub>cat</sub>, k<sub>cat</sub>/K<sub>M</sub> | ODE, fits, 🎛️ |
| 14 | Cooperativity: glucokinase is a sensor | Hill equation, mutants | 🎛️ |
| 15 | Inhibitors and activators | types, K<sub>i</sub>, Dixon, IC₅₀ | 🎛️ |
| 16 | Temperature and pH | Arrhenius, Eyring, pKa | 🎛️ |
| | **🏁 Wrap-up** | | |
| 17 | Summary, glossary and exercises | map of the course | |

> 💡 **Tip.** In Colab, open the **Table of contents** panel (☰ icon on the left) to jump between
> sections.

---

### ⚙️ Before you start

**Two run modes.** In `rapido` mode ("fast", the default) the notebook reads results that were
already computed and draws everything in a few minutes. In `completo` mode ("full") it reruns
the simulations (hours on the Colab CPU). Change `MODO` in the next cell if you want to recompute.

> 🎛️ **How to use the interactive plots.** Cells marked with 🎛️ show sliders.
> Move them and watch how the curve changes; underneath, a sentence translates the numbers.
> If the sliders don't appear (for example when viewing the notebook on GitHub), run the cell.
""")

code(r'''
# @title ⚙️ Set up the environment (run this cell first)
# --- Environment setup (run this cell first) --------------------------------------------
import os, sys, subprocess, pathlib

MODO = "rapido"  # @param ["rapido", "completo"]
# "rapido" (fast): uses precomputed results (minutes) | "completo" (full): recomputes everything (hours)

REPO = "https://github.com/juvenalyosa/cinetica-enzimatica.git"
if "google.colab" in sys.modules:
    # In Colab: download the course or, if it is already there from a previous session, update it to the latest version
    destino = pathlib.Path("/content/cinetica-enzimatica")
    if destino.exists():
        subprocess.run(["git", "-C", str(destino), "fetch", "-q", "origin"], check=True)
        subprocess.run(["git", "-C", str(destino), "reset", "-q", "--hard", "origin/main"], check=True)
    else:
        subprocess.run(["git", "clone", "-q", REPO, str(destino)], check=True)
    os.chdir(destino)
elif not pathlib.Path("enzimas").exists() and pathlib.Path("../enzimas").exists():
    os.chdir("..")                                     # run from notebooks/ in Jupyter
if os.getcwd() not in sys.path:
    sys.path.insert(0, os.getcwd())
# if the cell is run again, forget the version of the package loaded before
for nombre in [m for m in sys.modules if m == "enzimas" or m.startswith("enzimas.")]:
    del sys.modules[nombre]

from enzimas import colab_setup
entorno = colab_setup.instalar(MODO, idioma="en")
''')

code(r'''
# @title 📦 Load the course tools
# --- Imports and plot style -------------------------------------------------------------
%matplotlib inline
%config InlineBackend.figure_format = "retina"
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enzimas import kinetics as cin, viz, datos, interactivo, visor3d, textos
textos.usar("en")                           # language of plots, cards and viewers
viz.apply_style()

def leer(nombre):
    """Reads a text file from data/precalculado (decompresses .gz if needed)."""
    return pathlib.Path(datos.ruta(nombre)).read_text()

ref = cin.glucokinase_reference()          # experimental values from the literature, with their sources
def valor(clave, defecto=None):
    v = ref.get(clave, {})
    return v.get("value", defecto) if isinstance(v, dict) else (v if v is not None else defecto)

print("Ready. Mode:", MODO)
''')
