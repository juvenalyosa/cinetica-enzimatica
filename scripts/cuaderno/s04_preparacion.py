"""Sección 4. preparación."""
from .celdas import code, md

# ============================================================================ 4. preparación
md(r"""
## 4. Preparar la enzima para simularla

> 🎯 **En esta sección** verás por qué una estructura cristalina no se puede simular tal cual y
> qué cuatro arreglos la convierten en un modelo listo para la dinámica molecular.

---

### 💡 La analogía: de la foto a la maqueta

Una foto cristalográfica es como la foto de un edificio: sirve para saber cómo es, pero para
estudiar cómo se mueve con el viento hay que construir una **maqueta completa**. A la foto le
faltan cosas:

* los **hidrógenos** (los rayos X no los ven),
* **trozos de cadena** que estaban desordenados en el cristal y no se ven,
* las **reglas** de cómo se mueve cada átomo (resortes y cargas),
* y el **agua** con iones que la rodea en la célula.

[[fig:preparacion_pasos | Cuatro pasos de preparación: añadir hidrógenos, completar bucles, poner parámetros y sumergir en agua con iones]]

---

### 🔬 Los datos reales: qué hizo el script `scripts/01_preparar_sistema.py`

| Paso | Qué se hizo | Herramienta |
|---|---|---|
| 1. Completar la cadena | Cadena A de 3FGU, residuos 5–458; tres bucles internos faltantes modelados | PDBFixer |
| 2. Recuperar el ATP | AMP‑PNP → ATP (N3B → O3B) | edición del PDB |
| 3. Poner parámetros | glucosa y ATP con GAFF2 y cargas AM1‑BCC; proteína con Amber ff14SB; Mg²⁺ y K⁺ con parámetros de Li–Merz | antechamber, tleap |
| 4. Agua e iones | aguas cristalográficas conservadas; caja de agua TIP3P de 10 Å y 21 Na⁺ para neutralizar | tleap |

> ⚠️ **Ojo.** En modo `completo` esa preparación requiere AmberTools (antechamber, tleap), que no
> está en Colab; por eso los archivos preparados viajan con el repositorio.

Veamos el resultado:
""")

code(r'''
prep = datos.json_("sistema/preparacion.json")
solv = prep["sistema_solvatado"]
print("Fuente:", prep["fuente"])
print("Bucles modelados (índice, residuos):", prep["bucles_modelados"])
print("Campo de fuerza:", prep["campo_de_fuerza"])
print(f"\nSistema solvatado: {solv['n_atoms']:,} átomos, carga total {solv['total_charge']:+.3f} e")
print("Cargas de los ligandos:", {k: round(v, 2) for k, v in solv["ligand_charges"].items()})
print("Residuos:", {k: v for k, v in solv["residue_counts"].items() if k in ("WAT", "Na+", "GLC", "ATP", "MG", "K+")})
''')

md(r"""
> ✅ **Para llevar.** Una estructura cristalina es el punto de partida, no el modelo: hay que
> añadir hidrógenos, completar los trozos que faltan, darle a cada átomo sus resortes y cargas y
> sumergirla en agua con iones. Solo entonces se puede simular.
""")
