#!/usr/bin/env python
"""Genera notebooks/Cinetica_Enzimatica_Glucoquinasa.ipynb a partir de las celdas definidas aquí.

Mantener el notebook como código facilita revisarlo, versionarlo y regenerarlo:
    python scripts/construir_notebook.py
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "notebooks" / "Cinetica_Enzimatica_Glucoquinasa.ipynb"

CELLS = []


def md(text):
    cell = nbf.v4.new_markdown_cell(text.strip("\n"))
    cell["id"] = f"celda-{len(CELLS):03d}"  # id determinista: el notebook se regenera byte a byte igual
    CELLS.append(cell)


def code(text):
    cell = nbf.v4.new_code_cell(text.strip("\n"))
    cell["id"] = f"celda-{len(CELLS):03d}"
    CELLS.append(cell)


# ============================================================================ 0. portada
md(r"""
# Cinética enzimática, paso a paso, con simulaciones QM/MM

**Enzima de estudio: la glucoquinasa humana** (hexoquinasa IV, EC 2.7.1.2)

> glucosa + ATP  ⟶  glucosa‑6‑fosfato + ADP

Este cuaderno explica la cinética enzimática **desde cero**, para estudiantes de pregrado.
Cada sección sigue el mismo orden: primero **una analogía o una idea simple**, después **la
ecuación con el significado de cada término**, luego **qué pasa si un término sube o baja**
(con gráficas interactivas que puedes mover), y por último **los datos reales**: una
estructura cristalina, simulaciones de dinámica molecular y de mecánica cuántica (QM/MM), y
los números medidos en el laboratorio para esta enzima, con sus fuentes.

**Índice**

| # | Sección | Qué aprendes | Herramienta |
|---|---|---|---|
| 1 | ¿Qué es una enzima y por qué la glucoquinasa? | vocabulario básico, datos experimentales | lectura |
| 2 | ¿Qué es la velocidad de una reacción? | v₀, unidades, cómo se mide | simulación |
| 3 | La estructura 3D del complejo catalítico | sitio activo | PDB 3FGU, visor 3D |
| 4 | Preparar la enzima para simularla | de la foto al modelo | PDBFixer, Amber |
| 5 | La enzima se mueve: dinámica molecular | conformaciones reactivas | OpenMM |
| 6 | Colinas de energía y velocidad | ecuación de Eyring | 🎛️ interactivo |
| 7 | Mirar la reacción con mecánica cuántica: QM/MM | regiones QM y MM | MOPAC PM7 + Amber |
| 8 | Reactivo, producto y perfil de energía | coordenada de reacción | escaneo |
| 9 | El estado de transición | NEB, dímero, frecuencia imaginaria | ASE |
| 10 | De ΔE‡ a ΔG‡: mejoras del cálculo | termoquímica, promedios, método | ASE, MOPAC |
| 11 | La reacción fuera de la enzima | SADDLE (QST2), TS, FORCETS | MOPAC nativo |
| 12 | De la barrera a *k*<sub>cat</sub> | comparación con el experimento | Eyring |
| 13 | Michaelis–Menten desde el mecanismo | K<sub>M</sub>, V<sub>max</sub>, k<sub>cat</sub>, k<sub>cat</sub>/K<sub>M</sub> | ODE, ajustes, 🎛️ |
| 14 | Cooperatividad: la glucoquinasa es un sensor | ecuación de Hill, mutantes | 🎛️ |
| 15 | Inhibidores y activadores | tipos, K<sub>i</sub>, Dixon, IC₅₀ | 🎛️ |
| 16 | Temperatura y pH | Arrhenius, Eyring, pKa | 🎛️ |
| 17 | Resumen, glosario y ejercicios | | |

**Dos modos de ejecución.** En modo `rapido` (por defecto) el cuaderno lee resultados ya
calculados y dibuja todo en unos minutos. En modo `completo` reejecuta las simulaciones
(horas en la CPU de Colab). Cambia `MODO` en la celda siguiente si quieres recalcular.

> **Cómo usar las gráficas interactivas.** Las celdas marcadas con 🎛️ muestran deslizadores.
> Muévelos y observa cómo cambia la curva; debajo aparece una frase que traduce los números.
> Si los deslizadores no aparecen (por ejemplo al ver el cuaderno en GitHub), ejecuta la celda.
""")

code(r'''
# --- Preparación del entorno (ejecuta esta celda primero) -------------------------------
import os, sys, subprocess, pathlib

MODO = "rapido"          # "rapido": usa resultados precalculados | "completo": recalcula todo

if not pathlib.Path("enzimas").exists():
    if pathlib.Path("../enzimas").exists():            # ejecutado desde notebooks/
        os.chdir("..")
    else:                                              # Google Colab: clonar el repositorio
        subprocess.run(["git", "clone", "-q", "https://github.com/juvenalyosa/cinetica-enzimatica.git"], check=True)
        os.chdir("cinetica-enzimatica")
sys.path.insert(0, os.getcwd())

from enzimas import colab_setup
entorno = colab_setup.instalar(MODO)
''')

code(r'''
# --- Importaciones y estilo gráfico -----------------------------------------------------
%matplotlib inline
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enzimas import kinetics as cin, viz, datos, interactivo
viz.apply_style()

def leer(nombre):
    """Lee un archivo de texto de data/precalculado (descomprime .gz si hace falta)."""
    return pathlib.Path(datos.ruta(nombre)).read_text()

ref = cin.glucokinase_reference()          # valores experimentales de la literatura, con sus fuentes
def valor(clave, defecto=None):
    v = ref.get(clave, {})
    return v.get("value", defecto) if isinstance(v, dict) else (v if v is not None else defecto)

print("Listo. Modo:", MODO)
''')

# ============================================================================ 1. enzimas
md(r"""
## 1. ¿Qué es una enzima y por qué la glucoquinasa?

**La analogía.** Imagina que quieres cruzar una cordillera. Sin guía, buscarías el paso a
ciegas y tardarías semanas. Un buen guía conoce el **paso más bajo** y te lleva por él en
horas. La montaña sigue ahí, pero el camino es otro. Una **enzima** es ese guía: no cambia
el punto de partida ni el de llegada de una reacción química, pero la lleva por un camino
mucho más fácil, y por eso ocurre miles o millones de veces más rápido.

**Vocabulario que usaremos todo el tiempo.**

* **Sustrato (S)**: la molécula que la enzima transforma. Aquí, glucosa y ATP.
* **Producto (P)**: lo que sale. Aquí, glucosa‑6‑fosfato y ADP.
* **Sitio activo**: el hueco de la enzima donde encajan los sustratos y ocurre la química.
* **Complejo enzima‑sustrato (ES)**: la enzima con el sustrato ya dentro.
* **Estado de transición (TS)**: el "paso de montaña", el momento más difícil de la reacción.
  Casi todo este cuaderno gira en torno a él.

**Nuestra enzima.** La **glucoquinasa** toma una glucosa y le pega un grupo fosfato que
viene del ATP. El producto, glucosa‑6‑fosfato, ya no puede salir de la célula: es el primer
paso para "atrapar" y usar la glucosa. La glucoquinasa vive en el hígado y en las células β
del páncreas, donde funciona como **sensor de glucosa**: su actividad decide cuánta insulina
se libera. Mutaciones en su gen causan una diabetes hereditaria (GCK‑MODY, antes MODY2) y,
si la vuelven hiperactiva, hipoglucemia congénita.

Tres cosas la hacen ideal para aprender cinética:

1. Su reacción es una **transferencia de fosforilo** clásica, con Mg²⁺ y una base catalítica (Asp205).
2. Su curva de velocidad **no** es la hipérbola de Michaelis–Menten sino una **sigmoide** (sección 14).
3. Tiene inhibidores, una proteína reguladora y activadores farmacológicos bien estudiados (sección 15).

**Números medidos en el laboratorio.** Todo lo que sigue se compara con estos valores, tomados
de artículos originales (columna *fuente*). Guárdalos: son el "resultado experimental" que la
teoría debe explicar. Fíjate en que distintos laboratorios obtienen k<sub>cat</sub> entre 38 y 66 s⁻¹
para la misma enzima: en bioquímica, un número siempre viene con sus condiciones.
""")

code(r'''
filas = []
for k, v in ref.items():
    if isinstance(v, dict) and "value" in v:
        filas.append((v.get("label", k), v["value"], v.get("unit", ""), v.get("conditions", ""), v.get("source", "")))
tabla = pd.DataFrame(filas, columns=["parámetro", "valor", "unidad", "condiciones", "fuente"])
pd.set_option("display.max_colwidth", 90)
display(tabla)
print(ref.get("notes", ""))
''')

# ============================================================================ 2. velocidad
md(r"""
## 2. ¿Qué es la velocidad de una reacción?

**La idea simple.** La velocidad de una reacción es **cuánto producto aparece por segundo**
(o cuánto sustrato desaparece). Si en un tubo con enzima la glucosa‑6‑fosfato pasa de 0 a
10 µM en 10 s, la velocidad es 1 µM/s.

**Por qué "velocidad inicial".** Al principio del experimento el sustrato apenas se ha gastado
y el producto apenas se ha acumulado, así que la velocidad es constante y refleja solo la
enzima y la concentración de sustrato que pusimos. Después la curva se dobla (el sustrato se
agota). Por eso los bioquímicos miden la **velocidad inicial, v₀**: la pendiente de la curva
de producto al comienzo.

**Cómo se mide en la glucoquinasa.** No se ve la glucosa‑6‑fosfato directamente: se acopla
otra enzima (glucosa‑6‑fosfato deshidrogenasa) que la convierte y, al hacerlo, produce NADPH,
que absorbe luz a 340 nm. El espectrofotómetro registra la absorbancia frente al tiempo, y la
pendiente inicial es v₀ (el ensayo "acoplado a G6PDH" de casi todos los artículos de la tabla
anterior). Simulemos esa curva:
""")

code(r'''
# Curva de progreso simulada: producto frente a tiempo para tres concentraciones de sustrato
E0, k1, k_1, k2 = 0.05, 1.0, 50.0, 60.0          # µM, µM⁻¹s⁻¹, s⁻¹, s⁻¹
fig, ax = viz.figure(7.5, 4.2)
for S0, color in zip((50, 200, 1000), viz.sequential_blue(3)):
    sim = cin.simulate_mechanism(e0=E0, s0=S0, k1=k1, k_minus1=k_1, k2=k2, t_end=40.0)
    ax.plot(sim["t"], sim["P"], color=color, lw=2, label=f"[S]₀ = {S0} µM")
    ventana = (sim["t"] > 0.02) & (sim["t"] < 1.0)
    pend = cin.linear_fit(sim["t"][ventana], sim["P"][ventana])["slope"]
    ax.plot([0, 8], [0, 8 * pend], color=viz.INK_MUTED, ls="--", lw=1)
    ax.annotate(f"v₀ = {pend:.2f} µM/s", (8, 8 * pend), xytext=(6, 0), textcoords="offset points", fontsize=9, color=viz.INK_SECONDARY)
ax.set_xlabel("Tiempo (s)"); ax.set_ylabel("Producto (µM)"); ax.legend()
ax.set_title("Curva de progreso: la pendiente inicial es la velocidad v₀", loc="left")
fig;
print("Con más sustrato la pendiente inicial es mayor, pero no proporcionalmente. Ese 'no proporcionalmente' es toda la cinética enzimática.")
''')

# ============================================================================ 3. estructura
md(r"""
## 3. La estructura 3D del complejo catalítico

**La idea simple.** Para simular una enzima necesitamos saber dónde está cada átomo. Eso lo
dan los cristalógrafos: la estructura **3FGU** del Protein Data Bank es una "foto" de la
glucoquinasa humana con **glucosa**, un análogo de ATP (AMP‑PNP, que no puede reaccionar y por
eso se deja fotografiar) y el ion **Mg²⁺**, todos dentro del sitio activo, justo antes de la
reacción.

**Más detalle.** La glucoquinasa tiene dos dominios (grande y pequeño) que se cierran sobre la
glucosa como una almeja. En la forma cerrada el fósforo terminal del ATP (Pγ) queda a ~2.7 Å
del oxígeno O6 de la glucosa, que recibirá el fosfato. El carboxilato de **Asp205** está a
2.5 Å del hidroxilo O6–H: es la **base catalítica**, que se llevará el protón. **Lys169** y el
Mg²⁺ neutralizan las cargas negativas del fosfato. Gira la molécula con el ratón:
""")

code(r'''
pdb_cristal = pathlib.Path("data/raw/3FGU.pdb").read_text()
vista = viz.view_complex(pdb_cristal, ligand_resnames=("BGC", "ANP"), ion_resnames=("MG", "K"),
                         highlight_residues=(("ASP", 205), ("LYS", 169), ("THR", 228), ("SER", 151)),
                         label_map={205: "Asp205", 169: "Lys169", 228: "Thr228", 151: "Ser151"})
vista
''')

code(r'''
# Distancias clave en el cristal (Å): ¿está todo listo para reaccionar?
def coords_pdb(texto, resn, name, resi=None):
    for l in texto.splitlines():
        if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == resn and l[12:16].strip() == name and l[16] in " A":
            if resi is None or int(l[22:26]) == resi:
                return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])

PG, O6, OD1 = coords_pdb(pdb_cristal, "ANP", "PG"), coords_pdb(pdb_cristal, "BGC", "O6"), coords_pdb(pdb_cristal, "ASP", "OD1", 205)
print(f"Pγ ··· O6 (glucosa)  = {np.linalg.norm(PG - O6):.2f} Å   (un enlace P–O mide 1.6 Å: están muy cerca, pero aún no unidos)")
print(f"O6 ··· OD1 (Asp205)  = {np.linalg.norm(O6 - OD1):.2f} Å   (enlace de hidrógeno corto: Asp205 está listo para tomar el protón)")
''')

md(r"""
> **Para pensar.** El cristal usa AMP‑PNP, en el que el oxígeno entre los fósforos β y γ se
> cambió por un nitrógeno (N3B). Ese cambio impide la reacción y permite capturar el complejo.
> Para simular la reacción real, nosotros volveremos a poner el oxígeno: ATP verdadero.
""")

# ============================================================================ 4. preparación
md(r"""
## 4. Preparar la enzima para simularla

**La idea simple.** Una foto cristalográfica no está lista para simular: le faltan los
hidrógenos (los rayos X no los ven), le faltan trozos de cadena que estaban desordenados y no
tiene agua alrededor. Hay que completarla y sumergirla en una caja de agua con iones, como
está en la célula.

**Qué hizo el script `scripts/01_preparar_sistema.py`.**

1. Cadena A de 3FGU, residuos 5–458; tres bucles internos faltantes modelados con PDBFixer.
2. AMP‑PNP → ATP (N3B → O3B). Glucosa y ATP parametrizados con GAFF2 y cargas AM1‑BCC.
3. Proteína con el campo de fuerza Amber ff14SB; Mg²⁺ y K⁺ con parámetros de Li–Merz.
4. Aguas cristalográficas conservadas; caja de agua TIP3P de 10 Å y 21 Na⁺ para neutralizar.

En modo `completo` esa preparación requiere AmberTools (antechamber, tleap), que no está en
Colab; por eso los archivos preparados viajan con el repositorio.
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

# ============================================================================ 5. MD
md(r"""
## 5. La enzima se mueve: dinámica molecular

**La analogía.** La estructura cristalina es una foto; la **dinámica molecular (MD)** es la
película. Cada átomo es una bolita unida a sus vecinas por resortes (los enlaces) y que siente
las cargas y los choques de las demás. Con las leyes de Newton calculamos hacia dónde se mueve
cada una en pasos de 2 femtosegundos (0.000000000000002 s). Un nanosegundo de película son
500 000 pasos.

**Qué se hizo.** Con OpenMM y Amber: minimización, calentamiento de 0 a 300 K, equilibración
a 1 bar y 1 ns de producción (script `02_dinamica_molecular.py`). Aquí no hay química (los
enlaces no se rompen): la MD sirve para ver **cuánto tiempo pasa el sitio activo en una
geometría reactiva**, con el Pγ del ATP cerca del O6 de la glucosa. Esas geometrías se llaman
**conformaciones de ataque cercano** (*near‑attack conformations*, NAC).

**Cómo leer las gráficas.** El RMSD mide cuánto se ha alejado la proteína de la estructura
inicial (1 Å es poco: la enzima es estable). Las distancias muestran que el fosfato y la
glucosa "vibran" pero siguen apuntándose.
""")

code(r'''
if MODO == "completo":
    subprocess.run([sys.executable, "scripts/02_dinamica_molecular.py"], check=True)

md_df = datos.csv("md/analisis.csv")
md_info = datos.json_("md/md.json")
print(f"Producción: {md_info['produccion_ns']} ns, {md_info['n_atomos']:,} átomos, plataforma {md_info['plataforma']}")
fig = viz.plot_md_summary(md_df, crystal_value=2.68, subtitle="RMSD del esqueleto y distancias del sitio activo durante la producción")
fig;
''')

code(r'''
umbral = 3.5  # Å: definición práctica de conformación de ataque cercano
fraccion = (md_df["d_PG_O6_A"] < umbral).mean()
print(f"d(Pγ–O6) media = {md_df['d_PG_O6_A'].mean():.2f} Å; fracción de cuadros con d < {umbral} Å = {fraccion:.0%}")
print("Idea clave: la enzima mantiene los reactivos 'apuntándose' la mayor parte del tiempo. Sin enzima, en agua, ese encuentro sería un golpe de suerte.")
''')

# ============================================================================ 6. Eyring
md(r"""
## 6. Colinas de energía y velocidad de reacción

**La analogía del paso de montaña.** Una pelota en un valle quiere pasar a otro valle más bajo,
pero entre los dos hay una colina. Cuanto más alta, más raro es que la pelota la cruce. En
química la colina es la **barrera de activación** y su cima es el **estado de transición**.

**La ecuación (Eyring, 1935) y qué significa cada término:**

$$k \;=\; \kappa\,\frac{k_B T}{h}\;\exp\!\left(-\frac{\Delta G^{\ddagger}}{RT}\right)$$

| Término | Qué es | Valor típico |
|---|---|---|
| $k$ | constante de velocidad: cuántas veces por segundo cruza la barrera cada molécula lista para hacerlo | ≈ 60 s⁻¹ para la glucoquinasa |
| $k_B T/h$ | "frecuencia de intentos": las veces por segundo que la molécula se lanza contra la colina | 6.2 × 10¹² s⁻¹ a 25 °C |
| $\Delta G^{\ddagger}$ | altura de la colina (energía libre de activación) | 10–25 kcal/mol |
| $R T$ | energía térmica disponible | 0.59 kcal/mol a 25 °C |
| $\kappa$ | coeficiente de transmisión: fracción de cruces que no "rebotan" | ≈ 1 |

**Qué pasa si…**

* **ΔG‡ sube 1.36 kcal/mol** → la exponencial cae a la décima parte: la reacción es **10 veces
  más lenta**. Es la regla del pulgar más útil de toda la cinética.
* **T sube de 25 a 37 °C** → RT crece un 4 %: para una barrera de 15 kcal/mol la velocidad se
  duplica, aproximadamente.
* **La enzima baja ΔG‡ 10 kcal/mol** → la reacción va 10⁷ veces más rápido. Así de potente es
  bajar la colina.

Muévelo tú mismo:
""")

code(r'''
# 🎛️ Explora la ecuación de Eyring: barrera y temperatura
interactivo.explorar_eyring()
''')

code(r'''
# 🎛️ Perfil de energía: barrera (cuesta) y energía de reacción (desnivel entre valles)
interactivo.explorar_perfil_energia()
''')

code(r'''
barreras = np.linspace(5, 30, 200)
fig, ax = viz.figure(7, 4)
ax.semilogy(barreras, cin.eyring_rate(barreras), color=viz.COLORS["reactivo"], lw=2)
kcat = valor("kcat_s", 60.0)
for dg, texto in [(cin.barrier_from_rate(kcat), f"k_cat de la glucoquinasa (≈{kcat:.0f} s⁻¹)"),
                  (cin.barrier_from_rate(1e-8), "una reacción que tarda años")]:
    ax.axvline(dg, color=viz.INK_MUTED, ls="--", lw=1)
    ax.annotate(f"{texto}\nΔG‡ ≈ {dg:.1f} kcal/mol", (dg, cin.eyring_rate(dg)), xytext=(8, 10), textcoords="offset points",
                fontsize=9, color=viz.INK_SECONDARY)
ax.set_xlabel("ΔG‡ (kcal/mol)"); ax.set_ylabel("k (s⁻¹)")
ax.set_title("Cada 1.36 kcal/mol de barrera, la velocidad cambia 10 veces", loc="left")
fig;
''')

# ============================================================================ 7. QM/MM
md(r"""
## 7. Mirar la reacción con mecánica cuántica: QM/MM

**La analogía.** Para fotografiar un partido no necesitas enfocar a los 50 000 espectadores:
enfocas el balón y los jugadores cercanos, y el estadio queda como fondo. En la enzima, "el
balón" es el fosfato que salta; el estadio, el resto de la proteína y el agua.

**Por qué hace falta.** La dinámica molecular clásica no puede romper ni formar enlaces: sus
resortes no saben de electrones. Para ver cómo el fosfato salta del ATP a la glucosa hay que
usar **mecánica cuántica (QM)**, que sí describe los electrones. Pero es carísima: imposible
para 57 000 átomos.

**La solución: QM/MM.** Se trata con QM solo el pedacito donde ocurre la química (la **región
QM**, unas decenas de átomos) y el resto con el campo de fuerza clásico (**región MM**). Las
dos regiones "se ven":

* las cargas de la región MM crean un potencial eléctrico que entra en el cálculo cuántico
  (**embedding electrostático**), y
* la región MM empuja a la QM con fuerzas de van der Waals (**Lennard‑Jones**): sin esto la
  región QM se "hundiría" en la proteína, porque la electrostática sola atrae.

**Nuestra receta (herramientas de la suite Leonardo, MOPAC PM7 + Amber):**

* **Región QM** (78 átomos, carga −2): glucosa, el fragmento C5'–trifosfato del ATP, las cadenas
  laterales de Asp205, Lys169 y Thr228, el Mg²⁺ y su agua coordinada.
* Donde un enlace covalente cruza la frontera se pone un **átomo de enlace** (un H).
* El **entorno MM** (2300 cargas dentro de 16 Å) se mantiene **fijo** (aproximación de
  entorno rígido).
* Energía total: $E = E_{\mathrm{PM7}}[\text{QM en el potencial de la enzima}] + E_{\mathrm{LJ}}(\text{QM–MM})$.
""")

code(r'''
particion = datos.json_("qmmm/particion.json")
print(f"Región QM: {particion['n_qm']} átomos + {particion['n_link']} H de enlace; libres: {particion['n_free']}")
print(f"Entorno MM: {particion['n_mm']} cargas puntuales; carga QM = {particion['qm_charge']}; carga MM = {particion['mm_charge']:.2f} e")
print("Átomos QM por residuo:", pd.Series([a["resname"] for a in particion["qm_atoms"]]).value_counts().to_dict())
viz.view_qm_region(leer("qmmm/region_qm_inicial.pdb"), leer("qmmm/entorno_mm_8A.pdb"), link_atom_indices=particion["fixed"])
''')

md(r"""
**¿Cómo sabemos que el programa calcula bien las fuerzas?** Comparando el gradiente analítico
con el numérico (mover un átomo 0.005 Å y ver cuánto cambia la energía). Es una comprobación
que todo código de simulación debería pasar antes de usarse:
""")

code(r'''
grad = datos.csv("qmmm/comprobacion_gradiente.csv")
grad["error relativo"] = (grad["diff"].abs() / grad["analytic"].abs().clip(lower=1e-6)).map(lambda x: f"{x:.1%}")
display(grad.round(3))
''')

# ============================================================================ 8. R, P y escaneo
md(r"""
## 8. Reactivo, producto y el perfil de energía

**La idea simple.** Primero relajamos la región QM en su punto de partida (el **reactivo**:
glucosa + ATP). Luego construimos el **producto** (glucosa‑6‑fosfato + ADP, con el protón de
O6 ya en Asp205) y también lo relajamos. La diferencia de energía dice si la reacción es
"cuesta abajo" o "cuesta arriba" dentro de la enzima.

**La coordenada de reacción.** Para seguir la transferencia usamos un solo número que cambia
de forma continua entre reactivo y producto:

$$\xi = d(\mathrm{P_\gamma\!-\!O_{3\beta}}) - d(\mathrm{P_\gamma\!-\!O_6})$$

| Situación | d(Pγ–O3β) | d(Pγ–O6) | ξ |
|---|---|---|---|
| reactivo (fosfato en el ATP) | corta (1.6 Å) | larga (3.2 Å) | negativa |
| estado de transición | intermedia | intermedia | ≈ 0 |
| producto (fosfato en la glucosa) | larga | corta (1.7 Å) | positiva |

Un **escaneo relajado** fija ξ en valores intermedios (con una restricción armónica) y deja
que todo lo demás se acomode: el resultado es un **perfil de energía** aproximado, y su punto
más alto es una primera estimación del estado de transición.
""")

code(r'''
if MODO == "completo":
    subprocess.run([sys.executable, "scripts/03_qmmm_reaccion.py"], check=True)

res = datos.json_("qmmm/resumen.json")
R, P = res["etapas"]["reactivo"], res["etapas"]["producto"]
print("Reactivo:  d(Pγ–O6) = %.2f Å, d(Pγ–O3β) = %.2f Å, ξ = %+.2f Å" % (R["d_PG_O6"], R["d_PG_O3B"], R["xi"]))
print("Producto:  d(Pγ–O6) = %.2f Å, d(Pγ–O3β) = %.2f Å, ξ = %+.2f Å" % (P["d_PG_O6"], P["d_PG_O3B"], P["xi"]))
print("Protón del O6: en el reactivo d(O6–H) = %.2f Å; en el producto d(OD1–H) = %.2f Å (Asp205 protonado)" % (R["d_O6_H"], P["d_OD1_H"]))
print("ΔE(reacción) = %+.1f kcal/mol dentro de la enzima (entorno fijo, PM7)" % P["dE_reaccion_kcal"])
''')

code(r'''
escaneo = datos.csv("qmmm/escaneo.csv")
xi_esc = np.r_[R["xi"], escaneo["xi"].values]                   # el perfil arranca en el reactivo relajado (E = 0)
e_esc = np.r_[0.0, escaneo["energia_rel_kcal"].values]
fig = viz.plot_energy_profile(xi_esc, e_esc, relative=False, xlabel="ξ = d(Pγ–O3β) − d(Pγ–O6)  (Å)",
                              ts_index=int(np.argmax(e_esc)), title="Escaneo relajado de la transferencia de fosforilo",
                              subtitle="PM7 en el campo de la enzima; el máximo es una primera estimación del TS")
fig;
''')

md(r"""
> **Un detalle que llama la atención.** El perfil arranca con un salto de ~4 kcal/mol entre el
> reactivo (ξ = −1.56 Å, E = 0) y el primer punto del escaneo. La optimización inicial del
> reactivo se detuvo en un *hombro* de la superficie; el verdadero mínimo apareció después, al
> descender desde el estado de transición (sección 9). Lección práctica: un optimizador
> "convergido" no garantiza el mínimo más bajo; por eso los caminos de reacción se verifican en
> las dos direcciones.
""")

code(r'''
# Animación del escaneo: el fosfato viaja del ATP a la glucosa
simbolos, cuadros, comentarios = datos.leer_xyz_multiple("qmmm/escaneo.xyz")
viz.view_frames(cuadros, simbolos, interval_ms=300)
''')

# ============================================================================ 9. TS
md(r"""
## 9. El estado de transición

**La analogía del paso de montaña, otra vez.** El escaneo nos dio una colina aproximada, pero
forzando una sola distancia. El verdadero estado de transición es un **punto de silla**: un
máximo a lo largo del camino y un mínimo en todas las demás direcciones, exactamente como el
paso entre dos montañas (si te sales del camino, subes).

**Cómo se localiza (tres pasos):**

1. **NEB con imagen trepadora** (*nudged elastic band*): una cadena de geometrías entre reactivo
   y producto, unidas por muelles, que se relaja hasta dibujar el camino de mínima energía. La
   imagen más alta "trepa" hasta la cima. Es el análogo moderno del método **QST2** (dos
   extremos → TS) de Gaussian y MOPAC.
2. **Método del dímero**: refina la cima usando solo gradientes, hasta que la fuerza es cero.
3. **Frecuencias**: en un punto de silla verdadero hay **exactamente una** frecuencia imaginaria,
   y su vector describe el movimiento de la reacción. Es la prueba definitiva.
""")

code(r'''
neb = datos.csv("qmmm/neb.csv")
perfiles = [(xi_esc, e_esc, "escaneo restringido"),
            (np.r_[R["xi"], neb["xi"].values[1:]], np.r_[0.0, neb["energia_rel_kcal"].values[1:]], "NEB (imagen trepadora)")]
fig = viz.plot_energy_profiles(perfiles, xlabel="ξ (Å)", relative=False, title="Del escaneo al camino de mínima energía",
                               subtitle="El NEB relaja todas las coordenadas a la vez; la imagen trepadora sube a la cima")
fig;
''')

code(r'''
ts = res["etapas"].get("dimero")
if ts is None or not ts.get("convergido", True):
    print("(el dímero no convergió: se usa la imagen trepadora del NEB como TS)")
    ts = ts or {**res["etapas"]["neb"], "xi": neb.loc[res["etapas"]["neb"]["imagen_ts"], "xi"], "d_PG_O6": float("nan"), "d_PG_O3B": float("nan")}
fr = res["etapas"].get("frecuencias", {})
print("TS refinado (dímero): barrera ΔE‡ = %.1f kcal/mol; ξ = %+.2f Å; d(Pγ–O6) = %.2f Å; d(Pγ–O3β) = %.2f Å"
      % (ts["barrera_kcal"], ts["xi"], ts["d_PG_O6"], ts["d_PG_O3B"]))
print("Frecuencias más bajas (cm⁻¹, negativas = imaginarias):", np.round(fr["frecuencias_mas_bajas"], 1))
val = fr["validacion"]
print("Modos imaginarios:", val["imaginary_mode_count"], "→", "punto de silla de primer orden ✔" if val["ok"] else "revisar ✘")
''')

md(r"""
**Qué nos dice la geometría del estado de transición.** En el TS el fósforo está a la misma
distancia de los dos oxígenos (≈ 2.1–2.2 Å de O6 y de O3β): es un TS **concertado y "en
línea"**; el fosfato pasa de un oxígeno al otro invirtiendo sus tres oxígenos como un paraguas.
El protón del O6, en cambio, **todavía no se ha movido** (sigue a ~1.07 Å de O6). Solo después
de cruzar la cima, cuando el enlace P–O6 ya está formado, el protón salta a Asp205: el camino
baja por una "meseta" hasta ξ ≈ 1 Å y allí termina de caer. En este modelo Asp205 actúa de base
**tras** la transferencia del fosforilo. Un estudio QM/MM independiente (Zhang et al., 2009)
describe a Asp205 como base general y a Lys169 como ácido general que protona el fosfato; que
dos modelos discrepen en el *orden* de los pasos, pero coincidan en la barrera (sección 12), es
un buen ejemplo de por qué los mecanismos se estudian con simulaciones y con mutantes.

**El movimiento del estado de transición.** La frecuencia imaginaria corresponde a una vibración
que no oscila sino que "cae" hacia reactivo o hacia producto. Animándola vemos la química:
""")

code(r'''
vibra = np.load(datos.ruta("qmmm/frecuencias_ts.npz"))
cuadros_modo = viz.mode_animation_frames(vibra["ts_xyz"], vibra["modo_imaginario"], n_frames=24, amplitude=0.6)
viz.view_frames(cuadros_modo, list(vibra["simbolos"]), interval_ms=60, bonds_from=0)
''')

md(r"""
**Bajar de la cima.** Si empujamos el TS un poquito a cada lado a lo largo del modo imaginario
y dejamos que la geometría relaje, descendemos por el camino de mínima energía hasta el
reactivo y hasta el producto. Esto confirma que **este** TS conecta **estos** reactivos con
**estos** productos (la idea del IRC, *intrinsic reaction coordinate*).
""")

code(r'''
camino = datos.csv("qmmm/camino_descenso.csv")
fig = viz.plot_energy_profile(camino["xi"].values, camino["energia_rel_kcal"].values, relative=False, xlabel="ξ (Å)",
                              ts_index=int(camino["energia_rel_kcal"].idxmax()), smooth=False,
                              title="Camino de mínima energía descendente desde el TS",
                              subtitle="Cada punto es una geometría relajada; el TS está en el centro")
fig;
''')

code(r'''
niveles = [("E·S (reactivo)", 0.0), ("TS", ts["barrera_kcal"]), ("E·P (producto)", P["dE_reaccion_kcal"])]
fig = viz.plot_energy_levels(niveles, ts_indices=[1], title="Diagrama de energía de la reacción en la enzima",
                             subtitle="PM7/Amber, entorno fijo: energías electrónicas relativas (ΔE), todavía no energías libres")
fig;
''')

# ============================================================================ 10. mejoras
md(r"""
## 10. De ΔE‡ a ΔG‡: tres mejoras del cálculo

Hasta aquí tenemos una **energía electrónica** ΔE‡, calculada en **una** conformación de la
enzima con **un** método semiempírico. La ecuación de Eyring pide una **energía libre** ΔG‡. Un
buen científico computacional se hace tres preguntas, y las responde con cálculos adicionales:

**1. ¿Y las vibraciones?** Los átomos nunca están quietos: incluso a 0 K tienen la **energía
de punto cero** (ZPE), y a 300 K, energía y entropía vibracionales. Con las frecuencias del
reactivo y del TS (aproximación armónica) se obtiene

$$\Delta G^{\ddagger} = \Delta E^{\ddagger} + \Delta \mathrm{ZPE} + \Delta H_{\mathrm{vib}}(T) - T\,\Delta S_{\mathrm{vib}}(T).$$

* **ΔZPE**: el enlace P–O que se rompe pierde su vibración de tensión en el TS, así que suele
  ser **negativa** (baja la barrera).
* **−TΔS‡**: si el TS es más "rígido" que el reactivo (ΔS‡ < 0), este término **sube** la barrera.

**2. ¿Y si la enzima estuviera en otra postura?** Como el entorno está fijo, cada instantánea
de la MD da una barrera distinta. Repetimos el cálculo en varias instantáneas (250, 500, 750 y
1000 ps) y miramos la **media y la dispersión**. Un buen resultado no es un número: es un
número con su incertidumbre.

**3. ¿Y si el método se equivoca?** PM7 es semiempírico. Recalculamos las energías con otro
hamiltoniano (PM6‑D3H4) en las mismas geometrías: si la barrera cambia poco, el resultado es
robusto; si cambia mucho, sabemos cuánto (des)confiar.
""")

code(r'''
termo = res["etapas"].get("termoquimica", {})
T_ref = "298.15"
if T_ref in termo:
    t = termo[T_ref]
    print("Termoquímica armónica a 298.15 K (átomos libres de la región QM):")
    print(f"  ΔE‡ = {t['dE_kcal']:.1f}   ΔZPE = {t['dZPE_kcal']:+.1f}   ΔH‡ = {t['dH_kcal']:.1f} kcal/mol   ΔS‡(vib) = {t['dS_vib_cal']:+.1f} cal/mol/K   →   ΔG‡ = {t['dG_kcal']:.1f} kcal/mol")
    for Tk in ("303.15", "310.15"):
        if Tk in termo:
            print(f"  a {float(Tk) - 273.15:.0f} °C: ΔG‡ = {termo[Tk]['dG_kcal']:.1f} kcal/mol")
    print("  ", termo.get("nota", ""))
else:
    print("(termoquímica no disponible en estos datos)")
''')

code(r'''
try:
    inst = datos.json_("qmmm/instantaneas.json")
    filas = [(r["instantanea"], r.get("barrera_kcal", np.nan), r.get("ts_d_PG_O6", np.nan), r.get("ts_d_PG_O3B", np.nan),
              r.get("ts_origen", "dímero" if "barrera_kcal" in r else "falló"))
             for r in inst["instantaneas"]]
    df_inst = pd.DataFrame(filas, columns=["instantánea", "ΔE‡ (kcal/mol)", "d(Pγ–O6) TS (Å)", "d(Pγ–O3β) TS (Å)", "TS obtenido por"])
    display(df_inst.round(2))
    ok = df_inst.dropna(subset=["ΔE‡ (kcal/mol)"])
    media, sd = ok["ΔE‡ (kcal/mol)"].mean(), ok["ΔE‡ (kcal/mol)"].std(ddof=1) if len(ok) > 1 else 0.0
    if len(ok) > 1:
        print(f"Barrera media = {media:.1f} ± {sd:.1f} kcal/mol (n = {len(ok)})")
        fig, ax = viz.figure(6.5, 3.8)
        ax.bar(range(len(ok)), ok["ΔE‡ (kcal/mol)"], color=viz.COLORS["ts"], width=0.55)
        ax.axhline(media, color=viz.INK, lw=1, ls="--")
        ax.set_xticks(range(len(ok))); ax.set_xticklabels(ok["instantánea"], rotation=15)
        ax.set_ylabel("ΔE‡ (kcal/mol)"); ax.set_title("La barrera depende de la conformación de la enzima", loc="left")
        fig;
    else:
        print(f"Solo {len(ok)} conformación con camino completo (barrera {media:.1f} kcal/mol): las demás no llegan a un producto estable (ver abajo).")
    sin_prod = [r for r in inst["instantaneas"] if r.get("sin_producto_estable")]
    if sin_prod:
        print("Instantáneas sin producto estable (la energía sube hasta ξ = 2 Å sin máximo):")
        for r in sin_prod:
            print(f"  {r['instantanea']}: E(ξ = {r['xi_max']:.1f}) = {r['escaneo_max_rel_kcal']:.0f} kcal/mol sobre R; d(Asp205 OD1···H–O6) en el reactivo = {r['d_OD1_H_reactivo']:.2f} Å")
except FileNotFoundError:
    print("(promedio sobre instantáneas no disponible en estos datos)")
''')

md(r"""
**Lo que enseñan las instantáneas de la MD.** En la conformación cristalina minimizada, el
hidroxilo O6–H de la glucosa está unido por puente de hidrógeno a Asp205 (1.7 Å) y el camino
tiene un producto estable. En varias instantáneas de la dinámica clásica, en cambio, Asp205 se
ha alejado del O6–H (en la MD la distancia O6···OD1 promedia 4.3 Å, sección 5): al forzar la
transferencia del fosforilo el protón **no encuentra a la base** y la energía sube sin parar,
sin mínimo de producto. Dos lecciones: (1) la reacción solo puede ocurrir desde las
**conformaciones de ataque cercano** correctas, que son una fracción del tiempo; (2) el campo de
fuerza clásico no siempre mantiene la geometría catalítica (aquí, la interacción
hidroxilo–carboxilato junto al Mg²⁺), y por eso el modelador debe comparar la MD con el cristal
antes de elegir la estructura de partida del QM/MM. Una barrera "promedio" solo tiene sentido
sobre conformaciones reactivas.
""")

code(r'''
met = res["etapas"].get("metodos", {})
for nombre, d in met.items():
    if isinstance(d, dict):
        print(f"{nombre:>9s}//PM7: barrera {d['barrera_kcal']:.1f} kcal/mol, ΔE(reacción) {d['dE_reaccion_kcal']:+.1f} kcal/mol")
print(met.get("nota", ""))
''')

md(r"""
**Lo que dicen estos números.** La corrección vibracional es pequeña (ΔZPE ≈ −0.2 y −TΔS‡ ≈ +1.3
kcal/mol se compensan en parte): ΔG‡ ≈ ΔE‡ + 0.6 kcal/mol. El promedio sobre instantáneas dice
cuánto "respira" la barrera con la conformación de la enzima. Y el cambio de hamiltoniano es la
sorpresa: PM6‑D3H4 da energías muy distintas en las mismas geometrías. Los métodos semiempíricos
se parametrizan con conjuntos de moléculas concretos y el trío fosfato + Mg²⁺ + carboxilato es
un caso difícil para ellos; PM7 se ajustó con muchos más compuestos de fósforo que PM6, y por eso
lo usamos. La lección: **antes de creer una barrera, hay que saber para qué química fue
calibrado el método** y, si es posible, comparar con un cálculo de más nivel (DFT) o con el
experimento (sección 12).
""")

# ============================================================================ 11. agua
md(r"""
## 11. La misma reacción **fuera** de la enzima

**La idea simple.** ¿Cuánto de la catálisis lo hace el *resto* de la enzima? Sacamos el sitio
activo (el mismo clúster QM, con sus anclajes fijos) y lo ponemos en **agua** (disolvente
implícito COSMO), sin el campo eléctrico de los otros 2300 átomos. Si la barrera cambia, ese
cambio es obra del entorno proteico.

**Los métodos nativos de MOPAC.** Fuera de la enzima no hay potencial externo que "congelar",
así que aquí sí usamos los buscadores de estado de transición nativos de MOPAC, con los nombres
con que los orquesta la suite Leonardo:

| Nombre en Leonardo | Palabra clave MOPAC | Qué hace |
|---|---|---|
| QST2 | `SADDLE` | Busca el TS a partir de reactivo y producto (dos extremos) |
| QST3 / TS | `TS` | Refina el punto de silla siguiendo el modo de curvatura negativa |
| Validación | `FORCETS` | Frecuencias sobre las coordenadas libres: debe haber una imaginaria |
| Camino | `IRC=1*` | Sigue el camino intrínseco de reacción en ambas direcciones |

Un detalle práctico: `SADDLE` da una *estimación* del TS que hay que refinar. Aquí el
refinamiento se hizo con el método del dímero y después con la palabra clave `TS` de MOPAC,
que llegaron al mismo punto de silla.
""")

code(r'''
agua = res["etapas"]["agua"]
print("1) SADDLE (QST2): estimación del TS %.1f kcal/mol por encima del reactivo (ξ = %+.2f Å)" % (agua["barrera_saddle_qst2_kcal"], agua["saddle_xi"]))
print("2) Refinamiento (%s): barrera = %.1f kcal/mol; ΔE(reacción) = %+.1f kcal/mol" % (agua["metodo_ts"], agua["barrera_kcal"], agua["dE_reaccion"]))
print("   geometría del TS en agua: d(Pγ–O6) = %.2f Å, d(Pγ–O3β) = %.2f Å, d(O6–H) = %.2f Å" % (agua["ts_d_PG_O6"], agua["ts_d_PG_O3B"], agua["ts_d_O6_H"]))
print("3) FORCETS de MOPAC, frecuencias más bajas (cm⁻¹):", np.round(agua["frecuencias_mas_bajas"], 1))
if "frecuencias_ase_mas_bajas" in agua:
    print("   Hessiano numérico (ASE), frecuencias más bajas (cm⁻¹):", np.round(agua["frecuencias_ase_mas_bajas"], 1),
          "→ modos imaginarios:", agua["validacion_ts_ase"]["imaginary_mode_count"])
''')

md(r"""
**Una lección sobre validar.** El Hessiano de `FORCETS` sobre la superficie COSMO muestra varias
frecuencias "imaginarias" pequeñas (|ν| < 40 cm⁻¹): son ruido numérico de la cavidad del
disolvente, no movimientos reales; por eso conviene mirar la *magnitud* de las frecuencias y no
solo contarlas. El Hessiano por diferencias finitas más finas (ASE) deja un solo modo claro
(≈ −160 cm⁻¹, el fósforo que salta) y un residuo de −37 cm⁻¹ atribuible a la cavidad. El `IRC`
de MOPAC abortó por un fallo interno del programa en esta superficie con disolvente (también
eso pasa en la práctica); el camino se obtuvo con el mismo descenso desde el TS de la sección 9.
""")

code(r'''
try:
    cam_agua = datos.csv("qmmm/agua_camino_descenso.csv")
    fig = viz.plot_energy_profile(cam_agua["xi"].values, cam_agua["energia_kcal"].values - agua["E_reactivo"], relative=False,
                                  xlabel="ξ (Å)", smooth=False, ts_index=int(np.argmax(cam_agua["energia_kcal"].values)),
                                  title="Camino de reacción del sitio activo en agua (COSMO)",
                                  subtitle="Energías relativas al reactivo en agua; el máximo es el TS refinado")
    fig;
except FileNotFoundError:
    print("(camino en agua no disponible)")
''')

code(r'''
con = [("E·S", 0.0), ("TS", ts["barrera_kcal"]), ("E·P", P["dE_reaccion_kcal"])]
sin = [("R", 0.0), ("TS", agua["barrera_kcal"]), ("P", agua["dE_reaccion"])]
fig = viz.plot_energy_levels(con, compare=sin, ts_indices=[1], label="dentro de la enzima", compare_label="sitio activo en agua",
                             title="Dentro de la enzima frente a fuera de ella",
                             subtitle="La diferencia entre las dos colinas es el efecto del resto de la proteína")
fig;
print("Efecto del entorno proteico sobre la barrera: %+.1f kcal/mol" % (ts["barrera_kcal"] - agua["barrera_kcal"]))
''')

md(r"""
**Cómo leer esta comparación.** El clúster "en agua" **no** es la reacción sin catalizar:
todavía contiene la base catalítica (Asp205), la carga positiva de Lys169 y el Mg²⁺, es decir,
la maquinaria química esencial. Por eso las dos colinas se parecen: la diferencia (unas pocas
kcal/mol) mide lo que aporta el *resto* de la proteína, sobre todo su campo eléctrico y el
hecho de mantener los reactivos alineados. La reacción de verdad sin enzima, glucosa y ATP
solos en agua, es extraordinariamente lenta: la hidrólisis espontánea de un fosfato dianión
tiene una vida media de ~10¹² años (Lad, Williams y Wolfenden, 2003; ΔG‡ ≈ 44 kcal/mol),
porque ninguna base, ningún catión ni ningún campo eléctrico están ahí para ayudar. Entre esa
reacción y k<sub>cat</sub> ≈ 60 s⁻¹ hay ~10²¹ veces: ese es el tamaño real de la catálisis.
""")

# ============================================================================ 12. barrera -> kcat
md(r"""
## 12. De la barrera a *k*<sub>cat</sub>

**La idea simple.** Con la altura de la colina y la ecuación de Eyring estimamos cuántas veces
por segundo la enzima completa la reacción: ese número es **k**<sub>cat</sub>, el "número de
recambio". Lo comparamos con el valor medido en el laboratorio y con otro cálculo QM/MM
publicado.

**Con honestidad.** Nuestra ΔG‡ viene de un método semiempírico (PM7), de un entorno congelado y
de una aproximación armónica. Un error de 1.4 kcal/mol ya es un factor 10 en velocidad. Por eso
comparamos órdenes de magnitud y, sobre todo, *tendencias*, no decimales.
""")

code(r'''
T = 298.15
kcat_exp = valor("kcat_s", 60.0)
dG_qmmm = termo.get(T_ref, {}).get("dG_kcal", ts["barrera_kcal"])
filas = [("ΔE‡ PM7, una conformación (este cuaderno)", ts["barrera_kcal"]),
         ("ΔG‡ PM7 + termoquímica armónica (este cuaderno)", dG_qmmm)]
try:
    if len(ok) > 1:
        filas.append(("ΔE‡ promedio sobre instantáneas de MD (este cuaderno)", media))
except NameError:
    pass
filas.append(("ΔE‡ QM/MM publicada, Zhang et al. 2009", valor("qmmm_barrier_literature_kcal", np.nan)))
filas.append(("sitio activo en agua, COSMO (este cuaderno)", agua["barrera_kcal"]))
filas.append(("experimento: ΔG‡ que implica k_cat", cin.barrier_from_rate(kcat_exp, T)))
tabla = pd.DataFrame([(n, b, cin.eyring_rate(b, T)) for n, b in filas], columns=["estimación", "barrera (kcal/mol)", "k (s⁻¹) por Eyring"])
tabla["k (s⁻¹) por Eyring"] = tabla["k (s⁻¹) por Eyring"].map(lambda x: f"{x:.2e}")
display(tabla.round(1))
print(f"k_cat experimental ≈ {kcat_exp:.0f} s⁻¹  ({ref['kcat_s']['source']})")
print("Nuestra barrera y la de Zhang et al. (18.3 kcal/mol, otro programa y otra partición QM/MM) coinciden dentro de 1 kcal/mol;")
print("ambas sobrestiman la experimental (≈15) en 3-4 kcal/mol, es decir, un factor ~10²-10³ en k. Es lo esperable de un método semiempírico con entorno fijo.")
''')

# ============================================================================ 13. MM
md(r"""
## 13. Michaelis–Menten desde el mecanismo

**La analogía del cajero.** Una caja de supermercado (la enzima) atiende clientes (el sustrato).
Con pocos clientes, cada uno que llega es atendido enseguida: la "velocidad" de la caja crece
con el número de clientes. Con muchos clientes se forma cola: la caja trabaja a tope y da igual
que lleguen más. La velocidad se **satura**. Ese comportamiento tiene forma de hipérbola y se
describe con dos números:

* **V**<sub>max</sub>: la velocidad con la caja siempre ocupada (toda la enzima en forma ES).
* **K**<sub>M</sub>: la concentración de clientes a la que la caja trabaja a la **mitad** de su
  máximo. Una K<sub>M</sub> pequeña significa que la enzima "se llena" con poco sustrato.

**El mecanismo que hay detrás (Michaelis y Menten 1913; Briggs y Haldane 1925):**

$$\mathrm{E + S \;\underset{k_{-1}}{\overset{k_1}{\rightleftharpoons}}\; ES \;\overset{k_2}{\longrightarrow}\; E + P}$$

| Constante | Qué describe | Unidades |
|---|---|---|
| $k_1$ | el sustrato entra al sitio activo (choque productivo) | M⁻¹ s⁻¹ |
| $k_{-1}$ | el sustrato se suelta sin reaccionar | s⁻¹ |
| $k_2$ (= k<sub>cat</sub>) | el complejo ES cruza la barrera y suelta el producto | s⁻¹ |

**La derivación, paso a paso.** (Sigue cada línea: solo es álgebra.)

1. Velocidad de formación de producto: $v_0 = k_2[\mathrm{ES}]$.
2. **Hipótesis del estado estacionario**: ES se forma tan rápido como se consume, así que su
   concentración casi no cambia: $k_1[\mathrm{E}][\mathrm{S}] = (k_{-1} + k_2)[\mathrm{ES}]$.
3. La enzima total se reparte: $[\mathrm{E}]_0 = [\mathrm{E}] + [\mathrm{ES}]$, luego $[\mathrm{E}] = [\mathrm{E}]_0 - [\mathrm{ES}]$.
4. Sustituyendo (3) en (2) y despejando: $[\mathrm{ES}] = \dfrac{[\mathrm{E}]_0[\mathrm{S}]}{K_M + [\mathrm{S}]}$ con
   $K_M \equiv \dfrac{k_{-1}+k_2}{k_1}$.
5. Metiendo (4) en (1): $\boxed{v_0 = \dfrac{V_{\max}[\mathrm{S}]}{K_M + [\mathrm{S}]}}$ con $V_{\max} = k_2[\mathrm{E}]_0$.

**Qué significa cada término y qué pasa si cambia.**

* $[\mathrm{S}] \ll K_M$: $v_0 \approx (V_{\max}/K_M)[\mathrm{S}]$, una recta: la enzima está casi vacía y
  cada molécula de sustrato cuenta. Pendiente = $k_{\mathrm{cat}}/K_M \cdot [\mathrm{E}]_0$.
* $[\mathrm{S}] = K_M$: $v_0 = V_{\max}/2$.
* $[\mathrm{S}] \gg K_M$: $v_0 \to V_{\max}$: saturación.
* **Si K<sub>M</sub> sube** (la enzima "agarra" peor el sustrato o lo suelta más rápido), hace falta
  más sustrato para la misma velocidad: la curva se estira a la derecha.
* **Si V<sub>max</sub> sube** (más enzima, o una enzima más rápida), toda la curva escala hacia arriba.
* $K_M$ solo es igual a la constante de disociación $K_d = k_{-1}/k_1$ cuando $k_2 \ll k_{-1}$.

En lugar de creer la fórmula, **resolvamos el mecanismo numéricamente** y veamos cómo la
hipérbola aparece sola. Primero la película de las concentraciones:
""")

code(r'''
# 🎛️ Mueve k1, k-1, k2, [E]0 y [S]0 y observa el pre-estado estacionario y la formación de producto
interactivo.explorar_mecanismo()
''')

code(r'''
E0, k1, k_1, k2 = 0.05, 1.0, 50.0, 60.0          # µM, µM⁻¹s⁻¹, s⁻¹, s⁻¹ (k2 = k_cat ≈ 60 s⁻¹)
S = np.array([5, 10, 20, 40, 80, 150, 300, 600, 1200, 2400])       # µM
v0_sim = cin.initial_rates_from_simulation(E0, S, k1, k_1, k2, t_window=(0.01, 0.1))
ajuste = cin.fit_michaelis_menten(S, v0_sim)
teoria = cin.steady_state_parameters(k1, k_1, k2)
fig = viz.plot_initial_rates_from_ode(S, v0_sim, fit=ajuste)
fig;
print(f"K_M ajustada a las velocidades simuladas = {ajuste['km']:.1f} µM   frente a   (k₋₁ + k₂)/k₁ = {teoria['km']:.1f} µM")
print(f"V_max ajustada = {ajuste['vmax']:.3f} µM/s   frente a   k₂·[E]₀ = {k2 * E0:.3f} µM/s   →   la derivación funciona.")
''')

code(r'''
# 🎛️ La hipérbola: mueve Vmax y Km
interactivo.explorar_michaelis_menten()
''')

md(r"""
**Tres números que salen de la curva, y cómo derivarlos de los datos.**

| Constante | Cómo se obtiene | Qué mide | Glucoquinasa (literatura) |
|---|---|---|---|
| **K<sub>M</sub>** (o S<sub>0.5</sub>) | ajuste de v₀ frente a [S] | afinidad aparente; a qué [S] la enzima está medio llena | ≈ 7.7 mM de glucosa |
| **k<sub>cat</sub>** = V<sub>max</sub>/[E]<sub>0</sub> | V<sub>max</sub> del ajuste ÷ concentración de enzima | recambio: reacciones por segundo por molécula de enzima | ≈ 62–66 s⁻¹ |
| **k<sub>cat</sub>/K<sub>M</sub>** | cociente de los dos anteriores | eficiencia a baja [S]; tope físico ≈ 10⁸–10⁹ M⁻¹s⁻¹ (difusión) | ≈ 8 × 10³ M⁻¹s⁻¹ |

Ahora un "experimento": datos simulados con los parámetros de la literatura y un 4 % de ruido,
tratados como lo haría un bioquímico: **ajuste no lineal** y, para comparar, las tres
linealizaciones clásicas. Verás que **Lineweaver–Burk** amplifica el error de los puntos a baja
concentración: hoy se usa para *visualizar*, y el ajuste no lineal para *cuantificar*.
""")

code(r'''
rng = np.random.default_rng(7)
S_mM = np.array([0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 40, 60])
E0_uM = 0.2                                     # enzima en el tubo (µM)
Km_real = valor("s_half_mm", 7.5)               # tratamos la enzima como hiperbólica en este ejercicio
Vmax_real = valor("kcat_s", 60.0) * E0_uM       # µM/s
v_obs = cin.michaelis_menten(S_mM, Vmax_real, Km_real) * (1 + 0.04 * rng.standard_normal(S_mM.size))
aj = cin.fit_michaelis_menten(S_mM, v_obs)
fig = viz.plot_michaelis_menten(S_mM, v_obs, fit=aj, title="Curva de saturación de Michaelis–Menten",
                                subtitle="Datos simulados con ruido del 4 %; línea: ajuste no lineal")
fig;
print(f"V_max = {aj['vmax']:.2f} ± {aj['vmax_err']:.2f} µM/s;  K_M = {aj['km']:.2f} ± {aj['km_err']:.2f} mM;  R² = {aj['r2']:.4f}")
kc = cin.kcat_km_from_fit(aj, E0_uM)
print(f"k_cat = V_max/[E]₀ = {kc['kcat']:.0f} s⁻¹;  k_cat/K_M = {kc['kcat_over_km_molar']:.2e} M⁻¹s⁻¹;  tiempo por recambio = {cin.turnover_time(kc['kcat'])*1000:.0f} ms")
''')

code(r'''
fig = viz.plot_linearizations(S_mM, v_obs, subtitle="Las tres rectas dan V_max y K_M por los cortes con los ejes")
fig;
lb = cin.linear_fit(*cin.lineweaver_burk(S_mM, v_obs))
print(f"Lineweaver–Burk: V_max = {1/lb['intercept']:.2f} µM/s, K_M = {lb['slope']/lb['intercept']:.2f} mM   (ajuste no lineal: {aj['vmax']:.2f}, {aj['km']:.2f})")
''')

# ============================================================================ 14. Hill
md(r"""
## 14. Cooperatividad: la glucoquinasa es un sensor

**La analogía del interruptor.** Un regulador de luz (hipérbola) sube la intensidad poco a poco
desde el principio. Un interruptor con umbral (sigmoide) casi no hace nada hasta cierto punto y
entonces se enciende de golpe. El páncreas necesita un interruptor: por debajo de ~5 mM de
glucosa no debe soltar insulina, y por encima sí. La glucoquinasa es ese interruptor.

**La ecuación de Hill y sus términos:**

$$v_0 = \frac{V_{\max}[\mathrm{S}]^{n}}{S_{0.5}^{\,n} + [\mathrm{S}]^{n}}$$

* $S_{0.5}$: concentración a la que se alcanza la mitad de V<sub>max</sub> (el análogo de K<sub>M</sub>).
* $n$ = **coeficiente de Hill**: cuán abrupto es el interruptor. $n = 1$ es Michaelis–Menten;
  $n > 1$ es cooperatividad positiva. Para la glucoquinasa $n ≈ 1.7$.

**Qué pasa si…** con $n = 1$ hace falta multiplicar [S] por **81** para pasar del 10 % al 90 % de
la actividad; con $n = 1.7$ basta multiplicarla por ~13; con $n = 4$ (hemoglobina) por 3.

**Lo sorprendente.** La glucoquinasa es un **monómero** con un solo sitio para glucosa: no puede
haber "comunicación entre subunidades" como en la hemoglobina. Su cooperatividad es
**cinética**: la enzima cambia entre una forma poco activa (abierta) y una activa (cerrada) a
una velocidad **comparable a la de la catálisis** (k<sub>ex</sub> ≈ 5–100 s⁻¹ frente a
k<sub>cat</sub> ≈ 60 s⁻¹, Larion et al. 2012). Con mucha glucosa la enzima no tiene tiempo de
relajarse a la forma lenta entre un ciclo y el siguiente, y trabaja más de lo que "debería". Es
el modelo **mnemónico** o de **transición lenta** (Storer y Cornish‑Bowden 1977; Neet y Ainslie;
Cárdenas 1984). Un detalle elegante: un inhibidor competitivo como la N‑acetilglucosamina
*suprime* la cooperatividad (n → 1), porque mantiene ocupado el sitio y rompe la "memoria". Las
simulaciones de MD de la sección 5 (dominios que se abren y cierran) son la imagen molecular de
ese cambio.
""")

code(r'''
# 🎛️ Mueve n y S0.5: compara con la hipérbola y mira la ventana 10 %–90 %
interactivo.explorar_hill()
''')

code(r'''
S_g = np.linspace(0.01, 30, 300)
s_half, n_h = valor("s_half_mm", 7.5), valor("hill_n", 1.7)
fig = viz.plot_hill_vs_mm(S_g, cin.michaelis_menten(S_g, 1.0, s_half), cin.hill(S_g, 1.0, s_half, n_h), n_hill=n_h, s_half=s_half,
                          ylabel="v₀ / V_max", title="Glucoquinasa: sigmoide (Hill) frente a hipérbola",
                          subtitle="Entre 4 y 10 mM de glucosa (el rango fisiológico) la sigmoide es mucho más sensible")
ax = fig.axes[0]; ax.axvspan(4, 7, color=viz.GRID, alpha=0.6, zorder=0)
ax.annotate("glucosa en sangre\nen ayunas (4–7 mM)", (5.5, 0.05), ha="center", fontsize=9, color=viz.INK_SECONDARY)
fig;
s10, s90 = cin.substrate_at_fraction(1.0, s_half, n_h, 0.1), cin.substrate_at_fraction(1.0, s_half, n_h, 0.9)
print(f"Con n = {n_h}: del 10 % al 90 % de V_max entre {s10:.1f} y {s90:.1f} mM (factor {s90/s10:.0f}); con n = 1 haría falta un factor 81.")
''')

code(r'''
# Ajuste de Hill a datos simulados y gráfico de Hill (la pendiente es n)
S_h = np.array([1, 2, 3, 4, 5, 6, 7.5, 9, 11, 14, 18, 25, 35, 50])
v_h = cin.hill(S_h, 10.0, s_half, n_h) * (1 + 0.03 * rng.standard_normal(S_h.size))
ajh = cin.fit_hill(S_h, v_h)
print(f"Ajuste de Hill: V_max = {ajh['vmax']:.2f}, S_0.5 = {ajh['s_half']:.2f} mM, n = {ajh['n']:.2f} ± {ajh['n_err']:.2f}")
fig = viz.plot_hill_plot(S_h, v_h, ajh["vmax"], subtitle="log[v/(V_max − v)] frente a log[S]: la pendiente es el coeficiente de Hill")
fig;
''')

md(r"""
**Mutaciones que mueven el interruptor: GCK‑MODY e hiperinsulinismo.** Una mutación que sube
S<sub>0.5</sub> (o baja k<sub>cat</sub>) hace que el páncreas "vea" menos glucosa de la que hay y
libere insulina tarde: glucemia alta desde el nacimiento, pero estable (GCK‑MODY). Una mutación
activadora baja S<sub>0.5</sub>: insulina de más, hipoglucemia congénita. Con los valores
publicados (Valentínová 2012; Sayed 2009) dibujamos qué hace cada mutante a 5 mM de glucosa:
""")

code(r'''
mut = ref.get("mutants", [])
if mut:
    S_g = np.linspace(0.01, 30, 300)
    elegidos = [m for m in mut if m["name"] in ("V244G", "G223S", "I110N", "W99L", "M197I")]
    curvas = [("silvestre", s_half, n_h, 1.0)] + [(m["name"], m["s_half_mm"], m.get("hill_n", n_h), m.get("kcat_rel", 1.0)) for m in elegidos]
    fig, ax = viz.figure(8, 5.2)
    for (nombre, s05, nn, krel), color in zip(curvas, [viz.INK] + list(viz.PALETTE[: len(elegidos)])):
        ax.plot(S_g, krel * cin.hill(S_g, 1.0, s05, nn), color=color, lw=2, label=f"{nombre}: S₀.₅ = {s05:g} mM, k_cat ×{krel:.2f}")
    ax.axvline(5.0, color=viz.INK_MUTED, ls="--", lw=1)
    ax.annotate("5 mM de glucosa", (5.0, 0.02), xytext=(4, 0), textcoords="offset points", fontsize=9, color=viz.INK_SECONDARY)
    ax.set_xlabel("[glucosa] (mM)"); ax.set_ylabel("actividad relativa a V_max silvestre")
    ax.legend(fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.18), frameon=False)
    ax.set_title("Mutantes de la glucoquinasa: el interruptor se desplaza", loc="left")
    fig;
    tabla = pd.DataFrame([(m["name"], m["kind"], m["s_half_mm"], m.get("hill_n"), m.get("kcat_s"),
                           round(m.get("kcat_rel", 1.0) * cin.hill(5.0, 1.0, m["s_half_mm"], m.get("hill_n", n_h)), 3)) for m in mut],
                         columns=["mutante", "fenotipo", "S₀.₅ (mM)", "n", "k_cat (s⁻¹)", "actividad a 5 mM (rel.)"])
    tabla.loc[len(tabla)] = ["silvestre", "-", s_half, n_h, valor("kcat_s"), round(cin.hill(5.0, 1.0, s_half, n_h), 3)]
    display(tabla)
    print("Fuente:", ref.get("mutants_source", ""))
else:
    print("(sin datos de mutantes en la tabla de referencia)")
''')

# ============================================================================ 15. inhibición
md(r"""
## 15. Inhibidores y activadores

**La analogía del cajero, ampliada.** Hay tres maneras de frenar la caja:

* **Competitivo**: alguien se pone en la cola *sin comprar nada* y ocupa el turno. Con muchos
  clientes reales, el intruso pasa desapercibido: V<sub>max</sub> no cambia, pero hace falta más
  sustrato para llegar a ella (K<sub>M</sub> aparente sube).
* **Acompetitivo**: alguien traba la caja *solo cuando ya hay un cliente*: baja V<sub>max</sub>
  y, curiosamente, también K<sub>M</sub> (atrapa el complejo ES).
* **No competitivo**: alguien apaga la luz de la caja, haya o no cliente: la caja trabaja más
  lento (baja V<sub>max</sub>) pero los clientes entran igual (K<sub>M</sub> no cambia).

**La derivación para el competitivo (las otras son análogas).** El inhibidor I se une a la
enzima libre: $\mathrm{E + I \rightleftharpoons EI}$ con $K_i = [\mathrm{E}][\mathrm{I}]/[\mathrm{EI}]$.
Ahora la enzima total es $[\mathrm{E}]_0 = [\mathrm{E}] + [\mathrm{ES}] + [\mathrm{EI}]$ y, repitiendo el
álgebra de la sección 13,

$$v_0 = \frac{V_{\max}[\mathrm{S}]}{K_M\left(1 + \dfrac{[\mathrm{I}]}{K_i}\right) + [\mathrm{S}]}
\qquad\Rightarrow\qquad K_M^{\mathrm{ap}} = K_M\left(1+\frac{[\mathrm{I}]}{K_i}\right),\; V_{\max}^{\mathrm{ap}} = V_{\max}.$$

| Tipo | Se une a | V<sub>max</sub> aparente | K<sub>M</sub> aparente | Huella en Lineweaver–Burk |
|---|---|---|---|---|
| competitivo | E | igual | × (1 + [I]/K<sub>i</sub>) | rectas que se cruzan en el eje y |
| acompetitivo | ES | ÷ (1 + [I]/K<sub>i</sub>′) | ÷ (1 + [I]/K<sub>i</sub>′) | rectas paralelas |
| no competitivo | E y ES por igual | ÷ (1 + [I]/K<sub>i</sub>) | igual | rectas que se cruzan en el eje x |
| mixto | E y ES, distinto | ÷ (1 + [I]/K<sub>i</sub>′) | × (1+[I]/K<sub>i</sub>)/(1+[I]/K<sub>i</sub>′) | se cruzan a la izquierda del eje y |

**Qué pasa si…** [I] = K<sub>i</sub> duplica la K<sub>M</sub> aparente de un competitivo; [I] = 10 K<sub>i</sub>
la multiplica por 11. Y un competitivo siempre puede "vencerse" con mucho sustrato; un no
competitivo, nunca.
""")

code(r'''
# 🎛️ Elige el tipo de inhibidor y mueve [I] y Ki: mira la curva y su huella en Lineweaver–Burk
interactivo.explorar_inhibicion()
''')

md(r"""
**Cómo se obtiene K<sub>i</sub> de los datos (tres caminos que deben coincidir).**

1. **Ajuste global**: todos los puntos (varias [S] y varias [I]) contra la ecuación completa
   → V<sub>max</sub>, K<sub>M</sub> y K<sub>i</sub> con sus errores. Es el método recomendado hoy.
2. **Gráfico secundario**: se ajusta cada [I] por separado y se dibuja K<sub>M</sub><sup>ap</sup>
   frente a [I]: es una recta de pendiente K<sub>M</sub>/K<sub>i</sub>.
3. **Gráfico de Dixon**: 1/v₀ frente a [I] para varias [S]; en un competitivo las rectas se
   cruzan en [I] = −K<sub>i</sub>.

Y una conversión útil en farmacología: la **IC₅₀** (concentración que reduce la actividad a la
mitad) depende de [S]; la ecuación de **Cheng–Prusoff** la convierte en K<sub>i</sub>:
$K_i = \mathrm{IC}_{50} / (1 + [\mathrm{S}]/K_M)$ para un competitivo.
""")

code(r'''
Ki_real = 3.0                                   # mM (inhibidor competitivo hipotético)
S_exp = np.array([1, 2, 4, 8, 16, 32, 64]); I_exp = (0, 2, 5, 10)
filas = []
for I in I_exp:
    v = cin.competitive_inhibition(S_exp, Vmax_real, Km_real, I, Ki_real) * (1 + 0.03 * rng.standard_normal(S_exp.size))
    filas += [(s, vv, I) for s, vv in zip(S_exp, v)]
df_i = pd.DataFrame(filas, columns=["S", "v", "I"])
# 1) ajuste global
aji = cin.fit_inhibition(df_i.S.values, df_i.v.values, df_i.I.values, kind="competitive")
print(f"1) Ajuste global:      K_i = {aji['ki']:.2f} ± {aji['ki_err']:.2f} mM   (K_M = {aji['km']:.2f} mM, V_max = {aji['vmax']:.2f} µM/s)")
# 2) gráfico secundario: Km aparente frente a [I]
km_ap = np.array([cin.fit_michaelis_menten(df_i[df_i.I == I].S.values, df_i[df_i.I == I].v.values)["km"] for I in I_exp])
sec = cin.secondary_plot_competitive(np.array(I_exp), km_ap, aji["km"])
print(f"2) Gráfico secundario: K_i = {sec['ki']:.2f} mM (R² = {sec['r2']:.3f})")
# 3) Dixon
dx = cin.dixon_plot(df_i.S.values, df_i.v.values, df_i.I.values)
kd = cin.ki_from_dixon(sorted(dx["lines"]), dx, aji["km"], aji["vmax"])
print(f"3) Dixon:              K_i = {kd['ki']:.2f} mM (por el cruce de las rectas)")
fig, axes = viz.figure(11, 4, ncols=2)
for (s_val, (ii, inv_v)), color in zip(sorted(dx["lines"].items()), viz.sequential_blue(len(dx["lines"]))):
    fit = dx["fits"][s_val]
    xx = np.linspace(-1.3 * Ki_real, max(I_exp), 50)
    axes[0].plot(xx, fit["intercept"] + fit["slope"] * xx, color=color, lw=1.5)
    axes[0].plot(ii, inv_v, "o", color=color, ms=5, label=f"[S] = {s_val:g} mM")
axes[0].axvline(-Ki_real, color=viz.INK_MUTED, ls="--", lw=1); axes[0].set_xlabel("[I] (mM)"); axes[0].set_ylabel("1/v₀ (s/µM)")
axes[0].set_title("Gráfico de Dixon: las rectas se cruzan en [I] = −Kᵢ", loc="left", fontsize=11); axes[0].legend(fontsize=8)
axes[1].plot(I_exp, km_ap, "o", color=viz.COLORS["datos"], ms=6)
xx = np.linspace(0, max(I_exp), 20); axes[1].plot(xx, aji["km"] * (1 + xx / aji["ki"]), color=viz.INK_MUTED, lw=1.5)
axes[1].set_xlabel("[I] (mM)"); axes[1].set_ylabel("K_M aparente (mM)"); axes[1].set_title("Gráfico secundario: pendiente = K_M/Kᵢ", loc="left", fontsize=11)
fig;
# 4) IC50 -> Ki (Cheng–Prusoff)
S_ensayo = 8.0
ic50 = Ki_real * (1 + S_ensayo / Km_real)
print(f"4) Cheng–Prusoff: con [S] = {S_ensayo} mM la IC₅₀ sería {ic50:.1f} mM y devuelve K_i = {cin.cheng_prusoff(ic50, S_ensayo, Km_real):.2f} mM")
''')

md(r"""
**Los inhibidores y reguladores reales de la glucoquinasa** (tabla con fuentes). Tres cosas que
sorprenden:

* La glucosa‑6‑fosfato, que frena a las hexoquinasas I–III, **no** inhibe a la glucoquinasa: por
  eso el hígado sigue fosforilando glucosa aunque el producto se acumule (Viñuela, Salas y Sols, 1963).
* La **manoheptulosa**, que muchos libros llaman "competitiva", resulta de tipo **mixto** cuando se
  mide con cuidado (Scruel 1998), y su K<sub>i</sub> varía 100 veces entre fuentes secundarias: por
  eso no damos un número.
* En el hígado, la **proteína reguladora GKRP** actúa como inhibidor competitivo frente a la glucosa
  (sube S<sub>0.5</sub> sin cambiar V; IC₅₀ ≈ 13 mM de glucosa) y secuestra a la enzima en el núcleo
  cuando la glucosa baja; la fructosa‑6‑fosfato refuerza ese secuestro y la fructosa‑1‑fosfato lo
  deshace.
""")

code(r'''
inh = ref.get("inhibitors", [])
display(pd.DataFrame([(d.get("name"), d.get("kind"), d.get("versus", ""), d.get("ki", ""), d.get("note", ""), d.get("source", "")) for d in inh],
                     columns=["inhibidor / regulador", "tipo", "frente a", "Kᵢ / efecto", "nota", "fuente"]))
act = ref.get("activators", [])
display(pd.DataFrame([(d.get("name"), d.get("fold", ""), d.get("ec50", ""), d.get("note", ""), d.get("source", "")) for d in act],
                     columns=["activador", "activación (veces)", "EC₅₀", "nota", "fuente"]))
''')

md(r"""
**Activadores alostéricos (GKA).** Se unen a un sitio a ~20 Å del sitio de glucosa (Kamata 2004),
bajan S<sub>0.5</sub> y suben V<sub>max</sub>: lo contrario de un inhibidor. El compuesto de
referencia RO‑28‑1675 activa a la enzima silvestre ~16 veces con EC₅₀ ≈ 7 µM (Sayed 2009) y se
ensayó como fármaco para la diabetes tipo 2 (Grimsby 2003). ¿Qué pasa con la actividad a 5 mM de
glucosa cuando bajas S<sub>0.5</sub>? Muévelo:
""")

code(r'''
# 🎛️ Un activador alostérico: baja S0.5 y sube Vmax. ¿Qué pasa con la actividad a 5 mM de glucosa?
interactivo.explorar_activador()
''')

# ============================================================================ 16. T y pH
md(r"""
## 16. Temperatura y pH

**La idea simple.** Calentar acelera casi todas las reacciones (hasta que la enzima se
desnaturaliza). Y cada enzima tiene un pH óptimo: por encima o por debajo, los grupos que hacen
la química pierden la carga que necesitan.

**Arrhenius**: $\ln k = \ln A - E_a/RT$. En un gráfico de ln k frente a 1/T la pendiente es
$-E_a/R$: la **energía de activación** es "cuánto sube la velocidad al calentar".

**Eyring**: $\ln(k/T) = \ln(k_B/h) + \Delta S^{\ddagger}/R - \Delta H^{\ddagger}/RT$. Separa la barrera
en **entalpía** (energía que hay que aportar) y **entropía** (orden que hay que imponer: un
ΔS‡ negativo significa que el estado de transición es más ordenado que los reactivos).

**Qué pasa si…** con E<sub>a</sub> = 12 kcal/mol, subir de 25 a 37 °C multiplica k por ~2.2;
con E<sub>a</sub> = 6 kcal/mol, solo por 1.5. Cuanto mayor la barrera, más sensible a la temperatura.

**pH**: si un grupo debe estar desprotonado (Asp205, la base catalítica) y otro protonado
(Lys169, que estabiliza el fosfato), la actividad tiene forma de campana:

$$v = \frac{v_{\max}}{1 + 10^{pK_1 - \mathrm{pH}} + 10^{\mathrm{pH} - pK_2}}$$

El primer término del denominador "apaga" la enzima a pH bajo (Asp protonado); el segundo la
apaga a pH alto (Lys desprotonada). El máximo está en (pK₁ + pK₂)/2. Para la glucoquinasa
humana el óptimo medido es pH 8.5–8.7 (Šimčíková y Heneberg 2019); durante décadas se creyó más
bajo porque el ATP acidifica los tampones del ensayo: un ejemplo de cómo un detalle técnico
cambia un "hecho" de libro de texto.
""")

code(r'''
# 🎛️ Entalpía y entropía de activación: mira el gráfico de Eyring y k(T)
interactivo.explorar_temperatura()
''')

code(r'''
T_K = np.array([283.15, 288.15, 293.15, 298.15, 303.15, 308.15, 313.15])
dH, dS = 12.0, -8.0                       # kcal/mol, cal/mol/K (ilustrativos)
k_T = np.array([cin.eyring_rate(dH - T * dS / 1000, T) for T in T_K]) * (1 + 0.03 * rng.standard_normal(T_K.size))
aj_arr, aj_eyr = cin.fit_arrhenius(T_K, k_T), cin.fit_eyring(T_K, k_T)
fig, axes = viz.figure(11, 4, ncols=2)
viz.plot_arrhenius(T_K, k_T, fit=aj_arr, ax=axes[0]); viz.plot_eyring(T_K, k_T, fit=aj_eyr, ax=axes[1])
fig;
print(f"Arrhenius: E_a = {aj_arr['ea_kcal']:.1f} kcal/mol.   Eyring: ΔH‡ = {aj_eyr['delta_h_kcal']:.1f} kcal/mol, ΔS‡ = {aj_eyr['delta_s_cal']:.1f} cal/mol/K, ΔG‡(25 °C) = {aj_eyr['delta_g_kcal']:.1f} kcal/mol")
print("Relación entre ambas: E_a ≈ ΔH‡ + RT (0.6 kcal/mol a 25 °C).")
''')

code(r'''
# 🎛️ pKa de la base y del ácido: la campana de pH
interactivo.explorar_ph()
''')

# ============================================================================ 17. resumen
md(r"""
## 17. Resumen, glosario y ejercicios

**Lo que hemos visto, en una frase cada cosa:**

1. Una enzima acelera una reacción **bajando la barrera** del estado de transición, como un guía
   que conoce el paso más bajo.
2. La estructura cristalina (3FGU) muestra a los reactivos ya alineados; la MD muestra que la
   enzima **los mantiene alineados**.
3. QM/MM permite **ver la química** dentro de la enzima; hay que cuidar la partición, los átomos
   de enlace, la electrostática y la repulsión QM–MM.
4. El estado de transición es un **punto de silla**: se localiza (escaneo → NEB → dímero) y se
   **valida** con una única frecuencia imaginaria y con el camino que baja hacia R y P.
5. De ΔE‡ a ΔG‡: vibraciones, promedio sobre conformaciones y sensibilidad al método. Un
   resultado computacional lleva **incertidumbre**, y nuestra barrera coincide con otro cálculo
   publicado y sobrestima el experimento en ~3–4 kcal/mol.
6. Eyring convierte la barrera en *k*<sub>cat</sub>; 1.4 kcal/mol son un factor 10.
7. Michaelis–Menten **emerge** del mecanismo E + S ⇌ ES → E + P; K<sub>M</sub>, V<sub>max</sub>,
   k<sub>cat</sub> y k<sub>cat</sub>/K<sub>M</sub> se derivan de los datos por ajuste no lineal.
8. La glucoquinasa es **sigmoide** (n ≈ 1.7) sin tener varias subunidades: cooperatividad
   cinética; sus mutantes desplazan el interruptor de la insulina.
9. Cada tipo de inhibidor deja una huella distinta; K<sub>i</sub> se obtiene por ajuste global,
   gráfico secundario o Dixon; los activadores hacen lo contrario.
10. Temperatura y pH modulan k a través de ΔH‡, ΔS‡ y los pK<sub>a</sub> del sitio activo.

**Glosario rápido.** *Sustrato*: molécula transformada. *Sitio activo*: donde ocurre la
química. *ES*: complejo enzima‑sustrato. *Estado de transición*: cima de la barrera. *ΔG‡*:
altura de la barrera. *k<sub>cat</sub>*: reacciones por segundo por enzima. *K<sub>M</sub>*: [S] a
media velocidad. *k<sub>cat</sub>/K<sub>M</sub>*: eficiencia a baja [S]. *S<sub>0.5</sub>, n*: parámetros de
Hill. *K<sub>i</sub>*: constante de disociación del inhibidor. *QM/MM*: mecánica cuántica para la
región reactiva y campo de fuerza para el resto. *NEB, dímero*: métodos para hallar el TS.
*Frecuencia imaginaria*: la firma de un punto de silla.

**Ejercicios**

1. Cambia `umbral` en la sección 5 a 3.0 y 4.0 Å. ¿Cómo cambia la fracción de conformaciones de
   ataque cercano? ¿Qué pasaría con *k*<sub>cat</sub> si la enzima no cerrara sus dominios?
2. Con el explorador de Eyring, averigua cuánto tendría que bajar la barrera para multiplicar
   *k*<sub>cat</sub> por 1000. Compáralo con la diferencia entre el sitio activo en agua y en la enzima.
3. En la sección 13, cambia `k2` a 6 s⁻¹ y a 600 s⁻¹. ¿Cómo cambian K<sub>M</sub> y V<sub>max</sub>?
   ¿Cuándo K<sub>M</sub> ≈ K<sub>d</sub> = k<sub>−1</sub>/k<sub>1</sub>?
4. Simula datos con n = 1.0, 1.4 y 2.0 en la sección 14 y ajústalos con Michaelis–Menten.
   ¿Qué error cometes si ignoras la cooperatividad? ¿Qué mutante de la tabla dejaría de liberar
   insulina a 5 mM?
5. Diseña un experimento (concentraciones de S e I) que distinga un inhibidor competitivo de uno
   mixto con K<sub>i</sub>′ = 3K<sub>i</sub>. Usa `cin.fit_inhibition` y el gráfico de Dixon para comprobarlo.
6. Con la termoquímica de la sección 10, ¿cuánto cambia ΔG‡ entre 25 y 37 °C? ¿Qué factor en
   *k* supone? Compáralo con lo que predice el explorador de temperatura.
7. (Avanzado, modo `completo`) En `enzimas/glucoquinasa.py` saca Lys169 de `QM_RESIDUES` (quedará
   como cargas MM) y recalcula el escaneo. ¿Sube o baja la barrera? Es un "mutante computacional";
   Zhang et al. (2009) obtuvieron 32 kcal/mol para K169A frente a 18 para la silvestre.

**Referencias**
""")

code(r'''
for r in ref.get("references", []):
    print("•", r)
print("\nHerramientas QM/MM adaptadas de la suite Leonardo (Juvenal Yosa, MIT): https://github.com/juvenalyosa/Leonardo")
''')


def main():
    nb = nbf.v4.new_notebook()
    nb["cells"] = CELLS
    nb["metadata"] = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
        "colab": {"name": "Cinética enzimática con QM/MM: glucoquinasa", "provenance": []},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(nb, str(OUT))
    print("escrito", OUT, "con", len(CELLS), "celdas")


if __name__ == "__main__":
    main()
