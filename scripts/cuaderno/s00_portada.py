"""Sección 0. portada."""
from .celdas import code, md

# ============================================================================ 0. portada
md(r"""
# Cinética enzimática, paso a paso, con simulaciones QM/MM

**Enzima de estudio: la glucoquinasa humana** (hexoquinasa IV, EC 2.7.1.2)

$$\text{glucosa} + \text{ATP} \;\xrightarrow{\;\text{glucoquinasa}\;}\; \text{glucosa‑6‑fosfato} + \text{ADP}$$

[[fig:portada_banner | La reacción de la glucoquinasa y el recorrido del curso: estructura, dinámica molecular, química cuántica y cinética]]

Este cuaderno explica la cinética enzimática **desde cero**, para estudiantes de pregrado. No
necesitas saber programar: basta con ejecutar las celdas en orden (▶ o `Shift + Enter`) y leer.
Usamos una enzima real y todo lo que se puede usar con ella: una estructura cristalina,
simulaciones de dinámica molecular y de mecánica cuántica (QM/MM), y los números medidos en el
laboratorio para esta enzima, con sus fuentes.

---

### 🧭 Cómo está organizada cada sección

Todas las secciones siguen el mismo orden, marcado siempre con los mismos símbolos:

| Símbolo | Bloque | Qué encontrarás |
|:---:|---|---|
| 🎯 | **En esta sección** | la pregunta que vamos a responder |
| 💡 | **La analogía** | una imagen de la vida diaria, con un dibujo |
| 🧮 | **La ecuación en palabras** | por qué necesitamos una ecuación y qué dice, sin símbolos |
| 📐 | **La ecuación completa** | la versión formal, término a término, con los mismos colores del dibujo |
| 🎛️ | **Qué pasa si…** | gráficas con deslizadores: mueve un término y mira qué cambia |
| 🔬 | **Los datos reales** | simulaciones y valores medidos en el laboratorio |
| ✅ | **Para llevar** | la idea clave en dos líneas |

Además, los recuadros **⚠️ Ojo** avisan de una trampa frecuente y los **✍️ Para pensar** te dejan una pregunta.

---

### 🗺️ Índice

| # | Sección | Qué aprendes | Herramienta |
|---|---|---|---|
| | **🧱 Bloque 1 · Las bases** | | |
| 1 | ¿Qué es una enzima y por qué la glucoquinasa? | vocabulario básico, datos experimentales | lectura |
| 2 | ¿Qué es la velocidad de una reacción? | v₀, unidades, cómo se mide | simulación |
| | **🔭 Bloque 2 · Ver la enzima** | | |
| 3 | La estructura 3D del complejo catalítico | sitio activo | PDB 3FGU, visor 3D |
| 4 | Preparar la enzima para simularla | de la foto al modelo | PDBFixer, Amber |
| 5 | La enzima se mueve: dinámica molecular | conformaciones reactivas | OpenMM |
| | **⚛️ Bloque 3 · La química, con mecánica cuántica** | | |
| 6 | Colinas de energía y velocidad | ecuación de Eyring | 🎛️ interactivo |
| 7 | Mirar la reacción con mecánica cuántica: QM/MM | regiones QM y MM | MOPAC PM7 + Amber |
| 8 | Reactivo, producto y perfil de energía | coordenada de reacción | escaneo |
| 9 | El estado de transición | NEB, dímero, frecuencia imaginaria | ASE |
| 10 | De ΔE‡ a ΔG‡: mejoras del cálculo | termoquímica, promedios, método | ASE, MOPAC |
| 11 | La reacción fuera de la enzima | SADDLE (QST2), TS, FORCETS | MOPAC nativo |
| 12 | De la barrera a *k*<sub>cat</sub> | comparación con el experimento | Eyring |
| | **📈 Bloque 4 · La cinética del laboratorio** | | |
| 13 | Michaelis–Menten desde el mecanismo | K<sub>M</sub>, V<sub>max</sub>, k<sub>cat</sub>, k<sub>cat</sub>/K<sub>M</sub> | ODE, ajustes, 🎛️ |
| 14 | Cooperatividad: la glucoquinasa es un sensor | ecuación de Hill, mutantes | 🎛️ |
| 15 | Inhibidores y activadores | tipos, K<sub>i</sub>, Dixon, IC₅₀ | 🎛️ |
| 16 | Temperatura y pH | Arrhenius, Eyring, pKa | 🎛️ |
| | **🏁 Cierre** | | |
| 17 | Resumen, glosario y ejercicios | mapa del curso | |

> 💡 **Consejo.** En Colab, abre el panel **Índice** (icono ☰ a la izquierda) para saltar entre
> secciones.

---

### ⚙️ Antes de empezar

**Dos modos de ejecución.** En modo `rapido` (por defecto) el cuaderno lee resultados ya
calculados y dibuja todo en unos minutos. En modo `completo` reejecuta las simulaciones
(horas en la CPU de Colab). Cambia `MODO` en la celda siguiente si quieres recalcular.

> 🎛️ **Cómo usar las gráficas interactivas.** Las celdas marcadas con 🎛️ muestran deslizadores.
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
