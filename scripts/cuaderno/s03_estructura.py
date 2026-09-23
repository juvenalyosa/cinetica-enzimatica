"""Sección 3. estructura."""
from .celdas import code, md

# ============================================================================ 3. estructura
md(r"""
## 3. La estructura 3D del complejo catalítico

> 🎯 **En esta sección** verás dónde está cada átomo de la glucoquinasa justo antes de la
> reacción, qué hace cada pieza del sitio activo y por qué los reactivos ya están "apuntándose".

---

### 💡 La analogía: una almeja que se cierra sobre la glucosa

Para simular una enzima necesitamos saber dónde está cada átomo. Eso lo dan los cristalógrafos:
la estructura **3FGU** del Protein Data Bank es una "foto" de la glucoquinasa humana con
**glucosa**, un análogo de ATP y el ion **Mg²⁺**, todos dentro del sitio activo, justo antes de la
reacción.

La glucoquinasa tiene dos partes, el **dominio grande** y el **dominio pequeño**, unidas por una
bisagra. Sin glucosa la enzima está **abierta**, como una almeja entreabierta. Cuando la glucosa
entra, el dominio pequeño se cierra sobre ella: los reactivos quedan encerrados, sin agua
alrededor y alineados para reaccionar.

[[fig:almeja_dominios | La glucoquinasa abierta y cerrada: el dominio pequeño se cierra sobre la glucosa como una almeja]]

---

### 🔬 Los datos reales: quién es quién en el sitio activo

En la forma cerrada, cada pieza tiene un papel:

| Pieza | Dónde está | Qué hace | En la analogía |
|---|---|---|---|
| **O6 de la glucosa** | a ~2.7 Å del fósforo Pγ del ATP | recibirá el fosfato | la mano que atrapa la pelota |
| **Pγ del ATP** | al final del trifosfato | es el fósforo que salta | la pelota |
| **Asp205** | a 2.5 Å del hidroxilo O6–H | **base catalítica**: se llevará el protón | quien libera la mano para atrapar |
| **Lys169** y **Mg²⁺** | junto a los fosfatos | neutralizan sus cargas negativas | amortiguadores |

Un enlace P–O mide 1.6 Å: a 2.7 Å el O6 y el Pγ están **muy cerca, pero aún no unidos**.

[[fig:sitio_activo | Esquema del sitio activo de la glucoquinasa en el cristal 3FGU con las distancias clave]]

Explora la enzima real, átomo por átomo: **arrastra** para girarla, usa la **rueda** para acercarte
y los botones para ver el **sitio activo** o la **superficie**:
""")

code(r'''
# @title 🧬 La glucoquinasa en 3D (gírala con el ratón)
pdb_cristal = pathlib.Path("data/raw/3FGU.pdb").read_text()
visor3d.complejo_cristal(pdb_cristal)
''')

md(r"""
Y comprobemos las dos distancias clave directamente en el archivo del cristal:
""")

code(r'''
# @title 📏 Distancias clave en el cristal
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
> ✍️ **Para pensar.** El cristal usa AMP‑PNP, en el que el oxígeno entre los fósforos β y γ se
> cambió por un nitrógeno (N3B). Ese cambio impide la reacción y permite capturar el complejo.
> Para simular la reacción real, nosotros volveremos a poner el oxígeno: ATP verdadero.

[[fig:ampnp_atp | AMP-PNP con un puente de nitrógeno que no se rompe frente al ATP con un puente de oxígeno que sí se rompe]]

> ✅ **Para llevar.** La glucoquinasa se cierra como una almeja sobre la glucosa y deja al O6 a
> 2.7 Å del fósforo que debe atacar, con Asp205 listo para tomar su protón y Mg²⁺ y Lys169
> sujetando las cargas del fosfato. La foto ya muestra una enzima **preparada para reaccionar**.
""")
