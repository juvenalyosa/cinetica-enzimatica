"""Sección 7. QM/MM."""
from .celdas import code, md

# ============================================================================ 7. QM/MM

md(r"""
## 7. Mirar la reacción con mecánica cuántica: QM/MM

> 🎯 **En esta sección** entenderás por qué hace falta mecánica cuántica para ver cómo el
> fosfato salta del ATP a la glucosa, por qué no se puede usar para toda la enzima y cómo el
> método **QM/MM** combina lo mejor de los dos mundos.

---

### 💡 La analogía: fotografiar un partido

Para fotografiar un partido no necesitas enfocar a los 50 000 espectadores: enfocas el balón y
los jugadores cercanos, y el estadio queda como fondo. En la enzima, "el balón" es el fosfato que
salta; el estadio, el resto de la proteína y el agua.

[[fig:qmmm_estadio | Analogía del estadio: se enfoca el balón (región QM) y las gradas quedan de fondo (región MM)]]

**¿Por qué no basta la dinámica molecular?** La dinámica molecular clásica de la sección 5 no
puede romper ni formar enlaces: sus resortes no saben de electrones. Para ver cómo el fosfato
salta del ATP a la glucosa hay que usar **mecánica cuántica (QM)**, que sí describe los
electrones. Pero es carísima: imposible para 57 000 átomos.

**La solución: QM/MM.** Se trata con QM solo el pedacito donde ocurre la química (la **región
QM**, unas decenas de átomos) y el resto con el campo de fuerza clásico (**región MM**). Las
dos regiones "se ven":

* las cargas de la región MM crean un potencial eléctrico que entra en el cálculo cuántico
  (**embedding electrostático**), y
* la región MM empuja a la QM con fuerzas de van der Waals (**Lennard‑Jones**): sin esto la
  región QM se "hundiría" en la proteína, porque la electrostática sola atrae.

---

### 🧮 La ecuación en palabras

**¿Para qué queremos una ecuación?** Para buscar el camino de la reacción, el programa necesita
la energía (y las fuerzas) de cada geometría que prueba. En QM/MM esa energía tiene dos partes:

$$\text{energía} \;=\; \underbrace{\text{energía cuántica de la región QM (en el campo de la enzima)}}_{\text{la foto enfocada}} \;+\; \underbrace{\text{los «codos» entre QM y MM}}_{\text{no atravesar al público}}$$

Vuelve al estadio: los jugadores que enfocamos (la región QM) oyen los gritos de las gradas (las
cargas de la enzima, que atraen o repelen a sus electrones), pero no pueden meterse *dentro* del
público: los codos de la gente de alrededor se lo impiden (la repulsión de corto alcance).

[[fig:qmmm_regiones | Región QM rodeada por la región MM con cargas puntuales; átomo de enlace en la frontera; electrostática y Lennard-Jones]]

---

### 📐 La ecuación completa, término a término

$$E \;=\; {\color{#eb6834}{E_{\mathrm{PM7}}\big[\text{QM en el potencial de la enzima}\big]}} \;+\; {\color{#4a3aa7}{E_{\mathrm{LJ}}(\text{QM–MM})}}$$

| Término | Qué es | En la analogía | En nuestro modelo |
|---|---|---|---|
| ${\color{#eb6834}{E_{\mathrm{PM7}}[\dots]}}$ | energía de electrones y núcleos de la región QM, calculada con el método semiempírico PM7 | el grupo del centro, bien enfocado | 78 átomos, carga −2 |
| potencial de la enzima | el campo eléctrico de las cargas MM, incluido *dentro* del cálculo cuántico (embedding) | la música y los gritos del público | ~2300 cargas Amber dentro de 16 Å, fijas |
| ${\color{#4a3aa7}{E_{\mathrm{LJ}}}}$ | Lennard‑Jones entre átomos QM y MM: atracción débil de lejos, repulsión fuerte de cerca | los codos del público | parámetros Amber |
| átomo de enlace | un H que "tapa" cada enlace covalente cortado por la frontera | — | cadenas laterales cortadas en la frontera |

**Nuestra receta (herramientas de la suite Leonardo, MOPAC PM7 + Amber):**

* **Región QM** (78 átomos, carga −2): glucosa, el fragmento C5'–trifosfato del ATP, las cadenas
  laterales de Asp205, Lys169 y Thr228, el Mg²⁺ y su agua coordinada.
* Donde un enlace covalente cruza la frontera se pone un **átomo de enlace** (un H).
* El **entorno MM** (2300 cargas dentro de 16 Å) se mantiene **fijo** (aproximación de
  entorno rígido).
* Energía total: $E = E_{\mathrm{PM7}}[\text{QM en el potencial de la enzima}] + E_{\mathrm{LJ}}(\text{QM–MM})$.

---

### 🎛️ Qué pasa si…

| Si… | entonces… |
|---|---|
| **quitamos el potencial de la enzima** | la región QM ya no "ve" a la proteína: es como sacarla de ella (la sección 11 lo hace, en agua) |
| **quitamos el término Lennard‑Jones** | la electrostática, que solo atrae, "hunde" la región QM en la proteína y la geometría deja de tener sentido |
| **agrandamos la región QM** | el modelo es más fiel, pero cada cálculo cuesta mucho más |
| **dejamos que el entorno MM se mueva** | el cálculo sería más realista y mucho más caro; con el entorno fijo, la barrera depende de la conformación elegida (sección 10) |

---

### 🔬 Los datos reales: nuestra partición
""")

code(r'''
# @title 🔍 La región cuántica dentro de la enzima
particion = datos.json_("qmmm/particion.json")
print(f"Región QM: {particion['n_qm']} átomos + {particion['n_link']} H de enlace; libres: {particion['n_free']}")
print(f"Entorno MM: {particion['n_mm']} cargas puntuales; carga QM = {particion['qm_charge']}; carga MM = {particion['mm_charge']:.2f} e")
print("Átomos QM por residuo:", pd.Series([a["resname"] for a in particion["qm_atoms"]]).value_counts().to_dict())
display(visor3d.region_qm())
''')

md(r"""
**¿Cómo sabemos que el programa calcula bien las fuerzas?** Comparando el gradiente analítico
con el numérico (mover un átomo 0.005 Å y ver cuánto cambia la energía). Es una comprobación
que todo código de simulación debería pasar antes de usarse:
""")

code(r'''
# @title 🧪 Control de calidad: fuerzas analíticas frente a numéricas
grad = datos.csv("qmmm/comprobacion_gradiente.csv")
grad["error relativo"] = (grad["diff"].abs() / grad["analytic"].abs().clip(lower=1e-6)).map(lambda x: f"{x:.1%}")
display(grad.round(3))
''')

md(r"""
> ✅ **Para llevar.** QM/MM enfoca la mecánica cuántica donde ocurre la química (78 átomos) y
> trata el resto de la enzima como cargas y "codos" clásicos. La energía tiene dos partes: la
> energía cuántica de la región QM **dentro del campo eléctrico de la enzima** y la repulsión
> Lennard‑Jones que le impide hundirse en ella.
""")
