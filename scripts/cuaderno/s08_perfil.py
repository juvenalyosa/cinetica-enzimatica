"""Sección 8. R, P y escaneo."""
from .celdas import code, md

# ============================================================================ 8. R, P y escaneo

md(r"""
## 8. Reactivo, producto y el perfil de energía

> 🎯 **En esta sección** describirás toda la reacción con **un solo número**, la coordenada de
> reacción ξ, y obtendrás el primer perfil de energía de la transferencia del fosfato dentro de
> la enzima.

---

### 💡 La analogía: el pañuelo del tira y afloja

En un tira y afloja, el pañuelo atado al centro de la cuerda dice quién va ganando: si está del
lado de un equipo, gana ese equipo; si está en el centro, la cosa está igualada. Aquí la cuerda
es el grupo fosfato (el Pγ), y los dos "equipos" son el oxígeno O3β del ATP (que lo suelta) y el
oxígeno O6 de la glucosa (que lo recibe).

Primero relajamos la región QM en su punto de partida (el **reactivo**: glucosa + ATP). Luego
construimos el **producto** (glucosa‑6‑fosfato + ADP, con el protón de O6 ya en Asp205) y
también lo relajamos. La diferencia de energía dice si la reacción es "cuesta abajo" o "cuesta
arriba" dentro de la enzima.

[[fig:xi_coordenada | Tres viñetas del fosfato pasando del O3β del ATP al O6 de la glucosa, con las dos distancias y el valor de ξ]]

---

### 🧮 La ecuación en palabras

**¿Para qué queremos una ecuación?** Una reacción mueve decenas de átomos a la vez; para dibujar
un perfil de energía necesitamos un eje horizontal, un número que avance de forma continua
desde el reactivo hasta el producto.

$$\xi \;=\; \text{cuánto se ha alejado el fósforo del ADP} \;-\; \text{cuánto le falta para llegar a la glucosa}$$

---

### 📐 La ecuación completa, término a término

$$\xi \;=\; {\color{#2a78d6}{d(\mathrm{P_\gamma\!-\!O_{3\beta}})}} \;-\; {\color{#1baf7a}{d(\mathrm{P_\gamma\!-\!O_6})}}$$

El término **azul** es el enlace que se rompe (lado del reactivo); el **aguamarina**, el que se
forma (lado del producto).

| Situación | d(Pγ–O3β) | d(Pγ–O6) | ξ | En la analogía | En nuestro cálculo |
|---|---|---|---|---|---|
| reactivo (fosfato en el ATP) | corta (1.6 Å) | larga (3.2 Å) | negativa | el pañuelo en el campo del ATP | ξ = −1.56 Å |
| estado de transición | intermedia | intermedia | ≈ 0 | el pañuelo en el centro | ξ = +0.12 Å |
| producto (fosfato en la glucosa) | larga | corta (1.7 Å) | positiva | el pañuelo en el campo de la glucosa | ξ = +1.75 Å |

---

### 🎛️ Qué pasa si…

| Si… | entonces… | porque… |
|---|---|---|
| **el fósforo avanza 0.1 Å** hacia la glucosa, en línea | ξ sube unos **0.2 Å** | una distancia crece y la otra se acorta a la vez |
| **ξ = 0** | el fósforo está a igual distancia de los dos oxígenos, pero eso **no garantiza** que sea la cima | la cima es el máximo de energía, que se busca aparte (sección 9) |
| **fijamos ξ** en un valor intermedio y dejamos relajar lo demás | obtenemos un punto del **escaneo relajado** | la restricción sostiene la geometría "a medio camino" |

Un **escaneo relajado** fija ξ en valores intermedios (con una restricción armónica) y deja que
todo lo demás se acomode: el resultado es un **perfil de energía** aproximado, y su punto más
alto es una primera estimación del estado de transición.

---

### 🔬 Los datos reales
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

md(r"""
> ✅ **Para llevar.** ξ = d(Pγ–O3β) − d(Pγ–O6) resume la reacción en un solo número: negativa
> con el fosfato en el ATP, cerca de cero en la cima, positiva con el fosfato en la glucosa. Un
> escaneo relajado a lo largo de ξ da el primer perfil de energía y una primera estimación de la
> barrera, que la sección 9 afina.
""")
