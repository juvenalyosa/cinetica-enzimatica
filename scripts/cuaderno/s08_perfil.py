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
# @title ⚗️ Reactivo y producto dentro de la enzima
if MODO == "completo":
    subprocess.run([sys.executable, "scripts/03_qmmm_reaccion.py"], check=True)

res = datos.json_("qmmm/resumen.json")
R, P = res["etapas"]["reactivo"], res["etapas"]["producto"]
viz.tarjetas(
    [("Reactivo · ξ", f"{R['xi']:+.2f}".replace("-", "−"), "Å",
      f"Pγ–O3β {R['d_PG_O3B']:.2f} Å (enlazado) · Pγ–O6 {R['d_PG_O6']:.2f} Å (libre)", "azul"),
     ("Producto · ξ", f"{P['xi']:+.2f}", "Å",
      f"Pγ–O6 {P['d_PG_O6']:.2f} Å (enlazado) · Pγ–O3β {P['d_PG_O3B']:.2f} Å (libre)", "agua"),
     ("El protón del O6", f"{P['d_OD1_H']:.2f}", "Å",
      f"en el producto está en Asp205 (en el reactivo, O6–H = {R['d_O6_H']:.2f} Å)", "violeta"),
     ("ΔE de reacción", f"{P['dE_reaccion_kcal']:+.1f}", "kcal/mol", "cuesta arriba dentro de la enzima (PM7, entorno fijo)", "naranja")],
    titulo="Los dos extremos de la reacción, relajados en la enzima")
''')

code(r'''
# @title 📈 El perfil de energía de la reacción
escaneo = datos.csv("qmmm/escaneo.csv")
xi_esc = np.r_[R["xi"], escaneo["xi"].values]                   # el perfil arranca en el reactivo relajado (E = 0)
e_esc = np.r_[0.0, escaneo["energia_rel_kcal"].values]
fig = viz.plot_energy_profile(xi_esc, e_esc, relative=False, xlabel="ξ = d(Pγ–O3β) − d(Pγ–O6)  (Å)",
                              ts_index=int(np.argmax(e_esc)), annotate_states=True,
                              state_names=("reactivo\nfosfato en el ATP", "cima: 1.ª estimación del TS", "producto\nfosfato en la glucosa"),
                              title=f"Para pasar el fosfato hay que subir una colina de {e_esc.max():.0f} kcal/mol",
                              subtitle="Escaneo relajado: se fija ξ en cada punto y se relaja todo lo demás (PM7 en el campo de la enzima)")
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

md(r"""
### 🎬 Mírala ocurrir

Ahora la película: los átomos a la izquierda y, a la derecha, el **perfil de energía con un punto
amarillo que sigue a la película**. Pulsa ▶ o arrastra la barra. Fíjate en tres cosas:

1. El enlace **rosa** (Pγ–O3β, del ATP) se adelgaza y se vuelve discontinuo: **se está rompiendo**.
2. El enlace **verde agua** (Pγ–O6, hacia la glucosa) aparece y se engrosa: **se está formando**.
3. En la cima de la gráfica (el estado de transición) los dos están a medias: el fósforo está
   **entre** los dos oxígenos.

Solo se dibujan enlaces reales: los de la topología del modelo, más estos cuatro que cambian, con un
grosor proporcional a su orden de enlace (enteros si están formados, discontinuos si están a medio
formar o romper). Las líneas verdes finas son la coordinación del Mg²⁺, que no es un enlace covalente.
La película sigue el camino de mínima energía que se calcula en la sección 9, desde el reactivo
relajado hasta el producto.
""")

code(r'''
# @title 🎬 La reacción en movimiento, sincronizada con su energía
visor3d.pelicula_reaccion("enzima")
''')

md(r"""
> ✅ **Para llevar.** ξ = d(Pγ–O3β) − d(Pγ–O6) resume la reacción en un solo número: negativa
> con el fosfato en el ATP, cerca de cero en la cima, positiva con el fosfato en la glucosa. Un
> escaneo relajado a lo largo de ξ da el primer perfil de energía y una primera estimación de la
> barrera, que la sección 9 afina.
""")
