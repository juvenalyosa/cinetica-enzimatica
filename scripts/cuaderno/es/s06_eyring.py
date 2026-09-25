"""Sección 6. Eyring."""
from ..celdas import code, md

# ============================================================================ 6. Eyring
md(r"""
## 6. Colinas de energía y velocidad de reacción

> 🎯 **En esta sección** entenderás por qué una colina de energía más baja significa una reacción
> más rápida, y cómo la ecuación de Eyring convierte la altura de esa colina en un número de
> reacciones por segundo.

---

### 💡 La analogía: una pelota que intenta cruzar una colina

Imagina una pelota en el fondo de un valle. Entre ese valle (los **reactivos**: glucosa + ATP) y
el valle vecino (los **productos**: glucosa‑6‑fosfato + ADP) hay una colina. La pelota no está
quieta: el calor la sacude sin parar y la lanza cuesta arriba una y otra vez. Casi siempre se
queda corta y rueda de vuelta. Muy de vez en cuando, un empujón especialmente fuerte la lleva
hasta la cima, y entonces cae al otro lado: **la reacción ocurrió**.

[[fig:colina_intentos | Una pelota en un valle intenta cruzar una colina; la mayoría de los intentos falla y uno llega a la cima]]

En química la colina se llama **barrera de activación** y su cima, **estado de transición**.
Una enzima no empuja más fuerte a la pelota: **rebaja la colina**.

---

### 🧮 La ecuación en palabras

**¿Para qué queremos una ecuación?** Para pasar de "la colina es alta" a "la reacción ocurre
60 veces por segundo". Con un número podemos comparar la simulación con el experimento.

La idea es la de cualquier juego de azar: si tiras un dado muchas veces, el número de "seises"
por minuto es *(tiradas por minuto)* × *(probabilidad de sacar seis)*. Aquí:

$$\text{velocidad} \;=\; \underbrace{\text{intentos por segundo}}_{\text{muchísimos}} \;\times\; \underbrace{\text{probabilidad de llegar a la cima}}_{\text{diminuta}}$$

[[fig:eyring_anatomia | La ecuación de Eyring dividida en dos factores: intentos por segundo y probabilidad de éxito]]

---

### 📐 La ecuación completa, término a término (Eyring, 1935)

$$k \;=\; {\color{#2a78d6}{\kappa\,\frac{k_B T}{h}}}\;\times\;{\color{#eb6834}{\exp\!\left(-\frac{\Delta G^{\ddagger}}{RT}\right)}}$$

El factor **azul** son los intentos; el **naranja**, la probabilidad de que un intento llegue a la
cima. Cada símbolo:

| Término | Qué es | En la analogía | Valor típico |
|---|---|---|---|
| $k$ | constante de velocidad | cuántas pelotas cruzan por segundo | ≈ 60 s⁻¹ para la glucoquinasa |
| $k_B T/h$ | frecuencia de intentos | cuántas veces por segundo la pelota se lanza cuesta arriba | 6.2 × 10¹² s⁻¹ a 25 °C |
| $\Delta G^{\ddagger}$ | energía libre de activación | la altura de la colina | 10–25 kcal/mol |
| $R T$ | energía térmica disponible | la fuerza típica de las sacudidas | 0.59 kcal/mol a 25 °C |
| $\kappa$ | coeficiente de transmisión | fracción de pelotas que, ya en la cima, no se arrepienten y vuelven | ≈ 1 |

¿Por qué una **exponencial**? Porque llegar alto exige muchas sacudidas favorables seguidas, y
la probabilidad de que ocurran a la vez se multiplica: cada escalón extra de altura divide la
probabilidad por el mismo factor. Por eso lo que importa es el cociente $\Delta G^{\ddagger}/RT$:
cuántas "sacudidas típicas" mide la colina.

---

### 🎛️ Qué pasa si…

| Si… | entonces… | porque… |
|---|---|---|
| **ΔG‡ sube 1.36 kcal/mol** | la reacción es **10 veces más lenta** | la exponencial cae a la décima parte (la regla de oro) |
| **T sube de 25 a 37 °C** | para una barrera de 15 kcal/mol, la velocidad **casi se duplica** | las sacudidas (RT) son un 4 % más fuertes y la exponencial lo amplifica |
| **la enzima baja ΔG‡ 10 kcal/mol** | la reacción va **10⁷ veces** más rápido | 10 / 1.36 ≈ 7 escalones de ×10 |

Muévelo tú mismo con los deslizadores:
""")

code(r'''
# @title 🎛️ Explora la ecuación de Eyring
# 🎛️ Explora la ecuación de Eyring: barrera y temperatura
interactivo.explorar_eyring()
''')

code(r'''
# @title 🎛️ Dibuja tu propia colina de energía
# 🎛️ Perfil de energía: barrera (cuesta) y energía de reacción (desnivel entre valles)
interactivo.explorar_perfil_energia()
''')

code(r'''
# @title 📉 Cada 1.36 kcal/mol, un factor 10
barreras = np.linspace(5, 30, 200)
kcat = valor("kcat_s", 60.0)
fig, ax = viz.figure(10.5, 5.8)
k_curva = cin.eyring_rate(barreras)
ax.fill_between(barreras, k_curva.min() / 1e3, k_curva, color=viz._tint(viz.COLORS["reactivo"], 0.9), zorder=1, lw=0)
ax.semilogy(barreras, k_curva, color=viz.COLORS["reactivo"], lw=2.8, zorder=3)
puntos = [(cin.barrier_from_rate(kcat), f"la glucoquinasa\nk_cat ≈ {kcat:.0f} s⁻¹", viz.COLORS["ts"], (12, 6), "left"),
          (cin.barrier_from_rate(1e-8), "una reacción que\ntarda años", viz.INK_SECONDARY, (-16, -10), "right")]
for dg, texto, color, off, ha in puntos:
    viz._guide(ax, "v", dg, start=k_curva.min() / 1e3, end=cin.eyring_rate(dg))
    viz._keypoint(ax, dg, cin.eyring_rate(dg), color, size=11)
    ax.annotate(f"{texto}\nΔG‡ ≈ {dg:.1f} kcal/mol", (dg, cin.eyring_rate(dg)), xytext=off, textcoords="offset points",
                ha=ha, va="bottom" if off[1] > 0 else "top", fontsize=11, color=viz.INK)
# la escalera: 3 escalones de 1.36 kcal/mol desde la barrera de la glucoquinasa
dg0 = cin.barrier_from_rate(kcat)
for j in range(3):
    a, b = dg0 + j * 1.364, dg0 + (j + 1) * 1.364
    ax.plot([a, b, b], [cin.eyring_rate(a), cin.eyring_rate(a), cin.eyring_rate(b)], color=viz.INK_MUTED, lw=1.2, zorder=4)
    ax.annotate("÷10", (b, np.sqrt(cin.eyring_rate(a) * cin.eyring_rate(b))), xytext=(4, 0), textcoords="offset points",
                ha="left", va="center", fontsize=10, color=viz.INK_SECONDARY)
ax.set_xlim(5, 30); ax.set_ylim(k_curva.min() / 1e3, k_curva.max() * 10)
viz._finish(ax, "ΔG‡: altura de la colina (kcal/mol)", "k (s⁻¹, escala log)",
            "Cada 1.36 kcal/mol de barrera, la velocidad cambia 10 veces",
            "Cada escalón gris sube la colina 1.36 kcal/mol y divide la velocidad entre 10.")
dg_tab = np.r_[np.arange(10, 30, 2.0), dg0]
dg_tab.sort()
k_tab = cin.eyring_rate(dg_tab)
viz.mostrar(fig, viz.datos(
    pd.DataFrame({"ΔG‡": dg_tab, "k": k_tab, "tiempo por reacción": [viz.duracion(1 / k) for k in k_tab]}),
    "la ecuación de Eyring",
    "Cada fila es un punto de la recta: se elige una altura de colina, se calcula k con la ecuación de Eyring a 25 °C "
    "y, dándole la vuelta, cuánto tarda en promedio una reacción.",
    x="ΔG‡", y="k", calculadas={"k": "(k_B·T/h) · exp(−ΔG‡/RT), con k_B·T/h = 6.2 × 10¹² s⁻¹ y RT = 0.59 kcal/mol",
                               "tiempo por reacción": "1 / k"},
    unidades={"ΔG‡": "kcal/mol", "k": "s⁻¹"}, formatos={"ΔG‡": "{:.1f}"},
    resaltar={int(np.argmin(np.abs(dg_tab - dg0))): ("k_cat", "naranja")}))
''')

md(r"""
> ✅ **Para llevar.** Velocidad = intentos × probabilidad de éxito. Los intentos son casi fijos
> (~10¹³ por segundo); lo que decide todo es la altura de la colina, y cada 1.36 kcal/mol es un
> factor 10. Una enzima acelera la reacción **bajando la colina**, no empujando más fuerte.
""")
