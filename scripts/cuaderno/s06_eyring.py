"""Sección 6. Eyring."""
from .celdas import code, md

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

md(r"""
> ✅ **Para llevar.** Velocidad = intentos × probabilidad de éxito. Los intentos son casi fijos
> (~10¹³ por segundo); lo que decide todo es la altura de la colina, y cada 1.36 kcal/mol es un
> factor 10. Una enzima acelera la reacción **bajando la colina**, no empujando más fuerte.
""")
