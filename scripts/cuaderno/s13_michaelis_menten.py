"""Sección 13. Michaelis–Menten desde el mecanismo."""
from .celdas import code, md

# ============================================================================ 13. MM
md(r"""
## 13. Michaelis–Menten desde el mecanismo

> 🎯 **En esta sección** verás por qué la velocidad de una enzima **se satura** cuando hay mucho
> sustrato, de dónde sale la ecuación de Michaelis–Menten (sin creerla: la derivaremos y la
> simularemos) y cómo se obtienen K<sub>M</sub>, V<sub>max</sub>, k<sub>cat</sub> y
> k<sub>cat</sub>/K<sub>M</sub> a partir de datos.

---

### 💡 La analogía: el cajero del supermercado

Una caja de supermercado (la **enzima**) atiende clientes (el **sustrato**). Con pocos clientes,
cada uno que llega es atendido enseguida: la "velocidad" de la caja crece con el número de
clientes. Con muchos clientes se forma cola: la caja trabaja a tope y da igual que lleguen más.
La velocidad se **satura**.

[[fig:cajero_saturacion | Tres viñetas de una caja de supermercado con pocos, algunos y muchos clientes, y los tres tramos correspondientes de la curva de saturación]]

Ese comportamiento tiene forma de hipérbola y se describe con dos números:

* **V**<sub>max</sub>: la velocidad con la caja siempre ocupada (toda la enzima en forma ES).
* **K**<sub>M</sub>: la concentración de clientes a la que la caja trabaja a la **mitad** de su
  máximo. Una K<sub>M</sub> pequeña significa que la enzima "se llena" con poco sustrato.

---

### 🧮 La ecuación en palabras

**¿Para qué queremos una ecuación?** Para predecir cuánto producto sale por segundo con cualquier
cantidad de sustrato, y para resumir una enzima entera en dos números que se pueden comparar entre
laboratorios, mutantes y fármacos.

La idea es sencilla: la velocidad es la **velocidad máxima** multiplicada por la **fracción del
tiempo que la caja está ocupada**:

$$\text{velocidad} \;=\; \underbrace{\text{velocidad con la caja siempre ocupada}}_{V_{\max}} \;\times\; \underbrace{\text{fracción de tiempo ocupada}}_{\text{de 0 a 1}}$$

[[fig:mm_anatomia | La ecuación de Michaelis–Menten como producto de la velocidad máxima por la fracción de enzima ocupada]]

La fracción ocupada, $[\mathrm{S}]/(K_M + [\mathrm{S}])$, vale 0 sin sustrato, 1/2 cuando
$[\mathrm{S}] = K_M$ y se acerca a 1 (sin llegar nunca) cuando hay muchísimo. Eso es todo.

---

### 📐 La ecuación completa, término a término

**El mecanismo que hay detrás (Michaelis y Menten 1913; Briggs y Haldane 1925):**

$$\mathrm{E + S \;\underset{k_{-1}}{\overset{k_1}{\rightleftharpoons}}\; ES \;\overset{k_2}{\longrightarrow}\; E + P}$$

[[fig:mecanismo_mm | Esquema del mecanismo: la enzima une el sustrato (k1), lo puede soltar (k−1) o transformarlo en producto (k2)]]

| Constante | Qué describe | En la analogía | Unidades |
|---|---|---|---|
| $k_1$ | el sustrato entra al sitio activo (choque productivo) | el cliente llega a la caja | M⁻¹ s⁻¹ |
| $k_{-1}$ | el sustrato se suelta sin reaccionar | se va sin comprar | s⁻¹ |
| $k_2$ (= k<sub>cat</sub>) | el complejo ES cruza la barrera y suelta el producto | paga y sale con su compra | s⁻¹ |

**La derivación, paso a paso.** (Sigue cada línea: solo es álgebra.)

1. Velocidad de formación de producto: $v_0 = k_2[\mathrm{ES}]$ — *cuántos clientes pagan por
   segundo = rapidez del cajero × clientes siendo atendidos*.
2. **Hipótesis del estado estacionario**: ES se forma tan rápido como se consume, así que su
   concentración casi no cambia: $k_1[\mathrm{E}][\mathrm{S}] = (k_{-1} + k_2)[\mathrm{ES}]$ — *llegan
   a la caja tantos clientes como se van (comprando o no)*.
3. La enzima total se reparte: $[\mathrm{E}]_0 = [\mathrm{E}] + [\mathrm{ES}]$, luego $[\mathrm{E}] = [\mathrm{E}]_0 - [\mathrm{ES}]$ —
   *cada caja está libre u ocupada*.
4. Sustituyendo (3) en (2) y despejando: $[\mathrm{ES}] = \dfrac{[\mathrm{E}]_0[\mathrm{S}]}{K_M + [\mathrm{S}]}$ con
   $K_M \equiv \dfrac{k_{-1}+k_2}{k_1}$.
5. Metiendo (4) en (1):

$$\boxed{v_0 \;=\; {\color{#2a78d6}{V_{\max}}}\;\times\;{\color{#1baf7a}{\frac{[\mathrm{S}]}{K_M + [\mathrm{S}]}}}} \qquad\text{con}\qquad {\color{#2a78d6}{V_{\max} = k_2[\mathrm{E}]_0}}$$

El factor **azul** es la velocidad máxima; el **verde**, la fracción de enzima ocupada (los mismos
colores que en el dibujo).

| Término | Qué es | En la analogía |
|---|---|---|
| $v_0$ | velocidad inicial de formación de producto | clientes que salen con su compra por segundo |
| ${\color{#2a78d6}{V_{\max}}}$ | velocidad con toda la enzima ocupada: k<sub>cat</sub>·[E]<sub>0</sub> | todas las cajas atendiendo sin pausa |
| [S] | concentración de sustrato | cuántos clientes hay en la tienda |
| $K_M$ | [S] a la que la enzima está medio ocupada: (k<sub>−1</sub> + k<sub>2</sub>)/k<sub>1</sub> | cuántos clientes hacen falta para que la caja esté ocupada la mitad del tiempo |

> ⚠️ **Ojo.** $K_M$ solo es igual a la constante de disociación $K_d = k_{-1}/k_1$ cuando
> $k_2 \ll k_{-1}$ (el cliente se va sin comprar muchas más veces de las que compra). En general
> K<sub>M</sub> no es una "afinidad" pura.

---

### 🎛️ Qué pasa si…

| Si… | entonces… | porque… |
|---|---|---|
| [S] ≪ K<sub>M</sub> | v₀ ≈ (V<sub>max</sub>/K<sub>M</sub>)·[S]: **una recta** | la enzima está casi vacía y cada molécula de sustrato cuenta; pendiente = (k<sub>cat</sub>/K<sub>M</sub>)·[E]<sub>0</sub> |
| [S] = K<sub>M</sub> | v₀ = V<sub>max</sub>/2 | la caja está ocupada la mitad del tiempo |
| [S] ≫ K<sub>M</sub> | v₀ → V<sub>max</sub>: **saturación** | siempre hay cola: más clientes no aceleran al cajero |
| **K<sub>M</sub> sube** | la curva se estira a la derecha | la enzima "agarra" peor el sustrato o lo suelta más rápido: hace falta más para la misma velocidad |
| **V<sub>max</sub> sube** | toda la curva escala hacia arriba | más enzima, o una enzima más rápida |

---

### 🔬 Simulemos el mecanismo

En lugar de creer la fórmula, **resolvamos el mecanismo numéricamente** y veamos cómo la
hipérbola aparece sola. Primero la película de las concentraciones:
""")

code(r'''
# @title 🎛️ El mecanismo en acción: E + S ⇌ ES → E + P
# 🎛️ Mueve k1, k-1, k2, [E]0 y [S]0 y observa el pre-estado estacionario y la formación de producto
interactivo.explorar_mecanismo()
''')

code(r'''
# @title 🔁 De la simulación a la hipérbola
k2 = 60  # @param {type:"slider", min:6, max:600, step:6}
# k2 = k_cat (s⁻¹); prueba 6 y 600 (ejercicio 3)
E0, k1, k_1 = 0.05, 1.0, 50.0                   # µM, µM⁻¹s⁻¹, s⁻¹
S = np.array([5, 10, 20, 40, 80, 150, 300, 600, 1200, 2400])       # µM
v0_sim = cin.initial_rates_from_simulation(E0, S, k1, k_1, k2, t_window=(0.01, 0.1))
ajuste = cin.fit_michaelis_menten(S, v0_sim)
teoria = cin.steady_state_parameters(k1, k_1, k2)
fig = viz.plot_initial_rates_from_ode(S, v0_sim, fit=ajuste)
viz.mostrar(fig, viz.tarjetas([
    ("K_M del ajuste", f"{ajuste['km']:.1f}", "µM", "ajustando la hipérbola a las v₀ simuladas", "azul"),
    ("K_M de la fórmula", f"{teoria['km']:.1f}", "µM", "(k₋₁ + k₂)/k₁", "agua"),
    ("V_max del ajuste", f"{ajuste['vmax']:.3f}", "µM/s", "meseta de la curva", "azul"),
    ("V_max de la fórmula", f"{k2 * E0:.3f}", "µM/s", "k₂·[E]₀", "agua"),
], titulo="¿Coinciden la simulación y la derivación en papel?",
   nota="Si las parejas coinciden, la derivación funciona: la hipérbola sale sola del mecanismo."))
''')

code(r'''
# @title 🎛️ La hipérbola: mueve V_max y K_M
# 🎛️ La hipérbola: mueve Vmax y Km
interactivo.explorar_michaelis_menten()
''')

md(r"""
### 🔬 Los datos reales: tres números que salen de la curva

| Constante | Cómo se obtiene | Qué mide | En la analogía | Glucoquinasa (literatura) |
|---|---|---|---|---|
| **K<sub>M</sub>** (o S<sub>0.5</sub>) | ajuste de v₀ frente a [S] | afinidad aparente; a qué [S] la enzima está medio llena | clientes necesarios para media ocupación | ≈ 7.7 mM de glucosa |
| **k<sub>cat</sub>** = V<sub>max</sub>/[E]<sub>0</sub> | V<sub>max</sub> del ajuste ÷ concentración de enzima | recambio: reacciones por segundo por molécula de enzima | clientes por segundo que despacha **una** caja a tope | ≈ 62–66 s⁻¹ |
| **k<sub>cat</sub>/K<sub>M</sub>** | cociente de los dos anteriores | eficiencia a baja [S]; tope físico ≈ 10⁸–10⁹ M⁻¹s⁻¹ (difusión) | lo bien que la caja aprovecha una tienda casi vacía | ≈ 8 × 10³ M⁻¹s⁻¹ |

Ahora un "experimento": datos simulados con los parámetros de la literatura y un 4 % de ruido,
tratados como lo haría un bioquímico: **ajuste no lineal** y, para comparar, las tres
linealizaciones clásicas.

> ⚠️ **Ojo.** **Lineweaver–Burk** (1/v₀ frente a 1/[S]) amplifica el error de los puntos a baja
> concentración: hoy se usa para *visualizar*, y el ajuste no lineal para *cuantificar*.
""")

code(r'''
# @title 🔬 Un experimento simulado y su ajuste
rng = np.random.default_rng(7)
S_mM = np.array([0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 40, 60])
E0_uM = 0.2                                     # enzima en el tubo (µM)
Km_real = valor("s_half_mm", 7.5)               # tratamos la enzima como hiperbólica en este ejercicio
Vmax_real = valor("kcat_s", 60.0) * E0_uM       # µM/s
v_obs = cin.michaelis_menten(S_mM, Vmax_real, Km_real) * (1 + 0.04 * rng.standard_normal(S_mM.size))
aj = cin.fit_michaelis_menten(S_mM, v_obs)
fig = viz.plot_michaelis_menten(S_mM, v_obs, fit=aj, title=f"La enzima llega a la mitad de su máximo con {aj['km']:.1f} mM de sustrato",
                                subtitle="Datos simulados con ruido del 4 % (puntos) y ajuste no lineal de Michaelis–Menten (línea)")
kc = cin.kcat_km_from_fit(aj, E0_uM)
viz.mostrar(fig, viz.tarjetas([
    ("V_max", f"{aj['vmax']:.2f}", f"± {aj['vmax_err']:.2f} µM/s", "la meseta: toda la enzima ocupada", "azul"),
    ("K_M", f"{aj['km']:.2f}", f"± {aj['km_err']:.2f} mM", f"[S] a media velocidad (R² = {aj['r2']:.4f})", "naranja"),
    ("k_cat = V_max/[E]₀", f"{kc['kcat']:.0f}", "s⁻¹", f"una reacción cada {cin.turnover_time(kc['kcat'])*1000:.0f} ms por enzima", "agua"),
    ("k_cat/K_M", f"{kc['kcat_over_km_molar']:.2e}", "M⁻¹s⁻¹", "eficiencia cuando hay poco sustrato", "violeta"),
], titulo="Los cuatro números que salen de la curva"))
''')

code(r'''
# @title 📐 Las tres linealizaciones clásicas
fig = viz.plot_linearizations(S_mM, v_obs, title="Tres maneras de convertir la hipérbola en una recta",
                              subtitle="Los cortes con los ejes (puntos grises) dan V_max y K_M; círculos naranjas: los puntos que más distorsionan Lineweaver–Burk")
lb = cin.linear_fit(*cin.lineweaver_burk(S_mM, v_obs))
viz.mostrar(fig, viz.tarjetas([
    ("V_max por Lineweaver–Burk", f"{1/lb['intercept']:.2f}", "µM/s", f"ajuste no lineal: {aj['vmax']:.2f}", "gris"),
    ("K_M por Lineweaver–Burk", f"{lb['slope']/lb['intercept']:.2f}", "mM", f"ajuste no lineal: {aj['km']:.2f}", "gris"),
], titulo="Lineweaver–Burk frente al ajuste no lineal",
   nota="Las rectas sirven para visualizar el tipo de comportamiento; para cuantificar se usa el ajuste no lineal."))
''')

md(r"""
> ✅ **Para llevar.** Velocidad = velocidad máxima × fracción de enzima ocupada. La hipérbola
> **emerge** del mecanismo E + S ⇌ ES → E + P; K<sub>M</sub> es la [S] de media ocupación,
> k<sub>cat</sub> las reacciones por segundo de cada enzima y k<sub>cat</sub>/K<sub>M</sub> su
> eficiencia cuando el sustrato escasea. Se miden con un ajuste no lineal.
""")
