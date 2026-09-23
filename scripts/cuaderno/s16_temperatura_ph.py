"""Sección 16. Temperatura y pH."""
from .celdas import code, md

# ============================================================================ 16. T y pH
md(r"""
## 16. Temperatura y pH

> 🎯 **En esta sección** verás por qué calentar acelera las reacciones (y cuánto), cómo separar la
> barrera en **entalpía** y **entropía** con las ecuaciones de Arrhenius y Eyring, y por qué la
> actividad de una enzima frente al pH tiene forma de **campana**.

---

### 💡 La analogía: calentar = sacudidas más fuertes

Vuelve a la pelota de la sección 6. La temperatura es la fuerza de las sacudidas: al calentar,
todas las moléculas se mueven más y **crece la fracción** que tiene energía suficiente para
cruzar la colina. Como esa fracción es la "cola" de una distribución, un cambio pequeño de
temperatura la agranda mucho más de lo que uno esperaría.

[[fig:temperatura_sacudidas | Distribución de energías de las moléculas a 25 y 37 °C; la cola que supera la barrera es mayor a 37 °C]]

(Hasta que la enzima se desnaturaliza: calentar demasiado desarma la proteína.)

---

### 🧮 La ecuación en palabras

**¿Para qué queremos una ecuación?** Para leer la **altura de la barrera** en un experimento
sencillo: medir k a varias temperaturas. Cuánto sube la velocidad al calentar dice cuán alta es
la colina.

$$\text{ln(velocidad)} \;=\; \text{constante} \;-\; \frac{\text{altura de la barrera}}{\text{energía de las sacudidas}}$$

Si dibujas ln k frente a 1/T sale una **recta**, y su pendiente es la barrera: cuanto más empinada,
más alta la colina y más sensible la reacción a la temperatura.

---

### 📐 Las ecuaciones completas, término a término

**Arrhenius**: $\ln k = \ln A - {\color{#eb6834}{E_a}}/RT$. En un gráfico de ln k frente a 1/T la
pendiente es $-{\color{#eb6834}{E_a}}/R$: la **energía de activación** es "cuánto sube la velocidad
al calentar".

**Eyring**: $\ln(k/T) = \ln(k_B/h) + {\color{#4a3aa7}{\Delta S^{\ddagger}}}/R - {\color{#eb6834}{\Delta H^{\ddagger}}}/RT$.
Separa la barrera en **entalpía** (energía que hay que aportar) y **entropía** (orden que hay que
imponer: un ΔS‡ negativo significa que el estado de transición es más ordenado que los reactivos).

| Término | Qué es | En la analogía |
|---|---|---|
| ${\color{#eb6834}{E_a,\ \Delta H^{\ddagger}}}$ | energía de activación / entalpía de activación | la altura de la colina |
| ${\color{#4a3aa7}{\Delta S^{\ddagger}}}$ | entropía de activación | lo estrecho del paso: si hay que "enhebrar" la cima con precisión, ΔS‡ < 0 |
| $A$ | factor preexponencial | los intentos por segundo |
| $RT$ | energía térmica | la fuerza típica de las sacudidas |

---

### 🎛️ Qué pasa si…

| Si… | al subir de 25 a 37 °C, k se multiplica por… | porque… |
|---|---|---|
| E<sub>a</sub> = 12 kcal/mol | **~2.2** | la cola sobre una colina alta crece mucho con las sacudidas |
| E<sub>a</sub> = 6 kcal/mol | **~1.5** | sobre una colina baja ya pasan muchas: calentar añade menos |

Cuanto mayor la barrera, más sensible a la temperatura.
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

md(r"""
---

### 💡 La analogía del pH: dos grupos, dos condiciones

Cada enzima tiene un pH óptimo: por encima o por debajo, los grupos que hacen la química pierden
la carga que necesitan. En la glucoquinasa, **Asp205** (la base catalítica) debe estar
**desprotonado** para poder aceptar el protón de la glucosa, y **Lys169** debe estar **protonada**
para estabilizar el fosfato. Es como una puerta con dos cerraduras: solo se abre cuando las dos
están en la posición correcta.

[[fig:ph_campana | Curva de actividad frente a pH en forma de campana, con los estados de protonación de Asp205 y Lys169 a pH bajo, óptimo y alto]]

### 🧮 La ecuación en palabras

$$\text{actividad} \;=\; \frac{\text{actividad máxima}}{1 + (\text{penalización si falta la base}) + (\text{penalización si falta el ácido})}$$

Cada penalización es casi cero en el pH correcto y crece 10 veces por cada unidad de pH que te
alejas: de ahí la campana.

### 📐 La ecuación completa, término a término

$$v = \frac{v_{\max}}{1 + {\color{#eb6834}{10^{\,pK_1 - \mathrm{pH}}}} + {\color{#4a3aa7}{10^{\,\mathrm{pH} - pK_2}}}}$$

| Término | Qué es | Cuándo pesa |
|---|---|---|
| ${\color{#eb6834}{10^{\,pK_1 - \mathrm{pH}}}}$ | cociente Asp205 **protonado** (inútil como base) / desprotonado | a pH bajo: "apaga" la enzima |
| ${\color{#4a3aa7}{10^{\,\mathrm{pH} - pK_2}}}$ | cociente Lys169 **desprotonada** (sin carga positiva) / protonada | a pH alto: "apaga" la enzima |
| (pK<sub>1</sub> + pK<sub>2</sub>)/2 | posición del máximo | el pH óptimo |

> 🔬 **Los datos reales.** Para la glucoquinasa humana el óptimo medido es pH 8.5–8.7 (Šimčíková y
> Heneberg 2019); durante décadas se creyó más bajo porque el ATP acidifica los tampones del ensayo:
> un ejemplo de cómo un detalle técnico cambia un "hecho" de libro de texto.
""")

code(r'''
# 🎛️ pKa de la base y del ácido: la campana de pH
interactivo.explorar_ph()
''')

md(r"""
> ✅ **Para llevar.** Calentar agranda la cola de moléculas que superan la barrera: cuanto más alta
> la colina, más gana la reacción con la temperatura (Arrhenius mide E<sub>a</sub>; Eyring la separa
> en ΔH‡ y ΔS‡). El pH decide si los grupos catalíticos tienen la carga correcta: dos pK<sub>a</sub>
> dibujan una campana con el máximo entre ellos.
""")
