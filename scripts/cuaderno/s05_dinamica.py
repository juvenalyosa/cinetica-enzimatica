"""Sección 5. MD."""
from .celdas import code, md

# ============================================================================ 5. MD
md(r"""
## 5. La enzima se mueve: dinámica molecular

> 🎯 **En esta sección** verás cómo se "filma" una proteína con las leyes de Newton y para qué
> sirve esa película: medir cuánto tiempo pasan los reactivos bien colocados para reaccionar.

---

### 💡 La analogía: la foto y la película

La estructura cristalina es una **foto**: una sola pose. Pero una proteína a 37 °C nunca está
quieta: se dobla, respira y vibra. La **dinámica molecular (MD)** es la **película**: una
sucesión de fotogramas separados por 2 femtosegundos (0.000000000000002 s). Un nanosegundo de
película son 500 000 fotogramas.

[[fig:pelicula_md | Una foto del cristal frente a una tira de película de dinámica molecular con fotogramas cada 2 fs]]

---

### 🧮 La ecuación en palabras

**¿Para qué queremos una ecuación?** Para saber dónde estará cada átomo en el fotograma
siguiente. La receta es la de una mesa de billar: si sabes con qué fuerza empujan a una bola,
sabes cómo se acelera y dónde estará un instante después.

Cada átomo es una **bolita** con masa, unida a sus vecinas por **resortes** (los enlaces), y
además siente las **cargas** y los **choques** de las demás. Con eso, cada fotograma es una
vuelta del mismo ciclo:

$$\text{fuerzas} \;\rightarrow\; \text{aceleración} \;\rightarrow\; \text{posición nueva} \;\rightarrow\; \text{fuerzas} \;\rightarrow\; \cdots$$

[[fig:newton_ciclo | Modelo de bolitas y resortes y el ciclo de Newton: posiciones, fuerzas, aceleración y avance de 2 fs]]

---

### 📐 La ecuación completa, término a término

La segunda ley de Newton para cada átomo *i*, y un paso de tiempo:

$${\color{#eb6834}{\mathbf{F}_i}} = m_i\,{\color{#4a3aa7}{\mathbf{a}_i}}
\qquad\Rightarrow\qquad
\mathbf{v}_i \leftarrow \mathbf{v}_i + {\color{#4a3aa7}{\mathbf{a}_i}}\,{\color{#1baf7a}{\Delta t}},
\qquad
{\color{#2a78d6}{\mathbf{r}_i}} \leftarrow {\color{#2a78d6}{\mathbf{r}_i}} + \mathbf{v}_i\,{\color{#1baf7a}{\Delta t}}$$

| Término | Qué es | En la analogía | Valor aquí |
|---|---|---|---|
| ${\color{#2a78d6}{\mathbf{r}_i}}$ | posición del átomo *i* | dónde está la bolita | ~57 000 átomos |
| ${\color{#eb6834}{\mathbf{F}_i}}$ | fuerza total sobre el átomo | tirones de los resortes + atracción/repulsión de las cargas + choques | campo de fuerza Amber ff14SB |
| $m_i$ | masa del átomo | lo pesada que es la bolita | H = 1, C = 12, O = 16 g/mol |
| ${\color{#4a3aa7}{\mathbf{a}_i}}$ | aceleración | cuánto cambia su velocidad | $\mathbf{F}_i/m_i$ |
| $\mathbf{v}_i$ | velocidad | hacia dónde y qué tan rápido va | ~ la que da la temperatura (300 K) |
| ${\color{#1baf7a}{\Delta t}}$ | paso de tiempo | tiempo entre fotogramas | 2 fs |

Las fuerzas salen de una **energía** $U$ que suma resortes, cargas y choques: cada átomo es
empujado "cuesta abajo" en esa energía, $\mathbf{F}_i = -\partial U/\partial \mathbf{r}_i$. El
programa (OpenMM) usa un integrador de Langevin, una versión más cuidadosa de este mismo ciclo que
además mantiene la temperatura en 300 K.

---

### 🎛️ Qué pasa si…

| Si… | entonces… | porque… |
|---|---|---|
| **Δt fuera mucho mayor** | la simulación "explota" | los átomos saltarían demasiado entre fotogramas y los resortes se estirarían sin control; 2 fs es posible porque se fijan las longitudes de los enlaces con H |
| **simulamos más tiempo** | vemos movimientos más lentos (dominios que se abren y cierran) | 1 ns son 500 000 vueltas del ciclo; los cambios grandes pueden tardar mucho más |
| **quisiéramos romper un enlace** | la MD clásica no puede | sus resortes no saben de electrones: para eso hace falta QM/MM (sección 7) |

---

### 🔬 Los datos reales

**Qué se hizo.** Con OpenMM y Amber: minimización, calentamiento hasta 300 K, equilibración a
1 bar y 1 ns de producción (script `02_dinamica_molecular.py`). Aquí no hay química (los enlaces
no se rompen): la MD sirve para ver **cuánto tiempo pasa el sitio activo en una geometría
reactiva**, con el Pγ del ATP cerca del O6 de la glucosa. Esas geometrías se llaman
**conformaciones de ataque cercano** (*near‑attack conformations*, NAC).

[[fig:ataque_cercano | Conformación de ataque cercano con el O6 a menos de 3.5 Å del fósforo frente a una conformación no reactiva]]

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

md(r"""
> ✅ **Para llevar.** La MD es una película hecha con las leyes de Newton, un fotograma cada 2 fs.
> No rompe enlaces, pero dice algo clave: la enzima **mantiene a los reactivos apuntándose** la
> mayor parte del tiempo, y solo desde esas conformaciones de ataque cercano puede ocurrir la química.
""")
