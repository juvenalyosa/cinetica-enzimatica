"""Sección 14. Cooperatividad (Hill)."""
from .celdas import code, md

# ============================================================================ 14. Hill
md(r"""
## 14. Cooperatividad: la glucoquinasa es un sensor

> 🎯 **En esta sección** verás por qué la curva de la glucoquinasa no es una hipérbola sino una
> **sigmoide**, qué mide el coeficiente de Hill *n* y cómo una enzima con **un solo sitio** puede
> comportarse como un interruptor que decide cuándo el páncreas libera insulina.

---

### 💡 La analogía: regulador de luz frente a interruptor

Un **regulador de luz** (hipérbola) sube la intensidad poco a poco desde el principio. Un
**interruptor con umbral** (sigmoide) casi no hace nada hasta cierto punto y entonces se enciende
de golpe. El páncreas necesita un interruptor: por debajo de ~5 mM de glucosa no debe soltar
insulina, y por encima sí. La glucoquinasa es ese interruptor.

[[fig:regulador_interruptor | Comparación de una hipérbola (regulador de luz) y una sigmoide (interruptor): la sigmoide pasa del 10 al 90 por ciento con un cambio de sustrato mucho menor]]

Mira las franjas sombreadas: para pasar del 10 % al 90 % de la actividad, la hipérbola necesita
multiplicar la glucosa **por 81**; la sigmoide de la glucoquinasa, solo **por 13**.

---

### 🧮 La ecuación en palabras

**¿Para qué queremos otra ecuación?** Porque Michaelis–Menten solo sabe dibujar hipérbolas y los
datos de la glucoquinasa no lo son. Necesitamos un número que diga **cuán abrupto** es el
interruptor.

La idea es la misma de la sección 13 (velocidad = máxima × fracción ocupada), pero la
concentración de sustrato entra **elevada a una potencia n**:

$$\text{velocidad} \;=\; V_{\max} \times \frac{(\text{sustrato})^{\,n}}{(\text{umbral})^{\,n} + (\text{sustrato})^{\,n}}$$

Elevar a *n* > 1 hace que las concentraciones **por debajo** del umbral cuenten todavía menos y las
de **por encima** cuenten todavía más: la curva se aplana al principio y se empina en el centro.
Es el interruptor.

---

### 📐 La ecuación completa, término a término (Hill)

$$v_0 = \frac{V_{\max}[\mathrm{S}]^{\color{#eb6834}{n}}}{{\color{#2a78d6}{S_{0.5}}}^{\,\color{#eb6834}{n}} + [\mathrm{S}]^{\color{#eb6834}{n}}}$$

| Término | Qué es | En la analogía | Glucoquinasa |
|---|---|---|---|
| $V_{\max}$ | velocidad máxima | la luz a tope | k<sub>cat</sub>·[E]<sub>0</sub> |
| ${\color{#2a78d6}{S_{0.5}}}$ | concentración a la que se alcanza la mitad de V<sub>max</sub> (el análogo de K<sub>M</sub>) | dónde está el umbral del interruptor | ≈ 7.7 mM de glucosa |
| ${\color{#eb6834}{n}}$ | **coeficiente de Hill**: cuán abrupto es el cambio | cuán "de golpe" se enciende la luz | ≈ 1.7 |

$n = 1$ es Michaelis–Menten; $n > 1$ es cooperatividad positiva.

---

### 🎛️ Qué pasa si…

| Si… | para pasar del 10 % al 90 % de la actividad hay que multiplicar [S] por… | ejemplo |
|---|---|---|
| n = 1 | **81** | una enzima michaeliana |
| n = 1.7 | **~13** | la glucoquinasa |
| n = 4 | **3** | la hemoglobina |

---

### 💡 Lo sorprendente: una enzima con memoria

La glucoquinasa es un **monómero** con un solo sitio para glucosa: no puede haber "comunicación
entre subunidades" como en la hemoglobina. Su cooperatividad es **cinética**: la enzima cambia
entre una forma poco activa (abierta) y una activa (cerrada) a una velocidad **comparable a la de
la catálisis** (k<sub>ex</sub> ≈ 5–100 s⁻¹ frente a k<sub>cat</sub> ≈ 60 s⁻¹, Larion et al. 2012).

[[fig:cooperatividad_cinetica | La glucoquinasa alterna entre una forma abierta lenta y una cerrada rápida; con mucha glucosa no le da tiempo a relajarse y se mantiene en la forma rápida]]

Con mucha glucosa la enzima no tiene tiempo de relajarse a la forma lenta entre un ciclo y el
siguiente, y trabaja más de lo que "debería". Es el modelo **mnemónico** o de **transición lenta**
(Storer y Cornish‑Bowden 1977; Neet y Ainslie; Cárdenas 1984).

> ✍️ **Para pensar.** Un inhibidor competitivo como la N‑acetilglucosamina *suprime* la
> cooperatividad (n → 1), porque mantiene ocupado el sitio y rompe la "memoria". Las simulaciones
> de MD de la sección 5 (dominios que se abren y cierran) son la imagen molecular de ese cambio.

---

### 🔬 Los datos reales

Mueve *n* y S<sub>0.5</sub>, y compara con la hipérbola:
""")

code(r'''
# 🎛️ Mueve n y S0.5: compara con la hipérbola y mira la ventana 10 %–90 %
interactivo.explorar_hill()
''')

code(r'''
S_g = np.linspace(0.01, 30, 300)
s_half, n_h = valor("s_half_mm", 7.5), valor("hill_n", 1.7)
fig = viz.plot_hill_vs_mm(S_g, cin.michaelis_menten(S_g, 1.0, s_half), cin.hill(S_g, 1.0, s_half, n_h), n_hill=n_h, s_half=s_half,
                          ylabel="v₀ / V_max", title="Glucoquinasa: sigmoide (Hill) frente a hipérbola",
                          subtitle="Entre 4 y 10 mM de glucosa (el rango fisiológico) la sigmoide es mucho más sensible")
ax = fig.axes[0]; ax.axvspan(4, 7, color=viz.GRID, alpha=0.6, zorder=0)
ax.annotate("glucosa en sangre\nen ayunas (4–7 mM)", (5.5, 0.05), ha="center", fontsize=9, color=viz.INK_SECONDARY)
fig;
s10, s90 = cin.substrate_at_fraction(1.0, s_half, n_h, 0.1), cin.substrate_at_fraction(1.0, s_half, n_h, 0.9)
print(f"Con n = {n_h}: del 10 % al 90 % de V_max entre {s10:.1f} y {s90:.1f} mM (factor {s90/s10:.0f}); con n = 1 haría falta un factor 81.")
''')

code(r'''
# Ajuste de Hill a datos simulados y gráfico de Hill (la pendiente es n)
S_h = np.array([1, 2, 3, 4, 5, 6, 7.5, 9, 11, 14, 18, 25, 35, 50])
v_h = cin.hill(S_h, 10.0, s_half, n_h) * (1 + 0.03 * rng.standard_normal(S_h.size))
ajh = cin.fit_hill(S_h, v_h)
print(f"Ajuste de Hill: V_max = {ajh['vmax']:.2f}, S_0.5 = {ajh['s_half']:.2f} mM, n = {ajh['n']:.2f} ± {ajh['n_err']:.2f}")
fig = viz.plot_hill_plot(S_h, v_h, ajh["vmax"], subtitle="log[v/(V_max − v)] frente a log[S]: la pendiente es el coeficiente de Hill")
fig;
''')

md(r"""
### 🔬 Mutaciones que mueven el interruptor: GCK‑MODY e hiperinsulinismo

Una mutación que sube S<sub>0.5</sub> (o baja k<sub>cat</sub>) hace que el páncreas "vea" menos
glucosa de la que hay y libere insulina tarde: glucemia alta desde el nacimiento, pero estable
(GCK‑MODY). Una mutación activadora baja S<sub>0.5</sub>: insulina de más, hipoglucemia congénita.
En la analogía: el interruptor se desplaza a la derecha o a la izquierda. Con los valores publicados
(Valentínová 2012; Sayed 2009) dibujamos qué hace cada mutante a 5 mM de glucosa:
""")

code(r'''
mut = ref.get("mutants", [])
if mut:
    S_g = np.linspace(0.01, 30, 300)
    elegidos = [m for m in mut if m["name"] in ("V244G", "G223S", "I110N", "W99L", "M197I")]
    curvas = [("silvestre", s_half, n_h, 1.0)] + [(m["name"], m["s_half_mm"], m.get("hill_n", n_h), m.get("kcat_rel", 1.0)) for m in elegidos]
    fig, ax = viz.figure(8, 5.2)
    for (nombre, s05, nn, krel), color in zip(curvas, [viz.INK] + list(viz.PALETTE[: len(elegidos)])):
        ax.plot(S_g, krel * cin.hill(S_g, 1.0, s05, nn), color=color, lw=2, label=f"{nombre}: S₀.₅ = {s05:g} mM, k_cat ×{krel:.2f}")
    ax.axvline(5.0, color=viz.INK_MUTED, ls="--", lw=1)
    ax.annotate("5 mM de glucosa", (5.0, 0.02), xytext=(4, 0), textcoords="offset points", fontsize=9, color=viz.INK_SECONDARY)
    ax.set_xlabel("[glucosa] (mM)"); ax.set_ylabel("actividad relativa a V_max silvestre")
    ax.legend(fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.18), frameon=False)
    ax.set_title("Mutantes de la glucoquinasa: el interruptor se desplaza", loc="left")
    fig;
    tabla = pd.DataFrame([(m["name"], m["kind"], m["s_half_mm"], m.get("hill_n"), m.get("kcat_s"),
                           round(m.get("kcat_rel", 1.0) * cin.hill(5.0, 1.0, m["s_half_mm"], m.get("hill_n", n_h)), 3)) for m in mut],
                         columns=["mutante", "fenotipo", "S₀.₅ (mM)", "n", "k_cat (s⁻¹)", "actividad a 5 mM (rel.)"])
    tabla.loc[len(tabla)] = ["silvestre", "-", s_half, n_h, valor("kcat_s"), round(cin.hill(5.0, 1.0, s_half, n_h), 3)]
    display(tabla)
    print("Fuente:", ref.get("mutants_source", ""))
else:
    print("(sin datos de mutantes en la tabla de referencia)")
''')

md(r"""
> ✅ **Para llevar.** El coeficiente de Hill *n* mide cuán abrupto es el interruptor: con n ≈ 1.7
> la glucoquinasa pasa de apagada a encendida con un cambio de glucosa ~6 veces menor que una enzima
> michaeliana. Su cooperatividad no viene de varias subunidades sino de una **memoria cinética**, y
> las mutaciones que desplazan S<sub>0.5</sub> desplazan el umbral de la insulina.
""")
