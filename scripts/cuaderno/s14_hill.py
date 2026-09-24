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
# @title 🎛️ Sigmoide frente a hipérbola
# 🎛️ Mueve n y S0.5: compara con la hipérbola y mira la ventana 10 %–90 %
interactivo.explorar_hill()
''')

code(r'''
# @title 🩸 La glucoquinasa en el rango de glucosa de la sangre
S_g = np.linspace(0.01, 30, 300)
s_half, n_h = valor("s_half_mm", 7.5), valor("hill_n", 1.7)
fig = viz.plot_hill_vs_mm(S_g, cin.michaelis_menten(S_g, 1.0, s_half), cin.hill(S_g, 1.0, s_half, n_h), n_hill=n_h, s_half=s_half,
                          ylabel="v₀ / V_max", title="En el rango de la sangre, la glucoquinasa responde como un interruptor",
                          subtitle="Entre 4 y 10 mM de glucosa (el rango fisiológico) la sigmoide es mucho más sensible que la hipérbola")
fig.set_size_inches(10.5, 5.8)
ax = fig.axes[0]
viz._kband(ax, 4, 7, "glucosa en sangre\nen ayunas (4–7 mM)", color=viz._tint(viz.COLORS["ts"], 0.9), y_text=0.97)
s10, s90 = cin.substrate_at_fraction(1.0, s_half, n_h, 0.1), cin.substrate_at_fraction(1.0, s_half, n_h, 0.9)
S_tab = np.array(sorted({1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, round(s_half, 2)}), dtype=float)
mm_tab, hill_tab = cin.michaelis_menten(S_tab, 1.0, s_half), cin.hill(S_tab, 1.0, s_half, n_h)
tabla_hill = pd.DataFrame({"[glucosa]": S_tab, "[S]/S₀.₅": S_tab / s_half, "hipérbola (n = 1)": mm_tab,
                           f"sigmoide (n = {n_h:g})": hill_tab, "diferencia": hill_tab - mm_tab})
resalte = {int(np.argmin(np.abs(S_tab - s_half))): ("S₀.₅", "naranja")}
resalte.update({i: ("ayuno", "agua") for i, sv in enumerate(S_tab) if 4 <= sv <= 7 and i not in resalte})
datos_hill = viz.datos(
    tabla_hill, "sigmoide frente a hipérbola",
    "Las dos curvas evaluadas en concentraciones redondas de glucosa, como fracción de V_max. Con poca glucosa la "
    "sigmoide va por debajo (la enzima casi no responde); cerca de S₀.₅ se cruzan; por encima, la sigmoide sube más rápido.",
    x="[glucosa]", y=["hipérbola (n = 1)", f"sigmoide (n = {n_h:g})"],
    calculadas={"[S]/S₀.₅": "[glucosa] ÷ S₀.₅",
                "hipérbola (n = 1)": "[S] / (S₀.₅ + [S])  (Michaelis–Menten con K_M = S₀.₅)",
                f"sigmoide (n = {n_h:g})": "[S]ⁿ / (S₀.₅ⁿ + [S]ⁿ)  (Hill)",
                "diferencia": "sigmoide − hipérbola"},
    unidades={"[glucosa]": "mM", "hipérbola (n = 1)": "v₀/V_max", f"sigmoide (n = {n_h:g})": "v₀/V_max", "diferencia": "v₀/V_max"},
    formatos={"[glucosa]": "{:g}", "[S]/S₀.₅": "{:.2f}", "hipérbola (n = 1)": "{:.0%}", f"sigmoide (n = {n_h:g})": "{:.0%}",
              "diferencia": "{:+.0%}"},
    resaltar=resalte, barra=f"sigmoide (n = {n_h:g})")
viz.mostrar(fig, datos_hill, viz.tarjetas([
    ("S₀.₅", f"{s_half:g}", "mM", "glucosa a media actividad", "naranja"),
    ("coeficiente de Hill", f"{n_h:g}", "", "n > 1: sigmoide", "naranja"),
    ("del 10 % al 90 %", f"{s10:.1f} → {s90:.1f}", "mM", f"[S] × {s90/s10:.0f}", "naranja"),
    ("con n = 1", "× 81", "", "la hipérbola necesitaría mucho más", "azul"),
], titulo="El interruptor de la glucoquinasa en números"))
''')

code(r'''
# @title 📐 Ajuste de Hill y gráfico de Hill
# Ajuste de Hill a datos simulados y gráfico de Hill (la pendiente es n)
S_h = np.array([1, 2, 3, 4, 5, 6, 7.5, 9, 11, 14, 18, 25, 35, 50])
v_h = cin.hill(S_h, 10.0, s_half, n_h) * (1 + 0.03 * rng.standard_normal(S_h.size))
ajh = cin.fit_hill(S_h, v_h)
fig = viz.plot_hill_plot(S_h, v_h, ajh["vmax"], title=f"Los datos suben más empinados que la referencia: n = {ajh['n']:.2f}",
                         subtitle="log[v/(V_max − v)] frente a log[S]: la pendiente es el coeficiente de Hill; la línea discontinua es n = 1")
frac = v_h / (ajh["vmax"] - v_h)
tabla_hp = pd.DataFrame({"[S]": S_h, "v₀": v_h, "log₁₀[S]": np.log10(S_h), "v₀/(V_max − v₀)": frac,
                         "log₁₀[v₀/(V_max − v₀)]": np.log10(frac)})
datos_hp = viz.datos(
    tabla_hp, "el gráfico de Hill",
    "Los 14 tubos del experimento simulado. Para el gráfico de Hill se transforman las dos columnas: el eje x es "
    "log₁₀[S] y el eje y es log₁₀[v₀/(V_max − v₀)], con el V_max del ajuste. Si la enzima sigue la ecuación de Hill, "
    "los puntos caen en una recta de pendiente n; cruzan el 0 del eje y justo en [S] = S₀.₅ (mitad ocupada).",
    x="log₁₀[S]", y="log₁₀[v₀/(V_max − v₀)]",
    calculadas={"log₁₀[S]": "logaritmo decimal de [S]",
                "v₀/(V_max − v₀)": "enzima «encendida» ÷ enzima «apagada»",
                "log₁₀[v₀/(V_max − v₀)]": "logaritmo de la columna anterior (= n·log[S] − n·log S₀.₅)"},
    unidades={"[S]": "mM", "v₀": "µM/s"},
    formatos={"[S]": "{:g}", "v₀": "{:.2f}", "log₁₀[S]": "{:.3f}", "v₀/(V_max − v₀)": "{:.3f}", "log₁₀[v₀/(V_max − v₀)]": "{:+.3f}"},
    resaltar={int(np.argmin(np.abs(np.log10(frac)))): ("≈ S₀.₅", "naranja")})
viz.mostrar(fig, datos_hp, viz.tarjetas([
    ("V_max", f"{ajh['vmax']:.2f}", "", "meseta del ajuste", "azul"),
    ("S₀.₅", f"{ajh['s_half']:.2f}", "mM", "glucosa a media actividad", "naranja"),
    ("n (Hill)", f"{ajh['n']:.2f}", f"± {ajh['n_err']:.2f}", "n > 1: cooperatividad positiva", "naranja"),
], titulo="Ajuste de Hill a datos simulados con 3 % de ruido"))
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
# @title 🧬 Mutantes GCK‑MODY: el interruptor se desplaza
mut = ref.get("mutants", [])
if mut:
    S_g = np.linspace(0.01, 30, 300)
    elegidos = [m for m in mut if m["name"] in ("V244G", "G223S", "I110N", "W99L", "M197I")]
    curvas = [("silvestre", s_half, n_h, 1.0)] + [(m["name"], m["s_half_mm"], m.get("hill_n", n_h), m.get("kcat_rel", 1.0)) for m in elegidos]
    fig, ax = viz.figure(11.5, 6.4)
    viz._kband(ax, 4, 7, color=viz._tint(viz.COLORS["ts"], 0.9))
    etiquetas = []
    for (nombre, s05, nn, krel), color in zip(curvas, [viz.INK] + list(viz.PALETTE[: len(elegidos)])):
        y = krel * cin.hill(S_g, 1.0, s05, nn)
        ax.plot(S_g, y, color=color, lw=3.2 if nombre == "silvestre" else 2.4, zorder=4 if nombre == "silvestre" else 3,
                label=f"{nombre}: S₀.₅ = {s05:g} mM, k_cat ×{krel:.2f}")
        y5 = krel * cin.hill(5.0, 1.0, s05, nn)
        ax.plot([5.0], [y5], "o", ms=8, mfc=color, mec="white", mew=1.6, zorder=6)
        etiquetas.append([y[-1], nombre])
    # etiquetas directas al final de cada curva, separadas para que no se pisen
    etiquetas.sort()
    y_max = max(e[0] for e in etiquetas)
    for j in range(1, len(etiquetas)):
        etiquetas[j][0] = max(etiquetas[j][0], etiquetas[j - 1][0] + 0.06 * y_max)
    for y_lab, nombre in etiquetas:
        ax.annotate(nombre, (S_g[-1], y_lab), xytext=(6, 0), textcoords="offset points", ha="left", va="center",
                    fontsize=10.5, color=viz.INK, fontweight="semibold" if nombre == "silvestre" else "normal",
                    annotation_clip=False)
    viz._guide(ax, "v", 5.0, color=viz.INK_SECONDARY)
    ax.annotate("5 mM (franja: glucosa en ayunas, 4–7 mM)\ncada punto: la actividad con la que\nel páncreas «ve» la glucosa", (5.0, 0.0), xytext=(96, 14),
                textcoords="offset points", ha="left", va="bottom", fontsize=10.5, color=viz.INK_SECONDARY)
    ax.set_xlim(0, 30); ax.set_ylim(bottom=0)
    viz._finish(ax, "[glucosa] (mM)", "actividad relativa a V_max silvestre",
                "Cada mutación desplaza el interruptor de la insulina",
                "Curvas de Hill con los valores publicados: a la derecha o más abajo = el páncreas «ve» menos glucosa (GCK‑MODY)")
    ax.legend(fontsize=10, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.13), frameon=False)
    tabla = pd.DataFrame([(m["name"], m["kind"], m["s_half_mm"], m.get("hill_n"), m.get("kcat_s"),
                           round(m.get("kcat_rel", 1.0) * cin.hill(5.0, 1.0, m["s_half_mm"], m.get("hill_n", n_h)), 3)) for m in mut],
                         columns=["mutante", "fenotipo", "S₀.₅ (mM)", "n", "k_cat (s⁻¹)", "actividad a 5 mM (rel.)"])
    tabla.loc[len(tabla)] = ["silvestre", "-", s_half, n_h, valor("kcat_s"), round(cin.hill(5.0, 1.0, s_half, n_h), 3)]
    viz.mostrar(fig, viz.tabla(tabla, titulo="Los mutantes, uno por uno", nota="Fuente: " + ref.get("mutants_source", "")))
else:
    viz.mostrar(viz.mensaje("Sin datos de mutantes en la tabla de referencia.", tipo="ojo"))
''')

md(r"""
> ✅ **Para llevar.** El coeficiente de Hill *n* mide cuán abrupto es el interruptor: con n ≈ 1.7
> la glucoquinasa pasa de apagada a encendida con un cambio de glucosa ~6 veces menor que una enzima
> michaeliana. Su cooperatividad no viene de varias subunidades sino de una **memoria cinética**, y
> las mutaciones que desplazan S<sub>0.5</sub> desplazan el umbral de la insulina.
""")
