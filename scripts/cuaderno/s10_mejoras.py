"""Sección 10. mejoras."""
from .celdas import code, md

# ============================================================================ 10. mejoras

md(r"""
## 10. De ΔE‡ a ΔG‡: tres mejoras del cálculo

> 🎯 **En esta sección** convertirás la barrera electrónica ΔE‡ en la energía libre ΔG‡ que pide
> la ecuación de Eyring, y aprenderás a preguntarte cuánto fiarte de un número calculado.

---

### 💡 La analogía: el topógrafo cuidadoso

Un topógrafo que ha medido la altura de una montaña **una vez, un día, con un instrumento** se
hace tres preguntas antes de publicar el número: ¿el terreno se mueve? ¿otro día habría medido
otra cosa? ¿mi instrumento está bien calibrado?

Hasta aquí tenemos una **energía electrónica** ΔE‡, calculada en **una** conformación de la
enzima con **un** método semiempírico. La ecuación de Eyring pide una **energía libre** ΔG‡. Un
buen científico computacional se hace las mismas tres preguntas, y las responde con cálculos
adicionales:

[[fig:tres_mejoras | Tres preguntas: vibraciones, otras posturas de la enzima y otro método; la suma que lleva de ΔE‡ a ΔG‡]]

---

### 🧮 La ecuación en palabras

**¿Para qué queremos una ecuación?** Porque la altura "a secas" (ΔE‡) ignora que los átomos
vibran siempre. La energía libre suma lo que cuesta, o lo que se ahorra, en vibración y en orden:

$$\text{barrera libre} \;=\; \text{barrera electrónica} \;+\; \text{vibración en reposo} \;+\; \text{vibración térmica} \;-\; T\times\text{(cambio de desorden)}$$

---

### 📐 La ecuación completa, término a término

$$\Delta G^{\ddagger} = \Delta E^{\ddagger} + {\color{#1baf7a}{\Delta \mathrm{ZPE}}} + {\color{#4a3aa7}{\Delta H_{\mathrm{vib}}(T)}} - {\color{#e87ba4}{T\,\Delta S_{\mathrm{vib}}(T)}}$$

| Término | Qué es | En la analogía | Nuestro valor (25 °C) |
|---|---|---|---|
| **ΔE‡** | altura electrónica de la colina, átomos quietos | la medida "de un día" | 19.1 kcal/mol |
| **ΔZPE** | cambio de la **energía de punto cero**: incluso a 0 K los átomos vibran | el terreno tiembla aunque no haya viento | −0.2 kcal/mol |
| **ΔH<sub>vib</sub>(T)** | energía vibracional extra a la temperatura T | el temblor extra de un día caluroso | −0.4 kcal/mol |
| **−TΔS<sub>vib</sub>** | el precio del orden: un TS más rígido que el reactivo tiene ΔS‡ negativa | ordenar un cuarto cuesta trabajo | +1.3 kcal/mol |

* **ΔZPE**: el enlace P–O que se rompe pierde su vibración de tensión en el TS, así que suele
  ser **negativa** (baja la barrera).
* **−TΔS‡**: si el TS es más "rígido" que el reactivo (ΔS‡ < 0), este término **sube** la barrera.

**Las otras dos preguntas no cambian la fórmula, sino la confianza en el número:**

**2. ¿Y si la enzima estuviera en otra postura?** Como el entorno está fijo, cada instantánea
de la MD puede dar una barrera distinta. Repetimos el cálculo en varias instantáneas (250, 500,
750 y 1000 ps) para ver **si la barrera cambia y por qué**. Un buen resultado no es un número:
es un número con su incertidumbre, y a veces la respuesta es que la reacción ni siquiera es
posible desde algunas conformaciones.

**3. ¿Y si el método se equivoca?** PM7 es semiempírico. Recalculamos las energías con otro
hamiltoniano (PM6‑D3H4) en las mismas geometrías: si la barrera cambia poco, el resultado es
robusto; si cambia mucho, sabemos cuánto (des)confiar.

---

### 🎛️ Qué pasa si…

| Si… | entonces… | porque… |
|---|---|---|
| **la temperatura sube de 25 a 37 °C** | ΔG‡ apenas sube (de 19.75 a 19.80 kcal/mol) | con ΔS‡ ≈ −4.3 cal/mol/K, 12 K más solo añaden 12 × 0.0043 ≈ 0.05 kcal/mol |
| **el TS fuera más flexible que el reactivo** (ΔS‡ > 0) | el término entrópico **bajaría** la barrera | −TΔS‡ cambiaría de signo |
| **cambiamos de método** (PM7 → PM6‑D3H4) | la barrera pasa de 19.1 a 4.5 kcal/mol | los métodos semiempíricos se calibran con moléculas concretas |

---

### 🔬 Los datos reales
""")

code(r'''
termo = res["etapas"].get("termoquimica", {})
T_ref = "298.15"
if T_ref in termo:
    t = termo[T_ref]
    print("Termoquímica armónica a 298.15 K (átomos libres de la región QM):")
    print(f"  ΔE‡ = {t['dE_kcal']:.1f}   ΔZPE = {t['dZPE_kcal']:+.1f}   ΔH‡ = {t['dH_kcal']:.1f} kcal/mol   ΔS‡(vib) = {t['dS_vib_cal']:+.1f} cal/mol/K   →   ΔG‡ = {t['dG_kcal']:.1f} kcal/mol")
    for Tk in ("303.15", "310.15"):
        if Tk in termo:
            print(f"  a {float(Tk) - 273.15:.0f} °C: ΔG‡ = {termo[Tk]['dG_kcal']:.1f} kcal/mol")
    print("  ", termo.get("nota", ""))
else:
    print("(termoquímica no disponible en estos datos)")
''')

code(r'''
try:
    inst = datos.json_("qmmm/instantaneas.json")
    filas = [(r["instantanea"], r.get("barrera_kcal", np.nan), r.get("ts_d_PG_O6", np.nan), r.get("ts_d_PG_O3B", np.nan),
              r.get("ts_origen", "dímero" if "barrera_kcal" in r else "falló"))
             for r in inst["instantaneas"]]
    df_inst = pd.DataFrame(filas, columns=["instantánea", "ΔE‡ (kcal/mol)", "d(Pγ–O6) TS (Å)", "d(Pγ–O3β) TS (Å)", "TS obtenido por"])
    display(df_inst.round(2))
    ok = df_inst.dropna(subset=["ΔE‡ (kcal/mol)"])
    media, sd = ok["ΔE‡ (kcal/mol)"].mean(), ok["ΔE‡ (kcal/mol)"].std(ddof=1) if len(ok) > 1 else 0.0
    if len(ok) > 1:
        print(f"Barrera media = {media:.1f} ± {sd:.1f} kcal/mol (n = {len(ok)})")
        fig, ax = viz.figure(6.5, 3.8)
        ax.bar(range(len(ok)), ok["ΔE‡ (kcal/mol)"], color=viz.COLORS["ts"], width=0.55)
        ax.axhline(media, color=viz.INK, lw=1, ls="--")
        ax.set_xticks(range(len(ok))); ax.set_xticklabels(ok["instantánea"], rotation=15)
        ax.set_ylabel("ΔE‡ (kcal/mol)"); ax.set_title("La barrera depende de la conformación de la enzima", loc="left")
        fig;
    else:
        print(f"Solo {len(ok)} conformación con camino completo (barrera {media:.1f} kcal/mol): las demás no llegan a un producto estable (ver abajo).")
    sin_prod = [r for r in inst["instantaneas"] if r.get("sin_producto_estable")]
    if sin_prod:
        print("Instantáneas sin producto estable (la energía sube hasta ξ = 2 Å sin máximo):")
        for r in sin_prod:
            print(f"  {r['instantanea']}: E(ξ = {r['xi_max']:.1f}) = {r['escaneo_max_rel_kcal']:.0f} kcal/mol sobre R; d(Asp205 OD1···H–O6) en el reactivo = {r['d_OD1_H_reactivo']:.2f} Å")
except FileNotFoundError:
    print("(promedio sobre instantáneas no disponible en estos datos)")
''')

md(r"""
**Lo que enseñan las instantáneas de la MD.** En la conformación cristalina minimizada, el
hidroxilo O6–H de la glucosa está unido por puente de hidrógeno a Asp205 (1.7 Å) y el camino
tiene un producto estable. En varias instantáneas de la dinámica clásica, en cambio, Asp205 se
ha alejado del O6–H (en la MD la distancia O6···OD1 promedia 4.3 Å, sección 5): al forzar la
transferencia del fosforilo el protón **no encuentra a la base** y la energía sube sin parar,
sin mínimo de producto. Dos lecciones: (1) la reacción solo puede ocurrir desde las
**conformaciones de ataque cercano** correctas, que son una fracción del tiempo; (2) el campo de
fuerza clásico no siempre mantiene la geometría catalítica (aquí, la interacción
hidroxilo–carboxilato junto al Mg²⁺), y por eso el modelador debe comparar la MD con el cristal
antes de elegir la estructura de partida del QM/MM. Una barrera "promedio" solo tiene sentido
sobre conformaciones reactivas.
""")

code(r'''
met = res["etapas"].get("metodos", {})
for nombre, d in met.items():
    if isinstance(d, dict):
        print(f"{nombre:>9s}//PM7: barrera {d['barrera_kcal']:.1f} kcal/mol, ΔE(reacción) {d['dE_reaccion_kcal']:+.1f} kcal/mol")
print(met.get("nota", ""))
''')

md(r"""
**Lo que dicen estos números.** La corrección vibracional es pequeña (ΔZPE ≈ −0.2 y −TΔS‡ ≈ +1.3
kcal/mol se compensan en parte): ΔG‡ ≈ ΔE‡ + 0.6 kcal/mol. El promedio sobre instantáneas dice
cuánto "respira" la barrera con la conformación de la enzima. Y el cambio de hamiltoniano es la
sorpresa: PM6‑D3H4 da energías muy distintas en las mismas geometrías. Los métodos semiempíricos
se parametrizan frente a conjuntos concretos de moléculas, y el trío fosfato + Mg²⁺ + carboxilato
es un caso difícil para ellos (fósforo hipervalente, cargas altas, un catión divalente). La
lección: **antes de creer una barrera, hay que saber para qué química fue calibrado el método**
y, si es posible, comparar con un cálculo de más nivel (DFT) o con el experimento (sección 12);
aquí PM7 resulta coherente con otro cálculo QM/MM publicado y con el experimento, y PM6‑D3H4 no.
""")

md(r"""
> ✅ **Para llevar.** ΔG‡ = ΔE‡ + ΔZPE + ΔH_vib − TΔS_vib: las vibraciones corrigen poco
> (19.1 → 19.75 kcal/mol, porque ΔZPE y −TΔS‡ casi se compensan). Las otras dos preguntas no
> cambian el número sino la **confianza** en él: la barrera solo tiene sentido desde
> conformaciones reactivas, y un método mal calibrado para fosfato + Mg²⁺ puede dar un valor muy
> distinto.
""")
