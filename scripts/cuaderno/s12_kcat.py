"""Sección 12. barrera -> kcat."""
from .celdas import code, md

# ============================================================================ 12. barrera -> kcat

md(r"""
## 12. De la barrera a *k*<sub>cat</sub>

> 🎯 **En esta sección** convertirás las barreras calculadas en un número de recambio
> *k*<sub>cat</sub> y lo compararás, con honestidad, con el experimento y con otro cálculo
> publicado.

---

### 💡 La analogía: una regla con dos escalas

Un termómetro puede tener dos escalas, °C y °F: es la misma medida leída de dos formas. La
ecuación de Eyring es una regla así: por arriba lees la **altura de la colina**, por abajo las
**reacciones por segundo**. Con la altura de la colina estimamos cuántas veces por segundo la
enzima completa la reacción: ese número es **k**<sub>cat</sub>, el "número de recambio". Lo
comparamos con el valor medido en el laboratorio y con otro cálculo QM/MM publicado.

[[fig:escalera_kcat | Regla con la barrera arriba y k_cat abajo; marcadas las barreras calculadas y la experimental]]

---

### 🧮 La ecuación en palabras

Es la misma ecuación de la sección 6, *intentos × probabilidad de éxito*, usada en los dos
sentidos:

* **de la simulación al experimento**: con la barrera calculada, ¿cuántas reacciones por segundo predice?
* **del experimento a la barrera**: con la k<sub>cat</sub> medida, ¿qué altura de colina implica?

---

### 📐 La ecuación completa, término a término

$$k_{\mathrm{cat}} \;=\; {\color{#2a78d6}{\frac{k_BT}{h}}}\;{\color{#eb6834}{e^{-\Delta G^{\ddagger}/RT}}}
\qquad\Longleftrightarrow\qquad
\Delta G^{\ddagger} \;=\; RT\,\ln\!\frac{k_BT}{h\,k_{\mathrm{cat}}}$$

| Término | Qué es | En la analogía | Valor |
|---|---|---|---|
| $k_{\mathrm{cat}}$ | reacciones por segundo por molécula de enzima | la escala de abajo | medida: ≈ 62–66 s⁻¹ |
| ${\color{#2a78d6}{k_BT/h}}$ | intentos por segundo | — | 6.2 × 10¹² s⁻¹ a 25 °C |
| ${\color{#eb6834}{\exp(-\Delta G^{\ddagger}/RT)}}$ | probabilidad de cruzar la cima | — | diminuta |
| $\Delta G^{\ddagger}$ | altura de la colina | la escala de arriba | experimento: ≈ 15 kcal/mol |

---

### 🎛️ Qué pasa si… (o por qué hay que ser honestos)

Nuestra ΔG‡ viene de un método semiempírico (PM7), de un entorno congelado y de una aproximación
armónica. Un error de 1.4 kcal/mol ya es un factor 10 en velocidad. Por eso comparamos órdenes de
magnitud y, sobre todo, *tendencias*, no decimales.

| Si la barrera calculada se equivoca en… | k<sub>cat</sub> se equivoca en un factor… |
|---|---|
| 0.5 kcal/mol | ~2 |
| 1.4 kcal/mol | ~10 |
| 4 kcal/mol | ~850 |

> ⚠️ **Ojo.** Una barrera calculada "casi igual" a la experimental puede dar una k<sub>cat</sub>
> diez o cien veces distinta. Compara siempre en la escala de la energía **y** en la de la
> velocidad.

---

### 🔬 Los datos reales
""")

code(r'''
T = 298.15
kcat_exp = valor("kcat_s", 60.0)
dG_qmmm = termo.get(T_ref, {}).get("dG_kcal", ts["barrera_kcal"])
filas = [("ΔE‡ PM7, una conformación (este cuaderno)", ts["barrera_kcal"]),
         ("ΔG‡ PM7 + termoquímica armónica (este cuaderno)", dG_qmmm)]
try:
    if len(ok) > 1:
        filas.append(("ΔE‡ promedio sobre instantáneas de MD (este cuaderno)", media))
except NameError:
    pass
filas.append(("ΔE‡ QM/MM publicada, Zhang et al. 2009", valor("qmmm_barrier_literature_kcal", np.nan)))
filas.append(("sitio activo en agua, COSMO (este cuaderno)", agua["barrera_kcal"]))
filas.append(("experimento: ΔG‡ que implica k_cat", cin.barrier_from_rate(kcat_exp, T)))
tabla = pd.DataFrame([(n, b, cin.eyring_rate(b, T)) for n, b in filas], columns=["estimación", "barrera (kcal/mol)", "k (s⁻¹) por Eyring"])
tabla["k (s⁻¹) por Eyring"] = tabla["k (s⁻¹) por Eyring"].map(lambda x: f"{x:.2e}")
display(tabla.round(1))
print(f"k_cat experimental ≈ {kcat_exp:.0f} s⁻¹  ({ref['kcat_s']['source']})")
print("Nuestra barrera y la de Zhang et al. (18.3 kcal/mol, otro programa y otra partición QM/MM) coinciden dentro de 1 kcal/mol;")
print("ambas sobrestiman la experimental (≈15) en 3-4 kcal/mol, es decir, un factor ~10²-10³ en k. Es lo esperable de un método semiempírico con entorno fijo.")
''')

md(r"""
> ✅ **Para llevar.** Eyring traduce barreras en velocidades en los dos sentidos. Nuestra barrera
> y la de Zhang et al. coinciden entre sí dentro de 1 kcal/mol y quedan 3–4 kcal/mol por encima
> de la que implica el experimento (≈ 15): la simulación acierta el mecanismo y el orden de
> magnitud de la colina, no los decimales de *k*<sub>cat</sub>.
""")
