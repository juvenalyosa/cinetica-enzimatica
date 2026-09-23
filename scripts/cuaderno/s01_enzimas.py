"""Sección 1. enzimas."""
from .celdas import code, md

# ============================================================================ 1. enzimas
md(r"""
## 1. ¿Qué es una enzima y por qué la glucoquinasa?

> 🎯 **En esta sección** aprenderás qué hace una enzima (y qué **no** hace), el vocabulario que
> usaremos todo el curso, qué reacción cataliza la glucoquinasa y por qué es tan importante para
> tu cuerpo.

---

### 💡 La analogía: un guía que conoce el paso más bajo

Imagina que quieres cruzar una cordillera. Sin guía, buscarías el paso a ciegas y tardarías
semanas. Un buen guía conoce el **paso más bajo** y te lleva por él en horas. La montaña sigue
ahí, pero el camino es otro. Una **enzima** es ese guía: no cambia el punto de partida ni el de
llegada de una reacción química, pero la lleva por un camino mucho más fácil, y por eso ocurre
miles o millones de veces más rápido.

[[fig:guia_paso_bajo | Dos perfiles de energía entre los mismos reactivos y productos: sin enzima la colina es alta; con enzima es mucho más baja]]

Fíjate en lo que el dibujo **no** cambia: los dos valles están a la misma altura con y sin
enzima. Una enzima acelera la reacción, pero no decide hacia dónde va ni cuánto producto habrá
al final.

---

### 📖 Vocabulario que usaremos todo el tiempo

| Término | Qué es | En la analogía | En la glucoquinasa |
|---|---|---|---|
| **Sustrato (S)** | la molécula que la enzima transforma | el viajero en el valle de partida | glucosa y ATP |
| **Producto (P)** | lo que sale | el viajero ya en el valle de llegada | glucosa‑6‑fosfato y ADP |
| **Sitio activo** | el hueco de la enzima donde encajan los sustratos y ocurre la química | el punto de encuentro con el guía | el bolsillo con Asp205, Lys169 y Mg²⁺ |
| **Complejo ES** | la enzima con el sustrato ya dentro | viajero y guía, juntos | glucoquinasa·glucosa·ATP |
| **Estado de transición (TS)** | el momento más difícil de la reacción | la cima del paso | el fosfato a medio camino |

Casi todo este cuaderno gira en torno al **estado de transición**.

---

### 🔬 Nuestra enzima: la glucoquinasa

La **glucoquinasa** toma una glucosa y le pega un grupo fosfato que viene del ATP. El
producto, glucosa‑6‑fosfato, ya no puede salir de la célula: es el primer paso para "atrapar" y
usar la glucosa.

[[fig:reaccion_glucoquinasa | Esquema de la reacción: el fosfato gamma del ATP salta al oxígeno O6 de la glucosa, con Mg2+, Lys169 y Asp205 ayudando]]

La glucoquinasa vive en el hígado y en las células β del páncreas, donde funciona como **sensor
de glucosa**: su actividad decide cuánta insulina se libera.

[[fig:sensor_glucosa | La glucosa de la sangre entra a la célula beta, la glucoquinasa la detecta y se libera insulina]]

Mutaciones en su gen causan una diabetes hereditaria (GCK‑MODY, antes MODY2) y, si la vuelven
hiperactiva, hipoglucemia congénita.

**Tres cosas la hacen ideal para aprender cinética:**

1. Su reacción es una **transferencia de fosforilo** clásica, con Mg²⁺ y una base catalítica (Asp205).
2. Su curva de velocidad **no** es la hipérbola de Michaelis–Menten sino una **sigmoide** (sección 14).
3. Tiene inhibidores, una proteína reguladora y activadores farmacológicos bien estudiados (sección 15).

---

### 🔬 Los datos reales: números medidos en el laboratorio

Todo lo que sigue se compara con estos valores, tomados de artículos originales (columna
*fuente*). Guárdalos: son el "resultado experimental" que la teoría debe explicar.

> ⚠️ **Ojo.** Distintos laboratorios obtienen k<sub>cat</sub> entre 38 y 66 s⁻¹ para la misma
> enzima: en bioquímica, un número siempre viene con sus condiciones (temperatura, pH, tampón).
""")

code(r'''
filas = []
for k, v in ref.items():
    if isinstance(v, dict) and "value" in v:
        filas.append((v.get("label", k), v["value"], v.get("unit", ""), v.get("conditions", ""), v.get("source", "")))
tabla = pd.DataFrame(filas, columns=["parámetro", "valor", "unidad", "condiciones", "fuente"])
pd.set_option("display.max_colwidth", 90)
display(tabla)
print(ref.get("notes", ""))
''')

md(r"""
> ✅ **Para llevar.** Una enzima es un guía: baja la colina entre reactivos y productos sin
> cambiar ni el punto de partida ni el de llegada. La glucoquinasa mueve un fosfato del ATP a la
> glucosa, y como sensor de glucosa decide cuánta insulina libera el páncreas.
""")
