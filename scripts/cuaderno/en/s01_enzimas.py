"""Section 1. Enzymes."""
from ..celdas import code, md

# ============================================================================ 1. enzymes
md(r"""
## 1. What is an enzyme, and why glucokinase?

> 🎯 **In this section** you will learn what an enzyme does (and what it does **not** do), the
> vocabulary we will use throughout the course, which reaction glucokinase catalyzes and why it
> matters so much to your body.

---

### 💡 The analogy: a guide who knows the lowest pass

Imagine you want to cross a mountain range. Without a guide, you would search for a pass blindly
and it would take you weeks. A good guide knows the **lowest pass** and takes you through it in
hours. The mountain is still there, but the route is different. An **enzyme** is that guide: it
does not change the starting point or the destination of a chemical reaction, but it takes it along
a much easier route, and that is why it happens thousands or millions of times faster.

[[fig:guia_paso_bajo | Two energy profiles between the same reactants and products: without the enzyme the hill is high; with the enzyme it is much lower]]

Notice what the drawing does **not** change: the two valleys are at the same height with and
without the enzyme. An enzyme speeds up the reaction, but it does not decide which way it goes or
how much product there will be at the end.

---

### 📖 Vocabulary we will use all the time

| Term | What it is | In the analogy | In glucokinase |
|---|---|---|---|
| **Substrate (S)** | the molecule the enzyme transforms | the traveler in the starting valley | glucose and ATP |
| **Product (P)** | what comes out | the traveler already in the destination valley | glucose‑6‑phosphate and ADP |
| **Active site** | the pocket of the enzyme where the substrates fit and the chemistry happens | the meeting point with the guide | the pocket with Asp205, Lys169 and Mg²⁺ |
| **ES complex** | the enzyme with the substrate already inside | traveler and guide, together | glucokinase·glucose·ATP |
| **Transition state (TS)** | the hardest moment of the reaction | the top of the pass | the phosphate halfway across |

Almost all of this notebook revolves around the **transition state**.

---

### 🔬 Our enzyme: glucokinase

**Glucokinase** takes a glucose and attaches to it a phosphate group that comes from ATP. The
product, glucose‑6‑phosphate, can no longer leave the cell: it is the first step to "trap" and
use glucose.

[[fig:reaccion_glucoquinasa | Reaction scheme: the gamma phosphate of ATP jumps to the O6 oxygen of glucose, with Mg2+, Lys169 and Asp205 helping]]

Glucokinase lives in the liver and in the β cells of the pancreas, where it works as a **glucose
sensor**: its activity decides how much insulin is released.

[[fig:sensor_glucosa | Blood glucose enters the beta cell, glucokinase senses it and insulin is released]]

Mutations in its gene cause a hereditary form of diabetes (GCK‑MODY, formerly MODY2) and, if they
make it hyperactive, congenital hypoglycemia.

**Three things make it ideal for learning kinetics:**

1. Its reaction is a classic **phosphoryl transfer**, with Mg²⁺ and a catalytic base (Asp205).
2. Its rate curve is **not** the Michaelis–Menten hyperbola but a **sigmoid** (section 14).
3. It has well-studied inhibitors, a regulatory protein and pharmacological activators (section 15).

---

### 🔬 The real data: numbers measured in the laboratory

Everything that follows is compared with these values, taken from original papers (*source*
column). Keep them in mind: they are the "experimental result" that the theory must explain.

> ⚠️ **Watch out.** Different laboratories obtain k<sub>cat</sub> between 38 and 66 s⁻¹ for the same
> enzyme: in biochemistry, a number always comes with its conditions (temperature, pH, buffer).
""")

code(r'''
# @title 📋 Glucokinase numbers measured in the laboratory
filas = []
for k, v in ref.items():
    if isinstance(v, dict) and "value" in v:
        valor_txt = v["value"]
        if isinstance(valor_txt, (int, float)) and abs(valor_txt) >= 1e5:  # 1.1e12 → 1.1 × 10¹²
            m, e = f"{valor_txt:.1e}".split("e")
            valor_txt = f"{m} × 10" + str(int(e)).translate(str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))
        filas.append((v.get("label", k), valor_txt, v.get("unit", ""), v.get("conditions", ""), v.get("source", "")))
tabla = pd.DataFrame(filas, columns=["parameter", "value", "unit", "conditions", "source"])
viz.mostrar(
    viz.tabla(tabla, titulo="Human glucokinase in the laboratory",
              nota="Every number comes with its conditions and its source: different laboratories measure somewhat different values."),
    viz.mensaje(ref.get("notes", ""), tipo="dato", titulo="Note on the data") if ref.get("notes") else None,
)
''')

md(r"""
> ✅ **Takeaway.** An enzyme is a guide: it lowers the hill between reactants and products without
> changing either the starting point or the destination. Glucokinase moves a phosphate from ATP to
> glucose, and as a glucose sensor it decides how much insulin the pancreas releases.
""")
