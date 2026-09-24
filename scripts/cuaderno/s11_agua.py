"""Sección 11. agua."""
from .celdas import code, md

# ============================================================================ 11. agua

md(r"""
## 11. La misma reacción **fuera** de la enzima

> 🎯 **En esta sección** medirás cuánto de la catálisis hace *el resto* de la proteína, repitiendo
> la reacción con el sitio activo sumergido en agua, y conocerás los buscadores de estado de
> transición nativos de MOPAC.

---

### 💡 La analogía: la misma canción, fuera del estadio

¿Cuánto del sonido de un concierto se debe a la banda y cuánto a la acústica del estadio? Para
saberlo, la banda toca la misma canción en un estudio. Aquí la "banda" es el sitio activo y el
"estadio", el resto de la proteína.

Sacamos el sitio activo (el mismo clúster QM, con sus anclajes fijos) y lo ponemos en **agua**
(disolvente implícito COSMO), sin el campo eléctrico de los otros 2300 átomos. Si la barrera
cambia, ese cambio es obra del entorno proteico.

[[fig:dentro_fuera | Perfiles de energía dentro de la enzima y del clúster en agua; qué conserva y qué pierde el clúster]]

---

### 🧮 La ecuación en palabras

**¿Para qué queremos una ecuación?** Para poner un número al "aporte del estadio":

$$\text{efecto del resto de la proteína} \;=\; \text{barrera dentro de la enzima} \;-\; \text{barrera del clúster en agua}$$

$$\Delta\Delta E^{\ddagger}_{\text{entorno}} \;=\; {\color{#2a78d6}{\Delta E^{\ddagger}_{\text{enzima}}}} \;-\; {\color{#898781}{\Delta E^{\ddagger}_{\text{agua}}}} \;=\; 19.1 - 21.6 \;\approx\; -2.5\ \text{kcal/mol}$$

Un número **negativo** significa que el entorno proteico **baja** la colina. Con la regla de oro
(1.36 kcal/mol = ×10), 2.5 kcal/mol equivalen a que la reacción vaya unas **70 veces** más
rápido dentro de la enzima que en el clúster aislado.

---

### 🧰 Las herramientas: los métodos nativos de MOPAC

Fuera de la enzima no hay potencial externo que "congelar", así que aquí sí usamos los
buscadores de estado de transición nativos de MOPAC, con los nombres con que los orquesta la
suite Leonardo:

| Nombre en Leonardo | Palabra clave MOPAC | Qué hace | En la analogía del paso de montaña |
|---|---|---|---|
| QST2 | `SADDLE` | Busca el TS a partir de reactivo y producto (dos extremos) | caminar desde los dos valles hasta encontrarse |
| QST3 / TS | `TS` | Refina el punto de silla siguiendo el modo de curvatura negativa | afinar la posición exacta del collado |
| Validación | `FORCETS` | Frecuencias sobre las coordenadas libres: debe haber una imaginaria | comprobar que solo una dirección baja |
| Camino | `IRC=1*` | Sigue el camino intrínseco de reacción en ambas direcciones | bajar del collado a los dos valles |

Un detalle práctico: `SADDLE` da una *estimación* del TS que hay que refinar. Aquí el
refinamiento se hizo con el método del dímero y después con la palabra clave `TS` de MOPAC,
que llegaron al mismo punto de silla.

---

### 🔬 Los datos reales
""")

code(r'''
# @title 💧 El sitio activo en agua: resultados de MOPAC
agua = res["etapas"]["agua"]
display(viz.tarjetas(
    [("1 · SADDLE (QST2)", f"{agua['barrera_saddle_qst2_kcal']:.1f}", "kcal/mol", f"primera estimación del TS, ξ = {agua['saddle_xi']:+.2f} Å", "gris"),
     ("2 · TS refinado", f"{agua['barrera_kcal']:.1f}", "kcal/mol", agua["metodo_ts"], "naranja"),
     ("Geometría del TS", f"{agua['ts_d_PG_O6']:.2f} / {agua['ts_d_PG_O3B']:.2f}", "Å",
      f"Pγ···O6 / Pγ···O3β; el O6–H sigue en {agua['ts_d_O6_H']:.2f} Å", "azul"),
     ("ΔE de reacción", f"{agua['dE_reaccion']:+.1f}", "kcal/mol", "en agua la reacción es cuesta abajo", "agua")],
    titulo="La misma reacción, fuera de la enzima (PM7 + agua implícita COSMO)"))
conjuntos = [("3 · FORCETS (MOPAC)", agua["frecuencias_mas_bajas"])]
if "frecuencias_ase_mas_bajas" in agua:
    conjuntos.append(("Hessiano más fino (ASE)", agua["frecuencias_ase_mas_bajas"]))
fig = viz.plot_frequencies(conjuntos, title="Mirar la magnitud, no solo contar las imaginarias",
                           subtitle="Gris: ruido de la cavidad del disolvente. Naranja: movimientos que «caen». El Hessiano más fino deja uno solo: la reacción")
n_f = max(len(f) for _, f in conjuntos)
tabla_f = pd.DataFrame({"modo": np.arange(1, n_f + 1),
                        **{nombre.split(" (")[0].replace("3 · ", ""): np.r_[np.asarray(f, float), [np.nan] * (n_f - len(f))]
                           for nombre, f in conjuntos}})
viz.mostrar(fig, viz.datos(
    tabla_f, "las frecuencias más bajas del TS en agua",
    "Cada columna es un cálculo del Hessiano; cada fila, uno de los modos más lentos ordenados de menor a mayor. "
    "Los valores entre −40 y +40 cm⁻¹ son ruido numérico de la cavidad del disolvente.",
    y=list(tabla_f.columns[1:]), unidades={c: "cm⁻¹" for c in tabla_f.columns[1:]}, formatos={c: "{:.1f}" for c in tabla_f.columns[1:]},
    resaltar={0: ("la reacción", "naranja")}))
''')

md(r"""
**Una lección sobre validar.** El Hessiano de `FORCETS` sobre la superficie COSMO muestra varias
frecuencias "imaginarias" pequeñas (|ν| < 40 cm⁻¹): son ruido numérico de la cavidad del
disolvente, no movimientos reales; por eso conviene mirar la *magnitud* de las frecuencias y no
solo contarlas. El Hessiano por diferencias finitas más finas (ASE) deja un solo modo claro
(≈ −160 cm⁻¹, el fósforo que salta) y un residuo de −37 cm⁻¹ atribuible a la cavidad. El `IRC`
de MOPAC abortó por un fallo interno del programa en esta superficie con disolvente (también
eso pasa en la práctica); el camino se obtuvo con el mismo descenso desde el TS de la sección 9.
""")

code(r'''
# @title 📈 El camino de reacción en agua
try:
    cam_agua = datos.csv("qmmm/agua_camino_descenso.csv")
    fig = viz.plot_energy_profile(cam_agua["xi"].values, cam_agua["energia_kcal"].values - agua["E_reactivo"], relative=False,
                                  xlabel="ξ (Å)", smooth=False, sort=False, ts_index=int(np.argmax(cam_agua["energia_kcal"].values)),
                                  annotate_states=True, annotate_reaction=True, state_names=("reactivo", "TS", "producto"),
                                  color=viz.PALETTE[6],
                                  title=f"En agua la colina mide {agua['barrera_kcal']:.1f} kcal/mol",
                                  subtitle="Camino de descenso desde el TS refinado; energías relativas al reactivo en agua")
    e_rel = cam_agua["energia_kcal"].values - agua["E_reactivo"]
    viz.mostrar(fig, viz.datos(
        pd.DataFrame({"paso": cam_agua["cuadro"], "ξ": cam_agua["xi"], "E": e_rel}), "el camino en agua",
        "Cada fila es un paso del descenso desde el TS en agua, en el orden del camino.",
        x="ξ", y="E", calculadas={"E": "energía del paso − energía del reactivo en agua"},
        unidades={"ξ": "Å", "E": "kcal/mol"}, formatos={"ξ": "{:+.2f}", "E": "{:.1f}"},
        resaltar={int(np.argmax(e_rel)): ("TS", "naranja")}, barra="E"))
except FileNotFoundError:
    display(viz.mensaje("Camino en agua no disponible.", "ojo"))
''')

md(r"""
La misma película que en la sección 8, pero ahora **sin el resto de la proteína**: solo el sitio
activo en agua. Compara la altura de la colina en la gráfica con la de la enzima.
""")

code(r'''
# @title 🎬 La reacción en agua, sincronizada con su energía
visor3d.pelicula_reaccion("agua")
''')

code(r'''
# @title ⚖️ Dentro de la enzima frente a fuera de ella
con = [("E·S", 0.0), ("TS", ts["barrera_kcal"]), ("E·P", P["dE_reaccion_kcal"])]
sin = [("R", 0.0), ("TS", agua["barrera_kcal"]), ("P", agua["dE_reaccion"])]
efecto = ts["barrera_kcal"] - agua["barrera_kcal"]
fig = viz.plot_energy_levels(con, compare=sin, ts_indices=[1], label="dentro de la enzima", compare_label="sitio activo en agua",
                             title=f"El resto de la proteína baja la colina {abs(efecto):.1f} kcal/mol",
                             subtitle="Azul/naranja: dentro de la enzima. Gris: el mismo sitio activo en agua")
viz.mostrar(fig, viz.tarjetas(
    [("Efecto del entorno proteico", f"{efecto:+.1f}", "kcal/mol", "barrera en la enzima − barrera en agua", "azul"),
     ("En velocidad", f"× {10 ** (-efecto / 1.364):.0f}", "", "cada 1.36 kcal/mol es un factor 10", "naranja")]),
    viz.datos(pd.DataFrame({"estado": ["reactivo", "estado de transición", "producto"],
                            "en la enzima": [e for _, e in con], "en agua": [e for _, e in sin],
                            "diferencia": [a - b for (_, a), (_, b) in zip(con, sin)]}),
              "los dos diagramas",
              "Las alturas de los dos diagramas, lado a lado. La fila del estado de transición es la que decide la velocidad.",
              y=["en la enzima", "en agua"], calculadas={"diferencia": "en la enzima − en agua"},
              unidades={"en la enzima": "kcal/mol", "en agua": "kcal/mol", "diferencia": "kcal/mol"},
              formatos={"en la enzima": "{:.1f}", "en agua": "{:.1f}", "diferencia": "{:+.1f}"},
              resaltar={1: ("TS", "naranja")}))
''')

md(r"""
**Cómo leer esta comparación.** El clúster "en agua" **no** es la reacción sin catalizar:
todavía contiene la base catalítica (Asp205), la carga positiva de Lys169 y el Mg²⁺, es decir,
la maquinaria química esencial. Por eso las dos colinas se parecen: la diferencia (unas pocas
kcal/mol) mide lo que aporta el *resto* de la proteína, sobre todo su campo eléctrico y el
hecho de mantener los reactivos alineados. La reacción de verdad sin enzima, glucosa y ATP
solos en agua, es extraordinariamente lenta: la hidrólisis espontánea de un fosfato dianión
tiene una vida media de ~10¹² años (Lad, Williams y Wolfenden, 2003; ΔG‡ ≈ 44 kcal/mol),
porque ninguna base, ningún catión ni ningún campo eléctrico están ahí para ayudar. Entre esa
reacción y k<sub>cat</sub> ≈ 60 s⁻¹ hay ~10²¹ veces: ese es el tamaño real de la catálisis.
""")

md(r"""
> ✅ **Para llevar.** El sitio activo solo, en agua, ya tiene casi toda la maquinaria química
> (Asp205, Lys169, Mg²⁺), por eso su barrera (21.6) se parece a la de la enzima (19.1). El resto
> de la proteína aporta unas 2.5 kcal/mol más (un factor ~70). La catálisis completa, frente a
> la reacción sin ninguna ayuda, es de ~10²¹ veces.
""")
