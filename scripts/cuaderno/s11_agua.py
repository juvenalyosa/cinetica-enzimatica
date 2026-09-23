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
agua = res["etapas"]["agua"]
print("1) SADDLE (QST2): estimación del TS %.1f kcal/mol por encima del reactivo (ξ = %+.2f Å)" % (agua["barrera_saddle_qst2_kcal"], agua["saddle_xi"]))
print("2) Refinamiento (%s): barrera = %.1f kcal/mol; ΔE(reacción) = %+.1f kcal/mol" % (agua["metodo_ts"], agua["barrera_kcal"], agua["dE_reaccion"]))
print("   geometría del TS en agua: d(Pγ–O6) = %.2f Å, d(Pγ–O3β) = %.2f Å, d(O6–H) = %.2f Å" % (agua["ts_d_PG_O6"], agua["ts_d_PG_O3B"], agua["ts_d_O6_H"]))
print("3) FORCETS de MOPAC, frecuencias más bajas (cm⁻¹):", np.round(agua["frecuencias_mas_bajas"], 1))
if "frecuencias_ase_mas_bajas" in agua:
    print("   Hessiano numérico (ASE), frecuencias más bajas (cm⁻¹):", np.round(agua["frecuencias_ase_mas_bajas"], 1),
          "→ modos imaginarios:", agua["validacion_ts_ase"]["imaginary_mode_count"])
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
try:
    cam_agua = datos.csv("qmmm/agua_camino_descenso.csv")
    fig = viz.plot_energy_profile(cam_agua["xi"].values, cam_agua["energia_kcal"].values - agua["E_reactivo"], relative=False,
                                  xlabel="ξ (Å)", smooth=False, ts_index=int(np.argmax(cam_agua["energia_kcal"].values)),
                                  title="Camino de reacción del sitio activo en agua (COSMO)",
                                  subtitle="Energías relativas al reactivo en agua; el máximo es el TS refinado")
    fig;
except FileNotFoundError:
    print("(camino en agua no disponible)")
''')

code(r'''
con = [("E·S", 0.0), ("TS", ts["barrera_kcal"]), ("E·P", P["dE_reaccion_kcal"])]
sin = [("R", 0.0), ("TS", agua["barrera_kcal"]), ("P", agua["dE_reaccion"])]
fig = viz.plot_energy_levels(con, compare=sin, ts_indices=[1], label="dentro de la enzima", compare_label="sitio activo en agua",
                             title="Dentro de la enzima frente a fuera de ella",
                             subtitle="La diferencia entre las dos colinas es el efecto del resto de la proteína")
fig;
print("Efecto del entorno proteico sobre la barrera: %+.1f kcal/mol" % (ts["barrera_kcal"] - agua["barrera_kcal"]))
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
