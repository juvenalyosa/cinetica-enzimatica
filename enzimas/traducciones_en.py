"""Traducciones al inglés de los textos de la biblioteca (gráficas, tarjetas, visores 3D, exploradores).

Clave: el texto exacto en español que aparece dentro de ``t(...)``; valor: su versión en inglés.
Los textos con números son plantillas de ``str.format`` y conservan los mismos campos ``{...}``.
"""

EN: dict[str, str] = {
    # ------------------------------------------------------------------ comunes: ejes, estados, unidades
    "Coordenada de reacción": "Reaction coordinate",
    "Energía relativa (kcal/mol)": "Relative energy (kcal/mol)",
    "Energía (kcal/mol)": "Energy (kcal/mol)",
    "Energía relativa al reactivo (kcal/mol)": "Energy relative to the reactant (kcal/mol)",
    "Tiempo (ps)": "Time (ps)",
    "Tiempo (s)": "Time (s)",
    "Concentración (µM)": "Concentration (µM)",
    "Frecuencia (cm⁻¹)": "Frequency (cm⁻¹)",
    "Producto formado (µM)": "Product formed (µM)",
    "Barrera (kcal/mol)": "Barrier (kcal/mol)",
    "Temperatura (K)": "Temperature (K)",
    "RMSD Cα (Å)": "Cα RMSD (Å)",
    "d(Pγ–O6) (Å)": "d(Pγ–O6) (Å)",
    "d(Pγ–O3β) (Å)": "d(Pγ–O3β) (Å)",
    "d(O6–OD1 Asp205) (Å)": "d(O6–OD1 Asp205) (Å)",
    "d(Mg–Pγ) (Å)": "d(Mg–Pγ) (Å)",
    "d(NZ Lys169–O6) (Å)": "d(NZ Lys169–O6) (Å)",
    "reactivo": "reactant",
    "reactivos": "reactants",
    "producto": "product",
    "productos": "products",
    "estado de transición": "transition state",
    "cima: {v} kcal/mol": "summit: {v} kcal/mol",
    "con enzima": "with enzyme",
    "sin enzima": "without enzyme",
    "cristal": "crystal",
    "cristal minimizado": "minimized crystal",
    "experimento": "experiment",
    "días": "days",
    "años": "years",
    "tiempo": "time",
    "tiempo (ps)": "time (ps)",
    "energía (kcal/mol)": "energy (kcal/mol)",
    "Energía": "Energy",
    "glucosa": "glucose",
    "trifosfato": "triphosphate",
    "dominio grande": "large domain",
    "dominio pequeño": "small domain",
    "Michaelis–Menten": "Michaelis–Menten",
    "Hill": "Hill",
    "Eyring": "Eyring",
    "pH": "pH",
    "k_cat": "k_cat",
    "K_M (ATP)": "K_M (ATP)",
    "-": "-",
    "MD 250 ps": "MD 250 ps",
    "MD 500 ps": "MD 500 ps",
    "MD 750 ps": "MD 750 ps",
    "MD 1000 ps": "MD 1000 ps",

    # ------------------------------------------------------------------ viz: dinámica molecular
    "zona de ataque cercano  (d < {u:g} Å)": "near-attack zone  (d < {u:g} Å)",
    "Distancia entre el fósforo γ del ATP y el O6 de la glucosa": "Distance between the γ-phosphorus of ATP and O6 of glucose",
    "del tiempo": "of the time",
    "fotogramas": "frames",
    "sin columna «{c}»": "no column \"{c}\"",
    "La proteína es estable": "The protein is stable",
    "RMSD de los Cα (Å)": "Cα RMSD (Å)",
    "O6 ··· Asp205 (la base catalítica)": "O6 ··· Asp205 (the catalytic base)",
    "puente de H": "H-bond",
    "media {m:.2f} Å": "mean {m:.2f} Å",
    "La enzima pasa el {f:.0%} del tiempo lista para reaccionar": "The enzyme spends {f:.0%} of the time ready to react",
    "Resumen de la dinámica molecular": "Molecular dynamics summary",
    "Cada punto es un fotograma de la película (cada 10 ps); la línea gruesa es la media móvil":
        "Each point is a frame of the movie (every 10 ps); the thick line is the moving average",

    # ------------------------------------------------------------------ viz: frecuencias, progreso, instantáneas
    "ruido\nnumérico": "numerical\nnoise",
    "← imaginarias: la geometría «cae»": "← imaginary: the geometry \"falls\"",
    "reales: la geometría vibra →": "real: the geometry vibrates →",
    "Pendiente inicial (línea discontinua)": "Initial slope (dashed line)",
    "hay producto estable": "there is a stable product",
    "(camino completo)": "(full path)",
    "instantáneas de la MD (sin producto)": "MD snapshots (no product)",

    # ------------------------------------------------------------------ viz: cinética
    "casi lineal:\ncada molécula cuenta": "almost linear:\nevery molecule counts",
    "saturación: añadir más\nsustrato casi no ayuda": "saturation: adding more\nsubstrate barely helps",
    "por debajo de S₀.₅\nla sigmoide casi no responde": "below S₀.₅\nthe sigmoid barely responds",
    "Linealizaciones de Michaelis–Menten": "Michaelis–Menten linearizations",
    "[S] baja: pocos puntos\nque pesan demasiado": "low [S]: a few points\nthat weigh too much",
    " (control)": " (control)",
    "velocidad perdida\npor el inhibidor": "rate lost\nto the inhibitor",
    "Inhibición": "Inhibition",
    "Inhibición competitiva": "Competitive inhibition",
    "Inhibición acompetitiva": "Uncompetitive inhibition",
    "Inhibición no competitiva": "Noncompetitive inhibition",
    "Inhibición mixta": "Mixed inhibition",
    "las rectas se cruzan en el eje 1/v₀: Vmax no cambia, Km aparente aumenta":
        "the lines cross on the 1/v₀ axis: Vmax does not change, apparent Km increases",
    "rectas paralelas: Vmax y Km aparentes disminuyen en la misma proporción":
        "parallel lines: apparent Vmax and Km decrease in the same proportion",
    "las rectas se cruzan en el eje 1/[S]: Km no cambia, Vmax disminuye":
        "the lines cross on the 1/[S] axis: Km does not change, Vmax decreases",
    "las rectas se cruzan a la izquierda del eje 1/v₀ (fuera de los ejes)":
        "the lines cross to the left of the 1/v₀ axis (off the axes)",
    "aquí se cruzan": "they cross here",
    "rectas paralelas: nunca se cruzan": "parallel lines: they never cross",
    "Simulación del mecanismo": "Simulation of the mechanism",
    "S · sustrato": "S · substrate",
    "P · producto": "P · product",
    "El sustrato se gasta, el producto aparece": "Substrate is used up, product appears",
    "E · enzima libre": "E · free enzyme",
    "ES · complejo": "ES · complex",
    "En milisegundos, ES llega a una meseta": "Within milliseconds, ES reaches a plateau",
    "Enzima libre y complejo ES": "Free enzyme and ES complex",
    "Michaelis–Menten emerge del mecanismo": "Michaelis–Menten emerges from the mechanism",
    "v₀ de la simulación (EDO)": "v₀ from the simulation (ODE)",
    "ajuste de Michaelis–Menten": "Michaelis–Menten fit",
    "Gráfico de Arrhenius": "Arrhenius plot",
    "pendiente = −Ea/R": "slope = −Ea/R",
    "Gráfico de Eyring": "Eyring plot",
    "pendiente = −ΔH‡/R": "slope = −ΔH‡/R",
    "Perfil de pH": "pH profile",
    "Gráfico de Hill": "Hill plot",
    "referencia n = 1\n(sin cooperatividad)": "reference n = 1\n(no cooperativity)",
    "pendiente n = {n:.2f}": "slope n = {n:.2f}",

    # ------------------------------------------------------------------ viz: tabla de datos
    "eje x": "x axis",
    "eje y": "y axis",
    "calculada": "calculated",
    "{c} en barra": "{c} as a bar",
    "Se muestran {k} de {n} filas, repartidas a lo largo de la curva.": "Showing {k} of {n} rows, spread along the curve.",
    "Los datos de la gráfica:": "The data behind the plot:",
    "punto": "point",
    "puntos": "points",
    "clic para ocultar": "click to hide",
    "clic para ver": "click to show",

    # ------------------------------------------------------------------ exploradores (interactivo)
    "Mueve los deslizadores: la gráfica y la frase de abajo se actualizan.":
        "Move the sliders: the plot and the sentence below update.",
    "Qué dicen estos números": "What these numbers say",
    "saturación": "saturation",
    "mitad del máximo\nen [S] = Km = {km:g} mM": "half the maximum\nat [S] = Km = {km:g} mM",
    "Velocidad inicial frente a sustrato": "Initial rate versus substrate",
    "glucosa en sangre (5 mM)": "blood glucose (5 mM)",
    "enzima ocupada (%)": "enzyme occupied (%)",
    "¿Qué fracción de la enzima trabaja?": "What fraction of the enzyme is working?",
    "Con Km = {km:g} mM, la enzima va a media máquina cuando [S] = {km:g} mM":
        "With Km = {km:g} mM, the enzyme runs at half speed when [S] = {km:g} mM",
    "Izquierda: la hipérbola. Derecha: la misma idea como ocupación de la enzima, [S]/(Km + [S]).":
        "Left: the hyperbola. Right: the same idea as enzyme occupancy, [S]/(Km + [S]).",
    "Con Km = {km:g} mM, a [S] = {km:g} mM la enzima trabaja al 50 % de su máximo "
    "(v₀ = {v_half:g} µM/s); a 5 mM (glucosa en sangre) trabaja al {pct:.0f} % y "
    "necesita [S] = 9·Km = {km9:g} mM para llegar al 90 %.":
        "With Km = {km:g} mM, at [S] = {km:g} mM the enzyme works at 50 % of its maximum "
        "(v₀ = {v_half:g} µM/s); at 5 mM (blood glucose) it works at {pct:.0f} % and "
        "it needs [S] = 9·Km = {km9:g} mM to reach 90 %.",
    "Vmax: velocidad máxima (µM/s)": "Vmax: maximum rate (µM/s)",
    "Km: [S] a media velocidad (mM)": "Km: [S] at half rate (mM)",
    "La hipérbola de Michaelis–Menten": "The Michaelis–Menten hyperbola",
    "S (sustrato) se gasta": "S (substrate) is used up",
    "P (producto) se acumula": "P (product) builds up",
    "La reacción completa": "The whole reaction",
    "arranque": "start-up",
    "estado estacionario: ES casi constante": "steady state: ES almost constant",
    "E (enzima libre)": "E (free enzyme)",
    "ES (complejo)": "ES (complex)",
    "Zoom a los primeros instantes: la enzima": "Zoom into the first moments: the enzyme",
    "ES se estabiliza en {tp} y la enzima produce a v₀ = {v0:.3g} µM/s":
        "ES settles within {tp} and the enzyme produces at v₀ = {v0:.3g} µM/s",
    "La hipótesis del estado estacionario (ES constante) es la base de Michaelis–Menten.":
        "The steady-state assumption (constant ES) is the basis of Michaelis–Menten.",
    "Km = (k₋₁ + k₂)/k₁ = {km:.3g} µM; kcat = k₂ = {k2:.3g} s⁻¹; Vmax = kcat·[E]₀ = {vmax:.3g} µM/s. "
    "Con [S]₀ = {s0:g} µM ({ratio:.2g}·Km) la v₀ simulada es {v0:.3g} µM/s "
    "({pct:.0f} % de Vmax); el estado estacionario de ES se alcanza en ≈ {tp}.":
        "Km = (k₋₁ + k₂)/k₁ = {km:.3g} µM; kcat = k₂ = {k2:.3g} s⁻¹; Vmax = kcat·[E]₀ = {vmax:.3g} µM/s. "
        "With [S]₀ = {s0:g} µM ({ratio:.2g}·Km) the simulated v₀ is {v0:.3g} µM/s "
        "({pct:.0f} % of Vmax); the ES steady state is reached in ≈ {tp}.",
    "k₁: el sustrato entra (µM⁻¹·s⁻¹)": "k₁: substrate binds (µM⁻¹·s⁻¹)",
    "k₋₁: el sustrato se suelta (s⁻¹)": "k₋₁: substrate lets go (s⁻¹)",
    "k₂ = kcat: reacciona (s⁻¹)": "k₂ = kcat: it reacts (s⁻¹)",
    "[E]₀: enzima en el tubo (µM)": "[E]₀: enzyme in the tube (µM)",
    "[S]₀: sustrato inicial (µM)": "[S]₀: initial substrate (µM)",
    "El mecanismo E + S ⇌ ES → E + P, paso a paso": "The mechanism E + S ⇌ ES → E + P, step by step",
    "Con n = {n:.2g}, basta multiplicar [S] por {f:.1f} para pasar del 10 % al 90 %":
        "With n = {n:.2g}, multiplying [S] by {f:.1f} is enough to go from 10 % to 90 %",
    "Franja naranja: la ventana 10 %–90 % de la sigmoide. Con n = 1 (hipérbola) haría falta ×81.":
        "Orange band: the 10 %–90 % window of the sigmoid. With n = 1 (hyperbola) it would take ×81.",
    "Con S₀.₅ = {sh:g} mM y n = {n:.2g}, pasar del 10 % al 90 % de activación exige subir [S] "
    "de {s10:.3g} a {s90:.3g} mM (×{f:.1f}); con n = 1 (hipérbola) haría falta ×81.":
        "With S₀.₅ = {sh:g} mM and n = {n:.2g}, going from 10 % to 90 % activation requires raising [S] "
        "from {s10:.3g} to {s90:.3g} mM (×{f:.1f}); with n = 1 (hyperbola) it would take ×81.",
    "S₀.₅: [S] a media actividad (mM)": "S₀.₅: [S] at half activity (mM)",
    "n: coeficiente de Hill": "n: Hill coefficient",
    "Interruptor (sigmoide) frente a regulador (hipérbola)": "Threshold switch (sigmoid) versus dimmer (hyperbola)",
    "competitiva": "competitive",
    "acompetitiva": "uncompetitive",
    "no competitiva": "noncompetitive",
    "mixta": "mixed",
    "Vmax sin inhibidor": "Vmax without inhibitor",
    "Km aparente": "apparent Km",
    "v₀ frente a [S]: puntos = (Km, Vmax/2)": "v₀ versus [S]: points = (Km, Vmax/2)",
    "Su huella en Lineweaver–Burk": "Its fingerprint in Lineweaver–Burk",
    "{nombre}: Vmax ×{rv:.2f}, Km aparente ×{rk:.2f}": "{nombre}: Vmax ×{rv:.2f}, apparent Km ×{rk:.2f}",
    "En Lineweaver–Burk {patron}.": "In Lineweaver–Burk, {patron}.",
    "{nombre} con [I] = {i:g} mM y {ki_text}: "
    "Vmax_app = {va:.3g} µM/s (×{rv:.2f}), Km_app = {ka:.3g} mM (×{rk:.2f}). "
    "En Lineweaver–Burk {patron}.":
        "{nombre} with [I] = {i:g} mM and {ki_text}: "
        "Vmax_app = {va:.3g} µM/s (×{rv:.2f}), Km_app = {ka:.3g} mM (×{rk:.2f}). "
        "In Lineweaver–Burk, {patron}.",
    "Tipo de inhibidor": "Inhibitor type",
    "[I]: cantidad de inhibidor (mM)": "[I]: amount of inhibitor (mM)",
    "Ki: potencia (menor = más potente, mM)": "Ki: potency (smaller = more potent, mM)",
    "Ki′ (mM, solo en la mixta)": "Ki′ (mM, mixed type only)",
    "Cuatro formas de frenar una enzima": "Four ways to slow down an enzyme",
    "k_cat de la glucoquinasa (~60 s⁻¹)": "glucokinase k_cat (~60 s⁻¹)",
    "una vez por hora": "once per hour",
    "una vez por año": "once per year",
    "ΔG‡: altura de la colina (kcal/mol)": "ΔG‡: height of the hill (kcal/mol)",
    "k (s⁻¹, escala log)": "k (s⁻¹, log scale)",
    "k = (kB·T/h)·exp(−ΔG‡/RT) a {tk:g} K ({tc:.0f} °C)": "k = (kB·T/h)·exp(−ΔG‡/RT) at {tk:g} K ({tc:.0f} °C)",
    "(la tuya)": "(yours)",
    "barrera ΔG‡ (kcal/mol)": "barrier ΔG‡ (kcal/mol)",
    "±{d:.2f} kcal/mol = ×10 en k": "±{d:.2f} kcal/mol = ×10 in k",
    "Con ΔG‡ = {dg:g} kcal/mol, cada molécula reacciona en {tt}": "With ΔG‡ = {dg:g} kcal/mol, each molecule reacts in {tt}",
    "Mueve ΔG‡: fíjate en que la escala vertical es logarítmica (cada línea es un factor 10).":
        "Move ΔG‡: note that the vertical scale is logarithmic (each gridline is a factor of 10).",
    "ΔG‡ = {dg:g} kcal/mol a {tk:g} K → k = {k:.3g} s⁻¹; tiempo de recambio 1/k = {t1}, "
    "vida media t½ = ln2/k = {th}. A esta temperatura cada {d:.2f} kcal/mol "
    "adicionales de barrera dividen k entre 10 (1.36 kcal/mol a 298 K).":
        "ΔG‡ = {dg:g} kcal/mol at {tk:g} K → k = {k:.3g} s⁻¹; turnover time 1/k = {t1}, "
        "half-life t½ = ln2/k = {th}. At this temperature every extra {d:.2f} kcal/mol "
        "of barrier divides k by 10 (1.36 kcal/mol at 298 K).",
    "T: temperatura (K)": "T: temperature (K)",
    "La ecuación de Eyring: de la colina a la velocidad": "The Eyring equation: from the hill to the rate",
    "Gráfico de Eyring: la pendiente es −ΔH‡/R": "Eyring plot: the slope is −ΔH‡/R",
    "cuerpo\nhumano": "human\nbody",
    "k(T): calentar acelera, y más cuanto mayor es ΔH‡": "k(T): heating speeds things up, more so the larger ΔH‡ is",
    "De 25 a 37 °C la reacción va ×{f:.2f} más rápido": "From 25 to 37 °C the reaction runs ×{f:.2f} faster",
    "ΔH‡ = {dh:g} kcal/mol decide cuánto acelera el calor; ΔS‡ = {ds} cal/(mol·K) mueve todo arriba o abajo.":
        "ΔH‡ = {dh:g} kcal/mol decides how much heat speeds it up; ΔS‡ = {ds} cal/(mol·K) shifts everything up or down.",
    "ΔH‡ = {dh:g} kcal/mol y ΔS‡ = {ds} cal/(mol·K) → ΔG‡(298 K) = ΔH‡ − TΔS‡ = {dg:.2f} kcal/mol, "
    "k(25 °C) = {k25:.3g} s⁻¹ y k(37 °C) = {k37:.3g} s⁻¹ (×{f:.2f}); "
    "Q10 (25→35 °C) = {q10:.2f}.":
        "ΔH‡ = {dh:g} kcal/mol and ΔS‡ = {ds} cal/(mol·K) → ΔG‡(298 K) = ΔH‡ − TΔS‡ = {dg:.2f} kcal/mol, "
        "k(25 °C) = {k25:.3g} s⁻¹ and k(37 °C) = {k37:.3g} s⁻¹ (×{f:.2f}); "
        "Q10 (25→35 °C) = {q10:.2f}.",
    "ΔH‡: energía a aportar (kcal/mol)": "ΔH‡: energy to supply (kcal/mol)",
    "ΔS‡: orden a imponer (cal/(mol·K))": "ΔS‡: order to impose (cal/(mol·K))",
    "Temperatura: entalpía y entropía de activación": "Temperature: activation enthalpy and entropy",
    "El óptimo está a mitad de camino entre los dos pKa: pH {p:.2f}": "The optimum lies halfway between the two pKa values: pH {p:.2f}",
    "v = Vmax / (1 + 10^(pKa₁ − pH) + 10^(pH − pKa₂)) · zonas grises: un grupo catalítico está «apagado»":
        "v = Vmax / (1 + 10^(pKa₁ − pH) + 10^(pH − pKa₂)) · gray zones: a catalytic group is \"switched off\"",
    "la base (Asp205)\nya tiene protón": "the base (Asp205)\nalready has a proton",
    "la lisina (Lys169)\nperdió su carga +": "the lysine (Lys169)\nlost its + charge",
    "glucoquinasa medida: pH 8.5–8.7": "measured glucokinase: pH 8.5–8.7",
    "óptimo: pH {p:.2f} ({pct:.0f} % de Vmax)": "optimum: pH {p:.2f} ({pct:.0f} % of Vmax)",
    "Con pKa₁ = {p1:g} y pKa₂ = {p2:g} el óptimo está en pH = (pKa₁ + pKa₂)/2 = {po:.2f}, "
    "donde la enzima alcanza el {pct:.0f} % de Vmax":
        "With pKa₁ = {p1:g} and pKa₂ = {p2:g} the optimum is at pH = (pKa₁ + pKa₂)/2 = {po:.2f}, "
        "where the enzyme reaches {pct:.0f} % of Vmax",
    " (los pKa están tan próximos que nunca se llega al máximo teórico).":
        " (the pKa values are so close that the theoretical maximum is never reached).",
    "; la actividad cae a la mitad en pH ≈ {p1:g} y ≈ {p2:g}.": "; activity drops to half at pH ≈ {p1:g} and ≈ {p2:g}.",
    "pKa₁ (grupo que debe estar sin protón)": "pKa₁ (group that must be deprotonated)",
    "pKa₂ (grupo que debe tener protón)": "pKa₂ (group that must be protonated)",
    "La campana de pH: dos grupos, dos condiciones": "The pH bell: two groups, two conditions",
    "glucosa en sangre\nen ayunas (4–7 mM)": "fasting blood\nglucose (4–7 mM)",
    "sin activador (S₀.₅ = {sh} mM, n = {n})": "without activator (S₀.₅ = {sh} mM, n = {n})",
    "con activador (Vmax ×{vf:g}, S₀.₅ = {sh:.2g} mM)": "with activator (Vmax ×{vf:g}, S₀.₅ = {sh:.2g} mM)",
    "[glucosa] (mM)": "[glucose] (mM)",
    "actividad (% de la Vmax basal)": "activity (% of basal Vmax)",
    "A 5 mM de glucosa, el activador multiplica la actividad por {f:.1f}": "At 5 mM glucose, the activator multiplies activity by {f:.1f}",
    "El activador baja S₀.₅ y sube Vmax: la curva se desplaza hacia arriba y a la izquierda.":
        "The activator lowers S₀.₅ and raises Vmax: the curve shifts up and to the left.",
    "A 5 mM de glucosa: sin activador la glucoquinasa trabaja al {a:.0f} % de su Vmax basal; "
    "con activador (Vmax ×{vf:g}, S₀.₅ ×{sf:g} → {sh:.2g} mM) trabaja al "
    "{b:.0f} % (×{f:.1f}).":
        "At 5 mM glucose: without activator glucokinase works at {a:.0f} % of its basal Vmax; "
        "with activator (Vmax ×{vf:g}, S₀.₅ ×{sf:g} → {sh:.2g} mM) it works at "
        "{b:.0f} % (×{f:.1f}).",
    "Vmax se multiplica por": "Vmax is multiplied by",
    "S₀.₅ se multiplica por": "S₀.₅ is multiplied by",
    "Un activador alostérico de la glucoquinasa": "An allosteric activator of glucokinase",
    "exergónica (ΔE < 0)": "exergonic (ΔE < 0)",
    "endergónica (ΔE > 0)": "endergonic (ΔE > 0)",
    "termoneutra": "thermoneutral",
    "cuesta abajo": "downhill",
    "cuesta arriba": "uphill",
    "a nivel": "level",
    "La colina decide la velocidad (k ≈ {k:.3g} s⁻¹); el desnivel, hacia dónde va ({c})":
        "The hill sets the rate (k ≈ {k:.3g} s⁻¹); the drop sets the direction ({c})",
    "Altura de la cima → rapidez. Diferencia entre valles → equilibrio.":
        "Height of the summit → speed. Difference between valleys → equilibrium.",
    "ΔE‡ = {b:g} kcal/mol → k(298 K, Eyring, tratando ΔE‡ como ΔG‡) = {k:.3g} s⁻¹ "
    "(1/k = {t1}); reacción {kind}: barrera inversa {br:g} kcal/mol, "
    "k inversa = {kr:.3g} s⁻¹, K_eq = exp(−ΔE/RT) = {keq:.3g}.":
        "ΔE‡ = {b:g} kcal/mol → k(298 K, Eyring, treating ΔE‡ as ΔG‡) = {k:.3g} s⁻¹ "
        "(1/k = {t1}); {kind} reaction: reverse barrier {br:g} kcal/mol, "
        "reverse k = {kr:.3g} s⁻¹, K_eq = exp(−ΔE/RT) = {keq:.3g}.",
    " (La barrera se elevó a {b:g} kcal/mol para que el TS quede por encima de los productos.)":
        " (The barrier was raised to {b:g} kcal/mol so that the TS lies above the products.)",
    "ΔE‡: altura de la colina (kcal/mol)": "ΔE‡: height of the hill (kcal/mol)",
    "ΔE: desnivel entre valles (kcal/mol)": "ΔE: drop between valleys (kcal/mol)",
    "Dibuja tu propia colina de energía": "Draw your own energy hill",
    "Mecanismo (EDO)": "Mechanism (ODE)",
    "Temperatura": "Temperature",
    "Activador (GKA)": "Activator (GKA)",
    "Perfil de energía": "Energy profile",

    # ------------------------------------------------------------------ visores 3D (Python)
    "O6 (glucosa)": "O6 (glucose)",
    "Estado de transición: el fósforo está a medio camino entre los dos oxígenos":
        "Transition state: the phosphorus is halfway between the two oxygens",
    "Reactivo → cima: el fosfato se separa del ATP y se acerca a la glucosa":
        "Reactant → summit: the phosphate leaves ATP and approaches glucose",
    "Después de la cima: el fosfato ya está en la glucosa; el protón de O6 aún no se ha movido":
        "After the summit: the phosphate is already on glucose; the O6 proton has not moved yet",
    "Hacia el producto: el fosfato está en la glucosa y Asp205 empieza a tomar el protón de O6":
        "Toward the product: the phosphate is on glucose and Asp205 starts taking the O6 proton",
    "La reacción dentro de la enzima": "The reaction inside the enzyme",
    "La reacción en agua (sin el resto de la enzima)": "The reaction in water (without the rest of the enzyme)",
    "Arrastra para girar · rueda para acercar · ▶ para reproducir": "Drag to rotate · scroll to zoom · ▶ to play",
    "Arrastra para girar · rueda para acercar · botones para cambiar la vista":
        "Drag to rotate · scroll to zoom · buttons to change the view",
    "enlace que se rompe": "breaking bond",
    "enlace que se forma": "forming bond",
    "Mg²⁺ coordinado (no covalente)": "coordinated Mg²⁺ (non-covalent)",
    "El modo imaginario del estado de transición ({f:.0f} cm⁻¹)": "The imaginary mode of the transition state ({f:.0f} cm⁻¹)",
    "No es una vibración normal: hacia un lado cae al reactivo, hacia el otro al producto":
        "Not a normal vibration: one way it falls to the reactant, the other way to the product",
    "un ciclo de la vibración": "one vibration cycle",
    "desplazamiento en ξ (Å)": "displacement in ξ (Å)",
    "¿Hacia dónde empuja el modo?": "Which way does the mode push?",
    "hacia el producto": "toward the product",
    "hacia el reactivo": "toward the reactant",
    "ξ relativo": "relative ξ",
    "El fósforo se acerca al O6 de la glucosa: hacia el producto": "The phosphorus approaches glucose O6: toward the product",
    "El fósforo vuelve hacia el O3β del ATP: hacia el reactivo": "The phosphorus goes back toward ATP O3β: toward the reactant",
    "La película de la enzima: 1 ns de dinámica molecular": "The enzyme movie: 1 ns of molecular dynamics",
    "¿Se apuntan los reactivos?": "Are the reactants pointing at each other?",
    "ataque cercano (< {u} Å)": "near attack (< {u} Å)",
    "cristal {c} Å": "crystal {c} Å",
    "Lista para reaccionar: el fósforo del ATP apunta al O6 de la glucosa": "Ready to react: the ATP phosphorus points at glucose O6",
    "Un momento de separación: los reactivos se alejan un poco": "A moment apart: the reactants drift away a little",
    "La glucoquinasa humana (cristal 3FGU)": "Human glucokinase (crystal 3FGU)",
    "AMP‑PNP (análogo del ATP)": "AMP‑PNP (ATP analog)",
    "Asp205 · base": "Asp205 · base",
    "Pγ (ATP) ··· O6 (glucosa)": "Pγ (ATP) ··· O6 (glucose)",
    "O6 ··· OD1 (Asp205)": "O6 ··· OD1 (Asp205)",
    "QM/MM: la región cuántica dentro de la enzima": "QM/MM: the quantum region inside the enzyme",
    "Bolas y varillas = región QM (PM7) · esferas rojas/azules = cargas de la enzima (−/+)":
        "Balls and sticks = QM region (PM7) · red/blue spheres = enzyme charges (−/+)",
    "carga negativa de la enzima": "negative enzyme charge",
    "carga positiva de la enzima": "positive enzyme charge",

    # ------------------------------------------------------------------ visores 3D (plantilla JS)
    "No se pudo cargar 3Dmol.js (¿sin conexión a internet?). Vuelve a ejecutar la celda.":
        "Could not load 3Dmol.js (no internet connection?). Run the cell again.",
    "Toda la enzima": "Whole enzyme",
    "Sitio activo": "Active site",
    "Superficie": "Surface",
    "Girar": "Spin",
    "Distancias": "Distances",
    "<div style='font-weight:650;color:#e8ecf3;font-size:15px;margin-bottom:6px'>Qué estás viendo</div>":
        "<div style='font-weight:650;color:#e8ecf3;font-size:15px;margin-bottom:6px'>What you are seeing</div>",
    "▶  Reproducir": "▶  Play",
    "⏸  Pausa": "⏸  Pause",
    "Etiquetas": "Labels",
    "Proteína": "Protein",
    "Centrar": "Center",
    "Error en el visor 3D: ": "3D viewer error: ",
    "Las cintas son la cadena de la proteína: <b style='color:#9fb6de'>dominio grande</b> y "
    "<b style='color:#efb3d1'>dominio pequeño</b>. ":
        "The ribbons are the protein chain: <b style='color:#9fb6de'>large domain</b> and "
        "<b style='color:#efb3d1'>small domain</b>. ",
    "En la hendidura entre ambos están la <b style='color:#ffb38a'>glucosa</b> y el <b style='color:#ffd27a'>ATP</b> "
    "(aquí AMP‑PNP), con el <b style='color:#3ddc84'>Mg²⁺</b>.<br><br>":
        "In the cleft between them sit <b style='color:#ffb38a'>glucose</b> and <b style='color:#ffd27a'>ATP</b> "
        "(here AMP‑PNP), with <b style='color:#3ddc84'>Mg²⁺</b>.<br><br>",
    "Pulsa <b>Sitio activo</b> para acercarte: verás Asp205, la base que tomará el protón de la glucosa. ":
        "Press <b>Active site</b> to zoom in: you will see Asp205, the base that will take the glucose proton. ",
    "<b>Superficie</b> muestra la forma de la proteína: la glucosa queda casi enterrada.<br><br>":
        "<b>Surface</b> shows the shape of the protein: glucose ends up almost buried.<br><br>",
    "Las <b>líneas discontinuas</b> son las distancias medidas en el cristal (las mismas que calcula la celda siguiente): ":
        "The <b>dashed lines</b> are the distances measured in the crystal (the same ones the next cell calculates): ",
    "En <b>bolas y varillas</b>, los 78 átomos que se tratan con mecánica cuántica (PM7): glucosa, trifosfato, "
    "Mg²⁺ con su agua ":
        "In <b>balls and sticks</b>, the 78 atoms treated with quantum mechanics (PM7): glucose, triphosphate, "
        "Mg²⁺ with its water ",
    "y las cadenas laterales de Asp205, Lys169 y Thr228.<br><br>Las <b style='color:#ff8a8a'>esferas rojas</b> "
    "y <b style='color:#8ab4ff'>azules</b> son ":
        "and the side chains of Asp205, Lys169 and Thr228.<br><br>The <b style='color:#ff8a8a'>red</b> "
        "and <b style='color:#8ab4ff'>blue spheres</b> are ",
    "cargas negativas y positivas de la enzima (región MM): sus electrones no se calculan, pero su campo "
    "eléctrico sí entra en el cálculo cuántico.":
        "negative and positive charges of the enzyme (MM region): their electrons are not calculated, but their "
        "electric field does enter the quantum calculation.",

    # ------------------------------------------------------------------ instalación (colab_setup)
    "instalando:": "installing:",
    "descargando MOPAC": "downloading MOPAC",
    "entorno listo:": "environment ready:",

    # ------------------------------------------------------------------ valores de referencia (kinetics.glucokinase_reference)
    "glucoquinasa humana (hexoquinasa IV, GCK)": "human glucokinase (hexokinase IV, GCK)",
    "S₀.₅ (glucosa)": "S₀.₅ (glucose)",
    "silvestre, 0-100 mM glucosa, 5 mM ATP": "wild type, 0-100 mM glucose, 5 mM ATP",
    "coeficiente de Hill n (glucosa)": "Hill coefficient n (glucose)",
    "adimensional": "dimensionless",
    "silvestre": "wild type",
    "silvestre; 62.3 s⁻¹ en Sayed 2009; 38 s⁻¹ en Heredia 2006 (stopped-flow)":
        "wild type; 62.3 s⁻¹ in Sayed 2009; 38 s⁻¹ in Heredia 2006 (stopped-flow)",
    "ΔG‡ que implica k_cat (Eyring, κ = 1)": "ΔG‡ implied by k_cat (Eyring, κ = 1)",
    "derivado de k_cat = 62-66 s⁻¹ a 25 °C; 15.6 a 37 °C": "derived from k_cat = 62-66 s⁻¹ at 25 °C; 15.6 at 37 °C",
    "aritmética nuestra": "our arithmetic",
    "pH óptimo": "optimal pH",
    "8.5-8.7 tras corregir la acidificación por ATP": "8.5-8.7 after correcting for acidification by ATP",
    "IC₅₀ de GKRP (frente a glucosa)": "GKRP IC₅₀ (versus glucose)",
    "mM glucosa": "mM glucose",
    "k de intercambio conformacional (abierta ⇌ cerrada)": "conformational exchange k (open ⇌ closed)",
    "entre 5 y 100 s⁻¹ por RMN, escala de ms; comparable a k_cat": "between 5 and 100 s⁻¹ by NMR, ms timescale; comparable to k_cat",
    "ΔE‡ QM/MM publicada (silvestre)": "published QM/MM ΔE‡ (wild type)",
    "Asp205 base general; K169A: 32.1 kcal/mol": "Asp205 general base; K169A: 32.1 kcal/mol",
    "vida media de la hidrólisis no catalizada de un fosfato dianión": "half-life of the uncatalyzed hydrolysis of a phosphate dianion",
    "k = 2×10⁻²⁰ s⁻¹ a 25 °C (ΔG‡ ≈ 44 kcal/mol); hidrólisis, no transferencia a glucosa":
        "k = 2×10⁻²⁰ s⁻¹ at 25 °C (ΔG‡ ≈ 44 kcal/mol); hydrolysis, not transfer to glucose",
    "N-acetilglucosamina (GlcNAc)": "N-acetylglucosamine (GlcNAc)",
    "no localizada en fuente primaria": "not found in a primary source",
    "inhibidor de alta afinidad; competitivo y baja n hacia 1.0": "high-affinity inhibitor; competitive and lowers n toward 1.0",
    "manoheptulosa": "mannoheptulose",
    "discutida (0.25-20 mM según fuentes secundarias)": "disputed (0.25-20 mM depending on secondary sources)",
    "tipo mixto, predominantemente no competitivo/competitivo según Scruel 1998":
        "mixed type, predominantly noncompetitive/competitive according to Scruel 1998",
    "glucosamina": "glucosamine",
    "baja afinidad; también sustrato": "low affinity; also a substrate",
    "glucosa-6-fosfato (producto)": "glucose-6-phosphate (product)",
    "sin inhibición fisiológica (a diferencia de HK I-III)": "no physiological inhibition (unlike HK I-III)",
    "palmitoil-CoA": "palmitoyl-CoA",
    "15 µM → >90 % inactivación (horas)": "15 µM → >90 % inactivation (hours)",
    "proteína reguladora GKRP": "regulatory protein GKRP",
    "IC₅₀ 13.3 mM glucosa; K_i no localizada": "IC₅₀ 13.3 mM glucose; K_i not found",
    "sube S₀.₅ sin cambiar V; F6P refuerza (k_off ÷60), F1P antagoniza":
        "raises S₀.₅ without changing V; F6P strengthens it (k_off ÷60), F1P antagonizes it",
    "RO-28-1675 (RO0281675, GKA)": "RO-28-1675 (RO0281675, GKA)",
    "baja": "lowers",
    "sube": "raises",
    "activador mixto no esencial; sitio alostérico a ~20 Å del sitio de glucosa":
        "nonessential mixed activator; allosteric site ~20 Å from the glucose site",
    "activadora (hiperinsulinismo)": "activating (hyperinsulinism)",
    "GCK-MODY: Valentínová et al. 2012 (tabla 3); activadoras: Sayed et al. 2009 (tabla 2). "
    "k_cat relativo a la silvestre del mismo estudio.":
        "GCK-MODY: Valentínová et al. 2012 (table 3); activating: Sayed et al. 2009 (table 2). "
        "k_cat relative to the wild type of the same study.",
    "Valores de fuentes primarias (ver 'source' y 'references'); varían entre laboratorios y ensayos "
    "(acoplados a G6PDH, proteínas de fusión GST, temperatura no siempre indicada). Todos los números son "
    "aproximados (approx) y deben verificarse contra las fuentes primarias antes de usarse cuantitativamente.":
        "Values from primary sources (see 'source' and 'references'); they vary between laboratories and assays "
        "(G6PDH-coupled, GST fusion proteins, temperature not always stated). All numbers are "
        "approximate (approx) and should be checked against the primary sources before quantitative use.",
    # textos que vienen de los datos precalculados (resumen.json, instantaneas.json) y de las referencias
    "Aproximación armónica sobre los 70 átomos libres de la región QM con el entorno fijo; modos < 50 cm-1 elevados a 50 cm-1; sin contribuciones del entorno ni de la conformación de la proteína.":
        "Harmonic approximation over the 70 free atoms of the QM region with the environment fixed; modes < 50 cm-1 raised to 50 cm-1; no contributions from the environment or from the protein conformation.",
    "Energías de punto único en las geometrías PM7 (R, TS, P), mismo campo MM y mismo término LJ.":
        "Single-point energies at the PM7 geometries (R, TS, P), same MM field and same LJ term.",
    "dímero (ASE) + TS de MOPAC": "dimer (ASE) + MOPAC TS",
    "cristal minimizado": "minimized crystal",
    "punto": "point",
    "Senn H. M. & Thiel W. (2009) Angew. Chem. Int. Ed. 48:1198-1229 (métodos QM/MM)":
        "Senn H. M. & Thiel W. (2009) Angew. Chem. Int. Ed. 48:1198-1229 (QM/MM methods)",
}
