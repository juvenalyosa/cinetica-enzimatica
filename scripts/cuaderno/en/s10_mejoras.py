"""Section 10. Improvements."""
from ..celdas import code, md

# ============================================================================ 10. improvements

md(r"""
## 10. From ΔE‡ to ΔG‡: three improvements to the calculation

> 🎯 **In this section** you will turn the electronic barrier ΔE‡ into the free energy ΔG‡ that
> the Eyring equation asks for, and learn to ask yourself how much to trust a computed number.

---

### 💡 The analogy: the careful surveyor

A surveyor who has measured the height of a mountain **once, on one day, with one instrument**
asks three questions before publishing the number: does the ground move? would I have measured
something else on another day? is my instrument properly calibrated?

So far we have an **electronic energy** ΔE‡, computed in **one** conformation of the enzyme with
**one** semiempirical method. The Eyring equation asks for a **free energy** ΔG‡. A good
computational scientist asks the same three questions, and answers them with extra calculations:

[[fig:tres_mejoras | Three questions: vibrations, other poses of the enzyme and another method; the sum that takes ΔE‡ to ΔG‡]]

---

### 🧮 The equation in words

**Why do we want an equation?** Because the "bare" height (ΔE‡) ignores that atoms are always
vibrating. The free energy adds what it costs, or saves, in vibration and in order:

$$\text{free barrier} \;=\; \text{electronic barrier} \;+\; \text{vibration at rest} \;+\; \text{thermal vibration} \;-\; T\times\text{(change in disorder)}$$

---

### 📐 The full equation, term by term

$$\Delta G^{\ddagger} = \Delta E^{\ddagger} + {\color{#1baf7a}{\Delta \mathrm{ZPE}}} + {\color{#4a3aa7}{\Delta H_{\mathrm{vib}}(T)}} - {\color{#e87ba4}{T\,\Delta S_{\mathrm{vib}}(T)}}$$

| Term | What it is | In the analogy | Our value (25 °C) |
|---|---|---|---|
| **ΔE‡** | electronic height of the hill, atoms standing still | the "one-day" measurement | 19.1 kcal/mol |
| **ΔZPE** | change in the **zero-point energy**: even at 0 K atoms vibrate | the ground trembles even without wind | −0.2 kcal/mol |
| **ΔH<sub>vib</sub>(T)** | extra vibrational energy at temperature T | the extra trembling on a hot day | −0.4 kcal/mol |
| **−TΔS<sub>vib</sub>** | the price of order: a TS stiffer than the reactant has a negative ΔS‡ | tidying a room takes work | +1.3 kcal/mol |

* **ΔZPE**: the P–O bond that breaks loses its stretching vibration at the TS, so it is usually
  **negative** (it lowers the barrier).
* **−TΔS‡**: if the TS is "stiffer" than the reactant (ΔS‡ < 0), this term **raises** the barrier.

**The other two questions do not change the formula, but our confidence in the number:**

**2. What if the enzyme were in another pose?** Since the environment is fixed, each MD snapshot
can give a different barrier. We repeat the calculation on several snapshots (250, 500, 750 and
1000 ps) to see **whether the barrier changes and why**. A good result is not a number: it is a
number with its uncertainty, and sometimes the answer is that the reaction is not even possible
from some conformations.

**3. What if the method is wrong?** PM7 is semiempirical. We recompute the energies with another
Hamiltonian (PM6‑D3H4) on the same geometries: if the barrier changes little, the result is
robust; if it changes a lot, we know how much to (dis)trust it.

---

### 🎛️ What if…

| If… | then… | because… |
|---|---|---|
| **the temperature goes from 25 to 37 °C** | ΔG‡ barely rises (from 19.75 to 19.80 kcal/mol) | with ΔS‡ ≈ −4.3 cal/mol/K, 12 K more only add 12 × 0.0043 ≈ 0.05 kcal/mol |
| **the TS were more flexible than the reactant** (ΔS‡ > 0) | the entropic term would **lower** the barrier | −TΔS‡ would change sign |
| **we change the method** (PM7 → PM6‑D3H4) | the barrier goes from 19.1 to 4.5 kcal/mol | semiempirical methods are calibrated against specific molecules |

---

### 🔬 The real data
""")

code(r'''
# @title 🎻 The vibrations: from ΔE‡ to ΔG‡
termo = res["etapas"].get("termoquimica", {})
T_ref = "298.15"
if T_ref in termo:
    t = termo[T_ref]
    extra = [f"at {float(Tk) - 273.15:.0f} °C: ΔG‡ = {termo[Tk]['dG_kcal']:.2f} kcal/mol" for Tk in ("303.15", "310.15") if Tk in termo]
    display(viz.tarjetas(
        [("Electronic ΔE‡", f"{t['dE_kcal']:.2f}", "kcal/mol", "the hill with the atoms standing still", "gris"),
         ("ΔZPE", f"{t['dZPE_kcal']:+.2f}", "kcal/mol", "the breaking P–O spring stops vibrating", "agua"),
         ("Vibrational ΔS‡", f"{t['dS_vib_cal']:+.1f}", "cal/mol/K", f"stiffer TS → −TΔS‡ = {-298.15 * t['dS_vib_cal'] / 1000:+.2f} kcal/mol", "magenta"),
         ("ΔG‡ at 25 °C", f"{t['dG_kcal']:.2f}", "kcal/mol", f"ΔH‡ = {t['dH_kcal']:.2f} kcal/mol", "naranja")],
        titulo="Harmonic thermochemistry (free atoms of the QM region)",
        nota=" · ".join(extra) + ("  —  " + textos.t(termo["nota"]) if termo.get("nota") else "")))
else:
    display(viz.mensaje("Thermochemistry not available in these data.", "ojo"))
''')

code(r'''
# @title 📸 The barrier in different MD conformations
try:
    inst = datos.json_("qmmm/instantaneas.json")
    def _hay(v):
        return v is not None and np.isfinite(v)
    filas = []
    for r in inst["instantaneas"]:
        b = r.get("barrera_kcal")
        if _hay(b):
            que = f"TS found ({r.get('ts_origen', 'dimer')}): barrier {b:.1f} kcal/mol, stable product exists"
        else:
            que = (f"no stable product: the energy rises to {r['escaneo_max_rel_kcal']:.0f} kcal/mol at ξ = {r['xi_max']:.1f} Å; "
                   f"the O6 H is {r['d_OD1_H_reactivo']:.1f} Å away from Asp205")
        filas.append((textos.t(r["instantanea"]), b if _hay(b) else np.nan, r.get("ts_d_PG_O6", np.nan), r.get("ts_d_PG_O3B", np.nan), que))
    df_inst = pd.DataFrame(filas, columns=["conformation", "ΔE‡ (kcal/mol)", "d(Pγ–O6) TS (Å)", "d(Pγ–O3β) TS (Å)", "what happened?"])
    ok = df_inst.dropna(subset=["ΔE‡ (kcal/mol)"])
    media, sd = ok["ΔE‡ (kcal/mol)"].mean(), ok["ΔE‡ (kcal/mol)"].std(ddof=1) if len(ok) > 1 else 0.0

    esc_inst = datos.csv("qmmm/instantaneas_escaneos.csv")
    fallidas = [(g["xi"].values, g["energia_rel_kcal"].values, nombre) for nombre, g in esc_inst.groupby("instantanea", sort=False)]
    fig = viz.plot_snapshot_scans((xi_esc, e_esc), fallidas,       # the same crystal scan as in section 8
                                  title="From the wrong conformation, the reaction has nowhere to go",
                                  subtitle="ξ scans. Blue: the minimized crystal (summit and product valley). Red: four MD frames in which Asp205 moved away")
    viz.mostrar(fig, viz.tabla(df_inst, titulo="What happened in each conformation",
                               formatos={"ΔE‡ (kcal/mol)": "{:.1f}", "d(Pγ–O6) TS (Å)": "{:.2f}", "d(Pγ–O3β) TS (Å)": "{:.2f}"},
                               nota="\"—\": there is no transition state to measure, because without a stable product there is no summit between two valleys."),
                viz.datos(esc_inst.pivot_table(index="punto", columns="instantanea", values="energia_rel_kcal", sort=False)
                          .reset_index().rename(columns={"punto": textos.t("punto")}).assign(**{"ξ (mean)": esc_inst.groupby("punto")["xi"].mean().values})
                          [[textos.t("punto"), "ξ (mean)"] + list(dict.fromkeys(esc_inst["instantanea"]))],
                          "the red curves",
                          "Each row is a point of the scan (ξ fixed further and further to the right); each column, one MD "
                          "conformation. In no column does the energy reach a maximum and come back down: there is no summit and no product valley.",
                          x="ξ (mean)", y=list(dict.fromkeys(esc_inst["instantanea"])),
                          unidades={"ξ (mean)": "Å", **{c: "kcal/mol" for c in esc_inst["instantanea"].unique()}},
                          formatos={"ξ (mean)": "{:+.2f}", **{c: "{:.1f}" for c in esc_inst["instantanea"].unique()}},
                          nota="ξ varies by a few hundredths between conformations; its mean is shown."))
    if len(ok) > 1:
        display(viz.mensaje(f"Mean barrier = {media:.1f} ± {sd:.1f} kcal/mol (n = {len(ok)})", "dato"))
except FileNotFoundError:
    display(viz.mensaje("Average over snapshots not available in these data.", "ojo"))
''')

md(r"""
**What the MD snapshots teach us.** In the minimized crystal conformation, the O6–H hydroxyl of
glucose is hydrogen-bonded to Asp205 (1.7 Å) and the path has a stable product. In several
snapshots of the classical dynamics, however, Asp205 has moved away from O6–H (in the MD the
O6···OD1 distance averages 4.3 Å, section 5): when the phosphoryl transfer is forced, the proton
**does not find the base** and the energy keeps rising, with no product minimum. Two lessons:
(1) the reaction can only happen from the right **near-attack conformations**, which are a
fraction of the time; (2) the classical force field does not always preserve the catalytic
geometry (here, the hydroxyl–carboxylate interaction next to Mg²⁺), so the modeller must compare
the MD with the crystal before choosing the starting structure for QM/MM. An "average" barrier only
makes sense over reactive conformations.
""")

code(r'''
# @title ⚖️ And with another method?
met = res["etapas"].get("metodos", {})
colores = {"PM7": "azul", "PM6-D3H4": "magenta"}
items = [(f"{nombre}//PM7", f"{d['barrera_kcal']:.1f}", "kcal/mol barrier",
          f"ΔE(reaction) = {d['dE_reaccion_kcal']:+.1f} kcal/mol", colores.get(nombre, "gris"))
         for nombre, d in met.items() if isinstance(d, dict)]
viz.tarjetas(items, titulo="The same geometries, two semiempirical Hamiltonians", nota=textos.t(met.get("nota", "")))
''')

md(r"""
**What these numbers say.** The vibrational correction is small (ΔZPE ≈ −0.2 and −TΔS‡ ≈ +1.3
kcal/mol partly cancel): ΔG‡ ≈ ΔE‡ + 0.6 kcal/mol. The average over snapshots tells how much the
barrier "breathes" with the conformation of the enzyme. And changing the Hamiltonian is the
surprise: PM6‑D3H4 gives very different energies on the same geometries. Semiempirical methods are
parametrized against specific sets of molecules, and the trio phosphate + Mg²⁺ + carboxylate is a
hard case for them (hypervalent phosphorus, high charges, a divalent cation). The lesson:
**before believing a barrier, you need to know what chemistry the method was calibrated for**
and, if possible, compare with a higher-level calculation (DFT) or with experiment (section 12);
here PM7 turns out to be consistent with another published QM/MM calculation and with experiment,
and PM6‑D3H4 does not.
""")

md(r"""
> ✅ **Takeaway.** ΔG‡ = ΔE‡ + ΔZPE + ΔH_vib − TΔS_vib: the vibrations correct little
> (19.1 → 19.75 kcal/mol, because ΔZPE and −TΔS‡ nearly cancel). The other two questions do not
> change the number but our **confidence** in it: the barrier only makes sense from reactive
> conformations, and a method poorly calibrated for phosphate + Mg²⁺ can give a very different value.
""")
