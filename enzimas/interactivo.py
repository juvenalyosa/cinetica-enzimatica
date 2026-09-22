"""Exploradores interactivos (``ipywidgets``) para el cuaderno *Cinética Enzimática*.

Cada función ``explorar_*`` construye deslizadores en español, dibuja la figura
con el estilo de :mod:`enzimas.viz` y devuelve un ``ipywidgets.VBox`` (la celda
lo muestra si es la última expresión; también puede pasarse a ``display``).
Debajo de la figura se imprime una línea que explica qué significan los
números con los valores actuales.

Convenciones
------------
* Los identificadores son en inglés; etiquetas, descripciones y resúmenes en
  español.
* Cada explorador tiene una función de dibujo ``_dibujar_*(...)`` pura que
  devuelve ``(figura, resumen)`` y que las pruebas pueden llamar sin frontend.
* Semántica ``%matplotlib inline``: en cada actualización se crea una figura
  nueva, se muestra con ``display(fig)`` y se cierra con ``plt.close(fig)``.
  Con ``%matplotlib widget`` (ipympl) el cierre destruye el lienzo; use el
  backend ``inline``.
* Sin frontend (``nbconvert --execute``, pruebas): ``interactive_output``
  ejecuta la función una vez con los valores por defecto y captura la figura
  dentro del ``Output`` del widget; no se lanza ninguna excepción. En HTML
  estático la figura solo se ve si se guarda el estado de los widgets.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

try:
    import ipywidgets as widgets
except ImportError as exc:  # pragma: no cover - depende del entorno
    raise ImportError(
        "enzimas.interactivo necesita ipywidgets (pip install ipywidgets); "
        "en Google Colab ya viene instalado."
    ) from exc

from IPython import get_ipython
from IPython.display import display

from . import kinetics as kin
from . import viz
from .viz import AXIS, COLORS, INK, INK_MUTED, INK_SECONDARY, PALETTE, figure

__all__ = [
    "explorar_michaelis_menten",
    "explorar_mecanismo",
    "explorar_hill",
    "explorar_inhibicion",
    "explorar_eyring",
    "explorar_temperatura",
    "explorar_ph",
    "explorar_activador",
    "explorar_perfil_energia",
    "todo",
]

_T_REF = 298.15  # K
_GLUCOSE_BLOOD_MM = 5.0
_GK_S_HALF, _GK_N = 7.5, 1.7  # glucoquinasa humana (orden de magnitud)


# ---------------------------------------------------------------------------
# Infraestructura común
# ---------------------------------------------------------------------------
def _slider(description, value, lo, hi, step, fmt=".2f", log=False):
    """Deslizador de coma flotante con etiqueta completa y redibujado al soltar."""
    common = dict(
        value=value, description=description, continuous_update=False, readout_format=fmt,
        style={"description_width": "initial"}, layout=widgets.Layout(width="360px"),
    )
    if log:
        return widgets.FloatLogSlider(base=10, min=np.log10(lo), max=np.log10(hi), step=step, **common)
    return widgets.FloatSlider(min=lo, max=hi, step=step, **common)


def _mostrar(fig):
    """``display(fig)`` solo si hay un intérprete de IPython (evita imprimir el ``repr`` en consola)."""
    if get_ipython() is not None:
        display(fig)


def _explorer(update_fn, controls, titulo=None):
    """Conecta ``controls`` (``dict`` nombre → widget) con ``update_fn(**valores) -> (fig, resumen)``.

    Usa ``ipywidgets.interactive_output`` (que aporta el ``Output`` y el
    ``clear_output(wait=True)``): en cada cambio de valor se crea una figura
    nueva, se muestra, se imprime el resumen y se cierra la figura. Devuelve un
    ``VBox`` con atributos ``controles`` (el ``dict``), ``salida`` (``Output``),
    ``estado`` (último resumen y último error) y ``redibujar()``.
    """
    estado = {"error": None, "resumen": None}

    def _redibujar(**valores):
        fig = None
        try:
            fig, resumen = update_fn(**valores)
            _mostrar(fig)
            estado["resumen"] = resumen
            estado["error"] = None
            if resumen:
                print(resumen)
        except Exception as exc:
            estado["error"] = exc
            raise
        finally:
            if fig is not None:
                plt.close(fig)

    salida = widgets.interactive_output(_redibujar, controls)
    panel = widgets.HBox(list(controls.values()), layout=widgets.Layout(flex_flow="row wrap"))
    hijos = [panel, salida]
    if titulo:
        hijos.insert(0, widgets.HTML(f"<b>{titulo}</b>"))
    caja = widgets.VBox(hijos)
    caja.controles = controls
    caja.salida = salida
    caja.estado = estado
    caja.redibujar = lambda: _redibujar(**{name: w.value for name, w in controls.items()})
    return caja


def _tiempo(seconds):
    """Tiempo legible (ns … años) para tiempos de recambio y vidas medias."""
    s = float(seconds)
    if not np.isfinite(s):
        return "∞"
    scales = [
        (1e-9, "ns", 1e9), (1e-6, "µs", 1e6), (1e-3, "ms", 1e3), (1.0, "s", 1.0),
        (60.0, "min", 1 / 60.0), (3600.0, "h", 1 / 3600.0), (86400.0, "días", 1 / 86400.0),
        (3.15576e7, "años", 1 / 3.15576e7),
    ]
    chosen = scales[0]
    for threshold, unit, factor in scales:
        if s >= threshold:
            chosen = (threshold, unit, factor)
    _, unit, factor = chosen
    return f"{s * factor:.3g} {unit}"


def _num(value, decimals=1):
    return viz._num(value, decimals)


def _panel_title(ax, text):
    ax.set_title(text, loc="left", fontsize=11, color=INK_SECONDARY, fontweight="normal", pad=8)


# ---------------------------------------------------------------------------
# 1. Michaelis–Menten
# ---------------------------------------------------------------------------
def _dibujar_michaelis_menten(vmax, km, s_max=60.0):
    s = np.linspace(0.0, s_max, 400)
    v = kin.michaelis_menten(s, vmax, km)
    fig, (ax1, ax2) = figure(11.0, 4.2, ncols=2)

    ax1.plot(s, v, color=COLORS["michaelis_menten"], linewidth=2.0, zorder=2)
    ax1.set_xlim(0, s_max)
    ax1.set_ylim(0, vmax * 1.18)
    viz._guide(ax1, "h", vmax)
    ax1.annotate("Vmax", xy=(1.0, vmax), xycoords=("axes fraction", "data"), xytext=(0, 3),
                 textcoords="offset points", ha="right", va="bottom", color=INK_MUTED, fontsize=9.5)
    viz._guide(ax1, "h", vmax / 2, start=0, end=km)
    ax1.annotate("Vmax/2", xy=(km, vmax / 2), xytext=(8, -3), textcoords="offset points",
                 ha="left", va="top", color=INK_MUTED, fontsize=9.5)
    viz._guide(ax1, "v", km, start=0, end=vmax / 2)
    ax1.annotate("Km", xy=(km, 0), xytext=(4, 3), textcoords="offset points",
                 ha="left", va="bottom", color=INK_MUTED, fontsize=9.5)
    ax1.text(0.98, 0.06, f"Vmax = {vmax:g} µM/s\nKm = {km:g} mM", transform=ax1.transAxes,
             ha="right", va="bottom", color=INK_SECONDARY, fontsize=9.5)
    viz._finish(ax1, "[S] (mM)", "v₀ (µM/s)")
    _panel_title(ax1, "Velocidad inicial: v₀ = Vmax·[S]/(Km + [S])")

    frac = s / (km + s)
    ax2.plot(s, frac, color=COLORS["complejo_es"], linewidth=2.0, zorder=2)
    ax2.set_xlim(0, s_max)
    ax2.set_ylim(0, 1.05)
    viz._guide(ax2, "h", 0.5, start=0, end=km)
    viz._guide(ax2, "v", km, start=0, end=0.5)
    ax2.annotate("Km", xy=(km, 0), xytext=(4, 3), textcoords="offset points",
                 ha="left", va="bottom", color=INK_MUTED, fontsize=9.5)
    f_blood = _GLUCOSE_BLOOD_MM / (km + _GLUCOSE_BLOOD_MM)
    viz._markers(ax2, [_GLUCOSE_BLOOD_MM], [f_blood], COLORS["complejo_es"], size=8)
    ax2.annotate(f"[S] = 5 mM: {100 * f_blood:.0f} %", xy=(_GLUCOSE_BLOOD_MM, f_blood), xytext=(10, -4),
                 textcoords="offset points", ha="left", va="top", color=INK_SECONDARY, fontsize=9.5)
    ax2.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax2.set_yticklabels(["0 %", "25 %", "50 %", "75 %", "100 %"])
    viz._finish(ax2, "[S] (mM)", "fracción de enzima ocupada")
    _panel_title(ax2, "Ocupación: [ES]/[E]total = [S]/(Km + [S])")
    fig.tight_layout()

    resumen = (
        f"Con Km = {km:g} mM, a [S] = {km:g} mM la enzima trabaja al 50 % de su máximo "
        f"(v₀ = {vmax / 2:g} µM/s); a 5 mM (glucosa en sangre) trabaja al {100 * f_blood:.0f} % y "
        f"necesita [S] = 9·Km = {9 * km:g} mM para llegar al 90 %."
    )
    return fig, resumen


def explorar_michaelis_menten(vmax=8.0, km=5.0):
    """Explorador de Michaelis–Menten: deslizadores de ``Vmax`` (1–20 µM/s) y ``Km`` (0.5–30 mM).

    Panel izquierdo: ``v₀`` frente a ``[S]`` con guías de ``Vmax``, ``Vmax/2`` y
    ``Km``; panel derecho: fracción de enzima ocupada ``[S]/(Km + [S])``.
    """
    controls = {
        "vmax": _slider("Vmax (µM/s)", vmax, 1.0, 20.0, 0.5, ".1f"),
        "km": _slider("Km (mM)", km, 0.5, 30.0, 0.5, ".1f"),
    }
    return _explorer(_dibujar_michaelis_menten, controls)


# ---------------------------------------------------------------------------
# 2. Mecanismo E + S ⇌ ES → E + P
# ---------------------------------------------------------------------------
def _dibujar_mecanismo(k1, k_minus1, k2, e0, s0):
    ss = kin.steady_state_parameters(k1, k_minus1, k2)
    km, vmax = ss["km"], k2 * e0
    t_pre = 5.0 / (k1 * s0 + k_minus1 + k2)  # relajación del estado pre-estacionario
    t_end = max(1.2 * (s0 + km * np.log(50.0)) / vmax, 20.0 * t_pre)
    t_zoom = min(t_end, 10.0 * t_pre)
    sim = kin.simulate_mechanism(e0, s0, k1, k_minus1, k2, t_end, n_points=1200)
    zoom = kin.simulate_mechanism(e0, s0, k1, k_minus1, k2, t_zoom, n_points=600)

    fig, (ax1, ax2) = figure(11.0, 4.2, ncols=2)
    ax1.plot(sim["t"], sim["S"], color=COLORS["sustrato"], linewidth=2.0, label="S (sustrato)")
    ax1.plot(sim["t"], sim["P"], color=COLORS["producto_p"], linewidth=2.0, label="P (producto)")
    ax1.set_xlim(0, t_end)
    ax1.set_ylim(0, s0 * 1.08)
    viz._finish(ax1, "Tiempo (s)", "Concentración (µM)")
    _panel_title(ax1, "Sustrato y producto")
    viz._legend(ax1, loc="center right")

    ax2.plot(zoom["t"], zoom["E"], color=COLORS["enzima"], linewidth=2.0, label="E (enzima libre)")
    ax2.plot(zoom["t"], zoom["ES"], color=COLORS["complejo_es"], linewidth=2.0, label="ES (complejo)")
    ax2.set_xlim(0, t_zoom)
    ax2.set_ylim(0, e0 * 1.08)
    viz._finish(ax2, "Tiempo (s)", "Concentración (µM)")
    _panel_title(ax2, "Enzima libre y complejo ES (estado pre-estacionario)")
    viz._legend(ax2, loc="center right")
    fig.tight_layout()

    resumen = (
        f"Km = (k₋₁ + k₂)/k₁ = {km:.3g} µM; kcat = k₂ = {k2:.3g} s⁻¹; Vmax = kcat·[E]₀ = {vmax:.3g} µM/s. "
        f"Con [S]₀ = {s0:g} µM ({s0 / km:.2g}·Km) la v₀ simulada es {sim['v0']:.3g} µM/s "
        f"({100 * sim['v0'] / vmax:.0f} % de Vmax); el estado estacionario de ES se alcanza en ≈ {_tiempo(t_pre)}."
    )
    return fig, resumen


def explorar_mecanismo(k1=1.0, k_minus1=10.0, k2=5.0, e0=0.1, s0=50.0):
    """Explorador del mecanismo ``E + S ⇌ ES → E + P`` integrado numéricamente.

    Deslizadores logarítmicos de ``k₁`` (µM⁻¹·s⁻¹), ``k₋₁`` y ``k₂`` (s⁻¹) y
    lineales de ``[E]₀`` y ``[S]₀`` (µM). Muestra S y P (izquierda) y E y ES
    en el estado pre-estacionario (derecha) e imprime ``Km = (k₋₁ + k₂)/k₁`` y ``kcat``.
    """
    controls = {
        "k1": _slider("k₁ (µM⁻¹·s⁻¹)", k1, 0.01, 10.0, 0.1, ".3g", log=True),
        "k_minus1": _slider("k₋₁ (s⁻¹)", k_minus1, 0.1, 100.0, 0.1, ".3g", log=True),
        "k2": _slider("k₂ = kcat (s⁻¹)", k2, 0.1, 100.0, 0.1, ".3g", log=True),
        "e0": _slider("[E]₀ (µM)", e0, 0.01, 2.0, 0.01, ".2f"),
        "s0": _slider("[S]₀ (µM)", s0, 1.0, 200.0, 1.0, ".0f"),
    }
    return _explorer(_dibujar_mecanismo, controls)


# ---------------------------------------------------------------------------
# 3. Cooperatividad (Hill)
# ---------------------------------------------------------------------------
def _dibujar_hill(s_half, n, s_max=40.0):
    s = np.linspace(0.0, s_max, 400)
    v_mm = kin.michaelis_menten(s, 1.0, s_half)
    v_hill = kin.hill(s, 1.0, s_half, n)
    fig = viz.plot_hill_vs_mm(
        s, v_mm, v_hill, n_hill=n, s_half=s_half, ylabel="v₀ / Vmax",
        title="Cooperatividad: sigmoide de Hill frente a hipérbola",
    )
    ax = fig.axes[0]
    s10 = kin.substrate_at_fraction(1.0, s_half, n, 0.1)
    s90 = kin.substrate_at_fraction(1.0, s_half, n, 0.9)
    for frac, s_frac in ((0.1, s10), (0.9, s90)):
        viz._guide(ax, "h", frac, start=0, end=min(s_frac, s_max))
        if s_frac <= s_max:
            viz._guide(ax, "v", s_frac, start=0, end=frac)
            ax.annotate(f"{100 * frac:.0f} %: {s_frac:.3g} mM", xy=(s_frac, frac), xytext=(6, -4),
                        textcoords="offset points", ha="left", va="top", color=INK_MUTED, fontsize=9)
    ax.set_ylim(0, 1.08)
    resumen = (
        f"Con S₀.₅ = {s_half:g} mM y n = {n:.2g}, pasar del 10 % al 90 % de activación exige subir [S] "
        f"de {s10:.3g} a {s90:.3g} mM (×{s90 / s10:.1f}); con n = 1 (hipérbola) haría falta ×81."
    )
    return fig, resumen


def explorar_hill(s_half=7.5, n=1.7):
    """Explorador de cooperatividad: deslizadores de ``S₀.₅`` (1–20 mM) y ``n`` (1–4).

    Superpone la sigmoide de Hill y la hipérbola de Michaelis–Menten con el
    mismo ``S₀.₅`` e imprime el intervalo de ``[S]`` entre el 10 % y el 90 % de activación.
    """
    controls = {
        "s_half": _slider("S₀.₅ (mM)", s_half, 1.0, 20.0, 0.5, ".1f"),
        "n": _slider("Coeficiente de Hill n", n, 1.0, 4.0, 0.1, ".1f"),
    }
    return _explorer(_dibujar_hill, controls)


# ---------------------------------------------------------------------------
# 4. Inhibición
# ---------------------------------------------------------------------------
_INHIBITION_OPTIONS = [
    ("competitiva", "competitive"),
    ("acompetitiva", "uncompetitive"),
    ("no competitiva", "noncompetitive"),
    ("mixta", "mixed"),
]
_VMAX_INH, _KM_INH = 10.0, 5.0  # µM/s, mM (enzima de referencia del explorador)


def _dibujar_inhibicion(kind, i, ki, ki_prime, vmax=_VMAX_INH, km=_KM_INH):
    kind = str(kind).lower()
    kip = ki_prime if kind == "mixed" else None
    app = kin.apparent_parameters(kind, vmax, km, i, ki, kip)
    s_dense = np.linspace(0.0, 40.0, 300)
    s_points = np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0])
    curves_dense = [(0.0, kin.michaelis_menten(s_dense, vmax, km))]
    curves_points = [(0.0, kin.michaelis_menten(s_points, vmax, km))]
    if i > 0:
        curves_dense.append((i, kin.michaelis_menten(s_dense, app["vmax_app"], app["km_app"])))
        curves_points.append((i, kin.michaelis_menten(s_points, app["vmax_app"], app["km_app"])))

    fig, (ax1, ax2) = figure(12.0, 4.4, ncols=2)
    viz.plot_inhibition_family(s_dense, curves_dense, kind=kind, ax=ax1, title="")
    ax1.set_ylim(0, vmax * 1.12)
    viz._guide(ax1, "h", vmax)
    ax1.annotate("Vmax", xy=(1.0, vmax), xycoords=("axes fraction", "data"), xytext=(0, 3),
                 textcoords="offset points", ha="right", va="bottom", color=INK_MUTED, fontsize=9.5)
    _panel_title(ax1, f"{viz._INHIBITION_NAMES[kind]}: v₀ frente a [S]")
    viz.plot_inhibition_lineweaver(s_points, curves_points, kind=kind, ax=ax2, title="", subtitle="")
    _panel_title(ax2, "Lineweaver–Burk: 1/v₀ frente a 1/[S]")
    fig.tight_layout()

    ratio_v = app["vmax_app"] / vmax
    ratio_k = app["km_app"] / km
    nombre = next(name for name, key in _INHIBITION_OPTIONS if key == kind)
    ki_text = f"Ki = {ki:g} mM" + (f", Ki′ = {ki_prime:g} mM" if kind == "mixed" else "")
    resumen = (
        f"Inhibición {nombre} con [I] = {i:g} mM y {ki_text}: "
        f"Vmax_app = {app['vmax_app']:.3g} µM/s (×{ratio_v:.2f}), Km_app = {app['km_app']:.3g} mM (×{ratio_k:.2f}). "
        f"En Lineweaver–Burk {viz._INHIBITION_LB_PATTERN[kind]}."
    )
    return fig, resumen


def explorar_inhibicion(kind="competitive", i=2.0, ki=1.0, ki_prime=4.0):
    """Explorador de inhibición reversible (enzima de referencia: Vmax = 10 µM/s, Km = 5 mM).

    Desplegable del tipo (competitiva, acompetitiva, no competitiva, mixta) y
    deslizadores de ``[I]``, ``Ki`` y ``Ki′`` (solo activo en la mixta). Dos
    paneles: ``v₀`` frente a ``[S]`` con y sin inhibidor, y Lineweaver–Burk;
    imprime ``Vmax_app`` y ``Km_app``.
    """
    kind_widget = widgets.Dropdown(
        options=_INHIBITION_OPTIONS, value=kind, description="Tipo de inhibición",
        style={"description_width": "initial"}, layout=widgets.Layout(width="300px"),
    )
    controls = {
        "kind": kind_widget,
        "i": _slider("[I] (mM)", i, 0.0, 10.0, 0.25, ".2f"),
        "ki": _slider("Ki (mM)", ki, 0.1, 10.0, 0.1, ".1f"),
        "ki_prime": _slider("Ki′ (mM, solo mixta)", ki_prime, 0.1, 10.0, 0.1, ".1f"),
    }
    controls["ki_prime"].disabled = kind != "mixed"

    def _toggle(change):
        controls["ki_prime"].disabled = change["new"] != "mixed"

    kind_widget.observe(_toggle, "value")
    return _explorer(_dibujar_inhibicion, controls)


# ---------------------------------------------------------------------------
# 5. Eyring: barrera y temperatura
# ---------------------------------------------------------------------------
def _dibujar_eyring(delta_g, temperature):
    k = float(kin.eyring_rate(delta_g, temperature))
    decade = kin.R_KCAL * temperature * np.log(10.0)  # kcal/mol por factor 10 en k
    dg_grid = np.linspace(5.0, 30.0, 300)
    k_grid = kin.eyring_rate(dg_grid, temperature)

    fig, (ax1, ax2) = figure(11.0, 4.2, ncols=2)
    ax1.plot(dg_grid, k_grid, color=COLORS["datos"], linewidth=2.0, zorder=2)
    ax1.set_yscale("log")
    ax1.set_xlim(5.0, 30.0)
    viz._guide(ax1, "h", k, start=5.0, end=delta_g)
    viz._guide(ax1, "v", delta_g, start=k_grid.min(), end=k)
    viz._markers(ax1, [delta_g], [k], COLORS["ts"], size=8.5, zorder=5)
    ax1.annotate(f"k = {k:.3g} s⁻¹", xy=(delta_g, k), xytext=(8, 6), textcoords="offset points",
                 ha="left", va="bottom", color=INK, fontsize=10)
    viz._finish(ax1, "ΔG‡ (kcal/mol)", "k (s⁻¹, escala log)")
    _panel_title(ax1, f"k = (kB·T/h)·exp(−ΔG‡/RT) a {temperature:g} K")

    bars_dg = [delta_g - decade, delta_g, delta_g + decade]
    bars_k = [float(kin.eyring_rate(x, temperature)) for x in bars_dg]
    labels = [f"ΔG‡ − {decade:.2f}", f"ΔG‡ = {delta_g:g}", f"ΔG‡ + {decade:.2f}"]
    colors = [PALETTE[0], COLORS["ts"], PALETTE[0]]
    ax2.bar(range(3), bars_k, color=colors, width=0.6, zorder=2)
    ax2.set_yscale("log")
    ax2.set_xticks(range(3))
    ax2.set_xticklabels(labels)
    for x, kk in enumerate(bars_k):
        ax2.annotate(f"{kk:.3g} s⁻¹", xy=(x, kk), xytext=(0, 4), textcoords="offset points",
                     ha="center", va="bottom", color=INK_SECONDARY, fontsize=9.5)
    ax2.set_ylim(min(bars_k) / 30.0, max(bars_k) * 30.0)
    viz._finish(ax2, "barrera (kcal/mol)", "k (s⁻¹, escala log)")
    _panel_title(ax2, f"Cada {decade:.2f} kcal/mol cambian k en un factor 10")
    fig.tight_layout()

    resumen = (
        f"ΔG‡ = {delta_g:g} kcal/mol a {temperature:g} K → k = {k:.3g} s⁻¹; tiempo de recambio 1/k = {_tiempo(1 / k)}, "
        f"vida media t½ = ln2/k = {_tiempo(np.log(2) / k)}. A esta temperatura cada {decade:.2f} kcal/mol "
        f"adicionales de barrera dividen k entre 10 (1.36 kcal/mol a 298 K)."
    )
    return fig, resumen


def explorar_eyring(delta_g=15.0, temperature=298.15):
    """Explorador de la ecuación de Eyring: deslizadores de ``ΔG‡`` (5–30 kcal/mol) y ``T`` (273–323 K).

    Curva de ``k`` frente a ``ΔG‡`` en escala logarítmica y barras que
    ilustran la regla «cada ≈ 1.36 kcal/mol = ×10»; imprime ``1/k`` y ``t½``.
    """
    controls = {
        "delta_g": _slider("ΔG‡ (kcal/mol)", delta_g, 5.0, 30.0, 0.5, ".1f"),
        "temperature": _slider("T (K)", temperature, 273.0, 323.0, 1.0, ".0f"),
    }
    return _explorer(_dibujar_eyring, controls)


# ---------------------------------------------------------------------------
# 6. Dependencia con la temperatura (ΔH‡, ΔS‡)
# ---------------------------------------------------------------------------
def _dibujar_temperatura(delta_h, delta_s):
    temps = np.linspace(273.15, 333.15, 9)
    rates = kin.eyring_rate(delta_h - temps * delta_s / 1000.0, temps)
    delta_g = delta_h - _T_REF * delta_s / 1000.0
    fit = {"delta_h_kcal": delta_h, "delta_s_cal": delta_s, "delta_g_kcal": delta_g, "r2": 1.0}

    fig, (ax1, ax2) = figure(11.0, 4.2, ncols=2)
    viz.plot_eyring(temps, rates, fit=fit, ax=ax1, title="")
    _panel_title(ax1, "Gráfico de Eyring: ln(k/T) = ln(kB/h) + ΔS‡/R − (ΔH‡/R)·(1/T)")

    t_dense = np.linspace(273.15, 333.15, 300)
    k_dense = kin.eyring_rate(delta_h - t_dense * delta_s / 1000.0, t_dense)
    ax2.plot(t_dense, k_dense, color=COLORS["datos"], linewidth=2.0, zorder=2)
    k_ref = float(kin.eyring_rate(delta_g, _T_REF))
    k_body = float(kin.eyring_rate(delta_h - 310.15 * delta_s / 1000.0, 310.15))
    viz._markers(ax2, [_T_REF, 310.15], [k_ref, k_body], COLORS["datos"])
    ax2.annotate(f"25 °C: {k_ref:.3g} s⁻¹", xy=(_T_REF, k_ref), xytext=(-8, 4), textcoords="offset points",
                 ha="right", va="bottom", color=INK_SECONDARY, fontsize=9.5)
    ax2.annotate(f"37 °C: {k_body:.3g} s⁻¹", xy=(310.15, k_body), xytext=(8, -4), textcoords="offset points",
                 ha="left", va="top", color=INK_SECONDARY, fontsize=9.5)
    ax2.set_xlim(273.15, 333.15)
    ax2.set_ylim(bottom=0)
    viz._finish(ax2, "T (K)", "k (s⁻¹)")
    _panel_title(ax2, "k(T) = (kB·T/h)·exp(ΔS‡/R)·exp(−ΔH‡/RT)")
    fig.tight_layout()

    q10 = float(kin.eyring_rate(delta_h - 308.15 * delta_s / 1000.0, 308.15)) / k_ref
    resumen = (
        f"ΔH‡ = {delta_h:g} kcal/mol y ΔS‡ = {_num(delta_s, 0)} cal/(mol·K) → ΔG‡(298 K) = ΔH‡ − TΔS‡ = {delta_g:.2f} kcal/mol, "
        f"k(25 °C) = {k_ref:.3g} s⁻¹ y k(37 °C) = {k_body:.3g} s⁻¹ (×{k_body / k_ref:.2f}); "
        f"Q10 (25→35 °C) = {q10:.2f}."
    )
    return fig, resumen


def explorar_temperatura(delta_h=15.0, delta_s=-10.0):
    """Explorador de la dependencia con la temperatura: deslizadores de ``ΔH‡`` (kcal/mol) y ``ΔS‡`` (cal/mol/K).

    Gráfico de Eyring ``ln(k/T)`` frente a ``1/T`` (izquierda) y curva ``k(T)``
    (derecha) con marcas a 25 y 37 °C; imprime ``ΔG‡(298 K)`` y el ``Q10``.
    """
    controls = {
        "delta_h": _slider("ΔH‡ (kcal/mol)", delta_h, 5.0, 30.0, 0.5, ".1f"),
        "delta_s": _slider("ΔS‡ (cal/(mol·K))", delta_s, -30.0, 30.0, 1.0, ".0f"),
    }
    return _explorer(_dibujar_temperatura, controls)


# ---------------------------------------------------------------------------
# 7. Perfil de pH
# ---------------------------------------------------------------------------
def _dibujar_ph(pka1, pka2):
    ph = np.linspace(2.0, 12.0, 400)
    v = kin.bell_shaped_ph_profile(ph, 1.0, pka1, pka2)
    ph_opt = 0.5 * (pka1 + pka2)
    v_opt = float(kin.bell_shaped_ph_profile(ph_opt, 1.0, pka1, pka2))
    fig = viz.plot_ph_profile(
        ph, v, pkas=(pka1, pka2), ylabel="v₀ / Vmax", title="Perfil de pH en campana",
        subtitle="v = Vmax / (1 + 10^(pKa₁ − pH) + 10^(pH − pKa₂))",
    )
    ax = fig.axes[0]
    if pka2 - pka1 < 2.5:  # pKa próximos: la etiqueta de pKa₁ pasa a la izquierda de su guía
        for text in ax.texts:
            if text.get_text().startswith("pKa₁"):
                text.set_ha("right")
                text.xyann = (-4, -2)
    viz._markers(ax, [ph_opt], [v_opt], COLORS["ts"], size=8.5, zorder=5)
    ax.annotate(f"óptimo: pH {ph_opt:.2f} ({100 * v_opt:.0f} % de Vmax)", xy=(ph_opt, v_opt), xytext=(0, 8),
                textcoords="offset points", ha="center", va="bottom", color=INK, fontsize=9.5)
    ax.set_xlim(2.0, 12.0)
    ax.set_ylim(0, 1.18)
    resumen = (
        f"Con pKa₁ = {pka1:g} y pKa₂ = {pka2:g} el óptimo está en pH = (pKa₁ + pKa₂)/2 = {ph_opt:.2f}, "
        f"donde la enzima alcanza el {100 * v_opt:.0f} % de Vmax"
    )
    if pka2 - pka1 < 2.0:
        resumen += " (los pKa están tan próximos que nunca se llega al máximo teórico)."
    else:
        resumen += f"; la actividad cae a la mitad en pH ≈ {pka1:g} y ≈ {pka2:g}."
    return fig, resumen


def explorar_ph(pka1=6.0, pka2=9.0):
    """Explorador del perfil de pH en campana: deslizadores de ``pKa₁`` (3–8) y ``pKa₂`` (6–11).

    Dibuja ``v = Vmax/(1 + 10^(pKa₁ − pH) + 10^(pH − pKa₂))`` e imprime el pH óptimo.
    """
    controls = {
        "pka1": _slider("pKa₁ (grupo que debe estar desprotonado)", pka1, 3.0, 8.0, 0.1, ".1f"),
        "pka2": _slider("pKa₂ (grupo que debe estar protonado)", pka2, 6.0, 11.0, 0.1, ".1f"),
    }
    return _explorer(_dibujar_ph, controls)


# ---------------------------------------------------------------------------
# 8. Activador alostérico de la glucoquinasa (GKA)
# ---------------------------------------------------------------------------
def _dibujar_activador(vmax_factor, s_half_factor, s_max=30.0):
    s = np.linspace(0.0, s_max, 400)
    basal = kin.hill(s, 100.0, _GK_S_HALF, _GK_N)
    activated = kin.activator_effect(s, 100.0, _GK_S_HALF, _GK_N, vmax_factor, s_half_factor)
    a_blood = float(kin.hill(_GLUCOSE_BLOOD_MM, 100.0, _GK_S_HALF, _GK_N))
    b_blood = float(kin.activator_effect(_GLUCOSE_BLOOD_MM, 100.0, _GK_S_HALF, _GK_N, vmax_factor, s_half_factor))
    s_half_act = _GK_S_HALF * s_half_factor

    fig, ax = figure(8.0, 4.6)
    ax.plot(s, basal, color=COLORS["michaelis_menten"], linewidth=2.0, label=f"sin activador (S₀.₅ = {_GK_S_HALF} mM, n = {_GK_N})")
    ax.plot(s, activated, color=COLORS["hill"], linewidth=2.0, label=f"con activador (Vmax ×{vmax_factor:g}, S₀.₅ = {s_half_act:.2g} mM)")
    viz._guide(ax, "v", _GLUCOSE_BLOOD_MM)
    ax.annotate("glucosa en sangre (5 mM)", xy=(_GLUCOSE_BLOOD_MM, 1.0), xycoords=("data", "axes fraction"),
                xytext=(4, -2), textcoords="offset points", ha="left", va="top", color=INK_MUTED, fontsize=9.5)
    viz._markers(ax, [_GLUCOSE_BLOOD_MM], [a_blood], COLORS["michaelis_menten"], size=8)
    viz._markers(ax, [_GLUCOSE_BLOOD_MM], [b_blood], COLORS["hill"], size=8)
    ax.annotate(f"{a_blood:.0f} %", xy=(_GLUCOSE_BLOOD_MM, a_blood), xytext=(8, -4), textcoords="offset points",
                ha="left", va="top", color=INK, fontsize=10)
    ax.annotate(f"{b_blood:.0f} %", xy=(_GLUCOSE_BLOOD_MM, b_blood), xytext=(-8, 4), textcoords="offset points",
                ha="right", va="bottom", color=INK, fontsize=10)
    ax.set_xlim(0, s_max)
    ax.set_ylim(0, max(100.0 * vmax_factor, 100.0) * 1.12)
    viz._finish(ax, "[glucosa] (mM)", "actividad (% de la Vmax basal)",
                "Activador alostérico de la glucoquinasa (GKA)")
    viz._legend(ax, loc="lower right")

    resumen = (
        f"A 5 mM de glucosa: sin activador la glucoquinasa trabaja al {a_blood:.0f} % de su Vmax basal; "
        f"con activador (Vmax ×{vmax_factor:g}, S₀.₅ ×{s_half_factor:g} → {s_half_act:.2g} mM) trabaja al "
        f"{b_blood:.0f} % (×{b_blood / a_blood:.1f})."
    )
    return fig, resumen


def explorar_activador(vmax_factor=1.5, s_half_factor=0.4):
    """Explorador del activador de la glucoquinasa: deslizadores de ``factor Vmax`` (1–3) y ``factor S₀.₅`` (0.1–1).

    Curva de Hill de la glucoquinasa (S₀.₅ = 7.5 mM, n = 1.7) con y sin
    activador, con la glucemia (5 mM) marcada; imprime la actividad a 5 mM en ambos casos.
    """
    controls = {
        "vmax_factor": _slider("factor sobre Vmax", vmax_factor, 1.0, 3.0, 0.1, ".1f"),
        "s_half_factor": _slider("factor sobre S₀.₅", s_half_factor, 0.1, 1.0, 0.05, ".2f"),
    }
    return _explorer(_dibujar_activador, controls)


# ---------------------------------------------------------------------------
# 9. Perfil de energía de reacción (ilustrativo)
# ---------------------------------------------------------------------------
def _dibujar_perfil_energia(barrier, reaction_energy):
    from scipy.interpolate import CubicHermiteSpline

    barrier_eff = max(float(barrier), float(reaction_energy) + 1.0)
    nodes_x = np.array([0.0, 0.5, 1.0])
    nodes_e = np.array([0.0, barrier_eff, float(reaction_energy)])
    spline = CubicHermiteSpline(nodes_x, nodes_e, np.zeros(3))
    x = np.linspace(0.0, 1.0, 300)
    e = spline(x)

    fig, ax = figure(8.0, 4.6)
    ax.plot(x, e, color=COLORS["reactivo"], linewidth=2.0, zorder=2)
    for xi, ei, color, text in ((0.0, 0.0, COLORS["reactivo"], "R"), (0.5, barrier_eff, COLORS["ts"], "TS‡"),
                                (1.0, reaction_energy, COLORS["producto"], "P")):
        viz._markers(ax, [xi], [ei], color, size=8.5, zorder=5)
        ax.annotate(text, xy=(xi, ei), xytext=(0, 9), textcoords="offset points", ha="center", va="bottom",
                    color=INK_SECONDARY, fontsize=10)
    viz._guide(ax, "h", 0.0, start=0.0, end=0.62)
    viz._guide(ax, "h", barrier_eff, start=0.5, end=0.62, linestyle="-", linewidth=0.8, color=AXIS)
    viz._dimension_arrow(ax, 0.62, 0.0, barrier_eff, f"ΔE‡ = {_num(barrier_eff)} kcal/mol")
    viz._guide(ax, "h", 0.0, start=0.8, end=1.12)
    viz._guide(ax, "h", reaction_energy, start=1.0, end=1.12, linestyle="-", linewidth=0.8, color=AXIS)
    if abs(reaction_energy) > 0.05:
        viz._dimension_arrow(ax, 1.12, 0.0, reaction_energy, f"ΔE = {_num(reaction_energy)} kcal/mol")
    ax.set_xlim(-0.08, 1.55)
    lo = min(0.0, reaction_energy) - 4.0
    ax.set_ylim(lo, barrier_eff + 5.0)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_xticklabels(["reactivos", "estado de transición", "productos"])
    viz._finish(ax, "Coordenada de reacción", "Energía relativa (kcal/mol)", "Perfil de energía de la reacción",
                grid_axis="y")

    k = float(kin.eyring_rate(barrier_eff, _T_REF))
    k_rev = float(kin.eyring_rate(barrier_eff - reaction_energy, _T_REF))
    keq = float(np.exp(-reaction_energy / (kin.R_KCAL * _T_REF)))
    kind = "exergónica (ΔE < 0)" if reaction_energy < 0 else ("endergónica (ΔE > 0)" if reaction_energy > 0 else "termoneutra")
    resumen = (
        f"ΔE‡ = {barrier_eff:g} kcal/mol → k(298 K, Eyring, tratando ΔE‡ como ΔG‡) = {k:.3g} s⁻¹ "
        f"(1/k = {_tiempo(1 / k)}); reacción {kind}: barrera inversa {barrier_eff - reaction_energy:g} kcal/mol, "
        f"k inversa = {k_rev:.3g} s⁻¹, K_eq = exp(−ΔE/RT) = {keq:.3g}."
    )
    if barrier_eff != barrier:
        resumen += f" (La barrera se elevó a {barrier_eff:g} kcal/mol para que el TS quede por encima de los productos.)"
    return fig, resumen


def explorar_perfil_energia(barrier=15.0, reaction_energy=-5.0):
    """Explorador ilustrativo del perfil de energía: deslizadores de ``ΔE‡`` (5–30) y ``ΔE`` (−20…+10 kcal/mol).

    Dibuja un perfil suave (spline cúbica con pendiente nula en R, TS y P) e
    imprime ``k`` por Eyring y la constante de equilibrio ``exp(−ΔE/RT)``.
    """
    controls = {
        "barrier": _slider("ΔE‡ barrera (kcal/mol)", barrier, 5.0, 30.0, 0.5, ".1f"),
        "reaction_energy": _slider("ΔE de reacción (kcal/mol)", reaction_energy, -20.0, 10.0, 0.5, ".1f"),
    }
    return _explorer(_dibujar_perfil_energia, controls)


# ---------------------------------------------------------------------------
# Todos los exploradores en pestañas
# ---------------------------------------------------------------------------
_EXPLORADORES = [
    ("Michaelis–Menten", explorar_michaelis_menten),
    ("Mecanismo (EDO)", explorar_mecanismo),
    ("Hill", explorar_hill),
    ("Inhibición", explorar_inhibicion),
    ("Eyring", explorar_eyring),
    ("Temperatura", explorar_temperatura),
    ("pH", explorar_ph),
    ("Activador (GKA)", explorar_activador),
    ("Perfil de energía", explorar_perfil_energia),
]


def todo():
    """Pestañas (``ipywidgets.Tab``) con los nueve exploradores, construidos al llamar."""
    tab = widgets.Tab(children=[build() for _, build in _EXPLORADORES])
    for index, (name, _) in enumerate(_EXPLORADORES):
        tab.set_title(index, name)
    return tab
