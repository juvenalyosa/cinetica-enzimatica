"""Exploradores interactivos (``ipywidgets``) para el cuaderno *Cinética Enzimática*.

Cada función ``explorar_*`` construye deslizadores en español, dibuja la figura
con el estilo de :mod:`enzimas.viz` y devuelve un ``ipywidgets.VBox`` (la celda
lo muestra si es la última expresión; también puede pasarse a ``display``).
Debajo de la figura aparece un recuadro («Qué dicen estos números») que
traduce los valores actuales a palabras; el título de la figura es un mensaje
que también cambia con los deslizadores.

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
_FONT = "Figtree, 'Avenir Next', 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"


def _slider(description, value, lo, hi, step, fmt=".2f", log=False):
    """Deslizador con etiqueta completa en español, ancho cómodo y redibujado al soltar."""
    common = dict(
        value=value, description=description, continuous_update=False, readout_format=fmt,
        style={"description_width": "210px", "handle_color": viz.PALETTE[0]},
        layout=widgets.Layout(width="460px", margin="2px 18px 2px 0"),
    )
    if log:
        return widgets.FloatLogSlider(base=10, min=np.log10(lo), max=np.log10(hi), step=step, **common)
    return widgets.FloatSlider(min=lo, max=hi, step=step, **common)


def _mostrar(obj):
    """``display(obj)`` solo si hay un intérprete de IPython (evita imprimir el ``repr`` en consola)."""
    if get_ipython() is not None:
        display(obj)


def _cabecera(titulo, instruccion):
    return widgets.HTML(
        f'<div style="font-family:{_FONT};margin:0 0 6px 2px">'
        f'<div style="font-size:16px;font-weight:650;color:{INK}">🎛️ {titulo}</div>'
        f'<div style="font-size:13px;color:{INK_SECONDARY};margin-top:2px">{instruccion}</div></div>')


def _explorer(update_fn, controls, titulo=None, instruccion="Mueve los deslizadores: la gráfica y la frase de abajo se actualizan."):
    """Conecta ``controls`` (``dict`` nombre → widget) con ``update_fn(**valores) -> (fig, resumen)``.

    Usa ``ipywidgets.interactive_output`` (que aporta el ``Output`` y el
    ``clear_output(wait=True)``): en cada cambio de valor se crea una figura
    nueva, se muestra, debajo aparece el resumen en un recuadro con estilo y se
    cierra la figura. Devuelve un ``VBox`` con atributos ``controles`` (el ``dict``),
    ``salida`` (``Output``), ``estado`` (último resumen y último error) y ``redibujar()``.
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
                _mostrar(viz.mensaje(resumen, tipo="idea", titulo="Qué dicen estos números"))
        except Exception as exc:
            estado["error"] = exc
            raise
        finally:
            if fig is not None:
                plt.close(fig)

    salida = widgets.interactive_output(_redibujar, controls)
    panel = widgets.HBox(list(controls.values()), layout=widgets.Layout(flex_flow="row wrap", margin="0 0 4px 0"))
    hijos = [panel, salida]
    if titulo:
        hijos.insert(0, _cabecera(titulo, instruccion))
    caja = widgets.VBox(hijos, layout=widgets.Layout(border=f"1px solid {viz.GRID}", padding="12px 14px", margin="4px 0"))
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
    ax.set_title(text, loc="left", fontsize=12.5, color=INK_SECONDARY, fontweight="semibold", pad=10)


def _titulo(fig, titulo, subtitulo=None):
    """Título-mensaje de la figura (cambia con los deslizadores) y subtítulo de lectura."""
    viz._fig_title(fig, titulo, subtitulo)


# ---------------------------------------------------------------------------
# 1. Michaelis–Menten
# ---------------------------------------------------------------------------
def _dibujar_michaelis_menten(vmax, km, s_max=60.0):
    s = np.linspace(0.0, s_max, 400)
    v = kin.michaelis_menten(s, vmax, km)
    color = COLORS["michaelis_menten"]
    fig, (ax1, ax2) = figure(13.0, 5.4, ncols=2, gridspec_kw={"width_ratios": [1.55, 1.0]})

    ax1.fill_between(s, 0, v, color=viz._tint(color, 0.9), zorder=1, linewidth=0)
    ax1.plot(s, v, color=color, linewidth=2.8, zorder=3)
    ax1.set_xlim(0, s_max)
    ax1.set_ylim(0, 20.0 * 1.12)  # escala fija: al mover Vmax se ve crecer la curva
    if 4 * km < 0.85 * s_max:
        viz._kband(ax1, 4 * km, s_max, "saturación", y_text=0.97)
    viz._guide(ax1, "h", vmax, color=INK_SECONDARY)
    ax1.annotate(f"Vmax = {vmax:g}", xy=(1.0, vmax), xycoords=("axes fraction", "data"), xytext=(-4, 4),
                 textcoords="offset points", ha="right", va="bottom", color=INK, fontsize=11.5, fontweight="semibold")
    viz._guide(ax1, "h", vmax / 2, start=0, end=km)
    viz._guide(ax1, "v", km, start=0, end=vmax / 2)
    viz._keypoint(ax1, km, vmax / 2, COLORS["ts"])
    ax1.annotate(f"mitad del máximo\nen [S] = Km = {km:g} mM", xy=(km, vmax / 2), xytext=(12, -6),
                 textcoords="offset points", ha="left", va="top", color=INK, fontsize=11)
    viz._finish(ax1, "[S] (mM)", "v₀ (µM/s)")
    _panel_title(ax1, "Velocidad inicial frente a sustrato")

    marcas = [("Km/10", km / 10), ("Km", km), ("10·Km", 10 * km)]
    if abs(km - _GLUCOSE_BLOOD_MM) > 0.4:
        marcas.append(("glucosa en sangre (5 mM)", _GLUCOSE_BLOOD_MM))
    marcas = sorted(marcas, key=lambda m: m[1])
    y = np.arange(len(marcas))[::-1]
    fracs = [m[1] / (km + m[1]) for m in marcas]
    ax2.barh(y, [100] * len(marcas), color=viz.GRID, height=0.52, zorder=1)
    ax2.barh(y, [100 * f for f in fracs], color=COLORS["complejo_es"], height=0.52, zorder=2)
    for yy, (name, sv), f in zip(y, marcas, fracs):
        ax2.text(0, yy + 0.43, f"[S] = {name}" + ("" if "mM" in name else f"  ({sv:.3g} mM)"), ha="left", va="bottom",
                 color=INK_SECONDARY, fontsize=10.5)
        ax2.text(min(100 * f + 2, 86), yy, f"{100 * f:.0f} %", ha="left", va="center", color=INK, fontsize=11.5,
                 fontweight="semibold", zorder=3)
    ax2.set_xlim(0, 100)
    ax2.set_ylim(-0.6, len(marcas) - 0.1)
    ax2.set_yticks([])
    viz._finish(ax2, "enzima ocupada (%)", None, grid_axis=None)
    ax2.spines["left"].set_visible(False)
    _panel_title(ax2, "¿Qué fracción de la enzima trabaja?")
    _titulo(fig, f"Con Km = {km:g} mM, la enzima va a media máquina cuando [S] = {km:g} mM",
            "Izquierda: la hipérbola. Derecha: la misma idea como ocupación de la enzima, [S]/(Km + [S]).")

    f_blood = _GLUCOSE_BLOOD_MM / (km + _GLUCOSE_BLOOD_MM)
    resumen = (
        f"Con Km = {km:g} mM, a [S] = {km:g} mM la enzima trabaja al 50 % de su máximo "
        f"(v₀ = {vmax / 2:g} µM/s); a 5 mM (glucosa en sangre) trabaja al {100 * f_blood:.0f} % y "
        f"necesita [S] = 9·Km = {9 * km:g} mM para llegar al 90 %."
    )
    return fig, resumen


def explorar_michaelis_menten(vmax=8.0, km=5.0):
    """Explorador de Michaelis–Menten: deslizadores de ``Vmax`` (1–20 µM/s) y ``Km`` (0.5–30 mM).

    Panel izquierdo: ``v₀`` frente a ``[S]`` con ``Vmax`` y el punto ``(Km, Vmax/2)``;
    panel derecho: fracción de enzima ocupada ``[S]/(Km + [S])`` a varias ``[S]``.
    """
    controls = {
        "vmax": _slider("Vmax: velocidad máxima (µM/s)", vmax, 1.0, 20.0, 0.5, ".1f"),
        "km": _slider("Km: [S] a media velocidad (mM)", km, 0.5, 30.0, 0.5, ".1f"),
    }
    return _explorer(_dibujar_michaelis_menten, controls, "La hipérbola de Michaelis–Menten")


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

    fig, (ax1, ax2) = figure(13.0, 5.2, ncols=2)
    ax1.fill_between(sim["t"], 0, sim["P"], color=viz._tint(COLORS["producto_p"], 0.88), zorder=1, linewidth=0)
    ax1.plot(sim["t"], sim["S"], color=COLORS["sustrato"], linewidth=2.8, zorder=3)
    ax1.plot(sim["t"], sim["P"], color=COLORS["producto_p"], linewidth=2.8, zorder=3)
    k_lab = int(0.12 * len(sim["t"]))
    ax1.annotate("S (sustrato) se gasta", xy=(sim["t"][k_lab], sim["S"][k_lab]), xytext=(10, 6), textcoords="offset points",
                 ha="left", va="bottom", color=INK, fontsize=11)
    ax1.annotate("P (producto) se acumula", xy=(sim["t"][k_lab], sim["P"][k_lab]), xytext=(10, -6), textcoords="offset points",
                 ha="left", va="top", color=INK, fontsize=11)
    ax1.set_xlim(0, t_end)
    ax1.set_ylim(0, s0 * 1.1)
    viz._finish(ax1, "Tiempo (s)", "Concentración (µM)")
    _panel_title(ax1, "La reacción completa")

    ax2.axvspan(0, min(t_pre, t_zoom), color=viz._tint(COLORS["ts"], 0.88), zorder=0, linewidth=0)
    ax2.annotate("arranque", xy=(0, 0.97), xycoords=ax2.get_xaxis_transform(), xytext=(4, 0), textcoords="offset points",
                 ha="left", va="top", color=INK_MUTED, fontsize=10.5)
    ax2.text(0.5 * (min(t_pre, t_zoom) + t_zoom), 0.97, "estado estacionario: ES casi constante",
             transform=ax2.get_xaxis_transform(), ha="center", va="top", color=INK_MUTED, fontsize=10.5)
    ax2.plot(zoom["t"], zoom["E"], color=COLORS["enzima"], linewidth=2.8, label="E (enzima libre)", zorder=3)
    ax2.plot(zoom["t"], zoom["ES"], color=COLORS["complejo_es"], linewidth=2.8, label="ES (complejo)", zorder=3)
    ax2.set_xlim(0, t_zoom)
    ax2.set_ylim(0, e0 * 1.18)
    viz._finish(ax2, "Tiempo (s)", "Concentración (µM)")
    _panel_title(ax2, "Zoom a los primeros instantes: la enzima")
    viz._legend(ax2, loc="center right")
    _titulo(fig, f"ES se estabiliza en {_tiempo(t_pre)} y la enzima produce a v₀ = {sim['v0']:.3g} µM/s",
            "La hipótesis del estado estacionario (ES constante) es la base de Michaelis–Menten.")

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
    en el estado pre-estacionario (derecha) y resume ``Km = (k₋₁ + k₂)/k₁`` y ``kcat``.
    """
    controls = {
        "k1": _slider("k₁: el sustrato entra (µM⁻¹·s⁻¹)", k1, 0.01, 10.0, 0.1, ".3g", log=True),
        "k_minus1": _slider("k₋₁: el sustrato se suelta (s⁻¹)", k_minus1, 0.1, 100.0, 0.1, ".3g", log=True),
        "k2": _slider("k₂ = kcat: reacciona (s⁻¹)", k2, 0.1, 100.0, 0.1, ".3g", log=True),
        "e0": _slider("[E]₀: enzima en el tubo (µM)", e0, 0.01, 2.0, 0.01, ".2f"),
        "s0": _slider("[S]₀: sustrato inicial (µM)", s0, 1.0, 200.0, 1.0, ".0f"),
    }
    return _explorer(_dibujar_mecanismo, controls, "El mecanismo E + S ⇌ ES → E + P, paso a paso")


# ---------------------------------------------------------------------------
# 3. Cooperatividad (Hill)
# ---------------------------------------------------------------------------
def _dibujar_hill(s_half, n, s_max=40.0):
    s = np.linspace(0.0, s_max, 400)
    v_mm = kin.michaelis_menten(s, 1.0, s_half)
    v_hill = kin.hill(s, 1.0, s_half, n)
    s10 = kin.substrate_at_fraction(1.0, s_half, n, 0.1)
    s90 = kin.substrate_at_fraction(1.0, s_half, n, 0.9)
    fig = viz.plot_hill_vs_mm(
        s, v_mm, v_hill, n_hill=n, s_half=s_half, ylabel="v₀ / Vmax",
        title=f"Con n = {n:.2g}, basta multiplicar [S] por {s90 / s10:.1f} para pasar del 10 % al 90 %",
        subtitle="Franja naranja: la ventana 10 %–90 % de la sigmoide. Con n = 1 (hipérbola) haría falta ×81.",
    )
    fig.set_size_inches(10.5, 5.6)
    ax = fig.axes[0]
    if s10 < s_max:
        viz._kband(ax, s10, min(s90, s_max), color=viz._tint(COLORS["hill"], 0.9))
    for frac, s_frac in ((0.1, s10), (0.9, s90)):
        viz._guide(ax, "h", frac, start=0, end=min(s_frac, s_max))
        if s_frac <= s_max:
            viz._guide(ax, "v", s_frac, start=0, end=frac)
            viz._keypoint(ax, s_frac, frac, COLORS["hill"], size=7)
            ax.annotate(f"{100 * frac:.0f} %: {s_frac:.3g} mM", xy=(s_frac, frac), xytext=(8, -6),
                        textcoords="offset points", ha="left", va="top", color=INK, fontsize=10.5)
    ax.set_ylim(0, 1.08)
    resumen = (
        f"Con S₀.₅ = {s_half:g} mM y n = {n:.2g}, pasar del 10 % al 90 % de activación exige subir [S] "
        f"de {s10:.3g} a {s90:.3g} mM (×{s90 / s10:.1f}); con n = 1 (hipérbola) haría falta ×81."
    )
    return fig, resumen


def explorar_hill(s_half=7.5, n=1.7):
    """Explorador de cooperatividad: deslizadores de ``S₀.₅`` (1–20 mM) y ``n`` (1–4).

    Superpone la sigmoide de Hill y la hipérbola de Michaelis–Menten con el
    mismo ``S₀.₅`` y resume el intervalo de ``[S]`` entre el 10 % y el 90 % de activación.
    """
    controls = {
        "s_half": _slider("S₀.₅: [S] a media actividad (mM)", s_half, 1.0, 20.0, 0.5, ".1f"),
        "n": _slider("n: coeficiente de Hill", n, 1.0, 4.0, 0.1, ".1f"),
    }
    return _explorer(_dibujar_hill, controls, "Interruptor (sigmoide) frente a regulador (hipérbola)")


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

    fig, (ax1, ax2) = figure(13.5, 5.4, ncols=2)
    viz.plot_inhibition_family(s_dense, curves_dense, kind=kind, ax=ax1, title="")
    ax1.set_ylim(0, vmax * 1.15)
    viz._guide(ax1, "h", vmax, color=INK_SECONDARY)
    ax1.annotate("Vmax sin inhibidor", xy=(1.0, vmax), xycoords=("axes fraction", "data"), xytext=(-4, 4),
                 textcoords="offset points", ha="right", va="bottom", color=INK, fontsize=10.5)
    viz._keypoint(ax1, km, vmax / 2, viz.sequential_blue(2)[0], size=8)
    if i > 0:
        viz._keypoint(ax1, app["km_app"], app["vmax_app"] / 2, viz.sequential_blue(2)[1], size=8)
        if abs(app["km_app"] - km) > 0.3:
            ax1.annotate("", xy=(app["km_app"], app["vmax_app"] / 2), xytext=(km, vmax / 2),
                         arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=1.3, shrinkA=8, shrinkB=8), zorder=8)
        ax1.annotate("Km aparente", xy=(app["km_app"], app["vmax_app"] / 2), xytext=(10, -8), textcoords="offset points",
                     ha="left", va="top", color=INK, fontsize=10.5, fontweight="semibold")
    _panel_title(ax1, "v₀ frente a [S]: puntos = (Km, Vmax/2)")
    viz.plot_inhibition_lineweaver(s_points, curves_points, kind=kind, ax=ax2, title="", subtitle="")
    _panel_title(ax2, "Su huella en Lineweaver–Burk")

    ratio_v = app["vmax_app"] / vmax
    ratio_k = app["km_app"] / km
    nombre = next(name for name, key in _INHIBITION_OPTIONS if key == kind)
    _titulo(fig, f"Inhibición {nombre}: Vmax ×{ratio_v:.2f}, Km aparente ×{ratio_k:.2f}",
            f"En Lineweaver–Burk {viz._INHIBITION_LB_PATTERN[kind]}.")
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
    paneles: ``v₀`` frente a ``[S]`` con y sin inhibidor (con el desplazamiento de
    ``(Km, Vmax/2)``), y Lineweaver–Burk; resume ``Vmax_app`` y ``Km_app``.
    """
    kind_widget = widgets.Dropdown(
        options=_INHIBITION_OPTIONS, value=kind, description="Tipo de inhibidor",
        style={"description_width": "210px"}, layout=widgets.Layout(width="460px", margin="2px 18px 2px 0"),
    )
    controls = {
        "kind": kind_widget,
        "i": _slider("[I]: cantidad de inhibidor (mM)", i, 0.0, 10.0, 0.25, ".2f"),
        "ki": _slider("Ki: potencia (menor = más potente, mM)", ki, 0.1, 10.0, 0.1, ".1f"),
        "ki_prime": _slider("Ki′ (mM, solo en la mixta)", ki_prime, 0.1, 10.0, 0.1, ".1f"),
    }
    controls["ki_prime"].disabled = kind != "mixed"

    def _toggle(change):
        controls["ki_prime"].disabled = change["new"] != "mixed"

    kind_widget.observe(_toggle, "value")
    return _explorer(_dibujar_inhibicion, controls, "Cuatro formas de frenar una enzima")


# ---------------------------------------------------------------------------
# 5. Eyring: barrera y temperatura
# ---------------------------------------------------------------------------
_REFERENCIAS_K = [(60.0, "k_cat de la glucoquinasa (~60 s⁻¹)"), (1 / 3600.0, "una vez por hora"),
                  (1 / 3.15576e7, "una vez por año")]


def _dibujar_eyring(delta_g, temperature):
    k = float(kin.eyring_rate(delta_g, temperature))
    decade = kin.R_KCAL * temperature * np.log(10.0)  # kcal/mol por factor 10 en k
    dg_grid = np.linspace(5.0, 30.0, 300)
    k_grid = kin.eyring_rate(dg_grid, temperature)

    fig, (ax1, ax2) = figure(13.5, 5.4, ncols=2, gridspec_kw={"width_ratios": [1.6, 1.0]})
    ax1.plot(dg_grid, k_grid, color=COLORS["datos"], linewidth=2.8, zorder=3)
    ax1.fill_between(dg_grid, k_grid.min() / 10, k_grid, color=viz._tint(COLORS["datos"], 0.9), zorder=1, linewidth=0)
    ax1.set_yscale("log")
    ax1.set_xlim(5.0, 30.0)
    ax1.set_ylim(k_grid.min() / 10, k_grid.max() * 10)
    for k_ref, texto in _REFERENCIAS_K:
        ax1.axhline(k_ref, color=AXIS, linewidth=1.0, linestyle=(0, (2, 3)), zorder=2)
        ax1.annotate(texto, xy=(5.0, k_ref), xytext=(6, 3), textcoords="offset points", ha="left", va="bottom",
                     color=INK_MUTED, fontsize=10)
    viz._guide(ax1, "h", k, start=5.0, end=delta_g)
    viz._guide(ax1, "v", delta_g, start=k_grid.min() / 10, end=k)
    viz._keypoint(ax1, delta_g, k, COLORS["ts"], size=11)
    ax1.annotate(f"k = {k:.3g} s⁻¹", xy=(delta_g, k), xytext=(12, 8), textcoords="offset points",
                 ha="left", va="bottom", color=INK, fontsize=12, fontweight="semibold")
    viz._finish(ax1, "ΔG‡: altura de la colina (kcal/mol)", "k (s⁻¹, escala log)")
    _panel_title(ax1, f"k = (kB·T/h)·exp(−ΔG‡/RT) a {temperature:g} K ({temperature - 273.15:.0f} °C)")

    bars_dg = [delta_g - decade, delta_g, delta_g + decade]
    bars_k = [float(kin.eyring_rate(x, temperature)) for x in bars_dg]
    labels = [f"{_num(bars_dg[0], 2)}\n(−{decade:.2f})", f"{_num(delta_g, 2)}\n(la tuya)", f"{_num(bars_dg[2], 2)}\n(+{decade:.2f})"]
    colors = [viz._tint(COLORS["datos"], 0.35), COLORS["ts"], viz._tint(COLORS["datos"], 0.35)]
    ax2.bar(range(3), bars_k, color=colors, width=0.62, zorder=2)
    ax2.set_yscale("log")
    ax2.set_xticks(range(3))
    ax2.set_xticklabels(labels)
    for x, kk, tag in zip(range(3), bars_k, ("×10", "", "÷10")):
        ax2.annotate(f"{kk:.3g} s⁻¹" + (f"\n{tag}" if tag else ""), xy=(x, kk), xytext=(0, 5), textcoords="offset points",
                     ha="center", va="bottom", color=INK, fontsize=10.5, fontweight="semibold" if not tag else "normal")
    ax2.set_ylim(min(bars_k) / 30.0, max(bars_k) * 60.0)
    viz._finish(ax2, "barrera ΔG‡ (kcal/mol)", "k (s⁻¹, escala log)")
    _panel_title(ax2, f"±{decade:.2f} kcal/mol = ×10 en k")
    _titulo(fig, f"Con ΔG‡ = {delta_g:g} kcal/mol, cada molécula reacciona en {_tiempo(1 / k)}",
            "Mueve ΔG‡: fíjate en que la escala vertical es logarítmica (cada línea es un factor 10).")

    resumen = (
        f"ΔG‡ = {delta_g:g} kcal/mol a {temperature:g} K → k = {k:.3g} s⁻¹; tiempo de recambio 1/k = {_tiempo(1 / k)}, "
        f"vida media t½ = ln2/k = {_tiempo(np.log(2) / k)}. A esta temperatura cada {decade:.2f} kcal/mol "
        f"adicionales de barrera dividen k entre 10 (1.36 kcal/mol a 298 K)."
    )
    return fig, resumen


def explorar_eyring(delta_g=15.0, temperature=298.15):
    """Explorador de la ecuación de Eyring: deslizadores de ``ΔG‡`` (5–30 kcal/mol) y ``T`` (273–323 K).

    Curva de ``k`` frente a ``ΔG‡`` (escala log) con referencias cotidianas y
    barras que ilustran la regla «cada ≈ 1.36 kcal/mol = ×10»; resume ``1/k`` y ``t½``.
    """
    controls = {
        "delta_g": _slider("ΔG‡: altura de la colina (kcal/mol)", delta_g, 5.0, 30.0, 0.5, ".1f"),
        "temperature": _slider("T: temperatura (K)", temperature, 273.0, 323.0, 1.0, ".0f"),
    }
    return _explorer(_dibujar_eyring, controls, "La ecuación de Eyring: de la colina a la velocidad")


# ---------------------------------------------------------------------------
# 6. Dependencia con la temperatura (ΔH‡, ΔS‡)
# ---------------------------------------------------------------------------
def _dibujar_temperatura(delta_h, delta_s):
    temps = np.linspace(273.15, 333.15, 9)
    rates = kin.eyring_rate(delta_h - temps * delta_s / 1000.0, temps)
    delta_g = delta_h - _T_REF * delta_s / 1000.0
    fit = {"delta_h_kcal": delta_h, "delta_s_cal": delta_s, "delta_g_kcal": delta_g, "r2": 1.0}

    fig, (ax1, ax2) = figure(13.5, 5.6, ncols=2)
    viz.plot_eyring(temps, rates, fit=fit, ax=ax1, title="")
    _panel_title(ax1, "Gráfico de Eyring: la pendiente es −ΔH‡/R")

    t_dense = np.linspace(273.15, 333.15, 300)
    k_dense = kin.eyring_rate(delta_h - t_dense * delta_s / 1000.0, t_dense)
    ax2.fill_between(t_dense - 273.15, 0, k_dense, color=viz._tint(COLORS["datos"], 0.9), zorder=1, linewidth=0)
    ax2.plot(t_dense - 273.15, k_dense, color=COLORS["datos"], linewidth=2.8, zorder=3)
    viz._kband(ax2, 36.0, 38.0, "cuerpo\nhumano", color=viz._tint(COLORS["ts"], 0.85), y_text=0.97)
    k_ref = float(kin.eyring_rate(delta_g, _T_REF))
    k_body = float(kin.eyring_rate(delta_h - 310.15 * delta_s / 1000.0, 310.15))
    viz._keypoint(ax2, 25.0, k_ref, COLORS["datos"], size=8)
    viz._keypoint(ax2, 37.0, k_body, COLORS["ts"], size=9)
    ax2.annotate(f"25 °C: {k_ref:.3g} s⁻¹", xy=(25.0, k_ref), xytext=(-10, 6), textcoords="offset points",
                 ha="right", va="bottom", color=INK, fontsize=11)
    ax2.annotate(f"37 °C: {k_body:.3g} s⁻¹  (×{k_body / k_ref:.2f})", xy=(37.0, k_body), xytext=(-10, 8),
                 textcoords="offset points", ha="right", va="bottom", color=INK, fontsize=11, fontweight="semibold")
    ax2.set_xlim(0, 60)
    ax2.set_ylim(0, float(k_dense.max()) * 1.12)
    viz._finish(ax2, "T (°C)", "k (s⁻¹)")
    _panel_title(ax2, "k(T): calentar acelera, y más cuanto mayor es ΔH‡")
    _titulo(fig, f"De 25 a 37 °C la reacción va ×{k_body / k_ref:.2f} más rápido",
            f"ΔH‡ = {delta_h:g} kcal/mol decide cuánto acelera el calor; ΔS‡ = {_num(delta_s, 0)} cal/(mol·K) mueve todo arriba o abajo.")

    q10 = float(kin.eyring_rate(delta_h - 308.15 * delta_s / 1000.0, 308.15)) / k_ref
    resumen = (
        f"ΔH‡ = {delta_h:g} kcal/mol y ΔS‡ = {_num(delta_s, 0)} cal/(mol·K) → ΔG‡(298 K) = ΔH‡ − TΔS‡ = {delta_g:.2f} kcal/mol, "
        f"k(25 °C) = {k_ref:.3g} s⁻¹ y k(37 °C) = {k_body:.3g} s⁻¹ (×{k_body / k_ref:.2f}); "
        f"Q10 (25→35 °C) = {q10:.2f}."
    )
    return fig, resumen


def explorar_temperatura(delta_h=15.0, delta_s=-10.0):
    """Explorador de la dependencia con la temperatura: deslizadores de ``ΔH‡`` (kcal/mol) y ``ΔS‡`` (cal/mol/K).

    Gráfico de Eyring ``ln(k/T)`` frente a ``1000/T`` (izquierda) y curva ``k(T)``
    (derecha) con marcas a 25 y 37 °C; resume ``ΔG‡(298 K)`` y el ``Q10``.
    """
    controls = {
        "delta_h": _slider("ΔH‡: energía a aportar (kcal/mol)", delta_h, 5.0, 30.0, 0.5, ".1f"),
        "delta_s": _slider("ΔS‡: orden a imponer (cal/(mol·K))", delta_s, -30.0, 30.0, 1.0, ".0f"),
    }
    return _explorer(_dibujar_temperatura, controls, "Temperatura: entalpía y entropía de activación")


# ---------------------------------------------------------------------------
# 7. Perfil de pH
# ---------------------------------------------------------------------------
_GK_PH_OPT = (8.5, 8.7)  # óptimo medido para la glucoquinasa humana (Šimčíková y Heneberg 2019)


def _dibujar_ph(pka1, pka2):
    ph = np.linspace(2.0, 12.0, 400)
    v = kin.bell_shaped_ph_profile(ph, 1.0, pka1, pka2)
    ph_opt = 0.5 * (pka1 + pka2)
    v_opt = float(kin.bell_shaped_ph_profile(ph_opt, 1.0, pka1, pka2))
    fig = viz.plot_ph_profile(
        ph, v, pkas=(pka1, pka2), ylabel="v₀ / Vmax",
        title=f"El óptimo está a mitad de camino entre los dos pKa: pH {ph_opt:.2f}",
        subtitle="v = Vmax / (1 + 10^(pKa₁ − pH) + 10^(pH − pKa₂)) · zonas grises: un grupo catalítico está «apagado»",
        pka_labels=("la base (Asp205)\nya tiene protón", "la lisina (Lys169)\nperdió su carga +"),
    )
    fig.set_size_inches(10.5, 5.8)
    ax = fig.axes[0]
    if pka2 - pka1 < 2.5:  # pKa próximos: la etiqueta de pKa₁ pasa a la izquierda de su guía
        for text in ax.texts:
            if text.get_text().startswith("pKa₁"):
                text.set_ha("right")
                text.xyann = (-5, -3)
    ax.axvspan(*_GK_PH_OPT, color=viz._tint(COLORS["producto"], 0.55), alpha=0.45, zorder=2, linewidth=0)
    ax.text(0.5 * sum(_GK_PH_OPT), 0.03, "glucoquinasa medida: pH 8.5–8.7", rotation=90, ha="center", va="bottom",
            color=INK_SECONDARY, fontsize=9.5, zorder=2)
    viz._keypoint(ax, ph_opt, v_opt, COLORS["ts"], size=11)
    ax.annotate(f"óptimo: pH {ph_opt:.2f} ({100 * v_opt:.0f} % de Vmax)", xy=(ph_opt, v_opt), xytext=(0, 12),
                textcoords="offset points", ha="center", va="bottom", color=INK, fontsize=11.5, fontweight="semibold")
    ax.set_xlim(2.0, 12.0)
    ax.set_ylim(0, 1.22)
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

    Dibuja ``v = Vmax/(1 + 10^(pKa₁ − pH) + 10^(pH − pKa₂))`` y marca el pH óptimo.
    """
    controls = {
        "pka1": _slider("pKa₁ (grupo que debe estar sin protón)", pka1, 3.0, 8.0, 0.1, ".1f"),
        "pka2": _slider("pKa₂ (grupo que debe tener protón)", pka2, 6.0, 11.0, 0.1, ".1f"),
    }
    return _explorer(_dibujar_ph, controls, "La campana de pH: dos grupos, dos condiciones")


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

    fig, ax = figure(10.5, 5.8)
    viz._kband(ax, 4.0, 7.0, "glucosa en sangre\nen ayunas (4–7 mM)", color=viz._tint(COLORS["ts"], 0.88), y_text=0.97)
    ax.fill_between(s, basal, activated, color=viz._tint(COLORS["hill"], 0.85), zorder=1, linewidth=0)
    ax.plot(s, basal, color=COLORS["michaelis_menten"], linewidth=2.8, zorder=3,
            label=f"sin activador (S₀.₅ = {_GK_S_HALF} mM, n = {_GK_N})")
    ax.plot(s, activated, color=COLORS["hill"], linewidth=2.8, zorder=3,
            label=f"con activador (Vmax ×{vmax_factor:g}, S₀.₅ = {s_half_act:.2g} mM)")
    viz._keypoint(ax, _GLUCOSE_BLOOD_MM, a_blood, COLORS["michaelis_menten"], size=9)
    viz._keypoint(ax, _GLUCOSE_BLOOD_MM, b_blood, COLORS["hill"], size=9)
    ax.annotate("", xy=(_GLUCOSE_BLOOD_MM, b_blood), xytext=(_GLUCOSE_BLOOD_MM, a_blood),
                arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=1.4, shrinkA=9, shrinkB=9), zorder=8)
    ax.annotate(f"{a_blood:.0f} %", xy=(_GLUCOSE_BLOOD_MM, a_blood), xytext=(12, -4), textcoords="offset points",
                ha="left", va="top", color=INK, fontsize=11.5)
    ax.annotate(f"{b_blood:.0f} %  (×{b_blood / a_blood:.1f})", xy=(_GLUCOSE_BLOOD_MM, b_blood), xytext=(12, 4),
                textcoords="offset points", ha="left", va="bottom", color=INK, fontsize=11.5, fontweight="semibold")
    ax.set_xlim(0, s_max)
    ax.set_ylim(0, max(100.0 * vmax_factor, 100.0) * 1.15)
    viz._finish(ax, "[glucosa] (mM)", "actividad (% de la Vmax basal)",
                f"A 5 mM de glucosa, el activador multiplica la actividad por {b_blood / a_blood:.1f}",
                "El activador baja S₀.₅ y sube Vmax: la curva se desplaza hacia arriba y a la izquierda.")
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
    activador, con la glucemia marcada; resume la actividad a 5 mM en ambos casos.
    """
    controls = {
        "vmax_factor": _slider("Vmax se multiplica por", vmax_factor, 1.0, 3.0, 0.1, ".1f"),
        "s_half_factor": _slider("S₀.₅ se multiplica por", s_half_factor, 0.1, 1.0, 0.05, ".2f"),
    }
    return _explorer(_dibujar_activador, controls, "Un activador alostérico de la glucoquinasa")


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
    lo = min(0.0, reaction_energy) - 4.0

    fig, ax = figure(10.5, 5.8)
    ax.fill_between(x, lo, e, color=viz._tint(INK_MUTED, 0.9), zorder=1, linewidth=0)
    ax.plot(x, e, color=INK_SECONDARY, linewidth=2.8, zorder=3)
    for xi, ei, color, text in ((0.0, 0.0, COLORS["reactivo"], "reactivos"), (0.5, barrier_eff, COLORS["ts"], "estado de transición"),
                                (1.0, reaction_energy, COLORS["producto"], "productos")):
        viz._keypoint(ax, xi, ei, color, size=11)
        ax.annotate(text, xy=(xi, ei), xytext=(0, 12), textcoords="offset points", ha="center", va="bottom",
                    color=INK, fontsize=11, fontweight="semibold")
    viz._guide(ax, "h", 0.0, start=0.0, end=0.62)
    viz._guide(ax, "h", barrier_eff, start=0.5, end=0.62, linestyle="-", linewidth=0.8, color=AXIS)
    viz._dimension_arrow(ax, 0.62, 0.0, barrier_eff, f"ΔE‡ = {_num(barrier_eff)} kcal/mol")
    viz._guide(ax, "h", 0.0, start=0.8, end=1.12)
    viz._guide(ax, "h", reaction_energy, start=1.0, end=1.12, linestyle="-", linewidth=0.8, color=AXIS)
    if abs(reaction_energy) > 0.05:
        viz._dimension_arrow(ax, 1.12, 0.0, reaction_energy, f"ΔE = {_num(reaction_energy)} kcal/mol")
    ax.set_xlim(-0.1, 1.55)
    ax.set_ylim(lo, barrier_eff + 6.0)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_xticklabels(["R", "TS‡", "P"])
    k = float(kin.eyring_rate(barrier_eff, _T_REF))
    k_rev = float(kin.eyring_rate(barrier_eff - reaction_energy, _T_REF))
    keq = float(np.exp(-reaction_energy / (kin.R_KCAL * _T_REF)))
    kind = "exergónica (ΔE < 0)" if reaction_energy < 0 else ("endergónica (ΔE > 0)" if reaction_energy > 0 else "termoneutra")
    cuesta = "cuesta abajo" if reaction_energy < 0 else ("cuesta arriba" if reaction_energy > 0 else "a nivel")
    viz._finish(ax, "Coordenada de reacción", "Energía relativa (kcal/mol)",
                f"La colina decide la velocidad (k ≈ {k:.3g} s⁻¹); el desnivel, hacia dónde va ({cuesta})",
                "Altura de la cima → rapidez. Diferencia entre valles → equilibrio.", grid_axis="y")

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

    Dibuja un perfil suave (spline cúbica con pendiente nula en R, TS y P) y
    resume ``k`` por Eyring y la constante de equilibrio ``exp(−ΔE/RT)``.
    """
    controls = {
        "barrier": _slider("ΔE‡: altura de la colina (kcal/mol)", barrier, 5.0, 30.0, 0.5, ".1f"),
        "reaction_energy": _slider("ΔE: desnivel entre valles (kcal/mol)", reaction_energy, -20.0, 10.0, 0.5, ".1f"),
    }
    return _explorer(_dibujar_perfil_energia, controls, "Dibuja tu propia colina de energía")


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
