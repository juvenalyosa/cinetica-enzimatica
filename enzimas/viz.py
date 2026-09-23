"""Visualización para el cuaderno *Cinética Enzimática* (glucoquinasa humana, QM/MM).

Módulo de gráficos consistente para Google Colab: perfiles de energía,
análisis de dinámica molecular, cinética (Michaelis–Menten, Hill,
linealizaciones, inhibición, simulación de mecanismos, Arrhenius/Eyring, pH) y
visores 3D con ``py3Dmol``.

Convenciones
------------
* Cada función ``plot_*`` devuelve la :class:`matplotlib.figure.Figure` (crea
  su propia figura salvo que se pase ``ax=``) y nunca llama a ``plt.show()``.
* Las funciones ``view_*`` devuelven un ``py3Dmol.view`` (la celda del cuaderno
  lo muestra si es la última expresión). ``py3Dmol`` se importa de forma
  perezosa dentro de esas funciones.
* Identificadores en inglés; textos, títulos y etiquetas de ejes en español.

Sistema de diseño
-----------------
Superficie ``#fcfcfb``; tinta primaria ``#0b0b0b``, secundaria ``#52514e``,
atenuada ``#898781``; rejilla ``#e1e0d9``; ejes ``#c3c2b7``. Los colores de
serie (:data:`PALETTE`) se asignan en orden fijo por entidad y nunca se
reciclan; la magnitud (p. ej. varias [I]) usa la rampa azul secuencial
(:func:`sequential_blue`) y la polaridad la rampa divergente azul↔rojo con
punto medio gris (:data:`DIVERGING_CMAP`). Los colores de estado
(:data:`STATUS`) están reservados y no se usan para series. El texto nunca
lleva el color de la serie: las anotaciones van en tinta y una marca de color a
su lado aporta la identidad.
"""

from __future__ import annotations

import os
from typing import Sequence

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from . import kinetics as _kin

__all__ = [
    # sistema de diseño
    "SURFACE", "INK", "INK_SECONDARY", "INK_MUTED", "GRID", "AXIS",
    "PALETTE", "COLORS", "STATUS", "SEQUENTIAL_BLUE", "SEQUENTIAL_CMAP",
    "DIVERGING_CMAP", "sequential_blue", "style_path", "apply_style", "figure",
    "save_figure",
    # energía / reacción
    "plot_energy_profile", "plot_energy_profiles", "plot_energy_levels",
    # dinámica molecular
    "plot_md_timeseries", "plot_md_summary", "rolling_mean", "plot_progress_curves", "plot_snapshot_scans",
    "plot_estimates", "plot_frequencies", "mostrar", "oscurecer",
    # cinética
    "plot_michaelis_menten", "plot_hill_vs_mm", "plot_linearizations",
    "plot_inhibition_family", "plot_inhibition_lineweaver", "plot_ode_simulation",
    "plot_initial_rates_from_ode", "plot_arrhenius", "plot_eyring",
    "plot_ph_profile", "plot_hill_plot",
    # 3D
    "view_complex", "view_qm_region", "view_frames", "mode_animation_frames",
    "xyz_block",
]

# ---------------------------------------------------------------------------
# Sistema de diseño
# ---------------------------------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

#: Colores categóricos en orden fijo (1 azul … 8 rojo). Nunca se reciclan.
PALETTE = [
    "#2a78d6",  # 1 azul
    "#eb6834",  # 2 naranja
    "#1baf7a",  # 3 aguamarina
    "#eda100",  # 4 amarillo
    "#e87ba4",  # 5 magenta
    "#008300",  # 6 verde
    "#4a3aa7",  # 7 violeta
    "#e34948",  # 8 rojo
]
_BLUE, _ORANGE, _AQUA, _YELLOW, _MAGENTA, _GREEN, _VIOLET, _RED = PALETTE

#: Roles semánticos → color de serie. Una entidad conserva su color en todas las figuras.
COLORS = {
    # perfil de reacción
    "reactivo": _BLUE,
    "ts": _ORANGE,
    "producto": _AQUA,
    "intermedio": _YELLOW,
    # mecanismo E + S ⇌ ES → E + P
    "enzima": _BLUE,
    "sustrato": _ORANGE,
    "complejo_es": _AQUA,
    "producto_p": _MAGENTA,
    # datos y modelos
    "datos": _BLUE,
    "michaelis_menten": _BLUE,
    "hill": _ORANGE,
    "ion_mg": _GREEN,
    "ion_k": _VIOLET,
    "sin_enzima": INK_MUTED,
}

#: Colores de estado (reservados: nunca para series).
STATUS = {"good": "#0ca30c", "warning": "#fab219", "critical": "#d03b3b"}

#: Rampa secuencial azul (claro = valor bajo, oscuro = valor alto).
SEQUENTIAL_BLUE = ("#cde2fb", "#0d366b")
SEQUENTIAL_CMAP = LinearSegmentedColormap.from_list("enzimas_azul", list(SEQUENTIAL_BLUE))
#: Rampa divergente azul ↔ rojo con punto medio gris neutro.
DIVERGING_CMAP = LinearSegmentedColormap.from_list("enzimas_divergente", [_BLUE, "#f0efec", _RED])

_STYLE_APPLIED = False


def sequential_blue(n, lo=0.3, hi=1.0):
    """``n`` colores de la rampa azul, de claro (valor bajo) a oscuro (valor alto).

    ``lo`` = 0.3 evita que el primer paso se confunda con la superficie.
    """
    n = int(n)
    if n <= 1:
        return [matplotlib.colors.to_hex(SEQUENTIAL_CMAP(hi))]
    return [matplotlib.colors.to_hex(SEQUENTIAL_CMAP(t)) for t in np.linspace(lo, hi, n)]


def style_path():
    """Ruta absoluta del archivo ``enzimas.mplstyle`` distribuido con el paquete."""
    try:
        from importlib import resources

        return str(resources.files("enzimas").joinpath("enzimas.mplstyle"))
    except Exception:  # pragma: no cover - versiones antiguas de Python
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "enzimas.mplstyle")


def _register_fonts():
    """Registra la tipografía Figtree (OFL) incluida en ``enzimas/fuentes``: se ve igual en Colab y en local.

    Figtree no trae letras griegas; matplotlib (≥ 3.6) las toma de DejaVu Sans, la siguiente de la lista.
    """
    from matplotlib import font_manager

    folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fuentes")
    for name in sorted(os.listdir(folder)) if os.path.isdir(folder) else []:
        if name.endswith(".ttf"):
            font_manager.fontManager.addfont(os.path.join(folder, name))


def apply_style():
    """Activa el estilo compartido del curso (``matplotlib.style.use``). Devuelve la ruta usada."""
    global _STYLE_APPLIED
    _register_fonts()
    path = style_path()
    matplotlib.style.use(path)
    _STYLE_APPLIED = True
    return path


def _ensure_style():
    if not _STYLE_APPLIED:
        apply_style()


def figure(width=7.0, height=4.2, nrows=1, ncols=1, **kwargs):
    """Crea ``(fig, ax)`` con el estilo del curso; ``width``/``height`` en pulgadas.

    Con ``nrows``/``ncols`` > 1 devuelve el arreglo de ejes de ``plt.subplots``.
    """
    _ensure_style()
    fig, ax = plt.subplots(nrows, ncols, figsize=(width, height), facecolor=SURFACE, **kwargs)
    return fig, ax


def save_figure(fig, path, dpi=200):
    """Guarda ``fig`` en ``path`` (crea la carpeta) con recorte ajustado y fondo de superficie."""
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=SURFACE, edgecolor=SURFACE)
    return path


# ---------------------------------------------------------------------------
# Utilidades internas
# ---------------------------------------------------------------------------
def _get_ax(ax, width=7.0, height=4.2):
    if ax is None:
        fig, ax = figure(width, height)
        return fig, ax
    _ensure_style()
    return ax.figure, ax


def _style_axes(ax, grid_axis="y"):
    """Aplica bordes, rejilla y colores aunque el ``ax`` venga de fuera del estilo."""
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.grid(False)
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID, linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_MUTED, labelcolor=INK_SECONDARY, width=0.8, length=3.5)
    ax.xaxis.label.set_color(INK_MUTED)
    ax.yaxis.label.set_color(INK_MUTED)


def _title(ax, title=None, subtitle=None):
    if title:
        ax.set_title(title, loc="left", pad=24 if subtitle else 12, color=INK, fontsize=13, fontweight="semibold")
    if subtitle:
        ax.annotate(
            subtitle, xy=(0, 1), xycoords="axes fraction", xytext=(0, 6), textcoords="offset points",
            ha="left", va="bottom", color=INK_SECONDARY, fontsize=10,
        )


def _fig_title(fig, title=None, subtitle=None, top=None):
    """Título y subtítulo alineados a la izquierda para figuras multipanel.

    ``top`` (fracción de la figura reservada a los paneles) se calcula a partir
    de la altura en pulgadas si no se indica.
    """
    height_in = float(fig.get_size_inches()[1])
    if top is None:
        reserve = 0.10 + (0.32 if title else 0.0) + (0.26 if subtitle else 0.0)
        top = 1.0 - reserve / height_in
    fig.tight_layout(rect=(0, 0, 1, top if (title or subtitle) else 1))
    if not (title or subtitle):
        return
    left = fig.subplotpars.left
    y = 0.985
    if title:
        fig.suptitle(title, x=left, y=y, ha="left", va="top", color=INK, fontsize=13, fontweight="semibold")
        y -= 0.30 / height_in
    if subtitle:
        fig.text(left, y, subtitle, ha="left", va="top", color=INK_SECONDARY, fontsize=10)


def _legend(ax, **kwargs):
    kwargs.setdefault("frameon", False)
    kwargs.setdefault("labelcolor", INK_SECONDARY)
    leg = ax.legend(**kwargs)
    if leg is not None and leg.get_title() is not None:
        leg.get_title().set_color(INK_SECONDARY)
    return leg


def _finish(ax, xlabel=None, ylabel=None, title=None, subtitle=None, grid_axis="y"):
    _style_axes(ax, grid_axis)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    _title(ax, title, subtitle)


def _markers(ax, x, y, color, label=None, size=6.5, zorder=4, hollow=False):
    return ax.plot(
        x, y, linestyle="none", marker="o", markersize=size, label=label, zorder=zorder,
        markerfacecolor=(SURFACE if hollow else color), markeredgecolor=(color if hollow else SURFACE),
        markeredgewidth=(1.4 if hollow else 1.2),
    )


def _guide(ax, orientation, value, start=None, end=None, color=INK_MUTED, **kwargs):
    kwargs.setdefault("linewidth", 1.0)
    kwargs.setdefault("linestyle", (0, (4, 3)))
    kwargs.setdefault("zorder", 1)
    if orientation == "h":
        if start is None and end is None:
            return ax.axhline(value, color=color, **kwargs)
        return ax.plot([start, end], [value, value], color=color, **kwargs)[0]
    if start is None and end is None:
        return ax.axvline(value, color=color, **kwargs)
    return ax.plot([value, value], [start, end], color=color, **kwargs)[0]


def _fmt(value, digits=3):
    return f"{float(value):.{digits}g}".replace("-", "−")


def _value_pm(value, err=None):
    """``valor ± error`` con decimales fijados por la magnitud del error (menos tipográfico)."""
    value = float(value)
    if err is None or not np.isfinite(err) or err <= 0:
        return _fmt(value)
    decimals = int(max(0, -np.floor(np.log10(err)) + 1))
    return f"{value:.{decimals}f} ± {float(err):.{decimals}f}".replace("-", "−")


def _num(value, decimals=1):
    """Número con ``decimals`` decimales y signo menos tipográfico (U+2212)."""
    return f"{float(value):.{decimals}f}".replace("-", "−")


def _smooth_curve(x, y, n=300):
    """Curva monótona por tramos (PCHIP) que no sobrepasa los puntos; segmentos si hay < 4 puntos."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 4:
        return x, y
    try:
        from scipy.interpolate import PchipInterpolator

        xs = np.linspace(x[0], x[-1], int(n))
        return xs, PchipInterpolator(x, y)(xs)
    except Exception:  # pragma: no cover
        return x, y


# ---------------------------------------------------------------------------
# 1. Energía y coordenada de reacción
# ---------------------------------------------------------------------------
def _dimension_arrow(ax, x, y0, y1, text, x_text_offset=6):
    """Flecha de cota vertical entre ``y0`` e ``y1`` con el valor en una píldora a la derecha."""
    ax.annotate(
        "", xy=(x, y1), xytext=(x, y0),
        arrowprops=dict(arrowstyle="<|-|>", color=INK, linewidth=1.2, shrinkA=0, shrinkB=0, mutation_scale=11),
        zorder=6,
    )
    ax.annotate(
        text, xy=(x, 0.5 * (y0 + y1)), xytext=(x_text_offset + 4, 0), textcoords="offset points",
        ha="left", va="center", color=INK, fontsize=12, fontweight="semibold", zorder=7,
        bbox=dict(boxstyle="round,pad=0.35,rounding_size=0.8", facecolor="#ffffff", edgecolor=GRID, linewidth=1.0),
    )


def _halo_point(ax, x, y, color, size=11, zorder=6, label=None):
    """Punto destacado con halo suave (estados R, TS, P)."""
    ax.plot([x], [y], "o", markersize=size * 2.1, color=color, alpha=0.16, markeredgewidth=0, zorder=zorder - 1)
    return ax.plot([x], [y], linestyle="none", marker="o", markersize=size, markerfacecolor=color,
                   markeredgecolor="#ffffff", markeredgewidth=2.0, zorder=zorder, label=label)[0]


def _state_label(ax, x, y, text, detail=None, dy=14, ha="center", va="bottom", color=INK):
    """Etiqueta directa de un estado: nombre en negrita y, debajo, un detalle en tinta secundaria."""
    s = text if not detail else f"{text}\n{detail}"
    t = ax.annotate(s, xy=(x, y), xytext=(0, dy if va == "bottom" else -dy), textcoords="offset points",
                    ha=ha, va=va, color=color, fontsize=11.5, fontweight="semibold", zorder=8, linespacing=1.25)
    return t


def _soft_fill(ax, x, y, base, color, alpha=0.10):
    """Relleno suave bajo una curva hasta ``base`` (sin borde)."""
    ax.fill_between(x, y, base, color=color, alpha=alpha, linewidth=0, zorder=1)


def plot_energy_profile(
    x, energy, xlabel="Coordenada de reacción", ax=None, ts_index=None, labels=None, title=None,
    subtitle=None, annotate_barrier=True, relative=True, smooth=True, color=None, label=None,
    ylabel="Energía relativa (kcal/mol)", annotate_states=False, state_names=("reactivo", "estado de transición", "producto"),
    annotate_reaction=False, sort=True,
):
    """Perfil de energía a lo largo de la coordenada de reacción (kcal/mol).

    ``energy`` se refiere al primer punto si ``relative=True``. El estado de transición (``ts_index``, por
    defecto el máximo) se resalta en naranja con halo y la barrera ``ΔE‡`` se anota con una flecha de cota
    desde el nivel de reactivos. Con ``annotate_states`` se rotulan directamente reactivo (primer punto,
    azul), TS (naranja) y producto (último punto, aguamarina); con ``annotate_reaction`` se acota también
    ``ΔE`` de reacción. ``labels`` es una lista opcional de textos por punto (``None`` para omitir).
    Con ``sort=False`` los puntos se dibujan en el orden dado (un camino de reacción cuyo ξ no es monótono:
    el primer punto es el reactivo y el último el producto).
    """
    x = np.asarray(x, dtype=float)
    e = np.asarray(energy, dtype=float)
    order = np.argsort(x, kind="stable") if sort else np.arange(x.size)
    x, e = x[order], e[order]
    if labels is not None:
        labels = [labels[i] for i in order]
    if ts_index is not None:
        ts_index = int(np.argsort(order)[int(ts_index)])
    if relative:
        e = e - e[0]
    if ts_index is None:
        ts_index = int(np.argmax(e))
    color = color or COLORS["reactivo"]

    fig, ax = _get_ax(ax, 9.6, 5.4)
    xs, ys = _smooth_curve(x, e) if (smooth and sort) else (x, e)
    base = min(float(np.min(e)), 0.0) - (0.16 if annotate_states else 0.08) * max(float(np.ptp(e)), 1.0)
    _soft_fill(ax, xs, ys, base, color, alpha=0.10)
    ax.plot(xs, ys, color=color, linewidth=2.8, label=label, zorder=3)
    if x.size <= 40:                       # en caminos muy densos los puntos ensucian: basta la línea
        ax.plot(x, e, linestyle="none", marker="o", markersize=5.5, markerfacecolor="#ffffff",
                markeredgecolor=color, markeredgewidth=1.6, zorder=4)
    ts_point = _markers(ax, [x[ts_index]], [e[ts_index]], COLORS["ts"], size=11, zorder=7)[0]
    ts_point.set_markeredgecolor("#ffffff")
    ts_point.set_markeredgewidth(2.0)
    ax.plot([x[ts_index]], [e[ts_index]], "o", markersize=24, color=COLORS["ts"], alpha=0.16,
            markeredgewidth=0, zorder=6)

    if labels is not None:
        for xi, ei, text in zip(x, e, labels):
            if text:
                ax.annotate(text, xy=(xi, ei), xytext=(0, 11), textcoords="offset points",
                            ha="center", va="bottom", color=INK_SECONDARY, fontsize=10.5)

    span = float(np.ptp(x)) if np.ptp(x) > 0 else 1.0
    erange = max(float(np.ptp(e)), 1.0)
    if annotate_states:
        _halo_point(ax, x[0], e[0], COLORS["reactivo"], size=10)
        _halo_point(ax, x[-1], e[-1], COLORS["producto"], size=10)
        ax.annotate(state_names[0], xy=(x[0], e[0]), xytext=(12, -10), textcoords="offset points", ha="left",
                    va="top", color=INK, fontsize=11.5, fontweight="semibold", zorder=8, linespacing=1.25)
        _state_label(ax, x[ts_index], e[ts_index], state_names[1], dy=16)
        below = x.size > 2 and float(np.max(e[max(0, x.size - 6):-1])) > e[-1]
        _state_label(ax, x[-1], e[-1], state_names[2], dy=16, ha="right", va="top" if below else "bottom")

    x_lo, x_hi = float(np.min(x)), float(np.max(x))
    right = x_hi + 0.03 * span
    if annotate_barrier:
        e_ref = e[0]
        x_ts = x[ts_index]
        x_arrow = x_ts + 0.05 * span
        _guide(ax, "h", e_ref, start=x[0], end=x_arrow)
        _guide(ax, "h", e[ts_index], start=x_ts, end=x_arrow, linestyle="-", linewidth=0.8, color=AXIS)
        _dimension_arrow(ax, x_arrow, e_ref, e[ts_index], f"ΔE‡ = {_num(e[ts_index] - e_ref)} kcal/mol")
        right = max(right, x_arrow + 0.32 * span)
    if annotate_reaction:
        x_r = x_hi + 0.05 * span
        _guide(ax, "h", e[0], start=x[0], end=x_r)
        ax.annotate("", xy=(x_r, e[-1]), xytext=(x_r, e[0]),
                    arrowprops=dict(arrowstyle="<|-|>", color=INK_SECONDARY, linewidth=1.0, shrinkA=0, shrinkB=0,
                                    mutation_scale=10), zorder=6)
        d_rxn = e[-1] - e[0]
        ax.annotate(f"ΔE = {'+' if d_rxn >= 0 else ''}{_num(d_rxn)} kcal/mol",
                    xy=(x_r, 0.5 * (e[0] + e[-1])), xytext=(7, 0), textcoords="offset points", ha="left",
                    va="center", color=INK_SECONDARY, fontsize=11)
        right = max(right, x_r + 0.22 * span)
    ax.set_xlim(x_lo - 0.04 * span, right)
    ax.set_ylim(base, float(np.max(e)) + 0.22 * erange)
    _finish(ax, xlabel, ylabel, title, subtitle)
    if label:
        _legend(ax)
    return fig


def plot_energy_profiles(
    profiles, xlabel="Coordenada de reacción", ax=None, title=None, subtitle=None, relative=True,
    smooth=True, mark_ts=True, ylabel="Energía relativa (kcal/mol)", colors=None, emphasis=None,
):
    """Varios perfiles ``(x, energía, etiqueta)`` superpuestos, rotulados al final de cada curva.

    ``emphasis`` (índice) dibuja ese perfil más grueso y con relleno; los demás quedan finos. La leyenda
    se mantiene arriba como índice de colores.
    """
    fig, ax = _get_ax(ax, 9.6, 5.4)
    colors = list(colors) if colors else PALETTE
    if len(profiles) > len(colors):
        raise ValueError("máximo 8 perfiles: agrupa o usa múltiplos pequeños")
    if emphasis is None:
        emphasis = len(profiles) - 1
    all_e = []
    ends = []
    for k, ((x, e, label), color) in enumerate(zip(profiles, colors)):
        x = np.asarray(x, dtype=float)
        e = np.asarray(e, dtype=float)
        order = np.argsort(x)
        x, e = x[order], e[order]
        if relative:
            e = e - e[0]
        all_e.append(e)
        xs, ys = _smooth_curve(x, e) if smooth else (x, e)
        strong = k == emphasis
        if strong:
            _soft_fill(ax, xs, ys, min(0.0, float(np.min(e))) - 1.0, color, alpha=0.09)
        ax.plot(xs, ys, color=color, linewidth=2.8 if strong else 1.8, label=label, zorder=3 if strong else 2,
                alpha=1.0 if strong else 0.9, linestyle="-" if strong else (0, (5, 3)))
        ax.plot(x, e, linestyle="none", marker="o", markersize=5 if strong else 4, markerfacecolor="#ffffff",
                markeredgecolor=color, markeredgewidth=1.5, zorder=4)
        if mark_ts:
            i = int(np.argmax(e))
            _markers(ax, [x[i]], [e[i]], color, size=10, hollow=not strong, zorder=6)
            if strong:
                ax.annotate(f"cima: {_num(e[i])} kcal/mol", xy=(x[i], e[i]), xytext=(0, 12), textcoords="offset points",
                            ha="center", va="bottom", fontsize=11, color=INK, fontweight="semibold")
        ends.append((float(x[-1]), float(e[-1]), label, color))
    lo = min(float(np.min(a)) for a in all_e)
    hi = max(float(np.max(a)) for a in all_e)
    ax.set_ylim(min(0.0, lo) - 1.0, hi + 0.25 * max(hi - lo, 1.0))
    x_max = max(float(np.max(p[0])) for p in profiles)
    x_min = min(float(np.min(p[0])) for p in profiles)
    ax.set_xlim(x_min - 0.03 * (x_max - x_min), x_max + 0.03 * (x_max - x_min))
    _finish(ax, xlabel, ylabel, title, subtitle)
    _legend(ax, loc="upper left", ncol=len(profiles), handlelength=2.2)
    return fig


def plot_energy_levels(
    levels, ax=None, title=None, subtitle=None, connect=True, compare=None, label="con enzima",
    compare_label="sin enzima", annotate_barrier=True, show_values=True, ts_indices=None,
    ylabel="Energía (kcal/mol)", width=0.56,
):
    """Diagrama de niveles de energía por etapas.

    ``levels`` = lista de ``(etiqueta, energía_kcal)``, p. ej.
    ``[("E + S", 0), ("ES", -3), ("TS", 15), ("EP", -5), ("E + P", -8)]``.
    Los estados de transición (etiqueta con "TS" o "‡", o ``ts_indices``) se dibujan en naranja; el primer
    nivel en azul (reactivo) y el último en aguamarina (producto). Una curva suave une los niveles como un
    perfil de reacción. La barrera se anota desde el mínimo inmediatamente anterior. ``compare`` es una
    segunda lista (p. ej. sin enzima) dibujada en gris discontinuo con las mismas posiciones.
    """
    names = [str(name) for name, _ in levels]
    e = np.asarray([float(val) for _, val in levels])
    n = len(names)
    if ts_indices is None:
        ts_indices = [i for i, name in enumerate(names) if "TS" in name.upper() or "‡" in name]
    ts_set = set(int(i) for i in ts_indices)
    fig, ax = _get_ax(ax, 9.6, 5.4)
    half = width / 2.0

    def level_color(i):
        if i in ts_set:
            return COLORS["ts"]
        if i == n - 1 and n > 1:
            return COLORS["producto"]
        return COLORS["reactivo"]

    def curve(values, color, lw, style, alpha, zorder):
        # perfil suave: tramos horizontales en cada nivel unidos por curvas en S
        xs, ys = [], []
        for i, val in enumerate(values):
            xs += [i - half, i + half]
            ys += [val, val]
        xs, ys = np.asarray(xs), np.asarray(ys)
        fine_x, fine_y = [], []
        for k in range(len(xs) - 1):
            t = np.linspace(0, 1, 30)
            if k % 2 == 0:
                fine_x.append(xs[k] + (xs[k + 1] - xs[k]) * t)
                fine_y.append(np.full_like(t, ys[k]))
            else:
                s = 0.5 - 0.5 * np.cos(np.pi * t)
                fine_x.append(xs[k] + (xs[k + 1] - xs[k]) * t)
                fine_y.append(ys[k] + (ys[k + 1] - ys[k]) * s)
        ax.plot(np.concatenate(fine_x), np.concatenate(fine_y), color=color, linewidth=lw, linestyle=style,
                alpha=alpha, zorder=zorder, solid_capstyle="round")

    if compare is not None:
        ce = np.asarray([float(val) for _, val in compare])
        if connect:
            curve(ce, COLORS["sin_enzima"], 1.6, (0, (5, 3)), 0.9, 2)
        for i, val in enumerate(ce):
            ax.hlines(val, i - half, i + half, color=COLORS["sin_enzima"], linewidth=3.0, zorder=3, alpha=0.8)
            if show_values and not (i < len(e) and abs(val - e[i]) < 0.05):
                ax.annotate(_num(val), xy=(i + half, val), xytext=(4, 0), textcoords="offset points",
                            ha="left", va="center", color=INK_MUTED, fontsize=10)
    if connect:
        curve(e, INK_MUTED, 1.2, (0, (2, 3)), 0.9, 2)
    for i, val in enumerate(e):
        col = level_color(i)
        ax.fill_between([i - half, i + half], [val, val], [val - 0.02 * max(np.ptp(e), 1)] * 2, color=col, alpha=0.0)
        ax.hlines(val, i - half, i + half, color=col, linewidth=5.0, zorder=5, capstyle="round")

    if show_values:
        for i, val in enumerate(e):
            # si el nivel de comparación queda justo encima, el valor va debajo de la barra
            crowded = compare is not None and i < len(ce) and 0 < ce[i] - val < 0.15 * max(float(np.ptp(e)), 1.0)
            ax.annotate(_num(val), xy=(i, val), xytext=(0, -10 if crowded else 8), textcoords="offset points",
                        ha="center", va="top" if crowded else "bottom", color=INK, fontsize=12, fontweight="semibold")

    if annotate_barrier and ts_set:
        for i_ts in sorted(ts_set):
            if i_ts == 0:
                continue
            i_ref = int(np.argmin(e[:i_ts]))
            e_ref = e[i_ref]
            x_arrow = i_ts + half + 0.14
            _guide(ax, "h", e_ref, start=i_ref + half, end=x_arrow)
            _guide(ax, "h", e[i_ts], start=i_ts + half, end=x_arrow, linestyle="-", linewidth=0.8, color=AXIS)
            _dimension_arrow(ax, x_arrow, e_ref, e[i_ts], f"ΔE‡ = {_num(e[i_ts] - e_ref)} kcal/mol", x_text_offset=5)

    ax.set_xticks(range(n), names)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, n - 0.5 + (0.35 if (annotate_barrier and ts_set) else 0.0))
    lo = float(np.min(e)) if compare is None else min(float(np.min(e)), float(np.min(ce)))
    hi = float(np.max(e)) if compare is None else max(float(np.max(e)), float(np.max(ce)))
    ax.set_ylim(lo - 0.12 * max(hi - lo, 1), hi + 0.18 * max(hi - lo, 1))
    _finish(ax, None, ylabel, title, subtitle)
    ax.tick_params(axis="x", length=0, labelcolor=INK, labelsize=12)
    if compare is not None:
        handles = [Line2D([], [], color=COLORS["reactivo"], linewidth=4, label=label),
                   Line2D([], [], color=COLORS["sin_enzima"], linewidth=3, linestyle=(0, (5, 3)), label=compare_label)]
        _legend(ax, handles=handles, loc="upper right")
    return fig


# ---------------------------------------------------------------------------
# 2. Dinámica molecular
# ---------------------------------------------------------------------------
_MD_LABELS = {
    "temperatura_K": "Temperatura (K)",
    "rmsd_CA_A": "RMSD Cα (Å)",
    "d_PG_O6_A": "d(Pγ–O6) (Å)",
    "d_PG_O3B_A": "d(Pγ–O3β) (Å)",
    "d_O6_OD1asp205_A": "d(O6–OD1 Asp205) (Å)",
    "d_Mg_PG_A": "d(Mg–Pγ) (Å)",
    "d_NZlys169_O6_A": "d(NZ Lys169–O6) (Å)",
}


def rolling_mean(y, window=7):
    """Media móvil centrada (bordes con ventana recortada): tendencia de una serie ruidosa de MD."""
    y = np.asarray(y, dtype=float)
    window = max(1, int(window))
    kernel = np.ones(window)
    num = np.convolve(y, kernel, mode="same")
    den = np.convolve(np.ones_like(y), kernel, mode="same")
    return num / den


def _band(ax, y0, y1, color, text=None, alpha=0.13, text_x=0.99, text_color=None):
    """Banda horizontal de referencia con rótulo interior a la derecha."""
    ax.axhspan(y0, y1, color=color, alpha=alpha, linewidth=0, zorder=0)
    if text:
        ax.annotate(text, xy=(text_x, 0.5 * (y0 + y1)), xycoords=("axes fraction", "data"), ha="right", va="center",
                    fontsize=10.5, color=text_color or INK_SECONDARY, fontweight="semibold", zorder=8)


def plot_md_timeseries(
    df, columns, labels=None, ylabel="", ax=None, time_column="tiempo_ps", xlabel="Tiempo (ps)",
    title=None, subtitle=None, colors=None, smooth_window=7,
):
    """Series temporales de una dinámica molecular a partir de un ``DataFrame`` con ``tiempo_ps``.

    Cada serie se dibuja tenue (los fotogramas) con su media móvil encima (la tendencia). ``columns`` es una
    lista de nombres de columna; ``labels`` (opcional) sus etiquetas en español.
    """
    if isinstance(columns, str):
        columns = [columns]
    columns = list(columns)
    if labels is None:
        labels = [_MD_LABELS.get(c, c) for c in columns]
    colors = list(colors) if colors else PALETTE
    if len(columns) > len(colors):
        raise ValueError("máximo 8 series por panel: usa múltiplos pequeños")
    fig, ax = _get_ax(ax, 9.6, 4.6)
    t = np.asarray(df[time_column], dtype=float)
    for col, lab, color in zip(columns, labels, colors):
        y = np.asarray(df[col], dtype=float)
        ax.plot(t, y, color=color, linewidth=1.0, alpha=0.35, zorder=2)
        ax.plot(t, rolling_mean(y, smooth_window), color=color, linewidth=2.6, label=lab, zorder=3)
    ax.set_xlim(t.min(), t.max())
    _finish(ax, xlabel, ylabel, title, subtitle)
    if len(columns) > 1:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax + 0.18 * (ymax - ymin))
        _legend(ax, loc="upper left", ncol=min(len(columns), 3))
    return fig


def plot_md_summary(
    df, crystal_value=None, title=None, subtitle=None,
    time_column="tiempo_ps", crystal_label="cristal", bins=30, nac_threshold=3.5,
    crystal_hbond=None, hbond_max=3.2, smooth_window=7,
):
    """Resumen de la MD en una sola historia: ¿está la enzima lista para reaccionar?

    Arriba, grande: d(Pγ–O6) frente al tiempo con la banda de ataque cercano (d < ``nac_threshold``), la
    media móvil y el valor del cristal, y su distribución pegada a la derecha. Abajo, pequeños: el RMSD
    (estabilidad de la proteína) y d(O6–OD1 Asp205) con la banda de puente de hidrógeno (< ``hbond_max``).
    Si ``title`` es None se escribe uno con el porcentaje de tiempo en ataque cercano. Devuelve 4 ejes.
    """
    _ensure_style()
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(10.4, 7.6), facecolor=SURFACE)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.45, 1.0], width_ratios=[4.2, 1.0], hspace=0.55, wspace=0.10)
    ax_d = fig.add_subplot(gs[0, 0])
    ax_h = fig.add_subplot(gs[0, 1], sharey=ax_d)
    gs_low = gs[1, :].subgridspec(1, 2, wspace=0.28)
    ax_r = fig.add_subplot(gs_low[0, 0])
    ax_a = fig.add_subplot(gs_low[0, 1])
    t = np.asarray(df[time_column], dtype=float)
    green = COLORS["producto"]

    frac = None
    if "d_PG_O6_A" in df.columns:
        d = np.asarray(df["d_PG_O6_A"], dtype=float)
        frac = float(np.mean(d < nac_threshold))
        lo = min(float(np.min(d)), crystal_value if crystal_value is not None else float(np.min(d))) - 0.15
        hi = float(np.max(d)) + 0.18
        ax_d.axhspan(lo - 1, nac_threshold, color=green, alpha=0.10, linewidth=0, zorder=0)
        ax_d.axhline(nac_threshold, color=green, linewidth=1.2, linestyle=(0, (4, 3)), zorder=1)
        ax_d.plot(t, d, color=COLORS["reactivo"], linewidth=1.0, alpha=0.35, zorder=2)
        ax_d.plot(t, d, linestyle="none", marker="o", markersize=3.2, color=COLORS["reactivo"], alpha=0.45, zorder=2)
        ax_d.plot(t, rolling_mean(d, smooth_window), color=COLORS["reactivo"], linewidth=2.8, zorder=3)
        y_lab = 0.5 * (nac_threshold + (crystal_value if crystal_value is not None else lo))
        ax_d.annotate(f"zona de ataque cercano  (d < {nac_threshold:g} Å)", xy=(0.012, y_lab),
                      xycoords=("axes fraction", "data"), ha="left", va="center", fontsize=11.5,
                      color=oscurecer(green), fontweight="semibold", zorder=8)
        if crystal_value is not None:
            ax_d.axhline(float(crystal_value), color=INK, linewidth=1.3, zorder=2)
            ax_d.annotate(f"{crystal_label}: {float(crystal_value):.2f} Å", xy=(0.012, float(crystal_value)),
                          xycoords=("axes fraction", "data"), xytext=(0, 6), textcoords="offset points",
                          ha="left", va="bottom", fontsize=10.5, color=INK, zorder=8)
        ax_d.set_ylim(lo, hi)
        ax_d.set_xlim(t.min(), t.max())
        _finish(ax_d, "Tiempo (ps)", "d(Pγ–O6)  (Å)")
        ax_d.set_title("Distancia entre el fósforo γ del ATP y el O6 de la glucosa", loc="left", fontsize=12,
                       color=INK_SECONDARY, fontweight="normal", pad=8)

        # distribución marginal, horizontal, en dos colores (dentro / fuera de la zona)
        edges = np.linspace(lo, hi, bins + 1)
        counts, _ = np.histogram(d, bins=edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        cols = [green if c < nac_threshold else INK_MUTED for c in centers]
        ax_h.barh(centers, counts, height=(edges[1] - edges[0]) * 0.86, color=cols, alpha=0.85, zorder=2)
        ax_h.axhline(nac_threshold, color=green, linewidth=1.2, linestyle=(0, (4, 3)))
        ax_h.annotate(f"{frac:.0%}", xy=(0.5, 0.5 * (lo + nac_threshold)), xycoords=("axes fraction", "data"),
                      ha="center", va="center", fontsize=20, fontweight="bold", color=oscurecer(green))
        ax_h.annotate("del tiempo", xy=(0.5, 0.5 * (lo + nac_threshold)), xycoords=("axes fraction", "data"),
                      xytext=(0, -18), textcoords="offset points", ha="center", va="center", fontsize=10,
                      color=INK_SECONDARY)
        ax_h.set_xlim(0, max(counts.max(), 1) * 1.15)
        _style_axes(ax_h, grid_axis=None)
        ax_h.spines["left"].set_visible(False)
        ax_h.tick_params(axis="y", left=False, labelleft=False)
        ax_h.set_xlabel("fotogramas")
    else:
        ax_d.text(0.5, 0.5, "sin columna «d_PG_O6_A»", transform=ax_d.transAxes, ha="center", va="center", color=INK_MUTED)
        _style_axes(ax_h, grid_axis=None)

    small = [
        (ax_r, "rmsd_CA_A", "La proteína es estable", "RMSD de los Cα (Å)", COLORS["enzima"], None),
        (ax_a, "d_O6_OD1asp205_A", "O6 ··· Asp205 (la base catalítica)", "d(O6–OD1)  (Å)", PALETTE[6], hbond_max),
    ]
    for ax, col, panel_title, ylab, color, band in small:
        if col in df.columns:
            y = np.asarray(df[col], dtype=float)
            if band is not None:
                lo_b = min(float(np.min(y)), crystal_hbond or band) - 0.3
                ax.axhspan(lo_b - 5, band, color=green, alpha=0.10, linewidth=0, zorder=0)
                ax.axhline(band, color=green, linewidth=1.0, linestyle=(0, (4, 3)))
                ax.annotate("puente de H", xy=(0.02, band), xycoords=("axes fraction", "data"), xytext=(0, -6),
                            textcoords="offset points", ha="left", va="top", fontsize=10, color=oscurecer(green),
                            fontweight="semibold")
                if crystal_hbond is not None:
                    ax.axhline(float(crystal_hbond), color=INK, linewidth=1.1)
                    ax.annotate(f"{crystal_label}: {float(crystal_hbond):.1f} Å", xy=(0.98, float(crystal_hbond)),
                                xycoords=("axes fraction", "data"), xytext=(0, 5), textcoords="offset points",
                                ha="right", va="bottom", fontsize=10, color=INK)
                ax.set_ylim(lo_b, float(np.max(y)) + 0.25)
            ax.plot(t, y, color=color, linewidth=1.0, alpha=0.35)
            ax.plot(t, rolling_mean(y, smooth_window), color=color, linewidth=2.4)
            m = float(np.mean(y))
            ax.annotate(f"media {m:.2f} Å", xy=(0.98, 0.95), xycoords="axes fraction", ha="right", va="top",
                        fontsize=10.5, color=INK, fontweight="semibold")
            if band is None:
                ax.set_ylim(0, max(float(np.max(y)) * 1.35, 1.0))
            ax.set_xlim(t.min(), t.max())
        else:
            ax.text(0.5, 0.5, f"sin columna «{col}»", transform=ax.transAxes, ha="center", va="center", color=INK_MUTED)
        _finish(ax, "Tiempo (ps)", ylab)
        ax.set_title(panel_title, loc="left", fontsize=12, color=INK_SECONDARY, fontweight="normal", pad=8)

    if title is None:
        title = (f"La enzima pasa el {frac:.0%} del tiempo lista para reaccionar" if frac is not None
                 else "Resumen de la dinámica molecular")
    if subtitle is None:
        subtitle = "Cada punto es un fotograma de la película (cada 10 ps); la línea gruesa es la media móvil"
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.08, top=0.855)
    fig.suptitle(title, x=0.08, y=0.975, ha="left", va="top", color=INK, fontsize=16, fontweight="semibold")
    fig.text(0.08, 0.935, subtitle, ha="left", va="top", color=INK_SECONDARY, fontsize=11.5)
    return fig


def plot_frequencies(sets, title=None, subtitle=None, noise=40.0, xlabel="Frecuencia (cm⁻¹)"):
    """Las frecuencias vibracionales más bajas como barras: las imaginarias (negativas) a la izquierda.

    ``sets`` = lista de ``(etiqueta, frecuencias)``; una fila por conjunto. Las imaginarias grandes
    (|ν| > ``noise``) van en naranja (el movimiento de la reacción); las pequeñas, en gris claro con la
    banda de «ruido numérico» sombreada.
    """
    if isinstance(sets, np.ndarray) or (sets and not isinstance(sets[0], (tuple, list))):
        sets = [("", sets)]
    n = len(sets)
    fig, ax = _get_ax(None, 9.6, 1.5 * n + 2.0)
    ax.axvspan(-noise, noise, color=INK_MUTED, alpha=0.10, linewidth=0, zorder=0)
    ax.annotate("ruido\nnumérico", xy=(0, 0.03), xycoords=("data", "axes fraction"), ha="center", va="bottom",
                fontsize=10, color=INK_MUTED, linespacing=1.1)
    ax.axvline(0, color=AXIS, linewidth=1.0, zorder=1)
    all_f = []
    for row, (lab, freqs) in enumerate(sets):
        y = n - 1 - row
        freqs = np.asarray(freqs, dtype=float)
        all_f += list(freqs)
        n_big = 0
        for f in sorted(freqs):
            big_imag = f < -noise
            color = COLORS["ts"] if big_imag else (AXIS if f < 0 else COLORS["reactivo"])
            ax.plot([f, f], [y - 0.28, y + 0.28], color=color, linewidth=6 if big_imag else 4, solid_capstyle="round",
                    alpha=1.0 if (big_imag or f > 0) else 0.9, zorder=3)
            if big_imag:
                up = n_big % 2 == 0
                ax.annotate(f"{_num(f, 0)} cm⁻¹", xy=(f, y + (0.30 if up else -0.30)), xytext=(0, 4 if up else -4),
                            textcoords="offset points", ha="center", va="bottom" if up else "top", fontsize=11,
                            fontweight="semibold", color=INK)
                n_big += 1
    ax.set_yticks(range(n), [lab for lab, _ in sets][::-1])
    ax.tick_params(axis="y", length=0, labelcolor=INK, labelsize=11.5)
    lo, hi = min(all_f + [-60]), max(all_f + [60])
    ax.set_xlim(lo - 0.12 * (hi - lo), hi + 0.08 * (hi - lo))
    ax.set_ylim(-0.7, n - 0.1)
    ax.annotate("← imaginarias: la geometría «cae»", xy=(0.0, -0.13), xycoords="axes fraction", ha="left",
                va="top", fontsize=10.5, color=INK_SECONDARY)
    ax.annotate("reales: la geometría vibra →", xy=(1.0, -0.13), xycoords="axes fraction", ha="right",
                va="top", fontsize=10.5, color=INK_SECONDARY)
    _finish(ax, None, None, title, subtitle, grid_axis="x")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelcolor=INK, labelsize=11.5)
    ax.set_xlabel(xlabel, labelpad=26)
    return fig


def oscurecer(color, t=0.35):
    """Versión más oscura de un color (para texto que acompaña a una banda del mismo tono)."""
    r, g, b = matplotlib.colors.to_rgb(color)
    return matplotlib.colors.to_hex((r * (1 - t), g * (1 - t), b * (1 - t)))


def plot_progress_curves(curves, t_tangent=None, title=None, subtitle=None, xlabel="Tiempo (s)",
                         ylabel="Producto formado (µM)"):
    """Curvas de progreso ``[P](t)`` con su tangente inicial (v₀) rotulada directamente.

    ``curves`` = lista de ``(t, P, v0, etiqueta)``, de menor a mayor [S]₀ (rampa azul). Cada tangente llega
    hasta el 70 % del producto final (como mucho ``t_tangent``, por defecto el 45 % del eje de tiempo), así
    se ve dónde la curva real se separa de ella.
    """
    fig, ax = _get_ax(None, 9.6, 5.4)
    colors = sequential_blue(len(curves), lo=0.45)
    y_top = 0.0
    t_cap = t_tangent if t_tangent is not None else 0.45 * float(np.max(curves[0][0]))
    handles = []
    for (t, p, v0, lab), color in zip(curves, colors):
        t = np.asarray(t, dtype=float)
        p = np.asarray(p, dtype=float)
        t_tan = min(0.7 * float(np.max(p)) / v0, t_cap) if v0 > 0 else t_cap
        ax.fill_between(t, p, 0, color=color, alpha=0.05, linewidth=0)
        ax.plot(t, p, color=color, linewidth=2.8, zorder=3)
        ax.plot([0, t_tan], [0, t_tan * v0], color=INK, linewidth=1.3, linestyle=(0, (4, 3)), zorder=4, alpha=0.75)
        handles.append(Line2D([], [], color=color, linewidth=3.0, label=f"{lab}:  v₀ = {v0:.2f} µM/s"))
        ax.annotate(lab, xy=(t[-1], p[-1]), xytext=(-4, 8), textcoords="offset points", ha="right", va="bottom",
                    fontsize=11, color=INK_SECONDARY)
        y_top = max(y_top, float(np.max(p)), t_tan * v0)
    ax.set_xlim(0, float(np.max(curves[0][0])))
    ax.set_ylim(0, y_top * 1.12)
    _finish(ax, xlabel, ylabel, title, subtitle)
    leg = _legend(ax, handles=handles[::-1], loc="upper left", title="Pendiente inicial (línea discontinua)",
                  handlelength=2.4, borderaxespad=0.8)
    leg.get_title().set_fontweight("semibold")
    leg.get_title().set_color(INK)
    leg._legend_box.align = "left"
    return fig


def plot_snapshot_scans(reference, scans, xlabel="ξ = d(Pγ–O3β) − d(Pγ–O6)  (Å)", title=None, subtitle=None,
                        reference_label="cristal minimizado", ts_index=None):
    """El camino con estado de transición (``reference`` = (ξ, E)) frente a escaneos que suben sin parar.

    ``scans`` = lista de ``(ξ, E, etiqueta)``; se dibujan en rojo claro con su valor final rotulado.
    """
    fig, ax = _get_ax(None, 9.6, 5.6)
    xr, er = (np.asarray(a, dtype=float) for a in reference)
    order = np.argsort(xr)
    xr, er = xr[order], er[order]
    reds = [matplotlib.colors.to_hex(c) for c in
            LinearSegmentedColormap.from_list("r", ["#f4a9a8", _RED])(np.linspace(0.15, 1.0, max(len(scans), 1)))]
    top = float(np.max(er))
    ends = []
    for (x, e, lab), color in zip(scans, reds):
        x = np.asarray(x, dtype=float)
        e = np.asarray(e, dtype=float)
        ax.plot(x, e, color=color, linewidth=2.0, zorder=3, marker="o", markersize=4.2,
                markerfacecolor="#ffffff", markeredgecolor=color, markeredgewidth=1.3)
        ends.append((float(x[-1]), float(e[-1]), lab, color))
        top = max(top, float(np.max(e)))
    # rótulos al final, separados verticalmente para que no choquen
    ends.sort(key=lambda r: r[1])
    min_gap = 0.055 * top
    placed = []
    for x_end, e_end, lab, color in ends:
        y = e_end if not placed else max(e_end, placed[-1] + min_gap)
        placed.append(y)
        ax.annotate(f"{lab}: {e_end:.0f}", xy=(x_end, e_end), xytext=(x_end + 0.08, y), textcoords="data",
                    ha="left", va="center", fontsize=10.5, color=INK_SECONDARY,
                    arrowprops=dict(arrowstyle="-", color=color, linewidth=0.8, shrinkA=0, shrinkB=2))
    xs, ys = _smooth_curve(xr, er)
    _soft_fill(ax, xs, ys, min(0.0, float(np.min(er))) - 1.5, COLORS["reactivo"], alpha=0.10)
    ax.plot(xs, ys, color=COLORS["reactivo"], linewidth=3.0, zorder=4)
    i_ts = int(np.argmax(er)) if ts_index is None else int(ts_index)
    _markers(ax, [xr[i_ts]], [er[i_ts]], COLORS["ts"], size=11, zorder=7)
    ax.plot([xr[i_ts]], [er[i_ts]], "o", markersize=24, color=COLORS["ts"], alpha=0.16, markeredgewidth=0, zorder=6)
    ax.annotate(f"cima: {er[i_ts]:.1f} kcal/mol", xy=(xr[i_ts], er[i_ts]), xytext=(0, -22), textcoords="offset points",
                ha="center", va="top", fontsize=11.5, fontweight="semibold", color=INK,
                bbox=dict(boxstyle="round,pad=0.3,rounding_size=0.7", facecolor="#ffffff", edgecolor=GRID))
    _halo_point(ax, xr[-1], er[-1], COLORS["producto"], size=10)
    ax.annotate(f"{reference_label}:\nhay producto estable", xy=(xr[-1], er[-1]), xytext=(0, -16),
                textcoords="offset points", ha="center", va="top", fontsize=10.5, color=INK, fontweight="semibold")
    ax.set_xlim(min(float(np.min(xr)), min(float(np.min(s[0])) for s in scans)) - 0.1,
                max(float(np.max(xr)), max(float(np.max(s[0])) for s in scans)) + 0.85)
    ax.set_ylim(min(0.0, float(np.min(er))) - 1.5, top * 1.12)
    _finish(ax, xlabel, "Energía relativa al reactivo (kcal/mol)", title, subtitle)
    handles = [Line2D([], [], color=COLORS["reactivo"], linewidth=3, label=f"{reference_label} (camino completo)"),
               Line2D([], [], color=_RED, linewidth=2, marker="o", markersize=4, markerfacecolor="#ffffff",
                      label="instantáneas de la MD (sin producto)")]
    _legend(ax, handles=handles, loc="upper left")
    return fig


def plot_estimates(rows, reference=None, reference_label="experimento", xlabel="Barrera (kcal/mol)",
                   title=None, subtitle=None, rate_fn=None):
    """Comparación de estimaciones de una barrera: una fila por estimación (punto + valor), el experimento
    como banda vertical, y (si se da ``rate_fn``) la velocidad k que implica cada valor, en texto a la derecha.

    ``rows`` = lista de ``(etiqueta, valor, color)``.
    """
    n = len(rows)
    fig, ax = _get_ax(None, 9.6, 0.62 * n + 2.2)
    vals = [float(r[1]) for r in rows]
    lo, hi = min(vals + ([reference] if reference is not None else [])), max(vals + ([reference] if reference is not None else []))
    pad = 0.12 * (hi - lo + 1)
    if reference is not None:
        ax.axvspan(reference - 0.5, reference + 0.5, color=INK_MUTED, alpha=0.12, linewidth=0, zorder=0)
        ax.axvline(reference, color=INK, linewidth=1.3, zorder=1)
        ax.annotate(f"{reference_label}: {reference:.1f}", xy=(reference, 0.0), xycoords=("data", "axes fraction"),
                    xytext=(8, 6), textcoords="offset points", ha="left", va="bottom", fontsize=11,
                    fontweight="semibold", color=INK)
    for i, (lab, val, color) in enumerate(rows):
        y = n - 1 - i
        if reference is not None:
            ax.plot([reference, val], [y, y], color=color, linewidth=2.0, alpha=0.5, zorder=2)
        _markers(ax, [val], [y], color, size=12, zorder=4)
        ax.annotate(f"{val:.1f}", xy=(val, y), xytext=(0, 11), textcoords="offset points", ha="center", va="bottom",
                    fontsize=11, fontweight="semibold", color=INK)
        if rate_fn is not None:
            k = rate_fn(val)
            ax.annotate(f"k ≈ {_sci(k)} s⁻¹", xy=(1.0, y), xycoords=("axes fraction", "data"), xytext=(10, 0),
                        textcoords="offset points", ha="left", va="center", fontsize=11, color=INK_SECONDARY,
                        annotation_clip=False)
    ax.set_yticks(range(n), [r[0] for r in rows][::-1])
    ax.tick_params(axis="y", length=0, labelcolor=INK, labelsize=11.5)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(lo - pad, hi + pad)
    _finish(ax, xlabel, None, title, subtitle, grid_axis="x")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelcolor=INK, labelsize=11.5)
    ax.set_ylim(-1.1, n - 0.3)
    fig.subplots_adjust(right=0.80)
    return fig


def _sci(x):
    """Número en notación científica legible: 6.2 × 10¹²."""
    if x == 0 or not np.isfinite(x):
        return _fmt(x)
    exp = int(np.floor(np.log10(abs(x))))
    if -2 <= exp <= 3:
        return f"{x:.3g}"
    sup = str(exp).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{x / 10 ** exp:.1f} × 10{sup}"


# ---------------------------------------------------------------------------
# 3. Cinética
# ---------------------------------------------------------------------------
def _tint(color, t=0.85):
    """Mezcla ``color`` con blanco (t = 1 → blanco): rellenos y bandas suaves del mismo tono que la serie."""
    rgb = np.array(matplotlib.colors.to_rgb(color))
    return matplotlib.colors.to_hex(rgb + (1.0 - rgb) * float(t))


def _kbox(ax, x, y, text, ha="right", va="bottom", transform=None, fontsize=11, color=None):
    """Recuadro blanco redondeado con texto en tinta (parámetros de un ajuste, lectura de una gráfica)."""
    return ax.text(
        x, y, text, transform=transform or ax.transAxes, ha=ha, va=va, color=color or INK_SECONDARY,
        fontsize=fontsize, linespacing=1.45, zorder=8,
        bbox=dict(boxstyle="round,pad=0.55,rounding_size=0.6", facecolor="white", edgecolor=GRID, linewidth=1.0),
    )


def _kband(ax, x0, x1, label=None, color=None, alpha=1.0, y_text=0.975, fontsize=10.5, anchor="center"):
    """Banda vertical sombreada con su rótulo arriba (zonas: saturación, rango fisiológico…).

    ``anchor="left"`` alinea el rótulo al borde izquierdo de la banda (bandas estrechas junto al eje).
    """
    ax.axvspan(x0, x1, color=color or _tint(INK_MUTED, 0.86), alpha=alpha, zorder=0, linewidth=0)
    if label:
        x = x0 if anchor == "left" else 0.5 * (x0 + x1)
        ax.annotate(label, xy=(x, y_text), xycoords=ax.get_xaxis_transform(), xytext=(7 if anchor == "left" else 0, 0),
                    textcoords="offset points", ha="left" if anchor == "left" else "center", va="top",
                    color=INK_MUTED, fontsize=fontsize, linespacing=1.3, zorder=1)


def _keypoint(ax, x, y, color, size=10, zorder=7):
    """Punto clave: marcador grande con aro blanco y halo del mismo color."""
    ax.plot([x], [y], "o", markersize=size * 2.1, color=color, alpha=0.16, markeredgewidth=0, zorder=zorder - 1)
    return ax.plot([x], [y], "o", markersize=size, markerfacecolor=color, markeredgecolor="white",
                   markeredgewidth=2.0, zorder=zorder)[0]


def plot_michaelis_menten(
    s, v, fit=None, ax=None, show_km=True, title=None, subtitle=None, xlabel="[S] (mM)",
    ylabel="v₀ (µM/s)", show_params=True, point_label=None, curve_label=None, color=None,
):
    """Datos ``(s, v)`` y curva de Michaelis–Menten ajustada, con las tres zonas de la hipérbola.

    ``fit`` es el diccionario de :func:`enzimas.kinetics.fit_michaelis_menten`
    (o ``None`` para dibujar solo los puntos). Con ``show_km`` se sombrean la zona
    casi lineal (``[S] < Km``) y la de saturación (``[S] > 4·Km``), se marcan
    ``Vmax`` y el punto ``(Km, Vmax/2)``.
    """
    s = np.asarray(s, dtype=float)
    v = np.asarray(v, dtype=float)
    color = color or COLORS["datos"]
    fig, ax = _get_ax(ax)
    if fit is not None:
        vmax, km = float(fit["vmax"]), float(fit["km"])
        s_dense = np.linspace(0.0, float(np.max(s)) * 1.05, 400)
        v_dense = _kin.michaelis_menten(s_dense, vmax, km)
        ax.fill_between(s_dense, 0, v_dense, color=_tint(color, 0.9), zorder=1, linewidth=0)
        ax.plot(s_dense, v_dense, color=color, linewidth=2.6, label=curve_label, zorder=3)
        ax.set_ylim(0, max(vmax * 1.22, float(np.max(v)) * 1.08))
        ax.set_xlim(0, s_dense[-1])
        if show_km:
            if km < 0.45 * s_dense[-1]:
                _kband(ax, 0, km, "casi lineal:\ncada molécula cuenta", color=_tint(color, 0.93), anchor="left")
            if 4 * km < 0.8 * s_dense[-1]:
                _kband(ax, 4 * km, s_dense[-1], "saturación: añadir más\nsustrato casi no ayuda")
            _guide(ax, "h", vmax, color=INK_SECONDARY)
            ax.annotate("Vmax", xy=(1.0, vmax), xycoords=("axes fraction", "data"), xytext=(-4, 4),
                        textcoords="offset points", ha="right", va="bottom", color=INK, fontsize=11.5, fontweight="semibold")
            _guide(ax, "h", vmax / 2, start=0, end=km)
            _guide(ax, "v", km, start=0, end=vmax / 2)
            _keypoint(ax, km, vmax / 2, COLORS["ts"])
            if km / s_dense[-1] > 0.1:  # hay sitio junto al eje
                ax.annotate("Vmax/2", xy=(0, vmax / 2), xytext=(5, 4), textcoords="offset points",
                            ha="left", va="bottom", color=INK_SECONDARY, fontsize=11)
            else:  # Km pequeño: a la derecha del punto, bajo la curva
                ax.annotate("Vmax/2", xy=(km, vmax / 2), xytext=(14, -6), textcoords="offset points",
                            ha="left", va="top", color=INK_SECONDARY, fontsize=11)
            ax.annotate("Km", xy=(km, 0), xytext=(5, 4), textcoords="offset points",
                        ha="left", va="bottom", color=INK, fontsize=11.5, fontweight="semibold")
        if show_params:
            text = f"Vmax = {_value_pm(vmax, fit.get('vmax_err'))}\nKm = {_value_pm(km, fit.get('km_err'))}"
            if "r2" in fit:
                text += f"\nR² = {fit['r2']:.4f}"
            _kbox(ax, 0.98, 0.06, text)
    _markers(ax, s, v, color, label=point_label, size=8.5, zorder=5)
    _finish(ax, xlabel, ylabel, title, subtitle)
    if point_label or curve_label:
        _legend(ax, loc="center right")
    return fig


def plot_hill_vs_mm(
    s, v_mm, v_hill, n_hill=None, ax=None, title=None, subtitle=None, xlabel="[S] (mM)",
    ylabel="v₀ (µM/s)", s_half=None, labels=None,
):
    """Hipérbola de Michaelis–Menten (n = 1) frente a sigmoide de Hill (n > 1), con la diferencia sombreada.

    La zona sombreada entre curvas muestra dónde el interruptor (sigmoide) responde distinto del
    regulador (hipérbola): por debajo de ``S₀.₅`` la sigmoide está «apagada»; por encima, más «encendida».
    """
    s = np.asarray(s, dtype=float)
    v_mm = np.asarray(v_mm, dtype=float)
    v_hill = np.asarray(v_hill, dtype=float)
    if labels is None:
        n_text = f"n = {float(n_hill):.2g}" if n_hill is not None else "n > 1"
        labels = ("Michaelis–Menten (n = 1)", f"Hill ({n_text})")
    fig, ax = _get_ax(ax)
    ax.fill_between(s, v_mm, v_hill, color=_tint(COLORS["hill"], 0.85), zorder=1, linewidth=0)
    ax.plot(s, v_mm, color=COLORS["michaelis_menten"], linewidth=2.6, label=labels[0], zorder=3)
    ax.plot(s, v_hill, color=COLORS["hill"], linewidth=2.6, label=labels[1], zorder=4)
    if s_half is not None:
        s_half = float(s_half)
        y_half = float(np.interp(s_half, s, v_hill))
        _guide(ax, "v", s_half, start=0, end=y_half)
        _keypoint(ax, s_half, y_half, COLORS["hill"], size=9)
        ax.annotate("S₀.₅", xy=(s_half, 0.0), xytext=(5, 4), textcoords="offset points",
                    ha="left", va="bottom", color=INK, fontsize=11.5, fontweight="semibold")
        top = float(max(np.max(v_mm), np.max(v_hill)))
        below = s < s_half
        if below.any():
            k = int(np.argmax(np.where(below, v_mm - v_hill, -np.inf)))
            if v_mm[k] - v_hill[k] > 0.06 * top:
                ax.annotate("por debajo de S₀.₅\nla sigmoide casi no responde", xy=(s[k], 0.5 * (v_mm[k] + v_hill[k])),
                            xytext=(0.36, 0.2), textcoords="axes fraction", ha="left", va="center", color=INK_SECONDARY,
                            fontsize=10.5, arrowprops=dict(arrowstyle="-", color=INK_MUTED, linewidth=0.9, shrinkB=4))
    ax.set_xlim(0, s.max())
    ax.set_ylim(0, float(max(np.max(v_mm), np.max(v_hill))) * 1.12)
    _finish(ax, xlabel, ylabel, title, subtitle)
    _legend(ax, loc="lower right")
    return fig


_LINEARIZATIONS = (
    ("lineweaver_burk", "Lineweaver–Burk", "1/[S]", "1/v₀", "1/Vmax", "−1/Km"),
    ("eadie_hofstee", "Eadie–Hofstee", "v₀/[S]", "v₀", "Vmax", "Vmax/Km"),
    ("hanes_woolf", "Hanes–Woolf", "[S]", "[S]/v₀", "Km/Vmax", "−Km"),
)


def plot_linearizations(s, v, fits=None, title="Linealizaciones de Michaelis–Menten", subtitle=None, color=None):
    """Tres paneles: Lineweaver–Burk, Eadie–Hofstee y Hanes–Woolf.

    ``fits`` es un diccionario opcional ``{nombre: linear_fit_dict}`` con claves
    ``lineweaver_burk``, ``eadie_hofstee`` y ``hanes_woolf``; si falta alguna se
    calcula con :func:`enzimas.kinetics.linear_fit`. Las rectas se extrapolan
    (trazo discontinuo) hasta los cortes con los ejes, que se marcan y rotulan. En
    Lineweaver–Burk se señalan los puntos de [S] baja, que dominan (y distorsionan) el ajuste.
    """
    s = np.asarray(s, dtype=float)
    v = np.asarray(v, dtype=float)
    fits = dict(fits or {})
    color = color or COLORS["datos"]
    fig, axes = figure(14.0, 5.0, nrows=1, ncols=3)
    transforms = {"lineweaver_burk": _kin.lineweaver_burk, "eadie_hofstee": _kin.eadie_hofstee, "hanes_woolf": _kin.hanes_woolf}
    for ax, (key, name, xl, yl, y_int_label, x_int_label) in zip(axes, _LINEARIZATIONS):
        x, y = transforms[key](s, v)
        fit = fits.get(key) or _kin.linear_fit(x, y)
        slope, intercept = float(fit["slope"]), float(fit["intercept"])
        x_int = -intercept / slope if slope != 0 else np.nan
        lo = min(0.0, float(np.min(x)), x_int if np.isfinite(x_int) else 0.0)
        hi = max(0.0, float(np.max(x)), x_int if np.isfinite(x_int) else 0.0)
        pad = 0.14 * (hi - lo if hi > lo else 1.0)
        lo, hi = lo - pad, hi + pad
        ax.axhline(0, color=AXIS, linewidth=0.9, zorder=1)
        ax.axvline(0, color=AXIS, linewidth=0.9, zorder=1)
        xd = np.array([lo, hi])
        ax.plot(xd, slope * xd + intercept, color=color, linewidth=1.3, linestyle=(0, (4, 3)), zorder=2)
        xr = np.array([float(np.min(x)), float(np.max(x))])
        ax.plot(xr, slope * xr + intercept, color=color, linewidth=2.6, zorder=3)
        if key == "lineweaver_burk" and x.size >= 4:
            worst = np.argsort(x)[-2:]  # [S] más baja → 1/[S] más grande
            for k in worst:
                ax.plot([x[k]], [y[k]], "o", markersize=19, markerfacecolor="none", markeredgecolor=COLORS["ts"],
                        markeredgewidth=1.6, zorder=6)
            k = int(worst[-1])
            ax.annotate("[S] baja: pocos puntos\nque pesan demasiado", xy=(x[k], y[k]), xytext=(-12, 10),
                        textcoords="offset points", ha="right", va="bottom", color=INK_SECONDARY, fontsize=10)
        _markers(ax, x, y, color, size=8)
        _keypoint(ax, 0.0, intercept, INK_SECONDARY, size=7)
        ax.annotate(y_int_label, xy=(0.0, intercept), xytext=(-8, 6), textcoords="offset points",
                    ha="right", va="bottom", color=INK, fontsize=11, fontweight="semibold")
        if np.isfinite(x_int):
            _keypoint(ax, x_int, 0.0, INK_SECONDARY, size=7)
            ax.annotate(x_int_label, xy=(x_int, 0.0), xytext=(0, -10), textcoords="offset points",
                        ha="center", va="top", color=INK, fontsize=11, fontweight="semibold")
        ax.set_xlim(lo, hi)
        ax.margins(y=0.18)
        _finish(ax, xl, yl)
        ax.set_title(name, loc="left", fontsize=13, color=INK, fontweight="semibold", pad=22)
        ax.annotate(f"R² = {fit['r2']:.3f}", xy=(0, 1), xycoords="axes fraction", xytext=(0, 5),
                    textcoords="offset points", ha="left", va="bottom", color=INK_MUTED, fontsize=10.5)
    _fig_title(fig, title, subtitle)
    return fig


_INHIBITION_NAMES = {
    "competitive": "Inhibición competitiva",
    "uncompetitive": "Inhibición acompetitiva",
    "noncompetitive": "Inhibición no competitiva",
    "mixed": "Inhibición mixta",
}
_INHIBITION_LB_PATTERN = {
    "competitive": "las rectas se cruzan en el eje 1/v₀: Vmax no cambia, Km aparente aumenta",
    "uncompetitive": "rectas paralelas: Vmax y Km aparentes disminuyen en la misma proporción",
    "noncompetitive": "las rectas se cruzan en el eje 1/[S]: Km no cambia, Vmax disminuye",
    "mixed": "las rectas se cruzan a la izquierda del eje 1/v₀ (fuera de los ejes)",
}


def _inhibitor_label(i, unit):
    i = float(i)
    return f"[I] = {i:g} {unit}" + (" (control)" if i == 0 else "")


def plot_inhibition_family(
    s, curves, kind="competitive", ax=None, title=None, subtitle=None, xlabel="[S] (mM)",
    ylabel="v₀ (µM/s)", unit="mM", show_points=False,
):
    """Familia de curvas ``v(s)`` a varias concentraciones de inhibidor.

    ``curves`` = lista de ``(concentración_I, v_array)``; se colorean con la
    rampa azul secuencial (claro = poco inhibidor, oscuro = mucho). La velocidad
    que se pierde entre el control y la [I] más alta queda sombreada.
    """
    s = np.asarray(s, dtype=float)
    curves = sorted(((float(i), np.asarray(v, dtype=float)) for i, v in curves), key=lambda c: c[0])
    colors = sequential_blue(len(curves))
    fig, ax = _get_ax(ax)
    if len(curves) > 1:
        ax.fill_between(s, curves[-1][1], curves[0][1], color=_tint(_RED, 0.9),
                        zorder=1, linewidth=0)
        k = int(np.argmin(np.abs(s - s.max() * 0.62)))
        ax.annotate("velocidad perdida\npor el inhibidor", xy=(s[k], 0.5 * (curves[0][1][k] + curves[-1][1][k])),
                    ha="center", va="center", color=INK_SECONDARY, fontsize=10.5, zorder=6)
    for (i, v), color in zip(curves, colors):
        ax.plot(s, v, color=color, linewidth=2.6, label=_inhibitor_label(i, unit), zorder=3)
        if show_points:
            _markers(ax, s, v, color)
    ax.set_xlim(0, s.max())
    ax.set_ylim(bottom=0)
    if title is None:
        title = _INHIBITION_NAMES.get(str(kind).lower(), "Inhibición")
    _finish(ax, xlabel, ylabel, title, subtitle)
    _legend(ax, loc="lower right")
    return fig


def plot_inhibition_lineweaver(
    s, curves, kind="competitive", ax=None, title=None, subtitle=None, unit="mM", extend=1.25,
):
    """Lineweaver–Burk para una familia de inhibición: resalta el punto donde se cruzan las rectas.

    Cada recta se ajusta por mínimos cuadrados y se extrapola (trazo
    discontinuo) hasta más allá de su corte con el eje ``1/[S]``. El cruce entre
    el control y la [I] más alta (la «huella» del tipo de inhibición) se marca con un halo.
    """
    s = np.asarray(s, dtype=float)
    curves = sorted(((float(i), np.asarray(v, dtype=float)) for i, v in curves), key=lambda c: c[0])
    colors = sequential_blue(len(curves))
    fig, ax = _get_ax(ax)
    x_ints, x_max, y_max = [], 0.0, 0.0
    lines = []
    for (i, v), color in zip(curves, colors):
        x, y = _kin.lineweaver_burk(s, v)
        fit = _kin.linear_fit(x, y)
        lines.append((i, x, y, fit, color))
        if fit["slope"] != 0:
            x_ints.append(-fit["intercept"] / fit["slope"])
        x_max = max(x_max, float(np.max(x)))
        y_max = max(y_max, float(np.max(y)))
    x_min = min(x_ints) * extend if x_ints else -0.2 * x_max
    x_min = min(x_min, -0.15 * x_max)
    ax.axhline(0, color=AXIS, linewidth=0.9, zorder=1)
    ax.axvline(0, color=AXIS, linewidth=0.9, zorder=1)
    for i, x, y, fit, color in lines:
        xd = np.array([x_min, x_max * 1.05])
        ax.plot(xd, fit["slope"] * xd + fit["intercept"], color=color, linewidth=1.3, linestyle=(0, (4, 3)), zorder=2)
        xr = np.array([float(np.min(x)), float(np.max(x))])
        ax.plot(xr, fit["slope"] * xr + fit["intercept"], color=color, linewidth=2.6, zorder=3, label=_inhibitor_label(i, unit))
        _markers(ax, x, y, color, size=8)
    if len(lines) > 1:
        f0, f1 = lines[0][3], lines[-1][3]
        ds = f1["slope"] - f0["slope"]
        if abs(ds) > 1e-9 * max(abs(f0["slope"]), 1e-12):
            xc = (f0["intercept"] - f1["intercept"]) / ds
            yc = f0["slope"] * xc + f0["intercept"]
            if x_min <= xc <= x_max * 1.05:
                _keypoint(ax, xc, yc, COLORS["ts"], size=9)
                ax.annotate("aquí se cruzan", xy=(xc, yc), xytext=(-16, -30), textcoords="offset points",
                            ha="right", va="top", color=INK, fontsize=10.5, fontweight="semibold")
        else:
            ax.annotate("rectas paralelas: nunca se cruzan", xy=(0.03, 0.97), xycoords="axes fraction",
                        ha="left", va="top", color=INK, fontsize=10.5, fontweight="semibold")
    ax.set_xlim(x_min, x_max * 1.05)
    ax.set_ylim(top=y_max * 1.15)
    kind = str(kind).lower()
    if title is None:
        title = f"{_INHIBITION_NAMES.get(kind, 'Inhibición')}: Lineweaver–Burk"
    if subtitle is None:
        subtitle = _INHIBITION_LB_PATTERN.get(kind)
    _finish(ax, "1/[S]", "1/v₀", title, subtitle)
    _legend(ax, loc="upper left" if kind != "uncompetitive" else "lower right")
    return fig


def plot_ode_simulation(
    sim, title=r"Simulación del mecanismo E + S $\rightleftharpoons$ ES $\rightarrow$ E + P", subtitle=None, t_zoom=None,
    xlabel="Tiempo (s)", ylabel="Concentración (µM)",
):
    """Concentraciones E, S, ES y P frente al tiempo (diccionario de :func:`kinetics.simulate_mechanism`).

    Dos paneles con un solo eje cada uno: a la izquierda S y P (escala grande), rotulados sobre la curva;
    a la derecha E y ES (escala pequeña), con ``t_zoom`` para ampliar el estado pre-estacionario (límite
    superior del tiempo en el panel derecho), donde se marca la meseta de estado estacionario de ES.
    """
    t = np.asarray(sim["t"], dtype=float)
    fig, (ax1, ax2) = figure(11.6, 4.8, nrows=1, ncols=2)

    def direct(ax, x, y, text, color, where=0.62, dy=8):
        i = int(where * (len(x) - 1))
        ax.annotate(text, xy=(x[i], y[i]), xytext=(0, dy), textcoords="offset points", ha="center",
                    va="bottom" if dy > 0 else "top", fontsize=11.5, fontweight="semibold", color=INK,
                    bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.6", facecolor="#ffffff", edgecolor=color,
                              linewidth=1.2), zorder=6)

    S, P = np.asarray(sim["S"], dtype=float), np.asarray(sim["P"], dtype=float)
    ax1.fill_between(t, P, 0, color=COLORS["producto_p"], alpha=0.10, linewidth=0)
    ax1.plot(t, S, color=COLORS["sustrato"], linewidth=2.8)
    ax1.plot(t, P, color=COLORS["producto_p"], linewidth=2.8)
    direct(ax1, t, S, "S · sustrato", COLORS["sustrato"], where=0.30, dy=10)
    direct(ax1, t, P, "P · producto", COLORS["producto_p"], where=0.70, dy=-12)
    ax1.set_xlim(t.min(), t.max())
    ax1.set_ylim(0, max(float(S.max()), float(P.max())) * 1.1)
    _finish(ax1, xlabel, ylabel)
    ax1.set_title("El sustrato se gasta, el producto aparece", loc="left", fontsize=12.5, color=INK, pad=8)

    E, ES = np.asarray(sim["E"], dtype=float), np.asarray(sim["ES"], dtype=float)
    t_max = float(t_zoom) if t_zoom is not None else float(t.max())
    mask = t <= t_max
    ax2.fill_between(t[mask], ES[mask], 0, color=COLORS["complejo_es"], alpha=0.12, linewidth=0)
    ax2.plot(t, E, color=COLORS["enzima"], linewidth=2.8)
    ax2.plot(t, ES, color=COLORS["complejo_es"], linewidth=2.8)
    tm, Em, ESm = t[mask], E[mask], ES[mask]
    if len(tm) > 3:
        direct(ax2, tm, Em, "E · enzima libre", COLORS["enzima"], where=0.72, dy=8)
        direct(ax2, tm, ESm, "ES · complejo", COLORS["complejo_es"], where=0.72, dy=-12)
    if t_zoom is not None:
        ax2.set_xlim(0, float(t_zoom))
        ax2.set_title("En milisegundos, ES llega a una meseta", loc="left", fontsize=12.5, color=INK, pad=8)
    else:
        ax2.set_xlim(t.min(), t.max())
        ax2.set_title("Enzima libre y complejo ES", loc="left", fontsize=12.5, color=INK, pad=8)
    ax2.set_ylim(0, max(float(E.max()), float(ES.max())) * 1.18)
    _finish(ax2, xlabel, ylabel)
    _fig_title(fig, title, subtitle)
    return fig


def plot_initial_rates_from_ode(
    s_values, v0, fit=None, ax=None, title="Michaelis–Menten emerge del mecanismo", subtitle=None,
    xlabel="[S]₀ (µM)", ylabel="v₀ (µM/s)",
):
    """Velocidades iniciales obtenidas de la simulación (puntos) y ajuste de Michaelis–Menten (curva)."""
    return plot_michaelis_menten(
        s_values, v0, fit=fit, ax=ax, title=title, subtitle=subtitle, xlabel=xlabel, ylabel=ylabel,
        point_label="v₀ de la simulación (EDO)", curve_label="ajuste de Michaelis–Menten",
    )


def _slope_triangle(ax, x0, x1, slope, intercept, text, below=True):
    """Triángulo de pendiente entre ``x0`` y ``x1`` sobre la recta ``y = slope·x + intercept``."""
    y0, y1 = slope * x0 + intercept, slope * x1 + intercept
    corner = (x1, y0) if below else (x0, y1)
    ax.plot([x0, corner[0], x1], [y0, corner[1], y1], color=INK_MUTED, linewidth=1.1, zorder=4)
    ax.fill([x0, corner[0], x1], [y0, corner[1], y1], color=_tint(INK_MUTED, 0.9), zorder=1, linewidth=0)
    ax.annotate(text, xy=corner, xytext=(8, 0) if below else (-8, 0), textcoords="offset points",
                ha="left" if below else "right", va="center", color=INK_SECONDARY, fontsize=10.5)


def _celsius_axis(ax, transform, inverse, label="T (°C)"):
    """Eje superior en °C para gráficos de 1/T (Arrhenius, Eyring): el estudiante ve temperaturas reales."""
    sec = ax.secondary_xaxis("top", functions=(transform, inverse))
    sec.set_xlabel(label, color=INK_MUTED, fontsize=11)
    sec.tick_params(colors=INK_MUTED, labelcolor=INK_SECONDARY, width=0.8, length=3.5)
    sec.spines["top"].set_color(AXIS)
    return sec


def plot_arrhenius(temperatures_k, rates, fit=None, ax=None, title="Gráfico de Arrhenius", subtitle=None, color=None):
    """``ln k`` frente a ``1000/T``; ``fit`` es el diccionario de :func:`kinetics.fit_arrhenius`.

    La pendiente (−Ea/R) se dibuja como un triángulo y el eje superior muestra la temperatura en °C.
    """
    t = np.asarray(temperatures_k, dtype=float)
    k = np.asarray(rates, dtype=float)
    color = color or COLORS["datos"]
    x = 1000.0 / t
    fig, ax = _get_ax(ax)
    if fit is not None:
        pad = 0.04 * (x.max() - x.min())
        xd = np.linspace(x.min() - pad, x.max() + pad, 50)
        slope = -float(fit["ea_kcal"]) / _kin.R_KCAL / 1000.0
        intercept = np.log(float(fit["a"]))
        ax.plot(xd, slope * xd + intercept, color=color, linewidth=2.6, zorder=2)
        xa, xb = x.min() + 0.25 * (x.max() - x.min()), x.min() + 0.6 * (x.max() - x.min())
        _slope_triangle(ax, xa, xb, slope, intercept, "pendiente = −Ea/R", below=True)
        _kbox(ax, 0.97, 0.95, f"Ea = {_num(fit['ea_kcal'])} kcal/mol\nR² = {fit['r2']:.4f}", ha="right", va="top")
    _markers(ax, x, np.log(k), color, size=8.5, zorder=5)
    _finish(ax, "1000/T (K⁻¹)", "ln k", title, subtitle)
    _celsius_axis(ax, lambda u: 1000.0 / np.clip(u, 1e-9, None) - 273.15, lambda c: 1000.0 / (np.asarray(c) + 273.15))
    return fig


def plot_eyring(temperatures_k, rates, fit=None, ax=None, title="Gráfico de Eyring", subtitle=None, color=None):
    """``ln(k/T)`` frente a ``1/T``; ``fit`` es el diccionario de :func:`kinetics.fit_eyring`.

    La pendiente da ΔH‡ (−ΔH‡/R) y la ordenada en el origen ΔS‡; el eje superior muestra °C.
    """
    t = np.asarray(temperatures_k, dtype=float)
    k = np.asarray(rates, dtype=float)
    color = color or COLORS["datos"]
    x = 1000.0 / t
    fig, ax = _get_ax(ax)
    if fit is not None:
        pad = 0.04 * (x.max() - x.min())
        xd = np.linspace(x.min() - pad, x.max() + pad, 50)
        intercept = np.log(_kin.KB_J / _kin.H_J) + float(fit["delta_s_cal"]) / (_kin.R_KCAL * 1000.0)
        slope = -float(fit["delta_h_kcal"]) / _kin.R_KCAL / 1000.0
        ax.plot(xd, slope * xd + intercept, color=color, linewidth=2.6, zorder=2)
        xa, xb = x.min() + 0.25 * (x.max() - x.min()), x.min() + 0.6 * (x.max() - x.min())
        _slope_triangle(ax, xa, xb, slope, intercept, "pendiente = −ΔH‡/R", below=True)
        text = (f"ΔH‡ = {_num(fit['delta_h_kcal'])} kcal/mol\nΔS‡ = {_num(fit['delta_s_cal'])} cal/(mol·K)\n"
                f"ΔG‡(298 K) = {_num(fit['delta_g_kcal'])} kcal/mol")
        _kbox(ax, 0.97, 0.95, text, ha="right", va="top")
    _markers(ax, x, np.log(k / t), color, size=8.5, zorder=5)
    _finish(ax, "1000/T (K⁻¹)", "ln(k/T)", title, subtitle)
    _celsius_axis(ax, lambda u: 1000.0 / np.clip(u, 1e-9, None) - 273.15, lambda c: 1000.0 / (np.asarray(c) + 273.15))
    return fig


def plot_ph_profile(
    ph, v, ax=None, title="Perfil de pH", subtitle=None, pkas=None, xlabel="pH", ylabel="v₀ (µM/s)",
    color=None, fit_curve=None, pka_labels=None,
):
    """Actividad frente al pH (puntos; curva continua si se pasa ``fit_curve=(ph_denso, v_denso)``).

    ``pkas=(pKa₁, pKa₂)`` añade guías verticales rotuladas y sombrea los flancos
    ácido y básico; ``pka_labels=(texto_ácido, texto_básico)`` rotula esos flancos.
    """
    ph = np.asarray(ph, dtype=float)
    v = np.asarray(v, dtype=float)
    color = color or COLORS["datos"]
    fig, ax = _get_ax(ax)
    x_lo, x_hi = float(ph.min()), float(ph.max())
    if pkas:
        labels = pka_labels or (None, None)
        _kband(ax, x_lo, float(pkas[0]), labels[0], y_text=0.80)
        _kband(ax, float(pkas[-1]), x_hi, labels[-1], y_text=0.80)
    curve = None
    if fit_curve is not None:
        curve = (np.asarray(fit_curve[0], dtype=float), np.asarray(fit_curve[1], dtype=float))
    elif ph.size > 40:
        curve = (ph, v)
    if curve is not None:
        ax.fill_between(curve[0], 0, curve[1], color=_tint(color, 0.88), zorder=1, linewidth=0)
        ax.plot(curve[0], curve[1], color=color, linewidth=2.6, zorder=3)
    if ph.size <= 40:
        _markers(ax, ph, v, color, size=8.5, zorder=5)
    if pkas:
        for idx, pka in enumerate(pkas, start=1):
            _guide(ax, "v", float(pka))
            ax.annotate(f"pKa{'₁₂₃'[idx - 1]} = {float(pka):.1f}", xy=(float(pka), 1.0), xycoords=("data", "axes fraction"),
                        xytext=(5, -3), textcoords="offset points", ha="left", va="top", color=INK, fontsize=11,
                        fontweight="semibold")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(bottom=0)
    _finish(ax, xlabel, ylabel, title, subtitle)
    return fig


def plot_hill_plot(s, v, vmax, ax=None, title="Gráfico de Hill", subtitle=None, color=None, unit="mM"):
    """Gráfico de Hill ``log₁₀(v/(Vmax − v))`` frente a ``log₁₀[S]`` con la pendiente ``n`` anotada.

    Incluye, como referencia, la recta de pendiente 1 (sin cooperatividad) que pasa por el mismo ``S₀.₅``:
    cuanto más empinados los datos respecto de esa recta, más cooperativa la enzima.
    """
    x, y = _kin.hill_plot(s, v, float(vmax))
    fit = _kin.linear_fit(x, y)
    color = color or COLORS["datos"]
    fig, ax = _get_ax(ax)
    pad = 0.05 * (x.max() - x.min())
    xd = np.array([x.min() - pad, x.max() + pad])
    _guide(ax, "h", 0.0)
    if fit["slope"] != 0:
        log_s_half = -fit["intercept"] / fit["slope"]
        ax.plot(xd, xd - log_s_half, color=COLORS["michaelis_menten"], linewidth=1.6, linestyle=(0, (5, 4)), zorder=2,
                alpha=0.8)
        ax.annotate("referencia n = 1\n(sin cooperatividad)", xy=(xd[-1], xd[-1] - log_s_half), xytext=(-8, -12),
                    textcoords="offset points", ha="right", va="top", color=INK_SECONDARY, fontsize=10.5)
        _guide(ax, "v", log_s_half, start=float(np.min(y)), end=0.0)
        _keypoint(ax, log_s_half, 0.0, COLORS["ts"], size=9)
        ax.annotate(f"S₀.₅ ≈ {10 ** log_s_half:.3g} {unit}", xy=(log_s_half, 0.0), xytext=(10, -6),
                    textcoords="offset points", ha="left", va="top", color=INK, fontsize=11, fontweight="semibold")
    ax.plot(xd, fit["slope"] * xd + fit["intercept"], color=color, linewidth=2.6, zorder=3)
    _markers(ax, x, y, color, size=8.5, zorder=5)
    _kbox(ax, 0.03, 0.95, f"pendiente n = {fit['slope']:.2f}\nR² = {fit['r2']:.4f}", ha="left", va="top")
    ax.set_xlim(xd[0], xd[-1])
    _finish(ax, "log₁₀[S]", "log₁₀(v / (Vmax − v))", title, subtitle)
    return fig


# ---------------------------------------------------------------------------
# 4. Visores 3D (py3Dmol)
# ---------------------------------------------------------------------------
_LIGAND_CARBON_COLORS = {"GLC": _ORANGE, "BGC": _ORANGE, "ATP": _BLUE, "ADP": _BLUE, "ANP": _BLUE, "G6P": _AQUA}
_ION_COLORS = {"MG": _GREEN, "K": _VIOLET, "K+": _VIOLET, "NA": _VIOLET, "ZN": _VIOLET}
_COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05, "CL": 1.02,
    "MG": 1.41, "NA": 1.66, "K": 2.03, "ZN": 1.22, "CA": 1.76, "FE": 1.32, "MN": 1.39,
}


def _new_view(width, height):
    import py3Dmol  # importación perezosa

    view = py3Dmol.view(width=int(width), height=int(height))
    view.setBackgroundColor("white")
    return view


def view_complex(
    pdb_text, ligand_resnames=("GLC", "ATP"), ion_resnames=("MG", "K"),
    highlight_residues=(("ASP", 201), ("LYS", 165), ("THR", 224), ("SER", 147)), label_map=None,
    cartoon_color="#c3c2b7", width=700, height=450, ligand_colors=None, show_labels=True,
):
    """Complejo enzima–ligandos: cartoon gris, ligandos en barras, iones como esferas, residuos clave rotulados.

    Los carbonos de la glucosa van en naranja y los del ATP en azul (``ligand_colors``
    permite cambiarlos); ``label_map`` traduce la numeración del archivo a la del
    cristal en los rótulos, p. ej. ``{201: "Asp205"}``.
    """
    label_map = dict(label_map or {})
    carbon_colors = dict(_LIGAND_CARBON_COLORS)
    carbon_colors.update(ligand_colors or {})
    view = _new_view(width, height)
    view.addModel(pdb_text, "pdb")
    view.setStyle({}, {"cartoon": {"color": cartoon_color}})
    for name in ligand_resnames:
        view.addStyle({"resn": name}, {"stick": {"radius": 0.25}})
        view.addStyle({"resn": name, "elem": "C"}, {"stick": {"radius": 0.25, "color": carbon_colors.get(name, _AQUA)}})
    for name in ion_resnames:
        view.addStyle({"resn": name}, {"sphere": {"radius": 0.9, "color": _ION_COLORS.get(name, _VIOLET)}})
    for resn, resi in highlight_residues:
        sel = {"resn": str(resn), "resi": int(resi)}
        view.addStyle(sel, {"stick": {"radius": 0.14}})
        if show_labels:
            text = label_map.get(int(resi)) or f"{str(resn).capitalize()}{int(resi)}"
            view.addLabel(text, {
                "fontSize": 12, "fontColor": INK, "backgroundColor": "white", "backgroundOpacity": 0.75,
                "borderColor": AXIS, "borderThickness": 0.6, "inFront": True,
            }, sel)
    if ligand_resnames:
        view.zoomTo({"resn": list(ligand_resnames)})
    else:
        view.zoomTo()
    return view


def view_qm_region(qm_pdb_text, env_pdb_text=None, link_atom_indices=(), width=700, height=450, fmt="pdb", link_color=None):
    """Región QM en bolas y barras; entorno MM (opcional) en líneas finas; átomos de enlace H resaltados.

    ``link_atom_indices`` son índices 0-based de átomos del modelo QM.
    """
    link_color = link_color or COLORS["ts"]
    view = _new_view(width, height)
    view.addModel(qm_pdb_text, fmt)
    if env_pdb_text:
        view.addModel(env_pdb_text, fmt)
        view.setStyle({"model": 1}, {"line": {"linewidth": 1.0, "opacity": 0.8}})
    view.setStyle({"model": 0}, {"stick": {"radius": 0.16}, "sphere": {"scale": 0.24}})
    if link_atom_indices:
        view.addStyle({"model": 0, "index": [int(i) for i in link_atom_indices]},
                      {"sphere": {"scale": 0.36, "color": link_color}})
    view.zoomTo({"model": 0})
    return view


def xyz_block(symbols, xyz, comment=""):
    """Bloque XYZ (texto) para ``n`` átomos: ``symbols`` (n,) y ``xyz`` (n, 3) en Å."""
    xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)
    if len(symbols) != xyz.shape[0]:
        raise ValueError("symbols y xyz deben tener el mismo número de átomos")
    lines = [str(len(symbols)), str(comment).replace("\n", " ")]
    for sym, (x, y, z) in zip(symbols, xyz):
        lines.append(f"{sym:<2s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(lines) + "\n"


def _bond_list(symbols, xyz, factor=1.2):
    xyz = np.asarray(xyz, dtype=float)
    radii = np.array([_COVALENT_RADII.get(str(s).upper(), 0.76) for s in symbols])
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    cutoff = factor * (radii[:, None] + radii[None, :])
    i, j = np.where(np.triu(d < cutoff, k=1))
    return list(zip(i.tolist(), j.tolist()))


def _pdb_frames(frames, symbols, bonds):
    conect = {}
    for i, j in bonds:
        conect.setdefault(i, []).append(j)
        conect.setdefault(j, []).append(i)
    out = []
    for f, xyz in enumerate(frames, start=1):
        out.append(f"MODEL     {f:4d}")
        for k, (sym, (x, y, z)) in enumerate(zip(symbols, xyz), start=1):
            elem = str(sym).upper()[:2]
            name = (elem + str(k))[:4]
            out.append(f"HETATM{k:5d} {name:<4s} MOL A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}")
        for i in sorted(conect):
            partners = "".join(f"{j + 1:5d}" for j in conect[i])
            out.append(f"CONECT{i + 1:5d}{partners}")
        out.append("ENDMDL")
    return "\n".join(out) + "\nEND\n"


def view_frames(frames_xyz, symbols, interval_ms=80, loop=True, bonds_from=None, width=700, height=450, style=None):
    """Animación de una lista de geometrías ``(n_atomos, 3)`` en Å (p. ej. el modo imaginario).

    Por defecto se construye un XYZ multi-cuadro (``addModelsAsFrames``) y 3Dmol
    deduce los enlaces por distancias en cada cuadro. Con ``bonds_from`` (índice
    de un cuadro de referencia, o lista de pares ``(i, j)``) los enlaces se fijan
    mediante registros CONECT para que no parpadeen durante la vibración.
    """
    frames = [np.asarray(f, dtype=float).reshape(-1, 3) for f in frames_xyz]
    if not frames:
        raise ValueError("frames_xyz está vacío")
    view = _new_view(width, height)
    if bonds_from is None:
        text = "".join(xyz_block(symbols, f, f"cuadro {i}") for i, f in enumerate(frames))
        view.addModelsAsFrames(text, "xyz")
    else:
        if isinstance(bonds_from, (int, np.integer)):
            bonds = _bond_list(symbols, frames[int(bonds_from)])
        else:
            bonds = [(int(i), int(j)) for i, j in bonds_from]
        view.addModelsAsFrames(_pdb_frames(frames, symbols, bonds), "pdb")
    view.setStyle({}, style or {"stick": {"radius": 0.15}, "sphere": {"scale": 0.22}})
    view.zoomTo()
    view.animate({"loop": "backAndForth" if loop else "forward", "interval": int(interval_ms), "reps": 0 if loop else 1})
    return view


def mode_animation_frames(xyz, displacement, n_frames=24, amplitude=0.5):
    """Cuadros de una vibración: ``xyz + amplitude·sin(φ)·desplazamiento`` con ``φ`` en ``[0, 2π)``.

    El desplazamiento se normaliza para que ``amplitude`` sea la excursión
    máxima de un átomo en Å. Devuelve una lista de arreglos ``(n_atomos, 3)``.
    """
    xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)
    disp = np.asarray(displacement, dtype=float).reshape(xyz.shape)
    norm = float(np.max(np.linalg.norm(disp, axis=1)))
    if norm > 0:
        disp = disp / norm
    phases = np.linspace(0.0, 2.0 * np.pi, int(n_frames), endpoint=False)
    return [xyz + float(amplitude) * np.sin(phi) * disp for phi in phases]


# ---------------------------------------------------------------------------
# Resultados como tarjetas (HTML en la salida de la celda; Colab no lo sanea)
# ---------------------------------------------------------------------------
_ROLE_COLOR = {"azul": _BLUE, "naranja": _ORANGE, "agua": _AQUA, "amarillo": _YELLOW, "magenta": _MAGENTA,
               "verde": _GREEN, "violeta": _VIOLET, "rojo": _RED, "gris": INK_MUTED}
_CARD_FONT = "Figtree, 'Avenir Next', 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"


def _esc(text):
    import html

    return html.escape(str(text))


def tarjetas(items, titulo=None, nota=None, columnas=None):
    """Muestra resultados como tarjetas: ``items`` = lista de dicts o tuplas
    ``(etiqueta, valor, unidad, detalle, color)`` (los tres últimos opcionales; ``color`` es un nombre de
    ``_ROLE_COLOR`` o un hex). Devuelve un objeto ``IPython.display.HTML`` (última expresión de la celda).

    Sustituye a los ``print`` de resultados: el número grande, la unidad pequeña y una línea que lo traduce.
    """
    from IPython.display import HTML

    cards = []
    for it in items:
        if not isinstance(it, dict):
            it = dict(zip(("etiqueta", "valor", "unidad", "detalle", "color"), it))
        color = _ROLE_COLOR.get(it.get("color") or "azul", it.get("color") or _BLUE)
        unidad = f'<span style="font-size:15px;color:{INK_SECONDARY};font-weight:500;margin-left:4px">{_esc(it.get("unidad") or "")}</span>'
        detalle = (f'<div style="font-size:13px;color:{INK_SECONDARY};margin-top:6px;line-height:1.35">{_esc(it["detalle"])}</div>'
                   if it.get("detalle") else "")
        cards.append(
            f'<div style="background:#ffffff;border:1px solid {GRID};border-radius:14px;padding:14px 16px 14px 18px;'
            f'position:relative;overflow:hidden;min-width:0">'
            f'<div style="position:absolute;left:0;top:0;bottom:0;width:5px;background:{color}"></div>'
            f'<div style="font-size:13.5px;letter-spacing:.01em;color:{INK_MUTED};font-weight:600">{_esc(it.get("etiqueta", ""))}</div>'
            f'<div style="font-size:28px;font-weight:650;color:{INK};margin-top:4px;line-height:1.1">{_esc(it.get("valor", ""))}{unidad}</div>'
            f"{detalle}</div>")
    n = columnas or min(len(cards), 4)
    cabecera = (f'<div style="font-size:16px;font-weight:650;color:{INK};margin:0 0 10px 2px">{_esc(titulo)}</div>' if titulo else "")
    pie = (f'<div style="font-size:13.5px;color:{INK_SECONDARY};margin:12px 2px 0;line-height:1.45">{_esc(nota)}</div>' if nota else "")
    return HTML(
        f'<div style="font-family:{_CARD_FONT};background:{SURFACE};border:1px solid {GRID};border-radius:18px;padding:16px;max-width:960px">'
        f"{cabecera}<div style=\"display:grid;grid-template-columns:repeat({n},minmax(0,1fr));gap:12px\">{''.join(cards)}</div>{pie}</div>")


def mensaje(texto, tipo="idea", titulo=None):
    """Recuadro de mensaje para la salida de una celda. ``tipo``: idea, ok, ojo, dato."""
    from IPython.display import HTML

    estilos = {"idea": (_BLUE, "💡"), "ok": (STATUS["good"], "✅"), "ojo": (STATUS["warning"], "⚠️"), "dato": (_AQUA, "🔬")}
    color, icono = estilos.get(tipo, estilos["idea"])
    cab = f'<div style="font-weight:650;color:{INK};margin-bottom:4px">{icono} {_esc(titulo)}</div>' if titulo else ""
    cuerpo = _esc(texto) if titulo else f"{icono} {_esc(texto)}"
    return HTML(
        f'<div style="font-family:{_CARD_FONT};font-size:14.5px;line-height:1.5;color:{INK_SECONDARY};background:#ffffff;'
        f'border:1px solid {GRID};border-left:5px solid {color};border-radius:12px;padding:12px 16px;max-width:960px;margin:6px 0">'
        f"{cab}{cuerpo}</div>")


def tabla(df, titulo=None, nota=None, formatos=None):
    """Tabla con el estilo del curso (``pandas.DataFrame`` → HTML); sustituye a ``display(df)``.

    ``formatos`` = {columna: "{:.1f}"}; los NaN se muestran como «—».
    """
    from IPython.display import HTML

    df = df.copy()
    for col in df.columns:
        fmt = (formatos or {}).get(col)
        df[col] = [("—" if (isinstance(v, float) and not np.isfinite(v)) or v is None or v == ""
                    else (fmt.format(v) if fmt else v)) for v in df[col]]
    th = "".join(f'<th style="text-align:left;padding:8px 12px;color:{INK_MUTED};font-weight:600;font-size:12.5px;'
                 f'letter-spacing:.01em;border-bottom:1.5px solid {AXIS}">{_esc(c)}</th>' for c in df.columns)
    filas = []
    for i, row in enumerate(df.itertuples(index=False)):
        fondo = "#ffffff" if i % 2 == 0 else SURFACE
        tds = "".join(f'<td style="padding:8px 12px;color:{INK if j == 0 else INK_SECONDARY};font-weight:{600 if j == 0 else 400};'
                      f'border-bottom:1px solid {GRID};vertical-align:top">{_esc(v)}</td>' for j, v in enumerate(row))
        filas.append(f'<tr style="background:{fondo}">{tds}</tr>')
    cab = f'<div style="font-size:16px;font-weight:650;color:{INK};margin:0 0 10px 2px">{_esc(titulo)}</div>' if titulo else ""
    pie = f'<div style="font-size:13px;color:{INK_MUTED};margin:10px 2px 0;line-height:1.4">{_esc(nota)}</div>' if nota else ""
    return HTML(
        f'<div style="font-family:{_CARD_FONT};font-size:14px;background:{SURFACE};border:1px solid {GRID};border-radius:18px;'
        f'padding:16px;max-width:960px;overflow-x:auto">{cab}<table style="border-collapse:collapse;width:100%">'
        f"<thead><tr>{th}</tr></thead><tbody>{''.join(filas)}</tbody></table>{pie}</div>")


def mostrar(*objetos):
    """Muestra en orden figuras (y las cierra, para que no se dupliquen al final de la celda) y HTML.

    Con código oculto el orden importa: primero la gráfica, después la tarjeta que la resume.
    """
    from IPython.display import display

    for obj in objetos:
        if obj is None:
            continue
        display(obj)
        if isinstance(obj, Figure):
            plt.close(obj)
