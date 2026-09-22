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
    "plot_md_timeseries", "plot_md_summary",
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


def apply_style():
    """Activa el estilo compartido del curso (``matplotlib.style.use``). Devuelve la ruta usada."""
    global _STYLE_APPLIED
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
    ax.annotate(
        "", xy=(x, y1), xytext=(x, y0),
        arrowprops=dict(arrowstyle="<->", color=INK_SECONDARY, linewidth=1.0, shrinkA=0, shrinkB=0),
        zorder=5,
    )
    ax.annotate(
        text, xy=(x, 0.5 * (y0 + y1)), xytext=(x_text_offset, 0), textcoords="offset points",
        ha="left", va="center", color=INK, fontsize=10,
    )


def plot_energy_profile(
    x, energy, xlabel="Coordenada de reacción", ax=None, ts_index=None, labels=None, title=None,
    subtitle=None, annotate_barrier=True, relative=True, smooth=True, color=None, label=None,
    ylabel="Energía relativa (kcal/mol)",
):
    """Perfil de energía a lo largo de la coordenada de reacción (kcal/mol).

    ``energy`` se refiere al primer punto si ``relative=True``. El estado de
    transición (``ts_index``, por defecto el máximo) se resalta en naranja y la
    barrera ``ΔE‡`` se anota con una flecha de cota vertical desde el nivel de
    reactivos. ``labels`` es una lista opcional de textos por punto (``None``
    para omitir). Devuelve la figura.
    """
    x = np.asarray(x, dtype=float)
    e = np.asarray(energy, dtype=float)
    order = np.argsort(x, kind="stable")
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

    fig, ax = _get_ax(ax)
    xs, ys = _smooth_curve(x, e) if smooth else (x, e)
    ax.plot(xs, ys, color=color, linewidth=2.0, label=label, zorder=2)
    _markers(ax, x, e, color)
    _markers(ax, [x[ts_index]], [e[ts_index]], COLORS["ts"], size=8.5, zorder=5)

    if labels is not None:
        for xi, ei, text in zip(x, e, labels):
            if text:
                ax.annotate(text, xy=(xi, ei), xytext=(0, 9), textcoords="offset points",
                            ha="center", va="bottom", color=INK_SECONDARY, fontsize=9)

    if annotate_barrier:
        e_ref = e[0]
        x_ts = x[ts_index]
        span = float(np.ptp(x)) if np.ptp(x) > 0 else 1.0
        x_arrow = x_ts + 0.04 * span
        _guide(ax, "h", e_ref, start=x[0], end=x_arrow)
        _guide(ax, "h", e[ts_index], start=x_ts, end=x_arrow, linestyle="-", linewidth=0.8, color=AXIS)
        _dimension_arrow(ax, x_arrow, e_ref, e[ts_index], f"ΔE‡ = {_num(e[ts_index] - e_ref)} kcal/mol")
        ax.set_xlim(x[0] - 0.03 * span, max(x[-1] + 0.03 * span, x_arrow + 0.30 * span))

    ax.margins(y=0.12)
    _finish(ax, xlabel, ylabel, title, subtitle)
    if label:
        _legend(ax)
    return fig


def plot_energy_profiles(
    profiles, xlabel="Coordenada de reacción", ax=None, title=None, subtitle=None, relative=True,
    smooth=True, mark_ts=True, ylabel="Energía relativa (kcal/mol)", colors=None,
):
    """Varios perfiles ``(x, energía, etiqueta)`` superpuestos con leyenda (colores fijos por orden)."""
    fig, ax = _get_ax(ax)
    colors = list(colors) if colors else PALETTE
    if len(profiles) > len(colors):
        raise ValueError("máximo 8 perfiles: agrupa o usa múltiplos pequeños")
    for (x, e, label), color in zip(profiles, colors):
        x = np.asarray(x, dtype=float)
        e = np.asarray(e, dtype=float)
        order = np.argsort(x)
        x, e = x[order], e[order]
        if relative:
            e = e - e[0]
        xs, ys = _smooth_curve(x, e) if smooth else (x, e)
        ax.plot(xs, ys, color=color, linewidth=2.0, label=label)
        _markers(ax, x, e, color)
        if mark_ts:
            i = int(np.argmax(e))
            _markers(ax, [x[i]], [e[i]], color, size=8.5, hollow=True, zorder=5)
    _finish(ax, xlabel, ylabel, title, subtitle)
    _legend(ax)
    return fig


def plot_energy_levels(
    levels, ax=None, title=None, subtitle=None, connect=True, compare=None, label="con enzima",
    compare_label="sin enzima", annotate_barrier=True, show_values=True, ts_indices=None,
    ylabel="Energía (kcal/mol)", width=0.56,
):
    """Diagrama de niveles de energía por etapas.

    ``levels`` = lista de ``(etiqueta, energía_kcal)``, p. ej.
    ``[("E + S", 0), ("ES", -3), ("TS", 15), ("EP", -5), ("E + P", -8)]``.
    Los estados de transición (etiqueta con "TS" o "‡", o ``ts_indices``) se
    dibujan en naranja y la barrera se anota desde el mínimo inmediatamente
    anterior. ``compare`` es una segunda lista (p. ej. sin enzima) dibujada en
    gris; si tiene la misma longitud comparte las posiciones.
    """
    names = [str(name) for name, _ in levels]
    e = np.asarray([float(val) for _, val in levels])
    n = len(names)
    if ts_indices is None:
        ts_indices = [i for i, name in enumerate(names) if "TS" in name.upper() or "‡" in name]
    ts_set = set(int(i) for i in ts_indices)
    fig, ax = _get_ax(ax)
    half = width / 2.0

    def draw(values, colors, zorder, connector_color):
        for i, (val, col) in enumerate(zip(values, colors)):
            ax.hlines(val, i - half, i + half, color=col, linewidth=2.6, zorder=zorder)
        if connect:
            for i in range(len(values) - 1):
                ax.plot([i + half, i + 1 - half], [values[i], values[i + 1]], color=connector_color,
                        linewidth=1.0, linestyle=(0, (3, 3)), zorder=zorder - 1)

    if compare is not None:
        ce = np.asarray([float(val) for _, val in compare])
        draw(ce, [COLORS["sin_enzima"]] * len(ce), 2, AXIS)
    draw(e, [COLORS["ts"] if i in ts_set else COLORS["reactivo"] for i in range(n)], 3, INK_MUTED)

    if show_values:
        for i, val in enumerate(e):
            ax.annotate(_num(val), xy=(i, val), xytext=(0, 5), textcoords="offset points",
                        ha="center", va="bottom", color=INK_MUTED, fontsize=8.5)

    if annotate_barrier and ts_set:
        for i_ts in sorted(ts_set):
            if i_ts == 0:
                continue
            i_ref = int(np.argmin(e[:i_ts]))
            e_ref = e[i_ref]
            x_arrow = i_ts + half + 0.12
            _guide(ax, "h", e_ref, start=i_ref + half, end=x_arrow)
            _guide(ax, "h", e[i_ts], start=i_ts + half, end=x_arrow, linestyle="-", linewidth=0.8, color=AXIS)
            _dimension_arrow(ax, x_arrow, e_ref, e[i_ts], f"ΔE‡ = {_num(e[i_ts] - e_ref)} kcal/mol", x_text_offset=5)

    ax.set_xticks(range(n), names)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, n - 0.5 + (1.6 if (annotate_barrier and ts_set) else 0.0))
    ax.margins(y=0.15)
    _finish(ax, None, ylabel, title, subtitle)
    ax.tick_params(axis="x", length=0, labelcolor=INK_SECONDARY)
    if compare is not None:
        handles = [Line2D([], [], color=COLORS["reactivo"], linewidth=2.6, label=label),
                   Line2D([], [], color=COLORS["sin_enzima"], linewidth=2.6, label=compare_label)]
        _legend(ax, handles=handles, loc="upper right")
    return fig


# ---------------------------------------------------------------------------
# 2. Dinámica molecular
# ---------------------------------------------------------------------------
_MD_LABELS = {
    "temperatura_K": "Temperatura (K)",
    "rmsd_CA_A": "RMSD Cα (Å)",
    "d_PG_O6_A": "d(PG–O6) (Å)",
    "d_PG_O3B_A": "d(PG–O3B) (Å)",
    "d_O6_OD1asp205_A": "d(O6–OD1 Asp205) (Å)",
    "d_Mg_PG_A": "d(Mg–PG) (Å)",
    "d_NZlys169_O6_A": "d(NZ Lys169–O6) (Å)",
}


def plot_md_timeseries(
    df, columns, labels=None, ylabel="", ax=None, time_column="tiempo_ps", xlabel="Tiempo (ps)",
    title=None, subtitle=None, colors=None,
):
    """Series temporales de una dinámica molecular a partir de un ``DataFrame`` con ``tiempo_ps``.

    ``columns`` es una lista de nombres de columna; ``labels`` (opcional) sus
    etiquetas en español (por defecto se usan nombres conocidos o la columna).
    """
    if isinstance(columns, str):
        columns = [columns]
    columns = list(columns)
    if labels is None:
        labels = [_MD_LABELS.get(c, c) for c in columns]
    colors = list(colors) if colors else PALETTE
    if len(columns) > len(colors):
        raise ValueError("máximo 8 series por panel: usa múltiplos pequeños")
    fig, ax = _get_ax(ax)
    t = np.asarray(df[time_column], dtype=float)
    for col, lab, color in zip(columns, labels, colors):
        ax.plot(t, np.asarray(df[col], dtype=float), color=color, linewidth=1.6, label=lab)
    ax.set_xlim(t.min(), t.max())
    _finish(ax, xlabel, ylabel, title, subtitle)
    if len(columns) > 1:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax + 0.16 * (ymax - ymin))
        _legend(ax, loc="upper left", ncol=min(len(columns), 3))
    return fig


def plot_md_summary(
    df, crystal_value=None, title="Resumen de la dinámica molecular", subtitle=None,
    time_column="tiempo_ps", crystal_label="cristal", bins=30,
):
    """Múltiplos pequeños (2×2): RMSD, d(PG–O6), d(O6–OD1 Asp205) y el histograma de d(PG–O6).

    ``crystal_value`` dibuja una línea vertical de referencia (distancia en la
    estructura cristalina) sobre el histograma.
    """
    _ensure_style()
    fig, axes = figure(10.0, 6.4, nrows=2, ncols=2)
    axes = np.asarray(axes).ravel()
    t = np.asarray(df[time_column], dtype=float)
    panels = [
        ("rmsd_CA_A", "RMSD Cα", "RMSD (Å)"),
        ("d_PG_O6_A", "Distancia PG–O6 (ataque nucleófilo)", "d (Å)"),
        ("d_O6_OD1asp205_A", "Distancia O6–OD1 Asp205 (base catalítica)", "d (Å)"),
    ]
    for ax, (col, panel_title, ylab) in zip(axes[:3], panels):
        if col in df.columns:
            ax.plot(t, np.asarray(df[col], dtype=float), color=COLORS["datos"], linewidth=1.4)
            ax.set_xlim(t.min(), t.max())
        else:
            ax.text(0.5, 0.5, f"sin columna «{col}»", transform=ax.transAxes, ha="center", va="center",
                    color=INK_MUTED)
        _finish(ax, "Tiempo (ps)", ylab)
        ax.set_title(panel_title, loc="left", fontsize=11, color=INK_SECONDARY, fontweight="normal", pad=8)

    ax = axes[3]
    if "d_PG_O6_A" in df.columns:
        d = np.asarray(df["d_PG_O6_A"], dtype=float)
        ax.hist(d, bins=bins, color=COLORS["datos"], edgecolor=SURFACE, linewidth=1.0, zorder=2)
        ax.axvline(float(np.mean(d)), color=INK_SECONDARY, linewidth=1.0, linestyle=(0, (4, 3)), zorder=3)
        ax.annotate(f"media {np.mean(d):.2f} Å", xy=(float(np.mean(d)), 1.0), xycoords=("data", "axes fraction"),
                    xytext=(4, -2), textcoords="offset points", ha="left", va="top", color=INK_SECONDARY, fontsize=9)
        if crystal_value is not None:
            ax.axvline(float(crystal_value), color=INK, linewidth=1.2, zorder=3)
            # la línea suele quedar en el borde izquierdo: la etiqueta va a su derecha, un poco más abajo que la media
            ax.annotate(f"{crystal_label} {float(crystal_value):.2f} Å", xy=(float(crystal_value), 0.88),
                        xycoords=("data", "axes fraction"), xytext=(4, 0), textcoords="offset points",
                        ha="left", va="top", color=INK, fontsize=9)
    _finish(ax, "d(PG–O6) (Å)", "Frecuencia")
    ax.set_title("Distribución de d(PG–O6)", loc="left", fontsize=11, color=INK_SECONDARY, fontweight="normal", pad=8)
    _fig_title(fig, title, subtitle)
    return fig


# ---------------------------------------------------------------------------
# 3. Cinética
# ---------------------------------------------------------------------------
def plot_michaelis_menten(
    s, v, fit=None, ax=None, show_km=True, title=None, subtitle=None, xlabel="[S] (mM)",
    ylabel="v₀ (µM/s)", show_params=True, point_label=None, curve_label=None, color=None,
):
    """Datos ``(s, v)`` y curva de Michaelis–Menten ajustada.

    ``fit`` es el diccionario de :func:`enzimas.kinetics.fit_michaelis_menten`
    (o ``None`` para dibujar solo los puntos). Con ``show_km`` se añaden guías
    discontinuas para ``Vmax``, ``Vmax/2`` y ``Km``.
    """
    s = np.asarray(s, dtype=float)
    v = np.asarray(v, dtype=float)
    color = color or COLORS["datos"]
    fig, ax = _get_ax(ax)
    _markers(ax, s, v, color, label=point_label)
    if fit is not None:
        vmax, km = float(fit["vmax"]), float(fit["km"])
        s_dense = np.linspace(0.0, float(np.max(s)) * 1.05, 300)
        ax.plot(s_dense, _kin.michaelis_menten(s_dense, vmax, km), color=color, linewidth=2.0, label=curve_label, zorder=2)
        ax.set_ylim(0, max(vmax * 1.18, float(np.max(v)) * 1.05))
        ax.set_xlim(0, s_dense[-1])
        if show_km:
            _guide(ax, "h", vmax)
            ax.annotate("Vmax", xy=(1.0, vmax), xycoords=("axes fraction", "data"), xytext=(0, 3),
                        textcoords="offset points", ha="right", va="bottom", color=INK_MUTED, fontsize=9.5)
            _guide(ax, "h", vmax / 2, start=0, end=km)
            if km / s_dense[-1] > 0.12:
                ax.annotate("Vmax/2", xy=(km / 2, vmax / 2), xytext=(0, 3), textcoords="offset points",
                            ha="center", va="bottom", color=INK_MUTED, fontsize=9.5)
            else:  # Km pequeño: el hueco bajo la curva a la derecha de Km siempre está libre
                ax.annotate("Vmax/2", xy=(km, vmax / 2), xytext=(8, -3), textcoords="offset points",
                            ha="left", va="top", color=INK_MUTED, fontsize=9.5)
            _guide(ax, "v", km, start=0, end=vmax / 2)
            ax.annotate("Km", xy=(km, 0), xytext=(4, 3), textcoords="offset points",
                        ha="left", va="bottom", color=INK_MUTED, fontsize=9.5)
        if show_params:
            text = f"Vmax = {_value_pm(vmax, fit.get('vmax_err'))}\nKm = {_value_pm(km, fit.get('km_err'))}"
            if "r2" in fit:
                text += f"\nR² = {fit['r2']:.4f}"
            ax.text(0.98, 0.06, text, transform=ax.transAxes, ha="right", va="bottom", color=INK_SECONDARY, fontsize=9.5)
    _finish(ax, xlabel, ylabel, title, subtitle)
    if point_label or curve_label:
        _legend(ax, loc="center right")
    return fig


def plot_hill_vs_mm(
    s, v_mm, v_hill, n_hill=None, ax=None, title=None, subtitle=None, xlabel="[S] (mM)",
    ylabel="v₀ (µM/s)", s_half=None, labels=None,
):
    """Superpone la hipérbola de Michaelis–Menten (n = 1) y la sigmoide de Hill (n > 1)."""
    s = np.asarray(s, dtype=float)
    if labels is None:
        n_text = f"n = {float(n_hill):.2g}" if n_hill is not None else "n > 1"
        labels = ("Michaelis–Menten (n = 1)", f"Hill ({n_text})")
    fig, ax = _get_ax(ax)
    ax.plot(s, np.asarray(v_mm, dtype=float), color=COLORS["michaelis_menten"], linewidth=2.0, label=labels[0])
    ax.plot(s, np.asarray(v_hill, dtype=float), color=COLORS["hill"], linewidth=2.0, label=labels[1])
    if s_half is not None:
        _guide(ax, "v", float(s_half))
        ax.annotate("S₀.₅", xy=(float(s_half), 0.0), xytext=(4, 3), textcoords="offset points",
                    ha="left", va="bottom", color=INK_MUTED, fontsize=9.5)
    ax.set_xlim(0, s.max())
    ax.set_ylim(bottom=0)
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
    (trazo discontinuo) hasta los cortes con los ejes, que se marcan y rotulan.
    """
    s = np.asarray(s, dtype=float)
    v = np.asarray(v, dtype=float)
    fits = dict(fits or {})
    color = color or COLORS["datos"]
    fig, axes = figure(12.0, 3.9, nrows=1, ncols=3)
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
        ax.axhline(0, color=AXIS, linewidth=0.8, zorder=1)
        ax.axvline(0, color=AXIS, linewidth=0.8, zorder=1)
        xd = np.array([lo, hi])
        ax.plot(xd, slope * xd + intercept, color=color, linewidth=1.2, linestyle=(0, (4, 3)), zorder=2)
        xr = np.array([float(np.min(x)), float(np.max(x))])
        ax.plot(xr, slope * xr + intercept, color=color, linewidth=2.0, zorder=3)
        _markers(ax, x, y, color)
        _markers(ax, [0.0], [intercept], INK_SECONDARY, hollow=True, size=7, zorder=5)
        ax.annotate(y_int_label, xy=(0.0, intercept), xytext=(-6, 5), textcoords="offset points",
                    ha="right", va="bottom", color=INK_SECONDARY, fontsize=9.5)
        if np.isfinite(x_int):
            _markers(ax, [x_int], [0.0], INK_SECONDARY, hollow=True, size=7, zorder=5)
            ax.annotate(x_int_label, xy=(x_int, 0.0), xytext=(0, -8), textcoords="offset points",
                        ha="center", va="top", color=INK_SECONDARY, fontsize=9.5)
        ax.set_xlim(lo, hi)
        ax.margins(y=0.15)
        _finish(ax, xl, yl)
        ax.set_title(f"{name}  (R² = {fit['r2']:.3f})", loc="left", fontsize=11, color=INK_SECONDARY, fontweight="normal", pad=8)
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
    rampa azul secuencial (claro = poco inhibidor, oscuro = mucho).
    """
    s = np.asarray(s, dtype=float)
    curves = sorted(((float(i), np.asarray(v, dtype=float)) for i, v in curves), key=lambda c: c[0])
    colors = sequential_blue(len(curves))
    fig, ax = _get_ax(ax)
    for (i, v), color in zip(curves, colors):
        ax.plot(s, v, color=color, linewidth=2.0, label=_inhibitor_label(i, unit))
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
    """Lineweaver–Burk para una familia de inhibición: muestra el patrón de intersección.

    Cada recta se ajusta por mínimos cuadrados y se extrapola (trazo
    discontinuo) hasta más allá de su corte con el eje ``1/[S]``. El subtítulo
    por defecto describe el patrón característico de ``kind``.
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
    ax.axhline(0, color=AXIS, linewidth=0.8, zorder=1)
    ax.axvline(0, color=AXIS, linewidth=0.8, zorder=1)
    for i, x, y, fit, color in lines:
        xd = np.array([x_min, x_max * 1.05])
        ax.plot(xd, fit["slope"] * xd + fit["intercept"], color=color, linewidth=1.2, linestyle=(0, (4, 3)), zorder=2)
        xr = np.array([float(np.min(x)), float(np.max(x))])
        ax.plot(xr, fit["slope"] * xr + fit["intercept"], color=color, linewidth=2.0, zorder=3, label=_inhibitor_label(i, unit))
        _markers(ax, x, y, color)
    ax.set_xlim(x_min, x_max * 1.05)
    ax.set_ylim(top=y_max * 1.15)
    kind = str(kind).lower()
    if title is None:
        title = f"{_INHIBITION_NAMES.get(kind, 'Inhibición')}: Lineweaver–Burk"
    if subtitle is None:
        subtitle = _INHIBITION_LB_PATTERN.get(kind)
    _finish(ax, "1/[S]", "1/v₀", title, subtitle)
    _legend(ax, loc="upper left")
    return fig


def plot_ode_simulation(
    sim, title=r"Simulación del mecanismo E + S $\rightleftharpoons$ ES $\rightarrow$ E + P", subtitle=None, t_zoom=None,
    xlabel="Tiempo (s)", ylabel="Concentración (µM)",
):
    """Concentraciones E, S, ES y P frente al tiempo (diccionario de :func:`kinetics.simulate_mechanism`).

    Dos paneles con un solo eje cada uno: a la izquierda S y P (escala grande);
    a la derecha E y ES (escala pequeña), con ``t_zoom`` para ampliar el estado
    pre-estacionario (límite superior del tiempo en el panel derecho).
    """
    t = np.asarray(sim["t"], dtype=float)
    fig, (ax1, ax2) = figure(11.0, 4.0, nrows=1, ncols=2)
    ax1.plot(t, sim["S"], color=COLORS["sustrato"], linewidth=2.0, label="S (sustrato)")
    ax1.plot(t, sim["P"], color=COLORS["producto_p"], linewidth=2.0, label="P (producto)")
    ax1.set_xlim(t.min(), t.max())
    ax1.set_ylim(bottom=0)
    _finish(ax1, xlabel, ylabel)
    ax1.set_title("Sustrato y producto", loc="left", fontsize=11, color=INK_SECONDARY, fontweight="normal", pad=8)
    _legend(ax1, loc="center right")

    ax2.plot(t, sim["E"], color=COLORS["enzima"], linewidth=2.0, label="E (enzima libre)")
    ax2.plot(t, sim["ES"], color=COLORS["complejo_es"], linewidth=2.0, label="ES (complejo)")
    if t_zoom is not None:
        ax2.set_xlim(0, float(t_zoom))
        ax2.set_title("Enzima y complejo ES (estado pre-estacionario)", loc="left", fontsize=11,
                      color=INK_SECONDARY, fontweight="normal", pad=8)
    else:
        ax2.set_xlim(t.min(), t.max())
        ax2.set_title("Enzima libre y complejo ES", loc="left", fontsize=11, color=INK_SECONDARY, fontweight="normal", pad=8)
    ax2.set_ylim(bottom=0)
    _finish(ax2, xlabel, ylabel)
    _legend(ax2, loc="center right")
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


def plot_arrhenius(temperatures_k, rates, fit=None, ax=None, title="Gráfico de Arrhenius", subtitle=None, color=None):
    """``ln k`` frente a ``1000/T``; ``fit`` es el diccionario de :func:`kinetics.fit_arrhenius`."""
    t = np.asarray(temperatures_k, dtype=float)
    k = np.asarray(rates, dtype=float)
    color = color or COLORS["datos"]
    x = 1000.0 / t
    fig, ax = _get_ax(ax)
    if fit is not None:
        xd = np.linspace(x.min(), x.max(), 50)
        yd = np.log(float(fit["a"])) - float(fit["ea_kcal"]) / _kin.R_KCAL * (xd / 1000.0)
        ax.plot(xd, yd, color=color, linewidth=2.0, zorder=2)
        ax.text(0.98, 0.94, f"Ea = {_num(fit['ea_kcal'])} kcal/mol\nR² = {fit['r2']:.4f}", transform=ax.transAxes,
                ha="right", va="top", color=INK_SECONDARY, fontsize=9.5)
    _markers(ax, x, np.log(k), color)
    _finish(ax, "1000/T (K⁻¹)", "ln k", title, subtitle)
    return fig


def plot_eyring(temperatures_k, rates, fit=None, ax=None, title="Gráfico de Eyring", subtitle=None, color=None):
    """``ln(k/T)`` frente a ``1/T``; ``fit`` es el diccionario de :func:`kinetics.fit_eyring`."""
    t = np.asarray(temperatures_k, dtype=float)
    k = np.asarray(rates, dtype=float)
    color = color or COLORS["datos"]
    x = 1.0 / t
    fig, ax = _get_ax(ax)
    if fit is not None:
        xd = np.linspace(x.min(), x.max(), 50)
        intercept = np.log(_kin.KB_J / _kin.H_J) + float(fit["delta_s_cal"]) / (_kin.R_KCAL * 1000.0)
        slope = -float(fit["delta_h_kcal"]) / _kin.R_KCAL
        ax.plot(xd, slope * xd + intercept, color=color, linewidth=2.0, zorder=2)
        text = (f"ΔH‡ = {_num(fit['delta_h_kcal'])} kcal/mol\nΔS‡ = {_num(fit['delta_s_cal'])} cal/(mol·K)\n"
                f"ΔG‡(298 K) = {_num(fit['delta_g_kcal'])} kcal/mol")
        ax.text(0.98, 0.94, text, transform=ax.transAxes, ha="right", va="top", color=INK_SECONDARY, fontsize=9.5)
    _markers(ax, x, np.log(k / t), color)
    ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, -3), useMathText=True)
    ax.xaxis.get_offset_text().set_color(INK_MUTED)
    _finish(ax, "1/T (K⁻¹)", "ln(k/T)", title, subtitle)
    return fig


def plot_ph_profile(
    ph, v, ax=None, title="Perfil de pH", subtitle=None, pkas=None, xlabel="pH", ylabel="v₀ (µM/s)",
    color=None, fit_curve=None,
):
    """Actividad frente al pH (puntos; curva continua si se pasa ``fit_curve=(ph_denso, v_denso)``).

    ``pkas=(pKa₁, pKa₂)`` añade guías verticales rotuladas.
    """
    ph = np.asarray(ph, dtype=float)
    v = np.asarray(v, dtype=float)
    color = color or COLORS["datos"]
    fig, ax = _get_ax(ax)
    if fit_curve is not None:
        ax.plot(np.asarray(fit_curve[0], dtype=float), np.asarray(fit_curve[1], dtype=float), color=color, linewidth=2.0, zorder=2)
    elif ph.size > 40:
        ax.plot(ph, v, color=color, linewidth=2.0, zorder=2)
    if ph.size <= 40:
        _markers(ax, ph, v, color)
    if pkas:
        for idx, pka in enumerate(pkas, start=1):
            _guide(ax, "v", float(pka))
            ax.annotate(f"pKa{'₁₂₃'[idx - 1]} = {float(pka):.1f}", xy=(float(pka), 1.0), xycoords=("data", "axes fraction"),
                        xytext=(4, -2), textcoords="offset points", ha="left", va="top", color=INK_MUTED, fontsize=9.5)
    ax.set_ylim(bottom=0)
    _finish(ax, xlabel, ylabel, title, subtitle)
    return fig


def plot_hill_plot(s, v, vmax, ax=None, title="Gráfico de Hill", subtitle=None, color=None, unit="mM"):
    """Gráfico de Hill ``log₁₀(v/(Vmax − v))`` frente a ``log₁₀[S]`` con la pendiente ``n`` anotada."""
    x, y = _kin.hill_plot(s, v, float(vmax))
    fit = _kin.linear_fit(x, y)
    color = color or COLORS["datos"]
    fig, ax = _get_ax(ax)
    xd = np.array([x.min(), x.max()])
    ax.plot(xd, fit["slope"] * xd + fit["intercept"], color=color, linewidth=2.0, zorder=2)
    _markers(ax, x, y, color)
    _guide(ax, "h", 0.0)
    if fit["slope"] != 0:
        log_s_half = -fit["intercept"] / fit["slope"]
        _guide(ax, "v", log_s_half, start=float(np.min(y)), end=0.0)
        ax.annotate(f"S₀.₅ ≈ {10 ** log_s_half:.3g} {unit}", xy=(log_s_half, 0.0), xytext=(5, -4),
                    textcoords="offset points", ha="left", va="top", color=INK_MUTED, fontsize=9.5)
    ax.text(0.03, 0.94, f"pendiente n = {fit['slope']:.2f}\nR² = {fit['r2']:.4f}", transform=ax.transAxes,
            ha="left", va="top", color=INK_SECONDARY, fontsize=9.5)
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
