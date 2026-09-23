"""Mini-biblioteca SVG para las ilustraciones del cuaderno.

Todas las figuras comparten el sistema de diseño de ``enzimas/viz.py`` (superficie, tintas y
colores por entidad) para que dibujos y gráficas se lean como un solo sistema:

* reactivo / enzima = azul, estado de transición / sustrato = naranja, producto / ES = aguamarina,
  producto P = magenta, Mg²⁺ = verde, K⁺ = violeta, rojo = alerta/inhibidor.
* El texto va siempre en tinta (primaria, secundaria o atenuada), nunca en el color de la serie;
  una marca de color al lado aporta la identidad.
* Un solo grosor de trazo (2.5), esquinas redondeadas, mucho aire.

Cada figura es un :class:`Lienzo` de 1000 unidades de ancho. ``guardar`` escribe el SVG y lo
rasteriza a PNG (×1.6 → 1600 px) con ``rsvg-convert`` usando las fuentes del sistema, de modo que
en Colab la imagen se ve igual en cualquier navegador.
"""
from __future__ import annotations

import math
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parents[2]
SALIDA = ROOT / "assets" / "ilustraciones"

# ---------------------------------------------------------------- sistema de diseño (= viz.py)
SUPERFICIE = "#fcfcfb"
TINTA = "#0b0b0b"
TINTA2 = "#52514e"
TINTA3 = "#898781"
REJILLA = "#e1e0d9"
EJE = "#c3c2b7"

AZUL, NARANJA, AGUA, AMARILLO, MAGENTA, VERDE, VIOLETA, ROJO = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")

C = {
    "reactivo": AZUL, "ts": NARANJA, "producto": AGUA,
    "enzima": AZUL, "sustrato": NARANJA, "es": AGUA, "p": MAGENTA,
    "mg": VERDE, "k": VIOLETA, "inhibidor": ROJO, "neutro": TINTA3,
}

SANS = "'Avenir Next', 'Avenir', 'Helvetica Neue', Arial, sans-serif"
SERIF = "'STIX Two Text', 'STIXGeneral', 'Times New Roman', serif"
TRAZO = 2.5


def mezcla(color, t, fondo="#ffffff"):
    """Mezcla ``color`` con ``fondo``: t = 0 → color, t = 1 → fondo (tintes para rellenos)."""
    a = [int(color[i:i + 2], 16) for i in (1, 3, 5)]
    b = [int(fondo[i:i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(x + (y - x) * t):02x}" for x, y in zip(a, b))


def tinte(color, t=0.85):
    return mezcla(color, t)


def oscuro(color, t=0.35):
    return mezcla(color, t, "#000000")


def _num(v):
    return f"{v:.2f}".rstrip("0").rstrip(".") if isinstance(v, float) else str(v)


def _attrs(**kw):
    partes = []
    for k, v in kw.items():
        if v is None or v is False:
            continue
        k = k.rstrip("_").replace("_", "-")
        partes.append(f'{k}="{escape(_num(v), {chr(34): "&quot;"})}"')
    return " ".join(partes)


class Lienzo:
    """Una ilustración. Coordenadas en unidades SVG; ancho fijo 1000."""

    def __init__(self, alto, *, ancho=1000, titulo=None, etiqueta=None, subtitulo=None, marco=True):
        self.w, self.h = ancho, alto
        self.defs: list[str] = []
        self.cuerpo: list[str] = []
        self._flechas: set[str] = set()
        self._id = 0
        if marco:
            self.cuerpo.append(f'<rect x="1" y="1" width="{ancho - 2}" height="{alto - 2}" rx="22" '
                               f'fill="{SUPERFICIE}" stroke="{REJILLA}" stroke-width="2"/>')
        y = 52
        if etiqueta:
            self.pastilla(40, 30, etiqueta, color=etiqueta_color(etiqueta))
            y = 96
        if titulo:
            self.texto(40, y, titulo, size=27, weight=600)
            y += 30
        if subtitulo:
            self.texto(40, y, subtitulo, size=17, color=TINTA2)
        self.y0 = y + (28 if subtitulo else 10)   # primera coordenada libre bajo la cabecera

    # ------------------------------------------------------------------ básicos
    def nuevo_id(self, base="g"):
        self._id += 1
        return f"{base}{self._id}"

    def add(self, s):
        self.cuerpo.append(s)
        return self

    def grupo(self, contenido, **kw):
        self.cuerpo.append(f"<g {_attrs(**kw)}>{contenido}</g>")

    def rect(self, x, y, w, h, *, fill="none", stroke=None, sw=TRAZO, rx=10, dash=None, opacity=None):
        self.add(f"<rect {_attrs(x=x, y=y, width=w, height=h, rx=rx, fill=fill, stroke=stroke, stroke_width=sw if stroke else None, stroke_dasharray=dash, opacity=opacity)}/>")

    def circulo(self, cx, cy, r, *, fill="none", stroke=None, sw=TRAZO, dash=None, opacity=None):
        self.add(f"<circle {_attrs(cx=cx, cy=cy, r=r, fill=fill, stroke=stroke, stroke_width=sw if stroke else None, stroke_dasharray=dash, opacity=opacity)}/>")

    def elipse(self, cx, cy, rx, ry, *, fill="none", stroke=None, sw=TRAZO, opacity=None, dash=None):
        self.add(f"<ellipse {_attrs(cx=cx, cy=cy, rx=rx, ry=ry, fill=fill, stroke=stroke, stroke_width=sw if stroke else None, opacity=opacity, stroke_dasharray=dash)}/>")

    def linea(self, x1, y1, x2, y2, *, color=TINTA2, sw=TRAZO, dash=None, cap="round", opacity=None):
        self.add(f"<line {_attrs(x1=x1, y1=y1, x2=x2, y2=y2, stroke=color, stroke_width=sw, stroke_dasharray=dash, stroke_linecap=cap, opacity=opacity)}/>")

    def camino(self, d, *, fill="none", stroke=TINTA2, sw=TRAZO, dash=None, opacity=None, join="round", cap="round", flecha=False, flecha_ini=False):
        m_end = f"url(#{self._marcador(stroke)})" if flecha else None
        m_ini = f"url(#{self._marcador(stroke, inicio=True)})" if flecha_ini else None
        self.add(f"<path {_attrs(d=d, fill=fill, stroke=stroke, stroke_width=sw if stroke else None, stroke_dasharray=dash, opacity=opacity, stroke_linejoin=join, stroke_linecap=cap, marker_end=m_end, marker_start=m_ini)}/>")

    def poligono(self, puntos, **kw):
        d = "M " + " L ".join(f"{_num(float(x))} {_num(float(y))}" for x, y in puntos) + " Z"
        self.camino(d, **kw)

    def polilinea(self, puntos, **kw):
        d = "M " + " L ".join(f"{_num(float(x))} {_num(float(y))}" for x, y in puntos)
        self.camino(d, **kw)

    def curva(self, puntos, **kw):
        """Curva suave (Catmull-Rom → Bézier) por los puntos dados."""
        self.camino(suave(puntos), **kw)

    # ------------------------------------------------------------------ flechas
    def _marcador(self, color, inicio=False):
        mid = f"fl{'i' if inicio else ''}{color.strip('#')}"
        if mid not in self._flechas:
            self._flechas.add(mid)
            forma = "M 0 0 L 10 5 L 0 10 z" if not inicio else "M 10 0 L 0 5 L 10 10 z"
            self.defs.append(f'<marker id="{mid}" viewBox="0 0 10 10" refX="{8 if not inicio else 2}" refY="5" '
                             f'markerWidth="5.5" markerHeight="5.5" orient="auto-start-reverse" markerUnits="strokeWidth">'
                             f'<path d="{forma}" fill="{color}"/></marker>')
        return mid

    def flecha(self, x1, y1, x2, y2, *, color=TINTA2, sw=TRAZO, dash=None, doble=False, curvatura=0.0):
        """Flecha recta o curva (``curvatura`` = desplazamiento lateral relativo del punto de control)."""
        if curvatura:
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2
            dx, dy = x2 - x1, y2 - y1
            cx, cy = mx - dy * curvatura, my + dx * curvatura
            d = f"M {_num(float(x1))} {_num(float(y1))} Q {_num(float(cx))} {_num(float(cy))} {_num(float(x2))} {_num(float(y2))}"
        else:
            d = f"M {_num(float(x1))} {_num(float(y1))} L {_num(float(x2))} {_num(float(y2))}"
        self.camino(d, stroke=color, sw=sw, dash=dash, flecha=True, flecha_ini=doble)

    # ------------------------------------------------------------------ texto
    def texto(self, x, y, s, *, size=17, weight=400, color=TINTA, anchor="start", italic=False, font=SANS,
              lineas=1.3, rot=None, opacity=None):
        """Texto; ``s`` puede ser lista (varias líneas) o contener marcado ligero:
        ``_{sub}`` subíndice, ``^{sup}`` superíndice, ``*cursiva*`` y ``**negrita**``."""
        filas = s if isinstance(s, (list, tuple)) else [s]
        tr = f' transform="rotate({rot} {_num(float(x))} {_num(float(y))})"' if rot else ""
        for i, fila in enumerate(filas):
            yy = y + i * size * lineas
            self.add(f'<text {_attrs(x=x, y=yy, font_family=font, font_size=size, font_weight=weight, fill=color, text_anchor=anchor, font_style="italic" if italic else None, opacity=opacity)}{tr}>'
                     f"{_marcado(fila, size)}</text>")

    def mate(self, x, y, s, *, size=22, color=TINTA, anchor="start", weight=400):
        """Texto matemático: letras en cursiva serif con el mismo marcado que :meth:`texto`."""
        self.texto(x, y, s, size=size, color=color, anchor=anchor, font=SERIF, weight=weight)

    def pastilla(self, x, y, s, *, color=AZUL, size=14, anchor="start", relleno=None):
        """Etiqueta en forma de píldora: punto de color + texto en tinta."""
        ancho = ancho_texto(s.upper(), size - 1, bold=True) + 34
        x0 = x - ancho if anchor == "end" else (x - ancho / 2 if anchor == "middle" else x)
        self.rect(x0, y - size - 5, ancho, size + 16, rx=(size + 16) / 2, fill=relleno or tinte(color, 0.86))
        self.circulo(x0 + 15, y - size / 2 + 3, 5, fill=color)
        self.texto(x0 + 27, y + 1, s.upper(), size=size - 1, weight=600, color=TINTA2)
        return ancho

    def nota(self, x, y, w, lineas_txt, *, color=TINTA3, size=15.5, titulo=None, pad=16, fill=None):
        """Recuadro de nota con barra de acento a la izquierda. Devuelve la altura."""
        lineas_txt = lineas_txt if isinstance(lineas_txt, (list, tuple)) else [lineas_txt]
        h = pad * 2 + len(lineas_txt) * size * 1.35 + (size * 1.45 if titulo else 0) - size * 0.35
        self.rect(x, y, w, h, rx=12, fill=fill or tinte(color, 0.9))
        self.rect(x, y, 6, h, rx=3, fill=color)
        yy = y + pad + size * 0.85
        if titulo:
            self.texto(x + 22, yy, titulo, size=size, weight=600)
            yy += size * 1.45
        self.texto(x + 22, yy, list(lineas_txt), size=size, color=TINTA2, lineas=1.35)
        return h

    def numero(self, cx, cy, n, *, color=TINTA, r=15):
        """Círculo numerado (pasos)."""
        self.circulo(cx, cy, r, fill=color)
        self.texto(cx, cy + 6, str(n), size=r * 1.1, weight=600, color="#ffffff", anchor="middle")

    # ------------------------------------------------------------------ salida
    def svg(self):
        defs = f"<defs>{''.join(self.defs)}</defs>" if self.defs else ""
        return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}" height="{self.h}" '
                f'viewBox="0 0 {self.w} {self.h}">{defs}{"".join(self.cuerpo)}</svg>\n')

    def guardar(self, nombre, escala=1.6):
        SALIDA.mkdir(parents=True, exist_ok=True)
        svg = SALIDA / f"{nombre}.svg"
        png = SALIDA / f"{nombre}.png"
        svg.write_text(self.svg(), encoding="utf-8")
        subprocess.run(["rsvg-convert", "--zoom", str(escala), "-o", str(png), str(svg)], check=True)
        _comprimir(png)
        return png


def etiqueta_color(s):
    s = s.lower()
    if "analog" in s:
        return NARANJA
    if "ecuaci" in s or "términ" in s:
        return AZUL
    if "simula" in s or "cálcul" in s or "modelo" in s:
        return VIOLETA
    if "dato" in s or "experim" in s:
        return AGUA
    return TINTA3


def _comprimir(png):
    """Cuantiza a 256 colores (sin pérdida visible en arte vectorial plano) para aligerar el cuaderno."""
    try:
        from PIL import Image
    except ImportError:
        return
    im = Image.open(png).convert("RGBA")
    # MEDIANCUT conserva los tintes pálidos (FASTOCTREE los agrisaba) pero solo acepta RGB:
    # se cuantiza sobre blanco y las esquinas fuera de la tarjeta vuelven a ser transparentes.
    plano = Image.new("RGB", im.size, (255, 255, 255))
    plano.paste(im, mask=im.getchannel("A"))
    q = plano.quantize(colors=255, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    fuera = im.getchannel("A").point(lambda a: 255 if a == 0 else 0)
    pal = q.getpalette()[:255 * 3]
    q.putpalette(pal + [0] * (768 - len(pal)))
    q.paste(255, mask=fuera)
    q.save(png, transparency=255)


# ---------------------------------------------------------------- utilidades de texto
_ANCHOS = {"i": .28, "l": .28, "j": .28, "t": .36, "f": .34, "r": .38, "I": .3, " ": .28, ".": .28, ",": .28,
           ":": .28, ";": .28, "(": .34, ")": .34, "m": .86, "w": .78, "M": .9, "W": .98, "—": 1.0}


def ancho_texto(s, size, bold=False):
    """Estimación del ancho de un texto en Avenir Next (suficiente para cajas y píldoras)."""
    s = s.replace("_{", "").replace("^{", "").replace("}", "").replace("*", "")
    base = sum(_ANCHOS.get(ch, .68 if ch.isupper() else .55) for ch in s)
    return base * size * (1.06 if bold else 1.0)


def _marcado(s, size):
    """Convierte el marcado ligero en <tspan>."""
    out, i = [], 0
    while i < len(s):
        if s.startswith("**", i):
            j = s.index("**", i + 2)
            out.append(f'<tspan font-weight="700">{_marcado(s[i + 2:j], size)}</tspan>')
            i = j + 2
        elif s[i] == "*" and s.find("*", i + 1) > 0:
            j = s.index("*", i + 1)
            out.append(f'<tspan font-style="italic">{escape(s[i + 1:j])}</tspan>')
            i = j + 1
        elif s.startswith("_{", i) or s.startswith("^{", i):
            j = s.index("}", i)
            sub = s[i] == "_"
            dy = size * (0.28 if sub else -0.38)
            out.append(f'<tspan dy="{_num(dy)}" font-size="{_num(size * 0.68)}">{escape(s[i + 2:j])}</tspan>'
                       f'<tspan dy="{_num(-dy)}">​</tspan>')
            i = j + 1
        else:
            j = i
            while j < len(s) and not (s.startswith("**", j) or s[j] == "*" or s.startswith("_{", j) or s.startswith("^{", j)):
                j += 1
            out.append(escape(s[i:j]))
            i = j
    return "".join(out)


# ---------------------------------------------------------------- geometría
def suave(puntos, tension=0.5):
    """Camino SVG suave que pasa por todos los puntos (Catmull-Rom)."""
    p = [(float(x), float(y)) for x, y in puntos]
    d = f"M {_num(p[0][0])} {_num(p[0][1])}"
    for i in range(len(p) - 1):
        p0 = p[i - 1] if i > 0 else p[i]
        p1, p2 = p[i], p[i + 1]
        p3 = p[i + 2] if i + 2 < len(p) else p2
        c1 = (p1[0] + (p2[0] - p0[0]) * tension / 3, p1[1] + (p2[1] - p0[1]) * tension / 3)
        c2 = (p2[0] - (p3[0] - p1[0]) * tension / 3, p2[1] - (p3[1] - p1[1]) * tension / 3)
        d += f" C {_num(c1[0])} {_num(c1[1])} {_num(c2[0])} {_num(c2[1])} {_num(p2[0])} {_num(p2[1])}"
    return d


def muestrear(f, x0, x1, n=120):
    """Lista de puntos (x, f(x))."""
    return [(x0 + (x1 - x0) * i / (n - 1), f(x0 + (x1 - x0) * i / (n - 1))) for i in range(n)]


def perfil_barrera(x0, x1, y_base, y_ts, y_fin, *, centro=0.5, ancho=0.16):
    """Función y(x) de un perfil de reacción: valle - colina - valle (coordenadas SVG, y hacia abajo)."""
    def f(x):
        t = (x - x0) / (x1 - x0)
        s = 1 / (1 + math.exp(-(t - centro) / (ancho * 0.45)))          # escalón suave R → P
        base = y_base + (y_fin - y_base) * s
        altura = (y_ts - (y_base + (y_fin - y_base) * 0.5))
        return base + altura * math.exp(-((t - centro) / ancho) ** 2)
    return f
