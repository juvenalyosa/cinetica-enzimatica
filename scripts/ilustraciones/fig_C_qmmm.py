"""Secciones 7–12: QM/MM, coordenada de reacción, estado de transición, mejoras del cálculo,
la reacción fuera de la enzima y de la barrera a k_cat.

Los números que aparecen en los dibujos se leen de ``data/precalculado/qmmm/resumen.json`` (los
mismos que imprime el cuaderno), de modo que dibujo y cálculo no pueden discrepar.
"""
from __future__ import annotations

import json
import math
import random

import numpy as np

from .lienzo import (AGUA, AMARILLO, AZUL, EJE, MAGENTA, NARANJA, REJILLA, ROOT, TINTA, TINTA2, TINTA3,
                     VERDE, VIOLETA, Lienzo, muestrear, oscuro, perfil_barrera, suave, tinte)

# ---------------------------------------------------------------- datos del cálculo
_RES = json.loads((ROOT / "data" / "precalculado" / "qmmm" / "resumen.json").read_text())["etapas"]
_R, _P, _TS = _RES["reactivo"], _RES["producto"], _RES["dimero"]
_TERMO = _RES["termoquimica"]["298.15"]
_AGUA = _RES["agua"]
_MET = _RES["metodos"]

RT = 0.0019872036 * 298.15            # kcal/mol
KBT_H = 1.380649e-23 * 298.15 / 6.62607015e-34


def _k(dg):
    return KBT_H * math.exp(-dg / RT)


def _dg(k):
    return RT * math.log(KBT_H / k)


# ---------------------------------------------------------------- helpers locales
def _atomo(L, x, y, r, fill, etiqueta=None, *, size=15, dentro=True, dy=0, color_txt=TINTA, borde=None):
    L.circulo(x, y, r, fill=fill, stroke=borde or oscuro(fill, 0.25) if fill != "none" else borde, sw=1.5)
    if etiqueta:
        if dentro:
            L.texto(x, y + size * 0.36 + dy, etiqueta, size=size, weight=600, anchor="middle", color=color_txt)
        else:
            L.texto(x, y + r + size + 4 + dy, etiqueta, size=size, anchor="middle", color=TINTA2)


def _cota(L, x1, x2, y, texto, *, color=TINTA2, size=14):
    """Cota horizontal |←→| con texto encima."""
    L.linea(x1, y - 7, x1, y + 7, color=color, sw=1.5)
    L.linea(x2, y - 7, x2, y + 7, color=color, sw=1.5)
    L.linea(x1 + 4, y, x2 - 4, y, color=color, sw=1.5)
    L.texto((x1 + x2) / 2, y - 10, texto, size=size, anchor="middle", color=TINTA2)


def _tarjeta(L, x, y, w, h, *, fill="#ffffff"):
    L.rect(x, y, w, h, rx=16, fill=fill, stroke=REJILLA, sw=2)


def _muelle(x1, y1, x2, y2, n=7, amp=7):
    """Camino en zigzag (resorte) entre dos puntos."""
    dx, dy = x2 - x1, y2 - y1
    lon = math.hypot(dx, dy)
    ux, uy = dx / lon, dy / lon
    px, py = -uy, ux
    pts = [(x1, y1), (x1 + ux * lon * 0.12, y1 + uy * lon * 0.12)]
    for i in range(n):
        t = 0.12 + 0.76 * (i + 0.5) / n
        s = amp if i % 2 == 0 else -amp
        pts.append((x1 + ux * lon * t + px * s, y1 + uy * lon * t + py * s))
    pts += [(x1 + ux * lon * 0.88, y1 + uy * lon * 0.88), (x2, y2)]
    return "M " + " L ".join(f"{a:.1f} {b:.1f}" for a, b in pts)


# ================================================================ s07 · QM/MM
def ilustracion_qmmm_estadio():
    L = Lienzo(575, etiqueta="Analogía", titulo="QM/MM: enfocar el balón y dejar el estadio de fondo",
               subtitulo="La mecánica cuántica es cara: solo se usa donde ocurre la química")
    rnd = random.Random(7)

    # ---- izquierda: el estadio
    cx, cy = 255, 340
    L.elipse(cx, cy, 205, 150, fill=tinte(TINTA3, 0.82))
    for fila in range(5):                                    # gradas: espectadores desenfocados
        rx, ry = 190 - fila * 13, 137 - fila * 11
        n = 58 - fila * 5
        for i in range(n):
            a = 2 * math.pi * (i + 0.5 * (fila % 2)) / n
            L.circulo(cx + rx * math.cos(a), cy + ry * math.sin(a), 3.6, fill=tinte(TINTA3, 0.45 + 0.08 * rnd.random()))
    L.rect(cx - 118, cy - 70, 236, 140, rx=8, fill=tinte(VERDE, 0.88), stroke=tinte(VERDE, 0.6), sw=1.5)
    L.linea(cx, cy - 70, cx, cy + 70, color=tinte(VERDE, 0.55), sw=1.5)
    L.circulo(cx, cy, 22, stroke=tinte(VERDE, 0.55), sw=1.5)
    # jugadores cercanos al balón (enfocados) y lejanos (apagados)
    for x, y in [(cx - 90, cy - 40), (cx - 80, cy + 45), (cx + 95, cy - 30), (cx + 60, cy + 50)]:
        L.circulo(x, y, 6, fill=tinte(TINTA3, 0.5))
    bx, by = cx + 18, cy - 8
    for dx, dy in [(-30, -14), (26, -22), (-12, 28), (34, 18)]:
        L.circulo(bx + dx, by + dy, 8, fill=AZUL)
    L.circulo(bx, by, 9, fill=NARANJA)
    # visor de la cámara
    x0, y0, x1, y1 = bx - 58, by - 50, bx + 62, by + 50
    for (px, py, sx, sy) in [(x0, y0, 1, 1), (x1, y0, -1, 1), (x0, y1, 1, -1), (x1, y1, -1, -1)]:
        L.polilinea([(px, py + 16 * sy), (px, py), (px + 16 * sx, py)], stroke=TINTA, sw=3)
    L.texto(cx, 530, "Fotografía de un partido", size=17, weight=600, anchor="middle")
    L.texto(cx, 552, "enfocas el balón; las gradas quedan de fondo", size=14.5, color=TINTA2, anchor="middle")

    # ---- flecha de equivalencia
    L.texto(500, 348, "=", size=44, color=TINTA3, anchor="middle")

    # ---- derecha: la enzima
    ex, ey = 745, 340
    blob = [(ex + 205 * math.cos(a) * (1 + 0.07 * math.sin(3 * a)), ey + 150 * math.sin(a) * (1 + 0.06 * math.cos(4 * a)))
            for a in np.linspace(0, 2 * math.pi, 40, endpoint=False)]
    L.camino(suave(blob + blob[:3])[:], fill=tinte(AZUL, 0.9), stroke=tinte(AZUL, 0.6), sw=2)
    for _ in range(170):                                     # átomos MM: puntos pequeños
        a, r = rnd.uniform(0, 2 * math.pi), math.sqrt(rnd.uniform(0.12, 0.9))
        x, y = ex + 190 * r * math.cos(a), ey + 135 * r * math.sin(a)
        L.circulo(x, y, 3, fill=tinte(AZUL, 0.55))
    L.circulo(ex, ey, 62, fill="#ffffff", stroke=NARANJA, sw=2.5)
    # mini molécula en la región QM
    at = [(ex - 30, ey + 8), (ex - 8, ey - 14), (ex + 18, ey - 4), (ex + 32, ey + 20), (ex + 2, ey + 22)]
    for (a1, a2) in [(0, 1), (1, 2), (2, 3), (2, 4)]:
        L.linea(*at[a1], *at[a2], color=TINTA2, sw=3)
    for i, (x, y) in enumerate(at):
        L.circulo(x, y, 9 if i != 2 else 11, fill=NARANJA if i == 2 else tinte(NARANJA, 0.45), stroke="#ffffff", sw=1.5)
    L.texto(ex, 530, "La enzima en QM/MM", size=17, weight=600, anchor="middle")
    L.texto(ex, 552, "78 átomos cuánticos; ~57 000 clásicos", size=14.5, color=TINTA2, anchor="middle")
    # rótulos de las regiones
    L.linea(ex + 44, ey - 44, ex + 72, ey - 106, color=TINTA3, sw=1.5)
    L.texto(ex + 70, ey - 146, "región QM", size=15, weight=600)
    L.texto(ex + 70, ey - 128, "(mecánica cuántica)", size=13.5, color=TINTA2)
    L.linea(ex - 150, ey + 70, ex - 200, ey + 128, color=TINTA3, sw=1.5)
    L.texto(ex - 205, ey + 145, "región MM (campo de fuerza)", size=15, weight=600, anchor="middle")
    return L


def ilustracion_qmmm_regiones():
    L = Lienzo(665, etiqueta="La ecuación, término a término",
               titulo="Dos regiones que se ven: electricidad y «codos»",
               subtitulo="La región QM siente las cargas de la enzima y no puede atravesar sus átomos")
    rnd = random.Random(3)
    cx, cy = 250, 410
    # enlace que cruza la frontera (antes que las cargas, para dejar libre su zona)
    ax_, ay_ = cx - 22, cy - 18
    hx, hy = cx - 72, cy - 74
    mx, my = cx - 112, cy - 118
    # región MM: disco azul claro con cargas puntuales
    L.circulo(cx, cy, 188, fill=tinte(AZUL, 0.9), stroke=tinte(AZUL, 0.6), sw=2)
    puestos = []
    while len(puestos) < 30:
        a, r = rnd.uniform(0, 2 * math.pi), rnd.uniform(124, 170)
        x, y = cx + r * math.cos(a), cy + r * math.sin(a)
        libre = math.hypot(x - mx, y - my) > 36 and math.hypot(x - (hx + mx) / 2, y - (hy + my) / 2) > 30
        if libre and all(math.hypot(x - u, y - v) > 30 for u, v in puestos):
            puestos.append((x, y))
    for i, (x, y) in enumerate(puestos):
        L.circulo(x, y, 10, fill="#ffffff", stroke=tinte(AZUL, 0.3), sw=1.5)
        L.texto(x, y + 5, "+" if i % 3 == 0 else "−", size=15, weight=600, anchor="middle", color=TINTA2)
    # región QM
    L.circulo(cx, cy, 100, fill=tinte(NARANJA, 0.86), stroke=NARANJA, sw=2.5)
    at = [(cx - 50, cy + 14), (cx - 22, cy - 18), (cx + 14, cy - 4), (cx + 46, cy + 18), (cx + 8, cy + 30), (cx + 30, cy - 42)]
    for a1, a2 in [(0, 1), (1, 2), (2, 3), (2, 4), (2, 5)]:
        L.linea(*at[a1], *at[a2], color=TINTA2, sw=3)
    for i, (x, y) in enumerate(at):
        L.circulo(x, y, 12 if i == 2 else 9, fill=NARANJA if i == 2 else tinte(NARANJA, 0.4), stroke="#ffffff", sw=1.5)
    L.texto(cx, cy + 58, "región QM", size=15, weight=600, anchor="middle")
    L.texto(cx, cy + 75, "78 átomos, carga −2", size=12.5, color=TINTA2, anchor="middle")
    # átomo de enlace
    L.linea(ax_, ay_, hx, hy, color=TINTA2, sw=3)
    L.linea(hx, hy, mx, my, color=TINTA3, sw=3, dash="2 6")
    L.circulo(mx, my, 10, fill=tinte(AZUL, 0.5), stroke="#ffffff", sw=1.5)
    L.circulo(hx, hy, 9, fill="#ffffff", stroke=TINTA, sw=2)
    L.texto(hx, hy + 4.5, "H", size=11, weight=600, anchor="middle")
    L.linea(mx - 8, my - 10, 76, 202, color=TINTA3, sw=1.5)
    L.texto(40, 172, ["átomo de enlace: un H tapa el", "enlace cortado en la frontera"], size=13.5, color=TINTA2)
    L.texto(cx, cy + 222, "región MM: ~2300 cargas fijas", size=15, weight=600, anchor="middle")

    # tarjetas de las dos interacciones
    x, w = 505, 455
    _tarjeta(L, x, 180, w, 165)
    L.rect(x, 180, 6, 165, rx=3, fill=NARANJA)
    L.texto(x + 26, 212, "1 · Electrostática (embedding)", size=17, weight=600)
    L.texto(x + 26, 238, ["las cargas MM crean un potencial eléctrico que", "entra en el cálculo cuántico: los electrones",
                          "QM se desplazan hacia las cargas +"], size=14.5, color=TINTA2)
    # icono: nube que se deforma hacia una carga +
    ix, iy = x + 380, 318
    L.camino(f"M {ix - 58} {iy} C {ix - 58} {iy - 22} {ix - 10} {iy - 20} {ix + 4} {iy} C {ix - 10} {iy + 20} {ix - 58} {iy + 22} {ix - 58} {iy} Z",
             fill=tinte(NARANJA, 0.55), stroke=None)
    L.circulo(ix + 42, iy, 11, fill="#ffffff", stroke=tinte(AZUL, 0.3), sw=1.5)
    L.texto(ix + 42, iy + 5, "+", size=15, weight=600, anchor="middle", color=TINTA2)
    L.flecha(ix + 8, iy, ix + 26, iy, color=TINTA2, sw=1.8)
    L.texto(x + 26, 324, "atrae", size=14, weight=600, color=TINTA)

    _tarjeta(L, x, 365, w, 165)
    L.rect(x, 365, 6, 165, rx=3, fill=VIOLETA)
    L.texto(x + 26, 397, "2 · Lennard-Jones (los «codos»)", size=17, weight=600)
    L.texto(x + 26, 423, ["repulsión de corto alcance entre átomos QM y MM;", "sin ella la electrostática, que solo atrae,",
                          "«hundiría» la región QM en la proteína"], size=14.5, color=TINTA2)
    ix, iy = x + 380, 487
    L.circulo(ix - 18, iy, 16, fill=tinte(NARANJA, 0.4))
    L.circulo(ix + 18, iy, 16, fill=tinte(AZUL, 0.5))
    L.flecha(ix - 40, iy, ix - 62, iy, color=TINTA2, sw=2)
    L.flecha(ix + 40, iy, ix + 62, iy, color=TINTA2, sw=2)
    L.texto(x + 26, 510, "repele de cerca", size=14, weight=600)

    # ecuación
    y = 600
    L.mate(x + 12, y, "E  =", size=26)
    L.rect(x + 78, y - 26, 250, 38, rx=10, fill=tinte(NARANJA, 0.86))
    L.mate(x + 203, y, "E_{PM7}[QM en la enzima]", size=21, anchor="middle")
    L.mate(x + 344, y, "+", size=26, anchor="middle", color=TINTA2)
    L.rect(x + 362, y - 26, 90, 38, rx=10, fill=tinte(VIOLETA, 0.86))
    L.mate(x + 407, y, "E_{LJ}", size=21, anchor="middle")
    return L


# ================================================================ s08 · coordenada ξ
def ilustracion_xi_coordenada():
    L = Lienzo(640, etiqueta="La ecuación, término a término",
               titulo="ξ: un solo número que dice «dónde está el fosfato»",
               subtitulo="ξ = d(Pγ–O3β) − d(Pγ–O6): negativa en el reactivo, ≈ 0 en la cima, positiva en el producto")
    esc = 44                                                 # px por Å
    etapas = [("reactivo", AZUL, _R["d_PG_O3B"], _R["d_PG_O6"], _R["xi"], 1),
              ("estado de transición", NARANJA, _TS["d_PG_O3B"], _TS["d_PG_O6"], _TS["xi"], 0),
              ("producto", AGUA, _P["d_PG_O3B"], _P["d_PG_O6"], _P["xi"], -1)]
    for k, (nombre, color, d1, d2, xi, paraguas) in enumerate(etapas):
        px0 = 40 + k * 312
        _tarjeta(L, px0, 165, 296, 305)
        L.pastilla(px0 + 18, 198, nombre, color=color)
        cy = 300
        total = (d1 + d2) * esc
        xo3 = px0 + 148 - total / 2
        xp = xo3 + d1 * esc
        xo6 = xp + d2 * esc
        # enlaces: sólidos si están formados, discontinuos si son parciales o no existen
        for xa, d in [(xo3, d1), (xo6, d2)]:
            if d < 1.9:
                L.linea(xa, cy, xp, cy, color=TINTA2, sw=4)
            elif d < 2.6:
                L.linea(xa, cy, xp, cy, color=TINTA2, sw=3, dash="5 6")
        # paraguas PO3: tres oxígenos
        dx = 17 * paraguas
        tercero = (xp + 1.75 * dx, cy + 9, 7) if dx else (xp + 24, cy + 20, 7)
        for ox, oy, r in [tercero, (xp + dx, cy - 34, 8.5), (xp + dx, cy + 34, 8.5)]:
            L.linea(xp, cy, ox, oy, color=TINTA3, sw=2.5)
            L.circulo(ox, oy, r, fill=tinte(TINTA3, 0.55), stroke="#ffffff", sw=1.2)
        _atomo(L, xp, cy, 17, color, "P", size=15, color_txt="#ffffff")
        _atomo(L, xo3, cy, 12, tinte(TINTA3, 0.35))
        _atomo(L, xo6, cy, 12, tinte(TINTA3, 0.35))
        L.texto(xo3, cy - 22, "O3β", size=13.5, anchor="middle", color=TINTA2)
        L.texto(xo6, cy - 22, "O6", size=13.5, anchor="middle", color=TINTA2)
        L.texto(xo3, cy + 38, "ADP", size=12.5, anchor="middle", color=TINTA3)
        L.texto(xo6, cy + 38, "glucosa", size=12.5, anchor="middle", color=TINTA3)
        _cota(L, xo3, xp, 380, f"{d1:.2f} Å")
        _cota(L, xp, xo6, 380, f"{d2:.2f} Å")
        L.mate(px0 + 148, 440, f"ξ = {d1:.2f} − {d2:.2f} = {xi:+.2f} Å".replace("-", "−"), size=18, anchor="middle")

    # eje de ξ
    y = 545
    xa, xb = 110, 890
    lo, hi = -2.0, 2.0
    X = lambda v: xa + (v - lo) / (hi - lo) * (xb - xa)
    L.flecha(xa - 20, y, xb + 30, y, color=TINTA2, sw=2)
    for v in (-2, -1, 0, 1, 2):
        L.linea(X(v), y - 6, X(v), y + 6, color=TINTA2, sw=1.5)
        L.texto(X(v), y + 26, f"{v:+d}".replace("+0", "0").replace("-", "−"), size=13.5, anchor="middle", color=TINTA3)
    L.texto(xb + 34, y + 5, "ξ (Å)", size=14, color=TINTA2)
    for nombre, color, d1, d2, xi, _ in etapas:
        L.circulo(X(xi), y, 9, fill=color, stroke="#ffffff", sw=2)
    L.texto(X(-1.0), y - 18, "el fosfato sigue en el ATP", size=14, anchor="middle", color=TINTA2)
    L.texto(X(1.0), y - 18, "el fosfato ya está en la glucosa", size=14, anchor="middle", color=TINTA2)
    L.texto(500, 610, "Los enlaces sólidos están formados; los discontinuos, a medio formar o a medio romper.",
            size=14, anchor="middle", color=TINTA3)
    return L


# ================================================================ s09 · estado de transición
def _superficie(x, y):
    return (-1.00 * np.exp(-((x + 1.05) ** 2 + (y + 0.35) ** 2) / 0.42)
            - 0.80 * np.exp(-((x - 1.05) ** 2 + (y - 0.35) ** 2) / 0.42)
            + 0.22 * y ** 2 + 0.02 * x ** 2
            + 0.35 * np.exp(-((x) ** 2 + (y - 1.3) ** 2) / 0.35)
            + 0.30 * np.exp(-((x + 0.1) ** 2 + (y + 1.35) ** 2) / 0.35))


def _grad(x, y, h=1e-4):
    return ((_superficie(x + h, y) - _superficie(x - h, y)) / (2 * h),
            (_superficie(x, y + h) - _superficie(x, y - h)) / (2 * h))


def _camino_minimo():
    """Punto de silla por Newton y descenso de máxima pendiente hacia los dos valles."""
    p = np.array([0.0, 0.0])
    h = 1e-3
    for _ in range(50):
        g = np.array(_grad(*p))
        H = np.array([[(_grad(p[0] + h, p[1])[i] - _grad(p[0] - h, p[1])[i]) / (2 * h),
                       (_grad(p[0], p[1] + h)[i] - _grad(p[0], p[1] - h)[i]) / (2 * h)] for i in range(2)])
        p = p - np.linalg.solve(H, g)
    w, v = np.linalg.eigh(H)
    modo = v[:, 0]
    ramas = []
    for s in (+1, -1):
        q = p + s * 0.02 * modo
        pts = [q.copy()]
        for _ in range(4000):
            g = np.array(_grad(*q))
            if np.linalg.norm(g) < 1e-3:
                break
            q = q - 0.004 * g / max(np.linalg.norm(g), 0.05)
            pts.append(q.copy())
        ramas.append(pts)
    return p, ramas


def _curvas_nivel(L, X0, Y0, W, H, lim, niveles):
    import contourpy
    xs = np.linspace(-lim[0], lim[0], 220)
    ys = np.linspace(-lim[1], lim[1], 160)
    XX, YY = np.meshgrid(xs, ys)
    Z = _superficie(XX, YY)
    to = lambda x, y: (X0 + (x + lim[0]) / (2 * lim[0]) * W, Y0 + (lim[1] - y) / (2 * lim[1]) * H)
    gen = contourpy.contour_generator(xs, ys, Z, fill_type=contourpy.FillType.OuterCode,
                                      line_type=contourpy.LineType.Separate)
    rampa = [tinte(AZUL, t) for t in np.linspace(0.35, 0.97, len(niveles) - 1)]
    for (lo, hi), color in zip(zip(niveles[:-1], niveles[1:]), rampa):
        pts_list, codes_list = gen.filled(lo, hi)
        for pts, codes in zip(pts_list, codes_list):
            d = []
            for (x, y), c in zip(pts, codes):
                sx, sy = to(x, y)
                d.append(f"{'M' if c == 1 else 'L'} {sx:.1f} {sy:.1f}" + (" Z" if c == 79 else ""))
            L.camino(" ".join(d), fill=color, stroke=None)
    for nv in niveles[1:-1]:
        for linea in gen.lines(nv):
            L.polilinea([to(x, y) for x, y in linea], stroke="#ffffff", sw=1.2, opacity=0.8)
    return to


def ilustracion_punto_silla():
    L = Lienzo(655, etiqueta="Analogía", titulo="El estado de transición es un paso de montaña",
               subtitulo="Mapa de curvas de nivel: dos valles (azul oscuro) separados por un collado")
    X0, Y0, W, H = 40, 160, 560, 400
    L.rect(X0, Y0, W, H, rx=14, fill=tinte(AZUL, 0.97))
    niveles = list(np.linspace(-1.0, 0.62, 14))
    niveles = [-1.2] + niveles[1:] + [3.0]
    to = _curvas_nivel(L, X0, Y0, W, H, (2.2, 1.6), niveles)
    # marco con esquinas redondeadas por encima del relleno
    L.rect(X0, Y0, W, H, rx=14, stroke=REJILLA, sw=4)

    ts, ramas = _camino_minimo()
    ramas = sorted(ramas, key=lambda r: r[-1][0])          # la rama que acaba a la izquierda (x < 0) es el reactivo
    camino = ramas[0][::-1] + [ts] + ramas[1]
    pix = [to(*p) for p in camino]
    L.polilinea(pix, stroke=TINTA, sw=2.5, dash="1 7")
    # cuentas del NEB equiespaciadas a lo largo del camino
    lon = np.r_[0, np.cumsum([math.dist(pix[i], pix[i + 1]) for i in range(len(pix) - 1)])]
    n = 11
    cuentas = [pix[int(np.searchsorted(lon, t))] if t < lon[-1] else pix[-1] for t in np.linspace(0, lon[-1], n)]
    ts_pix = to(*ts)
    cuentas[n // 2] = ts_pix
    for a, b in zip(cuentas[:-1], cuentas[1:]):
        L.camino(_muelle(*a, *b, n=4, amp=4), stroke=TINTA2, sw=1.5)
    for i, (x, y) in enumerate(cuentas):
        if i == n // 2:
            continue
        L.circulo(x, y, 7, fill="#ffffff", stroke=TINTA, sw=2)
    L.circulo(*ts_pix, 11, fill=NARANJA, stroke="#ffffff", sw=2.5)
    rx_, ry_ = cuentas[0]
    px_, py_ = cuentas[-1]
    L.circulo(rx_, ry_, 11, fill=AZUL, stroke="#ffffff", sw=2.5)
    L.circulo(px_, py_, 11, fill=AGUA, stroke="#ffffff", sw=2.5)
    for (x, y, t, dy) in [(rx_, ry_, "reactivo", 40), (px_, py_, "producto", -30)]:
        L.rect(x - 44, y + dy - 18, 88, 26, rx=13, fill="#ffffff", opacity=0.92)
        L.texto(x, y + dy, t, size=15, weight=600, anchor="middle")
    L.linea(ts_pix[0] + 10, ts_pix[1] - 10, ts_pix[0] + 60, ts_pix[1] - 90, color=TINTA, sw=1.5)
    L.rect(ts_pix[0] + 8, ts_pix[1] - 138, 190, 50, rx=10, fill="#ffffff")
    L.texto(ts_pix[0] + 20, ts_pix[1] - 116, "punto de silla", size=15, weight=600)
    L.texto(ts_pix[0] + 20, ts_pix[1] - 97, "(imagen trepadora del NEB)", size=13, color=TINTA2)
    L.texto(X0 + 18, Y0 + H + 30, "Las cuentas unidas por muelles (el NEB) se asientan sobre el camino de mínima energía entre los dos valles.",
            size=13.5, color=TINTA3)

    # ---- derecha: los tres pasos
    x, w = 630, 330
    L.texto(x, 185, "Cómo se encuentra", size=18, weight=600)
    pasos = [
        ("NEB con imagen trepadora", ["una cadena de geometrías unidas por", "muelles se relaja sobre el camino;", "la cuenta central trepa a la cima"]),
        ("Método del dímero", ["dos puntos muy juntos «tantean» la", "curvatura y refinan la cima hasta", "que la fuerza es cero"]),
        ("Frecuencias", ["en una silla hay exactamente una", "frecuencia imaginaria: la dirección", "en la que la geometría «cae»"]),
    ]
    y = 225
    for i, (tit, txt) in enumerate(pasos, 1):
        L.numero(x + 16, y + 6, i, color=NARANJA if i == 1 else TINTA2)
        L.texto(x + 44, y + 12, tit, size=16, weight=600)
        L.texto(x + 44, y + 36, txt, size=14, color=TINTA2)
        y += 125
    return L


def ilustracion_frecuencia_imaginaria():
    L = Lienzo(520, etiqueta="La ecuación, término a término",
               titulo="La prueba del punto de silla: una sola frecuencia imaginaria",
               subtitulo="Cada frecuencia es un «resorte» del sistema; en la silla hay un resorte al revés")
    paneles = [
        ("En un valle (R o P)", "todas las direcciones suben", "todas las frecuencias son reales", AZUL, +1, +1),
        ("Silla, a lo ancho del paso", "subes hacia las montañas", "frecuencia real (vibra)", NARANJA, +1, None),
        ("Silla, a lo largo del camino", "bajas hacia R o hacia P", f"frecuencia imaginaria: {_RES['frecuencias']['frecuencias_mas_bajas'][0]:.0f} cm⁻¹".replace("-", "−"), NARANJA, -1, None),
    ]
    for k, (tit, sub, res, color, curv, _) in enumerate(paneles):
        x0 = 40 + k * 312
        _tarjeta(L, x0, 160, 296, 320)
        L.texto(x0 + 148, 196, tit, size=16, weight=600, anchor="middle")
        cx, cy = x0 + 148, 352
        f = lambda t, c=curv: cy - c * 90 * (t / 110) ** 2 + (0 if c > 0 else -40)
        pts = muestrear(lambda xx: f(xx - cx), cx - 110, cx + 110, 60)
        L.poligono(pts + [(cx + 110, 425), (cx - 110, 425)], fill=tinte(TINTA3, 0.88), stroke=None)
        L.curva(pts, stroke=TINTA2, sw=3)
        yb = f(0) - 14
        if curv > 0:
            L.circulo(cx, yb, 12, fill=color)
            L.flecha(cx - 18, yb - 26, cx - 58, yb - 20, color=TINTA3, sw=2)
            L.flecha(cx + 18, yb - 26, cx + 58, yb - 20, color=TINTA3, sw=2)
            L.texto(cx, yb - 42, "oscila", size=13.5, anchor="middle", color=TINTA2)
        else:
            L.circulo(cx, yb, 12, fill=color)
            L.flecha(cx - 20, yb + 2, cx - 72, f(-72) - 12, color=TINTA, sw=2)
            L.flecha(cx + 20, yb + 2, cx + 72, f(72) - 12, color=TINTA, sw=2)
            L.texto(cx - 112, f(-110) + 2, "R", size=15, weight=600, anchor="end")
            L.texto(cx + 112, f(110) + 2, "P", size=15, weight=600)
        L.texto(cx, 448, sub, size=14, color=TINTA2, anchor="middle")
        L.texto(cx, 470, res, size=14, weight=600, anchor="middle")
    L.texto(500, 505, "Animar el modo imaginario muestra la química: el fósforo va y viene entre O3β y O6.",
            size=14, anchor="middle", color=TINTA3)
    return L


# ================================================================ s10 · mejoras
def ilustracion_tres_mejoras():
    dE, dZPE = _TERMO["dE_kcal"], _TERMO["dZPE_kcal"]
    dHt = _TERMO["dH_vib_kcal"] - dZPE                       # parte térmica (sin ZPE)
    mTdS = -298.15 * _TERMO["dS_vib_cal"] / 1000
    dG = _TERMO["dG_kcal"]
    L = Lienzo(640, etiqueta="Modelo", titulo="De ΔE‡ a ΔG‡: tres preguntas de un buen científico",
               subtitulo="Cada pregunta es un cálculo extra; solo la primera cambia el número, las otras dicen cuánto fiarse")
    tarjetas = [
        ("¿Y las vibraciones?", ["los átomos nunca están quietos:", "energía de punto cero y entropía"], AGUA),
        ("¿Y otra postura de la enzima?", ["repetir en 4 instantáneas de", "la MD (250–1000 ps)"], AZUL),
        ("¿Y si el método se equivoca?", ["recalcular con PM6-D3H4 en", "las mismas geometrías"], MAGENTA),
    ]
    for k, (tit, txt, color) in enumerate(tarjetas):
        x0 = 40 + k * 312
        _tarjeta(L, x0, 160, 296, 250)
        L.numero(x0 + 30, 192, k + 1, color=color)
        L.texto(x0 + 54, 198, tit, size=16, weight=600)
        L.texto(x0 + 24, 368, txt, size=14, color=TINTA2)
        cx, cy = x0 + 148, 275
        if k == 0:                                           # dos átomos unidos por un resorte que vibra
            L.circulo(cx - 70, cy, 22, fill=tinte(color, 0.4))
            L.circulo(cx + 70, cy, 22, fill=tinte(color, 0.4))
            L.camino(_muelle(cx - 48, cy, cx + 48, cy, n=9, amp=10), stroke=TINTA2, sw=2.5)
            for s in (-1, 1):
                for j in range(3):
                    r = 30 + j * 9
                    ang = 0.6
                    x1 = cx + s * 70 + s * r * math.cos(ang)
                    L.camino(f"M {cx + s * 70 + s * r * math.cos(ang):.1f} {cy - r * math.sin(ang):.1f} "
                             f"A {r} {r} 0 0 {1 if s > 0 else 0} {cx + s * 70 + s * r * math.cos(ang):.1f} {cy + r * math.sin(ang):.1f}",
                             stroke=TINTA3, sw=1.8, opacity=1 - j * 0.25)
        elif k == 1:                                         # tres posturas de la «almeja»
            for j, (ang, op) in enumerate([(-18, .35), (-8, .6), (2, 1.0)]):
                a = math.radians(ang)
                L.camino(f"M {cx - 80} {cy + 30} Q {cx} {cy + 70} {cx + 80} {cy + 30} Z", fill=tinte(color, 0.7), stroke=None, opacity=op)
                tx, ty = cx - 80, cy + 20
                pts = [(cx - 80, cy + 20), (cx, cy - 30), (cx + 80, cy + 20)]
                rot = [(tx + (px - tx) * math.cos(a) - (py - ty) * math.sin(a), ty + (px - tx) * math.sin(a) + (py - ty) * math.cos(a)) for px, py in pts]
                L.camino(f"M {rot[0][0]:.1f} {rot[0][1]:.1f} Q {rot[1][0]:.1f} {rot[1][1]:.1f} {rot[2][0]:.1f} {rot[2][1]:.1f}",
                         stroke=tinte(color, 1 - op) if op < 1 else color, sw=3)
            L.circulo(cx, cy + 22, 9, fill=NARANJA)
        else:                                                # el mismo cálculo medido con dos «reglas»
            esc = 8.0
            for j, (nom, v, c) in enumerate([("PM7", _MET["PM7"]["barrera_kcal"], TINTA2),
                                             ("PM6-D3H4", _MET["PM6-D3H4"]["barrera_kcal"], color)]):
                yy = cy - 30 + j * 56
                L.texto(x0 + 24, yy - 10, nom, size=13.5, weight=600)
                L.rect(x0 + 24, yy, v * esc, 18, rx=5, fill=c)
                L.texto(x0 + 32 + v * esc, yy + 15, f"{v:.1f} kcal/mol", size=13.5, color=TINTA2)
    # barra de la suma
    y = 470
    L.texto(40, y, "Solo la pregunta 1 corrige la barrera (kcal/mol, 25 °C):", size=16, weight=600)
    # dos decimales: con uno solo la suma de los redondeos (19.8) no coincidiría con ΔG‡ (19.75)
    chips = [(f"ΔE‡ {dE:.2f}", TINTA3), (f"ΔZPE {dZPE:+.2f}", AGUA), (f"ΔH_{{vib}} {dHt:+.2f}", VIOLETA),
             (f"−TΔS‡ {mTdS:+.2f}", MAGENTA), (f"ΔG‡ {dG:.2f}", NARANJA)]
    x = 40
    for i, (txt, color) in enumerate(chips):
        txt = txt.replace("-", "−")
        w = 150 if i in (0, 4) else 150
        L.rect(x, y + 22, w, 50, rx=12, fill=tinte(color, 0.86))
        L.rect(x, y + 66, w, 6, rx=3, fill=color)
        L.mate(x + w / 2, y + 55, txt, size=20, anchor="middle")
        if i < 3:
            L.texto(x + w + 17, y + 56, "+", size=24, color=TINTA2, anchor="middle")
        elif i == 3:
            L.texto(x + w + 17, y + 56, "=", size=24, color=TINTA2, anchor="middle")
        x += w + 34
    L.texto(40, y + 120, ["El resorte que se rompe pierde su vibración (ΔZPE < 0 baja la barrera); el estado de transición es",
                         "más rígido que el reactivo (ΔS‡ < 0), así que −TΔS‡ la sube. Casi se compensan."],
            size=14, color=TINTA2)
    return L


# ================================================================ s11 · dentro vs fuera
def ilustracion_dentro_fuera():
    L = Lienzo(610, etiqueta="Analogía", titulo="¿Cuánto hace el resto de la enzima?",
               subtitulo="Sacamos el sitio activo de la proteína, lo ponemos en agua y volvemos a medir la colina")
    # ---- izquierda: dos perfiles (suaves, con el máximo exactamente en el TS)
    X0, X1, Yb = 110, 560, 395
    esc = 9.2                                                # px por kcal/mol
    xts = X0 + (X1 - X0) * 0.47

    def perfil(b, dr):
        def f(x):
            if x <= xts:
                u = (x - X0) / (xts - X0)
                return Yb - b * esc * (3 * u ** 2 - 2 * u ** 3)
            u = (x - xts) / (X1 - xts)
            return Yb - (b + (dr - b) * (3 * u ** 2 - 2 * u ** 3)) * esc
        return f

    L.linea(X0 - 20, Yb, X1 + 10, Yb, color=EJE, sw=1.5)
    for v in (-10, 0, 10, 20):
        L.linea(X0 - 26, Yb - v * esc, X0 - 20, Yb - v * esc, color=TINTA3, sw=1.5)
        L.texto(X0 - 32, Yb - v * esc + 5, str(v).replace("-", "−"), size=13, anchor="end", color=TINTA3)
    L.texto(X0 - 66, Yb - 60, "energía (kcal/mol)", size=13, color=TINTA3, rot=-90, anchor="middle")
    curvas = [("dentro de la enzima", _TS["barrera_kcal"], _P["dE_reaccion_kcal"], AZUL, None),
              ("sitio activo en agua", _AGUA["barrera_kcal"], _AGUA["dE_reaccion"], TINTA3, "7 6")]
    for nombre, b, dr, color, dash in curvas:
        L.curva(muestrear(perfil(b, dr), X0, X1, 90), stroke=color, sw=3.5 if not dash else 3, dash=dash)
        L.circulo(xts, Yb - b * esc, 6, fill=color)
    b1, b2 = _TS["barrera_kcal"], _AGUA["barrera_kcal"]
    y1, y2 = Yb - b1 * esc, Yb - b2 * esc
    L.camino(f"M {xts - 30} {y2} L {xts - 38} {y2} L {xts - 38} {y1} L {xts - 30} {y1}", stroke=TINTA, sw=1.8)
    L.linea(xts - 38, (y1 + y2) / 2, xts - 70, 168, color=TINTA3, sw=1.5)
    L.texto(xts - 74, 162, f"diferencia: {b2 - b1:.1f} kcal/mol", size=15, weight=600, anchor="end")
    L.texto(X0 - 4, Yb + 22, "R", size=14, weight=600, color=TINTA2)
    L.texto(X1 + 4, Yb - _P["dE_reaccion_kcal"] * esc - 10, "P", size=14, weight=600, color=TINTA2, anchor="end")
    # leyenda
    L.linea(X0, 530, X0 + 40, 530, color=AZUL, sw=3.5)
    L.texto(X0 + 50, 535, f"dentro de la enzima: barrera {b1:.1f} kcal/mol", size=14, color=TINTA2)
    L.linea(X0, 558, X0 + 40, 558, color=TINTA3, sw=3, dash="7 6")
    L.texto(X0 + 50, 563, f"el mismo clúster en agua (COSMO): {b2:.1f} kcal/mol", size=14, color=TINTA2)

    # ---- derecha: qué conserva el clúster
    x, w = 620, 340
    _tarjeta(L, x, 170, w, 390)
    L.texto(x + 22, 205, "El clúster en agua conserva…", size=16, weight=600)
    items = [("glucosa + trifosfato", NARANJA), ("Asp205, la base", AZUL), ("Lys169, la carga +", VIOLETA), ("Mg²⁺ y su agua", VERDE)]
    for i, (t, c) in enumerate(items):
        yy = 240 + i * 34
        L.circulo(x + 34, yy - 5, 8, fill=c)
        L.texto(x + 52, yy, t, size=15, color=TINTA2)
    L.linea(x + 22, 390, x + w - 22, 390, color=REJILLA, sw=1.5)
    L.texto(x + 22, 422, "…pero pierde", size=16, weight=600)
    L.circulo(x + 34, 452, 8, fill="none", stroke=TINTA3, sw=2, dash="3 3")
    L.texto(x + 52, 457, ["el campo eléctrico y el «molde»", "del resto de la proteína"], size=15, color=TINTA2)
    L.texto(x + 22, 535, "Por eso las dos colinas se parecen.", size=14, color=TINTA3)
    return L


# ================================================================ s12 · barrera → k_cat
def ilustracion_escalera_kcat():
    L = Lienzo(620, etiqueta="Datos reales", titulo="De la barrera al número de recambio",
               subtitulo="Una sola regla, dos escalas: cada 1.36 kcal/mol hacia la derecha, la reacción es 10 veces más lenta")
    xa, xb = 90, 910
    lo, hi = 13.0, 23.0
    X = lambda g: xa + (g - lo) / (hi - lo) * (xb - xa)
    y = 440
    # bandas de ×10
    k0 = math.floor(math.log10(_k(lo)))
    for i, e in enumerate(range(k0, -6, -1)):
        g1, g2 = _dg(10 ** e), _dg(10 ** (e - 1))
        if g1 > hi:
            break
        L.rect(X(max(g1, lo)), y - 22, X(min(g2, hi)) - X(max(g1, lo)), 44, rx=0, fill=tinte(NARANJA, 0.92 if i % 2 else 0.84))
    L.rect(xa, y - 22, xb - xa, 44, rx=6, stroke=REJILLA, sw=1.5)
    # escala superior: ΔG‡
    for g in range(13, 24):
        L.linea(X(g), y - 22, X(g), y - 30, color=TINTA2, sw=1.5)
        L.texto(X(g), y - 38, str(g), size=13.5, anchor="middle", color=TINTA2)
    L.texto(xa, y - 62, "barrera ΔG‡ (kcal/mol)", size=14, weight=600)
    # escala inferior: k
    for e in range(k0, -6, -1):
        g = _dg(10 ** e)
        if lo <= g <= hi:
            L.linea(X(g), y + 22, X(g), y + 30, color=TINTA2, sw=1.5)
            L.texto(X(g), y + 50, f"10^{{{e}}}".replace("-", "−"), size=14, anchor="middle", color=TINTA2)
    L.texto(xa, y + 82, "k = k_{cat} por Eyring (s⁻¹, 25 °C)", size=14, weight=600)

    exp_g = _dg(65.8)
    # (ΔG‡, nombre, valor, color, y de la etiqueta, anclaje)
    marcas = [
        (exp_g, "experimento", f"k_{{cat}} ≈ 66 s⁻¹ → ΔG‡ ≈ {exp_g:.1f}", TINTA, 200, "middle"),
        (18.3, "QM/MM publicada (Zhang 2009)", "18.3", VIOLETA, 262, "end"),
        (_TS["barrera_kcal"], "ΔE‡ de este cuaderno", f"{_TS['barrera_kcal']:.1f}", AZUL, 200, "start"),
        (_TERMO["dG_kcal"], "ΔG‡ con vibraciones", f"{_TERMO['dG_kcal']:.2f}", NARANJA, 300, "start"),
        (_AGUA["barrera_kcal"], "sitio activo en agua", f"{_AGUA['barrera_kcal']:.1f}", TINTA3, 262, "start"),
    ]
    for g, nombre, val, color, ytxt, anchor in marcas:
        xl = X(g) + {"start": -6, "end": 6, "middle": 0}[anchor]
        L.linea(X(g), ytxt + 30 if anchor == "middle" else ytxt - 16, X(g), y - 22, color=color, sw=2,
                dash="2 4" if color == TINTA3 else None)
        L.circulo(X(g), y, 8, fill=color, stroke="#ffffff", sw=2)
        dx = {"start": 10, "end": -10, "middle": 0}[anchor]
        L.texto(X(g) + dx, ytxt, nombre, size=14, weight=600, anchor=anchor)
        L.texto(X(g) + dx, ytxt + 20, val, size=13.5, anchor=anchor, color=TINTA2)
    # llave: la discrepancia
    g1, g2 = exp_g, _TS["barrera_kcal"]
    yb = y + 110
    L.camino(f"M {X(g1)} {yb - 10} L {X(g1)} {yb} L {X(g2)} {yb} L {X(g2)} {yb - 10}", stroke=TINTA2, sw=1.8)
    L.texto((X(g1) + X(g2)) / 2, yb + 24, "los cálculos QM/MM ponen 3–4 kcal/mol de más ⇒ una reacción ~10^{2}–10^{3} veces más lenta",
            size=14.5, anchor="middle", color=TINTA2)
    return L
