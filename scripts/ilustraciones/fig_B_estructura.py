"""Secciones 3–5: la estructura del complejo, su preparación y la dinámica molecular."""
from __future__ import annotations

import math

from .lienzo import (AGUA, AMARILLO, AZUL, EJE, NARANJA, REJILLA, ROJO, SANS, SERIF, TINTA, TINTA2, TINTA3,
                     VERDE, VIOLETA, Lienzo, _num, tinte)

FOSFATO = AMARILLO


# ---------------------------------------------------------------- helpers locales
def atomo(L, x, y, r, etiqueta=None, color=TINTA3, size=None, relleno=0.72, texto_color=TINTA):
    L.circulo(x, y, r, fill=tinte(color, relleno), stroke=color, sw=2)
    if etiqueta:
        griega = any(ch in etiqueta for ch in "αβγ")
        L.texto(x, y + (size or r * 0.8) * 0.36, etiqueta, size=(size or r * 0.8) * (1.12 if griega else 1),
                weight=600, anchor="middle", color=texto_color, font=SERIF if griega else SANS)


def enlace(L, x1, y1, x2, y2, sw=3.2, color=TINTA2):
    L.linea(x1, y1, x2, y2, color=color, sw=sw)


def hexagono(L, cx, cy, r, color, rot=0, sw=2.5, relleno=0.8):
    pts = [(cx + r * math.cos(math.radians(rot + 60 * i)), cy + r * math.sin(math.radians(rot + 60 * i)))
           for i in range(6)]
    L.poligono(pts, fill=tinte(color, relleno), stroke=color, sw=sw)
    return pts


def resorte(L, x1, y1, x2, y2, n=7, amp=7, color=TINTA2, sw=2.2):
    dx, dy = x2 - x1, y2 - y1
    lon = math.hypot(dx, dy)
    ux, uy = dx / lon, dy / lon
    px, py = -uy, ux
    pts = [(x1, y1), (x1 + ux * lon * 0.12, y1 + uy * lon * 0.12)]
    for i in range(1, 2 * n):
        t = 0.12 + 0.76 * i / (2 * n)
        s = amp if i % 2 else -amp
        pts.append((x1 + ux * lon * t + px * s, y1 + uy * lon * t + py * s))
    pts += [(x1 + ux * lon * 0.88, y1 + uy * lon * 0.88), (x2, y2)]
    L.polilinea(pts, stroke=color, sw=sw)


def elipse_rot(L, cx, cy, rx, ry, rot, fill, stroke, sw=2.5, pivote=None):
    px, py = pivote or (cx, cy)
    L.add(f'<ellipse cx="{_num(float(cx))}" cy="{_num(float(cy))}" rx="{rx}" ry="{ry}" fill="{fill}" stroke="{stroke}" '
          f'stroke-width="{sw}" transform="rotate({rot} {_num(float(px))} {_num(float(py))})"/>')


def candado(L, x, y, s=1.0, color=TINTA2):
    L.camino(f"M {x - 9 * s} {y} L {x - 9 * s} {y - 9 * s} A {9 * s} {9 * s} 0 0 1 {x + 9 * s} {y - 9 * s} L {x + 9 * s} {y}",
             stroke=color, sw=3 * s)
    L.rect(x - 14 * s, y - 2 * s, 28 * s, 22 * s, rx=4 * s, fill=color)
    L.circulo(x, y + 8 * s, 3 * s, fill="#ffffff")


def check(L, x, y, color="#0ca30c", r=15):
    L.circulo(x, y, r, fill=color)
    L.polilinea([(x - r * 0.45, y + r * 0.02), (x - r * 0.1, y + r * 0.38), (x + r * 0.5, y - r * 0.35)],
                stroke="#ffffff", sw=3)


def cruz(L, x, y, color="#d03b3b", r=15):
    L.circulo(x, y, r, fill=color)
    d = r * 0.38
    L.linea(x - d, y - d, x + d, y + d, color="#ffffff", sw=3)
    L.linea(x - d, y + d, x + d, y - d, color="#ffffff", sw=3)


def agua_molecula(L, x, y, s=1.0, rot=0):
    a = math.radians(rot)
    for ang in (-52, 52):
        b = a + math.radians(ang) + math.pi / 2
        hx, hy = x + 11 * s * math.cos(b), y + 11 * s * math.sin(b)
        L.linea(x, y, hx, hy, color=EJE, sw=2)
        L.circulo(hx, hy, 3.6 * s, fill="#ffffff", stroke=TINTA3, sw=1.2)
    L.circulo(x, y, 5.5 * s, fill=tinte(AZUL, 0.55), stroke=None)


# ---------------------------------------------------------------- sección 3
def ilustracion_almeja_dominios():
    L = Lienzo(520, etiqueta="Analogía", titulo="La glucoquinasa se cierra sobre la glucosa como una almeja",
               subtitulo="Dos dominios unidos por una bisagra: abierta atrapa poco; cerrada, alinea los reactivos")
    for cx, cerrada in ((250, False), (750, True)):
        cy = 330
        hx, hy = cx + 150, cy - 18                       # bisagra, a la derecha
        elipse_rot(L, cx, cy + 58, 175, 78, 0, tinte(AZUL, 0.62), AZUL)          # dominio grande
        elipse_rot(L, cx + 10, cy - 34, 140, 50, 0 if cerrada else 27, tinte(AZUL, 0.82), AZUL,
                   pivote=(hx, hy))                                              # dominio pequeño
        # sustratos sobre el dominio grande (el ATP y el Mg²⁺ ya están unidos)
        for i in range(3):
            atomo(L, cx + 22 + i * 26, cy + 12, 10, color=FOSFATO, relleno=0.5)
        L.circulo(cx + 48, cy + 40, 8, fill=VERDE)
        if cerrada:
            hexagono(L, cx - 30, cy + 10, 21, NARANJA, rot=30, relleno=0.55)
        else:
            hexagono(L, cx - 140, cy - 12, 21, NARANJA, rot=30, relleno=0.55)
            L.flecha(cx - 114, cy - 6, cx - 48, cy + 8, color=TINTA3, sw=2, dash="4 5")
        L.circulo(hx, hy, 9, fill=TINTA2)
    L.texto(250 - 140, 290, "glucosa", size=15, color=TINTA2, anchor="middle")
    L.texto(282, 240, "dominio pequeño", size=15, color=TINTA2, anchor="middle")
    L.texto(250, 430, "dominio grande", size=15, color=TINTA2, anchor="middle")
    L.texto(250 + 162, 318, "bisagra", size=14, color=TINTA2)
    L.texto(250, 495, "Abierta: poco activa", size=18, weight=600, anchor="middle")
    L.texto(750, 495, "Cerrada: lista para reaccionar", size=18, weight=600, anchor="middle")
    L.flecha(465, 250, 555, 250, color=TINTA, sw=3)
    L.texto(510, 232, "se cierra", size=14, color=TINTA2, anchor="middle")
    # leyenda
    lx, ly = 610, 190
    for i, (col, txt, forma) in enumerate([(NARANJA, "glucosa", "hex"), (FOSFATO, "fosfatos del ATP", "c"),
                                           (VERDE, "Mg²⁺", "m")]):
        x = lx + [0, 100, 250][i]
        if forma == "hex":
            hexagono(L, x, ly - 5, 8, col, rot=30, relleno=0.55, sw=1.8)
        elif forma == "c":
            atomo(L, x, ly - 5, 7, color=col, relleno=0.5)
        else:
            L.circulo(x, ly - 5, 6, fill=col)
        L.texto(x + 14, ly, txt, size=14, color=TINTA2)
    return L


SUPER = "#fcfcfb"


def ilustracion_sitio_activo():
    L = Lienzo(630, etiqueta="Datos reales", titulo="El sitio activo en el cristal (3FGU): todo listo para reaccionar",
               subtitulo="Esquema plano; las distancias son las medidas en la estructura cristalina")
    y = 330
    # --- glucosa
    gx, gy = 185, 360
    pts = hexagono(L, gx, gy, 46, NARANJA, rot=30, relleno=0.78)
    c6 = (300, 312)
    enlace(L, pts[5][0], pts[5][1], *c6)
    o6 = (400, y)
    enlace(L, *c6, *o6)
    L.texto(gx, gy + 7, "glucosa", size=15, weight=600, anchor="middle")
    L.texto(c6[0], c6[1] - 16, "C6", size=13, color=TINTA2, anchor="middle")
    hO6 = (418, 392)
    enlace(L, *o6, *hO6, sw=2.6)
    atomo(L, *hO6, 10, "H", color=TINTA3, relleno=0.95, size=12)
    atomo(L, *o6, 17, "O6", color=NARANJA, size=13)

    # --- Asp205
    od1 = (430, 465)
    cg = (482, 500)
    od2 = (482, 555)
    enlace(L, *od1, *cg)
    enlace(L, *cg, *od2)
    enlace(L, *cg, 560, 500)
    atomo(L, *od1, 15, "O", color=AZUL, size=13)
    atomo(L, *od2, 15, "O", color=AZUL, size=13)
    L.texto(570, 506, "Asp205", size=16, weight=600)
    L.texto(570, 526, "base catalítica (−)", size=14, color=TINTA2)
    # O6···OD1
    L.linea(o6[0] + 6, o6[1] + 16, od1[0] - 4, od1[1] - 15, color=AZUL, sw=2, dash="3 6")
    L.texto(365, 440, "2.5 Å", size=15, weight=600, anchor="end")
    L.texto(365, 458, "puente de H", size=13, color=TINTA2, anchor="end")

    # --- trifosfato
    pg, o3b, pb, o3a, pa = (560, y), (635, y), (710, y), (785, y), (860, y)
    for a, b in ((pg, o3b), (o3b, pb), (pb, o3a), (o3a, pa), (pa, (930, y))):
        enlace(L, *a, *b)
    terminales = {pg: [(560, 262), (522, 385), (598, 385)], pb: [(710, 262), (710, 398)], pa: [(860, 262), (860, 398)]}
    for p, os_ in terminales.items():
        for o in os_:
            enlace(L, *p, *o, sw=2.6)
            atomo(L, *o, 11, "O", color=FOSFATO, relleno=0.8, size=11)
    for p, lab in ((pg, "Pγ"), (pb, "Pβ"), (pa, "Pα")):
        atomo(L, *p, 22, lab, color=FOSFATO, relleno=0.45, size=16)
    for o in (o3b, o3a):
        atomo(L, *o, 12, "O", color=FOSFATO, relleno=0.8, size=11)
    L.texto(o3b[0], o3b[1] - 22, "O3β", size=13, color=TINTA2, anchor="middle")
    L.texto(935, y + 5, "…", size=18, color=TINTA2)
    L.texto(785, 205, "trifosfato del ATP", size=15, weight=600, anchor="middle")
    L.texto(785, 225, "(→ adenosina)", size=13, color=TINTA3, anchor="middle")

    # Pγ···O6
    L.linea(o6[0] + 18, y, pg[0] - 23, y, color=NARANJA, sw=2.4, dash="3 6")
    L.texto(480, y - 16, "2.7 Å", size=15, weight=600, anchor="middle")
    L.texto(480, y + 26, "aún sin", size=12, color=TINTA2, anchor="middle")
    L.texto(480, y + 40, "enlace", size=12, color=TINTA2, anchor="middle")

    # --- Mg2+
    mg = (655, 455)
    for o in ((598, 385), (710, 398)):
        L.linea(mg[0], mg[1], o[0], o[1], color=VERDE, sw=2, dash="2 5")
    L.circulo(*mg, 17, fill=VERDE)
    L.texto(mg[0], mg[1] + 5, "Mg", size=13, weight=600, color="#ffffff", anchor="middle")
    L.texto(mg[0] + 26, mg[1] + 32, "Mg²⁺ neutraliza", size=13.5, color=TINTA2)
    L.texto(mg[0] + 26, mg[1] + 49, "las cargas del fosfato", size=13.5, color=TINTA2)

    # --- Lys169
    lys = (560, 180)
    L.linea(lys[0], lys[1] + 17, 560, 262 - 12, color=AZUL, sw=2, dash="2 5")
    atomo(L, *lys, 17, "N⁺", color=AZUL, size=13)
    enlace(L, lys[0] - 17, lys[1], 470, lys[1])
    L.texto(460, lys[1] - 4, "Lys169", size=16, weight=600, anchor="end")
    L.texto(460, lys[1] + 15, "carga (+)", size=14, color=TINTA2, anchor="end")

    # --- flechas de lo que va a pasar
    L.flecha(410, y - 22, pg[0] - 20, y - 26, color=NARANJA, sw=2.4, curvatura=0.28)
    L.flecha(hO6[0] - 12, hO6[1] + 10, od1[0] - 16, od1[1] - 8, color=AZUL, sw=2.4, curvatura=0.35)

    # leyenda
    ly = 600
    L.linea(600, ly - 5, 630, ly - 5, color=TINTA2, sw=2.4, dash="3 6")
    L.texto(638, ly, "distancia medida", size=13.5, color=TINTA2)
    L.flecha(770, ly - 5, 805, ly - 5, color=TINTA2, sw=2.4)
    L.texto(813, ly, "lo que ocurrirá", size=13.5, color=TINTA2)
    return L


def ilustracion_ampnp_atp():
    L = Lienzo(360, etiqueta="Para pensar", titulo="El truco del cristalógrafo: un ATP que no puede reaccionar",
               subtitulo="Cambiar un solo átomo (O3B → N3B) congela el complejo justo antes de la química")
    for k, (x0, puente, titulo, sub, bloqueado) in enumerate([
            (60, "N", "AMP-PNP (en el cristal)", "el puente N no se rompe: no reacciona", True),
            (540, "O", "ATP (en nuestra simulación)", "el puente O sí se rompe: reacciona", False)]):
        y = 250
        L.rect(x0, 150, 400, 185, rx=16, fill="#ffffff", stroke=REJILLA, sw=2)
        L.texto(x0 + 24, 185, titulo, size=17, weight=600)
        L.texto(x0 + 24, 207, sub, size=14, color=TINTA2)
        pb, br, pg = (x0 + 90, y + 20), (x0 + 175, y + 20), (x0 + 260, y + 20)
        enlace(L, x0 + 30, y + 20, *pb)
        enlace(L, *pb, *br)
        enlace(L, *br, *pg)
        atomo(L, *pb, 22, "Pβ", color=FOSFATO, relleno=0.45, size=16)
        atomo(L, *pg, 22, "Pγ", color=FOSFATO, relleno=0.45, size=16)
        col = VIOLETA if bloqueado else FOSFATO
        atomo(L, *br, 16, puente, color=col, relleno=0.6, size=15)
        if bloqueado:
            candado(L, x0 + 345, y + 14, 1.2)
        else:
            L.linea(br[0] + 22, br[1] - 22, br[0] + 40, br[1] + 22, color=ROJO, sw=3)  # "tijeras": enlace que se corta
            L.flecha(pg[0] + 30, pg[1], pg[0] + 90, pg[1] - 30, color=NARANJA, sw=2.4)
            L.texto(pg[0] + 60, pg[1] + 28, "a la glucosa", size=13, color=TINTA2, anchor="middle")
    L.flecha(465, 262, 530, 262, color=TINTA, sw=3)
    L.texto(497, 247, "N → O", size=13, weight=600, anchor="middle")
    return L


# ---------------------------------------------------------------- sección 4
def ilustracion_preparacion_pasos():
    L = Lienzo(470, etiqueta="Modelo", titulo="De la foto del cristal a un modelo que se puede simular",
               subtitulo="Cuatro arreglos antes de encender la simulación (script 01_preparar_sistema.py)")
    xs = [40, 275, 510, 745]
    w, top = 215, 150
    textos = [
        ("Añadir hidrógenos", ["los rayos X no ven", "los átomos de H"]),
        ("Completar bucles", ["3 trozos desordenados", "se modelan (PDBFixer)"]),
        ("Poner parámetros", ["resortes y cargas:", "ff14SB, GAFF2, Li–Merz"]),
        ("Agua e iones", ["caja TIP3P de 10 Å y", "21 Na⁺ (puntos violeta)"]),
    ]
    for i, (x, (tit, lineas)) in enumerate(zip(xs, textos)):
        L.rect(x, top, w, 290, rx=16, fill="#ffffff", stroke=REJILLA, sw=2)
        L.numero(x + 28, top + 30, i + 1)
        L.texto(x + 52, top + 36, tit, size=17, weight=600)
        L.texto(x + 20, top + 240, lineas, size=14.5, color=TINTA2)
        cx, cy = x + w / 2, top + 135
        if i == 0:
            cad = [(cx - 70, cy + 10), (cx - 25, cy - 18), (cx + 20, cy + 10), (cx + 65, cy - 18)]
            for a, b in zip(cad, cad[1:]):
                enlace(L, *a, *b)
            for k, (ax, ay) in enumerate(cad):
                hy = ay + (30 if k % 2 == 0 else -30)
                enlace(L, ax, ay, ax, hy, sw=2.2, color=AGUA)
                L.circulo(ax, hy, 8, fill="#ffffff", stroke=AGUA, sw=2.4)
                L.texto(ax, hy + 4.5, "H", size=10, weight=600, anchor="middle", color=TINTA2)
                atomo(L, ax, ay, 13, color=TINTA3, relleno=0.6)
            L.texto(cx, cy + 75, "+ H", size=14, weight=600, anchor="middle", color=TINTA2)
        elif i == 1:
            cy += 18
            L.curva([(cx - 85, cy + 40), (cx - 60, cy - 10), (cx - 30, cy - 35)], stroke=AZUL, sw=6)
            L.curva([(cx + 30, cy - 35), (cx + 60, cy - 10), (cx + 85, cy + 40)], stroke=AZUL, sw=6)
            L.curva([(cx - 30, cy - 35), (cx - 5, cy - 70), (cx + 5, cy - 70), (cx + 30, cy - 35)],
                    stroke=AGUA, sw=6, dash="2 10")
            L.texto(cx, cy - 84, "bucle añadido", size=13, color=TINTA2, anchor="middle")
            L.texto(cx, cy + 58, "cadena del cristal", size=13, color=TINTA3, anchor="middle")
        elif i == 2:
            a, b, c = (cx - 60, cy + 20), (cx, cy - 25), (cx + 60, cy + 20)
            resorte(L, *a, *b, n=5)
            resorte(L, *b, *c, n=5)
            atomo(L, *a, 17, "δ+", color=AZUL, size=13, relleno=0.7)
            atomo(L, *b, 19, "δ−", color=ROJO, size=13, relleno=0.7)
            atomo(L, *c, 17, "δ+", color=AZUL, size=13, relleno=0.7)
            L.texto(cx, cy + 70, "enlaces = resortes", size=13, color=TINTA3, anchor="middle")
        else:
            L.rect(cx - 88, cy - 75, 176, 140, rx=10, fill=tinte(AZUL, 0.93), stroke=AZUL, sw=2, dash="6 5")
            L.elipse(cx, cy - 5, 40, 30, fill=tinte(AZUL, 0.6), stroke=AZUL, sw=2)
            for k, (px, py) in enumerate([(-65, -55), (-20, -60), (35, -58), (68, -40), (-70, -10), (70, 5),
                                          (-60, 40), (-15, 45), (30, 48), (65, 42), (-35, 20), (48, -30)]):
                agua_molecula(L, cx + px, cy + py, 0.85, rot=k * 47)
            for px, py in [(-52, -30), (58, 28), (8, 30)]:
                L.circulo(cx + px, cy + py, 7, fill=VIOLETA)
        if i < 3:
            L.flecha(x + w + 3, top + 145, x + w + 17, top + 145, color=TINTA3, sw=2.5)
    return L


# ---------------------------------------------------------------- sección 5
def _mini_molecula(L, cx, cy, fase, escala=1.0):
    """Pγ del ATP y O6 de la glucosa que vibran; ``fase`` cambia ligeramente las posiciones."""
    dx = 7 * math.sin(fase) * escala
    dy = 5 * math.cos(1.7 * fase) * escala
    sep = 27 * escala
    p = (cx - sep + dx, cy + dy)
    o = (cx + sep - dx * 0.6, cy - dy)
    rp, ro, rt, lr = 13 * escala, 11 * escala, 7 * escala, 20 * escala
    tp = (p[0] - rp - lr, p[1] + 9 * escala)
    to = (o[0] + ro + lr, o[1] + 9 * escala)
    resorte(L, tp[0] + rt, tp[1], p[0] - rp, p[1], n=3, amp=4 * escala, sw=1.8)
    resorte(L, o[0] + ro, o[1], to[0] - rt, to[1], n=3, amp=4 * escala, sw=1.8)
    atomo(L, *tp, rt, color=FOSFATO, relleno=0.6)
    L.circulo(*to, rt, fill=tinte(NARANJA, 0.6), stroke=NARANJA, sw=2)
    L.linea(p[0] + rp, p[1], o[0] - ro, o[1], color=TINTA3, sw=1.6, dash="2 4")
    atomo(L, *p, rp, "P", color=FOSFATO, relleno=0.45, size=12 * escala)
    atomo(L, *o, ro, "O", color=NARANJA, relleno=0.6, size=11 * escala)


def ilustracion_pelicula_md():
    L = Lienzo(460, etiqueta="Analogía", titulo="El cristal es una foto; la dinámica molecular es la película",
               subtitulo="Cada fotograma avanza 2 femtosegundos; un nanosegundo son 500 000 fotogramas")
    # foto
    fx, fy, fw, fh = 50, 205, 190, 150
    L.rect(fx, fy, fw, fh, rx=8, fill="#ffffff", stroke=TINTA2, sw=2.5)
    L.rect(fx + 10, fy + 10, fw - 20, fh - 45, rx=4, fill=tinte(AZUL, 0.92))
    _mini_molecula(L, fx + fw / 2, fy + 62, 0.0, escala=1.1)
    L.texto(fx + fw / 2, fy + fh - 13, "3FGU: una sola pose", size=13.5, color=TINTA2, anchor="middle")
    L.texto(fx + fw / 2, fy - 16, "Foto (cristal)", size=17, weight=600, anchor="middle")
    L.flecha(fx + fw + 14, fy + fh / 2, fx + fw + 64, fy + fh / 2, color=TINTA, sw=3)
    # tira de película
    sx, sy, n, cw = 320, 190, 4, 150
    ancho = n * cw + 20
    L.rect(sx, sy, ancho, 180, rx=10, fill=TINTA2)
    for k in range(int(ancho / 22)):
        L.rect(sx + 8 + k * 22, sy + 8, 11, 9, rx=2, fill=SUPER)
        L.rect(sx + 8 + k * 22, sy + 180 - 17, 11, 9, rx=2, fill=SUPER)
    for i in range(n):
        x = sx + 10 + i * cw
        L.rect(x + 5, sy + 26, cw - 10, 128, rx=4, fill=tinte(AZUL, 0.92))
        _mini_molecula(L, x + cw / 2, sy + 90, i * 1.3, escala=1.0)
        L.texto(x + cw / 2, sy + 205, f"t = {2 * i} fs", size=14, color=TINTA2, anchor="middle")
    L.texto(sx + ancho / 2, sy - 16, "Película (dinámica molecular)", size=17, weight=600, anchor="middle")
    L.texto(sx + ancho + 8, sy + 95, "…", size=26, color=TINTA2)
    L.texto(500, 440, "La distancia P···O cambia en cada fotograma: los reactivos «vibran», pero siguen apuntándose.",
            size=14.5, color=TINTA2, anchor="middle")
    return L


def ilustracion_newton_ciclo():
    L = Lienzo(545, etiqueta="La ecuación, en palabras",
               titulo="Cómo se hace cada fotograma: las leyes de Newton, una y otra vez",
               subtitulo="Átomos = bolitas con masa; enlaces = resortes; además cargas que se atraen o repelen")
    # modelo de bolitas y resortes (izquierda)
    a, b, c, d = (110, 330), (190, 275), (275, 320), (215, 405)
    resorte(L, *a, *b, n=5)
    resorte(L, *b, *c, n=5)
    resorte(L, *b, *d, n=5)
    atomo(L, *a, 18, "+", color=AZUL, size=18)
    atomo(L, *b, 24, "", color=TINTA3, relleno=0.6)
    atomo(L, *c, 18, "−", color=ROJO, size=20)
    atomo(L, *d, 16, "", color=TINTA3, relleno=0.75)
    L.flecha(b[0] + 5, b[1] - 30, b[0] + 45, b[1] - 75, color=NARANJA, sw=3)
    L.mate(b[0] + 52, b[1] - 78, "F", size=22)
    L.texto(190, 455, ["bolita = átomo (masa m)", "resorte = enlace"], size=14, color=TINTA2, anchor="middle")

    # ciclo (derecha)
    cx, cy, R = 650, 335, 118
    pasos = [
        (-90, "1. Posiciones", "dónde está cada átomo", AZUL),
        (0, "2. Fuerzas", "resortes + cargas + choques", NARANJA),
        (90, "3. Aceleración", "a = F / m", VIOLETA),
        (180, "4. Mover 2 fs", "posición nueva", AGUA),
    ]
    L.circulo(cx, cy, R, stroke=REJILLA, sw=3)
    for ang in (-45, 45, 135, 225):
        t = math.radians(ang)
        t2 = math.radians(ang + 16)
        x1, y1 = cx + R * math.cos(t - 0.14), cy + R * math.sin(t - 0.14)
        x2, y2 = cx + R * math.cos(t2), cy + R * math.sin(t2)
        L.flecha(x1, y1, x2, y2, color=TINTA3, sw=3, curvatura=-0.08)
    for ang, tit, sub, col in pasos:
        t = math.radians(ang)
        x, y = cx + R * math.cos(t), cy + R * math.sin(t)
        L.circulo(x, y, 13, fill=col)
        ax = {-90: "middle", 90: "middle", 0: "start", 180: "end"}[ang]
        ox = {0: 22, 180: -22}.get(ang, 0)
        oy = {-90: -40, 90: 42}.get(ang, -4)
        L.texto(x + ox, y + oy, tit, size=16, weight=600, anchor=ax)
        L.texto(x + ox, y + oy + 19, sub, size=13.5, color=TINTA2, anchor=ax)
    L.texto(cx, cy - 4, "×500 000", size=20, weight=600, anchor="middle")
    L.texto(cx, cy + 18, "vueltas = 1 ns", size=14, color=TINTA2, anchor="middle")
    return L


def ilustracion_ataque_cercano():
    L = Lienzo(440, etiqueta="Modelo", titulo="Conformación de ataque cercano: cuando los reactivos se apuntan",
               subtitulo="Solo desde esas posturas puede ocurrir la química; la MD mide qué fracción del tiempo pasa en ellas")
    for x0, cerca, tit in ((40, True, "Ataque cercano (NAC)"), (520, False, "No reactiva")):
        L.rect(x0, 160, 440, 250, rx=16, fill="#ffffff", stroke=REJILLA, sw=2)
        L.texto(x0 + 60, 197, tit, size=18, weight=600)
        (check if cerca else cruz)(L, x0 + 35, 191)
        y = 300
        pg = (x0 + 150, y)
        o6 = (x0 + (255 if cerca else 345), y + (0 if cerca else 35))
        # ATP a la izquierda
        enlace(L, x0 + 60, y, *pg)
        atomo(L, x0 + 70, y, 13, color=FOSFATO, relleno=0.6)
        atomo(L, *pg, 22, "Pγ", color=FOSFATO, relleno=0.45, size=15)
        # glucosa a la derecha
        hexagono(L, o6[0] + 62, o6[1] + 8, 26, NARANJA, rot=30, relleno=0.7)
        enlace(L, o6[0] + 14, o6[1] + 3, o6[0] + 38, o6[1] + 4)
        atomo(L, *o6, 15, "O6", color=NARANJA, size=12)
        L.linea(pg[0] + 24, pg[1] + (0 if cerca else 6), o6[0] - 17, o6[1] - (0 if cerca else 4),
                color=NARANJA if cerca else TINTA3, sw=2.4, dash="3 6")
        mx = (pg[0] + o6[0]) / 2
        L.texto(mx, y - (22 if cerca else 8), "d < 3.5 Å" if cerca else "d > 3.5 Å", size=16, weight=600,
                anchor="middle")
        L.texto(x0 + 30, 372, "el O6 apunta al fósforo y está cerca" if cerca else "demasiado lejos o mal orientado",
                size=14.5, color=TINTA2)
        L.texto(x0 + 30, 392, "→ la reacción es posible" if cerca else "→ hay que esperar a que se acerque",
                size=14.5, color=TINTA2)
    return L
