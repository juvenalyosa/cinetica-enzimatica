"""Portada y secciones 1, 2 y 17: recorrido del curso, la enzima como guía, la reacción de la
glucoquinasa, su papel de sensor, la velocidad inicial, el ensayo acoplado y el mapa del curso."""
from __future__ import annotations

import math

from .lienzo import (AGUA, AMARILLO, AZUL, EJE, MAGENTA, NARANJA, REJILLA, ROJO, TINTA, TINTA2, TINTA3,
                     VERDE, VIOLETA, Lienzo, muestrear, perfil_barrera, tinte)


# ---------------------------------------------------------------- helpers locales
def _hexagono(L, cx, cy, r, *, fill, stroke, sw=2.5):
    pts = [(cx + r * math.cos(math.radians(30 + 60 * k)), cy + r * math.sin(math.radians(30 + 60 * k))) for k in range(6)]
    L.poligono(pts, fill=fill, stroke=stroke, sw=sw)
    return pts


def _fosfato(L, cx, cy, letra, *, color=NARANJA, r=19):
    L.circulo(cx, cy, r, fill=tinte(color, 0.75), stroke=color, sw=2.2)
    L.texto(cx, cy + 6, "P", size=17, weight=600, anchor="middle")
    if letra:
        L.texto(cx, cy + r + 20, letra, size=14, color=TINTA2, anchor="middle", italic=True)


def _caja(L, x, y, w, h, titulo, lineas=None, *, color=AZUL, size=17):
    """Tarjeta blanca con barra de color arriba y texto en tinta."""
    L.rect(x, y, w, h, rx=14, fill="#ffffff", stroke=REJILLA, sw=2)
    L.rect(x, y, w, 7, rx=3.5, fill=color)
    L.texto(x + w / 2, y + 36, titulo, size=size, weight=600, anchor="middle")
    if lineas:
        L.texto(x + w / 2, y + 60, lineas, size=14, color=TINTA2, anchor="middle", lineas=1.3)


# ---------------------------------------------------------------- portada
def ilustracion_portada_banner():
    L = Lienzo(535, etiqueta="El curso en una imagen",
               titulo="De los átomos a la velocidad: la glucoquinasa, paso a paso",
               subtitulo="Una reacción real, cuatro miradas y los números medidos en el laboratorio")
    # --- reacción
    y = 215
    # enzima (almeja) en el centro
    cx = 500
    L.camino(f"M {cx - 70} {y + 30} Q {cx - 70} {y - 55} {cx} {y - 55} Q {cx + 70} {y - 55} {cx + 70} {y + 30} "
             f"L {cx + 40} {y + 30} Q {cx + 40} {y - 20} {cx} {y - 20} Q {cx - 40} {y - 20} {cx - 40} {y + 30} Z",
             fill=tinte(AZUL, 0.8), stroke=AZUL, sw=2.5)
    L.texto(cx, y + 58, "glucoquinasa", size=15, weight=600, anchor="middle")
    L.texto(cx, y + 77, "(la enzima)", size=13, color=TINTA2, anchor="middle")
    # sustratos
    _hexagono(L, 120, y, 30, fill=tinte(NARANJA, 0.8), stroke=NARANJA)
    L.texto(120, y + 58, "glucosa", size=15, weight=600, anchor="middle")
    L.texto(185, y + 7, "+", size=26, color=TINTA2, anchor="middle")
    for k in range(3):
        _fosfato(L, 235 + k * 40, y, None, r=16)
    L.texto(275, y + 58, "ATP", size=15, weight=600, anchor="middle")
    L.flecha(345, y, cx - 85, y, color=TINTA2)
    L.flecha(cx + 85, y, 655, y, color=TINTA2)
    # productos
    _hexagono(L, 705, y, 30, fill=tinte(AGUA, 0.8), stroke=AGUA)
    _fosfato(L, 752, y - 26, None, color=AGUA, r=14)
    L.texto(715, y + 58, "glucosa-6-fosfato", size=15, weight=600, anchor="middle")
    L.texto(800, y + 7, "+", size=26, color=TINTA2, anchor="middle")
    for k in range(2):
        _fosfato(L, 845 + k * 38, y, None, color=AGUA, r=16)
    L.texto(864, y + 58, "ADP", size=15, weight=600, anchor="middle")

    # --- recorrido del curso
    y2 = 395
    pasos = [("Estructura", "la foto 3D (3FGU)", AZUL, "§3–4"),
             ("Dinámica molecular", "la película de la enzima", AGUA, "§5"),
             ("Química cuántica", "QM/MM: la colina", NARANJA, "§6–12"),
             ("Cinética", "Michaelis–Menten, Hill, inhibidores", MAGENTA, "§13–16")]
    w, gap, x0 = 205, 33, 40
    L.linea(x0 + 20, y2 - 38, 960, y2 - 38, color=REJILLA, sw=2)
    for i, (t, sub, col, sec) in enumerate(pasos):
        x = x0 + i * (w + gap)
        L.circulo(x + w / 2, y2 - 38, 9, fill=col)
        L.rect(x, y2 - 10, w, 118, rx=14, fill="#ffffff", stroke=REJILLA, sw=2)
        L.texto(x + 18, y2 + 22, sec, size=13, weight=600, color=TINTA3)
        L.texto(x + 18, y2 + 50, t, size=18, weight=600)
        L.texto(x + 18, y2 + 76, sub if len(sub) < 26 else [sub[:sub.index(",") + 1], sub[sub.index(",") + 2:]],
                size=14, color=TINTA2)
        if i < 3:
            L.flecha(x + w + 5, y2 + 48, x + w + gap - 5, y2 + 48, color=TINTA3, sw=2)
    return L


# ---------------------------------------------------------------- sección 1
def ilustracion_guia_paso_bajo():
    L = Lienzo(590, etiqueta="Analogía", titulo="La enzima es un guía que conoce el paso más bajo",
               subtitulo="Mismo punto de partida, mismo destino: solo cambia la altura del camino")
    x0, x1 = 70, 930
    base, fin = 480, 520
    f_sin = perfil_barrera(x0, x1, base, 200, fin, centro=0.5, ancho=0.16)
    f_con = perfil_barrera(x0, x1, base, 350, fin, centro=0.5, ancho=0.19)
    sin = muestrear(f_sin, x0, x1, 160)
    con = muestrear(f_con, x0, x1, 160)
    L.curva(sin, stroke=TINTA3, sw=3, dash="8 7")
    L.curva(con, stroke=AZUL, sw=4)
    xm = (x0 + x1) / 2
    # cotas a la izquierda, sobre la llanura
    L.linea(95, 200, xm - 10, 200, color=EJE, sw=1.5, dash="4 5")
    L.linea(95, 350, xm - 10, 350, color=EJE, sw=1.5, dash="4 5")
    L.flecha(110, base - 4, 110, 206, color=TINTA3, sw=2.2, doble=True)
    L.flecha(190, base - 4, 190, 356, color=AZUL, sw=2.2, doble=True)
    L.mate(124, 290, "ΔG^{‡}", size=22)
    L.texto(124, 310, "sin enzima", size=14, color=TINTA2)
    L.mate(204, 418, "ΔG^{‡}", size=22)
    L.texto(204, 438, "con enzima", size=14, color=TINTA2)
    # rótulos de los caminos
    L.linea(xm + 30, 212, xm + 80, 190, color=TINTA3, sw=1.5)
    L.texto(xm + 86, 186, ["paso alto (sin enzima):", "casi nadie lo cruza"], size=15, color=TINTA2)
    L.linea(xm + 100, 380, xm + 150, 330, color=AZUL, sw=1.5)
    L.texto(xm + 156, 318, ["paso bajo (con enzima):", "el camino del guía"], size=15, weight=600)
    # valles
    L.circulo(290, base - 14, 14, fill=AZUL)
    L.circulo(xm + 400, f_con(xm + 400) - 14, 14, fill=AGUA)
    L.texto(290, 560, "reactivos", size=17, weight=600, anchor="middle")
    L.texto(xm + 400, 560, "productos", size=17, weight=600, anchor="middle")
    return L


def ilustracion_reaccion_glucoquinasa():
    L = Lienzo(625, etiqueta="La reacción", titulo="Qué hace la glucoquinasa: mover un fosfato",
               subtitulo="El fosfato γ del ATP salta al oxígeno O6 de la glucosa; Asp205 se lleva el protón")
    # --- ATP a la izquierda
    y = 330
    # adenosina
    L.rect(60, y - 32, 110, 64, rx=14, fill=tinte(TINTA3, 0.85), stroke=TINTA3, sw=2)
    L.texto(115, y + 6, "adenosina", size=15, weight=600, anchor="middle")
    xs = [215, 285, 355]
    L.linea(170, y, xs[0] - 19, y, color=TINTA2, sw=2.5)
    L.linea(xs[0] + 19, y, xs[1] - 19, y, color=TINTA2, sw=2.5)
    L.linea(xs[1] + 19, y, xs[2] - 19, y, color=TINTA2, sw=2.5, dash="5 5")
    for x, l in zip(xs, ["α", "β", "γ"]):
        _fosfato(L, x, y, l)
    L.texto(115, y + 70, "ATP", size=19, weight=600, anchor="middle")
    # enlace que se rompe
    L.texto(320, y - 30, "se rompe", size=13, color=TINTA2, anchor="middle")
    # Mg
    L.circulo(320, y + 70, 20, fill=tinte(VERDE, 0.7), stroke=VERDE, sw=2.2)
    L.texto(320, y + 76, "Mg²⁺", size=13, weight=600, anchor="middle")
    L.linea(305, y + 55, 290, y + 20, color=VERDE, sw=1.8, dash="3 4")
    L.linea(335, y + 55, 350, y + 20, color=VERDE, sw=1.8, dash="3 4")
    L.texto(320, y + 118, ["sujeta y calma", "las cargas del fosfato"], size=13, color=TINTA2, anchor="middle")
    # --- glucosa a la derecha
    gx, gy = 700, y
    _hexagono(L, gx, gy, 48, fill=tinte(NARANJA, 0.82), stroke=NARANJA)
    L.texto(gx, gy + 6, "glucosa", size=15, weight=600, anchor="middle")
    # brazo O6-H
    L.linea(gx - 42, gy - 24, gx - 90, gy - 70, color=TINTA2, sw=2.5)
    L.circulo(gx - 100, gy - 80, 16, fill="#ffffff", stroke=ROJO, sw=2.2)
    L.texto(gx - 100, gy - 74, "O6", size=13, weight=600, anchor="middle")
    L.linea(gx - 114, gy - 90, gx - 135, gy - 110, color=TINTA2, sw=2)
    L.circulo(gx - 141, gy - 116, 10, fill="#ffffff", stroke=TINTA3, sw=2)
    L.texto(gx - 141, gy - 111, "H", size=12, weight=600, anchor="middle")
    # salto del fosfato γ → O6
    L.flecha(xs[2] + 14, y - 22, gx - 118, gy - 86, color=NARANJA, sw=3, curvatura=-0.22)
    L.texto(300, y - 110, "el fosfato γ salta", size=16, weight=600, anchor="middle")
    L.texto(300, y - 90, "del ATP al O6", size=14, color=TINTA2, anchor="middle")
    # Asp205 toma el protón
    ax, ay = 860, 175
    L.rect(ax - 75, ay - 28, 150, 56, rx=12, fill=tinte(AZUL, 0.85), stroke=AZUL, sw=2)
    L.texto(ax, ay + 6, "Asp205 (–)", size=15, weight=600, anchor="middle")
    L.flecha(gx - 128, gy - 124, ax - 80, ay - 2, color=AZUL, sw=2.2, dash="6 5", curvatura=-0.18)
    L.texto(ax, ay + 50, ["la base: se lleva", "el protón del O6"], size=13, color=TINTA2, anchor="middle")
    # Lys169, junto al camino del fosfato
    lx, ly = 510, y + 105
    L.linea(lx, ly - 28, 478, y - 82, color=VIOLETA, sw=1.8, dash="3 4")
    L.rect(lx - 70, ly - 28, 140, 56, rx=12, fill=tinte(VIOLETA, 0.88), stroke=VIOLETA, sw=2)
    L.texto(lx, ly + 6, "Lys169 (+)", size=15, weight=600, anchor="middle")
    L.texto(lx, ly + 50, ["estabiliza el fosfato", "que viaja"], size=13, color=TINTA2, anchor="middle")
    # resultado
    L.nota(60, 540, 880, "Resultado: glucosa-6-fosfato + ADP. Con el fosfato pegado, la glucosa ya no puede salir de la célula.",
           color=AGUA, size=15)
    return L


def ilustracion_sensor_glucosa():
    L = Lienzo(470, etiqueta="Por qué importa", titulo="Un sensor de glucosa que decide cuánta insulina liberar",
               subtitulo="La glucoquinasa trabaja más cuanta más glucosa hay en la sangre")
    y = 250
    # sangre
    L.rect(50, y - 70, 190, 150, rx=16, fill=tinte(ROJO, 0.9), stroke=ROJO, sw=2)
    L.texto(145, y - 38, "sangre", size=16, weight=600, anchor="middle")
    for i, (dx, dy) in enumerate([(-50, 0), (0, -5), (48, 4), (-30, 40), (22, 38), (60, 45)]):
        _hexagono(L, 145 + dx, y + 10 + dy, 13, fill=tinte(NARANJA, 0.75), stroke=NARANJA, sw=1.8)
    L.texto(145, y + 110, "glucosa", size=14, color=TINTA2, anchor="middle")
    L.flecha(250, y, 330, y, color=TINTA2)
    # célula beta con GK
    L.elipse(470, y + 5, 125, 95, fill=tinte(AZUL, 0.9), stroke=AZUL, sw=2)
    L.texto(470, y - 55, "célula β del páncreas", size=15, weight=600, anchor="middle")
    L.rect(400, y - 25, 140, 50, rx=25, fill="#ffffff", stroke=AZUL, sw=2)
    L.texto(470, y + 6, "glucoquinasa", size=15, weight=600, anchor="middle")
    L.texto(470, y + 55, "más glucosa, más actividad", size=13, color=TINTA2, anchor="middle")
    L.flecha(600, y, 680, y, color=TINTA2)
    # insulina
    L.rect(690, y - 70, 250, 150, rx=16, fill=tinte(AGUA, 0.88), stroke=AGUA, sw=2)
    L.texto(815, y - 38, "insulina liberada", size=16, weight=600, anchor="middle")
    for i in range(7):
        L.circulo(735 + i * 30, y + 10 + (8 if i % 2 else -8), 9, fill=AGUA)
    L.texto(815, y + 60, "más glucosa → más insulina", size=13, color=TINTA2, anchor="middle")
    L.nota(50, 375, 890, ["Si una mutación frena a la glucoquinasa, el páncreas «ve» menos glucosa de la que hay: diabetes GCK-MODY.",
                          "Si la acelera, libera insulina de más: hipoglucemia congénita."], color=NARANJA, size=14.5)
    return L


# ---------------------------------------------------------------- sección 2
def ilustracion_velocidad_inicial():
    L = Lienzo(560, etiqueta="Analogía", titulo="Velocidad = cuánto producto aparece por segundo",
               subtitulo="Como contar las monedas que caen en un frasco: al principio el ritmo es constante")
    # tres instantáneas
    xs = [95, 245, 395]
    tiempos = ["t = 0 s", "t = 5 s", "t = 10 s"]
    cuentas = [0, 5, 10]
    ybase = 440
    for x, t, n in zip(xs, tiempos, cuentas):
        L.camino(f"M {x - 55} {ybase - 190} L {x - 45} {ybase} L {x + 45} {ybase} L {x + 55} {ybase - 190}",
                 stroke=TINTA2, sw=2.5)
        for k in range(n):
            fila, col = divmod(k, 4)
            L.circulo(x - 30 + col * 20 + (10 if fila % 2 else 0), ybase - 14 - fila * 19, 8.5, fill=MAGENTA)
        L.texto(x, ybase + 32, t, size=15, weight=600, anchor="middle")
        L.texto(x, ybase + 52, f"{n} µM de producto", size=13, color=TINTA2, anchor="middle")
    L.texto(245, 205, "1 µM por segundo ⇒ v₀ = 1 µM/s", size=16, weight=600, anchor="middle")
    L.linea(470, 190, 470, 500, color=REJILLA, sw=2)
    # curva de progreso
    gx0, gy0, gw, gh = 540, 470, 400, 260
    L.linea(gx0, gy0, gx0 + gw, gy0, color=EJE, sw=2)
    L.linea(gx0, gy0, gx0, gy0 - gh, color=EJE, sw=2)
    L.texto(gx0 + gw, gy0 + 28, "tiempo →", size=14, color=TINTA2, anchor="end")
    L.texto(gx0, gy0 - gh - 14, "producto ↑", size=14, color=TINTA2)
    f = lambda x: gy0 - gh * 0.88 * (1 - math.exp(-(x - gx0) / 150))
    L.curva(muestrear(f, gx0, gx0 + gw, 80), stroke=MAGENTA, sw=3.5)
    # tangente inicial
    pend = gh * 0.88 / 150
    L.linea(gx0, gy0, gx0 + 140, gy0 - pend * 140, color=TINTA, sw=2, dash="7 6")
    L.texto(gx0 + 150, gy0 - pend * 140 - 8, "pendiente inicial", size=15, weight=600)
    L.mate(gx0 + 150, gy0 - pend * 140 + 18, "= v_{0}", size=22)
    L.texto(gx0 + gw - 10, f(gx0 + gw) + 95, ["después la curva se dobla:", "el sustrato se agota"], size=13.5,
            color=TINTA2, anchor="end")
    return L


def ilustracion_ensayo_acoplado():
    L = Lienzo(560, etiqueta="Cómo se mide", titulo="El ensayo acoplado: convertir G6P en luz absorbida",
               subtitulo="Una segunda enzima produce NADPH, que absorbe luz ultravioleta de 340 nm")
    y = 200
    # cadena de reacciones
    items = [(95, "glucosa + ATP", NARANJA), (400, "G6P", AGUA), (705, "6-fosfoglucono-", TINTA3)]
    for x, t, c in items:
        L.rect(x - 85, y - 30, 170, 60, rx=14, fill=tinte(c, 0.85), stroke=c, sw=2)
    L.texto(95, y + 6, "glucosa + ATP", size=15, weight=600, anchor="middle")
    L.texto(400, y + 6, "G6P", size=16, weight=600, anchor="middle")
    L.texto(705, y - 3, "6-fosfoglucono-", size=14, weight=600, anchor="middle")
    L.texto(705, y + 15, "δ-lactona", size=14, weight=600, anchor="middle")
    L.flecha(185, y, 310, y, color=TINTA2)
    L.texto(248, y - 16, "glucoquinasa", size=14, weight=600, anchor="middle")
    L.texto(248, y + 26, "(la que medimos)", size=12.5, color=TINTA2, anchor="middle")
    L.flecha(490, y, 615, y, color=TINTA2)
    L.texto(553, y - 16, "G6PDH", size=14, weight=600, anchor="middle")
    L.texto(553, y + 26, "(en exceso)", size=12.5, color=TINTA2, anchor="middle")
    # NADP+ -> NADPH
    L.camino(f"M 510 {y + 75} Q 553 {y + 40} 596 {y + 75}", stroke=VIOLETA, sw=2.2, flecha=True)
    L.texto(500, y + 95, "NADP⁺", size=14, anchor="middle")
    L.rect(575, y + 75, 90, 32, rx=16, fill=tinte(VIOLETA, 0.8), stroke=VIOLETA, sw=2)
    L.texto(620, y + 97, "NADPH", size=14, weight=600, anchor="middle")
    L.texto(690, y + 87, ["1 G6P ⇒ 1 NADPH:", "contar NADPH es contar producto"], size=13.5, color=TINTA2)
    # espectrofotómetro
    sy = 425
    L.rect(60, sy - 40, 90, 80, rx=12, fill=tinte(AMARILLO, 0.8), stroke=AMARILLO, sw=2)
    L.texto(105, sy + 6, "lámpara", size=14, weight=600, anchor="middle")
    L.texto(105, sy + 62, "luz de 340 nm", size=13, color=TINTA2, anchor="middle")
    L.linea(150, sy, 270, sy, color=VIOLETA, sw=6, opacity=.85)
    L.rect(270, sy - 55, 60, 110, rx=6, fill="#ffffff", stroke=TINTA2, sw=2.5)
    L.rect(274, sy - 20, 52, 71, rx=3, fill=tinte(VIOLETA, 0.7))
    L.texto(300, sy + 78, "cubeta", size=13, color=TINTA2, anchor="middle")
    L.linea(330, sy, 450, sy, color=VIOLETA, sw=3, opacity=.5)
    L.rect(450, sy - 40, 90, 80, rx=12, fill=tinte(TINTA3, 0.85), stroke=TINTA3, sw=2)
    L.texto(495, sy + 6, "detector", size=14, weight=600, anchor="middle")
    L.texto(390, sy - 16, "sale menos luz", size=12.5, color=TINTA2, anchor="middle")
    L.flecha(550, sy, 610, sy, color=TINTA3, sw=2)
    # registro
    gx0, gy0 = 630, sy + 60
    L.linea(gx0, gy0, gx0 + 300, gy0, color=EJE, sw=2)
    L.linea(gx0, gy0, gx0, gy0 - 120, color=EJE, sw=2)
    L.linea(gx0, gy0 - 10, gx0 + 260, gy0 - 110, color=VIOLETA, sw=3)
    L.texto(gx0 + 300, gy0 + 22, "tiempo", size=13, color=TINTA2, anchor="end")
    L.texto(gx0 + 8, gy0 - 128, "absorbancia a 340 nm", size=13, color=TINTA2)
    L.texto(gx0 + 140, gy0 - 30, "pendiente ∝ v₀", size=14, weight=600)
    return L


# ---------------------------------------------------------------- sección 17
def ilustracion_mapa_curso():
    L = Lienzo(640, etiqueta="Resumen", titulo="El mapa del curso: de los átomos a la curva de velocidad",
               subtitulo="Cada flecha es una idea que conecta una sección con la siguiente")
    W, H = 200, 88
    nodos = {
        "est": (40, 170, "Estructura", ["3FGU: reactivos", "alineados"], AZUL, "§3–4"),
        "md": (280, 170, "Dinámica molecular", ["los mantiene", "alineados (NAC)"], AGUA, "§5"),
        "qm": (520, 170, "QM/MM", ["ve la química:", "R → TS → P"], NARANJA, "§7–9"),
        "dg": (760, 170, "ΔE‡ → ΔG‡", ["vibraciones, conf.,", "método (± error)"], NARANJA, "§10–11"),
        "ey": (760, 330, "Eyring", ["colina → k_{cat}", "1.36 kcal/mol = ×10"], NARANJA, "§6, §12"),
        "mm": (520, 330, "Michaelis–Menten", ["K_{M}, V_{max}, k_{cat}/K_{M}"], MAGENTA, "§13"),
        "hi": (280, 330, "Hill", ["sigmoide: sensor", "de glucosa (n ≈ 1.7)"], MAGENTA, "§14"),
        "in": (40, 330, "Inhibidores", ["y activadores: K_{i},", "huellas en L–B"], MAGENTA, "§15"),
        "tp": (400, 490, "Temperatura y pH", ["modulan k vía ΔH‡, ΔS‡ y los pK_{a}"], VIOLETA, "§16"),
    }
    for clave, (x, y, t, lin, col, sec) in nodos.items():
        w = W if clave != "tp" else 240
        L.rect(x, y, w, H + 12, rx=14, fill="#ffffff", stroke=REJILLA, sw=2)
        L.rect(x, y, 7, H + 12, rx=3.5, fill=col)
        L.texto(x + w - 14, y + H + 2, sec, size=12, weight=600, color=TINTA3, anchor="end")
        L.texto(x + 22, y + 30, t, size=16, weight=600)
        L.texto(x + 22, y + 55, lin, size=13, color=TINTA2)
    # flechas fila 1
    for a, b in [("est", "md"), ("md", "qm"), ("qm", "dg")]:
        xa, ya = nodos[a][0] + W, nodos[a][1] + 50
        L.flecha(xa + 6, ya, nodos[b][0] - 6, ya, color=TINTA3, sw=2)
    L.flecha(860, 170 + H + 16, 860, 330 - 6, color=TINTA3, sw=2)
    for a, b in [("ey", "mm"), ("mm", "hi"), ("hi", "in")]:
        L.flecha(nodos[a][0] - 6, nodos[a][1] + 50, nodos[b][0] + W + 6, nodos[b][1] + 50, color=TINTA3, sw=2)
    # T y pH afecta a la fila de cinética
    L.flecha(520, 490, 470, 330 + H + 18, color=VIOLETA, sw=2, dash="5 5")
    L.flecha(600, 490, 640, 330 + H + 18, color=VIOLETA, sw=2, dash="5 5")
    # leyenda de bloques
    ly = 610
    for i, (col, t) in enumerate([(AZUL, "estructura"), (AGUA, "dinámica"), (NARANJA, "barrera y velocidad"),
                                  (MAGENTA, "cinética de laboratorio"), (VIOLETA, "condiciones")]):
        x = 60 + i * 180
        L.circulo(x, ly - 5, 6, fill=col)
        L.texto(x + 14, ly, t, size=13.5, color=TINTA2)
    return L
