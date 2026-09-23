"""Sección 6: la pelota que intenta cruzar la colina y la anatomía de la ecuación de Eyring."""
from __future__ import annotations

from .lienzo import (AZUL, NARANJA, AGUA, TINTA, TINTA2, TINTA3, REJILLA, EJE, Lienzo, muestrear,
                     perfil_barrera, tinte)


def ilustracion_colina_intentos():
    L = Lienzo(615, etiqueta="Analogía", titulo="Una reacción es una pelota que intenta cruzar una colina",
               subtitulo="Rebota miles de millones de veces por segundo; solo muy pocos intentos llegan a la cima")
    x0, x1 = 70, 930
    base, cima, fin = 470, 215, 520
    f = perfil_barrera(x0, x1, base, cima, fin, centro=0.52, ancho=0.17)
    pts = muestrear(f, x0, x1, 160)
    # relleno del terreno
    L.poligono(pts + [(x1, 535), (x0, 535)], fill=tinte(TINTA3, 0.9), stroke=None)
    L.curva(pts, stroke=TINTA2, sw=3)

    # intentos fallidos: arcos grises que suben la ladera y vuelven
    xv = 205
    yv = f(xv)
    for alto, dx, op in [(45, 95, .55), (75, 140, .7), (105, 185, .85)]:
        xa = xv + 20
        L.camino(f"M {xa} {yv - 8} Q {xa + dx * 0.45} {f(xa + dx * 0.4) - alto * 1.9} {xa + dx} {f(xa + dx) - 12}",
                 stroke=TINTA3, sw=2, dash="3 7", opacity=op, flecha=True)
    # el intento con éxito
    xts = x0 + (x1 - x0) * 0.52
    L.camino(f"M {xv - 2} {yv - 34} C {xv + 40} {yv - 330} {xts + 30} {cima - 110} {xts + 190} {f(xts + 190) - 16}",
             stroke=NARANJA, sw=3, flecha=True)
    # pelota
    L.circulo(xv, yv - 16, 16, fill=AZUL)
    L.circulo(xv - 5, yv - 21, 5, fill="#ffffff", opacity=.55)
    L.circulo(xts + 195, f(xts + 195) - 16, 16, fill=AGUA, opacity=.95)

    # cima = estado de transición
    L.circulo(xts, cima, 7, fill=NARANJA)
    L.linea(xts + 6, cima - 6, xts + 40, cima - 40, color=TINTA3, sw=1.5)
    L.texto(xts + 46, cima - 58, "estado de transición", size=16, weight=600)
    L.texto(xts + 46, cima - 38, "(la cima de la colina)", size=14, color=TINTA2)

    # altura de la barrera
    xb = 120
    L.linea(xb - 20, base, xv + 60, base, color=EJE, sw=1.5, dash="4 5")
    L.linea(xb - 20, cima, xts - 20, cima, color=EJE, sw=1.5, dash="4 5")
    L.flecha(xb, base - 4, xb, cima + 6, color=TINTA, sw=2.2, doble=True)
    L.mate(xb + 14, (base + cima) / 2 - 4, "ΔG^{‡}", size=26)
    L.texto(xb + 14, (base + cima) / 2 + 20, ["altura de", "la colina"], size=14, color=TINTA2)

    L.texto(xv, 560, "reactivos", size=17, weight=600, anchor="middle")
    L.texto(xv, 580, "glucosa + ATP", size=14, color=TINTA2, anchor="middle")
    L.texto(800, 560, "productos", size=17, weight=600, anchor="middle")
    L.texto(800, 580, "G6P + ADP", size=14, color=TINTA2, anchor="middle")

    # leyenda de los arcos
    L.linea(360, 560, 400, 560, color=TINTA3, sw=2, dash="3 7")
    L.texto(410, 565, "intento fallido: vuelve al valle", size=15, color=TINTA2)
    L.linea(360, 588, 400, 588, color=NARANJA, sw=3)
    L.texto(410, 593, "intento con éxito: cruza la cima", size=15, color=TINTA2)
    return L


def ilustracion_eyring_anatomia():
    L = Lienzo(445, etiqueta="La ecuación, término a término",
               titulo="Eyring en palabras: velocidad = intentos × probabilidad de éxito")
    y = 205
    # ecuación grande con cajas de color bajo cada término
    L.mate(70, y, "k", size=46)
    L.mate(110, y, "=", size=40, color=TINTA2)
    cajas = [
        (165, 175, AZUL, "κ · k_{B}T / h", "intentos por segundo",
         ["la molécula «choca» contra", "la colina ~6 × 10^{12} veces/s"]),
        (395, 50, TINTA3, "×", None, None),
        (455, 250, NARANJA, "exp(−ΔG^{‡} / RT)", "probabilidad de llegar arriba",
         ["fracción de intentos con", "energía suficiente (muy pequeña)"]),
    ]
    for x, w, color, expr, nombre, expl in cajas:
        if nombre is None:
            L.mate(x + w / 2, y, expr, size=36, color=TINTA2, anchor="middle")
            continue
        L.rect(x, y - 52, w, 76, rx=14, fill=tinte(color, 0.88))
        L.mate(x + w / 2, y + 2, expr, size=32, anchor="middle")
        L.rect(x, y + 40, w, 5, rx=2.5, fill=color)
        L.texto(x + w / 2, y + 76, nombre, size=17, weight=600, anchor="middle")
        L.texto(x + w / 2, y + 100, expl, size=14.5, color=TINTA2, anchor="middle")
    # a la derecha: la regla de oro
    x = 745
    L.rect(x, 120, 215, 285, rx=16, fill="#ffffff", stroke=REJILLA, sw=2)
    L.texto(x + 20, 155, "Regla de oro", size=17, weight=600)
    L.texto(x + 20, 180, ["cada 1.36 kcal/mol más", "de colina ⇒ 10 veces", "más lenta"], size=15, color=TINTA2)
    # mini escalera
    for i, (h, lab) in enumerate([(0, "×1"), (1, "÷10"), (2, "÷100"), (3, "÷1000")]):
        xx, yy = x + 30 + i * 45, 370 - h * 28
        L.rect(xx, yy, 38, 370 - yy + 12, rx=5, fill=tinte(NARANJA, 0.55 + 0.1 * (3 - h)))
        L.texto(xx + 19, yy - 8, lab, size=13, color=TINTA2, anchor="middle", weight=600)
    L.texto(x + 107, 272, "más barrera →", size=13, color=TINTA3, anchor="middle")
    return L
