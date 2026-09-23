"""Secciones 13–16: Michaelis–Menten (el cajero), Hill (el interruptor), inhibidores, temperatura y pH."""
from __future__ import annotations

import math

from .lienzo import (AGUA, AZUL, EJE, MAGENTA, NARANJA, REJILLA, ROJO, TINTA, TINTA2, TINTA3, Lienzo,
                     tinte)

BLANCO = "#ffffff"
BUENO, MALO = "#0ca30c", "#d03b3b"          # colores de estado (viz.STATUS), siempre con ✓/✗


# ---------------------------------------------------------------- pictogramas locales
def _caja(L, x, y, *, w=70, h=46, apagada=False, etiqueta="caja"):
    """Mostrador de la caja (la enzima): rectángulo azul con una ranura de atención a la izquierda."""
    fill = tinte(TINTA3, 0.55) if apagada else AZUL
    L.rect(x, y, w, h, rx=9, fill=fill)
    L.rect(x + 10, y + 10, w - 20, 8, rx=4, fill=BLANCO, opacity=.45)
    if etiqueta:
        L.texto(x + w / 2, y + h + 20, etiqueta, size=13, color=TINTA3, anchor="middle")


def _cliente(L, cx, cy, *, color=NARANJA, r=11, opacity=None):
    L.circulo(cx, cy, r, fill=color, opacity=opacity)


def _panel(L, x, y, w, h):
    L.rect(x, y, w, h, rx=16, fill=BLANCO, stroke=REJILLA, sw=2)


def _barra(L, x, y, w, frac, color, *, h=10):
    L.rect(x, y, w, h, rx=h / 2, fill=tinte(TINTA3, 0.82))
    if frac > 0:
        L.rect(x, y, max(h, w * frac), h, rx=h / 2, fill=color)


def _pacman(L, cx, cy, r, *, fill=AZUL, boca=32, opacity=None):
    """Enzima: un disco con un hueco (el sitio activo) mirando a la derecha."""
    a = math.radians(boca)
    x1, y1 = cx + r * math.cos(a), cy - r * math.sin(a)
    x2, y2 = cx + r * math.cos(a), cy + r * math.sin(a)
    L.camino(f"M {cx + r * 0.25:.2f} {cy:.2f} L {x1:.2f} {y1:.2f} A {r} {r} 0 1 0 {x2:.2f} {y2:.2f} Z",
             fill=fill, stroke=None, opacity=opacity)


def _ejes(L, x0, y0, x1, y1, *, xlab=None, ylab=None, size=14):
    """Ejes en L: origen (x0, y0) abajo a la izquierda, hasta (x1, y1)."""
    L.linea(x0, y0, x1, y0, color=EJE, sw=2)
    L.linea(x0, y0, x0, y1, color=EJE, sw=2)
    if xlab:
        L.texto(x1, y0 + 24, xlab, size=size, color=TINTA2, anchor="end")
    if ylab:
        L.texto(x0 - 10, y1 + 4, ylab, size=size, color=TINTA2, anchor="end")


# ================================================================ s13 Michaelis–Menten
def ilustracion_cajero_saturacion():
    L = Lienzo(800, etiqueta="Analogía", titulo="La enzima es una caja de supermercado",
               subtitulo="Con pocos clientes la caja espera; con una cola larga trabaja a tope y ya no puede ir más rápido")
    paneles = [
        ("Pocos clientes", "[S] ≪ K_{M}", 0.12, ["la caja espera casi siempre:", "cada cliente nuevo cuenta"], 1),
        ("Algunos clientes", "[S] = K_{M}", 0.5, ["ocupada la mitad del tiempo:", "v_{0} = V_{max} / 2"], 2),
        ("Cola larga", "[S] ≫ K_{M}", 0.875, ["siempre ocupada: v_{0} → V_{max};", "da igual que lleguen más"], 3),
    ]
    pw, ph, py = 290, 280, 170
    for i, (tit, cond, frac, texto, n) in enumerate(paneles):
        x = 40 + i * (pw + 25)
        _panel(L, x, py, pw, ph)
        L.numero(x + 30, py + 32, n)
        L.texto(x + 54, py + 38, tit, size=17, weight=600)
        L.mate(x + pw - 18, py + 38, cond, size=17, color=TINTA2, anchor="end")
        # escena: caja a la derecha, cola a la izquierda
        cx, cy = x + pw - 100, py + 85
        _caja(L, cx, cy)
        suelo = cy + 46
        L.linea(x + 18, suelo, x + pw - 18, suelo, color=REJILLA, sw=2)
        if n >= 2:   # cliente atendido (complejo ES)
            L.circulo(cx - 20, suelo - 16, 17, fill=tinte(AGUA, 0.6))
            _cliente(L, cx - 20, suelo - 16)
        cola = {1: [x + 40], 2: [x + 70], 3: [cx - 52 - 27 * k for k in range(5)]}[n]
        for k, qx in enumerate(cola):
            _cliente(L, qx, suelo - 12, opacity=1 if n == 3 else .85)
        if n < 3:
            L.flecha(cola[0] + 18, suelo - 12, cola[0] + 48, suelo - 12, color=TINTA3, sw=2)
        if n == 1:
            L.texto(cx + 35, cy - 12, "libre", size=13, color=TINTA3, anchor="middle")
        # ocupación
        L.texto(x + 22, py + 190, "tiempo ocupada", size=13, color=TINTA3)
        L.texto(x + pw - 22, py + 190, f"{round(frac * 100)} %", size=13, weight=600, color=TINTA2, anchor="end")
        _barra(L, x + 22, py + 200, pw - 44, frac, AGUA)
        L.texto(x + 22, py + 240, texto, size=14.5, color=TINTA2)

    # la hipérbola con los tres tramos
    gx0, gx1, gy0, gy1 = 150, 900, 740, 520
    km_px = (gx1 - gx0) / 10
    def xs(s):
        return gx0 + s * km_px
    def yv(v):
        return gy0 - v * (gy0 - gy1)
    _ejes(L, gx0, gy0, gx1 + 20, gy1 - 30, xlab="[S]  (clientes que llegan)", ylab=None)
    L.texto(gx0 - 12, gy1 - 34, "v₀", size=15, weight=600, anchor="end")
    L.linea(gx0, yv(1), gx1, yv(1), color=EJE, sw=1.5, dash="5 6")
    L.mate(gx0 - 12, yv(1) + 6, "V_{max}", size=17, anchor="end")
    L.linea(gx0, yv(.5), xs(1), yv(.5), color=EJE, sw=1.5, dash="5 6")
    L.linea(xs(1), yv(.5), xs(1), gy0, color=EJE, sw=1.5, dash="5 6")
    L.mate(gx0 - 12, yv(.5) + 6, "V_{max}/2", size=17, anchor="end")
    L.mate(xs(1), gy0 + 24, "K_{M}", size=17, anchor="middle")
    L.curva([(xs(s), yv(s / (1 + s))) for s in [k * 0.125 for k in range(81)]],
            stroke=AZUL, sw=3)
    for s, n, dx, dy in [(0.14, 1, -26, -14), (1.0, 2, 24, 22), (7.0, 3, 0, 30)]:
        L.circulo(xs(s), yv(s / (1 + s)), 8, fill=AZUL, stroke=BLANCO, sw=2.5)
        L.numero(xs(s) + dx, yv(s / (1 + s)) + dy, n, r=12)
    L.texto(xs(3.6), yv(.9) - 38, "meseta: saturación", size=14, color=TINTA2, anchor="middle")
    L.texto(xs(1.25), yv(.1), "tramo casi recto (al principio)", size=14, color=TINTA2)
    return L


def ilustracion_mecanismo_mm():
    L = Lienzo(470, etiqueta="La ecuación, término a término", titulo="El mecanismo detrás de la hipérbola",
               subtitulo="Tres pasos elementales, tres constantes de velocidad")
    cy = 250
    # E + S
    _pacman(L, 115, cy, 42)
    L.texto(115, cy + 72, "E", size=18, weight=600, anchor="middle")
    L.texto(115, cy + 92, "enzima libre", size=13.5, color=TINTA2, anchor="middle")
    L.texto(185, cy + 8, "+", size=28, color=TINTA2, anchor="middle")
    L.circulo(235, cy, 16, fill=NARANJA)
    L.texto(235, cy + 72, "S", size=18, weight=600, anchor="middle")
    L.texto(235, cy + 92, "sustrato", size=13.5, color=TINTA2, anchor="middle")
    # ⇌
    L.flecha(285, cy - 12, 420, cy - 12, color=TINTA, sw=2.5)
    L.flecha(420, cy + 12, 285, cy + 12, color=TINTA3, sw=2.5)
    L.mate(352, cy - 26, "k_{1}", size=21, anchor="middle")
    L.mate(352, cy + 44, "k_{−1}", size=21, anchor="middle")
    # ES
    L.circulo(495, cy, 58, fill=tinte(AGUA, 0.72))
    _pacman(L, 488, cy, 42)
    L.circulo(518, cy, 14, fill=NARANJA)
    L.texto(495, cy + 72, "ES", size=18, weight=600, anchor="middle")
    L.texto(495, cy + 92, "complejo enzima-sustrato", size=13.5, color=TINTA2, anchor="middle")
    # →
    L.flecha(575, cy, 700, cy, color=TINTA, sw=2.5)
    L.mate(637, cy - 16, "k_{2} = k_{cat}", size=21, anchor="middle")
    # E + P
    _pacman(L, 765, cy, 42)
    L.texto(765, cy + 72, "E", size=18, weight=600, anchor="middle")
    L.texto(765, cy + 92, "vuelve libre", size=13.5, color=TINTA2, anchor="middle")
    L.texto(832, cy + 8, "+", size=28, color=TINTA2, anchor="middle")
    L.rect(866, cy - 16, 32, 32, rx=7, fill=MAGENTA)
    L.texto(882, cy + 72, "P", size=18, weight=600, anchor="middle")
    L.texto(882, cy + 92, "producto", size=13.5, color=TINTA2, anchor="middle")
    # explicaciones de cada constante
    notas = [("k_{1}", "el cliente llega a la caja", "(choque productivo)"),
             ("k_{−1}", "se va sin comprar", "(el sustrato se suelta)"),
             ("k_{2}", "paga y sale con la compra", "(cruza la barrera: química)")]
    for i, (k, a, b) in enumerate(notas):
        x = 60 + i * 300
        L.rect(x, 400, 270, 46, rx=10, fill=tinte(TINTA3, 0.9))
        L.mate(x + 16, 429, k, size=18)
        L.texto(x + 66, 421, a, size=14.5, weight=600)
        L.texto(x + 66, 439, b, size=13, color=TINTA2)
    L.texto(60, 390, "En la analogía del cajero:", size=14, color=TINTA3)
    return L


def ilustracion_mm_anatomia():
    L = Lienzo(445, etiqueta="La ecuación, término a término",
               titulo="Michaelis–Menten en palabras: velocidad = máxima × fracción ocupada")
    y = 205
    L.mate(60, y, "v_{0}", size=44)
    L.mate(128, y, "=", size=40, color=TINTA2)
    cajas = [(170, 170, AZUL, "V_{max}", "velocidad máxima", ["toda la enzima trabajando:", "k_{cat} · [E]_{0}"]),
             (383, 44, None, "×", None, None),
             (455, 270, AGUA, "[S] / (K_{M} + [S])", "fracción de enzima ocupada",
              ["de 0 (caja vacía) a 1", "(caja siempre ocupada)"])]
    for x, w, color, expr, nombre, expl in cajas:
        if nombre is None:
            L.mate(x + w / 2, y, expr, size=36, color=TINTA2, anchor="middle")
            continue
        L.rect(x, y - 52, w, 76, rx=14, fill=tinte(color, 0.84))
        L.mate(x + w / 2, y + 2, expr, size=26 if len(expr) > 10 else 30, anchor="middle")
        L.rect(x, y + 40, w, 5, rx=2.5, fill=color)
        L.texto(x + w / 2, y + 76, nombre, size=17, weight=600, anchor="middle")
        L.texto(x + w / 2, y + 100, expl, size=14.5, color=TINTA2, anchor="middle")
    # tres casos con barras de ocupación
    x = 755
    L.rect(x, 120, 205, 290, rx=16, fill=BLANCO, stroke=REJILLA, sw=2)
    L.texto(x + 20, 155, "La fracción ocupada", size=16, weight=600)
    casos = [("[S] = K_{M}/10", 1 / 11, "≈ 9 %"), ("[S] = K_{M}", 0.5, "50 %"), ("[S] = 10 K_{M}", 10 / 11, "≈ 91 %")]
    for i, (cond, fr, lab) in enumerate(casos):
        yy = 195 + i * 70
        L.mate(x + 20, yy, cond, size=16, color=TINTA2)
        L.texto(x + 185, yy, lab, size=14, weight=600, color=TINTA2, anchor="end")
        _barra(L, x + 20, yy + 14, 165, fr, AGUA)
    L.texto(x + 20, 395, "K_{M} = [S] a media ocupación", size=13, color=TINTA3)
    return L


# ================================================================ s14 Hill
def _curva_panel(L, x0, y0, w, h, n, color, titulo, icono, etiqueta_ventana):
    _panel(L, x0, y0, w, h)
    L.texto(x0 + 76, y0 + 40, titulo[0], size=17, weight=600)
    L.texto(x0 + 76, y0 + 62, titulo[1], size=14, color=TINTA2)
    icono(L, x0 + 40, y0 + 44)
    gx0, gx1, gy0, gy1 = x0 + 65, x0 + w - 30, y0 + h - 100, y0 + 110
    _ejes(L, gx0, gy0, gx1, gy1 - 10, size=13)
    L.texto(gx0 - 8, gy1 - 6, "v₀", size=14, weight=600, anchor="end")
    lmin, lmax = -2.0, 2.0          # eje x logarítmico: log10([S]/S0.5)
    def xs(s):
        return gx0 + (math.log10(s) - lmin) / (lmax - lmin) * (gx1 - gx0)
    def yv(v):
        return gy0 - v * (gy0 - gy1)
    for t, lab in [(0.1, "0.1"), (1, "1"), (10, "10")]:
        L.linea(xs(t), gy0, xs(t), gy0 + 6, color=EJE, sw=2)
        L.texto(xs(t), gy0 + 22, lab, size=12.5, color=TINTA3, anchor="middle")
    L.texto(gx1, gy0 + 22, "[S] / S₀.₅", size=12.5, color=TINTA3, anchor="end")
    s10, s90 = (1 / 9) ** (1 / n), 9 ** (1 / n)
    L.rect(xs(s10), gy1, xs(s90) - xs(s10), gy0 - gy1, rx=0, fill=tinte(color, 0.86))
    for sv in (s10, s90):
        L.linea(xs(sv), gy1, xs(sv), gy0, color=EJE, sw=1.5, dash="4 5")
        L.linea(xs(sv), gy0 + 32, xs(sv), gy0 + 52, color=EJE, sw=1.5)
    L.linea(gx0, yv(.1), xs(s10), yv(.1), color=EJE, sw=1.2, dash="3 5")
    L.linea(gx0, yv(.9), xs(s90), yv(.9), color=EJE, sw=1.2, dash="3 5")
    L.texto(gx0 - 8, yv(.1) + 5, "10 %", size=12.5, color=TINTA3, anchor="end")
    L.texto(gx0 - 8, yv(.9) + 5, "90 %", size=12.5, color=TINTA3, anchor="end")
    pts = [10 ** (lmin + k * (lmax - lmin) / 120) for k in range(121)]
    L.curva([(xs(sv), yv(sv ** n / (1 + sv ** n))) for sv in pts], stroke=color, sw=3)
    L.flecha(xs(s10) + 3, gy0 + 42, xs(s90) - 3, gy0 + 42, color=TINTA2, sw=2, doble=True)
    L.texto((xs(s10) + xs(s90)) / 2, gy0 + 72, etiqueta_ventana, size=16, weight=600, anchor="middle")


SUPERFICIE_BLANCA = BLANCO


def _icono_regulador(L, cx, cy):
    L.circulo(cx, cy, 20, fill=tinte(AZUL, 0.85), stroke=AZUL, sw=2.5)
    L.camino(f"M {cx - 13} {cy + 13} A 18.4 18.4 0 1 1 {cx + 13} {cy + 13}", stroke=AZUL, sw=3, fill="none")
    L.linea(cx, cy, cx + 10, cy - 11, color=TINTA, sw=3)


def _icono_interruptor(L, cx, cy):
    L.rect(cx - 13, cy - 22, 26, 44, rx=7, fill=tinte(NARANJA, 0.85), stroke=NARANJA, sw=2.5)
    L.rect(cx - 7, cy - 16, 14, 16, rx=4, fill=NARANJA)


def ilustracion_regulador_interruptor():
    L = Lienzo(560, etiqueta="Analogía", titulo="Regulador de luz frente a interruptor con umbral",
               subtitulo="Franja sombreada: cuánto hay que multiplicar [S] para pasar del 10 % al 90 % (eje x logarítmico)")
    _curva_panel(L, 40, 165, 450, 370, 1.0, AZUL, ("Regulador de luz", "hipérbola, n = 1 (Michaelis–Menten)"),
                 _icono_regulador, "[S] × 81")
    _curva_panel(L, 510, 165, 450, 370, 1.7, NARANJA, ("Interruptor con umbral", "sigmoide, n = 1.7 (glucoquinasa)"),
                 _icono_interruptor, "[S] × 13")
    return L


def ilustracion_cooperatividad_cinetica():
    L = Lienzo(615, etiqueta="Analogía", titulo="Cooperatividad cinética: una enzima con memoria",
               subtitulo="Un solo sitio activo, dos formas que se interconvierten tan despacio como la propia catálisis")
    # las dos formas
    def almeja(cx, cy, abierta, fill, fill_tapa):
        """Dos dominios articulados a la izquierda: el grande abajo y el pequeño (tapa) arriba."""
        L.camino(f"M {cx - 55} {cy} A 55 55 0 0 0 {cx + 55} {cy} Z", fill=fill, stroke=None)
        ang = -30 if abierta else 0
        L.add(f'<g transform="rotate({ang} {cx - 55} {cy - 3})">'
              f'<path d="M {cx - 55} {cy - 3} A 55 44 0 0 1 {cx + 55} {cy - 3} Z" fill="{fill_tapa}"/></g>')
        L.circulo(cx - 55, cy - 1, 6, fill=TINTA2)
    almeja(210, 280, True, tinte(AZUL, 0.45), tinte(AZUL, 0.62))
    L.texto(210, 362, "forma abierta", size=17, weight=600, anchor="middle")
    L.texto(210, 384, "poco activa (lenta)", size=14, color=TINTA2, anchor="middle")
    almeja(700, 280, False, AZUL, tinte(AZUL, 0.25))
    L.circulo(705, 274, 12, fill=NARANJA)
    L.texto(700, 362, "forma cerrada", size=17, weight=600, anchor="middle")
    L.texto(700, 384, "activa (rápida)", size=14, color=TINTA2, anchor="middle")
    L.flecha(300, 240, 610, 240, color=TINTA2, sw=2.5, curvatura=-0.12)
    L.flecha(610, 300, 300, 300, color=TINTA3, sw=2.5, curvatura=-0.12)
    L.texto(455, 205, "se cierra al unir glucosa", size=14.5, color=TINTA2, anchor="middle")
    L.texto(455, 345, "se relaja despacio si no llega otra glucosa", size=14.5, color=TINTA2, anchor="middle")
    L.mate(455, 283, "k_{ex} ≈ 5–100 s^{−1}  vs  k_{cat} ≈ 60 s^{−1}", size=16, color=TINTA2, anchor="middle")
    # dos líneas de tiempo
    filas = [("poca glucosa", [120, 400, 690], "entre una glucosa y la siguiente se relaja: casi siempre lenta"),
             ("mucha glucosa", [120, 190, 260, 330, 400, 470, 540, 610, 680, 750], "la siguiente llega antes de relajarse: sigue rápida")]
    for i, (nombre, llegadas, lema) in enumerate(filas):
        y = 440 + i * 80
        L.texto(40, y + 6, nombre, size=15, weight=600)
        x0, x1 = 175, 900
        L.rect(x0, y - 8, x1 - x0, 16, rx=8, fill=tinte(AZUL, 0.72))
        for a in llegadas:
            xa = x0 + (a - 100) * (x1 - x0) / 700
            rap = min(x1, xa + (95 if i == 0 else 80))
            L.rect(xa, y - 8, rap - xa, 16, rx=8, fill=AZUL)
            L.circulo(xa, y - 22, 7, fill=NARANJA)
        L.texto(x0, y + 32, lema, size=13.5, color=TINTA2)
    ly = 585
    L.rect(560, ly - 9, 18, 10, rx=5, fill=AZUL)
    L.texto(584, ly, "forma rápida", size=13, color=TINTA2)
    L.rect(690, ly - 9, 18, 10, rx=5, fill=tinte(AZUL, 0.72))
    L.texto(714, ly, "forma lenta", size=13, color=TINTA2)
    L.circulo(815, ly - 4, 6, fill=NARANJA)
    L.texto(827, ly, "llega una glucosa", size=13, color=TINTA2)
    return L


# ================================================================ s15 inhibición
def _lb_mini(L, x0, y0, w, h, rectas):
    """Mini Lineweaver–Burk: ejes cruzados en (0, 0); ``rectas`` = [(pendiente, corte, color)] en unidades Km=Vmax=1."""
    xmin, xmax, ymax = -1.6, 2.6, 6.2
    def X(u):
        return x0 + (u - xmin) / (xmax - xmin) * w
    def Y(v):
        return y0 + h - v / ymax * h
    L.linea(X(xmin), Y(0), X(xmax), Y(0), color=EJE, sw=1.8)
    L.linea(X(0), Y(0), X(0), Y(ymax), color=EJE, sw=1.8)
    L.texto(X(xmax), Y(0) + 20, "1/[S]", size=12.5, color=TINTA3, anchor="end")
    L.texto(X(0) + 6, Y(ymax) + 4, "1/v₀", size=12.5, color=TINTA3)
    for m, b, color in rectas:
        u0 = max(xmin, -b / m)
        L.linea(X(u0), Y(m * u0 + b), X(xmax), Y(m * xmax + b), color=color, sw=2.6)


def ilustracion_inhibidores_cajero():
    L = Lienzo(720, etiqueta="Analogía", titulo="Tres maneras de frenar la caja",
               subtitulo="Cada tipo de inhibidor deja una huella distinta en V_{max}, en K_{M} y en el gráfico de Lineweaver–Burk")
    tipos = [
        ("Competitivo", "se pone en la cola sin comprar", "ocupa el turno (se une a E)", "igual", "sube", [(1, 1, AZUL), (2, 1, ROJO)],
         "se cruzan en el eje vertical"),
        ("Acompetitivo", "traba la caja con un cliente", "atrapa al complejo ES", "baja", "baja", [(1, 1, AZUL), (1, 2, ROJO)],
         "rectas paralelas"),
        ("No competitivo", "apaga la luz de la caja", "se une a E y a ES por igual", "baja", "igual", [(1, 1, AZUL), (2, 2, ROJO)],
         "se cruzan en el eje horizontal"),
    ]
    pw = 290
    for i, (nom, analog, quimica, vmax, km, rectas, huella) in enumerate(tipos):
        x = 40 + i * (pw + 25)
        _panel(L, x, 165, pw, 530)
        L.texto(x + 22, 202, nom, size=18, weight=600)
        L.texto(x + 22, 225, analog, size=14.5, color=TINTA2)
        # escena
        cx, cy = x + pw - 105, 262
        suelo = cy + 46
        L.linea(x + 18, suelo, x + pw - 18, suelo, color=REJILLA, sw=2)
        _caja(L, cx, cy, apagada=(i == 2), etiqueta=None)
        if i == 0:
            L.circulo(cx - 20, suelo - 16, 15, fill=ROJO)
            L.texto(cx - 20, suelo - 11, "I", size=15, weight=700, color=BLANCO, anchor="middle")
            for k in range(3):
                _cliente(L, cx - 58 - 28 * k, suelo - 12)
        elif i == 1:
            L.circulo(cx - 20, suelo - 16, 17, fill=tinte(AGUA, 0.6))
            _cliente(L, cx - 20, suelo - 16)
            L.rect(cx - 42, cy - 16, 122, 14, rx=7, fill=ROJO)
            L.texto(cx + 19, cy - 5, "I", size=12, weight=700, color=BLANCO, anchor="middle")
            for k in range(2):
                _cliente(L, cx - 62 - 28 * k, suelo - 12)
        else:
            L.circulo(cx - 20, suelo - 16, 17, fill=tinte(AGUA, 0.6))
            _cliente(L, cx - 20, suelo - 16)
            bx, by = cx + 35, cy - 26
            L.circulo(bx, by, 11, fill=tinte(TINTA3, 0.6), stroke=ROJO, sw=2.5)
            L.linea(bx - 11, by + 11, bx + 11, by - 11, color=ROJO, sw=2.5)
            for k in range(2):
                _cliente(L, cx - 62 - 28 * k, suelo - 12)
        L.texto(x + 22, suelo + 30, quimica, size=14, color=TINTA2)
        # efectos
        for j, (lab, val) in enumerate([("V_{max} aparente", vmax), ("K_{M} aparente", km)]):
            yy = suelo + 70 + j * 32
            L.mate(x + 22, yy, lab, size=16, color=TINTA2)
            flecha = {"igual": "=  igual", "sube": "↑  sube", "baja": "↓  baja"}[val]
            L.texto(x + pw - 22, yy, flecha, size=15, weight=600, anchor="end")
        L.linea(x + 18, suelo + 100 + 20, x + pw - 18, suelo + 100 + 20, color=REJILLA, sw=1.5)
        _lb_mini(L, x + 30, suelo + 140, pw - 60, 175, rectas)
        L.texto(x + pw / 2, 680, huella, size=14, weight=600, color=TINTA2, anchor="middle")
    # leyenda
    L.linea(640, 145, 670, 145, color=AZUL, sw=2.6)
    L.texto(678, 150, "sin inhibidor", size=13.5, color=TINTA2)
    L.linea(790, 145, 820, 145, color=ROJO, sw=2.6)
    L.texto(828, 150, "con inhibidor", size=13.5, color=TINTA2)
    return L


# ================================================================ s16 temperatura y pH
def ilustracion_temperatura_sacudidas():
    L = Lienzo(520, etiqueta="Analogía", titulo="Calentar = sacudidas más fuertes",
               subtitulo="Cuántas moléculas tienen cada energía; solo la cola que supera la barrera reacciona")
    gx0, gx1, gy0, gy1 = 90, 700, 440, 190
    emax = 6.0
    def xs(e):
        return gx0 + e / emax * (gx1 - gx0)
    def f(e, kt):
        return 2 / math.sqrt(math.pi) * kt ** -1.5 * math.sqrt(max(e, 0)) * math.exp(-e / kt)
    esc = (gy0 - gy1) / 0.5
    def yv(v):
        return gy0 - v * esc
    barrera = 3.2
    temps = [(1.0, tinte(AZUL, 0.45), "25 °C"), (1.3, AZUL, "37 °C")]
    for kt, color, _ in temps:
        cola = [(xs(e), yv(f(e, kt))) for e in [barrera + k * 0.05 for k in range(int((emax - barrera) / 0.05) + 1)]]
        L.poligono(cola + [(xs(emax), gy0), (xs(barrera), gy0)], fill=tinte(NARANJA, 0.55 if kt > 1.1 else 0.3),
                   stroke=None, opacity=.85)
    for kt, color, lab in temps:
        L.curva([(xs(e), yv(f(e, kt))) for e in [k * 0.05 for k in range(int(emax / 0.05) + 1)]], stroke=color, sw=3)
    _ejes(L, gx0, gy0, gx1 + 10, gy1 - 20, xlab="energía de la molécula →")
    L.texto(gx0 - 10, gy1 - 16, ["nº de", "moléculas"], size=13, color=TINTA2, anchor="end")
    L.linea(xs(barrera), gy1 - 10, xs(barrera), gy0, color=TINTA, sw=2, dash="6 5")
    L.texto(xs(barrera) + 8, gy1 - 2, "barrera", size=15, weight=600)
    L.mate(xs(barrera) + 8, gy1 + 20, "E_{a}", size=18, color=TINTA2)
    L.flecha(xs(4.35), yv(0.12) - 10, xs(3.9), yv(0.05), color=TINTA2, sw=1.8)
    L.texto(xs(4.35) + 4, yv(0.12) - 16, "estas sí reaccionan", size=14, color=TINTA2)
    # leyenda
    for i, (kt, color, lab) in enumerate(temps):
        L.linea(xs(1.3), 215 + i * 28, xs(1.3) + 34, 215 + i * 28, color=color, sw=3)
        L.texto(xs(1.3) + 44, 220 + i * 28, lab, size=14.5, color=TINTA2)
    # nota
    L.nota(730, 190, 240, ["Con E_{a} = 12 kcal/mol,", "de 25 a 37 °C la cola", "crece ~2.2 veces:", "k se multiplica por 2.2."],
           color=NARANJA, titulo="En números", size=15)
    L.texto(730, 400, ["Dibujo esquemático: la", "diferencia entre las curvas", "está exagerada para verla."],
            size=13, color=TINTA3)
    return L


def ilustracion_ph_campana():
    L = Lienzo(640, etiqueta="Analogía", titulo="La campana de pH: dos grupos, dos condiciones",
               subtitulo="Asp205 debe estar sin protón (para aceptarlo) y Lys169 con protón (para sujetar el fosfato)")
    gx0, gx1, gy0, gy1 = 110, 900, 400, 185
    pk1, pk2 = 3.0, 7.0          # posiciones esquemáticas en un eje sin números
    def xs(p):
        return gx0 + p / 10 * (gx1 - gx0)
    def act(p):
        return 1 / (1 + 10 ** (pk1 - p) + 10 ** (p - pk2))
    amax = act(5.0)
    def yv(a):
        return gy0 - a / amax * (gy0 - gy1)
    L.poligono([(xs(p), yv(act(p))) for p in [k * 0.05 for k in range(201)]] + [(xs(10), gy0), (xs(0), gy0)],
               fill=tinte(AZUL, 0.9), stroke=None)
    L.curva([(xs(p), yv(act(p))) for p in [k * 0.05 for k in range(201)]], stroke=AZUL, sw=3)
    _ejes(L, gx0, gy0, gx1 + 10, gy1 - 25, xlab="pH →")
    L.texto(gx0 - 10, gy1 - 20, "actividad", size=13.5, color=TINTA2, anchor="end")
    for p, lab in [(pk1, "pK_{1}"), (pk2, "pK_{2}")]:
        L.linea(xs(p), yv(act(p)), xs(p), gy0, color=EJE, sw=1.5, dash="4 5")
        L.mate(xs(p), gy0 + 24, lab, size=17, anchor="middle")
    L.linea(xs(5), gy1, xs(5), gy0, color=TINTA3, sw=1.5, dash="2 5")
    L.texto(xs(5), gy1 - 12, "óptimo = (pK₁ + pK₂) / 2", size=14.5, weight=600, anchor="middle")
    L.texto(xs(5), gy0 + 24, "glucoquinasa: pH 8.5–8.7", size=13, color=TINTA3, anchor="middle")
    # tarjetas de estado
    tarjetas = [("pH bajo", [("Asp205", "COOH", False), ("Lys169", "NH₃⁺", True)], "la base ya tiene protón: no acepta otro"),
                ("pH óptimo", [("Asp205", "COO⁻", True), ("Lys169", "NH₃⁺", True)], "los dos grupos listos"),
                ("pH alto", [("Asp205", "COO⁻", True), ("Lys169", "NH₂", False)], "la lisina perdió su carga positiva")]
    tw = 290
    for i, (nom, grupos, lema) in enumerate(tarjetas):
        x = 40 + i * (tw + 25)
        y = 455
        L.rect(x, y, tw, 160, rx=14, fill=BLANCO, stroke=REJILLA, sw=2)
        L.texto(x + 20, y + 32, nom, size=16, weight=600)
        for j, (res, estado, ok) in enumerate(grupos):
            yy = y + 66 + j * 32
            L.texto(x + 20, yy, res, size=14.5, color=TINTA2)
            L.texto(x + 110, yy, estado, size=15, weight=600, font="'STIX Two Text', serif")
            L.circulo(x + tw - 34, yy - 5, 11, fill=BUENO if ok else MALO)
            L.texto(x + tw - 34, yy, "✓" if ok else "✗", size=13, weight=700, color=BLANCO, anchor="middle")
        L.texto(x + 20, y + 140, lema, size=13.5, color=TINTA2)
    return L
