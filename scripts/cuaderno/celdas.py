"""Constructores de celdas del cuaderno.

``md`` admite un marcador para insertar ilustraciones de ``assets/ilustraciones``:

    [[fig:colina_intentos | texto alternativo]]

que se sustituye por la imagen PNG incrustada en base64 (una data URI). Así el cuaderno es
autónomo: se ve igual en Colab, en Jupyter y sin conexión, sin depender de rutas relativas.
Las imágenes se generan con ``scripts/04_ilustraciones.py``; aquí solo se leen.
"""
from __future__ import annotations

import base64
import re
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[2]
ILUSTRACIONES = ROOT / "assets" / "ilustraciones"

CELLS: list = []
_FIG = re.compile(r"\[\[fig:\s*([\w\-]+)\s*(?:\|\s*([^\]]*?))?\s*\]\]")


def _imagen(m):
    nombre, alt = m.group(1), (m.group(2) or m.group(1)).strip()
    png = ILUSTRACIONES / f"{nombre}.png"
    datos = base64.b64encode(png.read_bytes()).decode("ascii")
    return f"![{alt}](data:image/png;base64,{datos})"


def md(text):
    cell = nbf.v4.new_markdown_cell(_FIG.sub(_imagen, text.strip("\n")))
    cell["id"] = f"celda-{len(CELLS):03d}"  # id determinista: el notebook se regenera byte a byte igual
    CELLS.append(cell)


def code(text):
    cell = nbf.v4.new_code_cell(text.strip("\n"))
    cell["id"] = f"celda-{len(CELLS):03d}"
    CELLS.append(cell)
