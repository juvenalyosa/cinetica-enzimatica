"""Acceso a los datos precalculados del repositorio (modo rápido del notebook).

Los archivos grandes viajan comprimidos (``.gz``); :func:`ruta` devuelve una copia
descomprimida en un caché local la primera vez que se pide.
"""
from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
PRECALCULADO = RAIZ / "data" / "precalculado"
CACHE = RAIZ / "work" / "cache"


def ruta(relativa, base=None):
    """Devuelve la ruta a ``relativa`` dentro de ``data/precalculado``; descomprime ``.gz`` si hace falta."""
    base = PRECALCULADO if base is None else Path(base)
    p = base / relativa
    if p.exists():
        return p
    gz = p.with_name(p.name + ".gz")
    if gz.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        destino = CACHE / p.name
        if not destino.exists() or destino.stat().st_mtime < gz.stat().st_mtime:
            with gzip.open(gz, "rb") as fi, destino.open("wb") as fo:
                shutil.copyfileobj(fi, fo)
        return destino
    raise FileNotFoundError(f"no existe {p} ni {gz}")


def json_(relativa):
    return json.loads(Path(ruta(relativa)).read_text(encoding="utf-8"))


def csv(relativa):
    import pandas as pd

    return pd.read_csv(ruta(relativa))


def leer_xyz_multiple(relativa):
    """Lee un archivo XYZ con varios cuadros -> (símbolos, lista de arrays (n,3), comentarios)."""
    lines = Path(ruta(relativa)).read_text().splitlines()
    frames, comments, symbols, k = [], [], None, 0
    while k < len(lines):
        if not lines[k].strip():
            k += 1
            continue
        n = int(lines[k].split()[0])
        comments.append(lines[k + 1])
        block = lines[k + 2: k + 2 + n]
        symbols = [b.split()[0] for b in block]
        frames.append(np.array([[float(v) for v in b.split()[1:4]] for b in block]))
        k += 2 + n
    return symbols, frames, comments


def energias_de_comentarios(comments):
    out = []
    for c in comments:
        val = np.nan
        for tok in c.split():
            if tok.startswith("E_kcal="):
                val = float(tok.split("=")[1])
        out.append(val)
    return np.array(out)
