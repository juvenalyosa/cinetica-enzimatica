#!/usr/bin/env python
"""Genera las ilustraciones del cuaderno en assets/ilustraciones/ (SVG fuente + PNG ×1.6).

Cada módulo ``scripts/ilustraciones/fig_*.py`` define funciones ``ilustracion_<nombre>()`` que
devuelven un :class:`Lienzo`; este script las descubre, las guarda y las rasteriza con
``rsvg-convert`` (hace falta localmente; el cuaderno solo lee los PNG ya generados).

    python scripts/04_ilustraciones.py            # todas
    python scripts/04_ilustraciones.py colina     # solo las que contienen "colina"
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))


def main(filtro=None):
    hechas = []
    for mod_path in sorted((AQUI / "ilustraciones").glob("fig_*.py")):
        mod = importlib.import_module(f"ilustraciones.{mod_path.stem}")
        for nombre in sorted(n for n in dir(mod) if n.startswith("ilustracion_")):
            corto = nombre.removeprefix("ilustracion_")
            if filtro and filtro not in corto:
                continue
            png = getattr(mod, nombre)().guardar(corto)
            hechas.append(f"{corto:32s} {png.stat().st_size / 1024:6.0f} KB")
    print("\n".join(hechas))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
