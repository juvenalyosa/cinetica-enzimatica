#!/usr/bin/env python
"""Genera notebooks/Cinetica_Enzimatica_Glucoquinasa.ipynb a partir de las celdas definidas en
``scripts/cuaderno/sNN_*.py`` (una sección por módulo, en orden de nombre).

Mantener el notebook como código facilita revisarlo, versionarlo y regenerarlo:
    python scripts/construir_notebook.py

Las ilustraciones se incrustan desde assets/ilustraciones/*.png (ver scripts/04_ilustraciones.py).
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "notebooks" / "Cinetica_Enzimatica_Glucoquinasa.ipynb"
sys.path.insert(0, str(Path(__file__).resolve().parent))

from cuaderno.celdas import CELLS  # noqa: E402


def main():
    for mod in sorted((ROOT / "scripts" / "cuaderno").glob("s[0-9][0-9]_*.py")):
        importlib.import_module(f"cuaderno.{mod.stem}")
    nb = nbf.v4.new_notebook()
    nb["cells"] = CELLS
    nb["metadata"] = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
        "colab": {"name": "Cinética enzimática con QM/MM: glucoquinasa", "provenance": [], "toc_visible": True},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(nb, str(OUT))
    print("escrito", OUT, "con", len(CELLS), "celdas")


if __name__ == "__main__":
    main()
