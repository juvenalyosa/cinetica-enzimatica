#!/usr/bin/env python
"""Genera los cuadernos del curso a partir de ``scripts/cuaderno/<idioma>/sNN_*.py`` (una sección por módulo):

* ``es`` → notebooks/Cinetica_Enzimatica_Glucoquinasa.ipynb
* ``en`` → notebooks/Enzyme_Kinetics_Glucokinase.ipynb

Los dos comparten las celdas de código (solo cambian los textos; lo comprueba tests/test_cuaderno.py), los datos
y el paquete ``enzimas``. Las ilustraciones se incrustan desde assets/ilustraciones/<idioma>/*.png.

    python scripts/construir_notebook.py            # los dos
    python scripts/construir_notebook.py es         # solo uno
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from cuaderno import celdas  # noqa: E402

CUADERNOS = {
    "es": ("Cinetica_Enzimatica_Glucoquinasa.ipynb", "Cinética enzimática con QM/MM: glucoquinasa"),
    "en": ("Enzyme_Kinetics_Glucokinase.ipynb", "Enzyme kinetics with QM/MM: glucokinase"),
}


def construir(idioma):
    archivo, nombre = CUADERNOS[idioma]
    celdas.reiniciar(idioma)
    for mod in sorted((ROOT / "scripts" / "cuaderno" / idioma).glob("s[0-9][0-9]_*.py")):
        nombre_mod = f"cuaderno.{idioma}.{mod.stem}"
        if nombre_mod in sys.modules:
            importlib.reload(sys.modules[nombre_mod])
        else:
            importlib.import_module(nombre_mod)
    nb = nbf.v4.new_notebook()
    nb["cells"] = list(celdas.CELLS)
    nb["metadata"] = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
        "colab": {"name": nombre, "provenance": [], "toc_visible": True},
    }
    out = ROOT / "notebooks" / archivo
    out.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(nb, str(out))
    print("escrito", out, "con", len(nb["cells"]), "celdas")
    return out


def main(idiomas):
    for idioma in idiomas:
        if not any((ROOT / "scripts" / "cuaderno" / idioma).glob("s[0-9][0-9]_*.py")):
            print(f"(sin secciones en scripts/cuaderno/{idioma}: se omite)")
            continue
        construir(idioma)


if __name__ == "__main__":
    main(sys.argv[1:] or ["es", "en"])
