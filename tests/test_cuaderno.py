"""El cuaderno se construye desde scripts/cuaderno y todas sus ilustraciones existen y pesan poco."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIG = re.compile(r"\[\[fig:\s*([\w\-]+)")


def _referencias():
    refs = {}
    for mod in sorted((ROOT / "scripts" / "cuaderno").glob("s[0-9][0-9]_*.py")):
        for nombre in FIG.findall(mod.read_text(encoding="utf-8")):
            refs.setdefault(nombre, mod.name)
    return refs


def test_ilustraciones_referenciadas_existen():
    faltan = {n: m for n, m in _referencias().items()
              if not (ROOT / "assets" / "ilustraciones" / f"{n}.png").exists()}
    assert not faltan, f"ilustraciones sin PNG: {faltan}"


def test_ilustraciones_ligeras():
    pesadas = [p.name for p in (ROOT / "assets" / "ilustraciones").glob("*.png") if p.stat().st_size > 200_000]
    assert not pesadas, f"PNG de más de 200 KB: {pesadas}"


def test_cuaderno_incrusta_imagenes():
    import nbformat
    nb = nbformat.read(str(ROOT / "notebooks" / "Cinetica_Enzimatica_Glucoquinasa.ipynb"), as_version=4)
    texto = "\n".join(c.source for c in nb.cells if c.cell_type == "markdown")
    assert "[[fig:" not in texto, "quedó un marcador de figura sin sustituir"
    assert texto.count("data:image/png;base64,") >= len(_referencias())
    secciones = [c.source.splitlines()[0] for c in nb.cells if c.cell_type == "markdown" and c.source.startswith("## ")]
    assert len(secciones) == 17


def test_celdas_de_codigo_ocultas_con_titulo():
    import nbformat
    nb = nbformat.read(str(ROOT / "notebooks" / "Cinetica_Enzimatica_Glucoquinasa.ipynb"), as_version=4)
    for c in nb.cells:
        if c.cell_type == "code":
            assert c.source.startswith("# @title "), c.source[:60]
            assert c.metadata.get("cellView") == "form"
            assert c.metadata.get("jupyter", {}).get("source_hidden") is True
