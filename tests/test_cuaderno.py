"""Los cuadernos (español e inglés) se construyen desde scripts/cuaderno/<idioma>, con sus ilustraciones, celdas
ocultas y el mismo código; el inglés no deja texto en español."""
from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIG = re.compile(r"\[\[fig:\s*([\w\-]+)")
CUADERNOS = {"es": "Cinetica_Enzimatica_Glucoquinasa.ipynb", "en": "Enzyme_Kinetics_Glucokinase.ipynb"}
IDIOMAS = [i for i in CUADERNOS if (ROOT / "notebooks" / CUADERNOS[i]).exists()
           and any((ROOT / "scripts" / "cuaderno" / i).glob("s[0-9][0-9]_*.py"))]


def _nb(idioma):
    import nbformat
    return nbformat.read(str(ROOT / "notebooks" / CUADERNOS[idioma]), as_version=4)


def _referencias(idioma):
    refs = {}
    for mod in sorted((ROOT / "scripts" / "cuaderno" / idioma).glob("s[0-9][0-9]_*.py")):
        for nombre in FIG.findall(mod.read_text(encoding="utf-8")):
            refs.setdefault(nombre, mod.name)
    return refs


@pytest.mark.parametrize("idioma", IDIOMAS)
def test_ilustraciones_referenciadas_existen(idioma):
    faltan = {n: m for n, m in _referencias(idioma).items()
              if not (ROOT / "assets" / "ilustraciones" / idioma / f"{n}.png").exists()}
    assert not faltan, f"ilustraciones sin PNG ({idioma}): {faltan}"


def test_ilustraciones_ligeras():
    pesadas = [str(p) for p in (ROOT / "assets" / "ilustraciones").rglob("*.png") if p.stat().st_size > 200_000]
    assert not pesadas, f"PNG de más de 200 KB: {pesadas}"


@pytest.mark.parametrize("idioma", IDIOMAS)
def test_cuaderno_incrusta_imagenes(idioma):
    nb = _nb(idioma)
    texto = "\n".join(c.source for c in nb.cells if c.cell_type == "markdown")
    assert "[[fig:" not in texto, "quedó un marcador de figura sin sustituir"
    assert texto.count("data:image/png;base64,") >= len(_referencias(idioma))
    secciones = [c.source.splitlines()[0] for c in nb.cells if c.cell_type == "markdown" and c.source.startswith("## ")]
    assert len(secciones) == 17


@pytest.mark.parametrize("idioma", IDIOMAS)
def test_celdas_de_codigo_ocultas_con_titulo(idioma):
    for c in _nb(idioma).cells:
        if c.cell_type == "code":
            assert c.source.startswith("# @title "), c.source[:60]
            assert c.metadata.get("cellView") == "form"
            assert c.metadata.get("jupyter", {}).get("source_hidden") is True


class _SinTextos(ast.NodeTransformer):
    """Sustituye cada texto por un marcador: dos celdas «iguales salvo el idioma» dan el mismo árbol."""

    def visit_Constant(self, node):
        if isinstance(node.value, str):
            return ast.copy_location(ast.Constant("·"), node)
        return node

    def visit_JoinedStr(self, node):
        return ast.copy_location(ast.Constant("·"), node)


def _esqueleto(fuente):
    lineas = [l for l in fuente.splitlines() if not l.lstrip().startswith(("%", "!"))]
    arbol = _SinTextos().visit(ast.parse("\n".join(lineas)))
    return ast.dump(arbol, annotate_fields=False)


@pytest.mark.skipif("en" not in IDIOMAS, reason="cuaderno en inglés aún no construido")
def test_mismo_codigo_en_los_dos_idiomas():
    es = [c.source for c in _nb("es").cells if c.cell_type == "code"]
    en = [c.source for c in _nb("en").cells if c.cell_type == "code"]
    assert len(es) == len(en)
    for k, (a, b) in enumerate(zip(es, en)):
        assert _esqueleto(a) == _esqueleto(b), f"la celda de código {k} difiere más allá de los textos:\n{a[:120]}\n---\n{b[:120]}"


_ESPANOL = re.compile(r"[áéíóúñ¿¡]|\b(el|la|los|las|del|con|para|una|que|cuando|reacción|enzima)\b", re.I)
_PERMITIDO = re.compile(r"Valentínová|Šimčíková|Cárdenas|Viñuela|Wolfenden|Glucoquinasa\.ipynb|cinetica-enzimatica|"
                        r"Cinetica_Enzimatica|enzimas|Leonardo|Juvenal|Yosa|Sols|Salas|Larion|Kamata|Grimsby|Sayed|Scruel|"
                        r"cin\.|datos\.|viz\.|visor3d\.|interactivo\.|rapido|completo|valor\(|ref\b")


@pytest.mark.skipif("en" not in IDIOMAS, reason="cuaderno en inglés aún no construido")
def test_ingles_sin_restos_de_espanol():
    restos = []
    for c in _nb("en").cells:
        texto = re.sub(r"data:image/png;base64,[A-Za-z0-9+/=]+", "", c.source)
        if c.cell_type == "code":          # en el código solo cuentan el título y los textos que ve el estudiante
            texto = "\n".join(l for l in texto.splitlines() if l.startswith("# @title") or '"' in l or "'" in l)
        for linea in texto.splitlines():
            limpia = _PERMITIDO.sub("", linea)
            if _ESPANOL.search(limpia):
                restos.append(linea.strip()[:140])
    assert not restos, "texto en español en el cuaderno en inglés:\n" + "\n".join(restos[:40])
