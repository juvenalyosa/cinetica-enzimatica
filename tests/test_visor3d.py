"""Los visores 3D se construyen con los datos precalculados y usan enlaces reales, no distancias."""
from __future__ import annotations

import json
import re

import numpy as np
import pytest

pytest.importorskip("IPython")
from enzimas import datos, visor3d  # noqa: E402


def _config(html):
    m = re.search(r"const cfg=(\{.*?\});\s*\nconst URL", html, re.S)
    return json.loads(m.group(1).replace("<\\/", "</"))


@pytest.mark.parametrize("fabrica", [visor3d.pelicula_reaccion, lambda: visor3d.pelicula_reaccion("agua"),
                                     visor3d.modo_imaginario, visor3d.pelicula_md, visor3d.complejo_cristal,
                                     visor3d.region_qm])
def test_visores_generan_html(fabrica):
    html = fabrica().data
    assert "3Dmol" in html and "const cfg=" in html


def test_reaccion_sincroniza_energia_y_fotogramas():
    cfg = _config(visor3d.pelicula_reaccion().data)
    n = len(cfg["grafica"]["x"])
    assert n == len(cfg["dinamicos"]) == len(cfg["fases"]) == len(cfg["lecturas"]["valores"])
    modelo = cfg["modelos"][0]
    import base64
    assert len(base64.b64decode(modelo["frames"])) == 2 * 3 * modelo["n"] * n
    # la barrera del perfil sincronizado es la del camino calculado (19.1 kcal/mol)
    assert max(cfg["grafica"]["y"]) == pytest.approx(19.1, abs=0.3)


def test_enlaces_de_la_topologia_no_por_distancia():
    topo = datos.json_("qmmm/topologia_qm.json")
    rompen = {tuple(sorted((d["i"], d["j"]))) for d in topo["dinamicos"]}
    fijos = {tuple(sorted(e)) for e in topo["enlaces"]}
    assert not rompen & fijos                         # los enlaces que cambian no están entre los fijos
    mg = topo["mg"]
    assert not any(mg in e for e in fijos)            # el Mg²⁺ no tiene enlaces covalentes
    cfg = _config(visor3d.pelicula_reaccion().data)
    visibles = [k for k, a in enumerate(topo["atomos"]) if a["nombre"] != "HL"]
    nuevo = {k: i for i, k in enumerate(visibles)}
    d = topo["destacados"]
    rompe = tuple(sorted((nuevo[d["PG"]], nuevo[d["O3B"]])))
    forma = tuple(sorted((nuevo[d["PG"]], nuevo[d["O6"]])))

    def orden(f):
        return {tuple(sorted(x[:2])): x[2] for x in cfg["dinamicos"][f] if x[2] >= 0}

    # reactivo: Pγ–O3β casi entero, Pγ–O6 aún no existe; estado de transición: los dos a medias
    assert orden(0)[rompe] > 0.85 and orden(0).get(forma, 0) < 0.2
    ts = orden(cfg["grafica"]["ts"])
    assert 0.25 < ts[rompe] < 0.65 and 0.25 < ts[forma] < 0.65
