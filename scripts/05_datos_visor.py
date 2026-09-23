#!/usr/bin/env python
"""Prepara los datos ligeros de los visores 3D del cuaderno (enzimas/visor3d.py).

1. ``qmmm/topologia_qm.json``: la conectividad REAL de la región QM, tomada de la topología de Amber
   (complejo_vacio.prmtop) y no adivinada por distancias. Separa
     * enlaces fijos (no cambian en la reacción),
     * enlaces que se rompen o se forman (Pγ–O3β, Pγ–O6, O6–H, H–OD Asp205): el visor los dibuja con un
       grosor que sigue a la distancia de cada fotograma,
     * coordinación del Mg²⁺ (no covalente: línea fina discontinua).
2. ``md/pelicula_md.json.gz``: la dinámica molecular (100 fotogramas, 10 ps entre ellos) alineada sobre los Cα
   al cristal minimizado (el mismo marco de coordenadas que el QM/MM), con el esqueleto completo de la proteína
   (para la caricatura), el sitio activo con todos sus átomos pesados y los ligandos completos.

    python scripts/05_datos_visor.py          (requiere parmed y mdtraj; en Colab solo se leen los resultados)
"""
from __future__ import annotations

import gzip
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRE = ROOT / "data" / "precalculado"
DESPLAZAMIENTO = 4          # numeración de tleap = cristal − 4 (Asp201 de tleap = Asp205 del cristal)


def _prmtop():
    import parmed

    tmp = tempfile.NamedTemporaryFile(suffix=".prmtop", delete=False)
    with gzip.open(PRE / "sistema" / "complejo_vacio.prmtop.gz", "rb") as f:
        shutil.copyfileobj(f, tmp)
    tmp.close()
    return parmed.load_file(tmp.name)


def _nombre_residuo(resname, resseq):
    if resname in ("GLC", "ATP", "MG", "K+"):
        return {"GLC": "glucosa", "ATP": "ATP", "MG": "Mg²⁺", "K+": "K⁺"}[resname]
    return f"{resname.capitalize()}{resseq + DESPLAZAMIENTO}"


def topologia_qm(prm):
    part = json.loads((PRE / "qmmm" / "particion.json").read_text())
    qm = part["qm_atoms"]
    glob = {a["global_index"]: a["k"] for a in qm if a.get("global_index") is not None and a["resname"] != "HL"}
    atomos = []
    for a in qm:
        atomos.append(dict(nombre=a["name"], elemento=a["element"], residuo=a["resname"], resseq=a["resseq"],
                           etiqueta=_nombre_residuo(a["resname"], a["resseq"])))
    enlaces = set()
    for b in prm.bonds:
        i, j = b.atom1.idx, b.atom2.idx
        if i in glob and j in glob:
            enlaces.add(tuple(sorted((glob[i], glob[j]))))
    # átomos de enlace (H que tapan el enlace cortado en la frontera QM/MM)
    hl = [a["k"] for a in qm if a["name"] == "HL"]
    for link, k in zip(part["links"], hl):
        enlaces.add(tuple(sorted((glob[link["qm"]], k))))

    def idx(res, name):
        return next(a["k"] for a in qm if a["resname"] == res and a["name"] == name)

    pg, o3b, o6 = idx("ATP", "PG"), idx("ATP", "O3B"), idx("GLC", "O6")
    vecinos_o6 = [j if i == o6 else i for i, j in enlaces if o6 in (i, j)]
    h6 = next(k for k in vecinos_o6 if qm[k]["element"] == "H")
    # el oxígeno de Asp205 que recibe el protón: el más cercano al H6 en el producto
    prod = _leer_pdb_xyz(PRE / "qmmm" / "producto.pdb")
    asp_o = [a["k"] for a in qm if a["resname"] == "ASP" and a["name"] in ("OD1", "OD2")]
    od = min(asp_o, key=lambda k: np.linalg.norm(prod[k] - prod[h6]))
    dinamicos = [
        dict(i=pg, j=o3b, nombre="Pγ–O3β", papel="se rompe", r0=1.63),
        dict(i=pg, j=o6, nombre="Pγ–O6", papel="se forma", r0=1.63),
        dict(i=o6, j=h6, nombre="O6–H", papel="se rompe", r0=0.97),
        dict(i=h6, j=od, nombre="H–O(Asp205)", papel="se forma", r0=0.97),
    ]
    for d in dinamicos:
        enlaces.discard(tuple(sorted((d["i"], d["j"]))))
    mg = idx("MG", "MG")
    reac = _leer_pdb_xyz(PRE / "qmmm" / "reactivo.pdb")
    coord = [k for k in range(len(qm)) if k != mg and qm[k]["element"] == "O" and np.linalg.norm(reac[k] - reac[mg]) < 2.5]
    destacados = dict(PG=pg, O3B=o3b, O6=o6, H6=h6, OD=od, MG=mg)
    return dict(atomos=atomos, enlaces=sorted(map(list, enlaces)), dinamicos=dinamicos, mg=mg, coordinacion_mg=coord,
                destacados=destacados,
                nota="Enlaces de la topología Amber (complejo_vacio.prmtop); los 'dinamicos' cambian en la reacción.")


def _leer_pdb_xyz(path):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in Path(path).read_text().splitlines() if l.startswith(("ATOM", "HETATM"))])


def pelicula_md(prm, n_max=100):
    import mdtraj as md

    t = md.load(str(PRE / "md" / "trayectoria_soluto.dcd"), top=str(PRE / "md" / "minimizado_soluto.pdb"))
    ref = md.load(str(PRE / "md" / "minimizado_soluto.pdb"))
    top = t.topology
    ca = top.select("protein and name CA")
    t.superpose(ref, atom_indices=ca)
    for a, b in zip(top.atoms, prm.atoms):            # misma numeración que la topología de Amber
        assert a.element.symbol == b.element_name and a.residue.index == b.residue.idx, (a, b)
    lig = top.select("resname GLC ATP MG")
    cerca = md.compute_neighbors(t[0], 0.65, lig)[0]
    res_sitio = sorted({top.atom(i).residue.index for i in cerca if top.atom(i).residue.is_protein})
    sel = set(top.select("protein and name N CA C O"))
    for r in res_sitio:
        sel |= {a.index for a in top.residue(r).atoms if a.element.symbol != "H" or a.name in ("HZ1", "HZ2", "HZ3", "HG1", "H")}
    sel |= {a.index for a in top.atoms if a.residue.name in ("GLC", "ATP", "MG")}
    sel = sorted(sel)
    pos = {g: k for k, g in enumerate(sel)}
    atomos = []
    for g in sel:
        a = top.atom(g)
        atomos.append(dict(n=a.name, e=a.element.symbol, r=a.residue.name, s=a.residue.resSeq + DESPLAZAMIENTO,
                           h=0 if a.residue.is_protein else 1, sitio=1 if a.residue.index in res_sitio or not a.residue.is_protein else 0))
    enlaces = []
    for b in prm.bonds:
        i, j = b.atom1.idx, b.atom2.idx
        if i in pos and j in pos:
            enlaces.append([pos[i], pos[j]])
    paso = max(1, t.n_frames // n_max)
    frames = [np.round(t.xyz[f, sel] * 10.0, 2).ravel().tolist() for f in range(0, t.n_frames, paso)]
    import pandas as pd

    # el DCD no guarda el tiempo real; analisis.csv tiene una fila por fotograma (cada 10 ps)
    ana = pd.read_csv(PRE / "md" / "analisis.csv")
    assert len(ana) == t.n_frames, (len(ana), t.n_frames)
    idx = list(range(0, t.n_frames, paso))
    tiempo = ana["tiempo_ps"].to_numpy()[idx].round(1).tolist()
    d = ana["d_PG_O6_A"].to_numpy()[idx].round(3).tolist()
    # comprobación: la distancia medida en estos fotogramas coincide con la del análisis
    i_pg = next(a.index for a in top.atoms if a.residue.name == "ATP" and a.name == "PG")
    i_o6 = next(a.index for a in top.atoms if a.residue.name == "GLC" and a.name == "O6")
    medida = np.linalg.norm(t.xyz[idx, i_pg] - t.xyz[idx, i_o6], axis=1) * 10
    assert np.allclose(medida, d, atol=0.05), np.abs(medida - d).max()
    return dict(atomos=atomos, enlaces=enlaces, frames=frames, tiempo_ps=tiempo, d_PG_O6=d,
                residuos_sitio=[_nombre_residuo(top.residue(r).name, top.residue(r).resSeq) for r in res_sitio],
                nota="Alineada sobre los Cα al cristal minimizado; coordenadas en Å con 2 decimales.")


def main():
    prm = _prmtop()
    topo = topologia_qm(prm)
    (PRE / "qmmm" / "topologia_qm.json").write_text(json.dumps(topo, ensure_ascii=False, indent=0) + "\n", encoding="utf-8")
    print("topología QM:", len(topo["enlaces"]), "enlaces fijos,", len(topo["dinamicos"]), "dinámicos; Mg coordinado a",
          len(topo["coordinacion_mg"]), "O")
    peli = pelicula_md(prm)
    salida = PRE / "md" / "pelicula_md.json.gz"
    with gzip.GzipFile(salida, "wb", mtime=0) as f:            # mtime=0: archivo reproducible byte a byte
        f.write(json.dumps(peli, separators=(",", ":"), ensure_ascii=False).encode("utf-8"))
    print("película MD:", len(peli["atomos"]), "átomos ×", len(peli["frames"]), "fotogramas →",
          f"{salida.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
