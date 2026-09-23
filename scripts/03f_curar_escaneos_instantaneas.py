#!/usr/bin/env python
"""Extrae los perfiles de energía de los escaneos de las instantáneas de MD (script 03e) a un CSV ligero.

03e solo guardó el último punto de cada escaneo en instantaneas.json; aquí se relee cada punto desde las
trayectorias de ASE en work/qmmm_instantaneas/<t>/escaneo/punto_NN/opt.traj y se resta la restricción
armónica (k = 300 kcal/mol/Å², centrada en el ξ objetivo) para obtener la energía sin sesgo, igual que
``GlucokinaseQMMM.scan``. Comprobación: el último punto reproduce ``escaneo_max_rel_kcal`` y ``xi_max``.

    python scripts/03f_curar_escaneos_instantaneas.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work" / "qmmm_instantaneas"
OUT = ROOT / "data" / "precalculado" / "qmmm"
EV_TO_KCAL = 23.060548
K_RESTRICCION = 300.0


def indices_reaccion():
    lineas = [l for l in (OUT / "region_qm_inicial.pdb").read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    nombres = {(l[17:20].strip(), l[12:16].strip()): i for i, l in enumerate(lineas)}
    return nombres[("ATP", "PG")], nombres[("ATP", "O3B")], nombres[("GLC", "O6")]


def main():
    pg, o3b, o6 = indices_reaccion()
    xi_de = lambda p: np.linalg.norm(p[pg] - p[o3b]) - np.linalg.norm(p[pg] - p[o6])
    inst = json.loads((OUT / "instantaneas.json").read_text(encoding="utf-8"))["instantaneas"]
    filas = []
    for r in inst:
        if not r.get("sin_producto_estable"):
            continue
        carpeta = WORK / f"t{int(r['t_ps']):04d}"
        reactivo = read(carpeta / "reactivo" / "opt.traj", -1)
        e_r, xi_r = reactivo.get_potential_energy(), xi_de(reactivo.positions)
        puntos = sorted((carpeta / "escaneo").glob("punto_*"))
        objetivos = np.linspace(xi_r, 2.0, len(puntos))
        for k, (d, xi0) in enumerate(zip(puntos, objetivos)):
            a = read(d / "opt.traj", -1)
            xi = xi_de(a.positions)
            e = (a.get_potential_energy() - e_r) * EV_TO_KCAL - 0.5 * K_RESTRICCION * (xi - xi0) ** 2
            filas.append(dict(instantanea=r["instantanea"], t_ps=r["t_ps"], punto=k, xi=round(float(xi), 4),
                              energia_rel_kcal=round(float(e), 3)))
        assert abs(filas[-1]["energia_rel_kcal"] - r["escaneo_max_rel_kcal"]) < 0.01, r["instantanea"]
        assert abs(filas[-1]["xi"] - r["xi_max"]) < 0.01, r["instantanea"]
    with open(OUT / "instantaneas_escaneos.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0]))
        w.writeheader()
        w.writerows(filas)
    print(f"{len(filas)} puntos → {OUT / 'instantaneas_escaneos.csv'}")


if __name__ == "__main__":
    main()
