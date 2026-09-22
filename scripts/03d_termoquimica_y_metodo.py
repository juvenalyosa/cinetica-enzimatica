#!/usr/bin/env python
"""Mejoras del cálculo QM/MM en la enzima:

1. Termoquímica armónica: frecuencias del reactivo (las del TS ya existen) -> ZPE, H, S, G
   vibracionales -> ΔG‡(298 K) y ΔH‡, ΔS‡ (subconjunto de átomos libres, entorno fijo).
2. Sensibilidad al hamiltoniano: energías de punto único PM6-D3H4 (y PM7) en las geometrías
   PM7 de R, TS y P, con el mismo campo y el mismo término LJ.
Actualiza data/precalculado/qmmm/resumen.json (claves "termoquimica" y "metodos").
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from enzimas.datos import ruta  # noqa: E402
from enzimas.glucoquinasa import ModeloGlucoquinasa  # noqa: E402

OUT = ROOT / "data" / "precalculado" / "qmmm"
WORK = ROOT / "work" / "qmmm"
MOPAC = "/opt/anaconda3/envs/protliginteract-full/bin/mopac"


def log(*a):
    print("[03d]", *a, flush=True)


def pdb_xyz(path):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in Path(path).read_text().splitlines() if l.startswith("HETATM")])


def main():
    t0 = time.time()
    m = ModeloGlucoquinasa(ruta("md/minimizado_completo.pdb"), ruta("sistema/complejo.prmtop"), mopac_exe=MOPAC, workdir=WORK, threads=2)
    summary = json.loads((OUT / "resumen.json").read_text(encoding="utf-8"))
    R, TS, P = pdb_xyz(OUT / "reactivo.pdb"), pdb_xyz(OUT / "ts.pdb"), pdb_xyz(OUT / "producto.pdb")
    # 1. termoquímica
    vib_ts = np.load(OUT / "frecuencias_ts.npz")["freq_cm_signed"]
    vib_r = m.vibrations(R, tag="frecuencias_reactivo")
    np.savez_compressed(OUT / "frecuencias_reactivo.npz", freq_cm_signed=vib_r["freq_cm_signed"], simbolos=np.array(m.symbols))
    termo = {}
    for T in (298.15, 303.15, 310.15):
        tr, tt = m.harmonic_thermo(vib_r["freq_cm_signed"], T), m.harmonic_thermo(vib_ts, T)
        dE = summary["etapas"]["dimero"]["barrera_kcal"]
        termo[str(T)] = dict(
            dE_kcal=dE, dZPE_kcal=tt["zpe_kcal"] - tr["zpe_kcal"], dH_vib_kcal=tt["h_vib_kcal"] - tr["h_vib_kcal"],
            dS_vib_cal=tt["s_vib_cal"] - tr["s_vib_cal"],
            dH_kcal=dE + tt["h_vib_kcal"] - tr["h_vib_kcal"],
            dG_kcal=dE + tt["g_vib_kcal"] - tr["g_vib_kcal"],
            reactivo=tr, ts=tt,
        )
        log("T = %.2f K: ΔE‡ %.2f, ΔZPE %.2f, ΔH‡ %.2f, ΔS‡ %.2f cal/mol/K, ΔG‡ %.2f kcal/mol" % (
            T, dE, termo[str(T)]["dZPE_kcal"], termo[str(T)]["dH_kcal"], termo[str(T)]["dS_vib_cal"], termo[str(T)]["dG_kcal"]))
    termo["nota"] = ("Aproximación armónica sobre los 70 átomos libres de la región QM con el entorno fijo; modos < 50 cm-1 "
                     "elevados a 50 cm-1; sin contribuciones del entorno ni de la conformación de la proteína.")
    termo["reactivo_n_imag"] = int((vib_r["freq_cm_signed"] < 0).sum())
    termo["reactivo_frecuencias_mas_bajas"] = [float(v) for v in np.sort(vib_r["freq_cm_signed"])[:4]]
    summary["etapas"]["termoquimica"] = termo
    # 2. sensibilidad al método
    metodos = {}
    for method in ("PM7", "PM6-D3H4"):
        m.method = method
        e = {}
        for name, xyz in (("R", R), ("TS", TS), ("P", P)):
            e[name] = m.energy_gradient(xyz, f"metodo_{method}/{name}")[0]
        metodos[method] = dict(barrera_kcal=e["TS"] - e["R"], dE_reaccion_kcal=e["P"] - e["R"], **{f"E_{k}": v for k, v in e.items()})
        log("%s//PM7: barrera %.2f, ΔE %.2f" % (method, e["TS"] - e["R"], e["P"] - e["R"]))
    m.method = "PM7"
    metodos["nota"] = "Energías de punto único en las geometrías PM7 (R, TS, P), mismo campo MM y mismo término LJ."
    summary["etapas"]["metodos"] = metodos
    (OUT / "resumen.json").write_text(json.dumps(summary, indent=1, default=float, ensure_ascii=False) + "\n", encoding="utf-8")
    log("listo en %.0f min" % ((time.time() - t0) / 60))


if __name__ == "__main__":
    main()
