#!/usr/bin/env python
"""Completa la etapa en agua: IRC de MOPAC ya calculado (job.xyz), frecuencias numéricas (ASE) y descenso desde el TS."""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from enzimas.datos import ruta  # noqa: E402
from enzimas.glucoquinasa import ModeloGlucoquinasa  # noqa: E402
from enzimas.qmmm_mopac import parse_mopac_xyz_trajectory  # noqa: E402

OUT = ROOT / "data" / "precalculado" / "qmmm"
WORK = ROOT / "work" / "qmmm"
MOPAC = "/opt/anaconda3/envs/protliginteract-full/bin/mopac"
EPS = 78.4


def log(*a):
    print("[03c]", *a, flush=True)


def save_frames(path, m, frames, energies):
    Path(path).write_text("".join(m.xyz_text(f, f"frame={k} E_kcal={e:.4f}") for k, (f, e) in enumerate(zip(frames, energies))))


def main():
    t0 = time.time()
    m = ModeloGlucoquinasa(ruta("md/minimizado_completo.pdb"), ruta("sistema/complejo.prmtop"), mopac_exe=MOPAC, workdir=WORK, threads=2)
    summary = json.loads((OUT / "resumen.json").read_text(encoding="utf-8"))
    agua = summary["etapas"]["agua"]
    Ra = agua["E_reactivo"]
    ts_xyz = np.array([[float(v) for v in l[30:54].split()] for l in (OUT / "agua_ts.pdb").read_text().splitlines() if l.startswith("HETATM")])
    assert ts_xyz.shape == (m.n_qm, 3)
    # 1. IRC nativo: cuadros escritos por MOPAC en job.xyz
    try:
        fr = parse_mopac_xyz_trajectory((WORK / "agua/irc/job.xyz").read_text(errors="replace"), expected_atoms=m.n_qm)
        coords = np.asarray(fr.get("coords_ang", fr.get("coords", [])))
        if "energy_hartree" in fr:
            en = np.asarray(fr["energy_hartree"], dtype=float) * 627.5094740631
        else:
            en = np.full(len(coords), np.nan)
        if coords.ndim == 3 and len(coords) > 1:
            xis = [m.reaction_coordinate(c)[0] for c in coords]
            pd.DataFrame(dict(cuadro=range(len(coords)), xi=xis, energia_kcal=en, energia_rel_kcal=en - Ra)).to_csv(OUT / "agua_irc.csv", index=False)
            save_frames(OUT / "agua_irc.xyz", m, coords, en)
            agua["irc_cuadros"] = int(len(coords))
            agua["irc_nota"] = "IRC=1* de MOPAC (COSMO); MOPAC terminó sin el mensaje final normal pero escribió los cuadros"
            log("IRC: %d cuadros, xi de %.2f a %.2f" % (len(coords), min(xis), max(xis)))
        else:
            agua["irc_cuadros"] = 0
            log("IRC sin cuadros; claves:", list(fr.keys()))
    except Exception:
        agua["irc_error"] = traceback.format_exc()[-300:]
    if "--solo-irc" in sys.argv:
        summary["etapas"]["agua"] = agua
        (OUT / "resumen.json").write_text(json.dumps(summary, indent=1, default=float, ensure_ascii=False) + "\n", encoding="utf-8")
        log("solo IRC: listo")
        return
    # 2. frecuencias numéricas (ASE) en agua
    try:
        vib = m.vibrations(ts_xyz, tag="agua/frecuencias_ase", embedding=False, eps=EPS)
        order = np.argsort(vib["freq_cm_signed"])
        agua["frecuencias_ase_mas_bajas"] = [float(v) for v in vib["freq_cm_signed"][order][:6]]
        agua["validacion_ts_ase"] = vib["validation"]
        mode = vib["modes"][int(order[0])]
        np.savez_compressed(OUT / "agua_frecuencias_ts_ase.npz", freq_cm_signed=vib["freq_cm_signed"], modo_imaginario=mode,
                            ts_xyz=ts_xyz, simbolos=np.array(m.symbols))
        # 3. descenso desde el TS en agua
        paths = {}
        for name, sign in (("a", -1), ("b", +1)):
            paths[name] = m.descend(ts_xyz, mode, f"agua/descenso_{name}", direction=sign, steps=400, embedding=False, eps=EPS)
            fin = m.optimize(paths[name]["frames"][-1], f"agua/descenso_{name}/minimo", fmax_kcal_a=0.5, steps=600, embedding=False, eps=EPS)
            paths[name]["frames"].append(fin["coords"])
            paths[name]["energies_kcal"] = np.append(paths[name]["energies_kcal"], fin["energy_kcal"])
            paths[name]["xi"].append(fin["xi"])
        a, b = paths["a"], paths["b"]
        if a["xi"][-1] > b["xi"][-1]:
            a, b = b, a
        frames = a["frames"][::-1] + [ts_xyz] + b["frames"]
        energies = list(a["energies_kcal"][::-1]) + [agua["E_ts"]] + list(b["energies_kcal"])
        xis = [m.reaction_coordinate(f)[0] for f in frames]
        pd.DataFrame(dict(cuadro=range(len(frames)), xi=xis, energia_kcal=energies, energia_rel_kcal=np.array(energies) - Ra)).to_csv(OUT / "agua_camino_descenso.csv", index=False)
        save_frames(OUT / "agua_camino_descenso.xyz", m, frames, energies)
        agua["descenso"] = dict(E_extremo_reactivo=float(energies[0]), E_extremo_producto=float(energies[-1]), xi_extremos=[float(xis[0]), float(xis[-1])])
        # referencia: el mínimo del lado del reactivo más bajo
        if energies[0] < Ra - 0.05:
            agua["E_reactivo_original"] = Ra
            agua["E_reactivo"] = float(energies[0])
            agua["barrera_kcal"] = agua["E_ts"] - agua["E_reactivo"]
            agua["dE_reaccion"] = agua["E_producto"] - agua["E_reactivo"]
            log("reactivo en agua más bajo por el descenso: barrera ahora %.2f" % agua["barrera_kcal"])
        log("frecuencias ASE: %s; imaginarias %d" % (np.round(agua["frecuencias_ase_mas_bajas"], 1), vib["validation"]["imaginary_mode_count"]))
    except Exception:
        agua["validacion_error"] = traceback.format_exc()[-400:]
        log("frecuencias/descenso fallaron")
    summary["etapas"]["agua"] = agua
    summary.get("errores", {}).pop("agua", None)
    (OUT / "resumen.json").write_text(json.dumps(summary, indent=1, default=float, ensure_ascii=False) + "\n", encoding="utf-8")
    log("listo en %.0f min" % ((time.time() - t0) / 60))


if __name__ == "__main__":
    main()
