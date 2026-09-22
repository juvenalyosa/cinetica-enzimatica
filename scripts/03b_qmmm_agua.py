#!/usr/bin/env python
"""Repite solo la etapa 8 (clúster catalítico en agua COSMO) reutilizando R, P y TS de la corrida QM/MM.

Estrategia robusta para el TS en agua:
  1. SADDLE (QST2) entre R y P  ->  estimación
  2. TS nativo desde la geometría SADDLE (RECALC=1); si falla, refinamiento con el método del dímero (ASE)
  3. FORCETS (una frecuencia imaginaria)  4. IRC=1*
Actualiza data/precalculado/qmmm/resumen.json (clave "agua") y los archivos agua_*.
"""
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

OUT = ROOT / "data" / "precalculado" / "qmmm"
WORK = ROOT / "work" / "qmmm"
MOPAC = "/opt/anaconda3/envs/protliginteract-full/bin/mopac"
EPS = 78.4


def log(*a):
    print("[03b]", *a, flush=True)


def save_frames(path, m, frames, energies):
    Path(path).write_text("".join(m.xyz_text(f, f"frame={k} E_kcal={e:.4f}") for k, (f, e) in enumerate(zip(frames, energies))))


def main():
    t0 = time.time()
    m = ModeloGlucoquinasa(ruta("md/minimizado_completo.pdb"), ruta("sistema/complejo.prmtop"), mopac_exe=MOPAC, workdir=WORK, threads=2)
    summary = json.loads((OUT / "resumen.json").read_text(encoding="utf-8"))
    R, P, TS = np.load(WORK / "R.npy"), np.load(WORK / "P.npy"), np.load(WORK / "TS.npy")
    if (WORK / "reactivo_desde_descenso" / "final").exists():
        pass
    agua = dict(eps=EPS, metodo="PM7/COSMO (agua implícita), anclajes fijos, sin el campo de la enzima")
    try:
        Ra = m.native_optimize(R, "agua/reactivo", eps=EPS)
        Pa = m.native_optimize(P, "agua/producto", eps=EPS)
        agua.update(E_reactivo=Ra["energy_kcal"], E_producto=Pa["energy_kcal"], dE_reaccion=Pa["energy_kcal"] - Ra["energy_kcal"])
        m.write_pdb(OUT / "agua_reactivo.pdb", Ra["coords"])
        m.write_pdb(OUT / "agua_producto.pdb", Pa["coords"])
        Sa = m.native_saddle(Ra["coords"], Pa["coords"], "agua/saddle_qst2", eps=EPS)
        agua["E_saddle_qst2"] = Sa["energy_kcal"]
        agua["barrera_saddle_qst2_kcal"] = Sa["energy_kcal"] - Ra["energy_kcal"]
        agua.update({f"saddle_{k}": v for k, v in m.reaction_coordinates(Sa["coords"]).items()})
        m.write_pdb(OUT / "agua_saddle_qst2.pdb", Sa["coords"])
        ts_xyz, ts_e, metodo_ts = None, None, None
        for guess_name, guess in (("saddle", Sa["coords"]), ("qmmm_ts", TS)):
            try:
                Ta = m.run(guess, f"agua/ts_desde_{guess_name}", "TS GNORM=1.0 CYCLES=2000 LET RECALC=1 DDMIN=0.0", embedding=False, eps=EPS)
                ts_xyz, ts_e, metodo_ts = Ta["coords"], Ta["energy_kcal"], f"MOPAC TS desde {guess_name}"
                log(f"TS nativo desde {guess_name}: E = {ts_e:.2f}")
                break
            except Exception as exc:
                log(f"TS nativo desde {guess_name} falló: {str(exc)[:120]}")
        if ts_xyz is None:
            tangent = Pa["coords"] - Ra["coords"]
            tangent[sorted(m.fixed)] = 0.0
            D = m.dimer(Sa["coords"], tag="agua/dimero", fmax_kcal_a=0.5, steps=400, mode_guess=tangent, embedding=False, eps=EPS)
            ts_xyz, ts_e, metodo_ts = D["coords"], D["energy_kcal"], "método del dímero (ASE) desde SADDLE"
        agua.update(E_ts=ts_e, barrera_kcal=ts_e - Ra["energy_kcal"], metodo_ts=metodo_ts)
        agua.update({f"ts_{k}": v for k, v in m.reaction_coordinates(ts_xyz).items()})
        m.write_pdb(OUT / "agua_ts.pdb", ts_xyz)
        Fa = m.native_frequencies(ts_xyz, "agua/forcets", eps=EPS)
        agua["validacion_ts"] = Fa["validation"]
        agua["frecuencias_mas_bajas"] = [float(v) for v in np.sort(Fa["freq_cm_signed"])[:6]]
        if Fa["normal_modes_cart"] is not None:
            np.savez_compressed(OUT / "agua_frecuencias_ts.npz", freq_cm_signed=Fa["freq_cm_signed"], modos=Fa["normal_modes_cart"],
                                ts_xyz=ts_xyz, simbolos=np.array(m.symbols))
        Ia = m.native_irc(ts_xyz, "agua/irc", eps=EPS)
        fr = Ia["irc_frames"]
        coords = np.asarray(fr.get("coords", [])) if isinstance(fr, dict) else np.zeros((0,))
        if coords.ndim == 3 and len(coords) > 1:
            en_key = next((k for k in ("energies_kcal", "energy_kcal", "energies", "heat_kcal") if k in fr), None)
            en = np.asarray(fr[en_key], dtype=float) if en_key else np.full(len(coords), np.nan)
            xis = [m.reaction_coordinate(c)[0] for c in coords]
            pd.DataFrame(dict(cuadro=range(len(coords)), xi=xis, energia_kcal=en, energia_rel_kcal=en - Ra["energy_kcal"])).to_csv(OUT / "agua_irc.csv", index=False)
            save_frames(OUT / "agua_irc.xyz", m, coords, en)
            agua["irc_cuadros"] = int(len(coords))
        else:
            agua["irc_cuadros"] = 0
            agua["irc_claves"] = list(fr.keys()) if isinstance(fr, dict) else str(type(fr))
        log("agua: barrera %.2f kcal/mol (%s); dE %.2f; imaginarias %d" % (agua["barrera_kcal"], metodo_ts, agua["dE_reaccion"], Fa["validation"]["imaginary_mode_count"]))
    except Exception as exc:
        summary.setdefault("errores", {})["agua"] = traceback.format_exc()
        log("falló:", exc)
    agua["segundos"] = time.time() - t0
    summary["etapas"]["agua"] = agua
    summary.get("errores", {}).pop("agua", None) if "barrera_kcal" in agua else None
    (OUT / "resumen.json").write_text(json.dumps(summary, indent=1, default=float, ensure_ascii=False) + "\n")
    log("listo en %.0f min" % ((time.time() - t0) / 60))


if __name__ == "__main__":
    main()
