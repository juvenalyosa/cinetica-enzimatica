#!/usr/bin/env python
"""Camino de reacción QM/MM (PM7/ff14SB) de la glucoquinasa: glucosa + ATP -> G6P + ADP.

Etapas (todas con ``enzimas.glucoquinasa.ModeloGlucoquinasa``):
  1. reactivo    optimización L-BFGS del complejo de Michaelis (anclajes fijos)
  2. producto    construcción geométrica + optimización restringida y libre
  3. escaneo     perfil relajado sobre xi = d(PG-O3B) - d(PG-O6) (restricción armónica)
  4. NEB         banda elástica con imagen trepadora entre R y P (imágenes del escaneo)
  5. dímero      refinamiento del punto de silla (solo gradientes)
  6. frecuencias Hessiano numérico sobre los átomos libres: UNA frecuencia imaginaria
  7. descenso    camino de mínima energía descendente desde el TS hacia R y P
  8. agua        el mismo clúster QM fuera de la enzima (COSMO), con MOPAC nativo:
                 SADDLE (= QST2), TS, FORCETS, IRC  -> barrera "sin el resto de la enzima"

Productos en data/precalculado/qmmm/ (ligeros, para el notebook en modo rápido).
Uso:  python scripts/03_qmmm_reaccion.py [--estructura PDB] [--rapido]
"""
from __future__ import annotations

import gzip
import json
import shutil
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
QUICK = "--rapido" in sys.argv
PRMTOP = ruta("sistema/complejo.prmtop")


def log(*a):
    print("[03]", *a, flush=True)


def structure_path():
    if "--estructura" in sys.argv:
        return Path(sys.argv[sys.argv.index("--estructura") + 1])
    return ruta("md/minimizado_completo.pdb")


def save_frames(path, m, frames, energies, extra=None):
    txt = ""
    for k, (f, e) in enumerate(zip(frames, energies)):
        txt += m.xyz_text(f, f"frame={k} E_kcal={e:.4f}" + ("" if extra is None else f" {extra[k]}"))
    Path(path).write_text(txt)


def main():
    t_start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    pdb = structure_path()
    m = ModeloGlucoquinasa(pdb, PRMTOP, mopac_exe=MOPAC, workdir=WORK, threads=2)
    summary = dict(estructura=str(pdb.name), particion=m.describe(), etapas={}, errores={})
    m.write_pdb(OUT / "region_qm_inicial.pdb", m.qm_xyz)
    m.write_env_pdb(OUT / "entorno_mm_8A.pdb", radius_a=8.0)
    (OUT / "particion.json").write_text(json.dumps(m.describe(), indent=1))

    # 0. comprobación del gradiente (documentada en el notebook)
    rows = m.check_gradient(m.qm_xyz, [m.O6, m.HO6, m.O1G, m.PG], h=0.005)
    pd.DataFrame(rows).to_csv(OUT / "comprobacion_gradiente.csv", index=False)
    log("gradiente: |dif| máx = %.4f kcal/mol/A" % max(abs(r["diff"]) for r in rows))

    # 1. reactivo
    t0 = time.time()
    R = m.optimize(m.qm_xyz, "reactivo", fmax_kcal_a=0.5, steps=600)
    summary["etapas"]["reactivo"] = dict(energia_kcal=R["energy_kcal"], convergido=R["converged"], pasos=R["steps"],
                                        segundos=time.time() - t0, **m.reaction_coordinates(R["coords"]),
                                        contactos=m.contacts(R["coords"]))
    np.save(WORK / "R.npy", R["coords"])
    m.write_pdb(OUT / "reactivo.pdb", R["coords"])

    # 2. producto
    t0 = time.time()
    P0 = m.build_product(R["coords"])
    P1 = m.optimize(P0, "producto_restringido", fmax_kcal_a=1.0, steps=400, restraint=(2.2, 100.0))
    P = m.optimize(P1["coords"], "producto", fmax_kcal_a=0.5, steps=600)
    summary["etapas"]["producto"] = dict(energia_kcal=P["energy_kcal"], convergido=P["converged"], pasos=P["steps"],
                                        segundos=time.time() - t0, **m.reaction_coordinates(P["coords"]),
                                        contactos=m.contacts(P["coords"]),
                                        dE_reaccion_kcal=P["energy_kcal"] - R["energy_kcal"])
    np.save(WORK / "P.npy", P["coords"])
    m.write_pdb(OUT / "producto.pdb", P["coords"])
    log("dE reacción = %.2f kcal/mol" % (P["energy_kcal"] - R["energy_kcal"]))

    # 3. escaneo relajado
    t0 = time.time()
    xi_r, xi_p = m.reaction_coordinate(R["coords"])[0], m.reaction_coordinate(P["coords"])[0]
    n_scan = 7 if QUICK else 15
    xi_values = np.linspace(xi_r, xi_p, n_scan)
    scan = m.scan(R["coords"], xi_values, tag="escaneo", k_kcal_a2=300.0, fmax_kcal_a=1.0, steps=300)
    df = pd.DataFrame([{k: v for k, v in p.items() if k != "coords"} for p in scan])
    df["energia_rel_kcal"] = df["energia_kcal"] - R["energy_kcal"]
    df.to_csv(OUT / "escaneo.csv", index=False)
    save_frames(OUT / "escaneo.xyz", m, [p["coords"] for p in scan], [p["energia_kcal"] for p in scan])
    summary["etapas"]["escaneo"] = dict(n_puntos=n_scan, segundos=time.time() - t0,
                                       barrera_aprox_kcal=float(df["energia_rel_kcal"].max()),
                                       xi_maximo=float(df.loc[df["energia_rel_kcal"].idxmax(), "xi"]))
    log("escaneo: barrera aproximada %.2f kcal/mol" % df["energia_rel_kcal"].max())

    # 4. NEB con imagen trepadora
    t0 = time.time()
    n_img = 5 if QUICK else 7
    pick = np.linspace(0, len(scan) - 1, n_img + 2)[1:-1].round().astype(int)
    images = [R["coords"]] + [scan[k]["coords"] for k in pick] + [P["coords"]]
    try:
        neb = m.neb(images, tag="neb", climb=True, fmax_kcal_a=1.0, steps=60 if QUICK else 250)
        k_ts = int(np.argmax(neb["energies_kcal"]))
        ts_guess = neb["coords"][k_ts]
        pd.DataFrame(dict(imagen=range(len(images)), xi=neb["xi"], energia_kcal=neb["energies_kcal"],
                          energia_rel_kcal=neb["energies_kcal"] - neb["energies_kcal"][0])).to_csv(OUT / "neb.csv", index=False)
        save_frames(OUT / "neb.xyz", m, neb["coords"], neb["energies_kcal"])
        summary["etapas"]["neb"] = dict(imagenes=len(images), convergido=neb["converged"], pasos=neb["steps"],
                                       segundos=time.time() - t0, imagen_ts=k_ts,
                                       barrera_kcal=float(neb["energies_kcal"][k_ts] - neb["energies_kcal"][0]))
    except Exception as exc:  # seguir con el máximo del escaneo
        summary["errores"]["neb"] = traceback.format_exc()
        log("NEB falló:", exc)
        k_max = int(df["energia_rel_kcal"].idxmax())
        ts_guess = scan[k_max]["coords"]
    np.save(WORK / "TS_guess.npy", ts_guess)

    # 5. dímero
    t0 = time.time()
    try:
        # modo inicial: tangente del camino en el máximo
        if "neb" in summary["etapas"]:
            tangent = neb["coords"][min(k_ts + 1, len(images) - 1)] - neb["coords"][max(k_ts - 1, 0)]
        else:
            tangent = scan[min(k_max + 1, len(scan) - 1)]["coords"] - scan[max(k_max - 1, 0)]["coords"]
        tangent[sorted(m.fixed)] = 0.0
        TS = m.dimer(ts_guess, tag="dimero", fmax_kcal_a=0.5, steps=100 if QUICK else 400, mode_guess=tangent)
        summary["etapas"]["dimero"] = dict(energia_kcal=TS["energy_kcal"], convergido=TS["converged"],
                                          curvatura_eV_A2=TS["curvature"], segundos=time.time() - t0,
                                          barrera_kcal=TS["energy_kcal"] - R["energy_kcal"],
                                          **m.reaction_coordinates(TS["coords"]))
        ts_xyz = TS["coords"]
        log("TS (dímero): barrera %.2f kcal/mol" % (TS["energy_kcal"] - R["energy_kcal"]))
    except Exception as exc:
        summary["errores"]["dimero"] = traceback.format_exc()
        log("dímero falló:", exc)
        ts_xyz = ts_guess
    np.save(WORK / "TS.npy", ts_xyz)
    m.write_pdb(OUT / "ts.pdb", ts_xyz)

    # 6. frecuencias
    t0 = time.time()
    try:
        vib = m.vibrations(ts_xyz, tag="frecuencias")
        order = np.argsort(vib["freq_cm_signed"])
        k_im = int(order[0])
        mode = vib["modes"][k_im]
        np.savez_compressed(OUT / "frecuencias_ts.npz", freq_cm_signed=vib["freq_cm_signed"], modo_imaginario=mode,
                            ts_xyz=ts_xyz, simbolos=np.array(m.symbols))
        summary["etapas"]["frecuencias"] = dict(segundos=time.time() - t0, validacion=vib["validation"],
                                               frecuencias_mas_bajas=[float(v) for v in vib["freq_cm_signed"][order][:6]])
        log("frecuencia imaginaria: %.1f cm-1; modos imaginarios: %d" % (vib["freq_cm_signed"][k_im], vib["validation"]["imaginary_mode_count"]))
    except Exception as exc:
        summary["errores"]["frecuencias"] = traceback.format_exc()
        log("frecuencias fallaron:", exc)
        mode = None

    # 7. descenso desde el TS (camino de mínima energía descendente)
    t0 = time.time()
    if mode is not None:
        try:
            paths = {}
            for name, sign in (("hacia_reactivo", -1), ("hacia_producto", +1)):
                # el signo se decide por la coordenada xi tras el desplazamiento
                paths[name] = m.descend(ts_xyz, mode, f"descenso_{name}", direction=sign, steps=150 if QUICK else 500)
                # el descenso con FIRE se detiene en zonas planas: se remata con L-BFGS hasta el mínimo
                fin = m.optimize(paths[name]["frames"][-1], f"descenso_{name}/minimo", fmax_kcal_a=0.5, steps=600)
                paths[name]["frames"].append(fin["coords"])
                paths[name]["energies_kcal"] = np.append(paths[name]["energies_kcal"], fin["energy_kcal"])
                paths[name]["xi"].append(fin["xi"])
            # ordenar: el que termina con xi menor es el lado del reactivo
            a, b = paths["hacia_reactivo"], paths["hacia_producto"]
            if a["xi"][-1] > b["xi"][-1]:
                a, b = b, a
            frames = a["frames"][::-1] + [ts_xyz] + b["frames"]
            energies = list(a["energies_kcal"][::-1]) + [summary["etapas"].get("dimero", {}).get("energia_kcal", np.nan)] + list(b["energies_kcal"])
            xis = [m.reaction_coordinate(f)[0] for f in frames]
            pd.DataFrame(dict(cuadro=range(len(frames)), xi=xis, energia_kcal=energies,
                              energia_rel_kcal=np.array(energies) - R["energy_kcal"])).to_csv(OUT / "camino_descenso.csv", index=False)
            save_frames(OUT / "camino_descenso.xyz", m, frames, energies)
            summary["etapas"]["descenso"] = dict(segundos=time.time() - t0, n_cuadros=len(frames),
                                                E_extremo_reactivo_rel=float(energies[0] - R["energy_kcal"]),
                                                E_extremo_producto_rel=float(energies[-1] - R["energy_kcal"]),
                                                xi_extremos=[float(xis[0]), float(xis[-1])])
            # Si el descenso encuentra un reactivo/producto más bajo que la optimización directa,
            # se reoptimiza desde ahí y se toma el mínimo más bajo como referencia.
            if energies[0] < R["energy_kcal"] - 0.05:
                R2 = m.optimize(frames[0], "reactivo_desde_descenso", fmax_kcal_a=0.5, steps=600)
                if R2["energy_kcal"] < R["energy_kcal"]:
                    log("reactivo más bajo hallado por el descenso: %.2f -> %.2f kcal/mol" % (R["energy_kcal"], R2["energy_kcal"]))
                    R = R2
                    m.write_pdb(OUT / "reactivo.pdb", R["coords"])
                    summary["etapas"]["reactivo"].update(energia_kcal=R["energy_kcal"], reoptimizado_desde_descenso=True,
                                                         **m.reaction_coordinates(R["coords"]))
            if energies[-1] < P["energy_kcal"] - 0.05:
                P2 = m.optimize(frames[-1], "producto_desde_descenso", fmax_kcal_a=0.5, steps=600)
                if P2["energy_kcal"] < P["energy_kcal"]:
                    log("producto más bajo hallado por el descenso: %.2f -> %.2f kcal/mol" % (P["energy_kcal"], P2["energy_kcal"]))
                    P = P2
                    m.write_pdb(OUT / "producto.pdb", P["coords"])
                    summary["etapas"]["producto"].update(energia_kcal=P["energy_kcal"], reoptimizado_desde_descenso=True,
                                                         **m.reaction_coordinates(P["coords"]))
            summary["etapas"]["producto"]["dE_reaccion_kcal"] = P["energy_kcal"] - R["energy_kcal"]
            if "dimero" in summary["etapas"]:
                summary["etapas"]["dimero"]["barrera_kcal"] = summary["etapas"]["dimero"]["energia_kcal"] - R["energy_kcal"]
            if "neb" in summary["etapas"]:
                summary["etapas"]["neb"]["barrera_kcal"] = float(neb["energies_kcal"][k_ts] - R["energy_kcal"])
            summary["etapas"]["escaneo"]["barrera_aprox_kcal"] = float(df["energia_kcal"].max() - R["energy_kcal"])
            # reescribir los CSV con la nueva referencia
            df["energia_rel_kcal"] = df["energia_kcal"] - R["energy_kcal"]
            df.to_csv(OUT / "escaneo.csv", index=False)
            pd.DataFrame(dict(cuadro=range(len(frames)), xi=xis, energia_kcal=energies,
                              energia_rel_kcal=np.array(energies) - R["energy_kcal"])).to_csv(OUT / "camino_descenso.csv", index=False)
            if "neb" in summary["etapas"]:
                pd.DataFrame(dict(imagen=range(len(images)), xi=neb["xi"], energia_kcal=neb["energies_kcal"],
                                  energia_rel_kcal=neb["energies_kcal"] - R["energy_kcal"])).to_csv(OUT / "neb.csv", index=False)
        except Exception as exc:
            summary["errores"]["descenso"] = traceback.format_exc()
            log("descenso falló:", exc)

    # 8. el mismo clúster fuera de la enzima (agua implícita COSMO): MOPAC nativo
    #    SADDLE (= QST2 de Leonardo) -> TS -> FORCETS -> IRC
    t0 = time.time()
    EPS = 78.4
    agua = dict(eps=EPS, metodo="PM7/COSMO, anclajes fijos, sin el campo de la enzima")
    try:
        Ra = m.native_optimize(R["coords"], "agua/reactivo", eps=EPS)
        Pa = m.native_optimize(P["coords"], "agua/producto", eps=EPS)
        agua.update(E_reactivo=Ra["energy_kcal"], E_producto=Pa["energy_kcal"], dE_reaccion=Pa["energy_kcal"] - Ra["energy_kcal"])
        m.write_pdb(OUT / "agua_reactivo.pdb", Ra["coords"])
        m.write_pdb(OUT / "agua_producto.pdb", Pa["coords"])
        try:
            Sa = m.native_saddle(Ra["coords"], Pa["coords"], "agua/saddle_qst2", eps=EPS)
            agua["E_saddle_qst2"] = Sa["energy_kcal"]
            agua["barrera_saddle_qst2_kcal"] = Sa["energy_kcal"] - Ra["energy_kcal"]
        except Exception as exc:
            agua["saddle_error"] = str(exc)[:300]
        # El refinamiento TS parte del punto de silla QM/MM (geometría más fiable que la de SADDLE)
        Ta = m.native_ts(ts_xyz, "agua/ts", eps=EPS)
        agua["E_ts"] = Ta["energy_kcal"]
        agua["barrera_kcal"] = Ta["energy_kcal"] - Ra["energy_kcal"]
        agua.update({f"ts_{k}": v for k, v in m.reaction_coordinates(Ta["coords"]).items()})
        m.write_pdb(OUT / "agua_ts.pdb", Ta["coords"])
        Fa = m.native_frequencies(Ta["coords"], "agua/forcets", eps=EPS)
        agua["validacion_ts"] = Fa["validation"]
        agua["frecuencias_mas_bajas"] = [float(v) for v in np.sort(Fa["freq_cm_signed"])[:6]]
        if Fa["normal_modes_cart"] is not None:
            np.savez_compressed(OUT / "agua_frecuencias_ts.npz", freq_cm_signed=Fa["freq_cm_signed"],
                                modos=Fa["normal_modes_cart"], ts_xyz=Ta["coords"], simbolos=np.array(m.symbols))
        Ia = m.native_irc(Ta["coords"], "agua/irc", eps=EPS)
        fr = Ia["irc_frames"]
        coords = np.asarray(fr.get("coords", [])) if isinstance(fr, dict) else np.zeros((0,))
        if coords.ndim == 3 and len(coords) > 1:
            en_key = next((k for k in ("energies_kcal", "energy_kcal", "energies", "heat_kcal") if k in fr), None)
            en = np.asarray(fr[en_key], dtype=float) if en_key else np.full(len(coords), np.nan)
            xis = [m.reaction_coordinate(c)[0] for c in coords]
            pd.DataFrame(dict(cuadro=range(len(coords)), xi=xis, energia_kcal=en,
                              energia_rel_kcal=en - Ra["energy_kcal"])).to_csv(OUT / "agua_irc.csv", index=False)
            save_frames(OUT / "agua_irc.xyz", m, coords, en)
            agua["irc_cuadros"] = int(len(coords))
        else:
            agua["irc_cuadros"] = 0
            agua["irc_claves"] = list(fr.keys()) if isinstance(fr, dict) else str(type(fr))
        log("agua (COSMO): barrera TS = %.2f kcal/mol; dE = %.2f" % (agua["barrera_kcal"], agua["dE_reaccion"]))
    except Exception as exc:
        summary["errores"]["agua"] = traceback.format_exc()
        log("modelo en agua falló:", exc)
    agua["segundos"] = time.time() - t0
    summary["etapas"]["agua"] = agua

    summary["tiempo_total_s"] = time.time() - t_start
    summary["llamadas_mopac"] = m.n_calls
    summary["modo_rapido"] = QUICK
    (OUT / "resumen.json").write_text(json.dumps(summary, indent=1, default=float, ensure_ascii=False) + "\n")
    log("listo en %.0f min, %d llamadas MOPAC" % ((time.time() - t_start) / 60, m.n_calls))


if __name__ == "__main__":
    main()
