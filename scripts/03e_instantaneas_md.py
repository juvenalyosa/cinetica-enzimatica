#!/usr/bin/env python
"""Barrera QM/MM promediada sobre conformaciones de la dinámica molecular.

La enzima no es rígida: cada instantánea de la MD da un entorno (y una barrera) algo distinto.
Para varias instantáneas (250, 500, 750 y 1000 ps) se hace:
  minimización clásica del sistema completo -> partición QM/MM -> reactivo (L-BFGS) ->
  escaneo corto de xi -> refinamiento del TS (dímero) -> descenso al reactivo (referencia más baja)
y se informa la barrera de cada una junto con la media y la desviación estándar, incluida la
instantánea original (cristal minimizado).  Resultados en data/precalculado/qmmm/instantaneas.json
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from enzimas.datos import ruta  # noqa: E402
from enzimas.glucoquinasa import ModeloGlucoquinasa  # noqa: E402

OUT = ROOT / "data" / "precalculado" / "qmmm"
WORK = ROOT / "work" / "qmmm_instantaneas"
MOPAC = "/opt/anaconda3/envs/protliginteract-full/bin/mopac"
FRAMES_PS = tuple(int(a) for a in sys.argv[1:]) or (250, 500, 750, 1000)


def log(*a):
    print("[03e]", *a, flush=True)


def minimized_snapshot(frame_index, tag):
    """Extrae un cuadro de la trayectoria completa y lo minimiza con OpenMM (sistema completo)."""
    import mdtraj as md
    import openmm
    from openmm import unit
    from openmm.app import AmberPrmtopFile, HBonds, PDBFile, PME, Simulation

    prmtop = AmberPrmtopFile(str(ruta("sistema/complejo.prmtop")))
    traj = md.load_frame(str(ROOT / "work" / "md" / "produccion_completa.dcd"), frame_index, top=md.Topology.from_openmm(prmtop.topology))
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0 * unit.nanometer, constraints=HBonds, rigidWater=True)
    platform = openmm.Platform.getPlatformByName("OpenCL") if any(openmm.Platform.getPlatform(i).getName() == "OpenCL" for i in range(openmm.Platform.getNumPlatforms())) else openmm.Platform.getPlatformByName("CPU")
    sim = Simulation(prmtop.topology, system, openmm.VerletIntegrator(0.001), platform)
    sim.context.setPositions(traj.xyz[0] * unit.nanometer)
    sim.context.setPeriodicBoxVectors(*(traj.unitcell_vectors[0] * unit.nanometer))
    sim.minimizeEnergy(maxIterations=3000)
    st = sim.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=False)
    path = WORK / f"{tag}_minimizado.pdb"
    with path.open("w") as fh:
        PDBFile.writeFile(prmtop.topology, st.getPositions(), fh)
    return path, st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def main():
    t0 = time.time()
    WORK.mkdir(parents=True, exist_ok=True)
    base = json.loads((OUT / "resumen.json").read_text(encoding="utf-8"))
    previos = []
    if (OUT / "instantaneas.json").exists():
        previos = [r for r in json.loads((OUT / "instantaneas.json").read_text(encoding="utf-8"))["instantaneas"]
                   if r.get("t_ps", 0) not in FRAMES_PS and r.get("t_ps", 0) != 0 and "error" not in r]
    results = [dict(instantanea="cristal minimizado", t_ps=0, barrera_kcal=base["etapas"]["dimero"]["barrera_kcal"],
                    dE_reaccion_kcal=base["etapas"]["producto"]["dE_reaccion_kcal"],
                    ts_d_PG_O6=base["etapas"]["dimero"]["d_PG_O6"], ts_d_PG_O3B=base["etapas"]["dimero"]["d_PG_O3B"],
                    frecuencia_imaginaria_cm=base["etapas"]["frecuencias"]["frecuencias_mas_bajas"][0])]
    for t_ps in FRAMES_PS:
        tag = f"t{t_ps:04d}"
        t1 = time.time()
        try:
            pdb, e_mm = minimized_snapshot(t_ps // 10 - 1, tag)
            m = ModeloGlucoquinasa(pdb, ruta("sistema/complejo.prmtop"), mopac_exe=MOPAC, workdir=WORK / tag, threads=2)
            R = m.optimize(m.qm_xyz, "reactivo", fmax_kcal_a=0.5, steps=600)
            xi_r = R["xi"]
            scan = m.scan(R["coords"], np.linspace(xi_r, 2.0, 12), tag="escaneo", k_kcal_a2=300.0, fmax_kcal_a=1.0, steps=300)
            energies = np.array([p["energia_kcal"] for p in scan])
            k_max = int(np.argmax(energies))
            d_asp = m.reaction_coordinates(R["coords"])["d_OD1_H"]
            if k_max == len(scan) - 1:  # sin máximo interior: el producto no es estable en esta conformación
                res = dict(instantanea=f"MD {t_ps} ps", t_ps=t_ps, barrera_kcal=float("nan"), sin_producto_estable=True,
                           escaneo_max_rel_kcal=float(energies[-1] - R["energy_kcal"]), xi_max=float(scan[-1]["xi"]),
                           d_OD1_H_reactivo=float(d_asp), n_qm=m.n_qm, e_mm_min_kcal=e_mm, segundos=time.time() - t1,
                           nota="La energía sube monótonamente hasta xi = 2.0: en esta conformación no hay un mínimo de producto "
                                "(el protón de O6 no encuentra a Asp205 cerca).")
                log("%s: sin máximo interior; E(xi=%.1f) = %.1f kcal/mol sobre R; d(OD1···H) = %.2f A" % (tag, scan[-1]["xi"], res["escaneo_max_rel_kcal"], d_asp))
                results.append(res)
                (OUT / "instantaneas.json").write_text(json.dumps(dict(instantaneas=results), indent=1, default=float, ensure_ascii=False) + "\n", encoding="utf-8")
                continue
            tangent = scan[min(k_max + 1, len(scan) - 1)]["coords"] - scan[max(k_max - 1, 0)]["coords"]
            tangent[sorted(m.fixed)] = 0.0
            TS = m.dimer(scan[k_max]["coords"], tag="dimero", fmax_kcal_a=0.5, steps=200, mode_guess=tangent)
            if (not TS["converged"]) or TS["energy_kcal"] > energies[k_max] + 5.0 or TS["curvature"] > 0:
                log(f"{tag}: el dímero no dio un punto de silla válido; se usa el máximo del escaneo")
                TS = dict(coords=scan[k_max]["coords"], energy_kcal=float(energies[k_max]), converged=False, curvature=float("nan"),
                          eigenmode=tangent, **m.reaction_coordinates(scan[k_max]["coords"]), origen="maximo_escaneo")
            # referencia: reactivo más bajo alcanzable desde el TS
            mode = TS["eigenmode"]
            desc = m.descend(TS["coords"], mode, "descenso", direction=-1, steps=400)
            fin = m.optimize(desc["frames"][-1], "descenso/minimo", fmax_kcal_a=0.5, steps=600)
            cand = [(R["energy_kcal"], "opt directa"), (fin["energy_kcal"], "descenso")]
            if fin["xi"] > 0.0:  # el descenso fue hacia el producto: probar el otro sentido
                desc2 = m.descend(TS["coords"], mode, "descenso_b", direction=+1, steps=400)
                fin2 = m.optimize(desc2["frames"][-1], "descenso_b/minimo", fmax_kcal_a=0.5, steps=600)
                cand.append((fin2["energy_kcal"], "descenso_b"))
                cand = [c for c in cand if not (c[1] == "descenso")]
            e_r, origen = min(cand)
            P_side = None
            res = dict(instantanea=f"MD {t_ps} ps", t_ps=t_ps, barrera_kcal=TS["energy_kcal"] - e_r, referencia_reactivo=origen,
                       ts_origen=TS.get("origen", "dimero"), d_OD1_H_reactivo=float(d_asp),
                       E_reactivo=e_r, E_ts=TS["energy_kcal"], ts_d_PG_O6=TS["d_PG_O6"], ts_d_PG_O3B=TS["d_PG_O3B"], ts_xi=TS["xi"],
                       dimero_convergido=TS["converged"], curvatura=TS["curvature"], n_qm=m.n_qm, e_mm_min_kcal=e_mm,
                       barrera_escaneo_kcal=scan[k_max]["energia_kcal"] - e_r, segundos=time.time() - t1)
            m.write_pdb(OUT / f"instantanea_{tag}_ts.pdb", TS["coords"])
            log("%s: barrera %.2f kcal/mol (TS xi %+.2f, d(P-O6) %.2f, d(P-O3B) %.2f) en %.0f s" % (
                tag, res["barrera_kcal"], TS["xi"], TS["d_PG_O6"], TS["d_PG_O3B"], res["segundos"]))
        except Exception:
            res = dict(instantanea=f"MD {t_ps} ps", t_ps=t_ps, error=traceback.format_exc()[-600:])
            log(tag, "falló")
        results.append(res)
        (OUT / "instantaneas.json").write_text(json.dumps(dict(instantaneas=results), indent=1, default=float, ensure_ascii=False) + "\n", encoding="utf-8")
    results = sorted(results + previos, key=lambda r: r.get("t_ps", 0))
    vals = np.array([r["barrera_kcal"] for r in results if "barrera_kcal" in r and np.isfinite(r["barrera_kcal"])])
    resumen = dict(instantaneas=results, barrera_media_kcal=float(vals.mean()), barrera_sd_kcal=float(vals.std(ddof=1)) if vals.size > 1 else 0.0,
                   n=int(vals.size), nota="Cada instantánea: minimización MM completa + camino QM/MM (PM7, entorno fijo).")
    (OUT / "instantaneas.json").write_text(json.dumps(resumen, indent=1, default=float, ensure_ascii=False) + "\n", encoding="utf-8")
    log("barrera media %.2f ± %.2f kcal/mol (n = %d); listo en %.0f min" % (resumen["barrera_media_kcal"], resumen["barrera_sd_kcal"], vals.size, (time.time() - t0) / 60))


if __name__ == "__main__":
    main()
