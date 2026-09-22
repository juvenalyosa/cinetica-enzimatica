#!/usr/bin/env python
"""Dinámica molecular clásica (OpenMM) del complejo glucoquinasa·glucosa·ATP·Mg2+.

Etapas (todas con Amber ff14SB/GAFF2/TIP3P, PME, 2 fs, enlaces con H restringidos):
  1. minimización         solvente libre, soluto pesado restringido; luego todo libre
  2. calentamiento NVT    0 -> 300 K en 100 ps, soluto pesado restringido (5 kcal/mol/A^2)
  3. equilibración NPT    200 ps, restricciones 5 -> 1 -> 0.2 kcal/mol/A^2
  4. producción NPT       1 ns libre, 300 K, 1 bar

Productos en data/precalculado/md/ (tamaños pequeños, pensados para GitHub y Colab):
  minimizado_soluto.pdb        estructura minimizada (proteína + ligandos + iones + aguas cristalinas)
  trayectoria_soluto.dcd       producción, solo soluto, un cuadro cada 10 ps
  analisis.csv                 tiempo, T, energías, RMSD, distancias catalíticas
  final_soluto.pdb             último cuadro (soluto)
  final_minimizado.pdb         último cuadro minimizado, sistema COMPLETO (punto de partida del QM/MM)
  md.json                      metadatos y tiempos de cómputo
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import openmm
from openmm import unit
from openmm.app import (DCDReporter, AmberInpcrdFile, AmberPrmtopFile, HBonds, PDBFile, PME, Simulation,
                        StateDataReporter)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from enzimas.datos import ruta  # noqa: E402  (descomprime .gz si hace falta)

SYS = ROOT / "data" / "precalculado" / "sistema"
OUT = ROOT / "data" / "precalculado" / "md"
WORK = ROOT / "work" / "md"
QUICK = "--rapido" in sys.argv  # prueba corta del protocolo

PROD_NS = 0.05 if QUICK else 1.0
HEAT_PS = 5 if QUICK else 100
EQ_PS = 6 if QUICK else 200


def log(*a):
    print("[md]", *a, flush=True)


def pick_platform():
    for name in ("CUDA", "OpenCL", "CPU"):
        try:
            p = openmm.Platform.getPlatformByName(name)
            props = {"Precision": "mixed"} if name == "CUDA" else {}
            return p, props
        except Exception:
            continue
    raise RuntimeError("sin plataforma OpenMM")


def solute_indices(topology):
    return [a.index for a in topology.atoms() if a.residue.name not in ("WAT", "HOH", "Na+", "Cl-")]


def crystal_water_count(prmtop_path):
    meta = json.loads((SYS / "preparacion.json").read_text(encoding="utf-8"))
    return int(meta["sistema_vacio"]["residue_counts"].get("WAT", 0))


def add_restraint(system, topology, positions, k_kcal, atom_filter):
    force = openmm.CustomExternalForce("0.5*k_res*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k_res", k_kcal * unit.kilocalories_per_mole / unit.angstrom**2)
    for n in ("x0", "y0", "z0"):
        force.addPerParticleParameter(n)
    for atom in topology.atoms():
        if atom_filter(atom):
            force.addParticle(atom.index, positions[atom.index])
    system.addForce(force)
    return force


def solute_heavy(atom):
    return atom.residue.name not in ("WAT", "Na+", "Cl-") and atom.element.symbol != "H"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    prmtop = AmberPrmtopFile(str(ruta("sistema/complejo.prmtop")))
    inpcrd = AmberInpcrdFile(str(ruta("sistema/complejo.inpcrd")))
    topology = prmtop.topology
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0 * unit.nanometer,
                                 constraints=HBonds, rigidWater=True)
    restraint = add_restraint(system, topology, inpcrd.positions, 5.0, solute_heavy)
    system.addForce(openmm.MonteCarloBarostat(1.0 * unit.bar, 300 * unit.kelvin, 25))
    barostat_index = system.getNumForces() - 1
    platform, props = pick_platform()
    log("plataforma", platform.getName())
    integrator = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    sim = Simulation(topology, system, integrator, platform, props)
    sim.context.setPositions(inpcrd.positions)
    sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    sim.context.setParameter("MonteCarloPressure", 0.0)  # sin barostato hasta NPT
    n_cryst_wat = crystal_water_count(SYS / "complejo.prmtop")
    sol_idx = solute_indices(topology)
    # las aguas cristalinas van justo después del soluto en el orden de tleap
    first_wat = max(sol_idx) + 1
    cryst_idx = list(range(first_wat, first_wat + 3 * n_cryst_wat))
    keep_idx = sol_idx + cryst_idx
    steps_per_ps = 500

    # 1. minimización
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setParameter("k_res", 0.0)
    sim.minimizeEnergy(maxIterations=2000)
    st = sim.context.getState(getEnergy=True, getPositions=True)
    e1 = st.getPotentialEnergy()
    log(f"minimización: {e0.value_in_unit(unit.kilocalorie_per_mole):.4g} -> {e1.value_in_unit(unit.kilocalorie_per_mole):.4g} kcal/mol")
    write_subset(topology, st.getPositions(), keep_idx, OUT / "minimizado_soluto.pdb")
    with (WORK / "minimizado_completo.pdb").open("w") as fh:
        PDBFile.writeFile(topology, st.getPositions(), fh)
    t_min = time.time() - t0

    # 2. calentamiento NVT con restricciones
    sim.context.setParameter("k_res", 5.0)
    sim.context.setVelocitiesToTemperature(5 * unit.kelvin)
    sim.reporters.append(StateDataReporter(str(WORK / "calentamiento.log"), 500, step=True, time=True,
                                           temperature=True, potentialEnergy=True, speed=True))
    n_heat = HEAT_PS * steps_per_ps
    for i in range(20):
        integrator.setTemperature((5 + (300 - 5) * (i + 1) / 20) * unit.kelvin)
        sim.step(n_heat // 20)
    log("calentamiento listo, T =", sim.context.getState(getEnergy=True).getKineticEnergy())
    sim.reporters.clear()

    # 3. equilibración NPT, restricciones decrecientes
    sim.context.setParameter("MonteCarloPressure", 1.0)
    integrator.setTemperature(300 * unit.kelvin)
    sim.reporters.append(StateDataReporter(str(WORK / "equilibracion.log"), 500, step=True, time=True,
                                           temperature=True, potentialEnergy=True, volume=True, speed=True))
    for k in (5.0, 1.0, 0.2):
        sim.context.setParameter("k_res", k)
        sim.step(EQ_PS * steps_per_ps // 3)
    sim.context.setParameter("k_res", 0.0)
    sim.reporters.clear()
    t_eq = time.time() - t0 - t_min

    # 4. producción
    frame_ps = 1 if QUICK else 10
    sim.reporters.append(DCDReporter(str(WORK / "produccion_completa.dcd"), frame_ps * steps_per_ps,
                                     enforcePeriodicBox=False))
    sim.reporters.append(StateDataReporter(str(WORK / "produccion.log"), frame_ps * steps_per_ps, step=True,
                                           time=True, temperature=True, potentialEnergy=True, kineticEnergy=True,
                                           volume=True, density=True, speed=True))
    sim.step(int(PROD_NS * 1000 * steps_per_ps))
    sim.reporters.clear()
    t_prod = time.time() - t0 - t_min - t_eq
    st = sim.context.getState(getPositions=True, enforcePeriodicBox=False)
    write_subset(topology, st.getPositions(), keep_idx, OUT / "final_soluto.pdb")

    # 5. último cuadro minimizado: punto de partida del QM/MM (sistema completo)
    sim.minimizeEnergy(maxIterations=2000)
    st = sim.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=False)
    with (OUT / "final_minimizado.pdb").open("w") as fh:
        PDBFile.writeFile(topology, st.getPositions(), fh)
    # 6. trayectoria del soluto + análisis
    analyze(topology, keep_idx, st.getPositions())
    meta = dict(
        plataforma=platform.getName(), produccion_ns=PROD_NS, calentamiento_ps=HEAT_PS, equilibracion_ps=EQ_PS,
        paso_fs=2, cuadro_ps=frame_ps, n_atomos=system.getNumParticles(), n_soluto=len(sol_idx),
        aguas_cristalinas=n_cryst_wat, tiempo_min_s=t_min, tiempo_eq_s=t_eq, tiempo_prod_s=t_prod,
        e_final_minimizada_kcal=st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole),
        protocolo=__doc__,
    )
    (OUT / "md.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
    log("listo en %.0f s" % (time.time() - t0))


def write_subset(topology, positions, indices, path):
    import mdtraj as md
    top = md.Topology.from_openmm(topology)
    xyz = np.array(positions.value_in_unit(unit.nanometer))[None, indices, :]
    traj = md.Trajectory(xyz, top.subset(indices))
    traj.save_pdb(str(path))


def analyze(topology, keep_idx, final_positions):
    import mdtraj as md
    import pandas as pd

    top = md.Topology.from_openmm(topology)
    traj = md.load(str(WORK / "produccion_completa.dcd"), top=top, atom_indices=keep_idx)
    ref = md.load(str(OUT / "minimizado_soluto.pdb"))
    traj.superpose(ref, atom_indices=ref.topology.select("backbone"))
    traj.save_dcd(str(OUT / "trayectoria_soluto.dcd"))
    ca = ref.topology.select("name CA")
    rmsd_ca = md.rmsd(traj, ref, atom_indices=ca) * 10.0
    sel = lambda q: int(ref.topology.select(q)[0])
    pg, o6 = sel("resname ATP and name PG"), sel("resname GLC and name O6")
    o3b = sel("resname ATP and name O3B")
    asp205 = sel("resname ASP and resSeq 201 and name OD1")
    asp205b = sel("resname ASP and resSeq 201 and name OD2")
    mg = sel("resname MG")
    lys169 = sel("resname LYS and resSeq 165 and name NZ")
    pairs = np.array([[pg, o6], [pg, o3b], [o6, asp205], [o6, asp205b], [mg, pg], [lys169, o6]])
    d = md.compute_distances(traj, pairs) * 10.0
    logdf = pd.read_csv(WORK / "produccion.log")
    n = min(len(logdf), traj.n_frames)
    df = pd.DataFrame({
        "tiempo_ps": logdf.iloc[:n, 1].values,
        "temperatura_K": logdf.iloc[:n]["Temperature (K)"].values,
        "energia_potencial_kJ": logdf.iloc[:n]["Potential Energy (kJ/mole)"].values,
        "densidad_g_mL": logdf.iloc[:n]["Density (g/mL)"].values,
        "rmsd_CA_A": rmsd_ca[:n],
        "d_PG_O6_A": d[:n, 0], "d_PG_O3B_A": d[:n, 1], "d_O6_OD1asp205_A": d[:n, 2],
        "d_O6_OD2asp205_A": d[:n, 3], "d_Mg_PG_A": d[:n, 4], "d_NZlys169_O6_A": d[:n, 5],
    })
    df.to_csv(OUT / "analisis.csv", index=False)
    log("análisis: RMSD CA medio %.2f A; d(PG-O6) media %.2f A" % (df.rmsd_CA_A.mean(), df.d_PG_O6_A.mean()))


if __name__ == "__main__":
    main()
