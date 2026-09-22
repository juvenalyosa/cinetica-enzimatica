#!/usr/bin/env python
"""Prepara el sistema glucoquinasa humana + glucosa + ATP + Mg2+ (PDB 3FGU) para MD y QM/MM.

Se ejecuta UNA vez en la máquina local (necesita AmberTools: antechamber, parmchk2, tleap).
Los productos que viajan al repositorio quedan en data/precalculado/sistema/:
  complejo.prmtop / complejo.inpcrd  topología y coordenadas Amber (OpenMM las lee sin conda)
  complejo_solv.pdb                  el mismo sistema en PDB (para visualizar)
  glc.mol2, glc.frcmod, atp.mol2, atp.frcmod  parámetros GAFF2/AM1-BCC de los ligandos
  preparacion.json                   metadatos: cargas, conteos, mapa de numeración

Decisiones de modelado (explicadas en el notebook):
  * cadena A, altloc A, residuos 5-458; los bucles faltantes internos se modelan con PDBFixer.
  * AMP-PNP (ANP) se convierte en ATP reemplazando el N3B puente por O3B.
  * Asp205 desprotonado (base catalítica); histidinas en estado por defecto de tleap (HIE).
  * Mg2+ y K+ con parámetros 12-6 de Li-Merz para TIP3P; aguas cristalográficas conservadas.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "precalculado" / "sistema"
WORK = ROOT / "work" / "prep"
PDB_ID = "3FGU"
ENV_BIN = Path(sys.executable).parent

GLUCOSE_SMILES = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"  # beta-D-glucopiranosa
ATP_SMILES = (
    "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)[C@@H](O)[C@H]1O"
)  # ATP4-


def log(*a):
    print("[prep]", *a, flush=True)


def run(cmd, cwd):
    log("$", " ".join(map(str, cmd)))
    cp = subprocess.run([str(c) for c in cmd], cwd=str(cwd), capture_output=True, text=True)
    (Path(cwd) / (Path(cmd[0]).name + ".log")).write_text(cp.stdout + "\n" + cp.stderr)
    if cp.returncode:
        raise RuntimeError(f"{cmd[0]} falló:\n{cp.stdout[-2000:]}\n{cp.stderr[-2000:]}")
    return cp


def download_pdb():
    RAW.mkdir(parents=True, exist_ok=True)
    target = RAW / f"{PDB_ID}.pdb"
    if not target.exists():
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{PDB_ID}.pdb", target)
    return target


def split_pdb(pdb_path):
    """Separa proteína (cadena A, altloc A), ligandos, iones y aguas.

    Las líneas SEQRES se conservan al inicio de la proteína: PDBFixer las necesita
    para detectar los residuos faltantes (bucles) que hay que modelar.
    """
    prot, bgc, anp, ions, waters = [], [], [], [], []
    for line in pdb_path.read_text().splitlines():
        if line.startswith("SEQRES") and line[11] == "A":
            prot.append(line)
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[21] != "A" or line[16] not in (" ", "A"):
            continue
        line = line[:16] + " " + line[17:]  # borrar altloc
        resn = line[17:20].strip()
        if line.startswith("ATOM"):
            prot.append(line)
        elif resn == "BGC":
            bgc.append(line)
        elif resn == "ANP":
            anp.append(line)
        elif resn in ("MG", "K"):
            ions.append(line)
        elif resn == "HOH":
            waters.append(line)
    return prot, bgc, anp, ions, waters


def fix_protein(prot_lines):
    """Modela bucles internos faltantes y átomos pesados faltantes con PDBFixer (sin hidrógenos)."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    src = WORK / "protein_raw.pdb"
    src.write_text("\n".join(prot_lines) + "\nEND\n")
    fixer = PDBFixer(filename=str(src))
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues):
        chain = chains[key[0]]
        n_res = len(list(chain.residues()))
        if key[1] == 0 or key[1] == n_res:  # no reconstruir extremos
            del fixer.missingResidues[key]
    modeled = {str(k): v for k, v in fixer.missingResidues.items()}
    fixer.findMissingAtoms()
    missing_atoms = {f"{r.name}{r.id}": [a.name for a in atoms] for r, atoms in fixer.missingAtoms.items()}
    fixer.addMissingAtoms()
    # Los bucles modelados por PDBFixer traen choques estéricos: se relajan con una
    # minimización en Amber14 manteniendo fijos (restringidos) los átomos pesados cristalinos.
    fixer.addMissingHydrogens(7.0)
    positions = relax_modeled_loops(fixer, modeled)
    out = WORK / "protein_heavy.pdb"
    with out.open("w") as fh:
        PDBFile.writeFile(fixer.topology, positions, fh, keepIds=True)
    # tleap añade sus propios hidrógenos: se escriben solo átomos pesados
    lines = [l for l in out.read_text().splitlines()
             if not (l.startswith("ATOM") and l[76:78].strip() == "H")]
    out.write_text("\n".join(lines) + "\n")
    # tleap: sin CONECT, HIS -> HIE, quitar OXT duplicado? (PDBFixer ya lo maneja)
    lines = [l for l in out.read_text().splitlines() if l.startswith(("ATOM", "TER", "END"))]
    lines = [l[:17] + "HIE" + l[20:] if l[17:20] == "HIS" else l for l in lines]
    out.write_text("\n".join(lines) + "\n")
    resids = []
    for l in lines:
        if l.startswith("ATOM"):
            key = (l[17:20], int(l[22:26]))
            if not resids or resids[-1] != key:
                resids.append(key)
    return out, modeled, missing_atoms, resids


def relax_modeled_loops(fixer, modeled):
    """Minimiza la proteína con restricciones fuertes sobre los átomos pesados cristalinos.

    Solo los bucles modelados y los hidrógenos se mueven libremente.  Se usa Amber14
    en vacío (sin corte) únicamente para eliminar choques; no es una simulación.
    """
    import openmm
    from openmm import unit
    from openmm.app import ForceField, NoCutoff

    ff = ForceField("amber14-all.xml")
    system = ff.createSystem(fixer.topology, nonbondedMethod=NoCutoff, constraints=None)
    modeled_keys = set()
    chains = list(fixer.topology.chains())
    for key, names in modeled.items():
        chain_index, start = eval(key)
        residues = list(chains[chain_index].residues())
        for k in range(len(names)):
            modeled_keys.add(residues[start + k].index)
    restraint = openmm.CustomExternalForce("0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
    restraint.addGlobalParameter("k", 100.0 * unit.kilocalories_per_mole / unit.angstrom**2)
    for name in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(name)
    n_restrained = 0
    for atom in fixer.topology.atoms():
        if atom.element.symbol == "H" or atom.residue.index in modeled_keys:
            continue
        restraint.addParticle(atom.index, fixer.positions[atom.index])
        n_restrained += 1
    system.addForce(restraint)
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("CPU"))
    context.setPositions(fixer.positions)
    e0 = context.getState(getEnergy=True).getPotentialEnergy()
    openmm.LocalEnergyMinimizer.minimize(context, tolerance=1.0, maxIterations=5000)
    state = context.getState(getEnergy=True, getPositions=True)
    log(f"relajación de bucles: {n_restrained} átomos restringidos, "
        f"E {e0.value_in_unit(unit.kilocalorie_per_mole):.3g} -> "
        f"{state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole):.4g} kcal/mol")
    return state.getPositions()


def ligand_from_pdb(lines, smiles, resname):
    """Reconstruye orden de enlace e hidrógenos de un ligando cristalográfico con RDKit."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    block = "\n".join(lines) + "\nEND\n"
    mol = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False, proximityBonding=True)
    template = Chem.MolFromSmiles(smiles)
    mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    Chem.SanitizeMol(mol)
    molh = Chem.AddHs(mol, addCoords=True)
    # asegurar que la estereoquímica proviene de las coordenadas 3D cristalinas
    Chem.AssignStereochemistryFrom3D(molh)
    path = WORK / f"{resname.lower()}.mol"
    Chem.MolToMolFile(molh, str(path))
    charge = sum(a.GetFormalCharge() for a in molh.GetAtoms())
    # nombres PDB de los átomos pesados (los H añadidos se numeran H1..Hn)
    names, n_h = [], 0
    for a in molh.GetAtoms():
        info = a.GetPDBResidueInfo()
        if a.GetAtomicNum() != 1 and info is not None and info.GetName().strip():
            names.append(info.GetName().strip())
        else:
            n_h += 1
            names.append(f"H{n_h}")
    return path, charge, names


def anp_to_atp(anp_lines):
    """AMP-PNP -> ATP: el N3B puente pasa a O3B; el resto de coordenadas se conserva."""
    out = []
    for l in anp_lines:
        name = l[12:16].strip()
        if name == "N3B":
            l = l[:12] + " O3B" + l[16:76] + " O"
        l = l[:17] + "ATP" + l[20:]
        out.append(l)
    return out


def rename_mol2_atoms(mol2_path, names):
    """antechamber renombra los átomos (P1, O1...); se restauran los nombres PDB (PG, O1G...)."""
    lines = mol2_path.read_text().splitlines()
    out, in_atoms, k = [], False, 0
    for l in lines:
        if l.startswith("@<TRIPOS>"):
            in_atoms = l.strip() == "@<TRIPOS>ATOM"
            out.append(l)
            continue
        if in_atoms and l.strip():
            parts = l.split()
            parts[1] = names[k]
            k += 1
            l = "%7s %-8s %10s %10s %10s %-6s %5s %-8s %10s" % tuple(parts[:9])
        out.append(l)
    assert k == len(names), (k, len(names))
    mol2_path.write_text("\n".join(out) + "\n")


def parametrize(mol_path, resname, net_charge, names):
    mol2 = WORK / f"{resname.lower()}.mol2"
    frcmod = WORK / f"{resname.lower()}.frcmod"
    run(
        [ENV_BIN / "antechamber", "-i", mol_path.name, "-fi", "mdl", "-o", mol2.name, "-fo", "mol2",
         "-c", "bcc", "-nc", str(net_charge), "-at", "gaff2", "-rn", resname, "-pf", "y", "-dr", "no"],
        WORK,
    )
    rename_mol2_atoms(mol2, names)
    run([ENV_BIN / "parmchk2", "-i", mol2.name, "-f", "mol2", "-o", frcmod.name, "-s", "gaff2"], WORK)
    return mol2, frcmod


def write_ions_waters(ions, waters, protein_pdb):
    """Escribe iones y aguas cristalinas; descarta aguas que chocan con los bucles modelados."""
    prot_xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                         for l in protein_pdb.read_text().splitlines() if l.startswith("ATOM")])
    kept, dropped = [], []
    for l in waters:
        o = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        if np.linalg.norm(prot_xyz - o, axis=1).min() < 2.2:
            dropped.append(l[22:26].strip())
        else:
            kept.append(l)
    log(f"aguas cristalinas: {len(kept)} conservadas, {len(dropped)} descartadas por choque {dropped}")
    waters = kept
    lines = []
    for l in ions:
        resn = l[17:20].strip()
        name = {"MG": "MG", "K": "K+"}[resn]
        resn_amber = {"MG": "MG", "K": "K+"}[resn]
        lines.append(l[:12] + f"{name:>4s}" + " " + f"{resn_amber:>3s}" + l[20:])
    for l in waters:
        lines.append(l[:12] + "   O " + "WAT" + l[20:])
    # TER entre residuos para que tleap no intente enlazarlos
    out = WORK / "ions_waters.pdb"
    txt = []
    for l in lines:
        txt.append(l)
        txt.append("TER")
    out.write_text("\n".join(txt) + "\nEND\n")
    return out


def build_tleap(protein_pdb, ions_pdb, glc, atp):
    script = f"""
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams frcmod.ionsjc_tip3p
loadamberparams frcmod.ions234lm_126_tip3p
loadamberparams {glc[1].name}
loadamberparams {atp[1].name}
GLC = loadmol2 {glc[0].name}
ATP = loadmol2 {atp[0].name}
prot = loadpdb {protein_pdb.name}
ions = loadpdb {ions_pdb.name}
complejo = combine {{prot GLC ATP ions}}
savepdb complejo complejo_vacio.pdb
saveamberparm complejo complejo_vacio.prmtop complejo_vacio.inpcrd
solvateBox complejo TIP3PBOX 10.0
addIonsRand complejo Na+ 0
addIonsRand complejo Cl- 0
savepdb complejo complejo_solv.pdb
saveamberparm complejo complejo.prmtop complejo.inpcrd
quit
"""
    (WORK / "tleap.in").write_text(script)
    run([ENV_BIN / "tleap", "-f", "tleap.in"], WORK)


def summarize(prmtop, inpcrd):
    from openmm import unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile
    from openmm import NonbondedForce

    top = AmberPrmtopFile(str(prmtop))
    crd = AmberInpcrdFile(str(inpcrd))
    system = top.createSystem()
    nb = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
    q = np.array([nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles())])
    residues = list(top.topology.residues())
    per_res = {}
    for r in residues:
        per_res.setdefault(r.name, 0)
        per_res[r.name] += 1
    qres = {}
    for r in residues:
        if r.name in ("GLC", "ATP", "MG", "K+"):
            qres[r.name] = float(sum(q[a.index] for a in r.atoms()))
    box = crd.boxVectors
    return dict(
        n_atoms=system.getNumParticles(),
        n_residues=len(residues),
        total_charge=float(q.sum()),
        ligand_charges=qres,
        residue_counts=per_res,
        box_nm=[float(v[i].value_in_unit(unit.nanometer)) for i, v in enumerate(box)] if box is not None else None,
    )


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    pdb = download_pdb()
    prot, bgc, anp, ions, waters = split_pdb(pdb)
    log(f"proteína {len(prot)} átomos, BGC {len(bgc)}, ANP {len(anp)}, iones {len(ions)}, aguas {len(waters)}")
    protein_pdb, modeled, missing_atoms, resids = fix_protein(prot)
    log("bucles modelados:", modeled)
    glc_mol, glc_q, glc_names = ligand_from_pdb(bgc, GLUCOSE_SMILES, "GLC")
    atp_mol, atp_q, atp_names = ligand_from_pdb(anp_to_atp(anp), ATP_SMILES, "ATP")
    log(f"glucosa: {len(glc_names)} átomos, carga {glc_q}; ATP: {len(atp_names)} átomos, carga {atp_q}")
    glc = parametrize(glc_mol, "GLC", glc_q, glc_names)
    atp = parametrize(atp_mol, "ATP", atp_q, atp_names)
    ions_pdb = write_ions_waters(ions, waters, protein_pdb)
    build_tleap(protein_pdb, ions_pdb, glc, atp)
    info = summarize(WORK / "complejo.prmtop", WORK / "complejo.inpcrd")
    info_vac = summarize(WORK / "complejo_vacio.prmtop", WORK / "complejo_vacio.inpcrd")
    for name in ("complejo.prmtop", "complejo.inpcrd", "complejo_solv.pdb", "complejo_vacio.prmtop",
                 "complejo_vacio.inpcrd", "complejo_vacio.pdb", "glc.mol2", "glc.frcmod", "atp.mol2", "atp.frcmod"):
        shutil.copy(WORK / name, OUT / name)
    shutil.copy(pdb, OUT / f"{PDB_ID}.pdb")
    # mapa de numeración: índice de residuo en tleap (1-based) -> residuo cristalino
    numbering = {i + 1: f"{resn}{resi}" for i, (resn, resi) in enumerate(resids)}
    meta = dict(
        fuente=f"RCSB PDB {PDB_ID} (glucoquinasa humana, glucosa, AMP-PNP, Mg2+, K+; 2.15 A)",
        cadena="A", altloc="A", residuos_cristal=f"{resids[0][1]}-{resids[-1][1]}",
        bucles_modelados=modeled, atomos_pesados_modelados=missing_atoms,
        conversion_ligando="ANP (AMP-PNP) -> ATP: N3B -> O3B",
        campo_de_fuerza=dict(proteina="ff14SB", ligandos="GAFF2 + AM1-BCC", agua="TIP3P",
                             iones="Li-Merz 12-6 (frcmod.ions234lm_126_tip3p) + Joung-Cheatham"),
        histidinas="HIE (por defecto de tleap)", asp205="desprotonado (base catalítica)",
        sistema_solvatado=info, sistema_vacio=info_vac,
        numeracion_tleap_a_cristal=numbering,
        offset_numeracion="residuo_cristal = residuo_tleap + %d (sin huecos: los bucles fueron modelados)" % (resids[0][1] - 1),
    )
    (OUT / "preparacion.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
    log(json.dumps({k: info[k] for k in ("n_atoms", "total_charge", "ligand_charges")}, indent=1))
    log("listo ->", OUT)


if __name__ == "__main__":
    main()
