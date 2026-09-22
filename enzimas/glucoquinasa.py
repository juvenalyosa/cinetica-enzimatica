"""Modelo QM/MM de la reacción de la glucoquinasa: glucosa + ATP -> glucosa-6-fosfato + ADP.

Energía QM/MM (embedding electrostático, entorno fijo):

    E(x_QM) = E_PM7[x_QM; V(x_QM)]  +  E_LJ(QM-MM)  (+ restricción opcional)

* ``E_PM7[x; V]`` es la energía MOPAC PM7 de la región QM en el potencial ``V`` que crean
  las cargas puntuales Amber del entorno (archivo ``mol.in``, palabra clave ``QMMM``).
  El potencial se **recalcula en cada evaluación** en las posiciones actuales.
* ``E_LJ`` es la repulsión/dispersión de Lennard-Jones entre átomos QM y MM con los
  parámetros del campo de fuerza (sin ella la región QM se hunde en la proteína,
  porque la electrostática sola es atractiva).
* Gradiente: MOPAC devuelve dE/dx tratando ``V`` como constante; se añade el término
  de Hellmann-Feynman ``q_i * grad V(x_i)`` (en NDDO la interacción es ``sum_i q_i V_i``)
  más el gradiente LJ.  ``check_gradient`` lo verifica por diferencias finitas.

La región QM: glucosa completa, fragmento C5'-trifosfato del ATP, cadenas laterales de
Asp205 (base catalítica), Lys169 y Thr228, Mg2+ y sus aguas.  Los cortes covalentes se
saturan con átomos de enlace (H a 1.09 A) que permanecen fijos junto con el átomo
QM al que están unidos (anclaje).  Numeración interna de tleap: cristal = tleap + 4.

Herramientas de camino de reacción (ASE): optimización L-BFGS, escaneo restringido de
la coordenada xi = d(PG-O3B) - d(PG-O6), NEB con imagen trepadora (los extremos R y P,
análogo moderno de QST2), refinamiento del TS con el método del dímero, frecuencias
numéricas y descenso desde el TS.  Los métodos nativos de MOPAC (SADDLE = QST2, TS,
FORCETS, IRC) se usan para el modelo **en vacío**, donde no hay potencial externo que
congelar.
"""
from __future__ import annotations

import json
import re
import subprocess
import time
from pathlib import Path

import numpy as np

from .qmmm_mopac import mopac_validate_ts_frequencies, parse_mopac_aux_text, parse_mopac_xyz_trajectory

COULOMB_KCAL_A = 332.0637  # kcal/mol * A / e^2
LINK_BOND_A = 1.09
KCAL_TO_EV = 0.0433641
EV_TO_KCAL = 1.0 / KCAL_TO_EV

QM_RESIDUES = {
    "GLC": dict(resseq=455, atoms="all"),
    "ATP": dict(resseq=456, atoms=["PG", "O1G", "O2G", "O3G", "PB", "O1B", "O2B", "O3B",
                                   "PA", "O1A", "O2A", "O3A", "O5'", "C5'"], hydrogens_of=["C5'"]),
    "ASP205": dict(resseq=201, atoms=["CB", "HB2", "HB3", "CG", "OD1", "OD2"]),
    "LYS169": dict(resseq=165, atoms=["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3",
                                      "CE", "HE2", "HE3", "NZ", "HZ1", "HZ2", "HZ3"]),
    "THR228": dict(resseq=224, atoms=["CB", "HB", "OG1", "HG1", "CG2", "HG21", "HG22", "HG23"]),
    "MG": dict(resseq=457, atoms="all"),
}
LINK_CUTS = [  # (resseq, átomo QM ancla, átomo MM sustituido por H, átomos MM cuya carga se elimina)
    (456, "C5'", "C4'", ["C4'"]),
    (201, "CB", "CA", ["CA", "HA"]),
    (165, "CB", "CA", ["CA", "HA"]),
    (224, "CB", "CA", ["CA", "HA"]),
]
QM_FORMAL_CHARGE = -2  # glucosa 0, trifosfato -4, Asp -1, Lys +1, Thr 0, Mg +2
MG_WATER_CUTOFF_A = 2.6
LJ_EXCLUSION_A = 2.8  # átomos MM 1-2/1-3 del corte: sin LJ con la región QM


def _log(*a):
    print("[qmmm]", *a, flush=True)


class ModeloGlucoquinasa:
    """Partición QM/MM, energía + gradiente QM/MM y herramientas de camino de reacción."""

    def __init__(self, pdb_full, prmtop, mopac_exe="mopac", workdir="work/qmmm", env_cutoff_a=16.0,
                 threads=2, precise=True):
        from openmm import NonbondedForce, unit
        from openmm.app import AmberPrmtopFile, PDBFile

        self.work = Path(workdir)
        self.work.mkdir(parents=True, exist_ok=True)
        self.mopac = str(mopac_exe)
        self.threads = int(threads)
        self.precise = bool(precise)
        self.method = "PM7"  # hamiltoniano semiempírico de MOPAC (PM7, PM6-D3H4, ...)
        top = AmberPrmtopFile(str(prmtop))
        system = top.createSystem()
        nb = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
        n = system.getNumParticles()
        q, sig, eps = np.zeros(n), np.zeros(n), np.zeros(n)
        for i in range(n):
            qi, si, ei = nb.getParticleParameters(i)
            q[i] = qi.value_in_unit(unit.elementary_charge)
            sig[i] = si.value_in_unit(unit.angstrom)
            eps[i] = ei.value_in_unit(unit.kilocalorie_per_mole)
        self.q_all, self.sig_all, self.eps_all = q, sig, eps
        pdb = PDBFile(str(pdb_full))
        self.xyz_all = np.array(pdb.positions.value_in_unit(unit.angstrom))
        self.atoms = [dict(index=a.index, name=a.name, resname=a.residue.name, resseq=int(a.residue.id),
                           element=a.element.symbol, z=a.element.atomic_number) for a in top.topology.atoms()]
        assert len(self.atoms) == len(self.xyz_all)
        self.n_calls = 0
        self._build_partition(env_cutoff_a)

    # ------------------------------------------------------------------ partición
    def _find(self, resseq, name):
        for a in self.atoms:
            if a["resseq"] == resseq and a["name"] == name:
                return a["index"]
        raise KeyError((resseq, name))

    def _bonded_h(self, heavy_index):
        c = self.xyz_all[heavy_index]
        res = self.atoms[heavy_index]["resseq"]
        return [a["index"] for a in self.atoms
                if a["resseq"] == res and a["element"] == "H" and np.linalg.norm(self.xyz_all[a["index"]] - c) < 1.2]

    def _build_partition(self, env_cutoff_a):
        qm = []
        for spec in QM_RESIDUES.values():
            rs = spec["resseq"]
            if spec["atoms"] == "all":
                qm += [a["index"] for a in self.atoms if a["resseq"] == rs]
            else:
                qm += [self._find(rs, n) for n in spec["atoms"]]
                for hn in spec.get("hydrogens_of", []):
                    qm += self._bonded_h(self._find(rs, hn))
        mg = self._find(QM_RESIDUES["MG"]["resseq"], "MG")
        for a in self.atoms:  # aguas coordinadas al Mg2+ (moléculas completas)
            if a["resname"] in ("WAT", "HOH") and a["element"] == "O":
                if np.linalg.norm(self.xyz_all[a["index"]] - self.xyz_all[mg]) < MG_WATER_CUTOFF_A:
                    qm += [b["index"] for b in self.atoms if b["resseq"] == a["resseq"] and b["resname"] == a["resname"]]
        qm = list(dict.fromkeys(qm))
        self.qm_global = qm
        link_xyz, link_meta, excluded, lj_excl = [], [], [], []
        for rs, qname, mname, drop in LINK_CUTS:
            iq, im = self._find(rs, qname), self._find(rs, mname)
            u = self.xyz_all[im] - self.xyz_all[iq]
            link_xyz.append(self.xyz_all[iq] + LINK_BOND_A * u / np.linalg.norm(u))
            link_meta.append(dict(index=-1, name="HL", resname=self.atoms[iq]["resname"], resseq=rs,
                                  element="H", z=1, qm_partner=iq, mm_partner=im))
            excluded += [self._find(rs, n) for n in drop]
            if rs == 456:
                excluded += self._bonded_h(im)
            # pares (QM, MM) vecinos 1-2/1-3 del corte: ambos a menos de LJ_EXCLUSION_A del ancla
            near_mm = [a["index"] for a in self.atoms
                       if a["resseq"] == rs and np.linalg.norm(self.xyz_all[a["index"]] - self.xyz_all[iq]) < LJ_EXCLUSION_A]
            near_qm = [g for g in qm if self.atoms[g]["resseq"] == rs
                       and np.linalg.norm(self.xyz_all[g] - self.xyz_all[iq]) < 1.7]  # unidos al ancla
            lj_excl.append((near_qm, near_mm))
        self.excluded = list(dict.fromkeys(excluded))
        self.qm_meta = [dict(self.atoms[i]) for i in qm] + link_meta
        self.qm_xyz = np.vstack([self.xyz_all[qm], np.array(link_xyz)])
        self.n_real, self.n_link = len(qm), len(link_xyz)
        self.n_qm = self.n_real + self.n_link
        center = self.qm_xyz[: self.n_real].mean(axis=0)
        dist = np.linalg.norm(self.xyz_all - center, axis=1)
        close_res = {(a["resname"], a["resseq"]) for a, d in zip(self.atoms, dist) if d < env_cutoff_a}
        qm_set, ex_set = set(qm), set(self.excluded)
        self.mm_global = [a["index"] for a in self.atoms
                          if (a["resname"], a["resseq"]) in close_res and a["index"] not in qm_set and a["index"] not in ex_set]
        self.mm_xyz = self.xyz_all[self.mm_global]
        self.mm_q = self.q_all[self.mm_global].copy()
        # conservación de carga total (región QM con carga formal entera)
        cluster = list(qm_set) + list(ex_set) + self.mm_global
        target_mm = float(self.q_all[cluster].sum()) - QM_FORMAL_CHARGE
        receivers = sorted({k for k, g in enumerate(self.mm_global)
                            if any(self.atoms[g]["resseq"] == rs for rs, *_ in LINK_CUTS) and self.atoms[g]["element"] != "H"})
        self.charge_correction = target_mm - float(self.mm_q.sum())
        self.mm_q[receivers] += self.charge_correction / len(receivers)
        self.receivers = receivers
        # Lennard-Jones QM-MM (Lorentz-Berthelot); sin LJ para H de enlace ni vecinos 1-2/1-3 del corte
        sig_qm = np.array([self.sig_all[i] for i in qm] + [0.0] * self.n_link)
        eps_qm = np.array([self.eps_all[i] for i in qm] + [0.0] * self.n_link)
        sig_mm, eps_mm = self.sig_all[self.mm_global], self.eps_all[self.mm_global]
        self.lj_sigma = 0.5 * (sig_qm[:, None] + sig_mm[None, :])
        self.lj_eps = np.sqrt(eps_qm[:, None] * eps_mm[None, :])
        mm_pos = {g: k for k, g in enumerate(self.mm_global)}
        n_excl = 0
        for near_qm, near_mm in lj_excl:
            for gq in near_qm:
                for gm in near_mm:
                    if gm in mm_pos:
                        self.lj_eps[qm.index(gq), mm_pos[gm]] = 0.0
                        n_excl += 1
        self.n_lj_excluded_pairs = n_excl
        self.symbols = [m["element"] for m in self.qm_meta]
        self.i = {(m["resseq"], m["name"]): k for k, m in enumerate(self.qm_meta) if m["index"] >= 0}
        self.PG, self.O6, self.O3B = self.i[(456, "PG")], self.i[(455, "O6")], self.i[(456, "O3B")]
        self.OD1, self.OD2, self.MG = self.i[(201, "OD1")], self.i[(201, "OD2")], self.i[(457, "MG")]
        self.O1G, self.O2G, self.O3G = self.i[(456, "O1G")], self.i[(456, "O2G")], self.i[(456, "O3G")]
        self.HO6 = min((k for k, m in enumerate(self.qm_meta) if m["resseq"] == 455 and m["element"] == "H"),
                       key=lambda k: np.linalg.norm(self.qm_xyz[k] - self.qm_xyz[self.O6]))
        self.fixed = set(range(self.n_real, self.n_qm))
        for m in link_meta:
            self.fixed.add(qm.index(m["qm_partner"]))
        self.free = [k for k in range(self.n_qm) if k not in self.fixed]
        _log(f"QM: {self.n_real} átomos + {self.n_link} H de enlace ({len(self.free)} libres); "
             f"MM: {len(self.mm_global)} cargas; corrección de carga {self.charge_correction:+.3f} e "
             f"sobre {len(receivers)} átomos de interfaz")

    def describe(self):
        return dict(
            n_qm=self.n_real, n_link=self.n_link, n_free=len(self.free), n_mm=len(self.mm_global),
            qm_charge=QM_FORMAL_CHARGE, mm_charge=float(self.mm_q.sum()), charge_correction_e=self.charge_correction,
            qm_atoms=[dict(k=k, name=m["name"], resname=m["resname"], resseq=m["resseq"], element=m["element"],
                           global_index=m["index"]) for k, m in enumerate(self.qm_meta)],
            fixed=sorted(self.fixed),
            links=[dict(qm=m["qm_partner"], mm=m["mm_partner"]) for m in self.qm_meta if "qm_partner" in m],
            residues_qm=list(QM_RESIDUES.keys()),
        )

    # ------------------------------------------------------------------ términos MM
    def potential(self, xyz):
        """Potencial electrostático (kcal/mol/e) del entorno sobre cada átomo QM y su gradiente."""
        r = xyz[:, None, :] - self.mm_xyz[None, :, :]
        d = np.linalg.norm(r, axis=2)
        inv = 1.0 / d
        v = COULOMB_KCAL_A * (self.mm_q[None, :] * inv).sum(axis=1)
        grad = -COULOMB_KCAL_A * ((self.mm_q[None, :] * inv**3)[:, :, None] * r).sum(axis=1)
        return v, grad

    def lennard_jones(self, xyz):
        r = xyz[:, None, :] - self.mm_xyz[None, :, :]
        d = np.linalg.norm(r, axis=2)
        sr6 = (self.lj_sigma / d) ** 6
        e_pair = 4.0 * self.lj_eps * (sr6**2 - sr6)
        de_dr = 4.0 * self.lj_eps * (-12.0 * sr6**2 + 6.0 * sr6) / d
        grad = ((de_dr / d)[:, :, None] * r).sum(axis=1)
        return float(e_pair.sum()), grad

    def reaction_coordinate(self, xyz):
        """xi = d(PG-O3B) - d(PG-O6): negativa en el reactivo, positiva en el producto."""
        xyz = np.asarray(xyz, dtype=float)
        v1, v2 = xyz[self.PG] - xyz[self.O3B], xyz[self.PG] - xyz[self.O6]
        d1, d2 = np.linalg.norm(v1), np.linalg.norm(v2)
        grad = np.zeros_like(xyz)
        grad[self.PG] = v1 / d1 - v2 / d2
        grad[self.O3B] = -v1 / d1
        grad[self.O6] = v2 / d2
        return float(d1 - d2), grad

    def reaction_coordinates(self, xyz):
        xyz = np.asarray(xyz, dtype=float)
        d = lambda a, b: float(np.linalg.norm(xyz[a] - xyz[b]))
        return dict(d_PG_O6=d(self.PG, self.O6), d_PG_O3B=d(self.PG, self.O3B), d_O6_H=d(self.O6, self.HO6),
                    d_OD1_H=d(self.OD1, self.HO6), d_Mg_O6=d(self.MG, self.O6), xi=d(self.PG, self.O3B) - d(self.PG, self.O6))

    # ------------------------------------------------------------------ MOPAC
    def _geometry_lines(self, xyz, fixed):
        return ["%-2s %14.8f %d %14.8f %d %14.8f %d" % (s, p[0], 0 if k in fixed else 1, p[1], 0 if k in fixed else 1,
                                                         p[2], 0 if k in fixed else 1)
                for k, (s, p) in enumerate(zip(self.symbols, xyz))]

    def run(self, xyz, tag, keywords, product_xyz=None, fixed=None, embedding=True, timeout_s=7200, eps=None):
        """Ejecuta MOPAC PM7 sobre la región QM.  Devuelve energía (kcal/mol), coordenadas y datos AUX."""
        fixed = self.fixed if fixed is None else set(fixed)
        work = self.work / tag
        work.mkdir(parents=True, exist_ok=True)
        xyz = np.asarray(xyz, dtype=float)
        phi = self.potential(xyz)[0] if embedding else np.zeros(len(xyz))
        mol_in = ["Potencial electrostatico Amber ff14SB/GAFF2/TIP3P (kcal/mol/e)", f"{self.n_real} {self.n_link}"]
        mol_in += ["%s %.10f %.10f %.10f %.12f" % (s, p[0], p[1], p[2], v) for s, p, v in zip(self.symbols, xyz, phi)]
        (work / "mol.in").write_text("\n".join(mol_in) + "\n")
        emb = " QMMM" if embedding else ""
        if eps is not None and not embedding:
            emb += f" EPS={float(eps):.1f}"  # disolvente implícito COSMO
        header = f"{getattr(self, 'method', 'PM7')} {keywords}{emb} AUX(9) XYZ CHARGE={QM_FORMAL_CHARGE} SINGLET GEO-OK THREADS={self.threads}"
        lines = [header, f"Glucoquinasa 3FGU: {tag}", ""] + self._geometry_lines(xyz, fixed)
        if product_xyz is not None:
            lines += [""] + self._geometry_lines(np.asarray(product_xyz, dtype=float), fixed)
        (work / "job.mop").write_text("\n".join(lines) + "\n")
        for stale in ("job.out", "job.aux", "job.arc", "job.xyz", "job.end"):
            if (work / stale).exists():
                (work / stale).unlink()
        t0 = time.perf_counter()
        cp = subprocess.run([self.mopac, "job.mop"], cwd=str(work), capture_output=True, text=True, timeout=timeout_s)
        seconds = time.perf_counter() - t0
        out = (work / "job.out").read_text(errors="replace") if (work / "job.out").exists() else ""
        hof = re.findall(r"FINAL HEAT OF FORMATION\s*=\s*([-+0-9.DEd]+)\s*KCAL", out)
        if not hof:  # los trabajos FORCE/FORCETS imprimen la energía con otro rótulo
            hof = re.findall(r"HEAT OF FORMATION\s*=\s*([-+0-9.DEd]+)\s*KCAL", out)
        aux = parse_mopac_aux_text((work / "job.aux").read_text(errors="replace")) if (work / "job.aux").exists() else {}
        if cp.returncode or (not hof and aux.get("heat_of_formation_ev") is None) or "JOB ENDED NORMALLY" not in out:
            (work / "console.log").write_text(cp.stdout + "\n" + cp.stderr)
            raise RuntimeError(f"MOPAC falló en {tag}:\n{out[-2500:]}\n{cp.stderr[-500:]}")
        energy = float(hof[-1].replace("D", "E")) if hof else float(aux["heat_of_formation_ev"]) * 23.060548
        coords = aux.get("optimized_coords_ang")
        if coords is None or np.asarray(coords).shape != xyz.shape:
            coords = xyz.copy()
        self.n_calls += 1
        return dict(tag=tag, energy_kcal=energy, coords=np.asarray(coords, dtype=float), aux=aux, seconds=seconds,
                    keywords=header, embedding=embedding, workdir=str(work))

    # ------------------------------------------------------------------ energía QM/MM completa
    def energy_gradient(self, xyz, tag="calc", restraint=None, embedding=True, eps=None):
        """E_QM/MM y gradiente (kcal/mol, kcal/mol/A).  ``restraint`` = (xi0, k) armónica sobre xi."""
        xyz = np.asarray(xyz, dtype=float)
        r = self.run(xyz, tag, "1SCF GRAD" + (" PRECISE" if self.precise else ""), embedding=embedding, eps=eps)
        g_raw = np.asarray(r["aux"]["gradients_kcal_mol_ang"], dtype=float).reshape(-1, 3)
        if g_raw.shape[0] == xyz.shape[0]:
            g = g_raw
        else:  # MOPAC solo imprime gradientes de las coordenadas optimizables (átomos libres)
            free = [k for k in range(self.n_qm) if k not in self.fixed]
            assert g_raw.shape[0] == len(free), (g_raw.shape, len(free))
            g = np.zeros_like(xyz)
            g[free] = g_raw
        e = r["energy_kcal"]
        terms = dict(e_qm_kcal=e, e_lj_kcal=0.0, e_restr_kcal=0.0)
        if embedding:
            q = np.asarray(r["aux"]["atom_charges"], dtype=float).reshape(-1)
            _, dphi = self.potential(xyz)
            g = g + q[:, None] * dphi
            e_lj, g_lj = self.lennard_jones(xyz)
            e, g = e + e_lj, g + g_lj
            terms["e_lj_kcal"] = e_lj
        if restraint is not None:
            xi0, k = restraint
            xi, dxi = self.reaction_coordinate(xyz)
            e_r = 0.5 * k * (xi - xi0) ** 2
            e, g = e + e_r, g + k * (xi - xi0) * dxi
            terms["e_restr_kcal"] = e_r
            terms["xi"] = xi
        g[sorted(self.fixed)] = 0.0
        terms.update(energy_kcal=e, seconds=r["seconds"])
        return e, g, terms

    def check_gradient(self, xyz, atom_indices, h=0.005, restraint=None):
        """Comprueba el gradiente analítico por diferencias finitas centradas."""
        xyz = np.asarray(xyz, dtype=float)
        _, g, _ = self.energy_gradient(xyz, "fd/ref", restraint=restraint)
        rows = []
        for k in atom_indices:
            for c in range(3):
                xp, xm = xyz.copy(), xyz.copy()
                xp[k, c] += h
                xm[k, c] -= h
                ep = self.energy_gradient(xp, "fd/p", restraint=restraint)[0]
                em = self.energy_gradient(xm, "fd/m", restraint=restraint)[0]
                num = (ep - em) / (2 * h)
                rows.append(dict(atom=int(k), name=self.qm_meta[k]["name"], comp="xyz"[c], analytic=float(g[k, c]),
                                 numeric=float(num), diff=float(g[k, c] - num)))
        return rows

    # ------------------------------------------------------------------ ASE
    def ase_atoms(self, xyz, restraint=None, embedding=True, tag="calc", eps=None):
        from ase import Atoms
        from ase.constraints import FixAtoms

        atoms = Atoms(symbols=self.symbols, positions=np.asarray(xyz, dtype=float))
        atoms.set_constraint(FixAtoms(indices=sorted(self.fixed)))
        atoms.calc = CalculadorQMMM(self, restraint=restraint, embedding=embedding, tag=tag, eps=eps)
        return atoms

    def optimize(self, xyz, tag, fmax_kcal_a=0.5, steps=600, restraint=None, embedding=True, logfile=None, eps=None):
        """Optimización L-BFGS (ASE) con los anclajes fijos.  fmax en kcal/mol/A."""
        from ase.optimize import LBFGS

        atoms = self.ase_atoms(xyz, restraint=restraint, embedding=embedding, tag=f"{tag}/calc", eps=eps)
        work = self.work / tag
        work.mkdir(parents=True, exist_ok=True)
        opt = LBFGS(atoms, logfile=str(work / "opt.log") if logfile is None else logfile, maxstep=0.15,
                    trajectory=str(work / "opt.traj"))
        t0 = time.perf_counter()
        converged = opt.run(fmax=fmax_kcal_a * KCAL_TO_EV, steps=steps)
        e, g, terms = self.energy_gradient(atoms.get_positions(), f"{tag}/final", restraint=restraint, embedding=embedding, eps=eps)
        rc = self.reaction_coordinates(atoms.get_positions())
        _log(f"{tag}: E = {terms['e_qm_kcal'] + terms['e_lj_kcal']:.2f} kcal/mol, xi = {rc['xi']:+.2f} A, "
             f"{opt.get_number_of_steps()} pasos, {'convergido' if converged else 'NO convergido'} "
             f"({time.perf_counter() - t0:.0f} s)")
        return dict(tag=tag, coords=atoms.get_positions().copy(), energy_kcal=terms["e_qm_kcal"] + terms["e_lj_kcal"],
                    terms=terms, converged=bool(converged), steps=opt.get_number_of_steps(), **rc)

    def scan(self, xyz_start, xi_values, tag="scan", k_kcal_a2=200.0, fmax_kcal_a=0.5, steps=400):
        """Escaneo relajado sobre xi con una restricción armónica; devuelve puntos con E sin la restricción."""
        points, xyz = [], np.asarray(xyz_start, dtype=float).copy()
        for n, xi0 in enumerate(xi_values):
            r = self.optimize(xyz, f"{tag}/punto_{n:02d}", fmax_kcal_a=fmax_kcal_a, steps=steps, restraint=(float(xi0), k_kcal_a2))
            xyz = r["coords"]
            points.append(dict(punto=n, xi_objetivo=float(xi0), energia_kcal=r["energy_kcal"], coords=xyz.copy(),
                               converged=r["converged"], **self.reaction_coordinates(xyz)))
        return points

    def neb(self, images_xyz, tag="neb", climb=True, fmax_kcal_a=1.0, steps=300, k_spring=0.1):
        """NEB con imagen trepadora entre reactivo y producto (imágenes intermedias del escaneo)."""
        from ase.mep import NEB
        from ase.optimize import FIRE

        images = [self.ase_atoms(x, tag=f"{tag}/img_{n:02d}") for n, x in enumerate(images_xyz)]
        neb = NEB(images, climb=climb, k=k_spring, method="improvedtangent", allow_shared_calculator=False)
        work = self.work / tag
        work.mkdir(parents=True, exist_ok=True)
        opt = FIRE(neb, logfile=str(work / "neb.log"), trajectory=str(work / "neb.traj"), maxstep=0.1)
        t0 = time.perf_counter()
        converged = opt.run(fmax=fmax_kcal_a * KCAL_TO_EV, steps=steps)
        energies = np.array([img.get_potential_energy() for img in images]) * EV_TO_KCAL
        coords = [img.get_positions().copy() for img in images]
        _log(f"{tag}: {len(images)} imágenes, E_max - E_0 = {energies.max() - energies[0]:.2f} kcal/mol, "
             f"{'convergido' if converged else 'NO convergido'} en {opt.get_number_of_steps()} pasos ({time.perf_counter() - t0:.0f} s)")
        return dict(energies_kcal=energies, coords=coords, converged=bool(converged), steps=opt.get_number_of_steps(),
                    xi=[self.reaction_coordinate(c)[0] for c in coords])

    def dimer(self, xyz_guess, tag="dimer", fmax_kcal_a=0.5, steps=400, mode_guess=None, embedding=True, eps=None):
        """Refinamiento del punto de silla con el método del dímero (solo gradientes)."""
        from ase.mep import DimerControl, MinModeAtoms, MinModeTranslate

        atoms = self.ase_atoms(xyz_guess, tag=f"{tag}/calc", embedding=embedding, eps=eps)
        work = self.work / tag
        work.mkdir(parents=True, exist_ok=True)
        mask = [k not in self.fixed for k in range(self.n_qm)]
        control = DimerControl(initial_eigenmode_method="displacement",
                               displacement_method="vector" if mode_guess is not None else "gauss",
                               logfile=str(work / "dimer.log"), mask=mask, dimer_separation=0.01,
                               max_num_rot=3, maximum_translation=0.1, trial_trans_step=0.01)
        d_atoms = MinModeAtoms(atoms, control)
        if mode_guess is not None:
            d_atoms.displace(displacement_vector=np.asarray(mode_guess, dtype=float) * 0.01)
        else:
            vec = np.zeros((self.n_qm, 3))
            for k in (self.PG, self.O6, self.O3B, self.O1G, self.O2G, self.O3G):
                vec[k] = np.random.default_rng(0).normal(size=3)
            d_atoms.displace(displacement_vector=vec / np.linalg.norm(vec) * 0.01)
        opt = MinModeTranslate(d_atoms, logfile=str(work / "translate.log"), trajectory=str(work / "dimer.traj"))
        t0 = time.perf_counter()
        converged = opt.run(fmax=fmax_kcal_a * KCAL_TO_EV, steps=steps)
        e, g, terms = self.energy_gradient(d_atoms.get_positions(), f"{tag}/final", embedding=embedding, eps=eps)
        rc = self.reaction_coordinates(d_atoms.get_positions())
        _log(f"{tag}: E = {e:.2f} kcal/mol, xi = {rc['xi']:+.2f}, curvatura {d_atoms.get_curvature():.4f} eV/A^2, "
             f"{'convergido' if converged else 'NO convergido'} ({time.perf_counter() - t0:.0f} s)")
        return dict(coords=d_atoms.get_positions().copy(), energy_kcal=e, converged=bool(converged),
                    curvature=float(d_atoms.get_curvature()), eigenmode=np.asarray(d_atoms.get_eigenmode()), **rc)

    def vibrations(self, xyz, tag="vib", indices=None, delta=0.01, nfree=2, embedding=True, eps=None):
        """Frecuencias numéricas (ASE) sobre los átomos libres.  Devuelve frecuencias con signo (cm-1) y modos."""
        from ase.vibrations import Vibrations

        atoms = self.ase_atoms(xyz, tag=f"{tag}/calc", embedding=embedding, eps=eps)
        indices = self.free if indices is None else list(indices)
        work = self.work / tag
        work.mkdir(parents=True, exist_ok=True)
        vib = Vibrations(atoms, indices=indices, name=str(work / "vib"), delta=delta, nfree=nfree)
        vib.clean()
        vib.run()
        vd = vib.get_vibrations()
        freqs = vd.get_frequencies()  # cm-1, complejas para modos imaginarios
        signed = np.array([-abs(f.imag) if abs(f.imag) > 1e-6 else f.real for f in freqs])
        modes = vd.get_modes(all_atoms=True)  # (n_modos, n_atoms, 3)
        validation = mopac_validate_ts_frequencies(freq_cm_signed=signed, threshold_cm=10.0)
        order = np.argsort(signed)
        _log(f"{tag}: frecuencias más bajas {np.round(signed[order][:4], 1)} cm-1; modos imaginarios: {validation['imaginary_mode_count']}")
        return dict(freq_cm_signed=signed, modes=np.asarray(modes), validation=validation, indices=indices)

    @staticmethod
    def harmonic_thermo(freq_cm_signed, temperature_k=298.15, min_cm=50.0):
        """Termoquímica armónica (ASE) a partir de frecuencias en cm-1.

        Se ignoran los modos imaginarios (el modo de reacción del TS) y los modos por debajo de
        ``min_cm`` se elevan a ``min_cm`` (corrección cuasi-armónica: los modos muy blandos
        dominan la entropía y son poco fiables).  Devuelve ZPE, H_vib, S_vib y G_vib en kcal/mol
        (S en cal/mol/K) para el subconjunto de átomos libres.
        """
        from ase.thermochemistry import HarmonicThermo
        from ase.units import invcm

        freqs = np.asarray(freq_cm_signed, dtype=float)
        real = np.clip(freqs[freqs > 0.0], min_cm, None)
        thermo = HarmonicThermo(vib_energies=real * invcm)  # eV
        zpe = thermo.get_ZPE_correction() * EV_TO_KCAL
        h = thermo.get_internal_energy(temperature_k, verbose=False) * EV_TO_KCAL
        s = thermo.get_entropy(temperature_k, verbose=False) * EV_TO_KCAL * 1000.0
        g = thermo.get_helmholtz_energy(temperature_k, verbose=False) * EV_TO_KCAL
        return dict(zpe_kcal=float(zpe), h_vib_kcal=float(h), s_vib_cal=float(s), g_vib_kcal=float(g),
                    n_modes=int(real.size), n_imag=int((freqs < 0).sum()), temperature_k=float(temperature_k))

    def descend(self, xyz_ts, mode, tag, direction=+1, displacement_a=0.15, fmax_kcal_a=0.5, steps=400, every=5,
                embedding=True, eps=None):
        """Desciende desde el TS a lo largo del modo imaginario (camino de mínima energía descendente)."""
        from ase.optimize import FIRE

        mode = np.asarray(mode, dtype=float)
        mode = mode / np.linalg.norm(mode)
        xyz0 = np.asarray(xyz_ts, dtype=float) + direction * displacement_a * mode
        atoms = self.ase_atoms(xyz0, tag=f"{tag}/calc", embedding=embedding, eps=eps)
        frames, energies = [], []

        def record():
            frames.append(atoms.get_positions().copy())
            energies.append(atoms.get_potential_energy() * EV_TO_KCAL)

        work = self.work / tag
        work.mkdir(parents=True, exist_ok=True)
        opt = FIRE(atoms, logfile=str(work / "descend.log"), maxstep=0.05, dt=0.05, dtmax=0.1)
        opt.attach(record, interval=every)
        record()
        converged = opt.run(fmax=fmax_kcal_a * KCAL_TO_EV, steps=steps)
        record()
        _log(f"{tag}: {len(frames)} cuadros, E final {energies[-1]:.2f} kcal/mol, xi = {self.reaction_coordinate(frames[-1])[0]:+.2f}, "
             f"{'convergido' if converged else 'NO convergido'}")
        return dict(frames=frames, energies_kcal=np.array(energies), converged=bool(converged),
                    xi=[self.reaction_coordinate(f)[0] for f in frames])

    # ------------------------------------------------------------------ MOPAC nativo (para el modelo en vacío)
    def native_optimize(self, xyz, tag, embedding=False, eps=None):
        r = self.run(xyz, tag, "EF GNORM=0.5 CYCLES=3000 LET", embedding=embedding, eps=eps)
        _log(f"{tag}: E = {r['energy_kcal']:.2f} kcal/mol ({r['seconds']:.0f} s)")
        return r

    def native_saddle(self, xyz_react, xyz_prod, tag, embedding=False, eps=None):
        """QST2 de Leonardo = SADDLE de MOPAC: busca el TS entre dos extremos."""
        r = self.run(xyz_react, tag, "SADDLE CYCLES=500 LET", product_xyz=xyz_prod, embedding=embedding, eps=eps)
        _log(f"{tag}: E = {r['energy_kcal']:.2f} kcal/mol ({r['seconds']:.0f} s)")
        return r

    def native_ts(self, xyz_guess, tag, embedding=False, eps=None):
        r = self.run(xyz_guess, tag, "TS GNORM=0.5 CYCLES=2000 LET RECALC=5", embedding=embedding, eps=eps)
        _log(f"{tag}: E = {r['energy_kcal']:.2f} kcal/mol ({r['seconds']:.0f} s)")
        return r

    def native_frequencies(self, xyz, tag, embedding=False, eps=None):
        r = self.run(xyz, tag, "FORCETS LET", embedding=embedding, eps=eps)
        signed = np.asarray(r["aux"].get("freq_cm_signed", []), dtype=float)
        r["freq_cm_signed"] = signed
        r["validation"] = mopac_validate_ts_frequencies(freq_cm_signed=signed, threshold_cm=10.0)
        modes = r["aux"].get("normal_modes_cart")
        r["normal_modes_cart"] = None if modes is None else np.asarray(modes, dtype=float)
        _log(f"{tag}: frecuencias más bajas {np.round(np.sort(signed)[:4], 1)} cm-1; imaginarias: {r['validation']['imaginary_mode_count']}")
        return r

    def native_irc(self, xyz_ts, tag, embedding=False, eps=None):
        r = self.run(xyz_ts, tag, "IRC=1* X-PRIORITY=0.05 LET", embedding=embedding, eps=eps)
        xyz_file = Path(r["workdir"]) / "job.xyz"
        r["irc_frames"] = parse_mopac_xyz_trajectory(xyz_file.read_text(errors="replace"), expected_atoms=self.n_qm) if xyz_file.exists() else {}
        return r

    # ------------------------------------------------------------------ química y salida
    def build_product(self, xyz_react, po_bond_a=1.62):
        """Geometría inicial del producto: H de O6 -> OD1(Asp205); PG sobre el eje O3B->O6 a 1.62 A de O6;
        los tres O gamma se reflejan a través del plano perpendicular al eje (inversión de paraguas)."""
        xyz = np.asarray(xyz_react, dtype=float).copy()
        o6, o3b, pg = xyz[self.O6], xyz[self.O3B], xyz[self.PG]
        u = (o6 - o3b) / np.linalg.norm(o6 - o3b)
        pg_new = o6 - po_bond_a * u
        for k in (self.O1G, self.O2G, self.O3G):
            v = xyz[k] - pg
            xyz[k] = pg_new + (v - 2.0 * np.dot(v, u) * u)
        xyz[self.PG] = pg_new
        w = o6 - xyz[self.OD1]
        xyz[self.HO6] = xyz[self.OD1] + 0.98 * w / np.linalg.norm(w)
        return xyz

    def contacts(self, xyz, cutoff=2.5):
        """Contactos pesados QM-MM por debajo de ``cutoff`` (diagnóstico de penetración)."""
        xyz = np.asarray(xyz, dtype=float)
        d = np.linalg.norm(xyz[:, None, :] - self.mm_xyz[None, :, :], axis=2)
        out = []
        for k in range(self.n_real):
            if self.symbols[k] == "H":
                continue
            for j in np.where(d[k] < cutoff)[0]:
                g = self.mm_global[j]
                if self.atoms[g]["element"] != "H":
                    out.append((f"{self.qm_meta[k]['resname']}:{self.qm_meta[k]['name']}",
                                f"{self.atoms[g]['resname']}{self.atoms[g]['resseq']}:{self.atoms[g]['name']}", round(float(d[k, j]), 2)))
        return out

    def write_pdb(self, path, xyz):
        lines = ["HETATM%5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s" % (
            k, m["name"][:4], m["resname"][:3], m["resseq"], p[0], p[1], p[2], m["element"])
            for k, (m, p) in enumerate(zip(self.qm_meta, xyz), 1)]
        Path(path).write_text("\n".join(lines) + "\nEND\n")

    def write_env_pdb(self, path, radius_a=8.0):
        center = self.qm_xyz[: self.n_real].mean(axis=0)
        lines = []
        for k, g in enumerate(self.mm_global, 1):
            m, p = self.atoms[g], self.xyz_all[g]
            if np.linalg.norm(p - center) > radius_a:
                continue
            lines.append("ATOM  %5d %-4s %3s B%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s" % (
                k % 100000, m["name"][:4], m["resname"][:3], m["resseq"] % 10000, p[0], p[1], p[2], self.mm_q[k - 1], m["element"]))
        Path(path).write_text("\n".join(lines) + "\nEND\n")

    def xyz_text(self, xyz, comment=""):
        return "\n".join([str(self.n_qm), comment] + ["%-2s %12.6f %12.6f %12.6f" % (s, *p) for s, p in zip(self.symbols, xyz)]) + "\n"


class CalculadorQMMM(__import__("ase.calculators.calculator", fromlist=["Calculator"]).Calculator):
    """Calculador ASE: energía en eV y fuerzas en eV/A a partir de ``ModeloGlucoquinasa``."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, modelo, restraint=None, embedding=True, tag="calc", eps=None):
        super().__init__()
        self.modelo, self.restraint, self.embedding, self.tag, self.eps = modelo, restraint, embedding, tag, eps
        self.terms = {}

    def calculate(self, atoms=None, properties=("energy",), system_changes=None):
        from ase.calculators.calculator import all_changes

        super().calculate(atoms, list(properties), all_changes if system_changes is None else system_changes)
        e, g, terms = self.modelo.energy_gradient(self.atoms.get_positions(), self.tag, restraint=self.restraint,
                                                  embedding=self.embedding, eps=self.eps)
        self.terms = terms
        self.results = dict(energy=e * KCAL_TO_EV, forces=-g * KCAL_TO_EV)
