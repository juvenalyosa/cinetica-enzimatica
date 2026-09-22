import math
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from enzimas.openmm_qmmm import classify_residue_kind
from enzimas.openmm_qmmm import classify_openmm_component_ff_class
from enzimas.openmm_qmmm import build_mopac_keywords_for_operation
from enzimas.openmm_qmmm import build_openmm_parameterization_audit
from enzimas.openmm_qmmm import compare_openmm_system_term_summaries
from enzimas.openmm_qmmm import compose_openmm_input_pdb
from enzimas.openmm_qmmm import dedupe_openmm_forcefield_xml_paths
from enzimas.openmm_qmmm import default_openmm_xml_selected_keys
from enzimas.openmm_qmmm import export_forcefield_xml_from_amber_prmtop
from enzimas.openmm_qmmm import merge_openmm_periodic_torsion_rows
from enzimas.openmm_qmmm import normalize_pdb_conect_records
from enzimas.openmm_qmmm import plan_openmm_staged_minimization
from enzimas.openmm_qmmm import assess_ligand_topology_signature_compatibility
from enzimas.openmm_qmmm import read_openmm_forcefield_residue_names
from enzimas.openmm_qmmm import rewrite_pdb_component_identity
from enzimas.openmm_qmmm import resolve_openmm_component_template_plan
from enzimas.openmm_qmmm import sanitize_pdb_conect_for_standalone_components
from enzimas.openmm_qmmm import split_openmm_parameterization_targets
from enzimas.openmm_qmmm import resolve_qmmm_regions
from enzimas.openmm_qmmm import run_openmm_qmmm_distance_scan_probe
from enzimas.openmm_qmmm import run_openmm_qmmm_single_point_probe


def test_classify_residue_kind_basics():
    assert classify_residue_kind("ALA") == "biopolymer"
    assert classify_residue_kind("HOH") == "water"
    assert classify_residue_kind("NA") == "ion"
    assert classify_residue_kind("LIG") == "hetero"


def test_classify_residue_kind_common_md_variants_are_not_misread_as_hetero():
    for resn in ("ASH", "GLH", "HSD", "HSE", "LYN", "CYM", "CYX", "MSE", "ACE", "NME"):
        assert classify_residue_kind(resn) == "biopolymer"


def test_classify_openmm_component_ff_class_treats_ca_ion_as_standard():
    comp = {"resn": "CA", "chain": "A", "resi": "480", "category": "ion", "uid": "A|480|CA|"}
    assert classify_openmm_component_ff_class(comp) == "standard"


def test_default_openmm_xml_selected_keys_excludes_standard_ions():
    components = [
        {"resn": "CA", "chain": "A", "resi": "480", "category": "ion", "uid": "A|480|CA|"},
        {"resn": "LIG", "chain": "B", "resi": "101", "category": "cofactor", "uid": "B|101|LIG|"},
    ]
    assert default_openmm_xml_selected_keys(components) == ["B|101|LIG|"]


def test_split_openmm_parameterization_targets_skips_standard_components():
    components = [
        {"resn": "CA", "chain": "A", "resi": "480", "category": "ion", "uid": "A|480|CA|"},
        {"resn": "LIG", "chain": "B", "resi": "101", "category": "cofactor", "uid": "B|101|LIG|"},
    ]
    plan = split_openmm_parameterization_targets(components)
    assert [row["uid"] for row in plan["parameterize"]] == ["B|101|LIG|"]
    assert [row["uid"] for row in plan["skipped_standard"]] == ["A|480|CA|"]


def test_read_openmm_forcefield_residue_names(tmp_path):
    xml_path = tmp_path / "ligand_ff.xml"
    xml_path.write_text(
        "<ForceField><Residues><Residue name='BEN'><Atom name='C1' type='c'/></Residue></Residues></ForceField>",
        encoding="utf-8",
    )
    assert read_openmm_forcefield_residue_names(str(xml_path)) == ["BEN"]


def test_dedupe_openmm_forcefield_xml_paths_skips_duplicate_residue_templates(tmp_path):
    ben_a = tmp_path / "ben_a.xml"
    ben_b = tmp_path / "ben_b.xml"
    ben_a.write_text(
        "<ForceField><Residues><Residue name='BEN'><Atom name='C1' type='c'/></Residue></Residues></ForceField>",
        encoding="utf-8",
    )
    ben_b.write_text(
        "<ForceField><Residues><Residue name='BEN'><Atom name='C1' type='c'/></Residue></Residues></ForceField>",
        encoding="utf-8",
    )
    plan = dedupe_openmm_forcefield_xml_paths([str(ben_a), str(ben_b)])
    assert plan["paths"] == [str(ben_a)]
    assert len(plan["skipped"]) == 1
    assert plan["skipped"][0]["residue_templates"] == ["BEN"]


def test_dedupe_openmm_forcefield_xml_paths_keeps_distinct_residue_templates(tmp_path):
    ben = tmp_path / "ben.xml"
    lig = tmp_path / "lig.xml"
    ben.write_text(
        "<ForceField><Residues><Residue name='BEN'><Atom name='C1' type='c'/></Residue></Residues></ForceField>",
        encoding="utf-8",
    )
    lig.write_text(
        "<ForceField><Residues><Residue name='LIG'><Atom name='C1' type='c'/></Residue></Residues></ForceField>",
        encoding="utf-8",
    )
    plan = dedupe_openmm_forcefield_xml_paths([str(ben), str(lig)])
    assert plan["paths"] == [str(ben), str(lig)]
    assert plan["skipped"] == []


def test_resolve_openmm_component_template_plan_uses_forcefield_xml(tmp_path):
    ffxml = tmp_path / "ligand_openmm_ff.xml"
    ffxml.write_text(
        "<ForceField><Residues><Residue name='BEN'><Atom name='C1' type='c'/></Residue></Residues></ForceField>",
        encoding="utf-8",
    )
    components = [{"resn": "BEN", "chain": "A", "resi": "1", "category": "cofactor", "uid": "A|1|BEN|"}]
    plan = resolve_openmm_component_template_plan(
        components=components,
        xml_map={"A|1|BEN|": str(ffxml)},
        xml_kind_fn=lambda p: "forcefield" if str(p).endswith("_ff.xml") else "",
    )
    assert plan["paths"] == [str(ffxml)]
    assert plan["missing"] == []


def test_resolve_openmm_component_template_plan_converts_system_xml_to_forcefield(tmp_path):
    sysxml = tmp_path / "ligand_openmm_system.xml"
    ffxml = tmp_path / "ligand_openmm_ff.xml"
    sysxml.write_text("<System/>", encoding="utf-8")
    ffxml.write_text(
        "<ForceField><Residues><Residue name='BEN'><Atom name='C1' type='c'/></Residue></Residues></ForceField>",
        encoding="utf-8",
    )
    components = [{"resn": "BEN", "chain": "A", "resi": "1", "category": "cofactor", "uid": "A|1|BEN|"}]
    plan = resolve_openmm_component_template_plan(
        components=components,
        xml_map={"A|1|BEN|": str(sysxml)},
        xml_kind_fn=lambda p: "system" if str(p).endswith("_system.xml") else ("forcefield" if str(p).endswith("_ff.xml") else ""),
        guess_forcefield_xml_from_system_fn=lambda _p: str(ffxml),
    )
    assert plan["paths"] == [str(ffxml)]
    assert plan["missing"] == []


def test_resolve_openmm_component_template_plan_reports_missing_nonstandard_component():
    components = [{"resn": "BEN", "chain": "A", "resi": "1", "category": "cofactor", "uid": "A|1|BEN|"}]
    plan = resolve_openmm_component_template_plan(
        components=components,
        xml_map={},
        xml_kind_fn=lambda _p: "",
    )
    assert plan["paths"] == []
    assert [row["uid"] for row in plan["missing"]] == ["A|1|BEN|"]


def test_sanitize_pdb_conect_for_standalone_components_removes_external_links():
    lines = [
        "HETATM    1  C1  BEN A   1      0.000   0.000   0.000  1.00 20.00           C",
        "HETATM    2  C2  BEN A   1      1.300   0.000   0.000  1.00 20.00           C",
        "ATOM      3  CA  ALA A   2      2.600   0.000   0.000  1.00 20.00           C",
        "CONECT    1    2    3",
        "CONECT    2    1",
        "CONECT    3    1",
    ]
    out = sanitize_pdb_conect_for_standalone_components(
        lines,
        [{"chain": "A", "resi": "1", "resn": "BEN"}],
    )
    assert out["removed_links"] == 2
    assert "CONECT    1    2" in out["lines"]
    assert "CONECT    3    1" not in out["lines"]


def test_rewrite_pdb_component_identity_relabels_component_atoms():
    lines = [
        "ATOM      1  C1  LIG X 999      0.000   0.000   0.000  1.00 20.00           C",
        "ATOM      2  H1  LIG X 999      0.000   0.000   1.000  1.00 20.00           H",
        "CONECT    1    2",
    ]
    out = rewrite_pdb_component_identity(lines, resn="BEN", chain="A", resi="223", hetero=True)
    assert out[0].startswith("HETATM")
    assert out[0][17:20] == "BEN"
    assert out[0][21] == "A"
    assert out[0][22:26].strip() == "223"
    assert out[2] == "CONECT    1    2"


def test_normalize_pdb_conect_records_deduplicates_targets():
    lines = [
        "ATOM      1  C1  BEN A   1      0.000   0.000   0.000  1.00 20.00           C",
        "ATOM      2  C2  BEN A   1      1.400   0.000   0.000  1.00 20.00           C",
        "ATOM      3  H1  BEN A   1      0.000   1.000   0.000  1.00 20.00           H",
        "CONECT    1    2    2    3",
        "CONECT    2    1    1",
    ]
    out = normalize_pdb_conect_records(lines)
    assert "CONECT    1    2    3" in out
    assert "CONECT    2    1" in out


def test_compose_openmm_input_pdb_renumbers_and_keeps_component_internal_conect():
    polymer = [
        "ATOM     10  N   ALA A   1      0.000   0.000   0.000  1.00 20.00           N",
        "ATOM     11  CA  ALA A   1      1.300   0.000   0.000  1.00 20.00           C",
    ]
    component = [
        "HETATM   21  C1  BEN A 480      2.600   0.000   0.000  1.00 20.00           C",
        "HETATM   22  C2  BEN A 480      3.900   0.000   0.000  1.00 20.00           C",
        "CONECT   21   22",
    ]
    out = compose_openmm_input_pdb(polymer, [component])
    lines = list(out["lines"])
    assert out["atom_count"] == 4
    assert out["conect_count"] == 1
    assert lines[0].startswith("ATOM      1")
    assert any(line.startswith("HETATM    3") for line in lines)
    assert "CONECT    3    4" in lines


def test_merge_openmm_periodic_torsion_rows_combines_duplicate_quadruplets():
    rows = [
        {"tag": "Proper", "type1": "A1", "type2": "A2", "type3": "A3", "type4": "A4", "periodicity1": "1", "phase1": "0.0", "k1": "0.5"},
        {"tag": "Proper", "type1": "A1", "type2": "A2", "type3": "A3", "type4": "A4", "periodicity1": "2", "phase1": "3.14", "k1": "1.0"},
    ]
    merged = merge_openmm_periodic_torsion_rows(rows)
    assert len(merged) == 1
    row = merged[0]
    assert row["periodicity1"] == "1"
    assert row["periodicity2"] == "2"
    assert row["k2"] == "1.0"


def test_build_openmm_parameterization_audit_reports_term_inventory_and_charge():
    audit = build_openmm_parameterization_audit(
        atom_records=[
            {"charge_e": -0.2, "sigma_nm": 0.3, "epsilon_kjmol": 0.1},
            {"charge_e": 0.1, "sigma_nm": 0.25, "epsilon_kjmol": 0.05},
            {"charge_e": 0.1, "sigma_nm": 0.28, "epsilon_kjmol": 0.0},
        ],
        residue_bond_pairs=[(0, 1), (1, 2)],
        bond_terms=[{"i": 0, "j": 1}, {"i": 1, "j": 2}],
        angle_terms=[{"i": 0, "j": 1, "k": 2}],
        torsion_terms=[{"i": 0, "j": 1, "k": 2, "l": 0, "improper": False}],
        torsion_rows=[{"tag": "Proper", "type1": "A1", "type2": "A2", "type3": "A3", "type4": "A1", "periodicity1": "1", "phase1": "0.0", "k1": "0.5"}],
    )
    assert audit["counts"]["atoms"] == 3
    assert audit["counts"]["bond_terms"] == 2
    assert audit["counts"]["angle_terms"] == 1
    assert audit["counts"]["torsion_terms"] == 1
    assert math.isclose(float(audit["charge"]["sum_e"]), 0.0, abs_tol=1.0e-12)
    assert audit["lj"]["zero_epsilon_atoms"] == 1
    assert audit["flags"]["has_bonded_terms"] is True
    assert audit["flags"]["has_lj_terms"] is True


def test_compare_openmm_system_term_summaries_accepts_full_parameter_match():
    src = {
        "particles": 2,
        "force_counts": {"NonbondedForce": 1, "HarmonicBondForce": 1, "HarmonicAngleForce": 1, "PeriodicTorsionForce": 1},
        "masses_da": [12.011, 1.008],
        "nonbonded_particles": [
            {"index": 0, "charge_e": -0.1, "sigma_nm": 0.33, "epsilon_kjmol": 0.12},
            {"index": 1, "charge_e": 0.1, "sigma_nm": 0.21, "epsilon_kjmol": 0.02},
        ],
        "nonbonded_exceptions": [
            {"a": 0, "b": 1, "chargeprod_e2": -0.008333333333, "sigma_nm": 0.25, "epsilon_kjmol": 0.01},
        ],
        "bonds": [{"a": 0, "b": 1, "length_nm": 0.109, "k_kjmol_nm2": 284512.0}],
        "angles": [{"a": 0, "b": 1, "c": 0, "angle_rad": 1.91, "k_kjmol_rad2": 350.0}],
        "torsions": [{"a": 0, "b": 1, "c": 0, "d": 1, "periodicity": 3, "phase_rad": 0.0, "k_kjmol": 1.5}],
    }
    out = compare_openmm_system_term_summaries(src, dict(src), tol=1.0e-10)
    assert out["ok"] is True
    assert out["issues"] == []
    assert math.isclose(float(out["max_abs_diff"]["nonbonded_particles"]), 0.0, abs_tol=1.0e-12)


def test_compare_openmm_system_term_summaries_flags_charge_mismatch():
    src = {
        "particles": 2,
        "force_counts": {"NonbondedForce": 1},
        "masses_da": [12.011, 1.008],
        "nonbonded_particles": [
            {"index": 0, "charge_e": -0.1, "sigma_nm": 0.33, "epsilon_kjmol": 0.12},
            {"index": 1, "charge_e": 0.1, "sigma_nm": 0.21, "epsilon_kjmol": 0.02},
        ],
        "nonbonded_exceptions": [],
        "bonds": [],
        "angles": [],
        "torsions": [],
    }
    dst = dict(src)
    dst["nonbonded_particles"] = [
        {"index": 0, "charge_e": -0.2, "sigma_nm": 0.33, "epsilon_kjmol": 0.12},
        {"index": 1, "charge_e": 0.2, "sigma_nm": 0.21, "epsilon_kjmol": 0.02},
    ]
    out = compare_openmm_system_term_summaries(src, dst, tol=1.0e-10)
    assert out["ok"] is False
    assert any("nonbonded_particles_mismatch" in issue for issue in out["issues"])


def test_assess_ligand_topology_signature_compatibility_accepts_matching_aromatic_signature():
    out = assess_ligand_topology_signature_compatibility(
        {"rotb": 0, "aromatic_bonds": 6},
        {"rotb": 0, "aromatic_bonds": 6},
    )
    assert out["ok"] is True
    assert out["issues"] == []


def test_assess_ligand_topology_signature_compatibility_rejects_lost_aromaticity():
    out = assess_ligand_topology_signature_compatibility(
        {"rotb": 0, "aromatic_bonds": 6},
        {"rotb": 2, "aromatic_bonds": 0},
    )
    assert out["ok"] is False
    assert any("aromatic_bonds_too_low" in issue for issue in out["issues"])
    assert any("unexpected_flexibility" in issue for issue in out["issues"])


def test_plan_openmm_staged_minimization_builds_heavy_atom_restraint_schedule():
    plan = plan_openmm_staged_minimization(
        atom_rows=[
            {"index": 0, "name": "N", "resn": "ALA", "element": "N"},
            {"index": 1, "name": "CA", "resn": "ALA", "element": "C"},
            {"index": 2, "name": "CB", "resn": "ALA", "element": "C"},
            {"index": 3, "name": "H", "resn": "ALA", "element": "H"},
            {"index": 4, "name": "C1", "resn": "BEN", "element": "C"},
            {"index": 5, "name": "H1", "resn": "BEN", "element": "H"},
            {"index": 6, "name": "O", "resn": "HOH", "element": "O"},
        ],
        min_iters=200,
        hydrogens_added=18,
    )
    assert plan["restrained_atom_count"] == 4
    assert plan["counts"]["polymer_backbone"] == 2
    assert plan["counts"]["polymer_sidechain"] == 1
    assert plan["counts"]["ligand_heavy"] == 1
    assert [stage["label"] for stage in plan["stages"]] == [
        "restrained_prep",
        "restrained_release",
        "unrestrained",
    ]
    assert math.isclose(float(plan["stages"][0]["scale"]), 1.0, abs_tol=1.0e-12)
    assert math.isclose(float(plan["stages"][1]["scale"]), 0.25, abs_tol=1.0e-12)
    assert math.isclose(float(plan["stages"][2]["scale"]), 0.0, abs_tol=1.0e-12)
    assert sum(int(stage["iterations"]) for stage in plan["stages"]) == 200


def test_plan_openmm_staged_minimization_handles_small_iteration_counts():
    plan = plan_openmm_staged_minimization(
        atom_rows=[
            {"index": 0, "name": "C1", "resn": "BEN", "element": "C"},
            {"index": 1, "name": "C2", "resn": "BEN", "element": "C"},
            {"index": 2, "name": "H1", "resn": "BEN", "element": "H"},
        ],
        min_iters=7,
        hydrogens_added=0,
    )
    assert plan["restrained_atom_count"] == 2
    assert sum(int(stage["iterations"]) for stage in plan["stages"]) == 7
    assert plan["stages"][0]["label"] == "restrained_prep"
    assert plan["stages"][-1]["label"] == "unrestrained"


def test_resolve_qmmm_regions_auto_hetero_ligand():
    atom_rows = [
        {"index": 0, "resn": "ALA", "chain": "A", "resi": "1"},
        {"index": 1, "resn": "ALA", "chain": "A", "resi": "1"},
        {"index": 2, "resn": "LIG", "chain": "L", "resi": "101"},
        {"index": 3, "resn": "LIG", "chain": "L", "resi": "101"},
        {"index": 4, "resn": "HOH", "chain": "W", "resi": "1"},
    ]
    reg = resolve_qmmm_regions(atom_rows=atom_rows, n_atoms=5, config={})
    assert reg["qm_indices"] == [2, 3]
    assert set(reg["mm_indices"]) == {0, 1, 4}
    assert reg["info"]["mode"] == "auto_largest_hetero_residue"


def test_resolve_qmmm_regions_explicit_indices():
    atom_rows = [{"index": i, "resn": "ALA", "chain": "A", "resi": "1"} for i in range(6)]
    reg = resolve_qmmm_regions(
        atom_rows=atom_rows,
        n_atoms=6,
        config={"qm_indices": [0, 1], "mm_indices": [2, 3, 4, 5]},
    )
    assert reg["qm_indices"] == [0, 1]
    assert reg["mm_indices"] == [2, 3, 4, 5]
    assert reg["info"]["mode"] == "explicit_qm_indices"


def test_resolve_qmmm_regions_residue_key_selection():
    atom_rows = [
        {"index": 0, "resn": "ALA", "chain": "A", "resi": "1"},
        {"index": 1, "resn": "LIG", "chain": "L", "resi": "10"},
        {"index": 2, "resn": "LIG", "chain": "L", "resi": "10"},
        {"index": 3, "resn": "HOH", "chain": "W", "resi": "5"},
    ]
    reg = resolve_qmmm_regions(
        atom_rows=atom_rows,
        n_atoms=4,
        config={"qm_residue_keys": [{"chain": "L", "resi": "10", "resn": "LIG"}]},
    )
    assert reg["qm_indices"] == [1, 2]
    assert reg["info"]["mode"] == "residue_key_selection"


def test_resolve_qmmm_regions_atom_descriptor_selection():
    atom_rows = [
        {"index": 0, "name": "CA", "resn": "ALA", "chain": "A", "resi": "1"},
        {"index": 1, "name": "CB", "resn": "ALA", "chain": "A", "resi": "1"},
        {"index": 2, "name": "C1", "resn": "LIG", "chain": "L", "resi": "10"},
        {"index": 3, "name": "O1", "resn": "LIG", "chain": "L", "resi": "10"},
        {"index": 4, "name": "ZN", "resn": "ZN", "chain": "Z", "resi": "1"},
    ]
    reg = resolve_qmmm_regions(
        atom_rows=atom_rows,
        n_atoms=5,
        config={
            "qm_atom_descriptors": [
                {"chain": "L", "resi": "10", "resn": "LIG", "name": "C1"},
                {"chain": "L", "resi": "10", "resn": "LIG", "name": "O1"},
            ],
            "mm_atom_descriptors": [
                {"chain": "A", "resi": "1", "resn": "ALA", "name": "CA"},
                {"chain": "A", "resi": "1", "resn": "ALA", "name": "CB"},
                {"chain": "Z", "resi": "1", "resn": "ZN", "name": "ZN"},
            ],
        },
    )
    assert reg["qm_indices"] == [2, 3]
    assert reg["mm_indices"] == [0, 1, 4]
    assert reg["info"]["mode"] == "atom_descriptor_selection"
    assert int(reg["info"]["qm_atom_count"]) == 2
    assert int(reg["info"]["mm_atom_count"]) == 3


def test_resolve_qmmm_regions_auto_ligand_pocket_uses_coordinates():
    atom_rows = [
        {"index": 0, "name": "CA", "resn": "ALA", "chain": "A", "resi": "1", "element": "C", "x": 8.0, "y": 0.0, "z": 0.0},
        {"index": 1, "name": "CB", "resn": "PHE", "chain": "A", "resi": "2", "element": "C", "x": 2.4, "y": 0.0, "z": 0.0},
        {"index": 2, "name": "C1", "resn": "LIG", "chain": "L", "resi": "10", "element": "C", "x": 0.0, "y": 0.0, "z": 0.0},
        {"index": 3, "name": "O1", "resn": "LIG", "chain": "L", "resi": "10", "element": "O", "x": 1.0, "y": 0.0, "z": 0.0},
    ]
    reg = resolve_qmmm_regions(
        atom_rows=atom_rows,
        n_atoms=4,
        config={"auto_include_pocket": True, "qm_pocket_cutoff_ang": 3.0, "operation": "optimization"},
    )
    assert set(reg["qm_indices"]) == {1, 2, 3}
    assert reg["mm_indices"] == [0]
    assert reg["info"]["mode"] == "auto_ligand_pocket"


def test_resolve_qmmm_regions_forces_reaction_atoms_into_auto_region():
    atom_rows = [
        {"index": 0, "name": "CA", "resn": "ALA", "chain": "A", "resi": "1"},
        {"index": 1, "name": "NZ", "resn": "LYS", "chain": "A", "resi": "50"},
        {"index": 2, "name": "C1", "resn": "LIG", "chain": "L", "resi": "10"},
        {"index": 3, "name": "O1", "resn": "LIG", "chain": "L", "resi": "10"},
    ]
    reg = resolve_qmmm_regions(
        atom_rows=atom_rows,
        n_atoms=4,
        config={
            "force_qm_atom_descriptors": [
                {"chain": "A", "resi": "50", "resn": "LYS", "name": "NZ"}
            ]
        },
    )
    assert set(reg["qm_indices"]) == {1, 2, 3}
    assert reg["mm_indices"] == [0]
    assert int(reg["info"]["forced_qm_atom_count"]) == 1


def test_resolve_qmmm_regions_forces_reaction_residue_into_auto_region():
    atom_rows = [
        {"index": 0, "name": "CA", "resn": "ALA", "chain": "A", "resi": "1"},
        {"index": 1, "name": "NZ", "resn": "LYS", "chain": "A", "resi": "50"},
        {"index": 2, "name": "CE", "resn": "LYS", "chain": "A", "resi": "50"},
        {"index": 3, "name": "C1", "resn": "LIG", "chain": "L", "resi": "10"},
        {"index": 4, "name": "O1", "resn": "LIG", "chain": "L", "resi": "10"},
    ]
    reg = resolve_qmmm_regions(
        atom_rows=atom_rows,
        n_atoms=5,
        config={
            "force_qm_residue_keys": [
                {"chain": "A", "resi": "50", "resn": "LYS"}
            ]
        },
    )
    assert set(reg["qm_indices"]) == {1, 2, 3, 4}
    assert reg["mm_indices"] == [0]
    assert int(reg["info"]["forced_qm_atom_count"]) == 2


def test_build_mopac_keywords_for_operation():
    kw_probe = build_mopac_keywords_for_operation("probe", "")
    kw_opt = build_mopac_keywords_for_operation("opt", "")
    kw_thermo = build_mopac_keywords_for_operation("thermo", "")
    kw_electro = build_mopac_keywords_for_operation("electrostatic", "")
    kw_md_coupled = build_mopac_keywords_for_operation("md_coupled", "")
    assert "QMMM" in kw_probe and "AUX" in kw_probe
    assert "EF" in kw_opt and "QMMM" in kw_opt
    assert "1SCF" not in kw_opt
    assert "THERMO" in kw_thermo and "QMMM" in kw_thermo
    assert "1SCF" not in kw_thermo
    assert "GRAD" in kw_electro and "QMMM" in kw_electro
    assert "GRAD" in kw_md_coupled and "QMMM" in kw_md_coupled
    kw_custom = build_mopac_keywords_for_operation("probe", "PM6 1SCF")
    assert "PM6" in kw_custom and "QMMM" in kw_custom and "AUX" in kw_custom


def test_build_mopac_keywords_switches_stale_default_to_selected_operation():
    kw_opt = build_mopac_keywords_for_operation("optimization", "PM7 1SCF GRAD AUX QMMM")
    assert "EF" in kw_opt
    assert "THERMO" not in kw_opt

    kw_thermo = build_mopac_keywords_for_operation("thermochemistry", "PM7 1SCF GRAD AUX QMMM")
    assert "FORCE" in kw_thermo
    assert "THERMO" in kw_thermo


def test_run_openmm_qmmm_single_point_probe_smoke(tmp_path):
    xyz_full = np.array(
        [
            [0.0, 0.0, 0.0],  # protein
            [1.0, 0.0, 0.0],  # protein
            [2.0, 0.0, 0.0],  # ligand
            [3.0, 0.0, 0.0],  # ligand
        ],
        dtype=float,
    )
    atomic_numbers = [6, 6, 6, 8]
    charges = np.array([0.2, -0.2, 0.0, 0.0], dtype=float)
    atom_rows = [
        {"index": 0, "name": "CA", "resn": "ALA", "chain": "A", "resi": "1"},
        {"index": 1, "name": "CB", "resn": "ALA", "chain": "A", "resi": "1"},
        {"index": 2, "name": "C1", "resn": "LIG", "chain": "L", "resi": "10"},
        {"index": 3, "name": "O1", "resn": "LIG", "chain": "L", "resi": "10"},
    ]

    def _runner(cmdline, cwd=None, capture_output=None, text=None, timeout=None):
        assert str(cmdline[0]) == "mopac"
        assert cwd is not None
        aux = """
TOTAL_ENERGY:EV= -35.0000
ATOM_CHARGES[2]=
 0.25000000 -0.25000000
GRADIENTS:KCAL/MOL/ANGSTROM[6]=
  1.00000000 0.00000000 0.00000000
 -1.00000000 0.00000000 0.00000000
"""
        (tmp_path / "step.aux").write_text(aux, encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    out = run_openmm_qmmm_single_point_probe(
        xyz_full_ang=xyz_full,
        atomic_numbers_full=atomic_numbers,
        mm_charges_e=charges,
        atom_rows=atom_rows,
        config={},
        workdir=str(tmp_path),
        runner=_runner,
    )
    assert out["ok"] is True
    assert out["qm_indices"] == [2, 3]
    assert out["mm_indices"] == [0, 1]
    qmmm = dict(out["qmmm"] or {})
    assert math.isclose(float(qmmm["energy_qm_hartree"]), -35.0 / 27.211386245988, rel_tol=1e-12)
    f_corr = np.asarray(out["force_correction_kj_mol_nm"], dtype=float)
    assert f_corr.shape == (4, 3)
    diag = dict(out["diagnostics"] or {})
    assert diag.get("third_law_ok") is True
    assert float(diag.get("net_electrostatic_force_norm", 1.0)) < 1.0e-6


def test_run_openmm_qmmm_distance_scan_probe_with_mock_step():
    xyz_full = np.array(
        [
            [0.0, 0.0, 0.0],  # protein
            [2.0, 0.0, 0.0],  # ligand atom A
            [3.0, 0.0, 0.0],  # ligand atom B (scanned)
        ],
        dtype=float,
    )
    atomic_numbers = [6, 6, 8]
    charges = np.array([0.1, 0.0, 0.0], dtype=float)
    atom_rows = [
        {"index": 0, "name": "CA", "resn": "ALA", "chain": "A", "resi": "1"},
        {"index": 1, "name": "C1", "resn": "LIG", "chain": "L", "resi": "10"},
        {"index": 2, "name": "O1", "resn": "LIG", "chain": "L", "resi": "10"},
    ]

    def _mock_step(
        xyz_full_ang,
        atomic_numbers_full,
        qm_indices,
        mm_indices,
        q_mm,
        mopac_executable="mopac",
        mopac_keywords="",
        qm_total_charge=0,
        qm_multiplicity=1,
        workdir="",
        dielectric=1.0,
        runner=None,
        timeout_s=120,
    ):
        xyz = np.asarray(xyz_full_ang, dtype=float)
        qidx = list(qm_indices)
        mmidx = list(mm_indices)
        # Energy minimum at distance 1.5 A between the 2 QM atoms.
        d = float(np.linalg.norm(xyz[qidx[1]] - xyz[qidx[0]]))
        e = float((d - 1.5) ** 2 - 10.0)
        return {
            "ok": True,
            "workdir": str(workdir or ""),
            "qm_indices": list(qidx),
            "mm_indices": list(mmidx),
            "phi_qm_au": np.zeros((len(qidx),), dtype=float),
            "q_qm": np.zeros((len(qidx),), dtype=float),
            "energy_qm_hartree": float(e),
            "forces_qm_internal_kj_mol_nm": np.zeros((len(qidx), 3), dtype=float),
            "energy_qmmm_elec_kj_mol": 0.0,
            "forces_qmmm_elec_qm_kj_mol_nm": np.zeros((len(qidx), 3), dtype=float),
            "forces_qmmm_elec_mm_kj_mol_nm": np.zeros((len(mmidx), 3), dtype=float),
        }

    out = run_openmm_qmmm_distance_scan_probe(
        xyz_full_ang=xyz_full,
        atomic_numbers_full=atomic_numbers,
        mm_charges_e=charges,
        atom_rows=atom_rows,
        scan_pair=(1, 2),
        scan_start_ang=1.0,
        scan_end_ang=2.0,
        n_points=5,
        config={},
        qmmm_step_fn=_mock_step,
    )
    assert out["ok"] is True
    assert int(out["n_points"]) == 5
    rows = list(out.get("rows", []) or [])
    assert len(rows) == 5
    best = dict(out.get("best_point", {}) or {})
    assert int(best.get("point", 0)) in (3, 2, 4)
    assert abs(float(best.get("distance_ang", 0.0)) - 1.5) <= 0.26
