from pathlib import Path

from metaboclip.core.atoms import parse_pdbqt_poses, read_structure_atoms
from metaboclip.core.workflow import run_single_pair
from metaboclip.core.config import load_yaml
from metaboclip.core.validation import validate_mechanism_file


def test_ligand_pdbqt_order_uses_line_order_not_file_serial(tmp_path):
    pdbqt = tmp_path / "lig.pdbqt"
    pdbqt.write_text(
        "MODEL 1\n"
        "ATOM     99  C1  LIG A   1       1.000   0.000   0.000  0.00  0.00      A\n"
        "ATOM      7  O1  LIG A   1       2.000   0.000   0.000  0.00  0.00      OA\n"
        "ENDMDL\n",
        encoding="utf-8",
    )
    poses = parse_pdbqt_poses(pdbqt)
    assert len(poses) == 1
    atoms = poses[0]["atoms"]
    assert atoms[0].serial == 1
    assert atoms[1].serial == 2
    assert atoms[0].extra["file_serial"] == 99
    assert atoms[1].extra["file_serial"] == 7


def test_linked_ligand_site_requires_same_group_and_instance(tmp_path):
    protein = tmp_path / "protein.pdbqt"
    protein.write_text(
        "ATOM      1  OG  SER A  10       0.000   0.000   0.000  0.00  0.00      OA\n"
        "ATOM      2  NE2 HIS A  20       0.000   3.000   0.000  0.00  0.00      NA\n",
        encoding="utf-8",
    )
    docked = tmp_path / "lig@protein.pdbqt"
    docked.write_text(
        "MODEL 1\n"
        "REMARK VINA RESULT: -6.0 0.0 0.0\n"
        "ATOM     50  C1  LIG A   1       3.000   0.000   0.000  0.00  0.00      C\n"
        "ATOM     51  O1  LIG A   1       4.000   0.000   0.000  0.00  0.00      OA\n"
        "ATOM     52  O2  LIG A   1       8.000   0.000   0.000  0.00  0.00      OA\n"
        "ENDMDL\n",
        encoding="utf-8",
    )
    role_table = tmp_path / "lig.role_table.csv"
    role_table.write_text(
        "ligand_id,group_id,instance_id,atom_label,atom_class,atom_role,source_atom_index,element,pdbqt_order,subtype,confidence\n"
        "lig,acetyl,grp1,acetyl.acyl_c,acyl_carbonyl_c,acyl_c,1,C,1,,1.0\n"
        "lig,benzoyl,grp1,benzoyl.carbonyl_o,carbonyl_oxygen,carbonyl_o,2,O,2,,1.0\n"
        "lig,acetyl,grp1,acetyl.carbonyl_o,carbonyl_oxygen,carbonyl_o,3,O,3,,1.0\n",
        encoding="utf-8",
    )
    mechanism = {
        "ligand_sites": {
            "carbonyl_c": {"atom_classes": ["acyl_carbonyl_c"], "required": True},
            "carbonyl_o": {"linked_to": "carbonyl_c", "atom_classes": ["carbonyl_oxygen"], "required": True},
        },
        "protein_roles": [
            {"role": "nucleophile", "from": "ligand.carbonyl_c", "radius": 5.0, "required": True, "residues": {"SER": ["OG"]}},
        ],
        "features": {
            "attack_distance": {
                "type": "distance",
                "atoms": ["protein.nucleophile", "ligand.carbonyl_c"],
                "gate": {"min": 0.0, "max": 5.0, "required": True},
                "score": {"level": 1, "transform": "distance_piecewise", "best": 3.2, "cutoff": 5.0, "weight": 1.0},
            },
            "attack_angle": {
                "type": "angle_3pt",
                "enabled": True,
                "a": "protein.nucleophile",
                "vertex": "ligand.carbonyl_c",
                "c": "ligand.carbonyl_o",
                "gate": {"min": 0.0, "max": 180.0, "required": True},
                "score": {"level": 1, "transform": "angle_gaussian", "target": 180.0, "sigma": 20.0, "weight": 1.0},
            },
        },
    }
    profile = {
        "scoring": {"affinity": {"full": -7.0, "zero": -3.0, "weight": 0.30}, "geometry": {"weight": 0.70, "level_weights": {1: 1.0}}},
        "aggregation": {"protein": {"total_conformations": 1, "coverage_threshold": 0.0, "coverage_target": 1.0, "alpha": 0.0}},
    }
    out = tmp_path / "out"
    result = run_single_pair(mechanism, profile, protein, docked, role_table, out)
    assert result["candidates"] == 1
    geom = (out / "geometry_features.csv").read_text(encoding="utf-8")
    assert "attack_angle" in geom
    # The same group and instance acetyl oxygen is at pdbqt_order 3, not the benzoyl atom at order 2.
    ligand_rows = (out / "resolved_ligand_sites.csv").read_text(encoding="utf-8")
    assert "acetyl.carbonyl_o" in ligand_rows
    assert "benzoyl.carbonyl_o" not in ligand_rows


def test_mechanism_validator_accepts_builtin_cyp450():
    root = Path(__file__).resolve().parents[1]
    report = validate_mechanism_file(root / "metaboclip" / "config" / "families" / "cyp450" / "mechanism.yaml")
    assert report.ok, report.errors


def test_optional_missing_support_level_is_not_scored(tmp_path):
    protein = tmp_path / "protein.pdbqt"
    protein.write_text(
        "ATOM      1  OG  SER A  10       0.000   0.000   0.000  0.00  0.00      OA\n",
        encoding="utf-8",
    )
    docked = tmp_path / "lig@protein.pdbqt"
    docked.write_text(
        "MODEL 1\n"
        "REMARK VINA RESULT: -7.0 0.0 0.0\n"
        "ATOM      1  C1  LIG A   1       3.200   0.000   0.000  0.00  0.00      C\n"
        "ENDMDL\n",
        encoding="utf-8",
    )
    role_table = tmp_path / "lig.role_table.csv"
    role_table.write_text(
        "ligand_id,group_id,instance_id,atom_label,atom_class,atom_role,source_atom_index,element,pdbqt_order,subtype,confidence\n"
        "lig,ester,e1,ester.carbonyl_c,ester_carbonyl_c,carbonyl_c,1,C,1,,1.0\n",
        encoding="utf-8",
    )
    mechanism = {
        "ligand_sites": {"target": {"atom_classes": ["ester_carbonyl_c"], "required": True}},
        "protein_roles": [
            {"role": "atom1", "from": "ligand.target", "radius": 6.0, "required": True, "residues": {"SER": ["OG"]}},
            {"role": "atom2", "from": "protein.atom1", "radius": 2.0, "required": False, "residues": {"HIS": ["NE2"]}},
        ],
        "features": {
            "primary_distance": {
                "type": "distance",
                "atoms": ["protein.atom1", "ligand.target"],
                "gate": {"min": 0.0, "max": 6.0, "required": True},
                "score": {"level": 1, "transform": "distance_piecewise", "best": 3.2, "cutoff": 6.0, "weight": 1.0},
            },
            "support_distance": {
                "type": "distance",
                "atoms": ["protein.atom1", "protein.atom2"],
                "gate": {"min": 0.0, "max": 6.0, "required": False},
                "score": {"level": 2, "transform": "distance_piecewise", "best": 3.2, "cutoff": 6.0, "weight": 1.0},
            },
        },
    }
    profile = {
        "scoring": {"affinity": {"full": -7.0, "zero": -3.0, "weight": 0.30}, "geometry": {"weight": 0.70, "level_weights": {1: 0.60, 2: 0.25}}},
        "aggregation": {"protein": {"total_conformations": 1, "coverage_threshold": 0.0, "coverage_target": 1.0, "alpha": 0.0}},
    }
    out = tmp_path / "out_optional"
    result = run_single_pair(mechanism, profile, protein, docked, role_table, out)
    assert result["candidates"] == 1
    text = (out / "candidate_scores.csv").read_text(encoding="utf-8")
    assert "level_2_score" not in text
    assert "100.0" in text or "100" in text


def test_collection_role_min_count_and_same_residue(tmp_path):
    protein = tmp_path / "protein.pdbqt"
    protein.write_text(
        "HETATM    1 FE   HEM A 501       0.000   0.000   0.000  0.00  0.00      Fe\n"
        "HETATM    2 NA   HEM A 501       1.000   0.000   0.000  0.00  0.00      N\n"
        "HETATM    3 NB   HEM A 501       0.000   1.000   0.000  0.00  0.00      N\n"
        "HETATM    4 NC   HEM A 501      -1.000   0.000   0.000  0.00  0.00      N\n"
        "HETATM    5 ND   HEM A 501       0.000  -1.000   0.000  0.00  0.00      N\n",
        encoding="utf-8",
    )
    docked = tmp_path / "lig@protein.pdbqt"
    docked.write_text(
        "MODEL 1\n"
        "REMARK VINA RESULT: -7.0 0.0 0.0\n"
        "ATOM      1  C1  LIG A   1       3.000   0.000   0.000  0.00  0.00      C\n"
        "ENDMDL\n",
        encoding="utf-8",
    )
    role_table = tmp_path / "lig.role_table.csv"
    role_table.write_text(
        "ligand_id,group_id,instance_id,atom_label,atom_class,atom_role,source_atom_index,element,pdbqt_order,subtype,confidence\n"
        "lig,ch,c1,c_h_site.reactive_c,hydrogen_bearing_carbon,reactive_c,1,C,1,,1.0\n",
        encoding="utf-8",
    )
    mechanism = {
        "ligand_sites": {"reactive_c": {"atom_classes": ["hydrogen_bearing_carbon"], "required": True}},
        "protein_roles": [
            {"role": "fe", "from": "ligand.reactive_c", "radius": 5.0, "required": True, "atoms": ["elem Fe"]},
            {"role": "heme_n", "from": "protein.fe", "radius": 2.0, "required": True, "collection": True, "same_residue": True, "min_count": 3, "preferred_count": 4, "residues": {"HEM": ["NA", "NB", "NC", "ND"]}},
        ],
        "geometry_refs": {
            "heme_plane": {"type": "best_fit_plane", "atoms_from_collection": "protein.heme_n", "min_atoms": 3},
            "distal_axis": {"type": "oriented_axis", "origin": "protein.fe", "direction": {"type": "plane_normal", "plane": "geometry.heme_plane"}},
        },
        "features": {
            "fe_c_distance": {"type": "distance", "atoms": ["protein.fe", "ligand.reactive_c"], "gate": {"min": 0.0, "max": 5.0, "required": True}, "score": {"level": 1, "transform": "distance_piecewise", "best": 3.0, "cutoff": 5.0, "weight": 1.0}},
            "axis": {"type": "axis_deviation", "axis": "geometry.distal_axis", "vector": {"from": "protein.fe", "to": "ligand.reactive_c"}, "gate": {"min": 0.0, "max": 180.0, "required": True}, "score": {"level": 1, "transform": "angle_deviation_linear", "cutoff": 180.0, "weight": 1.0}},
        },
    }
    profile = {"scoring": {"affinity": {"full": -7.0, "zero": -3.0, "weight": 0.30}, "geometry": {"weight": 0.70, "level_weights": {1: 1.0}}}, "aggregation": {"protein": {"total_conformations": 1, "coverage_threshold": 0.0, "coverage_target": 1.0, "alpha": 0.0}}}
    out = tmp_path / "out_collection"
    result = run_single_pair(mechanism, profile, protein, docked, role_table, out)
    assert result["candidates"] == 1
    rows = (out / "resolved_protein_roles.csv").read_text(encoding="utf-8")
    assert rows.count("heme_n") == 4
