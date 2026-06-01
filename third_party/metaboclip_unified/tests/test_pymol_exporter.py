from pathlib import Path
import csv

from metaboclip.pymol_exporter.export import export_pymol_script


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def test_export_writes_candidate_selections(tmp_path: Path):
    protein = tmp_path / "protein.pdbqt"
    protein.write_text(
        "ATOM      1  OG  SER A 171      10.000   8.000  -1.000  0.00  0.00           O\n"
        "ATOM      2  NE2 HIS A 340      11.000   8.000  -1.000  0.00  0.00           N\n",
        encoding="utf-8",
    )
    ligand = tmp_path / "ligand.pdbqt"
    ligand.write_text(
        "MODEL        1\n"
        "HETATM    7  C1  LIG L   1       1.000   2.000   3.000  0.00  0.00           C\n"
        "HETATM    8  O1  LIG L   1       2.000   2.000   3.000  0.00  0.00           O\n"
        "ENDMDL\n",
        encoding="utf-8",
    )
    candidate = tmp_path / "candidate_scores.csv"
    write_csv(candidate, [{
        "protein_id": "P1", "conformation_id": "0", "ligand_id": "L1", "pose_id": "1",
        "site_set_id": "pose_1_set_1", "candidate_score": "88.0",
    }])
    lig_rows = tmp_path / "resolved_ligand_sites.csv"
    write_csv(lig_rows, [{
        "protein_id": "P1", "conformation_id": "0", "ligand_id": "L1", "pose_id": "1",
        "site_set_id": "pose_1_set_1", "ligand_site": "carbonyl_c", "group_id": "ester",
        "instance_id": "ester_1", "atom_label": "ester.carbonyl_c", "atom_class": "ester_carbonyl_c",
        "atom_role": "carbonyl_c", "pdbqt_order": "1", "file_serial": "7", "element": "C",
        "x": "1.0", "y": "2.0", "z": "3.0",
    }])
    prot_rows = tmp_path / "resolved_protein_roles.csv"
    write_csv(prot_rows, [{
        "protein_id": "P1", "conformation_id": "0", "ligand_id": "L1", "pose_id": "1",
        "site_set_id": "pose_1_set_1", "protein_role": "nucleophile", "role_atom_index": "1",
        "chain": "A", "resi": "171", "resn": "SER", "atom_name": "OG", "element": "O",
        "x": "10.0", "y": "8.0", "z": "-1.0", "distance_to_anchor": "3.0",
    }])
    mechanism = tmp_path / "mechanism.yaml"
    mechanism.write_text(
        "features:\n"
        "  attack_distance:\n"
        "    type: distance\n"
        "    atoms:\n"
        "      - protein.nucleophile\n"
        "      - ligand.carbonyl_c\n",
        encoding="utf-8",
    )
    out_pml = tmp_path / "view.pml"
    export_pymol_script(
        protein_file=protein,
        docked_ligand_file=ligand,
        resolved_ligand_sites=lig_rows,
        resolved_protein_roles=prot_rows,
        candidate_scores=candidate,
        mechanism_file=mechanism,
        out_pml=out_pml,
        save_pse_name=str(tmp_path / "view.pse"),
    )
    text = out_pml.read_text(encoding="utf-8")
    assert "select LIGAND_SITE_carbonyl_c" in text
    assert "select PROTEIN_ROLE_nucleophile" in text
    assert "select CANDIDATE_LIGAND_ATOMS" in text
    assert "select CANDIDATE_PROTEIN_ATOMS" in text
    assert "create CANDIDATE_LIGAND_ATOMS_OBJECT" in text
    assert "distance DIST_attack_distance" in text
    assert "cmd.save(_pse_path)" in text
