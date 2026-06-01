from pathlib import Path

from metaboclip.core.workflow import run_directory


def test_generic_hydrolase_run(tmp_path):
    root = Path(__file__).resolve().parent
    result = run_directory(
        mechanism_path=Path(__file__).resolve().parents[1] / "examples" / "generic_hydrolase" / "mechanism.yaml",
        profile_path=Path(__file__).resolve().parents[1] / "metaboclip" / "config" / "profiles" / "default_profile.yaml",
        protein_dir=root / "data" / "proteins",
        docking_dir=root / "data" / "docking",
        role_table_dir=root / "data" / "ligand_roles",
        out_dir=tmp_path,
        ligand_id="lig001",
    )
    assert result["pairs"] == 1
    assert result["candidates"] >= 1
    assert (tmp_path / "lig001" / "ProteinA_conf1" / "candidate_scores.csv").exists()


def test_protein_score_aggregates_conformations(tmp_path):
    root = Path(__file__).resolve().parent
    src_protein = root / "data" / "proteins" / "ProteinA_conf1.pdbqt"
    src_pose = root / "data" / "docking" / "lig001@ProteinA_conf1.pdbqt"
    src_out = root / "data" / "docking" / "lig001_ProteinA_conf1_cavity_1.out"
    protein_dir = tmp_path / "proteins"
    docking_dir = tmp_path / "docking"
    role_dir = tmp_path / "roles"
    protein_dir.mkdir()
    docking_dir.mkdir()
    role_dir.mkdir()
    (role_dir / "lig001.role_table.csv").write_text((root / "data" / "ligand_roles" / "lig001.role_table.csv").read_text(), encoding="utf-8")
    for name in ["ProteinA", "ProteinA_1"]:
        (protein_dir / f"{name}.pdbqt").write_text(src_protein.read_text(), encoding="utf-8")
        (docking_dir / f"lig001@{name}.pdbqt").write_text(src_pose.read_text(), encoding="utf-8")
        (docking_dir / f"lig001_{name}_cavity_1.out").write_text(src_out.read_text(), encoding="utf-8")

    result = run_directory(
        mechanism_path=Path(__file__).resolve().parents[1] / "examples" / "generic_hydrolase" / "mechanism.yaml",
        profile_path=Path(__file__).resolve().parents[1] / "metaboclip" / "config" / "profiles" / "default_profile.yaml",
        protein_dir=protein_dir,
        docking_dir=docking_dir,
        role_table_dir=role_dir,
        out_dir=tmp_path / "out",
        ligand_id="lig001",
    )
    assert result["pairs"] == 2
    assert result["proteins"] == 1
    protein_scores = (tmp_path / "out" / "protein_scores.csv").read_text(encoding="utf-8")
    assert "protein_id,ligand_id,protein_score" in protein_scores
    assert "ProteinA,lig001" in protein_scores
    assert "observed_conformations" in protein_scores
