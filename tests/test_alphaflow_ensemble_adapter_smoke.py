from pathlib import Path

from DockingHub.ensemble import build_conformer_manifest
from DockingHub.manifest import ProteinRecord


def test_alphaflow_adapter_stages_csv_and_template(tmp_path: Path) -> None:
    pdb = tmp_path / "protein.pdb"
    pdb.write_text(
        "\n".join(
            [
                "ATOM      1  N   ALA A   1      11.104  13.207  12.011  1.00 20.00           N",
                "ATOM      2  CA  ALA A   1      12.000  13.000  12.500  1.00 20.00           C",
                "ATOM      3  N   GLY A   2      13.104  14.207  13.011  1.00 20.00           N",
                "ATOM      4  CA  GLY A   2      14.000  14.000  13.500  1.00 20.00           C",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    config = {
        "paths": {"aimd_root": str(tmp_path), "ensemble_dir": "ensemble"},
        "ensemble": {
            "enabled": True,
            "engine": "alphaflow",
            "fallback_to_input": True,
            "max_conformers_per_protein": 2,
            "alphaflow": {
                "project_dir": "",
                "run_msa": False,
                "run_prediction": False,
                "require_msa": False,
                "template_mode": "copy_input",
            },
            "align_to_reference": {"enabled": False},
        },
    }
    protein = ProteinRecord(
        protein_id="protein_A",
        cluster_id="C000001",
        batch_id="file_1",
        protein_path=pdb,
    )

    rows = build_conformer_manifest(config, [protein])

    assert len(rows) == 1
    assert rows[0]["status"] == "success"
    assert rows[0]["source"] == "input_protein"
    staged_csv = tmp_path / "ensemble" / "file_1" / "protein_A" / "input_csv" / "protein_A.csv"
    staged_template = tmp_path / "ensemble" / "file_1" / "protein_A" / "templates_dir" / "protein_A.pdb"
    assert staged_csv.read_text(encoding="utf-8").splitlines() == ["name,seqres", "protein_A,AG"]
    assert staged_template.exists()
