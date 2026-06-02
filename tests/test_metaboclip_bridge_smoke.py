from __future__ import annotations

import csv
import json
from pathlib import Path

from MetaBoClipBridge.bridge import _role_groups
from MetaBoClipBridge.main import run_metaboclip_bridge


DOCKING_FIELDS = [
    "job_id",
    "ligand_id",
    "protein_id",
    "cluster_id",
    "batch_id",
    "conformer_id",
    "pocket_id",
    "pocket_rank",
    "receptor_pdbqt_path",
    "ligand_pdbqt_path",
    "config_path",
    "out_pose_path",
    "log_path",
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
    "status",
    "return_code",
    "message",
    "best_affinity",
    "affinities",
    "n_affinities",
    "grid_size",
    "grid_space",
    "exhaustiveness",
    "random_seed",
    "pose_exists",
]


def _write_csv(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DOCKING_FIELDS)
        writer.writeheader()
        writer.writerow(row)


def _write_yaml(path: Path, config: dict) -> None:
    import yaml

    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")


def test_metaboclip_bridge_accepts_empty_role_group_list() -> None:
    assert _role_groups({"role_tables": {"groups": []}}) == []
    assert _role_groups({"role_tables": {"groups": ""}}) is None
    assert _role_groups({"role_tables": {"groups": ["hydroxyl", ""]}}) == ["hydroxyl"]


def test_metaboclip_bridge_runs_unified_backend_with_tiny_fixture(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    fixture_root = root / "metaboclip_unified" / "tests" / "data"
    manifest = tmp_path / "docking_result_manifest.csv"
    receptor = fixture_root / "proteins" / "ProteinA_conf1.pdbqt"
    pose = fixture_root / "docking" / "lig001@ProteinA_conf1.pdbqt"
    log = fixture_root / "docking" / "lig001_ProteinA_conf1_cavity_1.out"

    _write_csv(
        manifest,
        {
            "job_id": "smoke_lig001_ProteinA",
            "ligand_id": "lig001",
            "protein_id": "ProteinA",
            "cluster_id": "C000001",
            "batch_id": "file_1",
            "conformer_id": "conf_1",
            "pocket_id": "P1",
            "pocket_rank": "1",
            "receptor_pdbqt_path": str(receptor),
            "ligand_pdbqt_path": str(pose),
            "config_path": "",
            "out_pose_path": str(pose),
            "log_path": str(log),
            "center_x": "0",
            "center_y": "0",
            "center_z": "0",
            "size_x": "10",
            "size_y": "10",
            "size_z": "10",
            "status": "success",
            "return_code": "0",
            "message": "smoke fixture",
            "best_affinity": "-7.5",
            "affinities": "-7.5",
            "n_affinities": "1",
            "grid_size": "",
            "grid_space": "",
            "exhaustiveness": "",
            "random_seed": "1",
            "pose_exists": "true",
        },
    )

    config_path = tmp_path / "metaboclip_bridge.yaml"
    _write_yaml(
        config_path,
        {
            "module": "MetaBoClipBridge",
            "name": "Unified MetaboClip bridge smoke test",
            "backend": "unified",
            "paths": {
                "aimd_root": str(root),
                "refined_docking_manifest": str(manifest),
                "metaboclip_project_dir": "metaboclip_unified",
                "metaboclip_profile": "metaboclip_unified/metaboclip/config/profiles/default_profile.yaml",
                "staging_dir": str(tmp_path / "staging"),
                "unified_output_dir": str(tmp_path / "unified_runs"),
                "results_dir": str(tmp_path / "results"),
                "ligand_manifest": str(tmp_path / "missing_ligand_manifest.csv"),
                "ligand_source_manifest": str(tmp_path / "missing_ligand_source_manifest.csv"),
                "role_table_dir": "metaboclip_unified/tests/data/ligand_roles",
                "annotation_dir": str(tmp_path / "annotations"),
                "atom_map_dir": str(tmp_path / "atom_maps"),
            },
            "filtering": {
                "require_success_status": True,
                "require_existing_pose_and_log": True,
                "max_affinity_kcal": None,
            },
            "family_assignment": {
                "mode": "fixed",
                "fixed_families": ["hydrolase"],
                "all_families": ["hydrolase"],
            },
            "mechanisms": {
                "hydrolase": "metaboclip_unified/examples/generic_hydrolase/mechanism.yaml",
            },
            "role_tables": {
                "mode": "existing",
                "rules": "metaboclip_unified/rules/functional_groups.yaml",
                "groups": [],
                "ligand_source_column": "ligand_source_path",
                "prepared_pdbqt_column": "ligand_pdbqt_path",
                "max_dist": 0.5,
            },
            "execution": {
                "run_metaboclip": True,
                "continue_on_error": False,
            },
            "output": {
                "file_action": "copy",
                "overwrite": True,
            },
        },
    )

    final_path = run_metaboclip_bridge(config_path)
    results_dir = tmp_path / "results"

    assert final_path == results_dir / "metaboclip_final_ranking.csv"
    assert final_path.exists()
    assert (results_dir / "metaboclip_protein_scores_all.csv").exists()
    assert (results_dir / "metaboclip_conformation_scores_all.csv").exists()
    assert (results_dir / "metaboclip_pose_scores_all.csv").exists()
    assert (results_dir / "metaboclip_candidate_scores_all.csv").exists()
    assert (results_dir / "metaboclip_geometry_features_all.csv").exists()
    assert (results_dir / "metaboclip_resolved_ligand_sites_all.csv").exists()
    assert (results_dir / "metaboclip_resolved_protein_roles_all.csv").exists()

    rows = list(csv.DictReader(final_path.open("r", encoding="utf-8", newline="")))
    assert len(rows) == 1
    assert rows[0]["status"] == "success"
    assert rows[0]["protein_id"] == "ProteinA"
    assert rows[0]["ligand_id"] == "lig001"
    assert rows[0]["protein_score"]
    assert rows[0]["quality_score"]
    assert rows[0]["coverage"]
    assert "protein_score_norm" not in rows[0]
    assert "max_s_r" not in rows[0]

    report = json.loads((results_dir / "metaboclip_report.json").read_text(encoding="utf-8"))
    assert report["backend"] == "unified"
    assert report["scoring_executed"] is True
    assert report["n_success_rows"] == 1
