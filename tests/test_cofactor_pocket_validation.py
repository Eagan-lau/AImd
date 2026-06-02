from __future__ import annotations

import csv
from pathlib import Path

from DockingHub.cofactor import (
    compute_cofactor_pocket_validation,
    get_pocket_validation_config,
    build_cofactor_manifest,
)


def _atom(record: str, name: str, resn: str, chain: str, resi: int, coord: tuple[float, float, float]) -> dict[str, object]:
    return {
        "record": record,
        "name": name,
        "resn": resn,
        "chain": chain,
        "resi": str(resi),
        "icode": "",
        "element": "C" if name != "N" else "N",
        "coord": coord,
    }


def _template_atoms() -> list[dict[str, object]]:
    atoms = [_atom("HETATM", "C1", "H5L", "A", 1001, (0.0, 0.0, 0.0))]
    for index in range(1, 6):
        atoms.append(_atom("ATOM", "CA", "ALA", "A", index, (float(index), 0.0, 0.0)))
    return atoms


def _alignment_pairs(offset: float = 0.0, count: int = 5) -> list[dict[str, object]]:
    pairs: list[dict[str, object]] = []
    for index in range(1, count + 1):
        template_coord = (float(index), 0.0, 0.0)
        target_coord = (float(index) + offset, 0.0, 0.0)
        pairs.append(
            {
                "template": {
                    "name": "CA",
                    "resn": "ALA",
                    "chain": "A",
                    "resi": str(index),
                    "icode": "",
                    "coord": template_coord,
                },
                "target": {
                    "name": "CA",
                    "resn": "ALA",
                    "chain": "A",
                    "resi": str(index),
                    "icode": "",
                    "coord": target_coord,
                },
            }
        )
    return pairs


def _validation_config(**overrides: object) -> dict[str, object]:
    config = get_pocket_validation_config({"cofactor": {"alignment": {"pocket_validation": overrides}}})
    config["min_mapped_residues"] = int(overrides.get("min_mapped_residues", 5))
    return config


def test_pocket_validation_defaults_are_applied_when_missing() -> None:
    config = get_pocket_validation_config({"cofactor": {"alignment": {}}})
    assert config["enabled"] is True
    assert config["pocket_radius"] == 5.0
    assert config["rmsd_metric"] == "ca"
    assert config["local_rmsd_cutoff"] == 2.5
    assert config["min_mapped_residues"] == 5
    assert config["min_pocket_coverage"] == 0.70
    assert config["require_pass_before_transfer"] is True


def test_local_rmsd_pass_allows_transfer() -> None:
    result = compute_cofactor_pocket_validation(
        _template_atoms(),
        _alignment_pairs(offset=0.2),
        ["H5L"],
        _validation_config(),
    )
    assert result["cofactor_transfer_pass"] is True
    assert result["cofactor_validation_status"] == "success"
    assert result["mapped_pocket_residue_count"] == 5


def test_local_rmsd_fail_prevents_transfer() -> None:
    result = compute_cofactor_pocket_validation(
        _template_atoms(),
        _alignment_pairs(offset=3.0),
        ["H5L"],
        _validation_config(),
    )
    assert result["cofactor_transfer_pass"] is False
    assert result["cofactor_validation_status"] == "high_local_rmsd"
    assert "cutoff" in result["cofactor_validation_message"]


def test_low_pocket_coverage_prevents_transfer() -> None:
    result = compute_cofactor_pocket_validation(
        _template_atoms(),
        _alignment_pairs(offset=0.1, count=3),
        ["H5L"],
        _validation_config(min_mapped_residues=2, min_pocket_coverage=0.80),
    )
    assert result["cofactor_transfer_pass"] is False
    assert result["cofactor_validation_status"] == "low_pocket_coverage"
    assert result["mapped_pocket_residue_count"] == 3


def test_missing_cofactor_atoms_prevents_transfer_with_clear_message() -> None:
    result = compute_cofactor_pocket_validation(
        _template_atoms(),
        _alignment_pairs(offset=0.1),
        ["HEM"],
        _validation_config(),
    )
    assert result["cofactor_transfer_pass"] is False
    assert result["cofactor_validation_status"] == "missing_cofactor_atoms"
    assert "No configured cofactor heavy atoms" in result["cofactor_validation_message"]


def test_cofactor_manifest_contains_new_fields(tmp_path: Path) -> None:
    target = tmp_path / "target.pdb"
    target.write_text("ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00           C\nEND\n", encoding="utf-8")
    config = {
        "paths": {
            "aimd_root": str(tmp_path),
            "cofactor_mapped_dir": "cofactor_mapped",
            "cofactor_dir": "cofactor",
        },
        "cofactor": {"enabled": False},
    }
    rows = build_cofactor_manifest(
        config,
        [
            {
                "protein_id": "P1",
                "batch_id": "file_1",
                "conformer_id": "conf_0",
                "structure_path": str(target),
                "status": "success",
            }
        ],
    )
    assert rows[0]["transfer_mode"] == "skipped_cofactor_disabled"
    manifest = tmp_path / "cofactor_mapped" / "cofactor_manifest.csv"
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle))
    for field in [
        "cofactor_pocket_radius",
        "cofactor_pocket_residue_count",
        "mapped_pocket_residue_count",
        "cofactor_pocket_coverage",
        "cofactor_site_ca_rmsd",
        "cofactor_local_rmsd_cutoff",
        "cofactor_transfer_pass",
        "cofactor_validation_status",
        "cofactor_validation_message",
        "transfer_mode",
    ]:
        assert field in header
