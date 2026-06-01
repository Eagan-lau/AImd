from __future__ import annotations

from pathlib import Path
from metaboclip.core.config import package_path
from metaboclip.core.workflow import run_directory


def run(profile_path: Path, protein_dir: Path, docking_dir: Path, role_table_dir: Path, out_dir: Path, ligand_id: str | None = None) -> dict:
    mechanism_path = package_path("config", "families", "fe2og", "mechanism.yaml")
    return run_directory(mechanism_path, profile_path, protein_dir, docking_dir, role_table_dir, out_dir, ligand_id)
