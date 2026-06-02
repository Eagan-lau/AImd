#!/usr/bin/env python3
"""Create the canonical AImd data_input/data_output layout.

The script copies legacy local input folders into the clean layout and writes
input manifests that use project-root-relative paths. It does not move or
delete existing files unless explicitly requested by overwriting manifests.
"""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


PROTEIN_EXTENSIONS = {".pdb", ".cif", ".mmcif"}
LIGAND_EXTENSIONS = {".pdbqt"}

OUTPUT_DIRS = [
    "protein_structure",
    "cluster",
    "protein_batches",
    "ligand_transformation",
    "pocket",
    "ensemble",
    "cofactor_mapped",
    "receptor",
    "docking_configs",
    "docking_tasks",
    "docking_out",
    "scoring",
    "refinement",
    "refined",
    "metaboclip",
]


def _relative(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _copy_tree_contents(src: Path, dst: Path, overwrite: bool, dry_run: bool) -> None:
    if not src.exists():
        return
    for item in sorted(src.rglob("*")):
        rel = item.relative_to(src)
        if item.is_file() and item.name.endswith("_manifest.csv"):
            continue
        target = dst / rel
        if item.is_dir():
            if not dry_run:
                target.mkdir(parents=True, exist_ok=True)
            continue
        if target.exists() and not overwrite:
            continue
        print(f"[AImd data] copy {item} -> {target}")
        if not dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(item, target)


def _discover_files(base: Path, extensions: set[str]) -> list[Path]:
    if not base.exists():
        return []
    return sorted(
        path for path in base.rglob("*")
        if path.is_file() and path.suffix.lower() in extensions
    )


def _batch_id(path: Path, base: Path) -> str:
    rel_parts = path.relative_to(base).parts
    if len(rel_parts) > 1 and rel_parts[0].startswith("file_"):
        return rel_parts[0]
    return "file_1"


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]], overwrite: bool, dry_run: bool) -> None:
    if path.exists() and not overwrite:
        print(f"[AImd data] keep existing manifest: {path}")
        return
    print(f"[AImd data] write manifest: {path}")
    if dry_run:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_protein_manifest(root: Path, protein_dir: Path, overwrite: bool, dry_run: bool) -> None:
    rows = []
    for path in _discover_files(protein_dir, PROTEIN_EXTENSIONS):
        rows.append({
            "protein_id": path.stem,
            "pdb_path": _relative(path, root),
            "batch_id": _batch_id(path, protein_dir),
            "status": "success",
        })
    _write_csv(
        protein_dir / "structure_manifest.csv",
        ["protein_id", "pdb_path", "batch_id", "status"],
        rows,
        overwrite,
        dry_run,
    )


def _write_ligand_manifest(root: Path, ligand_dir: Path, overwrite: bool, dry_run: bool) -> None:
    rows = []
    for path in _discover_files(ligand_dir, LIGAND_EXTENSIONS):
        rel_path = _relative(path, root)
        rows.append({
            "ligand_id": path.stem,
            "batch_id": _batch_id(path, ligand_dir),
            "ligand_path": rel_path,
            "pdbqt_path": rel_path,
            "status": "success",
        })
    _write_csv(
        ligand_dir / "ligand_manifest.csv",
        ["ligand_id", "batch_id", "ligand_path", "pdbqt_path", "status"],
        rows,
        overwrite,
        dry_run,
    )


def prepare_layout(root: Path, overwrite: bool, dry_run: bool) -> None:
    data_input = root / "data" / "data_input"
    data_output = root / "data" / "data_output"

    for rel in [
        "protein/file_1",
        "ligand/file_1",
        "cofactor/file_1",
        "workflow",
    ]:
        path = data_input / rel
        print(f"[AImd data] ensure input dir: {path}")
        if not dry_run:
            path.mkdir(parents=True, exist_ok=True)

    for rel in OUTPUT_DIRS:
        path = data_output / rel
        print(f"[AImd data] ensure output dir: {path}")
        if not dry_run:
            path.mkdir(parents=True, exist_ok=True)

    _copy_tree_contents(root / "data" / "protein_input", data_input / "protein", overwrite, dry_run)
    _copy_tree_contents(root / "data" / "ligand_input", data_input / "ligand", overwrite, dry_run)
    _copy_tree_contents(root / "data" / "input", data_input / "workflow", overwrite, dry_run)
    _copy_tree_contents(root / "data" / "cofactor", data_input / "cofactor", overwrite, dry_run)

    _write_protein_manifest(root, data_input / "protein", overwrite=True, dry_run=dry_run)
    _write_ligand_manifest(root, data_input / "ligand", overwrite=True, dry_run=dry_run)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare the canonical AImd data layout")
    parser.add_argument("--root", default=".", help="AImd project root")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite copied files when targets already exist")
    parser.add_argument("--dry-run", action="store_true", help="Print actions without writing files")
    args = parser.parse_args()

    prepare_layout(Path(args.root).resolve(), overwrite=bool(args.overwrite), dry_run=bool(args.dry_run))


if __name__ == "__main__":
    main()
