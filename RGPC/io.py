#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Any

from .utils import resolve_path


@dataclass
class ProteinRecord:
    protein_id: str
    pdb_path: Path
    source: str = "structure_library"
    status: str = "success"


def scan_structure_dir(structure_dir: Path, extensions: List[str]) -> List[ProteinRecord]:
    records: List[ProteinRecord] = []
    ext_set = {e.lower() for e in extensions}
    for p in sorted(structure_dir.iterdir()):
        if not p.is_file():
            continue
        if p.suffix.lower() not in ext_set:
            continue
        protein_id = p.stem
        records.append(ProteinRecord(protein_id=protein_id, pdb_path=p.resolve()))
    return records


def read_structure_manifest(manifest_path: Path, root: Path, require_existing: bool = True) -> List[ProteinRecord]:
    records: List[ProteinRecord] = []
    with open(manifest_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        if "protein_id" not in reader.fieldnames or "pdb_path" not in reader.fieldnames:
            raise ValueError("structure_manifest must contain at least protein_id and pdb_path columns")
        for row in reader:
            status = row.get("status", "success")
            if status and status.lower() not in {"success", "ok", "true", "1"}:
                continue
            p = resolve_path(row["pdb_path"], root)
            if p is None:
                continue
            if require_existing and not p.exists():
                raise FileNotFoundError(f"PDB path not found for {row['protein_id']}: {p}")
            records.append(ProteinRecord(protein_id=row["protein_id"], pdb_path=p.resolve(), status="success"))
    return records


def load_protein_records(config: Dict[str, Any]) -> List[ProteinRecord]:
    paths = config.get("paths", {})
    input_cfg = config.get("input", {})
    root = Path(paths.get("aimd_root", ".")).resolve()
    structure_dir = resolve_path(paths.get("structure_dir"), root)
    manifest = resolve_path(paths.get("structure_manifest"), root)
    require_existing = bool(input_cfg.get("require_existing_files", True))
    extensions = input_cfg.get("file_extensions", [".pdb"])

    if manifest and manifest.exists():
        records = read_structure_manifest(manifest, root=root, require_existing=require_existing)
    else:
        if structure_dir is None or not structure_dir.exists():
            raise FileNotFoundError(f"structure_dir not found: {structure_dir}")
        records = scan_structure_dir(structure_dir, extensions)

    if not records:
        raise RuntimeError("No protein structure files found for RGPC")

    seen = set()
    unique = []
    for r in records:
        if r.protein_id in seen:
            raise ValueError(f"Duplicate protein_id detected: {r.protein_id}")
        seen.add(r.protein_id)
        unique.append(r)
    return unique


def write_input_manifest(records: List[ProteinRecord], path: Path) -> None:
    from .utils import write_csv
    rows = [
        {
            "protein_id": r.protein_id,
            "pdb_path": str(r.pdb_path),
            "source": r.source,
            "status": r.status,
        }
        for r in records
    ]
    write_csv(path, rows, ["protein_id", "pdb_path", "source", "status"])
