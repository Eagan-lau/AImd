#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .utils import as_bool, read_csv, resolve_path


@dataclass
class ProteinInput:
    protein_id: str
    cluster_id: str
    batch_id: str
    protein_path: Path
    source_pdb: str = ""
    is_representative: bool = False
    status: str = "success"


def _first_existing_key(row: dict[str, str], keys: list[str]) -> str:
    for key in keys:
        value = row.get(key, "")
        if value:
            return value
    return ""


def load_protein_manifest(config: dict[str, Any]) -> list[ProteinInput]:
    paths = config.get("paths", {})
    input_cfg = config.get("input", {})
    selection_cfg = config.get("selection", {})
    root = Path(paths.get("aimd_root", ".")).resolve()
    manifest_path = resolve_path(paths.get("protein_manifest", "data/data_output/protein_batches/protein_manifest.csv"), root)
    if manifest_path is None or not manifest_path.exists():
        raise FileNotFoundError(f"Protein manifest not found: {manifest_path}")

    required = {"protein_id", "batch_id"}
    rows = read_csv(manifest_path)
    if not rows:
        raise RuntimeError(f"Protein manifest is empty: {manifest_path}")
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"Protein manifest is missing required columns: {sorted(missing)}")

    require_existing = bool(input_cfg.get("require_existing_files", True))
    mode = str(selection_cfg.get("mode", "all"))
    max_proteins = selection_cfg.get("max_proteins")
    max_proteins = int(max_proteins) if max_proteins not in {None, "", 0, "0"} else None

    records: list[ProteinInput] = []
    for row in rows:
        status = str(row.get("status", "success")).strip().lower()
        if status and status not in {"success", "ok", "true", "1"}:
            continue
        is_rep = as_bool(row.get("is_representative", False))
        if mode == "representatives_only" and not is_rep:
            continue
        protein_rel = _first_existing_key(row, ["protein_path", "pdb_path", "source_pdb"])
        if not protein_rel:
            raise ValueError("Protein manifest must contain protein_path, pdb_path, or source_pdb")
        protein_path = resolve_path(protein_rel, root)
        if protein_path is None:
            continue
        if require_existing and not protein_path.exists():
            raise FileNotFoundError(f"Protein file not found for {row.get('protein_id')}: {protein_path}")
        records.append(
            ProteinInput(
                protein_id=row["protein_id"],
                cluster_id=row.get("cluster_id", ""),
                batch_id=row.get("batch_id", "file_1") or "file_1",
                protein_path=protein_path,
                source_pdb=row.get("source_pdb", ""),
                is_representative=is_rep,
                status="success",
            )
        )
        if max_proteins is not None and len(records) >= max_proteins:
            break

    if not records:
        raise RuntimeError("No protein entries selected for TApocket batch prediction")
    return records
