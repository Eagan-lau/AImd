#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .utils import as_bool, numeric, read_csv, resolve_path


@dataclass
class ProteinRecord:
    protein_id: str
    cluster_id: str
    batch_id: str
    protein_path: Path
    source_pdb: str = ""
    pdbqt_path: Path | None = None
    is_representative: bool = False


@dataclass
class LigandRecord:
    ligand_id: str
    batch_id: str
    ligand_path: Path
    smiles: str = ""
    name: str = ""


@dataclass
class PocketRecord:
    protein_id: str
    batch_id: str
    pocket_id: str
    pocket_rank: int
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    final_score: str = ""
    pocket_pdb_path: str = ""
    pocket_json_path: str = ""
    box_yaml_path: str = ""


def _first(row: dict[str, str], keys: list[str]) -> str:
    for key in keys:
        v = row.get(key, "")
        if v not in {None, ""}:
            return str(v)
    return ""


def load_proteins(config: dict[str, Any]) -> list[ProteinRecord]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    manifest_path = resolve_path(config.get("paths", {}).get("protein_manifest", "data/protein/protein_manifest.csv"), root)
    if manifest_path is None or not manifest_path.exists():
        raise FileNotFoundError(f"Protein manifest not found: {manifest_path}")
    rows = read_csv(manifest_path)
    if not rows:
        raise RuntimeError(f"Protein manifest is empty: {manifest_path}")

    selection = config.get("selection", {})
    mode = str(selection.get("protein_mode", "all"))
    max_n = selection.get("max_proteins")
    max_n = int(max_n) if max_n not in {None, "", 0, "0"} else None
    require_existing = bool(config.get("input", {}).get("require_existing_files", True))

    records: list[ProteinRecord] = []
    for row in rows:
        status = str(row.get("status", "success")).strip().lower()
        if status and status not in {"success", "ok", "true", "1"}:
            continue
        is_rep = as_bool(row.get("is_representative", False))
        if mode == "representatives_only" and not is_rep:
            continue
        protein_path_raw = _first(row, ["protein_path", "pdb_path", "source_pdb"])
        if not protein_path_raw:
            raise ValueError("protein_manifest must contain protein_path, pdb_path, or source_pdb")
        protein_path = resolve_path(protein_path_raw, root)
        if protein_path is None:
            continue
        if require_existing and not protein_path.exists():
            raise FileNotFoundError(f"Protein file not found for {row.get('protein_id')}: {protein_path}")
        pdbqt_raw = _first(row, ["pdbqt_path", "protein_pdbqt", "receptor_pdbqt"])
        pdbqt = resolve_path(pdbqt_raw, root) if pdbqt_raw else None
        if pdbqt is not None and require_existing and not pdbqt.exists():
            pdbqt = None
        records.append(
            ProteinRecord(
                protein_id=row.get("protein_id") or protein_path.stem,
                cluster_id=row.get("cluster_id", ""),
                batch_id=row.get("batch_id", "file_1") or "file_1",
                protein_path=protein_path,
                source_pdb=row.get("source_pdb", ""),
                pdbqt_path=pdbqt,
                is_representative=is_rep,
            )
        )
        if max_n is not None and len(records) >= max_n:
            break
    if not records:
        raise RuntimeError("No proteins selected for docking")
    return records


def load_ligands(config: dict[str, Any]) -> list[LigandRecord]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    ligand_manifest = resolve_path(config.get("paths", {}).get("ligand_manifest", "data/ligand/ligand_manifest.csv"), root)
    require_existing = bool(config.get("input", {}).get("require_existing_files", True))
    records: list[LigandRecord] = []
    if ligand_manifest is not None and ligand_manifest.exists():
        rows = read_csv(ligand_manifest)
        for row in rows:
            status = str(row.get("status", "success")).strip().lower()
            if status and status not in {"success", "ok", "true", "1"}:
                continue
            ligand_path_raw = _first(row, ["ligand_path", "pdbqt_path", "ligand_pdbqt", "pdbqt"])
            if not ligand_path_raw:
                continue
            ligand_path = resolve_path(ligand_path_raw, root)
            if ligand_path is None:
                continue
            if require_existing and not ligand_path.exists():
                raise FileNotFoundError(f"Ligand file not found for {row.get('ligand_id')}: {ligand_path}")
            records.append(
                LigandRecord(
                    ligand_id=row.get("ligand_id") or ligand_path.stem,
                    batch_id=row.get("batch_id", "file_1") or "file_1",
                    ligand_path=ligand_path,
                    smiles=row.get("smiles", ""),
                    name=row.get("name", ""),
                )
            )
    else:
        ligand_dir = resolve_path(config.get("paths", {}).get("ligand_dir", "data/ligand"), root)
        if ligand_dir and ligand_dir.exists():
            for batch_dir in sorted(ligand_dir.glob("file_*")):
                if not batch_dir.is_dir():
                    continue
                for fp in sorted(batch_dir.glob("*.pdbqt")):
                    records.append(LigandRecord(ligand_id=fp.stem, batch_id=batch_dir.name, ligand_path=fp))
    if not records:
        raise RuntimeError("No ligand PDBQT files found. Provide data/ligand/ligand_manifest.csv or data/ligand/file_*/*.pdbqt")
    return records


def load_pockets(config: dict[str, Any]) -> list[PocketRecord]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    pocket_manifest = resolve_path(config.get("paths", {}).get("pocket_manifest", "data/pocket/pocket_manifest.csv"), root)
    if pocket_manifest is None or not pocket_manifest.exists():
        raise FileNotFoundError(f"Pocket manifest not found: {pocket_manifest}")
    rows = read_csv(pocket_manifest)
    if not rows:
        raise RuntimeError(f"Pocket manifest is empty: {pocket_manifest}")
    records: list[PocketRecord] = []
    for row in rows:
        status = str(row.get("status", "success")).strip().lower()
        if status and status not in {"success", "ok", "true", "1"}:
            continue
        cx, cy, cz = numeric(row.get("center_x")), numeric(row.get("center_y")), numeric(row.get("center_z"))
        sx, sy, sz = numeric(row.get("size_x")), numeric(row.get("size_y")), numeric(row.get("size_z"))
        if None in {cx, cy, cz, sx, sy, sz}:
            continue
        rank = int(float(row.get("pocket_rank", row.get("rank", 1)) or 1))
        records.append(
            PocketRecord(
                protein_id=row.get("protein_id") or row.get("query_id") or "",
                batch_id=row.get("batch_id", "file_1") or "file_1",
                pocket_id=row.get("pocket_id", f"P{rank}"),
                pocket_rank=rank,
                center_x=float(cx), center_y=float(cy), center_z=float(cz),
                size_x=float(sx), size_y=float(sy), size_z=float(sz),
                final_score=row.get("final_score", row.get("pocket_score", "")),
                pocket_pdb_path=row.get("pocket_pdb_path", ""),
                pocket_json_path=row.get("pocket_json_path", ""),
                box_yaml_path=row.get("box_yaml_path", ""),
            )
        )
    if not records:
        raise RuntimeError("No valid pockets found in pocket_manifest.csv")
    return records


def group_pockets_by_protein(pockets: list[PocketRecord], top_n: int | None = None) -> dict[str, list[PocketRecord]]:
    grouped: dict[str, list[PocketRecord]] = {}
    for pocket in pockets:
        grouped.setdefault(pocket.protein_id, []).append(pocket)
    for pid in list(grouped):
        grouped[pid] = sorted(grouped[pid], key=lambda p: (p.pocket_rank, p.pocket_id))
        if top_n:
            grouped[pid] = grouped[pid][:top_n]
    return grouped
