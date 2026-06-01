#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Dict, List, Any

from .io import ProteinRecord
from .clustering import choose_representative
from .utils import ensure_dir, write_csv, resolve_path


def make_batches_by_cluster(clusters: Dict[str, List[str]], max_proteins: int) -> Dict[str, str]:
    """Assign each cluster to file_*, keeping clusters together.

    If a single cluster is larger than max_proteins, it is still kept in one batch.
    """
    cluster_to_batch: Dict[str, str] = {}
    batch_idx = 1
    current_count = 0
    for cid, members in sorted(clusters.items(), key=lambda x: (-len(x[1]), x[0])):
        n = len(members)
        if current_count > 0 and current_count + n > max_proteins:
            batch_idx += 1
            current_count = 0
        cluster_to_batch[cid] = f"file_{batch_idx}"
        current_count += n
    return cluster_to_batch


def copy_or_link(src: Path, dst: Path, action: str) -> None:
    ensure_dir(dst.parent)
    if dst.exists() or dst.is_symlink():
        return
    if action == "copy":
        shutil.copy2(src, dst)
    elif action == "symlink":
        os.symlink(src, dst)
    else:
        raise ValueError(f"Unsupported batching.file_action: {action}")


def write_protein_batches(config: Dict[str, Any], clusters: Dict[str, List[str]], records: List[ProteinRecord], edge_weights) -> Path | None:
    batch_cfg = config.get("batching", {})
    if not batch_cfg.get("enabled", True) or not batch_cfg.get("output_to_protein_dir", True):
        return None

    paths = config.get("paths", {})
    root = Path(paths.get("aimd_root", ".")).resolve()
    protein_out_dir = resolve_path(paths.get("protein_out_dir"), root)
    assert protein_out_dir is not None
    ensure_dir(protein_out_dir)

    max_n = int(batch_cfg.get("max_proteins_per_batch", 500))
    action = batch_cfg.get("file_action", "copy")
    include_only_reps = bool(batch_cfg.get("include_only_representatives", False))
    rep_method = config.get("representative", {}).get("method", "highest_mean_similarity")
    record_map = {r.protein_id: r for r in records}
    cluster_to_batch = make_batches_by_cluster(clusters, max_n)

    rows = []
    for cid, members in sorted(clusters.items()):
        batch = cluster_to_batch[cid]
        rep = choose_representative(members, edge_weights, rep_method)
        selected_members = [rep] if include_only_reps else sorted(members)
        for pid in selected_members:
            r = record_map[pid]
            dst = protein_out_dir / batch / r.pdb_path.name
            copy_or_link(r.pdb_path, dst, action)
            rows.append({
                "protein_id": pid,
                "cluster_id": cid,
                "batch_id": batch,
                "source_pdb": str(r.pdb_path),
                "protein_path": str(dst),
                "is_representative": str(pid == rep).lower(),
                "file_action": action,
                "status": "success",
            })

    manifest_path = protein_out_dir / "protein_manifest.csv"
    write_csv(
        manifest_path,
        rows,
        ["protein_id", "cluster_id", "batch_id", "source_pdb", "protein_path", "is_representative", "file_action", "status"],
    )
    return manifest_path
