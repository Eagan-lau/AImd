#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Any

from .io import ProteinRecord
from .utils import write_tsv, write_csv, resolve_path


def choose_representative(members: List[str], edge_weights: Dict[Tuple[str, str], float], method: str = "highest_mean_similarity") -> str:
    members = sorted(members)
    if len(members) == 1:
        return members[0]
    if method == "first":
        return members[0]
    if method != "highest_mean_similarity":
        raise ValueError(f"Unsupported representative.method: {method}")

    best_pid = members[0]
    best_mean = -1.0
    for pid in members:
        vals = []
        for other in members:
            if other == pid:
                continue
            vals.append(edge_weights.get(tuple(sorted([pid, other])), 0.0))
        mean_val = sum(vals) / len(vals) if vals else 0.0
        if mean_val > best_mean:
            best_mean = mean_val
            best_pid = pid
    return best_pid


def summarize_cluster(members: List[str], edge_weights: Dict[Tuple[str, str], float]) -> Dict[str, float | int | str]:
    vals = []
    for i, a in enumerate(members):
        for b in members[i + 1:]:
            if a == b:
                continue
            w = edge_weights.get(tuple(sorted([a, b])))
            if w is not None:
                vals.append(w)
    if vals:
        return {
            "n_edges_in_cluster": len(vals),
            "mean_weight": sum(vals) / len(vals),
            "max_weight": max(vals),
            "min_weight": min(vals),
        }
    return {"n_edges_in_cluster": 0, "mean_weight": "", "max_weight": "", "min_weight": ""}


def write_cluster_outputs(config: Dict[str, Any], clusters: Dict[str, List[str]], records: List[ProteinRecord], edge_weights: Dict[Tuple[str, str], float]) -> Tuple[Path, Path, Path]:
    paths = config.get("paths", {})
    root = Path(paths.get("aimd_root", ".")).resolve()
    out_dir = resolve_path(paths.get("out_dir"), root)
    assert out_dir is not None

    rep_method = config.get("representative", {}).get("method", "highest_mean_similarity")
    record_map = {r.protein_id: r for r in records}

    clusters_tsv = out_dir / "clusters.tsv"
    reps_tsv = out_dir / "representatives.tsv"
    summary_csv = out_dir / "cluster_summary.csv"
    manifest_csv = out_dir / "structure_cluster_manifest.csv"

    cluster_rows = []
    rep_rows = []
    summary_rows = []
    manifest_rows = []

    for cid, members in sorted(clusters.items()):
        members = sorted(members)
        rep = choose_representative(members, edge_weights, rep_method)
        stats = summarize_cluster(members, edge_weights)
        for pid in members:
            r = record_map[pid]
            cluster_rows.append({"cluster_id": cid, "protein_id": pid})
            manifest_rows.append({
                "protein_id": pid,
                "cluster_id": cid,
                "is_representative": str(pid == rep).lower(),
                "source_pdb": str(r.pdb_path),
                "status": "success",
            })
        rep_rows.append({
            "cluster_id": cid,
            "representative_id": rep,
            "n_proteins": len(members),
            "selection_method": "singleton" if len(members) == 1 else rep_method,
        })
        summary_rows.append({
            "cluster_id": cid,
            "n_proteins": len(members),
            "representative_id": rep,
            **stats,
        })

    write_tsv(clusters_tsv, cluster_rows, ["cluster_id", "protein_id"], header=True)
    write_tsv(reps_tsv, rep_rows, ["cluster_id", "representative_id", "n_proteins", "selection_method"], header=True)
    write_csv(summary_csv, summary_rows, ["cluster_id", "n_proteins", "representative_id", "n_edges_in_cluster", "mean_weight", "max_weight", "min_weight"])
    write_csv(manifest_csv, manifest_rows, ["protein_id", "cluster_id", "is_representative", "source_pdb", "status"])

    return clusters_tsv, reps_tsv, summary_csv
