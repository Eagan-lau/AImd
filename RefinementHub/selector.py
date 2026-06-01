#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path
from typing import Any

from .utils import as_bool, numeric, read_csv, resolve_path, write_csv


def _cluster_id(row: dict[str, str]) -> str:
    return str(row.get("cluster") or row.get("cluster_id") or row.get("family") or "").strip()


def load_top_clusters(config: dict[str, Any]) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    sel = config.get("selection", {})
    top_csv = resolve_path(config.get("paths", {}).get("clusterscore_top_csv", "data/scoring/ClusterScore/top10_clusters.csv"), root)
    all_csv = resolve_path(config.get("paths", {}).get("clusterscore_all_csv", "data/scoring/ClusterScore/cluster_binding_statistics.csv"), root)
    rows = read_csv(top_csv) if top_csv and top_csv.exists() else []
    if not rows:
        rows = read_csv(all_csv) if all_csv and all_csv.exists() else []
    if not rows:
        raise FileNotFoundError(f"No ClusterScore table found: {top_csv} or {all_csv}")

    explicit = [str(x) for x in sel.get("include_clusters", []) or []]
    if explicit:
        keep = {str(x) for x in explicit}
        rows = [r for r in rows if _cluster_id(r) in keep]
    else:
        sort_by = sel.get("sort_by", "rank")
        if sort_by in rows[0]:
            if str(sort_by).lower() == "rank":
                rows.sort(key=lambda r: numeric(r.get(sort_by), 10**12) or 10**12)
            else:
                asc = bool(sel.get("sort_ascending", False))
                rows.sort(key=lambda r: numeric(r.get(sort_by), 0.0) or 0.0, reverse=not asc)
        top_n = int(sel.get("top_n_clusters", 10) or 10)
        rows = rows[:top_n]
    out = []
    for i, row in enumerate(rows, start=1):
        cid = _cluster_id(row)
        if not cid:
            continue
        out.append({"selected_rank": i, "cluster_id": cid, **row})
    if not out:
        raise RuntimeError("ClusterScore table was found, but no cluster_id/cluster values could be selected")
    return out


def select_proteins_for_refinement(config: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    source_manifest = resolve_path(config.get("paths", {}).get("source_protein_manifest", "data/protein/protein_manifest.csv"), root)
    if source_manifest is None or not source_manifest.exists():
        raise FileNotFoundError(f"Protein manifest not found: {source_manifest}")
    protein_rows = read_csv(source_manifest)
    if not protein_rows:
        raise RuntimeError(f"Protein manifest is empty: {source_manifest}")
    cluster_rows = load_top_clusters(config)
    selected_clusters = [str(r["cluster_id"]) for r in cluster_rows]
    selected_set = set(selected_clusters)

    sel = config.get("selection", {})
    mode = str(sel.get("protein_mode", "all"))
    max_per_cluster = sel.get("max_proteins_per_cluster")
    max_per_cluster = int(max_per_cluster) if max_per_cluster not in {None, "", 0, "0"} else None
    keep_status = {"", "success", "ok", "true", "1"}
    by_cluster_count: dict[str, int] = {}
    selected: list[dict[str, Any]] = []
    cluster_rank = {str(r["cluster_id"]): str(r.get("selected_rank", "")) for r in cluster_rows}
    for row in protein_rows:
        status = str(row.get("status", "success")).strip().lower()
        if status not in keep_status:
            continue
        cid = str(row.get("cluster_id", "")).strip()
        if cid not in selected_set:
            continue
        if mode == "representatives_only" and not as_bool(row.get("is_representative", False)):
            continue
        if max_per_cluster is not None and by_cluster_count.get(cid, 0) >= max_per_cluster:
            continue
        new = dict(row)
        new["selected_for_refinement"] = "true"
        new["cluster_rank"] = cluster_rank.get(cid, "")
        selected.append(new)
        by_cluster_count[cid] = by_cluster_count.get(cid, 0) + 1
    if not selected:
        raise RuntimeError("No proteins were selected for refined docking from the selected ClusterScore clusters")
    return cluster_rows, selected


def write_refinement_selection(config: dict[str, Any]) -> tuple[Path, Path]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    out_manifest = resolve_path(config.get("paths", {}).get("selected_protein_manifest", "data/refinement/selected_protein_manifest.csv"), root)
    out_clusters = resolve_path(config.get("paths", {}).get("selected_clusters_csv", "data/refinement/selected_clusters.csv"), root)
    assert out_manifest is not None and out_clusters is not None
    cluster_rows, protein_rows = select_proteins_for_refinement(config)
    write_csv(out_clusters, cluster_rows)
    write_csv(out_manifest, protein_rows)
    return out_manifest, out_clusters
