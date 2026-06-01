#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Any, List, Tuple

from .io import ProteinRecord
from .utils import ensure_dir, resolve_path


RAW_FIELDS_5 = ["query", "target", "qtmscore", "ttmscore", "alntmscore"]
RAW_FIELDS_4 = ["query", "target", "qtmscore", "ttmscore"]


def _to_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def compute_weight(qtmscore: float, ttmscore: float, alntmscore: float, method: str) -> float:
    if method == "min_qtmscore_ttmscore":
        return min(qtmscore, ttmscore)
    if method == "mean_qtmscore_ttmscore":
        return (qtmscore + ttmscore) / 2.0
    if method == "geometric_mean":
        return math.sqrt(max(qtmscore, 0.0) * max(ttmscore, 0.0))
    if method == "alntmscore":
        return alntmscore
    raise ValueError(f"Unsupported weight_method: {method}")


def merge_values(values: List[float], method: str) -> float:
    if not values:
        return 0.0
    if method == "min":
        return min(values)
    if method == "max":
        return max(values)
    if method == "mean":
        return sum(values) / len(values)
    raise ValueError(f"Unsupported duplicate_edge_method: {method}")


def parse_foldseek_tsv(raw_tsv: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(raw_tsv, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                parts = line.split()
            if len(parts) >= 5:
                rows.append(dict(zip(RAW_FIELDS_5, parts[:5])))
            elif len(parts) == 4:
                row = dict(zip(RAW_FIELDS_4, parts[:4]))
                row["alntmscore"] = "nan"
                rows.append(row)
    return rows


def build_filtered_graph(config: Dict[str, Any], raw_tsv: Path, records: List[ProteinRecord]) -> Tuple[Path, Path, Dict[Tuple[str, str], float]]:
    paths = config.get("paths", {})
    root = Path(paths.get("aimd_root", ".")).resolve()
    out_dir = resolve_path(paths.get("out_dir"), root)
    assert out_dir is not None

    filter_cfg = config.get("filter", {})
    graph_cfg = config.get("graph", {})
    resume = config.get("resume", {})
    skip_existing = bool(resume.get("enabled", True) and resume.get("skip_existing", True))

    graph_dir = out_dir / "graph"
    ensure_dir(graph_dir)
    edges_tsv = graph_dir / "structure_edges.tsv"
    edges_abc = graph_dir / "structure_edges.abc"

    if skip_existing and edges_tsv.exists() and edges_abc.exists():
        edge_weights = load_edge_weights(edges_tsv)
        print(f"[RGPC] Skip graph construction: {edges_tsv}")
        return edges_tsv, edges_abc, edge_weights

    remove_self = bool(filter_cfg.get("remove_self_hits", True))
    min_q = float(filter_cfg.get("min_qtmscore", 0.5))
    min_t = float(filter_cfg.get("min_ttmscore", 0.5))
    weight_method = filter_cfg.get("weight_method", "min_qtmscore_ttmscore")
    make_symmetric = bool(graph_cfg.get("make_symmetric", True))
    duplicate_method = graph_cfg.get("duplicate_edge_method", "min")
    write_header = bool(graph_cfg.get("write_header_tsv", True))

    edge_values: Dict[Tuple[str, str], List[float]] = {}
    edge_detail: Dict[Tuple[str, str], Dict[str, float | str]] = {}

    raw_rows = parse_foldseek_tsv(raw_tsv)
    for row in raw_rows:
        q = row["query"]
        t = row["target"]
        if remove_self and q == t:
            continue
        qtm = _to_float(row["qtmscore"])
        ttm = _to_float(row["ttmscore"])
        atm = _to_float(row["alntmscore"])
        if math.isnan(qtm) or math.isnan(ttm):
            continue
        if qtm < min_q or ttm < min_t:
            continue
        weight = compute_weight(qtm, ttm, atm, weight_method)
        key = tuple(sorted([q, t])) if make_symmetric else (q, t)
        edge_values.setdefault(key, []).append(weight)
        prev = edge_detail.get(key)
        if prev is None or weight < float(prev["weight"]):
            edge_detail[key] = {
                "query": key[0],
                "target": key[1],
                "qtmscore": qtm,
                "ttmscore": ttm,
                "alntmscore": atm,
                "weight": weight,
            }

    merged: Dict[Tuple[str, str], float] = {k: merge_values(v, duplicate_method) for k, v in edge_values.items()}

    with open(edges_tsv, "w", encoding="utf-8") as f:
        if write_header:
            f.write("query\ttarget\tweight\n")
        for (q, t), w in sorted(merged.items()):
            f.write(f"{q}\t{t}\t{w:.6f}\n")

    with open(edges_abc, "w", encoding="utf-8") as f:
        for (q, t), w in sorted(merged.items()):
            f.write(f"{q}\t{t}\t{w:.6f}\n")

    summary_path = graph_dir / "structure_similarity_summary.csv"
    proteins = {r.protein_id for r in records}
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("n_proteins,n_edges,min_qtmscore,min_ttmscore,weight_method,make_symmetric,duplicate_edge_method\n")
        f.write(f"{len(proteins)},{len(merged)},{min_q},{min_t},{weight_method},{make_symmetric},{duplicate_method}\n")

    return edges_tsv, edges_abc, merged


def load_edge_weights(edges_tsv: Path) -> Dict[Tuple[str, str], float]:
    weights: Dict[Tuple[str, str], float] = {}
    with open(edges_tsv, "r", encoding="utf-8") as f:
        first = True
        for line in f:
            line = line.strip()
            if not line:
                continue
            if first and line.startswith("query\t"):
                first = False
                continue
            first = False
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            q, t, w = parts[:3]
            weights[tuple(sorted([q, t]))] = float(w)
    return weights
