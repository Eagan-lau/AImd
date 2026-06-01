#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path

from .utils import load_yaml, ensure_dir, resolve_path
from .io import load_protein_records, write_input_manifest
from .foldseek import run_foldseek
from .graph import build_filtered_graph
from .hipmcl import run_hipmcl, parse_hipmcl_output, create_singleton_clusters
from .clustering import write_cluster_outputs
from .batching import write_protein_batches


def run_rgpc(config_path: str | Path) -> None:
    config_path = Path(config_path).resolve()
    config = load_yaml(config_path)
    paths = config.get("paths", {})
    root = Path(paths.get("aimd_root", ".")).resolve()

    out_dir = resolve_path(paths.get("out_dir"), root)
    assert out_dir is not None
    ensure_dir(out_dir)

    print("[RGPC] Loading protein structure records...")
    records = load_protein_records(config)
    print(f"[RGPC] Loaded {len(records)} protein structures.")
    write_input_manifest(records, out_dir / "rgpc_input_manifest.csv")

    print("[RGPC] Running Foldseek...")
    raw_tsv = run_foldseek(config, records)

    print("[RGPC] Filtering Foldseek results and building graph...")
    edges_tsv, edges_abc, edge_weights = build_filtered_graph(config, raw_tsv, records)
    print(f"[RGPC] Graph edges: {len(edge_weights)}")

    if edge_weights:
        print("[RGPC] Running HipMCL...")
        hip_output = run_hipmcl(config, edges_abc)
        print("[RGPC] Parsing HipMCL clusters...")
        clusters = parse_hipmcl_output(hip_output, records)
    else:
        hip_output = out_dir / "hipmcl" / "hipmcl_output.txt"
        ensure_dir(hip_output.parent)
        hip_output.write_text("", encoding="utf-8")
        print("[RGPC] No graph edges found. Creating singleton clusters and skipping HipMCL.")
        clusters = create_singleton_clusters(records)
    print(f"[RGPC] Clusters: {len(clusters)}")

    print("[RGPC] Writing cluster outputs...")
    clusters_tsv, reps_tsv, summary_csv = write_cluster_outputs(config, clusters, records, edge_weights)

    print("[RGPC] Writing downstream protein batches...")
    protein_manifest = write_protein_batches(config, clusters, records, edge_weights)

    print("[RGPC] Done.")
    print(f"[RGPC] edges_tsv: {edges_tsv}")
    print(f"[RGPC] hipmcl_output: {hip_output}")
    print(f"[RGPC] clusters_tsv: {clusters_tsv}")
    print(f"[RGPC] representatives: {reps_tsv}")
    print(f"[RGPC] cluster_summary: {summary_csv}")
    if protein_manifest:
        print(f"[RGPC] protein_manifest: {protein_manifest}")


def main() -> None:
    parser = argparse.ArgumentParser(description="AImd/RGPC: Representative Graph-based Protein Clustering")
    parser.add_argument("--config", required=True, help="Path to configs/RGPC/rgpc.yaml")
    args = parser.parse_args()
    run_rgpc(args.config)


if __name__ == "__main__":
    main()
