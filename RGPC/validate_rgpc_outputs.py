#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path
import csv


def count_rows(path: Path, delimiter: str = "\t", header: bool = True) -> int:
    if not path.exists():
        return 0
    with open(path, "r", encoding="utf-8") as f:
        n = sum(1 for line in f if line.strip())
    return max(0, n - 1) if header else n


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate RGPC output files")
    parser.add_argument("--rgpc-dir", default="data/cluster/RGPC")
    parser.add_argument("--protein-manifest", default="data/protein/protein_manifest.csv")
    args = parser.parse_args()

    rgpc_dir = Path(args.rgpc_dir)
    protein_manifest = Path(args.protein_manifest)

    files = {
        "edges": rgpc_dir / "graph" / "structure_edges.tsv",
        "clusters": rgpc_dir / "clusters.tsv",
        "representatives": rgpc_dir / "representatives.tsv",
        "cluster_summary": rgpc_dir / "cluster_summary.csv",
        "protein_manifest": protein_manifest,
    }

    for name, path in files.items():
        print(f"{name}: {'OK' if path.exists() else 'MISSING'} -> {path}")

    print("n_edges:", count_rows(files["edges"], delimiter="\t", header=True))
    print("n_cluster_memberships:", count_rows(files["clusters"], delimiter="\t", header=True))
    print("n_representatives:", count_rows(files["representatives"], delimiter="\t", header=True))


if __name__ == "__main__":
    main()
