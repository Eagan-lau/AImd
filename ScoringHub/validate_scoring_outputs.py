#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate AImd ClusterScore outputs")
    parser.add_argument("--out-dir", default="data/data_output/scoring/ClusterScore")
    args = parser.parse_args()
    out = Path(args.out_dir)
    required = [
        "clusterscore_results.xlsx",
        "best_affinity_long.csv",
        "protein_ligand_matrix.csv",
        "protein_binding_statistics.csv",
        "cluster_binding_statistics.csv",
        "clusterscore_report.json",
    ]
    missing = [str(out / name) for name in required if not (out / name).exists()]
    if missing:
        raise SystemExit("Missing ClusterScore outputs:\n" + "\n".join(missing))
    print(f"ClusterScore outputs look complete: {out}")


if __name__ == "__main__":
    main()
