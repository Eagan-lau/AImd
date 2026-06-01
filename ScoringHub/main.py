#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path

from .clusterscore import run_clusterscore_core
from .utils import load_yaml, write_json


def run_clusterscore(config_path: str | Path) -> Path:
    config_path = Path(config_path).resolve()
    config = load_yaml(config_path)
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    print(f"[ClusterScore] AImd root: {root}")
    excel_path = run_clusterscore_core(config)
    out_dir = excel_path.parent
    write_json(out_dir / "run_config_snapshot.json", config)
    print(f"[ClusterScore] finished: {excel_path}")
    return excel_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run AImd ClusterScore: best-affinity extraction and cluster-level prioritization")
    parser.add_argument("--config", required=True, help="Path to configs/Scoring/cluster_score.yaml")
    args = parser.parse_args()
    run_clusterscore(args.config)


if __name__ == "__main__":
    main()
