#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate AImd RefinementHub outputs")
    ap.add_argument("--root", default=".")
    args = ap.parse_args()
    root = Path(args.root)
    required = [
        root / "data/data_output/refinement/selected_clusters.csv",
        root / "data/data_output/refinement/selected_protein_manifest.csv",
        root / "data/data_output/refined/docking_out/docking_result_manifest.csv",
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing RefinementHub outputs:\n" + "\n".join(missing))
    print("RefinementHub outputs look complete.")

if __name__ == "__main__":
    main()
