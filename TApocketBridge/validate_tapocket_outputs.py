#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path

from .utils import read_csv


def validate(pocket_manifest: str | Path, run_manifest: str | Path | None = None) -> dict[str, int]:
    pocket_manifest = Path(pocket_manifest).resolve()
    if not pocket_manifest.exists():
        raise FileNotFoundError(pocket_manifest)
    pocket_rows = read_csv(pocket_manifest)
    missing_files = 0
    checked_files = 0
    for row in pocket_rows:
        for key in ["pocket_pdb_path", "pocket_json_path", "box_tsv_path", "box_yaml_path", "residues_tsv_path", "summary_tsv_path"]:
            value = row.get(key, "")
            if not value:
                continue
            checked_files += 1
            if not Path(value).exists():
                missing_files += 1
                print(f"[missing] {key}: {value}")
    failed_runs = 0
    total_runs = 0
    if run_manifest:
        rp = Path(run_manifest).resolve()
        if rp.exists():
            runs = read_csv(rp)
            total_runs = len(runs)
            failed_runs = sum(1 for r in runs if r.get("status") != "success")
    result = {"pocket_rows": len(pocket_rows), "checked_files": checked_files, "missing_files": missing_files, "total_runs": total_runs, "failed_runs": failed_runs}
    print(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate AImd TApocketBridge outputs")
    parser.add_argument("--pocket-manifest", default="data/data_output/pocket/pocket_manifest.csv")
    parser.add_argument("--run-manifest", default="data/data_output/pocket/tapocket_run_manifest.csv")
    args = parser.parse_args()
    result = validate(args.pocket_manifest, args.run_manifest)
    if result["missing_files"] > 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
