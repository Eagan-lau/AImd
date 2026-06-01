#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path

from .utils import read_csv


def validate(root: str | Path = ".") -> int:
    root = Path(root).resolve()
    required = [
        root / "data/ensemble/conformer_manifest.csv",
        root / "data/cofactor_mapped/cofactor_manifest.csv",
        root / "data/receptor/receptor_manifest.csv",
        root / "data/docking_tasks/docking_task_manifest.csv",
        root / "data/docking_out/docking_result_manifest.csv",
    ]
    errors = []
    for fp in required:
        if not fp.exists():
            errors.append(f"missing: {fp}")
    if errors:
        print("[DockingHub validate] FAILED")
        for e in errors:
            print(" -", e)
        return 1
    results = read_csv(root / "data/docking_out/docking_result_manifest.csv")
    print(f"[DockingHub validate] OK: {len(results)} docking result rows")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate DockingHub outputs")
    parser.add_argument("--root", default=".", help="AImd root")
    args = parser.parse_args()
    raise SystemExit(validate(args.root))


if __name__ == "__main__":
    main()
