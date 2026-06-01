#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate AImd MetaBoClipBridge outputs")
    ap.add_argument("--out-dir", default="data/metaboclip/results")
    args = ap.parse_args()
    out = Path(args.out_dir)
    required = [
        out / "metaboclip_run_manifest.csv",
        out / "metaboclip_final_ranking.csv",
        out / "metaboclip_report.json",
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing MetaBoClipBridge outputs:\n" + "\n".join(missing))
    print(f"MetaBoClipBridge outputs look complete: {out}")

if __name__ == "__main__":
    main()
