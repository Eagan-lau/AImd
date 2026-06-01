#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path
from tool_manager import check_tool, load_tools_config

def main() -> None:
    parser = argparse.ArgumentParser(description="Check third-party tools configured for AImd.")
    parser.add_argument("--config", default="third_party/tools.yaml", help="Path to third_party/tools.yaml")
    parser.add_argument("--root", default=".", help="AImd root directory")
    parser.add_argument("--out", default="third_party/tool_check_report.csv", help="Output CSV report")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    cfg_path = Path(args.config)
    if not cfg_path.is_absolute():
        cfg_path = root / cfg_path
    cfg = load_tools_config(cfg_path, root)
    tool_names = sorted((cfg.get("tools") or {}).keys())
    rows = []
    failed_required = []
    for name in tool_names:
        row = check_tool(name, root, cfg_path)
        rows.append(row)
        prefix = "OK" if row["status"] == "ok" else ("SKIP" if row["status"] == "disabled" else "FAIL")
        print(f"[{prefix}] {name}: {row['executable']} | {row['status']} | {row['message']}")
        if row.get("required") and row.get("enabled") and row.get("status") != "ok":
            failed_required.append(name)
    out = Path(args.out)
    if not out.is_absolute():
        out = root / out
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["tool", "enabled", "required", "executable", "status", "return_code", "message"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"[AImd third_party] Wrote report: {out}")
    if failed_required:
        print("[AImd third_party] Missing or failed required tools: " + ", ".join(failed_required))
        sys.exit(2)
if __name__ == "__main__":
    main()
