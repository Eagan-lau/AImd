#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from .utils import ensure_dir, resolve_path, write_csv, write_json

AFF_RE = re.compile(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+")


def parse_vina_log(path: str | Path, num_affinities: int = 9) -> dict[str, Any]:
    path = Path(path)
    data: dict[str, Any] = {
        "affinities": [],
        "grid_size": "",
        "grid_space": "",
        "exhaustiveness": "",
        "random_seed": "",
    }
    if not path.exists():
        return data
    try:
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.startswith("Grid size"):
                parts = line.split(":", 1)[1].strip().split()
                data["grid_size"] = ",".join(parts[i] for i in range(1, len(parts), 2)) if len(parts) > 1 else ""
            elif line.startswith("Grid space"):
                data["grid_space"] = line.split(":", 1)[1].strip()
            elif line.startswith("Exhaustiveness"):
                data["exhaustiveness"] = line.split(":", 1)[1].strip()
            elif line.startswith("Performing docking") and "random seed" in line:
                try:
                    data["random_seed"] = line.split("(", 1)[1].split(")", 1)[0].replace("random seed:", "").strip()
                except Exception:
                    pass
            m = AFF_RE.match(line)
            if m and len(data["affinities"]) < num_affinities:
                data["affinities"].append(float(m.group(2)))
    except Exception as exc:
        data["parse_error"] = str(exc)
    return data


def summarize_results(config: dict[str, Any], run_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    out_dir = resolve_path(config.get("paths", {}).get("docking_out_dir", "data/data_output/docking_out"), root)
    assert out_dir is not None
    ensure_dir(out_dir)
    num_aff = int(config.get("postprocess", {}).get("num_affinities", 9) or 9)
    rows: list[dict[str, Any]] = []
    for row in run_rows:
        parsed = parse_vina_log(row.get("log_path", ""), num_aff)
        affinities = parsed.get("affinities", [])
        best = affinities[0] if affinities else ""
        out_pose = Path(row.get("out_pose_path", "")) if row.get("out_pose_path") else None
        rows.append({
            **row,
            "best_affinity": best,
            "affinities": ",".join(map(str, affinities)),
            "n_affinities": len(affinities),
            "grid_size": parsed.get("grid_size", ""),
            "grid_space": parsed.get("grid_space", ""),
            "exhaustiveness": parsed.get("exhaustiveness", ""),
            "random_seed": parsed.get("random_seed", ""),
            "pose_exists": str(bool(out_pose and out_pose.exists())).lower(),
        })
    write_csv(out_dir / "docking_result_manifest.csv", rows)

    success = sum(1 for r in rows if r.get("status") in {"success", "skipped"})
    failed = sum(1 for r in rows if r.get("status") == "failed")
    report = {
        "n_tasks": len(rows),
        "n_success_or_skipped": success,
        "n_failed": failed,
        "result_manifest": str(out_dir / "docking_result_manifest.csv"),
    }
    write_json(out_dir / "docking_report.json", report)
    return rows
