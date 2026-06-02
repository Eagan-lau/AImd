#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from .utils import ensure_dir, quote, which, write_csv, resolve_external_tool


def _run_one_vina(task: dict[str, Any], vina_cfg: dict[str, Any], overwrite: bool) -> dict[str, Any]:
    out_pose = Path(task["out_pose_path"])
    log_path = Path(task["log_path"])
    ensure_dir(out_pose.parent)
    if out_pose.exists() and log_path.exists() and not overwrite:
        return {**task, "status": "skipped", "return_code": "0", "message": "existing output"}
    root = Path(vina_cfg.get("aimd_root", ".")).resolve()
    executable = resolve_external_tool("vina", vina_cfg.get("executable", "auto"), root, {"third_party": vina_cfg.get("third_party", {})})
    if which(str(executable)) is None and Path(str(executable)).name == str(executable):
        return {**task, "status": "failed", "return_code": "127", "message": f"Vina executable not found: {executable}"}
    exhaustiveness = int(vina_cfg.get("exhaustiveness", 8) or 8)
    cpu = int(vina_cfg.get("cpu_per_job", 1) or 1)
    timeout = int(vina_cfg.get("timeout", 3600) or 3600)
    additional_params = str(vina_cfg.get("additional_params", "")).strip()
    command = f"{executable} --config {quote(task['config_path'])} --exhaustiveness {exhaustiveness} --cpu {cpu} {additional_params}".strip()
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"[AImd DockingHub] command: {command}\n")
        process = subprocess.Popen(command, shell=True, stdout=log, stderr=subprocess.STDOUT, text=True)
        try:
            process.communicate(timeout=timeout)
            rc = int(process.returncode or 0)
        except subprocess.TimeoutExpired:
            process.terminate()
            log.write(f"\n[AImd DockingHub] TIMEOUT after {timeout}s\n")
            rc = 124
    status = "success" if rc == 0 and out_pose.exists() else "failed"
    return {**task, "status": status, "return_code": str(rc), "message": "vina completed" if status == "success" else "vina failed"}


def run_vina_tasks(config: dict[str, Any], tasks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    run_cfg = dict(config.get("docking", {}).get("vina", {}))
    run_cfg["third_party"] = config.get("third_party", {})
    run_cfg["aimd_root"] = config.get("paths", {}).get("aimd_root", ".")
    execution = config.get("execution", {})
    overwrite = bool(config.get("output", {}).get("overwrite", False))
    if not bool(config.get("docking", {}).get("run", True)):
        rows = [{**task, "status": "not_run", "return_code": "", "message": "docking.run=false"} for task in tasks]
        return rows
    total_threads = int(execution.get("threads", min(64, os.cpu_count() or 1)) or 1)
    cpu = int(run_cfg.get("cpu_per_job", 1) or 1)
    workers = max(1, total_threads // max(1, cpu))
    workers = int(execution.get("workers", workers) or workers)
    rows: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_run_one_vina, task, run_cfg, overwrite) for task in tasks]
        for future in as_completed(futures):
            try:
                rows.append(future.result())
            except Exception as exc:
                rows.append({"job_id": "unknown", "status": "failed", "return_code": "", "message": str(exc)})
    task_dir = Path(config.get("paths", {}).get("aimd_root", ".")).resolve() / config.get("paths", {}).get("docking_task_dir", "data/data_output/docking_tasks")
    write_csv(task_dir / "docking_run_manifest.csv", rows)
    return rows
