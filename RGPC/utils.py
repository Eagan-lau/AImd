#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
import os
import shlex
import subprocess
from pathlib import Path
from typing import Iterable, List, Dict, Any, Optional

import yaml


def load_yaml(path: str | Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if data is None:
        return {}
    return data


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def resolve_path(path: str | Path | None, root: str | Path) -> Optional[Path]:
    if path is None or str(path).strip() == "":
        return None
    p = Path(path)
    if p.is_absolute():
        return p
    return Path(root) / p


def run_command(cmd: List[str] | str, log_path: str | Path | None = None, dry_run: bool = False) -> None:
    if isinstance(cmd, list):
        cmd_display = " ".join(shlex.quote(str(x)) for x in cmd)
    else:
        cmd_display = cmd

    print(f"[RGPC] Running: {cmd_display}")
    if dry_run:
        return

    if log_path is not None:
        ensure_dir(Path(log_path).parent)
        with open(log_path, "a", encoding="utf-8") as log:
            log.write(f"\n[COMMAND] {cmd_display}\n")
            proc = subprocess.run(cmd, shell=isinstance(cmd, str), stdout=log, stderr=subprocess.STDOUT, text=True)
    else:
        proc = subprocess.run(cmd, shell=isinstance(cmd, str))

    if proc.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {proc.returncode}: {cmd_display}")


def write_tsv(path: str | Path, rows: Iterable[Dict[str, Any]], fieldnames: List[str], header: bool = True) -> None:
    ensure_dir(Path(path).parent)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        if header:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_csv(path: str | Path, rows: Iterable[Dict[str, Any]], fieldnames: List[str]) -> None:
    ensure_dir(Path(path).parent)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_delimited(path: str | Path, delimiter: str = "\t", has_header: bool = True, fieldnames: List[str] | None = None) -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8") as f:
        if has_header:
            reader = csv.DictReader(f, delimiter=delimiter)
        else:
            reader = csv.DictReader(f, fieldnames=fieldnames, delimiter=delimiter)
        return [dict(r) for r in reader]


def shell_format(template: str, **kwargs: Any) -> str:
    return template.format(**{k: shlex.quote(str(v)) for k, v in kwargs.items()})

def resolve_external_tool(tool_name: str, configured_executable: str | None, root: str | Path, config: Dict[str, Any] | None = None) -> str:
    """Resolve a third-party executable through AImd/third_party/tools.yaml."""
    tools_config = None
    if config:
        tools_config = config.get("third_party", {}).get("tools_config")
    try:
        from third_party.tool_manager import resolve_tool
        return resolve_tool(tool_name, configured_executable, root, tools_config)
    except Exception:
        if configured_executable and str(configured_executable).strip().lower() not in {"", "auto", "null", "none"}:
            return str(configured_executable)
        return tool_name
