#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise RuntimeError("PyYAML is required. Install with: pip install PyYAML") from exc


def load_yaml(path: str | Path) -> dict[str, Any]:
    path = Path(path).resolve()
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def dump_yaml(path: str | Path, data: dict[str, Any]) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False, allow_unicode=True)


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_path(value: str | Path | None, root: str | Path) -> Path | None:
    if value in {None, ""}:
        return None
    p = Path(value).expanduser()
    return p.resolve() if p.is_absolute() else (Path(root).resolve() / p).resolve()


def read_csv(path: str | Path) -> list[dict[str, str]]:
    path = Path(path)
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: str | Path, rows: Iterable[dict[str, Any]], fieldnames: list[str] | None = None) -> Path:
    path = Path(path)
    ensure_dir(path.parent)
    rows = list(rows)
    if fieldnames is None:
        keys: list[str] = []
        for row in rows:
            for key in row.keys():
                if key not in keys:
                    keys.append(key)
        fieldnames = keys
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _stringify(row.get(k, "")) for k in fieldnames})
    return path


def write_json(path: str | Path, data: Any) -> Path:
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False)
    return path


def _stringify(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (list, tuple)):
        return ",".join(map(str, value))
    return str(value)


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on", "success", "ok"}


def safe_id(value: str) -> str:
    value = str(value).strip()
    value = re.sub(r"[\\/\s:;|]+", "_", value)
    value = re.sub(r"[^A-Za-z0-9_.@+=-]+", "_", value)
    return value.strip("_") or "NA"


def copy_or_link(src: Path, dst: Path, action: str = "copy", overwrite: bool = False) -> Path:
    ensure_dir(dst.parent)
    if dst.exists() or dst.is_symlink():
        if overwrite:
            if dst.is_dir() and not dst.is_symlink():
                shutil.rmtree(dst)
            else:
                dst.unlink()
        else:
            return dst
    if action == "copy":
        shutil.copy2(src, dst)
    elif action == "symlink":
        os.symlink(src, dst)
    elif action == "hardlink":
        os.link(src, dst)
    else:
        raise ValueError(f"Unsupported file action: {action}")
    return dst


def run_command(command: str, *, timeout: int | None = None, cwd: str | Path | None = None, log_path: str | Path | None = None) -> int:
    """Run a shell command and optionally tee stdout/stderr to a file."""
    if log_path:
        log_path = Path(log_path)
        ensure_dir(log_path.parent)
        with log_path.open("w", encoding="utf-8") as log:
            process = subprocess.Popen(command, shell=True, cwd=str(cwd) if cwd else None, stdout=log, stderr=subprocess.STDOUT, text=True)
            try:
                process.communicate(timeout=timeout)
                return int(process.returncode or 0)
            except subprocess.TimeoutExpired:
                process.terminate()
                log.write(f"\n[AImd DockingHub] TIMEOUT after {timeout}s: {command}\n")
                return 124
    process = subprocess.Popen(command, shell=True, cwd=str(cwd) if cwd else None)
    try:
        process.communicate(timeout=timeout)
        return int(process.returncode or 0)
    except subprocess.TimeoutExpired:
        process.terminate()
        return 124


def which(executable: str) -> str | None:
    return shutil.which(executable)


def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def quote(path: str | Path) -> str:
    return "'" + str(path).replace("'", "'\\''") + "'"


def is_probably_pdbqt(path: str | Path) -> bool:
    return Path(path).suffix.lower() == ".pdbqt"


def numeric(value: Any, default: float | None = None) -> float | None:
    try:
        if value in {None, ""}:
            return default
        return float(value)
    except Exception:
        return default

def resolve_external_tool(tool_name: str, configured_executable: str | None, root: str | Path, config: dict[str, Any] | None = None) -> str:
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
