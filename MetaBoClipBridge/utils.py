#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
import json
import os
import re
import shutil
from copy import deepcopy
from pathlib import Path
from typing import Any, Iterable

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise RuntimeError("PyYAML is required. Install with: pip install PyYAML") from exc


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def resolve_path(value: str | Path | None, root: str | Path) -> Path | None:
    if value in {None, ""}:
        return None
    p = Path(value).expanduser()
    return p.resolve() if p.is_absolute() else (Path(root).resolve() / p).resolve()


def load_yaml(path: str | Path) -> dict[str, Any]:
    p = Path(path).resolve()
    with p.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def dump_yaml(path: str | Path, data: dict[str, Any]) -> Path:
    p = Path(path)
    ensure_dir(p.parent)
    with p.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False, allow_unicode=True)
    return p


def write_json(path: str | Path, data: Any) -> Path:
    p = Path(path)
    ensure_dir(p.parent)
    with p.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False)
    return p


def read_csv(path: str | Path) -> list[dict[str, str]]:
    p = Path(path)
    if not p.exists():
        return []
    with p.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: str | Path, rows: Iterable[dict[str, Any]], fieldnames: list[str] | None = None) -> Path:
    p = Path(path)
    ensure_dir(p.parent)
    rows = list(rows)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row.keys():
                if key not in fieldnames:
                    fieldnames.append(key)
    with p.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _stringify(row.get(k, "")) for k in fieldnames})
    return p


def _stringify(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (list, tuple, set)):
        return ",".join(map(str, value))
    return str(value)


def safe_id(value: Any) -> str:
    value = str(value).strip()
    value = re.sub(r"[\\/\s:;|]+", "_", value)
    value = re.sub(r"[^A-Za-z0-9_.@+=-]+", "_", value)
    return value.strip("_") or "NA"


def copy_or_link(src: str | Path, dst: str | Path, action: str = "copy", overwrite: bool = False) -> Path:
    src = Path(src)
    dst = Path(dst)
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


def numeric(value: Any, default: float | None = None) -> float | None:
    try:
        if value in {None, "", "NA", "NaN", "nan"}:
            return default
        return float(value)
    except Exception:
        return default


def deep_merge(base: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy(base or {})
    for key, val in (overlay or {}).items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], val)
        else:
            out[key] = deepcopy(val)
    return out
