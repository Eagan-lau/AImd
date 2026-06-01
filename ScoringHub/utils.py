#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import ast
import csv
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise RuntimeError("PyYAML is required. Install with: pip install PyYAML") from exc


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_path(value: str | Path | None, root: str | Path) -> Path | None:
    if value in {None, ""}:
        return None
    p = Path(value).expanduser()
    return p.resolve() if p.is_absolute() else (Path(root).resolve() / p).resolve()


def load_yaml(path: str | Path) -> dict[str, Any]:
    path = Path(path).resolve()
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def write_json(path: str | Path, data: Any) -> Path:
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False)
    return path


def read_csv(path: str | Path) -> list[dict[str, str]]:
    path = Path(path)
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: str | Path, rows: Iterable[dict[str, Any]], fieldnames: list[str] | None = None) -> Path:
    rows = list(rows)
    path = Path(path)
    ensure_dir(path.parent)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row.keys():
                if key not in fieldnames:
                    fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _stringify(row.get(k, "")) for k in fieldnames})
    return path


def _stringify(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (list, tuple)):
        return ",".join(map(str, value))
    return str(value)


def numeric(value: Any, default: float | None = None) -> float | None:
    try:
        if value in {None, "", "NA", "NaN", "nan"}:
            return default
        if isinstance(value, float) and math.isnan(value):
            return default
        return float(value)
    except Exception:
        return default


AFFINITY_RE = re.compile(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+")


def parse_vina_log(path: str | Path, num_affinities: int = 9) -> dict[str, Any]:
    """Parse AutoDock Vina stdout/log and return affinities plus run parameters."""
    path = Path(path)
    data: dict[str, Any] = {
        "affinities": [],
        "grid_size": "",
        "grid_space": "",
        "exhaustiveness": "",
        "random_seed": "",
    }
    if not path.exists():
        data["parse_error"] = "log_not_found"
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
            match = AFFINITY_RE.match(line)
            if match and len(data["affinities"]) < num_affinities:
                data["affinities"].append(float(match.group(2)))
    except Exception as exc:
        data["parse_error"] = str(exc)
    return data


_ALLOWED_AST_NODES = {
    ast.Expression,
    ast.BinOp,
    ast.UnaryOp,
    ast.Num,
    ast.Constant,
    ast.Name,
    ast.Load,
    ast.Add,
    ast.Sub,
    ast.Mult,
    ast.Div,
    ast.Pow,
    ast.Mod,
    ast.USub,
    ast.UAdd,
    ast.Call,
    ast.Compare,
    ast.BoolOp,
    ast.And,
    ast.Or,
    ast.Eq,
    ast.NotEq,
    ast.Lt,
    ast.LtE,
    ast.Gt,
    ast.GtE,
    ast.IfExp,
}

_ALLOWED_FUNCS = {
    "min": min,
    "max": max,
    "abs": abs,
    "round": round,
    "sqrt": math.sqrt,
    "log": math.log,
    "log10": math.log10,
    "exp": math.exp,
}


def safe_eval_formula(expression: str, variables: dict[str, Any]) -> float:
    """Evaluate a simple numeric formula from YAML without exposing builtins."""
    tree = ast.parse(expression, mode="eval")
    for node in ast.walk(tree):
        if type(node) not in _ALLOWED_AST_NODES:
            raise ValueError(f"Unsupported formula syntax: {type(node).__name__}")
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name) or node.func.id not in _ALLOWED_FUNCS:
                raise ValueError("Only whitelisted numeric functions are allowed in formulas")
    names = {**_ALLOWED_FUNCS, **variables}
    return float(eval(compile(tree, "<ClusterScoreFormula>", "eval"), {"__builtins__": {}}, names))
