#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate the AImd project layout and basic configuration consistency.

This script does not run external tools. It checks:
1. required module directories;
2. required config files;
3. the active unified MetaboClip backend location;
4. YAML parseability;
5. Python importability of AImd wrapper modules;
6. canonical data_input/data_output layout;
7. availability of common external executables when present in PATH.
"""

from __future__ import annotations

import argparse
import importlib
import shlex
import shutil
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


REQUIRED_DIRS = [
    "MolLink", "RGPC", "TApocket", "TApocketBridge", "DockingHub",
    "ScoringHub", "RefinementHub", "MetaBoClipBridge",
    "orchestrator", "configs", "data", "docs", "third_party",
    "metaboclip_unified",
]

REQUIRED_DATA_DIRS = [
    "data/data_input",
    "data/data_input/protein",
    "data/data_input/ligand",
    "data/data_input/cofactor",
    "data/data_input/workflow",
    "data/data_output",
]

REQUIRED_CONFIGS = [
    "configs/RGPC/rgpc.yaml",
    "configs/MolLink/mollink.yaml",
    "configs/TApocket/tapocket_batch.yaml",
    "configs/Docking/docking.yaml",
    "configs/Scoring/cluster_score.yaml",
    "configs/Refinement/refine_from_clusterscore.yaml",
    "configs/MetaBoClip/metaboclip_bridge.yaml",
    "configs/workflows/full_iterative_metaboclip.yaml",
]

IMPORT_MODULES = [
    "MolLink.main",
    "RGPC.main",
    "TApocketBridge.main",
    "DockingHub.main",
    "ScoringHub.main",
    "RefinementHub.main",
    "MetaBoClipBridge.main",
    "orchestrator.full_iterative_metaboclip",
]

EXTERNAL_TOOLS = [
    "foldseek", "hipmcl", "prepare_ligand", "prepare_receptor", "vina",
]

OPTIONAL_EXTERNAL_TOOLS = [
    "pymol",
]

YAML_CHECK_DIRS = [
    "configs",
]


def check_layout(root: Path) -> list[str]:
    errors: list[str] = []
    for rel in REQUIRED_DIRS:
        if not (root / rel).is_dir():
            errors.append(f"missing directory: {rel}")
    for rel in REQUIRED_DATA_DIRS:
        if not (root / rel).is_dir():
            errors.append(f"missing data layout directory: {rel}")
    for rel in REQUIRED_CONFIGS:
        if not (root / rel).is_file():
            errors.append(f"missing config: {rel}")
    return errors


def check_yaml(root: Path) -> list[str]:
    if yaml is None:
        return ["PyYAML is not installed; cannot parse YAML configs"]
    errors: list[str] = []
    for rel_dir in YAML_CHECK_DIRS:
        base = root / rel_dir
        if not base.exists():
            continue
        for path in base.rglob("*.yaml"):
            try:
                with path.open("r", encoding="utf-8") as handle:
                    yaml.safe_load(handle)
            except Exception as exc:
                errors.append(f"YAML parse error: {path.relative_to(root)} :: {exc}")
    return errors


def check_imports(root: Path) -> list[str]:
    errors: list[str] = []
    sys.path.insert(0, str(root))
    unified_root = root / "metaboclip_unified"
    if str(unified_root) not in sys.path:
        sys.path.insert(0, str(unified_root))
    for mod in IMPORT_MODULES:
        try:
            importlib.import_module(mod)
        except Exception as exc:
            errors.append(f"import failed: {mod} :: {exc}")
    try:
        workflow = importlib.import_module("metaboclip.core.workflow")
        getattr(workflow, "run_directory")
        getattr(workflow, "run_single_pair")
    except Exception as exc:
        errors.append(f"unified MetaboClip backend import failed :: {exc}")
    return errors


def _load_tool_registry(root: Path) -> dict:
    if yaml is None:
        return {"tools": {}}
    path = root / "third_party" / "tools.yaml"
    if not path.is_file():
        return {"tools": {}}
    try:
        with path.open("r", encoding="utf-8") as handle:
            return yaml.safe_load(handle) or {"tools": {}}
    except Exception:
        return {"tools": {}}


def _resolve_registered_executable(root: Path, executable: str) -> str:
    executable = str(executable or "").strip()
    if not executable:
        return executable
    if any(ch.isspace() for ch in executable):
        parts = shlex.split(executable)
        if not parts:
            return executable
        first = _resolve_registered_executable(root, parts[0])
        return executable.replace(parts[0], first, 1)
    path = Path(executable).expanduser()
    if path.is_absolute() and path.exists():
        return str(path)
    if not path.is_absolute():
        root_candidate = root / path
        if root_candidate.exists():
            return str(root_candidate)
        local_bin = root / "third_party" / "bin" / executable
        if local_bin.exists():
            return str(local_bin)
    return shutil.which(executable) or executable


def _tool_available(root: Path, registry: dict, tool: str) -> bool:
    tool_cfg = (registry.get("tools") or {}).get(tool, {}) or {}
    executable = _resolve_registered_executable(root, str(tool_cfg.get("executable", tool)))
    first = shlex.split(executable)[0] if executable else tool
    return Path(first).exists() or shutil.which(first) is not None


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate AImd project layout")
    parser.add_argument("--root", default=".", help="AImd project root, default: current directory")
    parser.add_argument("--strict-tools", action="store_true", help="fail when external tools are missing from PATH")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    errors: list[str] = []
    errors.extend(check_layout(root))
    errors.extend(check_yaml(root))
    errors.extend(check_imports(root))

    registry = _load_tool_registry(root)
    missing_tools = [tool for tool in EXTERNAL_TOOLS if not _tool_available(root, registry, tool)]
    print(f"[AImd validate] root: {root}")
    if missing_tools:
        print("[AImd validate] external tools not found through PATH or third_party/tools.yaml:")
        for tool in missing_tools:
            print(f"  - {tool}")
        if args.strict_tools:
            errors.extend([f"missing external tool: {tool}" for tool in missing_tools])
    else:
        print("[AImd validate] all common external tools found through PATH or third_party/tools.yaml")

    missing_optional = [tool for tool in OPTIONAL_EXTERNAL_TOOLS if not _tool_available(root, registry, tool)]
    if missing_optional:
        print("[AImd validate] optional external tools not found through PATH or third_party/tools.yaml:")
        for tool in missing_optional:
            print(f"  - {tool}")

    if errors:
        print("[AImd validate] FAILED")
        for err in errors:
            print(f"  - {err}")
        raise SystemExit(1)

    print("[AImd validate] PASSED")


if __name__ == "__main__":
    main()
