#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate the AImd project layout and basic configuration consistency.

This script does not run external tools. It checks:
1. required module directories;
2. required config files;
3. key input/output folders;
4. YAML parseability;
5. Python importability of AImd wrapper modules;
6. availability of common external executables when present in PATH.
"""

from __future__ import annotations

import argparse
import importlib
import shutil
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


REQUIRED_DIRS = [
    "RGPC", "TApocket", "TApocketBridge", "DockingHub",
    "ScoringHub", "RefinementHub", "MetaBoClip", "MetaBoClipBridge",
    "orchestrator", "configs", "data", "docs",
]

REQUIRED_CONFIGS = [
    "configs/RGPC/rgpc.yaml",
    "configs/TApocket/tapocket_batch.yaml",
    "configs/Docking/docking.yaml",
    "configs/Scoring/cluster_score.yaml",
    "configs/Refinement/refine_from_clusterscore.yaml",
    "configs/MetaBoClip/metaboclip_bridge.yaml",
    "configs/workflows/full_iterative_metaboclip.yaml",
]

IMPORT_MODULES = [
    "RGPC.main",
    "TApocketBridge.main",
    "DockingHub.main",
    "ScoringHub.main",
    "RefinementHub.main",
    "MetaBoClipBridge.main",
    "orchestrator.full_iterative_metaboclip",
]

EXTERNAL_TOOLS = [
    "foldseek", "hipmcl", "pymol", "prepare_receptor", "vina",
]


def check_layout(root: Path) -> list[str]:
    errors: list[str] = []
    for rel in REQUIRED_DIRS:
        if not (root / rel).is_dir():
            errors.append(f"missing directory: {rel}")
    for rel in REQUIRED_CONFIGS:
        if not (root / rel).is_file():
            errors.append(f"missing config: {rel}")
    return errors


def check_yaml(root: Path) -> list[str]:
    if yaml is None:
        return ["PyYAML is not installed; cannot parse YAML configs"]
    errors: list[str] = []
    for path in root.rglob("*.yaml"):
        try:
            with path.open("r", encoding="utf-8") as handle:
                yaml.safe_load(handle)
        except Exception as exc:
            errors.append(f"YAML parse error: {path.relative_to(root)} :: {exc}")
    return errors


def check_imports(root: Path) -> list[str]:
    errors: list[str] = []
    sys.path.insert(0, str(root))
    for mod in IMPORT_MODULES:
        try:
            importlib.import_module(mod)
        except Exception as exc:
            errors.append(f"import failed: {mod} :: {exc}")
    return errors


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

    missing_tools = [tool for tool in EXTERNAL_TOOLS if shutil.which(tool) is None]
    print(f"[AImd validate] root: {root}")
    if missing_tools:
        print("[AImd validate] external tools not found in PATH:")
        for tool in missing_tools:
            print(f"  - {tool}")
        if args.strict_tools:
            errors.extend([f"missing external tool: {tool}" for tool in missing_tools])
    else:
        print("[AImd validate] all common external tools found in PATH")

    if errors:
        print("[AImd validate] FAILED")
        for err in errors:
            print(f"  - {err}")
        raise SystemExit(1)

    print("[AImd validate] PASSED")


if __name__ == "__main__":
    main()
