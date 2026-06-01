#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise RuntimeError("PyYAML is required. Install with: pip install PyYAML") from exc


def _load_yaml(path: str | Path) -> dict:
    path = Path(path).resolve()
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _resolve(value: str | Path, root: Path) -> Path:
    p = Path(value)
    return p.resolve() if p.is_absolute() else (root / p).resolve()


def run_workflow(config_path: str | Path) -> None:
    config_path = Path(config_path).resolve()
    config = _load_yaml(config_path)
    paths = config.get("paths", {})
    workflow = config.get("workflow", {})
    root = Path(paths.get("aimd_root", ".")).resolve()
    rgpc_config = _resolve(paths.get("rgpc_config", "configs/RGPC/rgpc.yaml"), root)
    tapocket_config = _resolve(paths.get("tapocket_batch_config", "configs/TApocket/tapocket_batch.yaml"), root)
    docking_config = _resolve(paths.get("docking_config", "configs/Docking/docking.yaml"), root)

    if bool(workflow.get("run_rgpc", True)):
        from RGPC.main import run_rgpc
        print(f"[AImdWorkflow] Running RGPC: {rgpc_config}")
        run_rgpc(rgpc_config)
    else:
        print("[AImdWorkflow] Skip RGPC")

    if bool(workflow.get("run_tapocket", True)):
        from TApocketBridge.runner import run_tapocket_batch
        print(f"[AImdWorkflow] Running TApocketBridge: {tapocket_config}")
        run_tapocket_batch(tapocket_config)
    else:
        print("[AImdWorkflow] Skip TApocketBridge")

    if bool(workflow.get("run_docking", True)):
        from DockingHub.main import run_dockinghub
        print(f"[AImdWorkflow] Running DockingHub: {docking_config}")
        run_dockinghub(docking_config)
    else:
        print("[AImdWorkflow] Skip DockingHub")

    print("[AImdWorkflow] RGPC → TApocket → Docking workflow finished.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run AImd RGPC + TApocket + Docking workflow")
    parser.add_argument("--config", required=True, help="Path to configs/workflows/rgpc_tapocket_docking.yaml")
    args = parser.parse_args()
    run_workflow(args.config)


if __name__ == "__main__":
    main()
