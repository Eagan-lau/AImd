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
    docking_config = _resolve(paths.get("broad_docking_config", "configs/Docking/docking.yaml"), root)
    clusterscore_config = _resolve(paths.get("clusterscore_config", "configs/Scoring/cluster_score.yaml"), root)
    refinement_config = _resolve(paths.get("refinement_config", "configs/Refinement/refine_from_clusterscore.yaml"), root)
    metaboclip_config = _resolve(paths.get("metaboclip_bridge_config", "configs/MetaBoClip/metaboclip_bridge.yaml"), root)

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

    if bool(workflow.get("run_broad_docking", True)):
        from DockingHub.main import run_dockinghub
        print(f"[AImdWorkflow] Running broad DockingHub: {docking_config}")
        run_dockinghub(docking_config)
    else:
        print("[AImdWorkflow] Skip broad DockingHub")

    if bool(workflow.get("run_clusterscore", True)):
        from ScoringHub.main import run_clusterscore
        print(f"[AImdWorkflow] Running ClusterScore: {clusterscore_config}")
        run_clusterscore(clusterscore_config)
    else:
        print("[AImdWorkflow] Skip ClusterScore")

    if bool(workflow.get("run_refinement", True)):
        from RefinementHub.main import run_refinement
        print(f"[AImdWorkflow] Running RefinementHub: {refinement_config}")
        run_refinement(refinement_config)
    else:
        print("[AImdWorkflow] Skip RefinementHub")

    if bool(workflow.get("run_metaboclip", True)):
        from MetaBoClipBridge.main import run_metaboclip_bridge
        print(f"[AImdWorkflow] Running MetaBoClipBridge: {metaboclip_config}")
        run_metaboclip_bridge(metaboclip_config)
    else:
        print("[AImdWorkflow] Skip MetaBoClipBridge")

    print("[AImdWorkflow] Full iterative MetaBoClip workflow finished.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run full AImd iterative workflow with refined docking and MetaBoClip scoring")
    parser.add_argument("--config", required=True, help="Path to configs/workflows/full_iterative_metaboclip.yaml")
    args = parser.parse_args()
    run_workflow(args.config)


if __name__ == "__main__":
    main()
