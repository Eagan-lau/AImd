#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from .selector import write_refinement_selection
from .utils import deep_merge, dump_yaml, load_yaml, resolve_path, write_json


def _build_refined_docking_config(config: dict[str, Any], selected_manifest: Path) -> Path:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    base_path = resolve_path(config.get("paths", {}).get("base_docking_config", "configs/Docking/docking.yaml"), root)
    out_config = resolve_path(config.get("paths", {}).get("generated_docking_config", "data/refinement/refined_docking.generated.yaml"), root)
    if base_path is None or not base_path.exists():
        raise FileNotFoundError(f"Base DockingHub config not found: {base_path}")
    assert out_config is not None
    base = load_yaml(base_path)
    default_overrides: dict[str, Any] = {
        "paths": {
            "aimd_root": str(root),
            "protein_manifest": str(selected_manifest),
            "pocket_manifest": "data/pocket/pocket_manifest.csv",
            "ligand_manifest": "data/ligand/ligand_manifest.csv",
            "cofactor_dir": "data/cofactor",
            "ensemble_dir": "data/refined/ensemble",
            "cofactor_mapped_dir": "data/refined/cofactor_mapped",
            "receptor_dir": "data/refined/receptor",
            "docking_config_dir": "data/refined/docking_configs",
            "docking_task_dir": "data/refined/docking_tasks",
            "docking_out_dir": "data/refined/docking_out",
        },
        "selection": {
            "protein_mode": "all",
            "top_n_pockets_per_protein": 1,
        },
        "ensemble": {
            "enabled": True,
            "fallback_to_input": True,
        },
        "cofactor": {
            "enabled": True,
            "use_foldseek": True,
            "fallback_to_first": True,
            "continue_without_cofactor_on_error": True,
        },
        "docking": {
            "engine": "vina",
            "run": True,
            "vina": {"exhaustiveness": 16, "cpu_per_job": 4},
        },
        "output": {"overwrite": False, "file_action": "copy"},
    }
    user_overrides = config.get("docking_overrides", {}) or {}
    merged = deep_merge(base, default_overrides)
    merged = deep_merge(merged, user_overrides)
    dump_yaml(out_config, merged)
    return out_config


def run_refinement(config_path: str | Path) -> Path:
    config_path = Path(config_path).resolve()
    config = load_yaml(config_path)
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    print(f"[RefinementHub] AImd root: {root}")
    selected_manifest, selected_clusters = write_refinement_selection(config)
    print(f"[RefinementHub] selected proteins: {selected_manifest}")
    print(f"[RefinementHub] selected clusters: {selected_clusters}")
    docking_config = _build_refined_docking_config(config, selected_manifest)
    print(f"[RefinementHub] generated refined DockingHub config: {docking_config}")

    run_cfg = config.get("run", {})
    if bool(run_cfg.get("run_dockinghub", True)):
        from DockingHub.main import run_dockinghub
        result = run_dockinghub(docking_config)
    else:
        result = Path(docking_config)
    out_report = resolve_path(config.get("paths", {}).get("report_json", "data/refinement/refinement_report.json"), root)
    assert out_report is not None
    write_json(out_report, {
        "selected_protein_manifest": str(selected_manifest),
        "selected_clusters_csv": str(selected_clusters),
        "generated_docking_config": str(docking_config),
        "refined_docking_result_or_config": str(result),
    })
    print(f"[RefinementHub] finished: {result}")
    return Path(result)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run AImd RefinementHub: select top ClusterScore clusters and re-run DockingHub")
    parser.add_argument("--config", required=True, help="Path to configs/Refinement/refine_from_clusterscore.yaml")
    args = parser.parse_args()
    run_refinement(args.config)


if __name__ == "__main__":
    main()
