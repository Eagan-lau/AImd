#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path

from .cofactor import build_cofactor_manifest
from .ensemble import build_conformer_manifest
from .manifest import load_ligands, load_pockets, load_proteins
from .postprocess import summarize_results
from .receptor import prepare_receptors
from .tasks import build_docking_tasks
from .utils import load_yaml, write_json
from .vina import run_vina_tasks


def run_dockinghub(config_path: str | Path) -> Path:
    config_path = Path(config_path).resolve()
    config = load_yaml(config_path)
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    print(f"[DockingHub] AImd root: {root}")
    proteins = load_proteins(config)
    ligands = load_ligands(config)
    pockets = load_pockets(config)
    print(f"[DockingHub] inputs: proteins={len(proteins)}, ligands={len(ligands)}, pockets={len(pockets)}")

    conformers = build_conformer_manifest(config, proteins)
    print(f"[DockingHub] conformers: {len(conformers)}")
    cofactors = build_cofactor_manifest(config, conformers)
    print(f"[DockingHub] cofactor/receptor structures: {len(cofactors)}")
    receptors = prepare_receptors(config, cofactors)
    ready = sum(1 for r in receptors if r.get("receptor_preparation_status") == "success")
    print(f"[DockingHub] prepared receptors: {ready}/{len(receptors)}")

    tasks = build_docking_tasks(config, receptors, ligands, pockets)
    print(f"[DockingHub] docking tasks: {len(tasks)}")
    run_rows = run_vina_tasks(config, tasks)
    result_rows = summarize_results(config, run_rows)
    result_manifest = root / config.get("paths", {}).get("docking_out_dir", "data/data_output/docking_out") / "docking_result_manifest.csv"
    write_json(root / config.get("paths", {}).get("docking_out_dir", "data/data_output/docking_out") / "run_config_snapshot.json", config)
    print(f"[DockingHub] finished: {result_manifest}")
    return result_manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Run AImd DockingHub: ensemble/cofactor-aware AutoDock Vina docking")
    parser.add_argument("--config", required=True, help="Path to configs/Docking/docking.yaml")
    args = parser.parse_args()
    run_dockinghub(args.config)


if __name__ == "__main__":
    main()
