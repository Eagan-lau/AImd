#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
from typing import Any

from .manifest import LigandRecord, PocketRecord, group_pockets_by_protein
from .utils import ensure_dir, read_csv, resolve_path, safe_id, write_csv


def _ligands_for_receptor(config: dict[str, Any], receptor_row: dict[str, Any], ligands: list[LigandRecord]) -> list[LigandRecord]:
    pairing = config.get("pairing", {})
    mode = pairing.get("mode", "all_vs_all")
    batch = receptor_row.get("batch_id", "file_1")
    if mode == "all_vs_all":
        return ligands
    if mode == "by_batch":
        same = [l for l in ligands if l.batch_id == batch]
        # Many projects store all ligands in file_1; use them if no same-batch ligand exists.
        return same or [l for l in ligands if l.batch_id == "file_1"] or ligands
    return ligands


def _explicit_pairs(config: dict[str, Any]) -> set[tuple[str, str]] | None:
    pairing = config.get("pairing", {})
    if pairing.get("mode") != "explicit_pairs":
        return None
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    pair_file = resolve_path(pairing.get("pairs_csv"), root)
    if pair_file is None or not pair_file.exists():
        raise FileNotFoundError(f"Explicit pair file not found: {pair_file}")
    rows = read_csv(pair_file)
    pairs: set[tuple[str, str]] = set()
    for row in rows:
        pairs.add((row.get("ligand_id", ""), row.get("protein_id", "")))
    return pairs


def _write_vina_config(path: Path, receptor: Path, ligand: Path, out_pose: Path, pocket: PocketRecord) -> None:
    ensure_dir(path.parent)
    lines = [
        f"receptor = {receptor}\n",
        f"ligand = {ligand}\n",
        f"out = {out_pose}\n",
        f"center_x = {pocket.center_x:.3f}\n",
        f"center_y = {pocket.center_y:.3f}\n",
        f"center_z = {pocket.center_z:.3f}\n",
        f"size_x = {pocket.size_x:.3f}\n",
        f"size_y = {pocket.size_y:.3f}\n",
        f"size_z = {pocket.size_z:.3f}\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")


def build_docking_tasks(
    config: dict[str, Any],
    receptors: list[dict[str, Any]],
    ligands: list[LigandRecord],
    pockets: list[PocketRecord],
) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    out_dir = resolve_path(config.get("paths", {}).get("docking_out_dir", "data/docking_out"), root)
    conf_dir = resolve_path(config.get("paths", {}).get("docking_config_dir", "data/docking_configs"), root)
    task_dir = resolve_path(config.get("paths", {}).get("docking_task_dir", "data/docking_tasks"), root)
    assert out_dir is not None and conf_dir is not None and task_dir is not None
    ensure_dir(out_dir); ensure_dir(conf_dir); ensure_dir(task_dir)
    pocket_top_n = config.get("selection", {}).get("top_n_pockets_per_protein", 1)
    pocket_top_n = int(pocket_top_n) if pocket_top_n not in {None, "", 0, "0"} else None
    pocket_by_protein = group_pockets_by_protein(pockets, pocket_top_n)
    explicit_pairs = _explicit_pairs(config)
    tasks: list[dict[str, Any]] = []

    for rec in receptors:
        if rec.get("receptor_preparation_status") != "success":
            continue
        receptor_path = Path(rec.get("receptor_pdbqt_path", ""))
        if not receptor_path.exists():
            continue
        pid = rec.get("protein_id", receptor_path.stem)
        batch = rec.get("batch_id", "file_1") or "file_1"
        conformer_id = rec.get("conformer_id", "conf_0")
        rec_pockets = pocket_by_protein.get(pid, [])
        if not rec_pockets:
            continue
        for ligand in _ligands_for_receptor(config, rec, ligands):
            if explicit_pairs is not None and (ligand.ligand_id, pid) not in explicit_pairs:
                continue
            for pocket in rec_pockets:
                job_core = f"{safe_id(ligand.ligand_id)}@{safe_id(pid)}@{safe_id(pocket.pocket_id)}@{safe_id(conformer_id)}"
                batch_out = out_dir / batch
                batch_conf = conf_dir / batch
                ensure_dir(batch_out); ensure_dir(batch_conf)
                out_pose = batch_out / f"{job_core}_out.pdbqt"
                log_path = batch_out / f"{job_core}.out"
                config_path = batch_conf / f"{job_core}.txt"
                _write_vina_config(config_path, receptor_path, ligand.ligand_path, out_pose, pocket)
                tasks.append({
                    "job_id": job_core,
                    "ligand_id": ligand.ligand_id,
                    "protein_id": pid,
                    "cluster_id": rec.get("cluster_id", ""),
                    "batch_id": batch,
                    "conformer_id": conformer_id,
                    "pocket_id": pocket.pocket_id,
                    "pocket_rank": pocket.pocket_rank,
                    "receptor_pdbqt_path": str(receptor_path),
                    "ligand_pdbqt_path": str(ligand.ligand_path),
                    "config_path": str(config_path),
                    "out_pose_path": str(out_pose),
                    "log_path": str(log_path),
                    "center_x": pocket.center_x,
                    "center_y": pocket.center_y,
                    "center_z": pocket.center_z,
                    "size_x": pocket.size_x,
                    "size_y": pocket.size_y,
                    "size_z": pocket.size_z,
                    "status": "pending",
                })
    write_csv(task_dir / "docking_task_manifest.csv", tasks)
    # Also write per-batch task files for scheduler-friendly operation.
    by_batch: dict[str, list[dict[str, Any]]] = {}
    for task in tasks:
        by_batch.setdefault(task["batch_id"], []).append(task)
    for batch, rows in by_batch.items():
        write_csv(task_dir / f"{batch}_tasks.csv", rows)
    return tasks
