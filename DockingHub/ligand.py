#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Any

from .utils import ensure_dir, quote, resolve_path, safe_id, write_csv, resolve_external_tool


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_source_rows(source_csv: Path, prep_cfg: dict[str, Any]) -> list[dict[str, Any]]:
    rows = _read_csv(source_csv)
    smiles_col = prep_cfg.get("smiles_column", "smiles")
    id_col = prep_cfg.get("id_column", "source_id")
    name_col = prep_cfg.get("name_column", "molecule_name")
    max_ligands = prep_cfg.get("max_ligands")
    max_ligands = int(max_ligands) if max_ligands not in {None, "", 0, "0"} else None
    selected: list[dict[str, Any]] = []
    for idx, row in enumerate(rows):
        smiles = str(row.get(smiles_col, "") or "").strip()
        if not smiles:
            continue
        ligand_id = str(row.get(id_col, "") or "").strip() or str(idx + 1)
        selected.append({
            "ligand_id": ligand_id,
            "name": str(row.get(name_col, "") or "").strip() or ligand_id,
            "smiles": smiles,
            "source_row": idx,
        })
        if max_ligands is not None and len(selected) >= max_ligands:
            break
    return selected


def _embed_and_optimize(row: dict[str, Any], sdf_path: Path, pdb_path: Path, prep_cfg: dict[str, Any]) -> dict[str, Any]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:
        return {"status": "failed", "message": f"RDKit import failed: {exc}"}

    smiles = str(row["smiles"])
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"status": "failed", "message": "RDKit could not parse SMILES"}

    mol = Chem.AddHs(mol)
    rdkit_cfg = prep_cfg.get("rdkit", {}) or {}
    seed = int(rdkit_cfg.get("random_seed", 61453))
    max_attempts = int(rdkit_cfg.get("max_attempts", 20))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.maxAttempts = max_attempts
    embed_status = int(AllChem.EmbedMolecule(mol, params))
    if embed_status != 0:
        return {"status": "failed", "message": f"RDKit embedding failed with code {embed_status}"}

    optimize = bool(rdkit_cfg.get("optimize", True))
    force_field = str(rdkit_cfg.get("force_field", "MMFF")).upper()
    max_iters = int(rdkit_cfg.get("max_iters", 500))
    optimization_status = "not_run"
    if optimize:
        try:
            if force_field == "MMFF" and AllChem.MMFFHasAllMoleculeParams(mol):
                rc = AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
                optimization_status = f"MMFF:{rc}"
            else:
                rc = AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
                optimization_status = f"UFF:{rc}"
        except Exception as exc:
            return {"status": "failed", "message": f"RDKit optimization failed: {exc}"}

    ensure_dir(sdf_path.parent)
    ensure_dir(pdb_path.parent)
    writer = Chem.SDWriter(str(sdf_path))
    mol.SetProp("_Name", str(row["ligand_id"]))
    mol.SetProp("smiles", smiles)
    writer.write(mol)
    writer.close()
    Chem.MolToPDBFile(mol, str(pdb_path))
    return {
        "status": "success",
        "message": "3D conformer generated",
        "embedding_status": str(embed_status),
        "optimization_status": optimization_status,
    }


def _prepare_pdbqt(input_pdb: Path, output_pdbqt: Path, prep_cfg: dict[str, Any], root: Path, config: dict[str, Any]) -> tuple[str, str, str]:
    pdbqt_cfg = prep_cfg.get("pdbqt", {}) or {}
    executable = resolve_external_tool("prepare_ligand", pdbqt_cfg.get("executable", "auto"), root, config)
    timeout = int(pdbqt_cfg.get("timeout", 3600) or 3600)
    extra_options = str(pdbqt_cfg.get("extra_options", "") or "").strip()
    command_template = pdbqt_cfg.get("command_template", "{executable} -l {input_pdb} -o {output_pdbqt} {extra_options}")
    command = command_template.format(
        executable=quote(executable),
        input_pdb=quote(input_pdb.name),
        output_pdbqt=quote(output_pdbqt),
        extra_options=extra_options,
    ).strip()
    ensure_dir(output_pdbqt.parent)
    log_path = output_pdbqt.with_suffix(".prepare_ligand.log")
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"[DockingHub] command: {command}\n")
        try:
            proc = subprocess.run(
                command,
                shell=True,
                cwd=str(input_pdb.parent),
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=timeout,
            )
            rc = int(proc.returncode)
        except subprocess.TimeoutExpired:
            log.write(f"\n[DockingHub] TIMEOUT after {timeout}s\n")
            return "failed", "124", "prepare_ligand timed out"
    if rc == 0 and output_pdbqt.exists():
        return "success", str(rc), "PDBQT prepared"
    return "failed", str(rc), f"prepare_ligand returned {rc}"


def prepare_ligands(config: dict[str, Any]) -> Path | None:
    prep_cfg = config.get("ligand_preparation", {}) or {}
    if not bool(prep_cfg.get("enabled", False)):
        return None

    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    source_csv = resolve_path(prep_cfg.get("source_csv", "data/data_input/ligand/taxane_molecules.csv"), root)
    output_dir = resolve_path(prep_cfg.get("output_dir", "data/data_output/ligand_preparation"), root)
    manifest_path = resolve_path(config.get("paths", {}).get("ligand_manifest", "data/data_output/ligand_preparation/ligand_manifest.csv"), root)
    if source_csv is None or not source_csv.exists():
        raise FileNotFoundError(f"Ligand source CSV not found: {source_csv}")
    if output_dir is None or manifest_path is None:
        raise ValueError("ligand_preparation.output_dir and paths.ligand_manifest are required")

    batch_id = str(prep_cfg.get("batch_id", "file_1") or "file_1")
    overwrite = bool(prep_cfg.get("overwrite", False))
    sdf_dir = output_dir / "sdf" / batch_id
    pdb_dir = output_dir / "pdb" / batch_id
    pdbqt_dir = output_dir / "pdbqt" / batch_id
    rows: list[dict[str, Any]] = []

    for src_row in _load_source_rows(source_csv, prep_cfg):
        ligand_id = safe_id(str(src_row["ligand_id"]))
        sdf_path = sdf_dir / f"{ligand_id}.sdf"
        pdb_path = pdb_dir / f"{ligand_id}.pdb"
        pdbqt_path = pdbqt_dir / f"{ligand_id}.pdbqt"
        status = "success"
        message = ""
        embed_info: dict[str, Any] = {
            "embedding_status": "",
            "optimization_status": "",
        }

        if overwrite or not (sdf_path.exists() and pdb_path.exists()):
            embed_info = _embed_and_optimize(src_row, sdf_path, pdb_path, prep_cfg)
            status = str(embed_info.get("status", "failed"))
            message = str(embed_info.get("message", ""))
        else:
            message = "3D conformer already exists"
            embed_info["embedding_status"] = "skipped_existing"
            embed_info["optimization_status"] = "skipped_existing"

        prepare_rc = ""
        if status == "success":
            if overwrite or not pdbqt_path.exists():
                pdbqt_status, prepare_rc, pdbqt_message = _prepare_pdbqt(pdb_path, pdbqt_path, prep_cfg, root, config)
                status = pdbqt_status
                message = pdbqt_message
            else:
                message = "PDBQT already exists"
                prepare_rc = "0"

        rows.append({
            "ligand_id": ligand_id,
            "batch_id": batch_id,
            "ligand_path": str(pdbqt_path) if pdbqt_path.exists() else "",
            "pdbqt_path": str(pdbqt_path) if pdbqt_path.exists() else "",
            "source_csv": str(source_csv),
            "source_row": src_row.get("source_row", ""),
            "smiles": src_row.get("smiles", ""),
            "name": src_row.get("name", ""),
            "sdf_path": str(sdf_path) if sdf_path.exists() else "",
            "pdb_path": str(pdb_path) if pdb_path.exists() else "",
            "embedding_status": embed_info.get("embedding_status", ""),
            "optimization_status": embed_info.get("optimization_status", ""),
            "prepare_ligand_return_code": prepare_rc,
            "status": status,
            "message": message,
        })

    fieldnames = [
        "ligand_id", "batch_id", "ligand_path", "pdbqt_path", "source_csv", "source_row",
        "smiles", "name", "sdf_path", "pdb_path", "embedding_status", "optimization_status",
        "prepare_ligand_return_code", "status", "message",
    ]
    write_csv(manifest_path, rows, fieldnames)
    success_count = sum(1 for row in rows if row.get("status") == "success")
    print(f"[DockingHub] ligand preparation: {success_count}/{len(rows)} successful -> {manifest_path}")
    return manifest_path
