#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import sys
import traceback
from pathlib import Path
from typing import Any

from .manifest import ProteinInput, load_protein_manifest
from .utils import copy_or_link, ensure_dir, load_yaml, read_json, read_tsv, resolve_path, safe_name, write_csv, write_json, write_yaml


FINAL_FILES = [
    "final_pockets.pdb",
    "final_pockets.json",
    "final_pocket_residues.tsv",
    "final_boxes.tsv",
    "summary.tsv",
    "run_summary.json",
    "provenance.json",
    "run.log",
]

RUN_MANIFEST_FIELDS = [
    "protein_id", "cluster_id", "batch_id", "is_representative", "protein_path",
    "tapocket_run_id", "tapocket_run_dir", "pocket_output_dir",
    "selected_template_count", "candidate_pocket_count", "final_pocket_count",
    "mcsa_enabled", "mcsa_action", "ai_enabled", "ai_used", "status", "error",
]

POCKET_MANIFEST_FIELDS = [
    "protein_id", "cluster_id", "batch_id", "pocket_id", "pocket_rank", "source",
    "center_x", "center_y", "center_z", "size_x", "size_y", "size_z", "padding_angstrom",
    "final_score", "query_residue_count", "mapping_coverage", "mapping_method",
    "mcsa_matched", "mcsa_match_basis", "mcsa_atoms_in_box", "mcsa_fraction_in_box", "ai_confidence",
    "protein_path", "pocket_pdb_path", "pocket_json_path", "box_tsv_path", "box_yaml_path",
    "residues_tsv_path", "summary_tsv_path", "tapocket_run_id", "tapocket_run_dir", "pocket_output_dir", "status",
]


def _install_tapocket_path(tapocket_project_dir: Path) -> None:
    """Make bundled AImd/TApocket importable without requiring pip install -e."""
    path = str(tapocket_project_dir.resolve())
    if path not in sys.path:
        sys.path.insert(0, path)


def _detect_final_dir(run_dir: Path) -> Path:
    if (run_dir / "final_boxes.tsv").exists() or (run_dir / "summary.tsv").exists():
        return run_dir
    if (run_dir / "final" / "final_boxes.tsv").exists() or (run_dir / "final" / "summary.tsv").exists():
        return run_dir / "final"
    return run_dir


def _copy_final_outputs(run_dir: Path, pocket_output_dir: Path, copy_mode: str, overwrite: bool) -> dict[str, Path]:
    final_dir = _detect_final_dir(run_dir)
    copied: dict[str, Path] = {}
    ensure_dir(pocket_output_dir)
    for name in FINAL_FILES:
        candidates = [final_dir / name, run_dir / name]
        if name == "run.log":
            candidates.append(run_dir / "logs" / "run.log")
        src = next((p for p in candidates if p.exists()), None)
        if src is None:
            continue
        copied[name] = copy_or_link(src, pocket_output_dir / name, mode=copy_mode, overwrite=overwrite)
    return copied


def _load_run_summary(pocket_output_dir: Path, fallback: dict[str, Any] | None = None) -> dict[str, Any]:
    summary_path = pocket_output_dir / "run_summary.json"
    if summary_path.exists():
        return read_json(summary_path)
    return fallback or {}


def _summary_by_pocket(summary_tsv: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(summary_tsv)
    return {row.get("pocket_id", ""): row for row in rows if row.get("pocket_id")}


def _as_float_or_text(value: Any) -> Any:
    text = "" if value is None else str(value).strip()
    if text == "":
        return ""
    try:
        return float(text)
    except Exception:
        return text


def _write_box_yaml(path: Path, protein: ProteinInput, run_id: str, box_row: dict[str, str]) -> None:
    data = {
        "protein_id": protein.protein_id,
        "cluster_id": protein.cluster_id,
        "batch_id": protein.batch_id,
        "pocket_id": box_row.get("pocket_id", ""),
        "source": box_row.get("source", ""),
        "tapocket_run_id": run_id,
        "docking_box": {
            "center_x": _as_float_or_text(box_row.get("center_x", "")),
            "center_y": _as_float_or_text(box_row.get("center_y", "")),
            "center_z": _as_float_or_text(box_row.get("center_z", "")),
            "size_x": _as_float_or_text(box_row.get("size_x", "")),
            "size_y": _as_float_or_text(box_row.get("size_y", "")),
            "size_z": _as_float_or_text(box_row.get("size_z", "")),
            "padding_angstrom": _as_float_or_text(box_row.get("padding_angstrom", "")),
        },
        "vina_box": {
            "center_x": _as_float_or_text(box_row.get("center_x", "")),
            "center_y": _as_float_or_text(box_row.get("center_y", "")),
            "center_z": _as_float_or_text(box_row.get("center_z", "")),
            "size_x": _as_float_or_text(box_row.get("size_x", "")),
            "size_y": _as_float_or_text(box_row.get("size_y", "")),
            "size_z": _as_float_or_text(box_row.get("size_z", "")),
        },
    }
    write_yaml(path, data)


def _collect_pocket_rows(protein: ProteinInput, run_id: str, run_dir: Path, pocket_output_dir: Path) -> list[dict[str, Any]]:
    boxes_path = pocket_output_dir / "final_boxes.tsv"
    summary_path = pocket_output_dir / "summary.tsv"
    residues_path = pocket_output_dir / "final_pocket_residues.tsv"
    pdb_path = pocket_output_dir / "final_pockets.pdb"
    json_path = pocket_output_dir / "final_pockets.json"
    boxes = read_tsv(boxes_path)
    summary_index = _summary_by_pocket(summary_path)
    rows: list[dict[str, Any]] = []
    for rank, box in enumerate(boxes, start=1):
        pocket_id = box.get("pocket_id", f"pocket_{rank}") or f"pocket_{rank}"
        safe_pocket = safe_name(pocket_id)
        box_yaml = pocket_output_dir / f"box_{safe_pocket}.yaml"
        _write_box_yaml(box_yaml, protein, run_id, box)
        srow = summary_index.get(pocket_id, {})
        rows.append({
            "protein_id": protein.protein_id,
            "cluster_id": protein.cluster_id,
            "batch_id": protein.batch_id,
            "pocket_id": pocket_id,
            "pocket_rank": rank,
            "source": box.get("source", srow.get("source", "")),
            "center_x": box.get("center_x", ""),
            "center_y": box.get("center_y", ""),
            "center_z": box.get("center_z", ""),
            "size_x": box.get("size_x", ""),
            "size_y": box.get("size_y", ""),
            "size_z": box.get("size_z", ""),
            "padding_angstrom": box.get("padding_angstrom", ""),
            "final_score": box.get("final_score", srow.get("final_score", "")),
            "query_residue_count": srow.get("query_residue_count", ""),
            "mapping_coverage": srow.get("mapping_coverage", ""),
            "mapping_method": srow.get("mapping_method", ""),
            "mcsa_matched": box.get("mcsa_matched", srow.get("mcsa_matched", "")),
            "mcsa_match_basis": srow.get("mcsa_match_basis", ""),
            "mcsa_atoms_in_box": srow.get("mcsa_atoms_in_box", ""),
            "mcsa_fraction_in_box": srow.get("mcsa_fraction_in_box", ""),
            "ai_confidence": box.get("ai_confidence", ""),
            "protein_path": str(protein.protein_path),
            "pocket_pdb_path": str(pdb_path) if pdb_path.exists() else "",
            "pocket_json_path": str(json_path) if json_path.exists() else "",
            "box_tsv_path": str(boxes_path) if boxes_path.exists() else "",
            "box_yaml_path": str(box_yaml),
            "residues_tsv_path": str(residues_path) if residues_path.exists() else "",
            "summary_tsv_path": str(summary_path) if summary_path.exists() else "",
            "tapocket_run_id": run_id,
            "tapocket_run_dir": str(run_dir),
            "pocket_output_dir": str(pocket_output_dir),
            "status": "success",
        })
    return rows


def _run_summary_row(protein: ProteinInput, run_id: str, run_dir: Path, pocket_output_dir: Path, summary: dict[str, Any], status: str, error: str = "") -> dict[str, Any]:
    return {
        "protein_id": protein.protein_id,
        "cluster_id": protein.cluster_id,
        "batch_id": protein.batch_id,
        "is_representative": str(protein.is_representative).lower(),
        "protein_path": str(protein.protein_path),
        "tapocket_run_id": run_id,
        "tapocket_run_dir": str(run_dir),
        "pocket_output_dir": str(pocket_output_dir),
        "selected_template_count": summary.get("selected_template_count", ""),
        "candidate_pocket_count": summary.get("candidate_pocket_count", ""),
        "final_pocket_count": summary.get("final_pocket_count", ""),
        "mcsa_enabled": summary.get("mcsa_enabled", ""),
        "mcsa_action": summary.get("mcsa_action", ""),
        "ai_enabled": summary.get("ai_enabled", ""),
        "ai_used": summary.get("ai_used", ""),
        "status": status,
        "error": error,
    }


def _write_batch_summary(out_dir: Path, pocket_rows: list[dict[str, Any]], run_rows: list[dict[str, Any]]) -> None:
    by_batch: dict[str, dict[str, Any]] = {}
    for row in run_rows:
        batch = str(row.get("batch_id", "file_1"))
        by_batch.setdefault(batch, {"batch_id": batch, "n_proteins": 0, "n_successful_runs": 0, "n_failed_runs": 0, "n_pockets": 0})
        by_batch[batch]["n_proteins"] += 1
        if row.get("status") == "success":
            by_batch[batch]["n_successful_runs"] += 1
        else:
            by_batch[batch]["n_failed_runs"] += 1
    for row in pocket_rows:
        batch = str(row.get("batch_id", "file_1"))
        by_batch.setdefault(batch, {"batch_id": batch, "n_proteins": 0, "n_successful_runs": 0, "n_failed_runs": 0, "n_pockets": 0})
        by_batch[batch]["n_pockets"] += 1
    write_csv(out_dir / "pocket_batch_summary.csv", by_batch.values(), ["batch_id", "n_proteins", "n_successful_runs", "n_failed_runs", "n_pockets"])


def prepare_tapocket_databases(config: dict[str, Any], tapocket_project_dir: Path, tapocket_config_path: Path, out_dir: Path) -> None:
    prep_cfg = config.get("prepare", {})
    if not any(bool(prep_cfg.get(k, False)) for k in ["check_layout", "build_manifest", "build_index"]):
        return
    _install_tapocket_path(tapocket_project_dir)
    from tapocket.core.config import load_config
    from tapocket.databases.manifest import check_database_layout, write_manifests

    tap_cfg = load_config(tapocket_config_path)
    report: dict[str, Any] = {}
    if bool(prep_cfg.get("check_layout", False)):
        report["check_layout"] = check_database_layout(tap_cfg)
        write_json(out_dir / "tapocket_check_layout.json", report["check_layout"])
    if bool(prep_cfg.get("build_manifest", False)):
        report["build_manifest"] = write_manifests(tap_cfg)
        write_json(out_dir / "tapocket_build_manifest.json", report["build_manifest"])
    if bool(prep_cfg.get("build_index", False)):
        from tapocket.databases.template_db import TemplateDB
        from tapocket.databases.mcsa_db import MCSADB
        from tapocket.retrievers.foldseek import build_foldseek_db_from_template_records, build_foldseek_db_from_mcsa_records
        binary = tap_cfg.get("foldseek", "binary", default="foldseek")
        db_mode = str(prep_cfg.get("build_index_db", "all"))
        force = bool(prep_cfg.get("force_index", False))
        create_index = bool(prep_cfg.get("create_index", False))
        index_results: list[dict[str, Any]] = []
        if db_mode in {"template", "all"}:
            template_manifest = tap_cfg.path("paths", "template_manifest")
            if not template_manifest.exists():
                write_manifests(tap_cfg)
            template_db = TemplateDB.from_manifest(template_manifest, tap_cfg.root)
            result = build_foldseek_db_from_template_records(
                records=template_db.records,
                root=tap_cfg.root,
                staging_dir=tap_cfg.path("paths", "staging_template_inputs"),
                output_db=tap_cfg.path("paths", "foldseek_template_db"),
                binary=binary,
                create_index=create_index,
                force=force,
            )
            result["kind"] = "template"
            index_results.append(result)
        if db_mode in {"mcsa", "all"}:
            mcsa_manifest = tap_cfg.path("paths", "mcsa_manifest")
            if not mcsa_manifest.exists():
                write_manifests(tap_cfg)
            mcsa_db = MCSADB.from_manifest(mcsa_manifest, tap_cfg.root)
            result = build_foldseek_db_from_mcsa_records(
                records=mcsa_db.records,
                root=tap_cfg.root,
                staging_dir=tap_cfg.path("paths", "staging_mcsa_inputs"),
                output_db=tap_cfg.path("paths", "foldseek_mcsa_db"),
                binary=binary,
                create_index=create_index,
                force=force,
            )
            result["kind"] = "mcsa"
            index_results.append(result)
        report["build_index"] = index_results
        write_json(out_dir / "tapocket_build_index.json", index_results)
    write_json(out_dir / "tapocket_prepare_report.json", report)


def run_tapocket_batch(config_path: str | Path) -> tuple[Path, Path]:
    config_path = Path(config_path).resolve()
    config = load_yaml(config_path)
    paths = config.get("paths", {})
    root = Path(paths.get("aimd_root", ".")).resolve()
    out_dir = resolve_path(paths.get("out_dir", "data/pocket"), root)
    tapocket_project_dir = resolve_path(paths.get("tapocket_project_dir", "TApocket"), root)
    tapocket_config_path = resolve_path(paths.get("tapocket_config", "TApocket/configs/tapocket_template_v1.yaml"), root)
    if out_dir is None or tapocket_project_dir is None or tapocket_config_path is None:
        raise ValueError("paths.out_dir, paths.tapocket_project_dir, and paths.tapocket_config are required")
    ensure_dir(out_dir)
    if not tapocket_project_dir.exists():
        raise FileNotFoundError(f"TApocket project directory not found: {tapocket_project_dir}")
    if not tapocket_config_path.exists():
        raise FileNotFoundError(f"TApocket config not found: {tapocket_config_path}")

    _install_tapocket_path(tapocket_project_dir)
    from tapocket.core.config import load_config
    from tapocket.core.pipeline_template import TApocketPipeline

    prepare_tapocket_databases(config, tapocket_project_dir, tapocket_config_path, out_dir)
    protein_records = load_protein_manifest(config)
    print(f"[TApocketBridge] Selected proteins: {len(protein_records)}")
    print(f"[TApocketBridge] Output pocket directory: {out_dir}")

    tap_cfg = load_config(tapocket_config_path)
    pipeline = TApocketPipeline(tap_cfg)
    execution_cfg = config.get("execution", {})
    copy_mode = str(config.get("output", {}).get("file_action", "copy"))
    overwrite = bool(config.get("output", {}).get("overwrite", True))
    resume = bool(config.get("resume", {}).get("enabled", True))
    skip_existing = bool(config.get("resume", {}).get("skip_existing", True))
    continue_on_error = bool(execution_cfg.get("continue_on_error", True))
    workers = int(execution_cfg.get("workers", 1))
    if workers != 1:
        print("[TApocketBridge][WARN] workers > 1 requested; running sequentially because the PyMOL mapper is safer in serial mode.")

    run_rows: list[dict[str, Any]] = []
    pocket_rows: list[dict[str, Any]] = []
    for idx, protein in enumerate(protein_records, start=1):
        run_id = f"{safe_name(protein.batch_id)}__{safe_name(protein.protein_id)}"
        pocket_output_dir = out_dir / protein.batch_id / safe_name(protein.protein_id)
        existing_summary = pocket_output_dir / "run_summary.json"
        print(f"[TApocketBridge] ({idx}/{len(protein_records)}) {protein.protein_id} -> {run_id}")
        try:
            if resume and skip_existing and existing_summary.exists():
                summary_dict = read_json(existing_summary)
                run_dir = Path(summary_dict.get("run_dir", "")) if summary_dict.get("run_dir") else tap_cfg.path("paths", "run_root") / run_id
                print(f"[TApocketBridge] Skip existing run: {protein.protein_id}")
            else:
                summary = pipeline.run(query=protein.protein_path, run_id=run_id)
                summary_dict = summary.to_dict()
                run_dir = Path(summary_dict["run_dir"]).resolve()
                _copy_final_outputs(run_dir, pocket_output_dir, copy_mode=copy_mode, overwrite=overwrite)
                if not (pocket_output_dir / "run_summary.json").exists():
                    write_json(pocket_output_dir / "run_summary.json", summary_dict)
            summary_dict = _load_run_summary(pocket_output_dir, fallback=summary_dict)
            run_dir = Path(summary_dict.get("run_dir", tap_cfg.path("paths", "run_root") / run_id)).resolve()
            run_rows.append(_run_summary_row(protein, run_id, run_dir, pocket_output_dir, summary_dict, status="success"))
            pocket_rows.extend(_collect_pocket_rows(protein, run_id, run_dir, pocket_output_dir))
        except Exception as exc:
            err = str(exc)
            print(f"[TApocketBridge][ERROR] {protein.protein_id}: {err}")
            if bool(execution_cfg.get("write_traceback", True)):
                ensure_dir(pocket_output_dir)
                (pocket_output_dir / "error.traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
            run_rows.append(_run_summary_row(protein, run_id, tap_cfg.path("paths", "run_root") / run_id, pocket_output_dir, {}, status="failed", error=err))
            if not continue_on_error:
                raise

    run_manifest = out_dir / "tapocket_run_manifest.csv"
    pocket_manifest = out_dir / "pocket_manifest.csv"
    write_csv(run_manifest, run_rows, RUN_MANIFEST_FIELDS)
    write_csv(pocket_manifest, pocket_rows, POCKET_MANIFEST_FIELDS)
    _write_batch_summary(out_dir, pocket_rows, run_rows)
    report = {
        "config": str(config_path),
        "aimd_root": str(root),
        "tapocket_project_dir": str(tapocket_project_dir),
        "tapocket_config": str(tapocket_config_path),
        "selected_proteins": len(protein_records),
        "successful_runs": sum(1 for r in run_rows if r.get("status") == "success"),
        "failed_runs": sum(1 for r in run_rows if r.get("status") != "success"),
        "final_pocket_rows": len(pocket_rows),
        "run_manifest": str(run_manifest),
        "pocket_manifest": str(pocket_manifest),
    }
    write_json(out_dir / "tapocket_batch_report.json", report)
    print("[TApocketBridge] Done.")
    print(f"[TApocketBridge] run_manifest: {run_manifest}")
    print(f"[TApocketBridge] pocket_manifest: {pocket_manifest}")
    return run_manifest, pocket_manifest
