#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

from .utils import copy_or_link, deep_merge, dump_yaml, ensure_dir, load_yaml, numeric, read_csv, resolve_path, safe_id, write_csv, write_json

DEFAULT_FAMILIES = ["cyp450", "fe2og", "ugt", "act", "ach"]
PROTEIN_DELIMITER = "__AIMD__"


def _truthy_status(status: str) -> bool:
    s = str(status or "success").strip().lower()
    return s in {"", "success", "skipped", "ok", "true", "1"}


def load_refined_docking_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    manifest = resolve_path(config.get("paths", {}).get("refined_docking_manifest", "data/refined/docking_out/docking_result_manifest.csv"), root)
    if manifest is None or not manifest.exists():
        raise FileNotFoundError(f"Refined docking result manifest not found: {manifest}")
    rows = read_csv(manifest)
    if not rows:
        raise RuntimeError(f"Refined docking result manifest is empty: {manifest}")
    filtering = config.get("filtering", {})
    require_success = bool(filtering.get("require_success_status", True))
    require_files = bool(filtering.get("require_existing_pose_and_log", True))
    min_affinity = filtering.get("max_affinity_kcal")
    min_affinity = float(min_affinity) if min_affinity not in {None, ""} else None
    out: list[dict[str, Any]] = []
    for row in rows:
        if require_success and not _truthy_status(row.get("status", "")):
            continue
        pose = Path(row.get("out_pose_path", "")) if row.get("out_pose_path") else None
        log = Path(row.get("log_path", "")) if row.get("log_path") else None
        receptor = Path(row.get("receptor_pdbqt_path", "")) if row.get("receptor_pdbqt_path") else None
        if require_files and not (pose and pose.exists() and log and log.exists() and receptor and receptor.exists()):
            continue
        aff = numeric(row.get("best_affinity"), None)
        if min_affinity is not None and aff is not None and aff > min_affinity:
            continue
        out.append(row)
    if not out:
        raise RuntimeError("No refined docking rows passed MetaBoClipBridge filtering")
    return out


def load_cluster_family_map(config: dict[str, Any]) -> dict[str, list[str]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    map_path = resolve_path(config.get("family_assignment", {}).get("cluster_family_map_csv", "data/input/cluster_family_map.csv"), root)
    mapping: dict[str, list[str]] = {}
    if map_path and map_path.exists():
        for row in read_csv(map_path):
            cid = str(row.get("cluster_id") or row.get("cluster") or "").strip()
            family_raw = str(row.get("family") or row.get("families") or "").strip()
            if not cid or not family_raw:
                continue
            fams = [f.strip().lower() for f in family_raw.replace(";", ",").split(",") if f.strip()]
            mapping[cid] = fams
    return mapping


def families_for_row(config: dict[str, Any], row: dict[str, Any], cluster_map: dict[str, list[str]]) -> list[str]:
    assign = config.get("family_assignment", {})
    mode = str(assign.get("mode", "run_all")).lower()
    if mode == "run_all":
        return [str(f).lower() for f in assign.get("all_families", DEFAULT_FAMILIES)]
    if mode == "fixed":
        return [str(f).lower() for f in assign.get("fixed_families", DEFAULT_FAMILIES)]
    if mode == "cluster_family_map":
        fams = cluster_map.get(str(row.get("cluster_id", "")), [])
        if fams:
            return fams
        if bool(assign.get("fallback_to_all_when_unmapped", False)):
            return [str(f).lower() for f in assign.get("all_families", DEFAULT_FAMILIES)]
        return []
    if mode == "column":
        col = str(assign.get("family_column", "family"))
        raw = str(row.get(col, "")).strip()
        return [f.strip().lower() for f in raw.replace(";", ",").split(",") if f.strip()]
    raise ValueError(f"Unsupported family_assignment.mode: {mode}")


def _protein_name(row: dict[str, Any]) -> str:
    protein = safe_id(row.get("protein_id", "protein"))
    pocket = safe_id(row.get("pocket_id", "P1"))
    conformer = safe_id(row.get("conformer_id", "conf_0"))
    return f"{protein}{PROTEIN_DELIMITER}{pocket}__{conformer}"


def stage_metaboclip_inputs(config: dict[str, Any]) -> tuple[Path, list[dict[str, Any]]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    staging_root = resolve_path(config.get("paths", {}).get("staging_dir", "data/metaboclip/staging"), root)
    assert staging_root is not None
    file_action = config.get("output", {}).get("file_action", "copy")
    overwrite = bool(config.get("output", {}).get("overwrite", False))
    cluster_map = load_cluster_family_map(config)
    refined_rows = load_refined_docking_rows(config)
    staged: list[dict[str, Any]] = []
    for row in refined_rows:
        families = families_for_row(config, row, cluster_map)
        if not families:
            continue
        ligand_id_raw = str(row.get("ligand_id", ""))
        ligand_id = safe_id(ligand_id_raw)
        protein_name = _protein_name(row)
        pose = Path(row["out_pose_path"])
        log = Path(row["log_path"])
        receptor = Path(row["receptor_pdbqt_path"])
        for family in families:
            family = safe_id(family.lower())
            unit_root = staging_root / family / ligand_id
            protein_dir = unit_root / "proteins"
            docking_dir = unit_root / "docking"
            receptor_dst = protein_dir / f"{protein_name}.pdbqt"
            ligand_dst = docking_dir / f"{ligand_id}@{protein_name}.pdbqt"
            log_dst = docking_dir / f"{ligand_id}_{protein_name}_cavity_1.out"
            copy_or_link(receptor, receptor_dst, file_action, overwrite=overwrite)
            copy_or_link(pose, ligand_dst, file_action, overwrite=overwrite)
            copy_or_link(log, log_dst, file_action, overwrite=overwrite)
            staged.append({
                "family": family,
                "ligand_id": ligand_id,
                "original_ligand_id": ligand_id_raw,
                "protein_id": row.get("protein_id", ""),
                "cluster_id": row.get("cluster_id", ""),
                "protein_name": protein_name,
                "pocket_id": row.get("pocket_id", ""),
                "conformer_id": row.get("conformer_id", ""),
                "job_id": row.get("job_id", ""),
                "best_affinity": row.get("best_affinity", ""),
                "protein_dir": str(protein_dir),
                "docking_dir": str(docking_dir),
                "staged_receptor_path": str(receptor_dst),
                "staged_pose_path": str(ligand_dst),
                "staged_log_path": str(log_dst),
                "source_receptor_pdbqt_path": str(receptor),
                "source_out_pose_path": str(pose),
                "source_log_path": str(log),
            })
    if not staged:
        raise RuntimeError("No rows were staged for MetaBoClip. Check family_assignment and refined docking results.")
    manifest_path = staging_root / "metaboclip_staging_manifest.csv"
    write_csv(manifest_path, staged)
    return manifest_path, staged


def _load_family_overlay(config: dict[str, Any], family: str, ligand_id: str) -> dict[str, Any]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    family_cfgs = config.get("family_instance_configs", {}) or {}
    overlay: dict[str, Any] = {}
    cfg_path = family_cfgs.get(family)
    if cfg_path:
        resolved = resolve_path(cfg_path, root)
        if resolved and resolved.exists():
            overlay = load_yaml(resolved)
    common_overlay = config.get("metaboclip_config_overlay", {}) or {}
    runtime_overlay = {
        "family": {"key": family, "plugin": family, "runner": family},
        "input": {
            "ligand_id": ligand_id,
            "protein_id_delimiter": PROTEIN_DELIMITER,
            "ligand_pattern": "{ligand_id}@{protein_name}{ligand_extension}",
            "out_pattern": "{ligand_id}_{protein_name}_cavity_1.out",
        },
    }
    return deep_merge(deep_merge(overlay, common_overlay), runtime_overlay)


def run_metaboclip_family_scoring(config: dict[str, Any], staged_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    meta_root = resolve_path(config.get("paths", {}).get("metaboclip_project_dir", "MetaBoClip"), root)
    results_root = resolve_path(config.get("paths", {}).get("results_dir", "data/metaboclip/results"), root)
    assert meta_root is not None and results_root is not None
    if not meta_root.exists():
        raise FileNotFoundError(f"MetaBoClip project directory not found: {meta_root}")
    if str(meta_root) not in sys.path:
        sys.path.insert(0, str(meta_root))
    from core.runtime import run_family  # type: ignore

    plugin_dirs = []
    for raw in config.get("plugin_dirs", []) or []:
        p = resolve_path(raw, root)
        if p:
            plugin_dirs.append(str(p))

    # Unique family-ligand runs.
    units: dict[tuple[str, str], dict[str, Any]] = {}
    for row in staged_rows:
        key = (row["family"], row["ligand_id"])
        units.setdefault(key, {
            "family": row["family"],
            "ligand_id": row["ligand_id"],
            "protein_dir": row["protein_dir"],
            "docking_dir": row["docking_dir"],
        })
    run_rows: list[dict[str, Any]] = []
    for (family, ligand_id), unit in sorted(units.items()):
        out_dir = results_root / family / ligand_id
        ensure_dir(out_dir)
        cfg = _load_family_overlay(config, family, ligand_id)
        effective_cfg_path = out_dir / "effective_metaboclip_config.yaml"
        dump_yaml(effective_cfg_path, cfg)
        try:
            result = run_family(
                cfg=cfg,
                protein_dir=Path(unit["protein_dir"]),
                docking_dir=Path(unit["docking_dir"]),
                out_dir=out_dir,
                plugin_dirs=plugin_dirs,
                explicit_family=family,
            )
            status = "success"
            message = "metaboclip completed"
        except Exception as exc:
            result = {}
            status = "failed"
            message = str(exc)
            if not bool(config.get("execution", {}).get("continue_on_error", True)):
                raise
        gating_csv = out_dir / f"file_{ligand_id}.gating.csv"
        pose_csv = out_dir / f"file_{ligand_id}.pose_scores.csv"
        conf_csv = out_dir / f"file_{ligand_id}.conformation_scores.csv"
        prot_csv = out_dir / f"file_{ligand_id}.protein_scores.csv"
        run_rows.append({
            "family": family,
            "ligand_id": ligand_id,
            "protein_dir": unit["protein_dir"],
            "docking_dir": unit["docking_dir"],
            "out_dir": str(out_dir),
            "effective_config": str(effective_cfg_path),
            "gating_csv": str(gating_csv) if gating_csv.exists() else "",
            "pose_scores_csv": str(pose_csv) if pose_csv.exists() else "",
            "conformation_scores_csv": str(conf_csv) if conf_csv.exists() else "",
            "protein_scores_csv": str(prot_csv) if prot_csv.exists() else "",
            "rows": result.get("rows", ""),
            "status": status,
            "message": message,
        })
    write_csv(results_root / "metaboclip_run_manifest.csv", run_rows)
    return run_rows


def _read_many_csv(rows: list[dict[str, Any]], csv_key: str, extra_keys: list[str]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for run in rows:
        path = Path(run.get(csv_key, "")) if run.get(csv_key) else None
        if path is None or not path.exists():
            continue
        for row in read_csv(path):
            enriched = {k: run.get(k, "") for k in extra_keys}
            enriched.update(row)
            out.append(enriched)
    return out


def collect_metaboclip_outputs(config: dict[str, Any], run_rows: list[dict[str, Any]]) -> Path:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    results_root = resolve_path(config.get("paths", {}).get("results_dir", "data/metaboclip/results"), root)
    assert results_root is not None
    extra = ["family", "ligand_id", "out_dir"]
    gating = _read_many_csv(run_rows, "gating_csv", extra)
    pose = _read_many_csv(run_rows, "pose_scores_csv", extra)
    conf = _read_many_csv(run_rows, "conformation_scores_csv", extra)
    prot = _read_many_csv(run_rows, "protein_scores_csv", extra)
    write_csv(results_root / "metaboclip_gating_all.csv", gating)
    write_csv(results_root / "metaboclip_pose_scores_all.csv", pose)
    write_csv(results_root / "metaboclip_conformation_scores_all.csv", conf)
    # final ranking from protein scores.
    for row in prot:
        row["protein_score_norm_float"] = numeric(row.get("protein_score_norm"), 0.0) or 0.0
        row["max_s_r_float"] = numeric(row.get("max_s_r"), 0.0) or 0.0
    prot.sort(key=lambda r: (float(r.get("protein_score_norm_float", 0.0)), float(r.get("max_s_r_float", 0.0))), reverse=True)
    for i, row in enumerate(prot, start=1):
        row["rank"] = i
    write_csv(results_root / "metaboclip_protein_scores_all.csv", prot)
    final_path = results_root / "metaboclip_final_ranking.csv"
    write_csv(final_path, prot)
    write_json(results_root / "metaboclip_report.json", {
        "n_family_ligand_runs": len(run_rows),
        "n_gating_rows": len(gating),
        "n_pose_score_rows": len(pose),
        "n_conformation_score_rows": len(conf),
        "n_protein_score_rows": len(prot),
        "final_ranking": str(final_path),
    })
    return final_path


def run_metaboclip_bridge_core(config: dict[str, Any]) -> Path:
    staging_manifest, staged = stage_metaboclip_inputs(config)
    print(f"[MetaBoClipBridge] staged inputs: {staging_manifest} ({len(staged)} rows)")
    if bool(config.get("execution", {}).get("run_metaboclip", True)):
        run_rows = run_metaboclip_family_scoring(config, staged)
    else:
        # Create dry-run units only.
        units = {}
        for row in staged:
            key = (row["family"], row["ligand_id"])
            units[key] = row
        run_rows = [{"family": k[0], "ligand_id": k[1], "protein_dir": v["protein_dir"], "docking_dir": v["docking_dir"], "status": "not_run", "message": "execution.run_metaboclip=false"} for k, v in units.items()]
        results_root = resolve_path(config.get("paths", {}).get("results_dir", "data/metaboclip/results"), Path(config.get("paths", {}).get("aimd_root", ".")).resolve())
        assert results_root is not None
        write_csv(results_root / "metaboclip_run_manifest.csv", run_rows)
    final_path = collect_metaboclip_outputs(config, run_rows)
    print(f"[MetaBoClipBridge] final ranking: {final_path}")
    return final_path
