#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

from .utils import copy_or_link, numeric, read_csv, resolve_path, safe_id, write_csv, write_json

DEFAULT_FAMILIES = ["cyp450", "fe2og", "ugt", "act", "ach"]
PROTEIN_DELIMITER = "__AIMD__"
SUPPORTED_BACKENDS = {"unified"}
PREPARE_ONLY_MESSAGE = "Unified backend inputs validated; scoring was not executed"

AIMD_DOCKING_METADATA_FIELDS = [
    "job_id",
    "ligand_id",
    "protein_id",
    "cluster_id",
    "batch_id",
    "conformer_id",
    "pocket_id",
    "pocket_rank",
    "receptor_pdbqt_path",
    "ligand_pdbqt_path",
    "config_path",
    "out_pose_path",
    "log_path",
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
    "return_code",
    "best_affinity",
    "affinities",
    "n_affinities",
    "grid_size",
    "grid_space",
    "exhaustiveness",
    "random_seed",
    "pose_exists",
]

RUN_METADATA_FIELDS = [
    "backend",
    "family",
    "ligand_id",
    "original_ligand_id",
    "protein_id",
    "cluster_id",
    "batch_id",
    "protein_name",
    "pocket_id",
    "pocket_rank",
    "conformer_id",
    "pose_id",
    "job_id",
    "best_affinity",
    "affinity_kcal",
    "protein_dir",
    "docking_dir",
    "out_dir",
    "source_receptor_pdbqt_path",
    "source_out_pose_path",
    "source_log_path",
    "staged_receptor_path",
    "staged_pose_path",
    "staged_log_path",
    "role_table_path",
    "annotation_json_path",
    "atom_map_json_path",
    "mechanism_path",
    "profile_path",
    "status",
    "message",
]

ROOT_OUTPUT_SPECS = {
    "protein_scores": ("protein_scores.csv", "metaboclip_protein_scores_all.csv"),
    "merged_conformation_scores": ("merged_conformation_scores.csv", "metaboclip_conformation_scores_all.csv"),
}

PAIR_OUTPUT_SPECS = {
    "pose_scores": ("pose_scores.csv", "metaboclip_pose_scores_all.csv"),
    "candidate_scores": ("candidate_scores.csv", "metaboclip_candidate_scores_all.csv"),
    "passing_candidates": ("passing_candidates.csv", "metaboclip_passing_candidates_all.csv"),
    "geometry_features": ("geometry_features.csv", "metaboclip_geometry_features_all.csv"),
    "resolved_ligand_sites": ("resolved_ligand_sites.csv", "metaboclip_resolved_ligand_sites_all.csv"),
    "resolved_protein_roles": ("resolved_protein_roles.csv", "metaboclip_resolved_protein_roles_all.csv"),
}


def _backend(config: dict[str, Any]) -> str:
    return str(config.get("backend", "unified")).strip().lower()


def _root(config: dict[str, Any]) -> Path:
    return Path(config.get("paths", {}).get("aimd_root", ".")).resolve()


def _require_unified_backend(config: dict[str, Any]) -> None:
    backend = _backend(config)
    if backend not in SUPPORTED_BACKENDS:
        raise RuntimeError(
            f"Unsupported MetaBoClipBridge backend: {backend}. "
            "The legacy MetaBoClip backend is disabled; set backend: unified."
        )


def _resolve_config_path(config: dict[str, Any], key: str, default: str) -> Path:
    path = resolve_path(config.get("paths", {}).get(key, default), _root(config))
    assert path is not None
    return path


def _unified_paths(config: dict[str, Any]) -> dict[str, Path]:
    return {
        "metaboclip_project_dir": _resolve_config_path(config, "metaboclip_project_dir", "metaboclip_unified"),
        "metaboclip_profile": _resolve_config_path(
            config,
            "metaboclip_profile",
            "metaboclip_unified/metaboclip/config/profiles/default_profile.yaml",
        ),
        "role_table_dir": _resolve_config_path(config, "role_table_dir", "data/data_output/metaboclip/ligand_roles/role_tables"),
        "annotation_dir": _resolve_config_path(config, "annotation_dir", "data/data_output/metaboclip/ligand_roles/annotations"),
        "atom_map_dir": _resolve_config_path(config, "atom_map_dir", "data/data_output/metaboclip/ligand_roles/atom_maps"),
        "ligand_manifest": _resolve_config_path(config, "ligand_manifest", "data/data_input/ligand/ligand_manifest.csv"),
        "ligand_source_manifest": _resolve_config_path(config, "ligand_source_manifest", "data/data_input/ligand/ligand_source_manifest.csv"),
        "unified_output_dir": _resolve_config_path(config, "unified_output_dir", "data/data_output/metaboclip/unified_runs"),
    }


def _ensure_unified_sys_path(config: dict[str, Any]) -> Path:
    meta_root = _unified_paths(config)["metaboclip_project_dir"]
    if not meta_root.exists():
        raise FileNotFoundError(f"Unified MetaboClip backend directory not found: {meta_root}")
    if str(meta_root) not in sys.path:
        sys.path.insert(0, str(meta_root))
    return meta_root


def _detect_unified_backend_api(config: dict[str, Any]) -> tuple[Any, Any]:
    meta_root = _ensure_unified_sys_path(config)
    try:
        from metaboclip.core.workflow import run_directory, run_single_pair  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "Could not import unified MetaboClip Python APIs "
            "metaboclip.core.workflow.run_directory and run_single_pair "
            f"from {meta_root}. Ensure metaboclip_unified is present "
            "and its Python dependencies are installed."
        ) from exc
    return run_directory, run_single_pair


def _detect_ligand_role_api(config: dict[str, Any]) -> dict[str, Any]:
    meta_root = _ensure_unified_sys_path(config)
    try:
        from metaboclip_ligand_roles.annotator_core import annotate_ligand, write_annotation  # type: ignore
        from metaboclip_ligand_roles.atom_map import recover_atom_map_by_coordinates, write_atom_map  # type: ignore
        from metaboclip_ligand_roles.role_table import annotation_to_role_rows, write_role_table  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "Could not import unified MetaboClip ligand role APIs from "
            f"{meta_root}. Role-table generation requires RDKit, NumPy, SciPy, and PyYAML."
        ) from exc
    return {
        "annotate_ligand": annotate_ligand,
        "write_annotation": write_annotation,
        "recover_atom_map_by_coordinates": recover_atom_map_by_coordinates,
        "write_atom_map": write_atom_map,
        "annotation_to_role_rows": annotation_to_role_rows,
        "write_role_table": write_role_table,
    }


def _families_from_config(config: dict[str, Any]) -> list[str]:
    assign = config.get("family_assignment", {})
    mode = str(assign.get("mode", "run_all")).lower()
    if mode == "fixed":
        values = assign.get("fixed_families", DEFAULT_FAMILIES)
    else:
        values = assign.get("all_families", DEFAULT_FAMILIES)
    return [safe_id(str(f).lower()) for f in values]


def _validate_mechanism_paths(config: dict[str, Any], families: list[str]) -> dict[str, Path]:
    mechanisms = config.get("mechanisms", {}) or {}
    missing: list[str] = []
    out: dict[str, Path] = {}
    root = _root(config)
    for family in sorted(set(families)):
        raw = mechanisms.get(family)
        if not raw:
            missing.append(f"{family}: no mechanism path configured under mechanisms.{family}")
            continue
        path = resolve_path(raw, root)
        if path is None or not path.exists():
            missing.append(f"{family}: {path}")
            continue
        out[family] = path
    if missing:
        raise FileNotFoundError(
            "Required unified MetaboClip mechanism YAML files are missing:\n"
            + "\n".join(f"- {item}" for item in missing)
        )
    return out


def _role_rules_path(config: dict[str, Any]) -> Path:
    role_cfg = config.get("role_tables", {}) or {}
    raw = role_cfg.get("rules", "metaboclip_unified/rules/functional_groups.yaml")
    path = resolve_path(raw, _root(config))
    assert path is not None
    return path


def _role_groups(config: dict[str, Any]) -> list[str] | None:
    groups = config.get("role_tables", {}).get("groups")
    if groups in {None, ""}:
        return None
    if isinstance(groups, str):
        return [g.strip() for g in groups.replace(";", ",").split(",") if g.strip()]
    return [str(g).strip() for g in groups if str(g).strip()]


def validate_unified_bridge_config(config: dict[str, Any], staged_rows: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    _require_unified_backend(config)
    paths = _unified_paths(config)
    _detect_unified_backend_api(config)
    if not paths["metaboclip_profile"].exists():
        raise FileNotFoundError(f"Unified MetaboClip profile YAML not found: {paths['metaboclip_profile']}")
    role_mode = str((config.get("role_tables", {}) or {}).get("mode", "existing")).strip().lower()
    if role_mode not in {"existing", "generate", "auto"}:
        raise ValueError(f"Unsupported role_tables.mode: {role_mode}. Use existing, generate, or auto.")
    if role_mode in {"generate", "auto"} and not _role_rules_path(config).exists():
        raise FileNotFoundError(f"Unified ligand role rules YAML not found: {_role_rules_path(config)}")
    families = [row["family"] for row in staged_rows] if staged_rows else _families_from_config(config)
    mechanisms = _validate_mechanism_paths(config, families)
    return {
        "backend": _backend(config),
        "metaboclip_project_dir": str(paths["metaboclip_project_dir"]),
        "metaboclip_profile": str(paths["metaboclip_profile"]),
        "mechanisms": {family: str(path) for family, path in mechanisms.items()},
    }


def _truthy_status(status: str) -> bool:
    s = str(status or "success").strip().lower()
    return s in {"", "success", "skipped", "ok", "true", "1"}


def _row_path(row: dict[str, Any], key: str, root: Path) -> Path | None:
    raw = str(row.get(key, "") or "").strip()
    return resolve_path(raw, root) if raw else None


def load_refined_docking_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    root = _root(config)
    manifest = resolve_path(config.get("paths", {}).get("refined_docking_manifest", "data/data_output/refined/docking_out/docking_result_manifest.csv"), root)
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
        pose = _row_path(row, "out_pose_path", root)
        log = _row_path(row, "log_path", root)
        receptor = _row_path(row, "receptor_pdbqt_path", root)
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
    root = _root(config)
    map_path = resolve_path(config.get("family_assignment", {}).get("cluster_family_map_csv", "data/data_input/workflow/cluster_family_map.csv"), root)
    mapping: dict[str, list[str]] = {}
    if map_path and map_path.exists():
        for row in read_csv(map_path):
            cid = str(row.get("cluster_id") or row.get("cluster") or "").strip()
            family_raw = str(row.get("family") or row.get("families") or "").strip()
            if not cid or not family_raw:
                continue
            fams = [safe_id(f.strip().lower()) for f in family_raw.replace(";", ",").split(",") if f.strip()]
            mapping[cid] = fams
    return mapping


def families_for_row(config: dict[str, Any], row: dict[str, Any], cluster_map: dict[str, list[str]]) -> list[str]:
    assign = config.get("family_assignment", {})
    mode = str(assign.get("mode", "run_all")).lower()
    if mode == "run_all":
        return [safe_id(str(f).lower()) for f in assign.get("all_families", DEFAULT_FAMILIES)]
    if mode == "fixed":
        return [safe_id(str(f).lower()) for f in assign.get("fixed_families", DEFAULT_FAMILIES)]
    if mode == "cluster_family_map":
        fams = cluster_map.get(str(row.get("cluster_id", "")), [])
        if fams:
            return fams
        if bool(assign.get("fallback_to_all_when_unmapped", False)):
            return [safe_id(str(f).lower()) for f in assign.get("all_families", DEFAULT_FAMILIES)]
        return []
    if mode == "column":
        col = str(assign.get("family_column", "family"))
        raw = str(row.get(col, "")).strip()
        return [safe_id(f.strip().lower()) for f in raw.replace(";", ",").split(",") if f.strip()]
    raise ValueError(f"Unsupported family_assignment.mode: {mode}")


def _protein_name(row: dict[str, Any]) -> str:
    protein = safe_id(row.get("protein_id", "protein"))
    pocket = safe_id(row.get("pocket_id", "P1"))
    conformer = safe_id(row.get("conformer_id", "conf_0"))
    match = re.search(r"(\d+)$", conformer)
    if match:
        base = re.sub(r"[_-]?\d+$", "", conformer).strip("_-") or "conf"
        return f"{protein}{PROTEIN_DELIMITER}{pocket}__{base}_{match.group(1)}"
    return f"{protein}{PROTEIN_DELIMITER}{pocket}__{conformer}"


def _split_unified_protein_name(protein_name: str) -> tuple[str, str, int]:
    match = re.match(r"^(?P<base>.+)_(?P<idx>\d+)$", str(protein_name))
    if match:
        idx = int(match.group("idx"))
        return match.group("base"), str(idx), idx
    return str(protein_name), "0", 0


def _unified_safe_name(value: Any) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))


def stage_metaboclip_inputs(config: dict[str, Any]) -> tuple[Path, list[dict[str, Any]]]:
    root = _root(config)
    staging_root = resolve_path(config.get("paths", {}).get("staging_dir", "data/data_output/metaboclip/staging"), root)
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
        unified_protein_id, unified_conf_id, unified_conf_index = _split_unified_protein_name(protein_name)
        pose = _row_path(row, "out_pose_path", root)
        log = _row_path(row, "log_path", root)
        receptor = _row_path(row, "receptor_pdbqt_path", root)
        if pose is None or log is None or receptor is None:
            raise FileNotFoundError(f"Missing refined docking pose, log, or receptor path for job_id={row.get('job_id', '')}")
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
            staged_row: dict[str, Any] = {
                "family": family,
                "ligand_id": ligand_id,
                "original_ligand_id": ligand_id_raw,
                "protein_id": row.get("protein_id", ""),
                "cluster_id": row.get("cluster_id", ""),
                "batch_id": row.get("batch_id", ""),
                "protein_name": protein_name,
                "unified_protein_id": unified_protein_id,
                "unified_conformation_id": unified_conf_id,
                "unified_conformation_index": unified_conf_index,
                "pocket_id": row.get("pocket_id", ""),
                "pocket_rank": row.get("pocket_rank", ""),
                "conformer_id": row.get("conformer_id", ""),
                "pose_id": row.get("pose_id") or row.get("conformer_id", ""),
                "job_id": row.get("job_id", ""),
                "best_affinity": row.get("best_affinity", ""),
                "affinity_kcal": row.get("best_affinity", ""),
                "protein_dir": str(protein_dir),
                "docking_dir": str(docking_dir),
                "staged_receptor_path": str(receptor_dst),
                "staged_pose_path": str(ligand_dst),
                "staged_log_path": str(log_dst),
                "receptor_pdbqt_path": str(receptor),
                "ligand_pdbqt_path": row.get("ligand_pdbqt_path", ""),
                "out_pose_path": str(pose),
                "log_path": str(log),
                "source_receptor_pdbqt_path": str(receptor),
                "source_out_pose_path": str(pose),
                "source_log_path": str(log),
                "source_status": row.get("status", ""),
                "source_message": row.get("message", ""),
            }
            for field in AIMD_DOCKING_METADATA_FIELDS:
                if field not in staged_row:
                    staged_row[field] = row.get(field, "")
            staged.append(staged_row)
    if not staged:
        raise RuntimeError("No rows were staged for MetaBoClip. Check family_assignment and refined docking results.")
    manifest_path = staging_root / "metaboclip_staging_manifest.csv"
    write_csv(manifest_path, staged)
    return manifest_path, staged


def _load_ligand_lookup(config: dict[str, Any]) -> dict[str, dict[str, str]]:
    paths = _unified_paths(config)
    lookup: dict[str, dict[str, str]] = {}
    for manifest in [paths["ligand_manifest"], paths["ligand_source_manifest"]]:
        if not manifest.exists():
            continue
        for row in read_csv(manifest):
            ligand_id = str(row.get("ligand_id") or "").strip()
            if not ligand_id:
                candidate = row.get("ligand_path") or row.get("ligand_source_path") or row.get("source_path")
                ligand_id = Path(candidate).stem if candidate else ""
            if ligand_id:
                lookup.setdefault(ligand_id, {}).update(row)
                lookup.setdefault(safe_id(ligand_id), {}).update(row)
    return lookup


def _first_existing_value(row: dict[str, Any], keys: list[str]) -> str:
    for key in keys:
        value = str(row.get(key, "") or "").strip()
        if value:
            return value
    return ""


def _resolve_optional_path(raw: str, root: Path) -> Path | None:
    if not raw:
        return None
    return resolve_path(raw, root)


def _ligand_source_path(config: dict[str, Any], staged_row: dict[str, Any], ligand_lookup: dict[str, dict[str, str]]) -> Path | None:
    root = _root(config)
    role_cfg = config.get("role_tables", {}) or {}
    ligand_id = str(staged_row.get("original_ligand_id") or staged_row.get("ligand_id") or "")
    lookup_row = ligand_lookup.get(ligand_id) or ligand_lookup.get(safe_id(ligand_id)) or {}
    column = str(role_cfg.get("ligand_source_column", "ligand_source_path"))
    raw = _first_existing_value(
        staged_row,
        [column, "ligand_source_path", "source_ligand_path", "source_path", "sdf_path", "mol_path", "mol2_path"],
    ) or _first_existing_value(
        lookup_row,
        [column, "ligand_source_path", "source_ligand_path", "source_path", "sdf_path", "mol_path", "mol2_path"],
    )
    return _resolve_optional_path(raw, root)


def _prepared_ligand_pdbqt_path(config: dict[str, Any], staged_row: dict[str, Any], ligand_lookup: dict[str, dict[str, str]]) -> Path | None:
    root = _root(config)
    role_cfg = config.get("role_tables", {}) or {}
    ligand_id = str(staged_row.get("original_ligand_id") or staged_row.get("ligand_id") or "")
    lookup_row = ligand_lookup.get(ligand_id) or ligand_lookup.get(safe_id(ligand_id)) or {}
    column = str(role_cfg.get("prepared_pdbqt_column", "ligand_pdbqt_path"))
    raw = _first_existing_value(
        staged_row,
        [column, "ligand_pdbqt_path", "prepared_ligand_pdbqt_path"],
    ) or _first_existing_value(
        lookup_row,
        [column, "ligand_pdbqt_path", "prepared_ligand_pdbqt_path", "ligand_path", "pdbqt_path"],
    )
    return _resolve_optional_path(raw, root)


def _role_paths_for_ligand(config: dict[str, Any], ligand_id: str) -> dict[str, Path]:
    paths = _unified_paths(config)
    clean_id = safe_id(ligand_id)
    return {
        "role_table": paths["role_table_dir"] / f"{clean_id}.role_table.csv",
        "annotation": paths["annotation_dir"] / f"{clean_id}.annotation.json",
        "atom_map": paths["atom_map_dir"] / f"{clean_id}.atom_map.json",
    }


def _generate_role_table(
    config: dict[str, Any],
    ligand_id: str,
    ligand_source: Path,
    prepared_pdbqt: Path,
    role_table: Path,
    annotation: Path,
    atom_map_path: Path,
) -> dict[str, Any]:
    api = _detect_ligand_role_api(config)
    rules = _role_rules_path(config)
    role_cfg = config.get("role_tables", {}) or {}
    max_dist = float(role_cfg.get("max_dist", 0.5))
    groups = _role_groups(config)
    annotation_obj = api["annotate_ligand"](ligand_source, ligand_id, rules, groups=groups)
    atom_map = api["recover_atom_map_by_coordinates"](ligand_source, prepared_pdbqt, ligand_id, max_dist=max_dist)
    role_rows = api["annotation_to_role_rows"](annotation_obj, atom_map=atom_map)
    api["write_role_table"](role_rows, role_table)
    api["write_annotation"](annotation_obj, annotation)
    api["write_atom_map"](atom_map, atom_map_path)
    return {
        "generated_role_rows": len(role_rows),
        "matched_heavy_atoms": atom_map.get("matched_atoms", ""),
        "atom_map_rmsd": atom_map.get("rmsd", ""),
    }


def _resolve_role_assets(
    config: dict[str, Any],
    staged_row: dict[str, Any],
    ligand_lookup: dict[str, dict[str, str]],
    allow_generation: bool,
) -> dict[str, str]:
    role_cfg = config.get("role_tables", {}) or {}
    mode = str(role_cfg.get("mode", "existing")).strip().lower()
    if mode not in {"existing", "generate", "auto"}:
        raise ValueError(f"Unsupported role_tables.mode: {mode}. Use existing, generate, or auto.")

    ligand_id = safe_id(staged_row.get("ligand_id", ""))
    role_paths = _role_paths_for_ligand(config, ligand_id)
    role_table = role_paths["role_table"]
    annotation = role_paths["annotation"]
    atom_map = role_paths["atom_map"]
    role_table_exists = role_table.exists()

    needs_generation = mode == "generate" or (mode == "auto" and not role_table_exists)
    if mode == "existing" and not role_table_exists:
        raise FileNotFoundError(
            f"Role table missing for ligand_id={ligand_id}: {role_table}. "
            "Provide the role table or set role_tables.mode to generate/auto with ligand source inputs."
        )

    ligand_source = _ligand_source_path(config, staged_row, ligand_lookup)
    prepared_pdbqt = _prepared_ligand_pdbqt_path(config, staged_row, ligand_lookup)
    generation_report: dict[str, Any] = {}
    if needs_generation:
        rules = _role_rules_path(config)
        if not rules.exists():
            raise FileNotFoundError(f"Functional-group rules YAML is missing: {rules}")
        if ligand_source is None:
            raise FileNotFoundError(
                f"Ligand source path missing for ligand_id={ligand_id}. "
                f"role_tables.mode={mode} requires a source path configured via "
                "paths.ligand_source_manifest or paths.ligand_manifest."
            )
        if not ligand_source.exists():
            raise FileNotFoundError(f"Ligand source path does not exist for ligand_id={ligand_id}: {ligand_source}")
        if prepared_pdbqt is None:
            raise FileNotFoundError(
                f"Prepared ligand PDBQT path missing for ligand_id={ligand_id}. "
                "Use ligand_pdbqt_path in the refined docking manifest or configure it in the ligand manifest."
            )
        if not prepared_pdbqt.exists():
            raise FileNotFoundError(f"Prepared ligand PDBQT path does not exist for ligand_id={ligand_id}: {prepared_pdbqt}")
        if allow_generation:
            generation_report = _generate_role_table(
                config,
                ligand_id,
                ligand_source,
                prepared_pdbqt,
                role_table,
                annotation,
                atom_map,
            )
            role_table_exists = role_table.exists()

    if not role_table_exists and allow_generation:
        raise FileNotFoundError(f"Role table was not created for ligand_id={ligand_id}: {role_table}")

    role_status = "existing" if role_table_exists and not needs_generation else "needs_generation"
    if allow_generation and needs_generation and role_table_exists:
        role_status = "generated"

    return {
        "role_table_path": str(role_table),
        "annotation_json_path": str(annotation),
        "atom_map_json_path": str(atom_map),
        "ligand_source_path": str(ligand_source or ""),
        "prepared_ligand_pdbqt_path": str(prepared_pdbqt or ""),
        "role_table_status": role_status,
        "annotation_status": "available" if annotation.exists() else "missing",
        "atom_map_status": "available" if atom_map.exists() else "missing",
        "generated_role_rows": str(generation_report.get("generated_role_rows", "")),
        "matched_heavy_atoms": str(generation_report.get("matched_heavy_atoms", "")),
        "atom_map_rmsd": str(generation_report.get("atom_map_rmsd", "")),
    }


def _base_run_row(
    config: dict[str, Any],
    row: dict[str, Any],
    mechanisms: dict[str, Path],
    role_assets: dict[str, str],
    status: str,
    message: str,
) -> dict[str, Any]:
    paths = _unified_paths(config)
    family = row["family"]
    ligand_id = row["ligand_id"]
    out_dir = paths["unified_output_dir"] / family / ligand_id
    pair_out = out_dir / _unified_safe_name(ligand_id) / _unified_safe_name(row.get("protein_name", ""))
    run_row: dict[str, Any] = {
        "backend": "unified",
        "family": family,
        "ligand_id": ligand_id,
        "original_ligand_id": row.get("original_ligand_id", ""),
        "protein_id": row.get("protein_id", ""),
        "cluster_id": row.get("cluster_id", ""),
        "batch_id": row.get("batch_id", ""),
        "protein_name": row.get("protein_name", ""),
        "unified_protein_id": row.get("unified_protein_id", ""),
        "unified_conformation_id": row.get("unified_conformation_id", ""),
        "unified_conformation_index": row.get("unified_conformation_index", ""),
        "pocket_id": row.get("pocket_id", ""),
        "pocket_rank": row.get("pocket_rank", ""),
        "conformer_id": row.get("conformer_id", ""),
        "pose_id": row.get("pose_id") or row.get("conformer_id", ""),
        "job_id": row.get("job_id", ""),
        "best_affinity": row.get("best_affinity", ""),
        "affinity_kcal": row.get("best_affinity", ""),
        "protein_dir": row.get("protein_dir", ""),
        "docking_dir": row.get("docking_dir", ""),
        "out_dir": str(out_dir),
        "pair_output_dir": str(pair_out),
        "mechanism_path": str(mechanisms.get(family, "")),
        "profile_path": str(paths["metaboclip_profile"]),
        "receptor_pdbqt_path": row.get("receptor_pdbqt_path", ""),
        "ligand_pdbqt_path": row.get("ligand_pdbqt_path", ""),
        "out_pose_path": row.get("out_pose_path", ""),
        "log_path": row.get("log_path", ""),
        "source_receptor_pdbqt_path": row.get("source_receptor_pdbqt_path", ""),
        "source_out_pose_path": row.get("source_out_pose_path", ""),
        "source_log_path": row.get("source_log_path", ""),
        "staged_receptor_path": row.get("staged_receptor_path", ""),
        "staged_pose_path": row.get("staged_pose_path", ""),
        "staged_log_path": row.get("staged_log_path", ""),
        "source_status": row.get("source_status", ""),
        "source_message": row.get("source_message", ""),
        **role_assets,
        "status": status,
        "message": message,
    }
    for field in AIMD_DOCKING_METADATA_FIELDS:
        if field not in run_row:
            run_row[field] = row.get(field, "")
    return run_row


def prepare_unified_metaboclip_runs(
    config: dict[str, Any],
    staged_rows: list[dict[str, Any]],
    allow_generation: bool = False,
) -> list[dict[str, Any]]:
    root = _root(config)
    results_root = resolve_path(config.get("paths", {}).get("results_dir", "data/data_output/metaboclip/results"), root)
    assert results_root is not None
    validation = validate_unified_bridge_config(config, staged_rows)
    mechanisms = {family: Path(path) for family, path in validation["mechanisms"].items()}
    ligand_lookup = _load_ligand_lookup(config)
    role_cache: dict[str, dict[str, str]] = {}
    run_rows: list[dict[str, Any]] = []
    for row in staged_rows:
        ligand_key = str(row.get("ligand_id", ""))
        try:
            if ligand_key not in role_cache:
                role_cache[ligand_key] = _resolve_role_assets(config, row, ligand_lookup, allow_generation=allow_generation)
            role_assets = role_cache[ligand_key]
            status = "prepared"
            message = PREPARE_ONLY_MESSAGE if not allow_generation else "Unified backend inputs are ready"
        except Exception as exc:
            role_assets = {
                "role_table_path": "",
                "annotation_json_path": "",
                "atom_map_json_path": "",
                "ligand_source_path": "",
                "prepared_ligand_pdbqt_path": "",
                "role_table_status": "",
                "annotation_status": "",
                "atom_map_status": "",
                "generated_role_rows": "",
                "matched_heavy_atoms": "",
                "atom_map_rmsd": "",
            }
            status = "failed"
            message = str(exc)
            if not bool(config.get("execution", {}).get("continue_on_error", True)):
                raise
        run_rows.append(_base_run_row(config, row, mechanisms, role_assets, status, message))
    write_csv(results_root / "metaboclip_run_manifest.csv", run_rows)
    return run_rows


def _group_ready_rows(run_rows: list[dict[str, Any]]) -> dict[tuple[str, str, str, str, str, str], list[dict[str, Any]]]:
    groups: dict[tuple[str, str, str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in run_rows:
        if row.get("status") == "failed":
            continue
        role_table = Path(str(row.get("role_table_path", "")))
        if not role_table.exists():
            row["status"] = "failed"
            row["message"] = f"Role table is not available for unified scoring: {role_table}"
            continue
        key = (
            str(row.get("family", "")),
            str(row.get("ligand_id", "")),
            str(row.get("protein_dir", "")),
            str(row.get("docking_dir", "")),
            str(role_table.parent),
            str(row.get("out_dir", "")),
        )
        groups[key].append(row)
    return groups


def execute_unified_metaboclip_runs(config: dict[str, Any], staged_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    run_directory, _run_single_pair = _detect_unified_backend_api(config)
    run_rows = prepare_unified_metaboclip_runs(config, staged_rows, allow_generation=True)
    continue_on_error = bool(config.get("execution", {}).get("continue_on_error", True))
    groups = _group_ready_rows(run_rows)
    for key, rows in groups.items():
        family, ligand_id, protein_dir, docking_dir, role_table_dir, out_dir = key
        mechanism_path = Path(str(rows[0].get("mechanism_path", "")))
        profile_path = Path(str(rows[0].get("profile_path", "")))
        try:
            result = run_directory(
                mechanism_path,
                profile_path,
                Path(protein_dir),
                Path(docking_dir),
                Path(role_table_dir),
                Path(out_dir),
                ligand_id,
            )
            pairs = int(result.get("pairs", 0))
            if pairs <= 0:
                raise RuntimeError(
                    "Unified backend found no matching protein, docking pose, and role-table pairs "
                    f"for family={family}, ligand_id={ligand_id}."
                )
            for row in rows:
                pair_out = Path(str(row.get("pair_output_dir", "")))
                row.update({
                    "status": "success",
                    "message": "Unified backend scoring completed",
                    "unified_pairs": result.get("pairs", ""),
                    "unified_candidates": result.get("candidates", ""),
                    "unified_proteins": result.get("proteins", ""),
                    "protein_scores_path": str(Path(out_dir) / "protein_scores.csv"),
                    "merged_conformation_scores_path": str(Path(out_dir) / "merged_conformation_scores.csv"),
                    "pose_scores_path": str(pair_out / "pose_scores.csv"),
                    "candidate_scores_path": str(pair_out / "candidate_scores.csv"),
                    "passing_candidates_path": str(pair_out / "passing_candidates.csv"),
                    "geometry_features_path": str(pair_out / "geometry_features.csv"),
                    "resolved_ligand_sites_path": str(pair_out / "resolved_ligand_sites.csv"),
                    "resolved_protein_roles_path": str(pair_out / "resolved_protein_roles.csv"),
                })
        except Exception as exc:
            for row in rows:
                row["status"] = "failed"
                row["message"] = str(exc)
            if not continue_on_error:
                raise
    results_root = resolve_path(config.get("paths", {}).get("results_dir", "data/data_output/metaboclip/results"), _root(config))
    assert results_root is not None
    write_csv(results_root / "metaboclip_run_manifest.csv", run_rows)
    return run_rows


def run_metaboclip_family_scoring(config: dict[str, Any], staged_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    _require_unified_backend(config)
    return execute_unified_metaboclip_runs(config, staged_rows)


def _metadata_for_run(row: dict[str, Any]) -> dict[str, Any]:
    metadata: dict[str, Any] = {}
    for field in RUN_METADATA_FIELDS:
        metadata[field] = row.get(field, "")
    for field in AIMD_DOCKING_METADATA_FIELDS:
        if field not in metadata:
            metadata[field] = row.get(field, "")
    metadata["source_status"] = row.get("source_status", "")
    metadata["source_message"] = row.get("source_message", "")
    return metadata


def _merge_metadata(unified_row: dict[str, Any], run_row: dict[str, Any]) -> dict[str, Any]:
    out = _metadata_for_run(run_row)
    for key, value in unified_row.items():
        if key in out and str(out.get(key, "")) not in {"", str(value)}:
            out[f"unified_{key}"] = value
        elif key in {"protein_id", "ligand_id", "pose_id", "conformation_id", "conformation_name", "conformation_index"}:
            out[f"unified_{key}"] = value
            out.setdefault(key, value)
        else:
            out[key] = value
    return out


def _find_run_for_unified_row(unified_row: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    conformation_name = str(unified_row.get("conformation_name", "") or "")
    protein_id = str(unified_row.get("protein_id", "") or "")
    ligand_id = str(unified_row.get("ligand_id", "") or "")
    for row in rows:
        if conformation_name and conformation_name == str(row.get("protein_name", "")):
            return row
    for row in rows:
        if protein_id and protein_id == str(row.get("unified_protein_id", "")):
            if not ligand_id or ligand_id == str(row.get("ligand_id", "")):
                return row
    for row in rows:
        if not ligand_id or ligand_id == str(row.get("ligand_id", "")):
            return row
    return rows[0]


def _successful_groups(run_rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in run_rows:
        if row.get("status") != "success":
            continue
        out_dir = str(row.get("out_dir", ""))
        if out_dir:
            groups[out_dir].append(row)
    return groups


def _write_aggregate(results_root: Path, filename: str, rows: list[dict[str, Any]], base_fields: list[str]) -> Path:
    return write_csv(results_root / filename, rows, fieldnames=None if rows else base_fields)


def _sort_final_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    def num(row: dict[str, Any], key: str) -> float:
        value = numeric(row.get(key), None)
        return float(value) if value is not None else float("-inf")

    return sorted(
        rows,
        key=lambda row: (
            -num(row, "protein_score"),
            -num(row, "quality_score"),
            -num(row, "coverage"),
            str(row.get("family", "")),
            str(row.get("ligand_id", "")),
            str(row.get("protein_id", "")),
        ),
    )


def collect_metaboclip_outputs(config: dict[str, Any], run_rows: list[dict[str, Any]]) -> Path:
    _require_unified_backend(config)
    root = _root(config)
    results_root = resolve_path(config.get("paths", {}).get("results_dir", "data/data_output/metaboclip/results"), root)
    assert results_root is not None
    write_csv(results_root / "metaboclip_run_manifest.csv", run_rows)

    aggregate_rows: dict[str, list[dict[str, Any]]] = {name: [] for name in ROOT_OUTPUT_SPECS | PAIR_OUTPUT_SPECS}
    output_paths: dict[str, str] = {}

    for out_dir, rows in _successful_groups(run_rows).items():
        for name, (source_name, output_name) in ROOT_OUTPUT_SPECS.items():
            source = Path(out_dir) / source_name
            if not source.exists():
                continue
            for unified_row in read_csv(source):
                run_row = _find_run_for_unified_row(unified_row, rows)
                aggregate_rows[name].append(_merge_metadata(unified_row, run_row))
            output_paths[name] = str(results_root / output_name)

    seen_pair_outputs: set[tuple[str, str]] = set()
    for row in run_rows:
        if row.get("status") != "success":
            continue
        pair_dir = Path(str(row.get("pair_output_dir", "")))
        for name, (source_name, output_name) in PAIR_OUTPUT_SPECS.items():
            key = (name, str(pair_dir))
            if key in seen_pair_outputs:
                continue
            seen_pair_outputs.add(key)
            source = pair_dir / source_name
            if not source.exists():
                continue
            for unified_row in read_csv(source):
                aggregate_rows[name].append(_merge_metadata(unified_row, row))
            output_paths[name] = str(results_root / output_name)

    base_fields = RUN_METADATA_FIELDS + ["unified_protein_id", "unified_ligand_id", "protein_score", "quality_score", "coverage"]
    for name, (_source_name, output_name) in ROOT_OUTPUT_SPECS.items():
        _write_aggregate(results_root, output_name, aggregate_rows[name], base_fields)
    for name, (_source_name, output_name) in PAIR_OUTPUT_SPECS.items():
        _write_aggregate(results_root, output_name, aggregate_rows[name], base_fields)

    final_rows = aggregate_rows["protein_scores"]
    if final_rows:
        final_rows = _sort_final_rows(final_rows)
    elif aggregate_rows["merged_conformation_scores"]:
        final_rows = list(aggregate_rows["merged_conformation_scores"])
    elif aggregate_rows["pose_scores"]:
        final_rows = list(aggregate_rows["pose_scores"])
    elif aggregate_rows["candidate_scores"]:
        final_rows = list(aggregate_rows["candidate_scores"])
    else:
        final_rows = list(run_rows)

    for rank, row in enumerate(final_rows, start=1):
        row["rank"] = rank
    final_path = results_root / "metaboclip_final_ranking.csv"
    write_csv(final_path, final_rows)

    scoring_executed = any(row.get("status") == "success" for row in run_rows)
    write_json(results_root / "metaboclip_report.json", {
        "backend": "unified",
        "scoring_executed": scoring_executed,
        "n_run_rows": len(run_rows),
        "n_success_rows": sum(1 for row in run_rows if row.get("status") == "success"),
        "n_failed_rows": sum(1 for row in run_rows if row.get("status") == "failed"),
        "n_final_rows": len(final_rows),
        "final_ranking": str(final_path),
        "outputs": output_paths,
        "legacy_columns_not_generated": ["protein_score_norm", "max_s_r"],
    })
    return final_path


def run_metaboclip_bridge_core(config: dict[str, Any]) -> Path:
    _require_unified_backend(config)
    staging_manifest, staged = stage_metaboclip_inputs(config)
    print(f"[MetaBoClipBridge] staged inputs: {staging_manifest} ({len(staged)} rows)")
    if bool(config.get("execution", {}).get("run_metaboclip", True)):
        run_rows = run_metaboclip_family_scoring(config, staged)
    else:
        run_rows = prepare_unified_metaboclip_runs(config, staged, allow_generation=False)
        for row in run_rows:
            if row.get("status") != "failed":
                row["status"] = "not_run"
                row["message"] = "execution.run_metaboclip=false"
        results_root = resolve_path(config.get("paths", {}).get("results_dir", "data/data_output/metaboclip/results"), _root(config))
        assert results_root is not None
        write_csv(results_root / "metaboclip_run_manifest.csv", run_rows)
    final_path = collect_metaboclip_outputs(config, run_rows)
    print(f"[MetaBoClipBridge] final ranking: {final_path}")
    return final_path
