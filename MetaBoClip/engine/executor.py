from __future__ import annotations

from pathlib import Path
from typing import Any

from core.output import write_rows_csv, write_score_tables
from core.pymol_utils import get_pymol, parse_vina_out, save_filtered_pse, split_ligand_states
from engine.operators import (
    attr_of,
    compute_feature,
    execute_object_op,
    format_selector,
    is_ref_object,
    is_missing,
    resolve_name,
)
from engine.scoring import score_rows


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _resolve_resource_path(raw_path: str, spec_path: Path | None) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    candidates = []
    if spec_path is not None:
        candidates.append(spec_path.parent / path)
        candidates.append(spec_path.parent.parent / path)
    candidates.append(PROJECT_ROOT / path)
    for cand in candidates:
        if cand.exists():
            return cand
    return candidates[0]


def _load_resources(cmd, cfg: dict, selector_vars: dict[str, Any], spec_path: Path | None) -> bool:
    for key, spec in (cfg.get("resources", {}) or {}).items():
        path = _resolve_resource_path(spec["path"], spec_path)
        obj_name = spec.get("object", key)
        if not path.exists():
            if spec.get("required", True):
                return False
            continue
        cmd.load(str(path), obj_name)
        selector_vars[f"{key}_obj"] = obj_name

        # Optional resource-to-protein alignment. This is used by ACH to make
        # the dAC-T template Ser171 anchor meaningful even when the loaded
        # protein is not already in the template coordinate frame. The ligand
        # is loaded later and remains in the protein coordinate frame, so only
        # the template/resource is moved.
        if spec.get("align_to_protein", False):
            mobile_tpl = spec.get("mobile_selector", f"{{{key}_obj}} and name CA")
            target_tpl = spec.get("target_selector", "{protein_obj} and name CA")
            mobile = format_selector(mobile_tpl, selector_vars)
            target = format_selector(target_tpl, selector_vars)
            alignment_required = bool(spec.get("alignment_required", True))

            try:
                if cmd.count_atoms(mobile) <= 0 or cmd.count_atoms(target) <= 0:
                    if alignment_required:
                        return False
                    continue

                method = str(spec.get("align_method", "super")).lower()
                if method == "align":
                    cmd.align(mobile, target)
                elif method == "cealign":
                    # PyMOL cealign takes target first, mobile second.
                    cmd.cealign(target, mobile)
                else:
                    # super is the robust default for homologous proteins.
                    cmd.super(mobile, target)
            except Exception:
                if alignment_required:
                    return False
    return True


def _resolve_named_objects(cmd, defs: dict[str, Any], named: dict[str, Any], selector_vars: dict[str, Any]) -> bool:
    for name, spec in (defs or {}).items():
        value = execute_object_op(cmd, spec["op"], spec, named, selector_vars)
        named[name] = value
        if spec.get("required", True) and is_missing(value):
            return False
    return True


def _candidate_rows(candidate_spec: dict[str, Any], named: dict[str, Any]) -> list[dict[str, Any]]:
    if not candidate_spec:
        return [{}]
    op = candidate_spec.get("op", "source")
    if op == "source":
        collection = resolve_name(candidate_spec["collection"], named)
        alias = candidate_spec.get("as", "item")
        out = []
        for item in collection:
            if is_ref_object(item):
                out.append({alias: item})
            elif isinstance(item, dict):
                out.append(dict(item))
            else:
                out.append({alias: item})
        return out
    if op == "cartesian_product":
        import itertools

        members = candidate_spec.get("members", {})
        buckets = []
        for alias, ref_name in members.items():
            buckets.append((alias, list(resolve_name(ref_name, named))))
        out = []
        for combo in itertools.product(*[vals for _, vals in buckets]):
            row = {}
            for (alias, _), value in zip(buckets, combo):
                if is_ref_object(value):
                    row[alias] = value
                elif isinstance(value, dict):
                    row.update(dict(value))
                else:
                    row[alias] = value
            out.append(row)
        return out
    raise RuntimeError(f"Unknown candidate operator: {op}")


def _rule_passes(features: dict[str, Any], rule: dict[str, Any]) -> bool:
    value = features[rule.get("feature") or rule.get("source")]
    op = rule["op"]
    if op == "between":
        return float(rule["min"]) <= float(value) <= float(rule["max"])
    if op == "<=":
        return float(value) <= float(rule["value"])
    if op == "<":
        return float(value) < float(rule["value"])
    if op == ">=":
        return float(value) >= float(rule["value"])
    if op == ">":
        return float(value) > float(rule["value"])
    if op == "==":
        return value == rule["value"]
    if op == "!=":
        return value != rule["value"]
    raise RuntimeError(f"Unknown gating operator: {op}")


def _passes_gating(features: dict[str, Any], gating_cfg: dict[str, Any]) -> bool:
    rules = gating_cfg.get("rules", [])
    if not rules:
        return True
    logic = gating_cfg.get("logic", "all")
    outcomes = [_rule_passes(features, rule) for rule in rules]
    return all(outcomes) if logic == "all" else any(outcomes)


def _output_fields(field_specs: list[dict[str, Any]], named: dict[str, Any]) -> dict[str, Any]:
    out = {}
    for spec in field_specs or []:
        name = spec["name"]
        if "value" in spec:
            out[name] = spec["value"]
            continue
        ref = resolve_name(spec["ref"], named)
        out[name] = attr_of(ref, spec["attr"]) if "attr" in spec else ref
    return out


def _persist_candidate_selection(cmd, named: dict[str, Any]) -> None:
    """Create the per-hit selection stored in selection_name for saved PSE sessions.

    The original ACH script creates selections such as ligand_0001_Csel_1 for
    passing carbonyl-carbon hits.  The YAML engine stores the same name in the
    output row; this helper materializes that selection only after the row has
    passed gating.
    """
    sel_name = named.get("selection_name")
    carbon = named.get("carbonyl_c")
    if not sel_name or not isinstance(carbon, dict):
        return
    obj_name = carbon.get("object")
    atom_index = carbon.get("pymol_index") or carbon.get("index")
    if not obj_name or atom_index is None:
        return
    try:
        cmd.select(str(sel_name), f"({obj_name} and index {int(atom_index)})")
    except Exception:
        pass


def _process_single_protein(cmd, cfg: dict, protein_path: Path, docking_dir: Path, ligand_id: str, out_dir: Path, spec_path: Path | None):
    rows = []
    passing_models = set()

    protein_name = protein_path.stem
    inp = cfg.get("input", {})
    ligand_ext = inp.get("ligand_extension", ".pdbqt")
    ligand_pattern = inp.get("ligand_pattern", "{ligand_id}@{protein_name}{ligand_extension}")
    out_pattern = inp.get("out_pattern", "{ligand_id}_{protein_name}_cavity_1.out")
    protein_obj = inp.get("protein_object", "protein")
    ligand_obj = inp.get("ligand_object", "ligand")

    ligand_file = docking_dir / ligand_pattern.format(ligand_id=ligand_id, protein_name=protein_name, ligand_extension=ligand_ext)
    out_file = docking_dir / out_pattern.format(ligand_id=ligand_id, protein_name=protein_name)
    if not protein_path.exists() or not ligand_file.exists() or not out_file.exists():
        return rows

    modes, mode2aff = parse_vina_out(str(out_file), float(inp.get("affinity_threshold", -3.0)))
    if not modes:
        return rows

    cmd.reinitialize()
    cmd.load(str(protein_path), protein_obj)
    selector_vars = {
        "protein_obj": protein_obj,
        "ligand_obj": ligand_obj,
        "protein_name": protein_name,
        "ligand_id": ligand_id,
    }
    if not _load_resources(cmd, cfg, selector_vars, spec_path):
        return rows

    global_named = {}
    if not _resolve_named_objects(cmd, cfg.get("globals", {}), global_named, selector_vars):
        return rows

    cmd.load(str(ligand_file), ligand_obj)
    all_models = split_ligand_states(cmd, ligand_obj, "ligand_")
    retained = [f"ligand_{str(mode).zfill(4)}" for mode in modes if f"ligand_{str(mode).zfill(4)}" in cmd.get_names("objects")]
    for model in all_models:
        if model not in retained:
            try:
                cmd.delete(model)
            except Exception:
                pass

    pose_cfg = cfg.get("pose", {}) or {}
    pose_defs = pose_cfg.get("objects", {}) or {}
    candidate_spec = cfg.get("candidates", {}) or {}
    feature_specs = cfg.get("features", []) or []
    output_cfg = cfg.get("output", {}) or {}
    field_specs = output_cfg.get("row_fields", []) or []

    for model in retained:
        selector_vars_pose = dict(selector_vars)
        selector_vars_pose["pose_obj"] = model
        if pose_cfg.get("add_hydrogens", False):
            try:
                cmd.h_add(model)
            except Exception:
                pass

        named = dict(global_named)
        if not _resolve_named_objects(cmd, pose_defs, named, selector_vars_pose):
            continue

        candidate_rows = _candidate_rows(candidate_spec, named)
        if not candidate_rows:
            continue

        mode_num = int(model.split("_")[1])
        affinity = mode2aff.get(mode_num, float("nan"))
        pose_any = False
        for candidate in candidate_rows:
            row_named = dict(named)
            row_named.update(candidate)
            features = {}
            for feature_spec in feature_specs:
                value = compute_feature(feature_spec["op"], feature_spec, row_named)
                features[feature_spec["name"]] = value
                row_named[feature_spec["name"]] = value
            if not _passes_gating(features, cfg.get("gating", {})):
                continue
            _persist_candidate_selection(cmd, row_named)
            delim = inp.get("protein_id_delimiter", None)
            if delim:
                base_protein_id = protein_name.split(str(delim))[0]
            elif "__" in protein_name:
                base_protein_id = protein_name.split("__")[0]
            else:
                base_protein_id = protein_name.split("_")[0]
            row = {
                "Ligand_id": str(ligand_id),
                "protein_id": base_protein_id,
                "conformation": protein_name,
                "mode": int(mode_num),
                "affinity_kcal": float(affinity),
            }
            row.update(features)
            row.update(_output_fields(field_specs, row_named))
            rows.append(row)
            pose_any = True
        if pose_any:
            passing_models.add(model)

    if passing_models and output_cfg.get("save_pse", True):
        sessions_dir = out_dir / "sessions"
        sessions_dir.mkdir(parents=True, exist_ok=True)
        save_filtered_pse(cmd, passing_models, sessions_dir / f"{protein_name}.pse")
    return rows


def run_family_spec(*, cfg: dict, protein_dir: Path, docking_dir: Path, out_dir: Path, spec_path: Path | None = None, registry: dict | None = None) -> dict[str, Any]:
    inp = cfg.get("input", {})
    ligand_id = str(inp.get("ligand_id") or Path(docking_dir).name.split("_")[-1])
    protein_ext = inp.get("protein_extension", ".pdbqt")
    proteins = sorted(Path(protein_dir).glob(f"*{protein_ext}"))
    rows = []
    cmd, pm = get_pymol()
    try:
        for protein_path in proteins:
            rows.extend(_process_single_protein(cmd, cfg, protein_path, Path(docking_dir), ligand_id, Path(out_dir), spec_path))
            try:
                cmd.delete("all")
            except Exception:
                pass
    finally:
        try:
            if pm is not None:
                pm.stop()
        except Exception:
            pass

    gating_csv = Path(out_dir) / f"file_{ligand_id}.gating.csv"
    write_rows_csv(rows, gating_csv)
    score_tables = score_rows(rows, cfg)
    write_score_tables(score_tables, Path(out_dir), f"file_{ligand_id}")
    return {
        "family": cfg.get("family", {}).get("key") or cfg.get("family", {}).get("plugin"),
        "rows": len(rows),
        "gating_csv": str(gating_csv),
    }
