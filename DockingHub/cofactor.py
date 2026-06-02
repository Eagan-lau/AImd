#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from typing import Any

from .utils import ensure_dir, quote, resolve_external_tool, resolve_path, run_command, safe_id, which, write_csv


POCKET_VALIDATION_DEFAULTS: dict[str, Any] = {
    "enabled": True,
    "pocket_radius": 5.0,
    "rmsd_metric": "ca",
    "local_rmsd_cutoff": 2.5,
    "min_mapped_residues": 5,
    "min_pocket_coverage": 0.70,
    "require_pass_before_transfer": True,
}

COFACTOR_VALIDATION_FIELDS = [
    "cofactor_pocket_radius",
    "cofactor_pocket_residue_count",
    "mapped_pocket_residue_count",
    "cofactor_pocket_coverage",
    "cofactor_site_ca_rmsd",
    "cofactor_local_rmsd_cutoff",
    "cofactor_transfer_pass",
    "cofactor_validation_status",
    "cofactor_validation_message",
    "transfer_mode",
]

SOLVENT_RESIDUES = {"HOH", "WAT", "DOD", "SOL"}


def _template_files(template_dir: Path) -> list[Path]:
    if not template_dir.exists():
        return []
    files: list[Path] = []
    for ext in ("*.pdb", "*.ent", "*.cif"):
        files.extend(template_dir.glob(ext))
    return sorted([p for p in files if p.is_file()])


def _parse_foldseek_best(tsv: Path) -> str | None:
    best_target = None
    best_score = -1.0
    if not tsv.exists():
        return None
    with tsv.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                parts = line.split()
            if len(parts) < 4:
                continue
            _query, target, qtmscore, ttmscore = parts[:4]
            try:
                score = min(float(qtmscore), float(ttmscore))
            except ValueError:
                continue
            if score > best_score:
                best_score = score
                best_target = target
    return best_target


def select_template_by_foldseek(config: dict[str, Any], target_structure: Path, template_dir: Path, work_dir: Path) -> Path | None:
    fold_cfg = config.get("cofactor", {}).get("foldseek", {})
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    executable = resolve_external_tool("foldseek", fold_cfg.get("executable", "auto"), root, config)
    if which(str(executable)) is None:
        return None
    templates = _template_files(template_dir)
    if not templates:
        return None
    ensure_dir(work_dir)
    template_db = work_dir / "template_db"
    query_db = work_dir / "query_db"
    result_db = work_dir / "result_db"
    tmp_dir = work_dir / "tmp"
    out_tsv = work_dir / "foldseek_hits.tsv"
    fmt = "query,target,qtmscore,ttmscore"
    commands = [
        f"{executable} createdb {quote(template_dir)} {quote(template_db)}",
        f"{executable} createdb {quote(target_structure)} {quote(query_db)}",
        f"{executable} search {quote(query_db)} {quote(template_db)} {quote(result_db)} {quote(tmp_dir)} -a --threads {int(fold_cfg.get('threads', 4) or 4)}",
        f"{executable} convertalis {quote(query_db)} {quote(template_db)} {quote(result_db)} {quote(out_tsv)} --format-output {quote(fmt)}",
    ]
    log_path = work_dir / "foldseek_template_select.log"
    with log_path.open("w", encoding="utf-8") as log:
        for command in commands:
            rc = run_command(command, timeout=int(fold_cfg.get("timeout", 86400) or 86400), log_path=work_dir / "last_cmd.log")
            log.write(f"$ {command}\nreturn_code={rc}\n")
            if rc != 0:
                return None
    best = _parse_foldseek_best(out_tsv)
    if not best:
        return None
    by_stem = {p.stem: p for p in templates}
    return by_stem.get(Path(best).stem) or by_stem.get(best) or templates[0]


def _as_bool(value: Any, default: bool) -> bool:
    if value in {None, ""}:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on", "success", "ok"}


def _as_float(value: Any, default: float) -> float:
    try:
        if value in {None, ""}:
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _as_int(value: Any, default: int) -> int:
    try:
        if value in {None, ""}:
            return default
        return int(value)
    except (TypeError, ValueError):
        return default


def get_pocket_validation_config(config: dict[str, Any]) -> dict[str, Any]:
    align_cfg = config.get("cofactor", {}).get("alignment", {}) or {}
    raw_cfg = align_cfg.get("pocket_validation", {}) or {}
    cfg = dict(POCKET_VALIDATION_DEFAULTS)
    cfg.update(raw_cfg)
    cfg["enabled"] = _as_bool(cfg.get("enabled"), bool(POCKET_VALIDATION_DEFAULTS["enabled"]))
    cfg["pocket_radius"] = _as_float(cfg.get("pocket_radius"), float(POCKET_VALIDATION_DEFAULTS["pocket_radius"]))
    cfg["rmsd_metric"] = str(cfg.get("rmsd_metric") or "ca").strip().lower()
    if cfg["rmsd_metric"] != "ca":
        cfg["rmsd_metric"] = "ca"
    cfg["local_rmsd_cutoff"] = _as_float(cfg.get("local_rmsd_cutoff"), float(POCKET_VALIDATION_DEFAULTS["local_rmsd_cutoff"]))
    cfg["min_mapped_residues"] = _as_int(cfg.get("min_mapped_residues"), int(POCKET_VALIDATION_DEFAULTS["min_mapped_residues"]))
    cfg["min_pocket_coverage"] = _as_float(cfg.get("min_pocket_coverage"), float(POCKET_VALIDATION_DEFAULTS["min_pocket_coverage"]))
    cfg["require_pass_before_transfer"] = _as_bool(
        cfg.get("require_pass_before_transfer"),
        bool(POCKET_VALIDATION_DEFAULTS["require_pass_before_transfer"]),
    )
    return cfg


def _validation_fields(
    validation_cfg: dict[str, Any],
    *,
    pocket_count: int = 0,
    mapped_count: int = 0,
    coverage: float = 0.0,
    rmsd: float | str | None = "",
    transfer_pass: bool = False,
    validation_status: str = "",
    validation_message: str = "",
    transfer_mode: str = "",
) -> dict[str, Any]:
    return {
        "cofactor_pocket_radius": validation_cfg["pocket_radius"],
        "cofactor_pocket_residue_count": pocket_count,
        "mapped_pocket_residue_count": mapped_count,
        "cofactor_pocket_coverage": coverage,
        "cofactor_site_ca_rmsd": "" if rmsd is None else rmsd,
        "cofactor_local_rmsd_cutoff": validation_cfg["local_rmsd_cutoff"],
        "cofactor_transfer_pass": transfer_pass,
        "cofactor_validation_status": validation_status,
        "cofactor_validation_message": validation_message,
        "transfer_mode": transfer_mode,
    }


def _atom_element(line: str, atom_name: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.upper()
    stripped = atom_name.strip()
    if not stripped:
        return ""
    if stripped[0].isdigit() and len(stripped) > 1:
        return stripped[1].upper()
    return stripped[0].upper()


def _is_hydrogen_atom(atom: dict[str, Any]) -> bool:
    element = str(atom.get("element", "")).strip().upper()
    name = str(atom.get("name", "")).strip().upper()
    if element == "H":
        return True
    if not name:
        return False
    if name[0].isdigit() and len(name) > 1:
        return name[1] == "H"
    return name.startswith("H")


def _parse_pdb_atoms(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"} or len(line) < 54:
                continue
            try:
                coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            atom_name = line[12:16].strip()
            atoms.append(
                {
                    "record": record,
                    "name": atom_name,
                    "resn": line[17:20].strip(),
                    "chain": line[21:22].strip(),
                    "resi": line[22:26].strip(),
                    "icode": line[26:27].strip(),
                    "element": _atom_element(line, atom_name),
                    "coord": coord,
                }
            )
    return atoms


def _residue_key(atom: dict[str, Any]) -> tuple[str, str, str, str]:
    return (
        str(atom.get("chain", "")).strip(),
        str(atom.get("resi", "")).strip(),
        str(atom.get("icode", "")).strip(),
        str(atom.get("resn", "")).strip(),
    )


def _payload_residue_key(payload: dict[str, Any]) -> tuple[str, str, str, str]:
    return (
        str(payload.get("chain", "")).strip(),
        str(payload.get("resi", "")).strip(),
        str(payload.get("icode", "")).strip(),
        str(payload.get("resn", "")).strip(),
    )


def _payload_coord(payload: dict[str, Any]) -> tuple[float, float, float] | None:
    coord = payload.get("coord")
    if not isinstance(coord, (list, tuple)) or len(coord) != 3:
        return None
    try:
        return float(coord[0]), float(coord[1]), float(coord[2])
    except (TypeError, ValueError):
        return None


def _distance_squared(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def _rmsd(coord_pairs: list[tuple[tuple[float, float, float], tuple[float, float, float]]]) -> float:
    if not coord_pairs:
        return math.nan
    total = 0.0
    for left, right in coord_pairs:
        total += _distance_squared(left, right)
    return math.sqrt(total / len(coord_pairs))


def _cofactor_heavy_atoms(template_atoms: list[dict[str, Any]], cofactor_residue_names: list[str]) -> list[dict[str, Any]]:
    names = {name.strip().upper() for name in cofactor_residue_names if str(name).strip()}
    atoms: list[dict[str, Any]] = []
    for atom in template_atoms:
        resn = str(atom.get("resn", "")).upper()
        if _is_hydrogen_atom(atom):
            continue
        if names:
            if resn in names:
                atoms.append(atom)
        elif atom.get("record") == "HETATM" and resn not in SOLVENT_RESIDUES:
            atoms.append(atom)
    return atoms


def _template_pocket_residues(
    template_atoms: list[dict[str, Any]],
    cofactor_atoms: list[dict[str, Any]],
    radius: float,
) -> set[tuple[str, str, str, str]]:
    radius_sq = radius * radius
    pocket: set[tuple[str, str, str, str]] = set()
    protein_atoms = [atom for atom in template_atoms if atom.get("record") == "ATOM" and not _is_hydrogen_atom(atom)]
    cofactor_coords = [atom["coord"] for atom in cofactor_atoms]
    for protein_atom in protein_atoms:
        coord = protein_atom["coord"]
        if any(_distance_squared(coord, cofactor_coord) <= radius_sq for cofactor_coord in cofactor_coords):
            pocket.add(_residue_key(protein_atom))
    return pocket


def compute_cofactor_pocket_validation(
    template_atoms: list[dict[str, Any]],
    alignment_pairs: list[dict[str, Any]],
    cofactor_residue_names: list[str],
    validation_cfg: dict[str, Any],
) -> dict[str, Any]:
    cofactor_atoms = _cofactor_heavy_atoms(template_atoms, cofactor_residue_names)
    if not cofactor_atoms:
        return _validation_fields(
            validation_cfg,
            validation_status="missing_cofactor_atoms",
            validation_message="No configured cofactor heavy atoms were found in the aligned template.",
        )

    pocket = _template_pocket_residues(template_atoms, cofactor_atoms, float(validation_cfg["pocket_radius"]))
    pocket_count = len(pocket)
    if pocket_count == 0:
        return _validation_fields(
            validation_cfg,
            validation_status="missing_pocket_residues",
            validation_message="No template protein pocket residues were found near the configured cofactor heavy atoms.",
        )

    if not alignment_pairs:
        return _validation_fields(
            validation_cfg,
            pocket_count=pocket_count,
            validation_status="alignment_unavailable",
            validation_message="No global alignment residue pairs were available for cofactor pocket validation.",
        )

    coord_pairs: list[tuple[tuple[float, float, float], tuple[float, float, float]]] = []
    seen_template_residues: set[tuple[str, str, str, str]] = set()
    for pair in alignment_pairs:
        template_atom = pair.get("template", {}) or {}
        target_atom = pair.get("target", {}) or {}
        if str(template_atom.get("name", "")).strip().upper() != "CA":
            continue
        if str(target_atom.get("name", "")).strip().upper() != "CA":
            continue
        template_key = _payload_residue_key(template_atom)
        if template_key not in pocket or template_key in seen_template_residues:
            continue
        template_coord = _payload_coord(template_atom)
        target_coord = _payload_coord(target_atom)
        if template_coord is None or target_coord is None:
            continue
        seen_template_residues.add(template_key)
        coord_pairs.append((template_coord, target_coord))

    mapped_count = len(coord_pairs)
    coverage = mapped_count / pocket_count if pocket_count else 0.0
    if mapped_count < int(validation_cfg["min_mapped_residues"]):
        return _validation_fields(
            validation_cfg,
            pocket_count=pocket_count,
            mapped_count=mapped_count,
            coverage=coverage,
            validation_status="too_few_mapped_residues",
            validation_message=(
                f"Only {mapped_count} cofactor pocket CA residue pairs were mapped; "
                f"minimum required is {validation_cfg['min_mapped_residues']}."
            ),
        )
    if coverage < float(validation_cfg["min_pocket_coverage"]):
        return _validation_fields(
            validation_cfg,
            pocket_count=pocket_count,
            mapped_count=mapped_count,
            coverage=coverage,
            validation_status="low_pocket_coverage",
            validation_message=(
                f"Cofactor pocket mapping coverage is {coverage:.3f}; "
                f"minimum required is {validation_cfg['min_pocket_coverage']:.3f}."
            ),
        )

    local_rmsd = _rmsd(coord_pairs)
    if math.isnan(local_rmsd):
        return _validation_fields(
            validation_cfg,
            pocket_count=pocket_count,
            mapped_count=mapped_count,
            coverage=coverage,
            validation_status="alignment_unavailable",
            validation_message="Cofactor pocket RMSD could not be computed from the mapped CA atoms.",
        )
    if local_rmsd > float(validation_cfg["local_rmsd_cutoff"]):
        return _validation_fields(
            validation_cfg,
            pocket_count=pocket_count,
            mapped_count=mapped_count,
            coverage=coverage,
            rmsd=local_rmsd,
            validation_status="high_local_rmsd",
            validation_message=(
                f"Cofactor pocket CA RMSD is {local_rmsd:.3f} Angstrom; "
                f"cutoff is {validation_cfg['local_rmsd_cutoff']:.3f} Angstrom."
            ),
        )

    return _validation_fields(
        validation_cfg,
        pocket_count=pocket_count,
        mapped_count=mapped_count,
        coverage=coverage,
        rmsd=local_rmsd,
        transfer_pass=True,
        validation_status="success",
        validation_message=(
            f"Cofactor pocket CA RMSD is {local_rmsd:.3f} Angstrom with "
            f"{coverage:.3f} mapped pocket coverage."
        ),
        transfer_mode="global_alignment_with_local_pocket_gate",
    )


def _write_pymol_mapping_script(script_path: Path) -> None:
    script_path.write_text(
        r'''
import json
import sys
from pymol import cmd

target, template, residues, out_cofactor_candidate, out_aligned_template, out_mapping_json = sys.argv[-6:]
residue_names = [r.strip() for r in residues.split(",") if r.strip()]


def write_mapping(payload, return_code=0):
    with open(out_mapping_json, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    if return_code:
        sys.exit(return_code)


def atom_payload(model_name, atom_index):
    model = cmd.get_model(f"{model_name} and index {atom_index}")
    if not model.atom:
        return None
    atom = model.atom[0]
    coord = atom.coord
    return {
        "chain": atom.chain,
        "resi": atom.resi,
        "icode": getattr(atom, "ins_code", "") or "",
        "resn": atom.resn,
        "name": atom.name,
        "coord": [float(coord[0]), float(coord[1]), float(coord[2])],
    }


cmd.reinitialize()
cmd.load(target, "target")
cmd.load(template, "template")

try:
    align_result = cmd.align(
        "template and polymer.protein and name CA",
        "target and polymer.protein and name CA",
        object="global_alignment",
    )
    raw_alignment = cmd.get_raw_alignment("global_alignment")
except Exception as exc:
    write_mapping({"alignment_pairs": [], "error": str(exc)}, return_code=2)

alignment_pairs = []
for column in raw_alignment:
    column_payload = {}
    for model_name, atom_index in column:
        object_name = str(model_name).split("/")[-1]
        if object_name not in {"target", "template"}:
            continue
        payload = atom_payload(object_name, atom_index)
        if payload is not None:
            column_payload[object_name] = payload
    if "template" in column_payload and "target" in column_payload:
        alignment_pairs.append({"template": column_payload["template"], "target": column_payload["target"]})

if residue_names:
    cofactor_selection = "template and resn " + "+".join(residue_names)
else:
    cofactor_selection = "template and hetatm and not solvent"
cmd.create("mapped_cofactor_candidate", cofactor_selection)
cmd.save(out_cofactor_candidate, "mapped_cofactor_candidate")
cmd.save(out_aligned_template, "template")
write_mapping({"alignment_pairs": alignment_pairs, "align_result": align_result})
'''.strip() + "\n",
        encoding="utf-8",
    )


def _load_alignment_pairs(mapping_json: Path) -> list[dict[str, Any]]:
    if not mapping_json.exists():
        return []
    try:
        with mapping_json.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return []
    pairs = payload.get("alignment_pairs", [])
    return pairs if isinstance(pairs, list) else []


def _remove_file(path: Path) -> None:
    try:
        if path.exists() or path.is_symlink():
            path.unlink()
    except OSError:
        pass


def _write_merged_receptor(target_structure: Path, mapped_cofactor: Path, out_merged: Path) -> None:
    ensure_dir(out_merged.parent)
    target_lines: list[str] = []
    with target_structure.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("END", "ENDMDL")):
                continue
            target_lines.append(line if line.endswith("\n") else line + "\n")
    cofactor_lines: list[str] = []
    with mapped_cofactor.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                cofactor_lines.append(line if line.endswith("\n") else line + "\n")
    with out_merged.open("w", encoding="utf-8") as handle:
        handle.writelines(target_lines)
        handle.writelines(cofactor_lines)
        handle.write("END\n")


def map_cofactor_with_local_gate(
    config: dict[str, Any],
    target_structure: Path,
    template: Path,
    out_cofactor: Path,
    out_merged: Path,
) -> dict[str, Any]:
    cof_cfg = config.get("cofactor", {})
    align_cfg = cof_cfg.get("alignment", {})
    validation_cfg = get_pocket_validation_config(config)
    residues = ",".join(cof_cfg.get("cofactor_residue_names", []))
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    executable = resolve_external_tool("pymol", align_cfg.get("pymol_executable", "auto"), root, config)
    if which(str(executable)) is None:
        return {
            "status": "failed",
            "message": f"PyMOL executable not found: {executable}",
            **_validation_fields(validation_cfg, validation_status="pymol_unavailable", validation_message=f"PyMOL executable not found: {executable}"),
        }

    ensure_dir(out_cofactor.parent)
    _remove_file(out_cofactor)
    _remove_file(out_merged)
    candidate_cofactor = out_cofactor.parent / f"{out_cofactor.stem}.candidate.pdb"
    aligned_template = out_cofactor.parent / f"{out_cofactor.stem}.aligned_template.pdb"
    mapping_json = out_cofactor.parent / f"{out_cofactor.stem}.alignment.json"
    for path in (candidate_cofactor, aligned_template, mapping_json):
        _remove_file(path)

    script_path = out_cofactor.parent / "map_cofactor_pymol.py"
    _write_pymol_mapping_script(script_path)
    cmd = (
        f"{executable} -cq {quote(script_path)} -- {quote(target_structure)} {quote(template)} "
        f"{quote(residues)} {quote(candidate_cofactor)} {quote(aligned_template)} {quote(mapping_json)}"
    )
    rc = run_command(cmd, timeout=int(align_cfg.get("timeout", 3600) or 3600), log_path=out_cofactor.parent / "map_cofactor.log")
    if rc != 0 or not candidate_cofactor.exists() or not aligned_template.exists():
        return {
            "status": "failed",
            "message": f"PyMOL mapping returned {rc}",
            **_validation_fields(
                validation_cfg,
                validation_status="pymol_alignment_failed",
                validation_message=f"PyMOL mapping returned {rc}.",
            ),
        }

    if not validation_cfg["enabled"]:
        shutil.copy2(candidate_cofactor, out_cofactor)
        _write_merged_receptor(target_structure, out_cofactor, out_merged)
        return {
            "status": "success",
            "message": "mapped_by_pymol_without_local_pocket_validation",
            **_validation_fields(
                validation_cfg,
                transfer_pass=True,
                validation_status="disabled",
                validation_message="Cofactor pocket validation is disabled.",
                transfer_mode="global_alignment_without_local_pocket_gate",
            ),
        }

    template_atoms = _parse_pdb_atoms(aligned_template)
    alignment_pairs = _load_alignment_pairs(mapping_json)
    validation = compute_cofactor_pocket_validation(
        template_atoms,
        alignment_pairs,
        list(cof_cfg.get("cofactor_residue_names", [])),
        validation_cfg,
    )
    require_pass = bool(validation_cfg["require_pass_before_transfer"])
    if require_pass and not bool(validation.get("cofactor_transfer_pass")):
        _remove_file(candidate_cofactor)
        return {
            "status": "failed",
            "message": str(validation.get("cofactor_validation_message", "Cofactor pocket validation failed.")),
            **validation,
            "transfer_mode": "skipped_validation_failed",
        }

    shutil.copy2(candidate_cofactor, out_cofactor)
    _write_merged_receptor(target_structure, out_cofactor, out_merged)
    validation["transfer_mode"] = validation.get("transfer_mode") or "global_alignment_with_local_pocket_gate"
    return {
        "status": "success",
        "message": "mapped_by_pymol_with_local_pocket_gate",
        **validation,
    }


def map_cofactor_with_pymol(config: dict[str, Any], target_structure: Path, template: Path, out_cofactor: Path, out_merged: Path) -> tuple[str, str]:
    result = map_cofactor_with_local_gate(config, target_structure, template, out_cofactor, out_merged)
    return str(result.get("status", "failed")), str(result.get("message", ""))


def _append_row(rows: list[dict[str, Any]], row: dict[str, Any], validation_cfg: dict[str, Any], validation: dict[str, Any] | None = None) -> None:
    merged = {**row}
    defaults = _validation_fields(validation_cfg)
    for field in COFACTOR_VALIDATION_FIELDS:
        merged[field] = defaults[field]
    if validation:
        for field in COFACTOR_VALIDATION_FIELDS:
            if field in validation:
                merged[field] = validation[field]
    rows.append(merged)


def build_cofactor_manifest(config: dict[str, Any], conformers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    cof_cfg = config.get("cofactor", {})
    validation_cfg = get_pocket_validation_config(config)
    out_dir = resolve_path(config.get("paths", {}).get("cofactor_mapped_dir", "data/data_output/cofactor_mapped"), root)
    template_root = resolve_path(config.get("paths", {}).get("cofactor_dir", "data/data_input/cofactor"), root)
    assert out_dir is not None and template_root is not None
    ensure_dir(out_dir)
    enabled = bool(cof_cfg.get("enabled", False))
    use_foldseek = bool(cof_cfg.get("use_foldseek", False))
    fallback_to_first = bool(cof_cfg.get("fallback_to_first", True))
    continue_without = bool(cof_cfg.get("continue_without_cofactor_on_error", True))
    rows: list[dict[str, Any]] = []

    success_like = {"", "success", "ok", "true", "1", "success_no_cofactor", "skipped"}

    for conf in conformers:
        target_raw = conf.get("structure_path", "")
        target = Path(target_raw) if target_raw else Path("")
        batch = conf.get("batch_id", "file_1") or "file_1"
        pid = conf.get("protein_id", target.stem if target_raw else "protein")
        cid = conf.get("conformer_id", "conf_0")
        out_sub = out_dir / batch / safe_id(pid) / safe_id(cid)
        ensure_dir(out_sub)

        conf_status = str(conf.get("status", "success")).strip().lower()
        if conf_status not in success_like or not target_raw or not target.exists():
            _append_row(
                rows,
                {
                    **conf,
                    "cofactor_enabled": str(bool(enabled)).lower(),
                    "template_path": "",
                    "mapped_cofactor_path": "",
                    "receptor_structure_path": "",
                    "status": "failed",
                    "message": f"conformer unavailable or failed before cofactor mapping: {conf.get('message', '')}",
                },
                validation_cfg,
                {
                    "cofactor_validation_status": "conformer_unavailable",
                    "cofactor_validation_message": "Conformer structure is unavailable before cofactor mapping.",
                    "transfer_mode": "skipped_validation_failed",
                },
            )
            continue

        if not enabled:
            _append_row(
                rows,
                {
                    **conf,
                    "cofactor_enabled": "false",
                    "template_path": "",
                    "mapped_cofactor_path": "",
                    "receptor_structure_path": str(target),
                    "status": "success",
                    "message": "cofactor disabled",
                },
                validation_cfg,
                {
                    "cofactor_validation_status": "disabled",
                    "cofactor_validation_message": "Cofactor transfer is disabled.",
                    "transfer_mode": "skipped_cofactor_disabled",
                },
            )
            continue
        template_dir = template_root / batch
        templates = _template_files(template_dir)
        selected: Path | None = None
        method = ""
        if use_foldseek:
            selected = select_template_by_foldseek(config, target, template_dir, out_sub / "foldseek")
            method = "foldseek"
        if selected is None and fallback_to_first and templates:
            selected = templates[0]
            method = "first_template"
        if selected is None:
            _append_row(
                rows,
                {
                    **conf,
                    "cofactor_enabled": "true",
                    "template_path": "",
                    "template_selection_method": "",
                    "mapped_cofactor_path": "",
                    "receptor_structure_path": str(target) if continue_without else "",
                    "status": "success_no_cofactor" if continue_without else "failed",
                    "message": f"No cofactor template found in {template_dir}",
                },
                validation_cfg,
                {
                    "cofactor_validation_status": "no_template",
                    "cofactor_validation_message": f"No cofactor template found in {template_dir}",
                    "transfer_mode": "skipped_no_template",
                },
            )
            continue
        out_cof = out_sub / f"{safe_id(pid)}_{safe_id(cid)}_mapped_cofactor.pdb"
        out_merged = out_sub / f"{safe_id(pid)}_{safe_id(cid)}_receptor_with_cofactor.pdb"
        result = map_cofactor_with_local_gate(config, target, selected, out_cof, out_merged)
        status = str(result.get("status", "failed"))
        message = str(result.get("message", ""))
        receptor_structure_path = str(out_merged) if status == "success" and out_merged.exists() else ""
        mapped_cofactor_path = str(out_cof) if status == "success" and out_cof.exists() else ""
        if status != "success" and continue_without:
            receptor_structure_path = str(target)
            message = message + "; continue_without_cofactor"
            status = "success_no_cofactor"
            result["transfer_mode"] = result.get("transfer_mode") or "skipped_validation_failed"
        _append_row(
            rows,
            {
                **conf,
                "cofactor_enabled": "true",
                "template_path": str(selected),
                "template_selection_method": method,
                "mapped_cofactor_path": mapped_cofactor_path,
                "receptor_structure_path": receptor_structure_path,
                "status": status,
                "message": message,
            },
            validation_cfg,
            result,
        )
    write_csv(out_dir / "cofactor_manifest.csv", rows)
    return rows
