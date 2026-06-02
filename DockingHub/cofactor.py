#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from .utils import ensure_dir, quote, resolve_path, run_command, safe_id, which, write_csv, resolve_external_tool


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
            query, target, qtmscore, ttmscore = parts[:4]
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
    # Foldseek target names usually correspond to file stems.
    by_stem = {p.stem: p for p in templates}
    return by_stem.get(Path(best).stem) or by_stem.get(best) or templates[0]


def _write_pymol_mapping_script(script_path: Path) -> None:
    script_path.write_text(
        r'''
import sys
from pathlib import Path
from pymol import cmd

target, template, residues, out_cofactor, out_merged = sys.argv[-5:]
residue_names = [r.strip() for r in residues.split(',') if r.strip()]
cmd.reinitialize()
cmd.load(target, 'target')
cmd.load(template, 'template')
# Align template protein onto target protein; mapped HETATM/cofactor moves with the template object.
cmd.align('template and polymer.protein', 'target and polymer.protein')
if residue_names:
    sel = 'template and resn ' + '+'.join(residue_names)
else:
    sel = 'template and hetatm and not solvent'
cmd.create('mapped_cofactor', sel)
cmd.save(out_cofactor, 'mapped_cofactor')
cmd.create('merged_receptor', 'target or mapped_cofactor')
cmd.save(out_merged, 'merged_receptor')
'''.strip() + "\n",
        encoding="utf-8",
    )


def map_cofactor_with_pymol(config: dict[str, Any], target_structure: Path, template: Path, out_cofactor: Path, out_merged: Path) -> tuple[str, str]:
    cof_cfg = config.get("cofactor", {})
    align_cfg = cof_cfg.get("alignment", {})
    residues = ",".join(cof_cfg.get("cofactor_residue_names", []))
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    executable = resolve_external_tool("pymol", align_cfg.get("pymol_executable", "auto"), root, config)
    if which(str(executable)) is None:
        return "failed", f"PyMOL executable not found: {executable}"
    script_path = out_cofactor.parent / "map_cofactor_pymol.py"
    _write_pymol_mapping_script(script_path)
    cmd = f"{executable} -cq {quote(script_path)} -- {quote(target_structure)} {quote(template)} {quote(residues)} {quote(out_cofactor)} {quote(out_merged)}"
    rc = run_command(cmd, timeout=int(align_cfg.get("timeout", 3600) or 3600), log_path=out_cofactor.parent / "map_cofactor.log")
    if rc != 0 or not out_cofactor.exists() or not out_merged.exists():
        return "failed", f"PyMOL mapping returned {rc}"
    return "success", "mapped_by_pymol"


def build_cofactor_manifest(config: dict[str, Any], conformers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    cof_cfg = config.get("cofactor", {})
    out_dir = resolve_path(config.get("paths", {}).get("cofactor_mapped_dir", "data/data_output/cofactor_mapped"), root)
    template_root = resolve_path(config.get("paths", {}).get("cofactor_dir", "data/data_input/cofactor"), root)
    assert out_dir is not None and template_root is not None
    ensure_dir(out_dir)
    enabled = bool(cof_cfg.get("enabled", False))
    use_foldseek = bool(cof_cfg.get("use_foldseek", False))
    fallback_to_first = bool(cof_cfg.get("fallback_to_first", True))
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
            rows.append({
                **conf,
                "cofactor_enabled": str(bool(enabled)).lower(),
                "template_path": "",
                "mapped_cofactor_path": "",
                "receptor_structure_path": "",
                "status": "failed",
                "message": f"conformer unavailable or failed before cofactor mapping: {conf.get('message', '')}",
            })
            continue

        if not enabled:
            rows.append({
                **conf,
                "cofactor_enabled": "false",
                "template_path": "",
                "mapped_cofactor_path": "",
                "receptor_structure_path": str(target),
                "status": "success",
                "message": "cofactor disabled",
            })
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
            rows.append({
                **conf,
                "cofactor_enabled": "true",
                "template_path": "",
                "mapped_cofactor_path": "",
                "receptor_structure_path": str(target),
                "status": "failed",
                "message": f"No cofactor template found in {template_dir}",
            })
            continue
        out_cof = out_sub / f"{safe_id(pid)}_{safe_id(cid)}_mapped_cofactor.pdb"
        out_merged = out_sub / f"{safe_id(pid)}_{safe_id(cid)}_receptor_with_cofactor.pdb"
        status, message = map_cofactor_with_pymol(config, target, selected, out_cof, out_merged)
        if status != "success" and bool(cof_cfg.get("continue_without_cofactor_on_error", True)):
            out_merged = target
            message = message + "; continue_without_cofactor"
            status = "success_no_cofactor"
        rows.append({
            **conf,
            "cofactor_enabled": "true",
            "template_path": str(selected),
            "template_selection_method": method,
            "mapped_cofactor_path": str(out_cof) if out_cof.exists() else "",
            "receptor_structure_path": str(out_merged),
            "status": status,
            "message": message,
        })
    write_csv(out_dir / "cofactor_manifest.csv", rows)
    return rows
