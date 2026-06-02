#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
from typing import Any

from .utils import copy_or_link, ensure_dir, is_probably_pdbqt, quote, resolve_path, run_command, safe_id, which, write_csv, resolve_external_tool


def prepare_receptors(config: dict[str, Any], cofactor_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    prep_cfg = config.get("receptor_preparation", {})
    out_dir = resolve_path(config.get("paths", {}).get("receptor_dir", "data/data_output/receptor"), root)
    assert out_dir is not None
    ensure_dir(out_dir)
    enabled = bool(prep_cfg.get("enabled", True))
    command_template = prep_cfg.get("command_template", "{executable} -r {input_pdb} -o {output_pdbqt} -A hydrogens {extra_options}")
    executable = resolve_external_tool("prepare_receptor", prep_cfg.get("executable", "auto"), root, config)
    extra_options = prep_cfg.get("extra_options", "")
    timeout = int(prep_cfg.get("timeout", 3600) or 3600)
    overwrite = bool(config.get("output", {}).get("overwrite", False))
    file_action = config.get("output", {}).get("file_action", "copy")
    rows: list[dict[str, Any]] = []
    success_like = {"", "success", "success_no_cofactor", "ok", "true", "1", "skipped"}

    for row in cofactor_rows:
        upstream_status = str(row.get("status", "success")).strip().lower()
        raw_src = row.get("receptor_structure_path") or row.get("structure_path") or ""
        src = Path(raw_src) if raw_src else Path("")
        batch = row.get("batch_id", "file_1") or "file_1"
        pid = row.get("protein_id", src.stem if raw_src else "protein")
        cid = row.get("conformer_id", "conf_0")
        out_sub = out_dir / batch
        ensure_dir(out_sub)
        out_pdbqt = out_sub / f"{safe_id(pid)}__{safe_id(cid)}.pdbqt"
        status = "success"
        message = ""

        if upstream_status not in success_like or not raw_src or not src.exists():
            rows.append({
                **row,
                "receptor_pdbqt_path": "",
                "receptor_preparation_status": "failed",
                "receptor_preparation_message": f"upstream receptor structure unavailable or failed: {row.get('message', '')}",
            })
            continue
        if is_probably_pdbqt(src):
            copy_or_link(src, out_pdbqt, file_action, overwrite=overwrite)
            message = "input already pdbqt"
        elif not enabled:
            status = "failed"
            message = "receptor_preparation disabled but input is not PDBQT"
        elif out_pdbqt.exists() and not overwrite:
            message = "skip_existing"
        else:
            if which(str(executable)) is None and "{executable}" in command_template:
                status = "failed"
                message = f"prepare receptor executable not found: {executable}"
            else:
                command = command_template.format(
                    executable=executable,
                    input_pdb=quote(src),
                    output_pdbqt=quote(out_pdbqt),
                    extra_options=extra_options,
                    protein_id=pid,
                    batch_id=batch,
                    conformer_id=cid,
                )
                rc = run_command(command, timeout=timeout, log_path=out_sub / f"{safe_id(pid)}__{safe_id(cid)}_prepare_receptor.log")
                if rc != 0 or not out_pdbqt.exists():
                    status = "failed"
                    message = f"prepare receptor returned {rc}"
                else:
                    message = "prepared"
        rows.append({
            **row,
            "receptor_pdbqt_path": str(out_pdbqt) if out_pdbqt.exists() else "",
            "receptor_preparation_status": status,
            "receptor_preparation_message": message,
        })
    write_csv(out_dir / "receptor_manifest.csv", rows)
    return rows
