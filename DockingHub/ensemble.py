#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
from typing import Any

from .manifest import ProteinRecord
from .utils import ensure_dir, quote, resolve_path, run_command, safe_id, which, write_csv, resolve_external_tool


def _collect_conformers(out_dir: Path, patterns: list[str]) -> list[Path]:
    files: list[Path] = []
    for pattern in patterns:
        files.extend(out_dir.glob(pattern))
    return sorted(set([p for p in files if p.is_file()]))


def _write_pymol_align_script(script_path: Path) -> None:
    script_path.write_text(
        r'''
import sys
from pymol import cmd

reference, mobile, output_pdb = sys.argv[-3:]
cmd.reinitialize()
cmd.load(reference, 'reference')
cmd.load(mobile, 'mobile')
# Align the sampled ensemble conformer back to the original structure frame.
# This keeps the original TApocket-derived pocket box coordinates valid.
cmd.align('mobile and polymer.protein', 'reference and polymer.protein')
cmd.save(output_pdb, 'mobile')
'''.strip() + "\n",
        encoding="utf-8",
    )


def _align_to_reference(config: dict[str, Any], mobile: Path, reference: Path, out_path: Path, log_path: Path) -> tuple[Path, str, str]:
    """Align one ensemble conformer to the original input structure coordinate frame.

    TApocket is intentionally run only on the original protein structure. Therefore,
    all AlphaFlow/ensemble conformers must be brought back to the original structure
    frame before docking so that the original pocket center/box remains meaningful.
    """
    align_cfg = config.get("ensemble", {}).get("align_to_reference", {})
    enabled = bool(align_cfg.get("enabled", True))
    fallback = bool(align_cfg.get("fallback_to_unaligned", True))
    if not enabled:
        return mobile, "skipped", "align_to_reference disabled"
    try:
        if mobile.resolve() == reference.resolve():
            return mobile, "success", "input structure; no alignment needed"
    except Exception:
        pass

    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    executable = resolve_external_tool("pymol", align_cfg.get("pymol_executable", "auto"), root, config)
    if which(str(executable)) is None:
        if fallback:
            return mobile, "success_unaligned", f"PyMOL not found: {executable}; using unaligned conformer"
        return mobile, "failed", f"PyMOL not found: {executable}"

    ensure_dir(out_path.parent)
    script_path = out_path.parent / "align_to_reference_pymol.py"
    _write_pymol_align_script(script_path)
    cmd = f"{executable} -cq {quote(script_path)} -- {quote(reference)} {quote(mobile)} {quote(out_path)}"
    rc = run_command(cmd, timeout=int(align_cfg.get("timeout", 3600) or 3600), log_path=log_path)
    if rc == 0 and out_path.exists():
        return out_path, "success", "aligned_to_original_structure"
    if fallback:
        return mobile, "success_unaligned", f"alignment failed with code {rc}; using unaligned conformer"
    return mobile, "failed", f"alignment failed with code {rc}"


def build_conformer_manifest(config: dict[str, Any], proteins: list[ProteinRecord]) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    ens_cfg = config.get("ensemble", {})
    out_dir = resolve_path(config.get("paths", {}).get("ensemble_dir", "data/data_output/ensemble"), root)
    assert out_dir is not None
    ensure_dir(out_dir)

    enabled = bool(ens_cfg.get("enabled", False))
    rows: list[dict[str, Any]] = []
    if not enabled:
        for protein in proteins:
            rows.append({
                "protein_id": protein.protein_id,
                "cluster_id": protein.cluster_id,
                "batch_id": protein.batch_id,
                "conformer_id": "conf_0",
                "structure_path": str(protein.protein_path),
                "original_structure_path": str(protein.protein_path),
                "unaligned_structure_path": str(protein.protein_path),
                "coordinate_frame": "original_structure",
                "alignment_status": "success",
                "alignment_message": "ensemble disabled; original structure used",
                "source": "input_protein",
                "status": "success",
                "message": "ensemble disabled",
            })
        write_csv(out_dir / "conformer_manifest.csv", rows)
        return rows

    command_template = str(ens_cfg.get("command_template", "")).strip()
    collect_patterns = ens_cfg.get("collect_patterns", ["*.pdb"])
    max_confs = int(ens_cfg.get("max_conformers_per_protein", 0) or 0)
    timeout = int(ens_cfg.get("timeout", 86400) or 86400)
    fallback = bool(ens_cfg.get("fallback_to_input", True))

    for protein in proteins:
        pid_safe = safe_id(protein.protein_id)
        p_out = out_dir / protein.batch_id / pid_safe
        ensure_dir(p_out)
        log_path = p_out / "alphaflow.log"
        status = "success"
        message = ""
        if command_template:
            command = command_template.format(
                input_pdb=quote(protein.protein_path),
                output_dir=quote(p_out),
                protein_id=protein.protein_id,
                batch_id=protein.batch_id,
                threads=int(ens_cfg.get("threads", 1) or 1),
            )
            rc = run_command(command, timeout=timeout, log_path=log_path)
            if rc != 0:
                status = "failed"
                message = f"AlphaFlow command returned {rc}"
        else:
            status = "skipped"
            message = "No ensemble.command_template provided"

        conformers = _collect_conformers(p_out, list(collect_patterns))
        # Avoid treating already aligned outputs from a previous run as new AlphaFlow conformers.
        aligned_root = p_out / "aligned_to_original"
        conformers = [p for p in conformers if aligned_root not in p.parents]
        if max_confs > 0:
            conformers = conformers[:max_confs]
        if not conformers and fallback:
            conformers = [protein.protein_path]
            status = "success"
            message = (message + "; " if message else "") + "fallback_to_input"

        for idx, raw_conf in enumerate(conformers):
            aligned_path = aligned_root / f"{pid_safe}__conf_{idx}.aligned_to_original.pdb"
            aligned_conf, align_status, align_message = _align_to_reference(
                config,
                Path(raw_conf),
                Path(protein.protein_path),
                aligned_path,
                aligned_root / f"{pid_safe}__conf_{idx}.align.log",
            )
            row_status = status
            if align_status == "failed":
                row_status = "failed"
            rows.append({
                "protein_id": protein.protein_id,
                "cluster_id": protein.cluster_id,
                "batch_id": protein.batch_id,
                "conformer_id": f"conf_{idx}",
                # Downstream modules must use this aligned structure path.
                "structure_path": str(aligned_conf),
                "original_structure_path": str(protein.protein_path),
                "unaligned_structure_path": str(raw_conf),
                "coordinate_frame": "original_structure" if align_status in {"success", "skipped"} else "unknown_or_unaligned",
                "alignment_status": align_status,
                "alignment_message": align_message,
                "source": "alphaflow" if Path(raw_conf) != protein.protein_path else "input_protein",
                "status": row_status,
                "message": message,
            })

    write_csv(out_dir / "conformer_manifest.csv", rows)
    return rows
