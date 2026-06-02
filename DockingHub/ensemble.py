#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
import os
import shutil
from pathlib import Path
from typing import Any

from .manifest import ProteinRecord
from .utils import ensure_dir, quote, resolve_path, run_command, safe_id, which, write_csv, resolve_external_tool

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z",
    "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H",
    "MSE": "M",
}


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


def _path_from_config(value: Any, root: Path, base: Path | None = None) -> Path | None:
    if value in {None, ""}:
        return None
    text = os.path.expandvars(str(value)).strip()
    if "$" in text:
        return None
    path = Path(text).expanduser()
    if path.is_absolute():
        return path.resolve()
    if base is not None:
        return (base / path).resolve()
    return (root / path).resolve()


def _command_executable(value: Any, env_name: str, default: str = "python") -> str:
    raw = str(value or "").strip() or os.environ.get(env_name, "") or default
    expanded = os.path.expandvars(raw)
    return expanded if "$" not in expanded else default


def _pdb_sequence(path: Path) -> str:
    seen: set[tuple[str, str, str]] = set()
    seq: list[str] = []
    with path.open("r", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM") or len(line) < 27:
                continue
            resname = line[17:20].strip().upper()
            aa = AA3_TO_1.get(resname)
            if not aa:
                continue
            key = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip())
            if key in seen:
                continue
            seen.add(key)
            seq.append(aa)
    return "".join(seq)


def _write_alphaflow_csv(path: Path, name: str, seqres: str) -> Path:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "seqres"])
        writer.writeheader()
        writer.writerow({"name": name, "seqres": seqres})
    return path


def _stage_template(src: Path, templates_dir: Path, name: str, mode: str) -> Path | None:
    if mode == "none":
        return None
    dst = templates_dir / f"{name}.pdb"
    if mode == "existing":
        return dst if dst.exists() else None
    ensure_dir(templates_dir)
    shutil.copy2(src, dst)
    return dst


def _msa_exists(msa_dir: Path, name: str) -> bool:
    return (msa_dir / name / "a3m" / f"{name}.a3m").is_file()


def _env_prefix(af_cfg: dict[str, Any]) -> str:
    assignments: list[str] = []
    cuda_raw = af_cfg.get("cuda_visible_devices", "") or os.environ.get("ALPHAFLOW_CUDA_VISIBLE_DEVICES", "")
    cuda_visible_devices = os.path.expandvars(str(cuda_raw)).strip()
    if cuda_visible_devices:
        assignments.append(f"CUDA_VISIBLE_DEVICES={quote(cuda_visible_devices)}")
    alloc_raw = af_cfg.get("pytorch_cuda_alloc_conf", "") or os.environ.get("ALPHAFLOW_PYTORCH_CUDA_ALLOC_CONF", "")
    alloc_conf = os.path.expandvars(str(alloc_raw)).strip()
    if alloc_conf:
        assignments.append(f"PYTORCH_CUDA_ALLOC_CONF={quote(alloc_conf)}")
    return " ".join(assignments) + (" " if assignments else "")


def _split_multimodel_pdb(pdb_path: Path, out_dir: Path, name: str) -> list[Path]:
    lines = pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if not any(line.startswith("MODEL") for line in lines):
        return [pdb_path]
    ensure_dir(out_dir)
    outputs: list[Path] = []
    current: list[str] = []
    in_model = False
    model_idx = 0
    for line in lines:
        if line.startswith("MODEL"):
            current = [line]
            in_model = True
            model_idx += 1
            continue
        if in_model:
            current.append(line)
            if line.startswith("ENDMDL"):
                out_path = out_dir / f"{name}__sample_{model_idx:03d}.pdb"
                out_path.write_text("\n".join(current) + "\n", encoding="utf-8")
                outputs.append(out_path)
                current = []
                in_model = False
    if in_model and current:
        out_path = out_dir / f"{name}__sample_{model_idx:03d}.pdb"
        out_path.write_text("\n".join(current) + "\nEND\n", encoding="utf-8")
        outputs.append(out_path)
    return outputs or [pdb_path]


def _collect_alphaflow_conformers(raw_out_dir: Path, split_dir: Path, name: str, patterns: list[str]) -> list[Path]:
    named = raw_out_dir / f"{name}.pdb"
    if named.exists():
        return _split_multimodel_pdb(named, split_dir, name)
    return _collect_conformers(raw_out_dir, patterns)


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


def _build_alphaflow_conformer_manifest(config: dict[str, Any], proteins: list[ProteinRecord], out_dir: Path) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    ens_cfg = config.get("ensemble", {})
    af_cfg = ens_cfg.get("alphaflow", {}) or {}
    project_dir = _path_from_config(af_cfg.get("project_dir") or os.environ.get("ALPHAFLOW_ROOT", ""), root)
    python_exec = _command_executable(af_cfg.get("python_executable"), "ALPHAFLOW_PYTHON", "python")
    run_msa = bool(af_cfg.get("run_msa", False))
    run_prediction = bool(af_cfg.get("run_prediction", True))
    require_msa = bool(af_cfg.get("require_msa", True))
    fallback = bool(ens_cfg.get("fallback_to_input", True))
    max_confs = int(ens_cfg.get("max_conformers_per_protein", 0) or 0)
    collect_patterns = list(ens_cfg.get("collect_patterns", ["*.pdb"]))
    timeout = int(ens_cfg.get("timeout", 86400) or 86400)
    template_mode = str(af_cfg.get("template_mode", "copy_input"))
    rows: list[dict[str, Any]] = []

    for protein in proteins:
        pid_safe = safe_id(protein.protein_id)
        p_out = out_dir / protein.batch_id / pid_safe
        input_csv_dir = p_out / "input_csv"
        templates_dir = _path_from_config((af_cfg.get("paths", {}) or {}).get("templates_dir"), root) or (p_out / "templates_dir")
        msa_dir = _path_from_config((af_cfg.get("paths", {}) or {}).get("msa_dir"), root) or (p_out / "msa_alignment")
        raw_out_dir = p_out / "raw_output"
        split_dir = p_out / "split_models"
        ensure_dir(p_out)

        status = "success"
        message = ""
        conformers: list[Path] = []
        seqres = _pdb_sequence(Path(protein.protein_path))
        if not seqres:
            status = "failed"
            message = "could not extract protein sequence from PDB"
        else:
            input_csv = _write_alphaflow_csv(input_csv_dir / f"{pid_safe}.csv", pid_safe, seqres)
            template_path = _stage_template(Path(protein.protein_path), templates_dir, pid_safe, template_mode)
            if template_mode != "none" and template_path is None:
                status = "failed"
                message = f"AlphaFlow template not found for {pid_safe}"

        if status == "success" and run_msa:
            if project_dir is None:
                status = "failed"
                message = "AlphaFlow project_dir is required for MSA generation"
            else:
                mmseqs_script = _path_from_config(af_cfg.get("mmseqs_script", "scripts/mmseqs_query.py"), root, project_dir)
                if mmseqs_script is None or not mmseqs_script.exists():
                    status = "failed"
                    message = f"AlphaFlow mmseqs script not found: {mmseqs_script}"
                else:
                    cmd = f"{_env_prefix(af_cfg)}{quote(python_exec)} {quote(mmseqs_script)} --split {quote(input_csv)} --outdir {quote(msa_dir)}"
                    rc = run_command(cmd, timeout=timeout, cwd=project_dir, log_path=p_out / "alphaflow_msa.log")
                    if rc != 0:
                        status = "failed"
                        message = f"AlphaFlow MSA command returned {rc}"

        if status == "success" and run_prediction:
            if project_dir is None:
                status = "failed"
                message = "AlphaFlow project_dir is required for prediction"
            elif require_msa and not _msa_exists(msa_dir, pid_safe):
                status = "failed"
                message = f"AlphaFlow MSA not found for {pid_safe}: {msa_dir / pid_safe / 'a3m' / (pid_safe + '.a3m')}"
            else:
                predict_script = _path_from_config(af_cfg.get("predict_script", "predict2.py"), root, project_dir)
                weights = _path_from_config(af_cfg.get("weights"), root, project_dir)
                if predict_script is None or not predict_script.exists():
                    status = "failed"
                    message = f"AlphaFlow predict script not found: {predict_script}"
                elif weights is None or not weights.exists():
                    status = "failed"
                    message = f"AlphaFlow weights not found: {weights}"
                else:
                    ensure_dir(raw_out_dir)
                    extra_args = str(af_cfg.get("extra_args", "") or "").strip()
                    cmd = (
                        f"{_env_prefix(af_cfg)}{quote(python_exec)} {quote(predict_script)} "
                        f"--mode {quote(str(af_cfg.get('mode', 'alphafold')))} "
                        f"--input_csv {quote(input_csv)} "
                        f"--msa_dir {quote(msa_dir)} "
                        f"--weights {quote(weights)} "
                        f"--samples {int(af_cfg.get('samples', 5) or 5)} "
                        f"--steps {int(af_cfg.get('steps', 10) or 10)} "
                        f"--outpdb {quote(raw_out_dir)} "
                        f"--templates_dir {quote(templates_dir)} "
                        f"{extra_args}"
                    ).strip()
                    rc = run_command(cmd, timeout=timeout, cwd=project_dir, log_path=p_out / "alphaflow_predict.log")
                    if rc != 0:
                        status = "failed"
                        message = f"AlphaFlow predict command returned {rc}"
                    else:
                        conformers = _collect_alphaflow_conformers(raw_out_dir, split_dir, pid_safe, collect_patterns)
                        if max_confs > 0:
                            conformers = conformers[:max_confs]
                        if not conformers:
                            status = "failed"
                            message = "AlphaFlow produced no PDB conformers"

        if not conformers and fallback:
            conformers = [Path(protein.protein_path)]
            if status != "success":
                message = (message + "; " if message else "") + "fallback_to_input"
            status = "success"

        for idx, raw_conf in enumerate(conformers):
            aligned_root = p_out / "aligned_to_original"
            aligned_path = aligned_root / f"{pid_safe}__conf_{idx}.aligned_to_original.pdb"
            aligned_conf, align_status, align_message = _align_to_reference(
                config,
                Path(raw_conf),
                Path(protein.protein_path),
                aligned_path,
                aligned_root / f"{pid_safe}__conf_{idx}.align.log",
            )
            row_status = status if align_status != "failed" else "failed"
            rows.append({
                "protein_id": protein.protein_id,
                "cluster_id": protein.cluster_id,
                "batch_id": protein.batch_id,
                "conformer_id": f"conf_{idx}",
                "structure_path": str(aligned_conf),
                "original_structure_path": str(protein.protein_path),
                "unaligned_structure_path": str(raw_conf),
                "coordinate_frame": "original_structure" if align_status in {"success", "skipped"} else "unknown_or_unaligned",
                "alignment_status": align_status,
                "alignment_message": align_message,
                "source": "alphaflow" if Path(raw_conf) != Path(protein.protein_path) else "input_protein",
                "status": row_status,
                "message": message,
                "alphaflow_name": pid_safe,
                "alphaflow_project_dir": str(project_dir) if project_dir else "",
                "alphaflow_msa_dir": str(msa_dir),
                "alphaflow_templates_dir": str(templates_dir),
            })

    write_csv(out_dir / "conformer_manifest.csv", rows)
    return rows


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

    engine = str(ens_cfg.get("engine", "command")).strip().lower()
    if engine == "alphaflow":
        return _build_alphaflow_conformer_manifest(config, proteins, out_dir)

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
