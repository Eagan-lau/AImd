#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import shlex
from pathlib import Path
from typing import Dict, Any, List, Set

from .io import ProteinRecord
from .utils import ensure_dir, resolve_path, run_command, resolve_external_tool


def run_hipmcl(config: Dict[str, Any], input_abc: Path) -> Path:
    paths = config.get("paths", {})
    root = Path(paths.get("aimd_root", ".")).resolve()
    hip_cfg = config.get("hipmcl", {})
    resume = config.get("resume", {})
    dry_run = bool(config.get("dry_run", False))
    skip_existing = bool(resume.get("enabled", True) and resume.get("skip_existing", True))

    output = resolve_path(hip_cfg.get("output"), root)
    log = resolve_path(hip_cfg.get("log"), root)
    assert output is not None
    ensure_dir(output.parent)

    if not hip_cfg.get("enabled", True):
        print("[RGPC] HipMCL disabled; expecting existing cluster output.")
        if not output.exists():
            raise FileNotFoundError(f"HipMCL disabled but output does not exist: {output}")
        return output

    if skip_existing and output.exists():
        print(f"[RGPC] Skip HipMCL: {output}")
        return output

    extra_args = hip_cfg.get("extra_args", "")
    if isinstance(extra_args, list):
        extra_args = " ".join(shlex.quote(str(arg)) for arg in extra_args)

    template = hip_cfg.get(
        "command_template",
        "OMP_NUM_THREADS={threads} mpirun -np {mpi_processes} {executable} -M {input} -I {inflation} "
        "-per-process-mem {per_process_mem_gb} -o {output} {extra_args}",
    )
    cmd = template.format(
        executable=resolve_external_tool("hipmcl", hip_cfg.get("executable", "auto"), root, config),
        input=str(input_abc),
        inflation=hip_cfg.get("inflation", 2.0),
        output=str(output),
        threads=hip_cfg.get("threads", 1),
        mpi_processes=hip_cfg.get("mpi_processes", 1),
        per_process_mem_gb=hip_cfg.get("per_process_mem_gb", 0),
        extra_args=str(extra_args).strip(),
    ).strip()
    run_command(cmd, log_path=log, dry_run=dry_run)
    return output


def _normalize_protein_label(label: str) -> str:
    """
    Normalize protein labels from Foldseek/HipMCL output.

    Examples:
        /path/to/TcChr01b00663.t1.pdb -> TcChr01b00663.t1
        TcChr01b00663.t1.pdb          -> TcChr01b00663.t1
        TcChr01b00663.t1.pdbqt        -> TcChr01b00663.t1
        TcChr01b00663.t1.cif          -> TcChr01b00663.t1
    """
    from pathlib import Path

    x = str(label).strip()
    x = x.strip(",;:")
    x = Path(x).name

    for suffix in [".pdbqt", ".pdb", ".cif", ".mmcif", ".ent"]:
        if x.endswith(suffix):
            x = x[: -len(suffix)]

    return x


def parse_hipmcl_output(path: Path, records: List[ProteinRecord]) -> Dict[str, List[str]]:
    all_ids: Set[str] = {r.protein_id for r in records}
    assigned: Set[str] = set()
    clusters: Dict[str, List[str]] = {}

    # Map normalized labels back to real protein_id.
    label_to_pid = {}
    for pid in all_ids:
        label_to_pid[pid] = pid
        label_to_pid[_normalize_protein_label(pid)] = pid

    unmatched_tokens = []

    with open(path, "r", encoding="utf-8") as f:
        idx = 1
        for line in f:
            line = line.strip()
            if not line:
                continue

            raw_members = [x for x in line.replace(",", "\t").split() if x]
            members = []

            for token in raw_members:
                candidates = [
                    token,
                    _normalize_protein_label(token),
                ]

                matched_pid = None
                for c in candidates:
                    if c in label_to_pid:
                        matched_pid = label_to_pid[c]
                        break

                if matched_pid is not None:
                    members.append(matched_pid)
                else:
                    unmatched_tokens.append(token)

            members = sorted(set(members))
            if not members:
                continue

            cid = f"C{idx:06d}"
            clusters[cid] = members
            assigned.update(members)
            idx += 1

    # Preserve singleton proteins not returned or not parsed from HipMCL.
    for pid in sorted(all_ids - assigned):
        cid = f"C{len(clusters) + 1:06d}"
        clusters[cid] = [pid]

    # Write a debug file next to hipmcl_output.txt.
    debug_path = Path(path).parent / "hipmcl_parse_debug.txt"
    with open(debug_path, "w", encoding="utf-8") as f:
        f.write(f"records_total\t{len(all_ids)}\n")
        f.write(f"clusters_parsed\t{len(clusters)}\n")
        f.write(f"assigned_proteins\t{len(assigned)}\n")
        f.write(f"unassigned_proteins\t{len(all_ids - assigned)}\n")
        f.write(f"unmatched_tokens\t{len(unmatched_tokens)}\n")
        f.write("\nExample unmatched tokens:\n")
        for x in unmatched_tokens[:50]:
            f.write(str(x) + "\n")

    return clusters

def create_singleton_clusters(records: List[ProteinRecord]) -> Dict[str, List[str]]:
    clusters: Dict[str, List[str]] = {}
    for idx, record in enumerate(records, start=1):
        clusters[f"C{idx:06d}"] = [record.protein_id]
    return clusters
