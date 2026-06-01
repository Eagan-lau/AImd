#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict, Any, List

from .io import ProteinRecord
from .utils import ensure_dir, resolve_path, run_command, resolve_external_tool


def _build_search_extra_args(foldseek_cfg: Dict[str, Any]) -> List[str]:
    args: List[str] = []
    search_cfg = foldseek_cfg.get("search", {})
    sensitivity = search_cfg.get("sensitivity")
    if sensitivity is not None:
        args.extend(["-s", str(sensitivity)])
    args.extend([str(x) for x in search_cfg.get("extra_args", [])])
    require_backtrace = bool(search_cfg.get("require_backtrace", True))
    if require_backtrace and "-a" not in args:
        args.insert(0, "-a")
    return args


def run_foldseek_all_vs_all(config: Dict[str, Any], records: List[ProteinRecord]) -> Path:
    paths = config.get("paths", {})
    root = Path(paths.get("aimd_root", ".")).resolve()
    foldseek_cfg = config.get("foldseek", {})
    resume = config.get("resume", {})
    dry_run = bool(config.get("dry_run", False))
    skip_existing = bool(resume.get("enabled", True) and resume.get("skip_existing", True))

    executable = resolve_external_tool("foldseek", foldseek_cfg.get("executable", "auto"), root, config)
    structure_dir = resolve_path(paths.get("structure_dir"), root)
    db_prefix = resolve_path(foldseek_cfg.get("createdb", {}).get("db_prefix"), root)
    result_db = resolve_path(foldseek_cfg.get("search", {}).get("result_db"), root)
    raw_tsv = resolve_path(foldseek_cfg.get("search", {}).get("raw_tsv"), root)
    tmp_dir = resolve_path(foldseek_cfg.get("tmp_dir"), root)
    threads = int(foldseek_cfg.get("threads", 1))
    fmt = foldseek_cfg.get("convertalis", {}).get("format_output", "query,target,qtmscore,ttmscore")

    assert structure_dir is not None and db_prefix is not None and result_db is not None and raw_tsv is not None and tmp_dir is not None
    ensure_dir(db_prefix.parent)
    ensure_dir(result_db.parent)
    ensure_dir(raw_tsv.parent)
    ensure_dir(tmp_dir)

    createdb_log = db_prefix.parent / "createdb.log"
    search_log = result_db.parent / "foldseek_search.log"
    convert_log = result_db.parent / "convertalis.log"

    db_marker = Path(str(db_prefix) + ".dbtype")
    if not skip_existing or not db_marker.exists():
        cmd = [executable, "createdb", str(structure_dir), str(db_prefix)]
        run_command(cmd, log_path=createdb_log, dry_run=dry_run)
    else:
        print(f"[RGPC] Skip Foldseek createdb: {db_prefix}")

    if not skip_existing or not Path(str(result_db) + ".dbtype").exists():
        cmd = [executable, "search", str(db_prefix), str(db_prefix), str(result_db), str(tmp_dir), "--threads", str(threads)]
        cmd.extend(_build_search_extra_args(foldseek_cfg))
        run_command(cmd, log_path=search_log, dry_run=dry_run)
    else:
        print(f"[RGPC] Skip Foldseek search: {result_db}")

    if not skip_existing or not raw_tsv.exists():
        cmd = [executable, "convertalis", str(db_prefix), str(db_prefix), str(result_db), str(raw_tsv), "--format-output", fmt]
        run_command(cmd, log_path=convert_log, dry_run=dry_run)
    else:
        print(f"[RGPC] Skip Foldseek convertalis: {raw_tsv}")

    return raw_tsv


def run_foldseek_per_query(config: Dict[str, Any], records: List[ProteinRecord]) -> Path:
    """Run Foldseek per query and merge TSV outputs.

    This mode is slower but easier to resume for very large jobs. Each query PDB is
    copied to a temporary single-query directory, converted to a query DB, searched
    against the whole target DB, converted to TSV, and finally merged.
    """
    paths = config.get("paths", {})
    root = Path(paths.get("aimd_root", ".")).resolve()
    foldseek_cfg = config.get("foldseek", {})
    resume = config.get("resume", {})
    dry_run = bool(config.get("dry_run", False))
    skip_existing = bool(resume.get("enabled", True) and resume.get("skip_existing", True))

    executable = resolve_external_tool("foldseek", foldseek_cfg.get("executable", "auto"), root, config)
    structure_dir = resolve_path(paths.get("structure_dir"), root)
    db_prefix = resolve_path(foldseek_cfg.get("createdb", {}).get("db_prefix"), root)
    raw_tsv = resolve_path(foldseek_cfg.get("search", {}).get("raw_tsv"), root)
    tmp_dir = resolve_path(foldseek_cfg.get("tmp_dir"), root)
    threads = int(foldseek_cfg.get("threads", 1))
    fmt = foldseek_cfg.get("convertalis", {}).get("format_output", "query,target,qtmscore,ttmscore")

    assert structure_dir is not None and db_prefix is not None and raw_tsv is not None and tmp_dir is not None
    ensure_dir(raw_tsv.parent)
    ensure_dir(tmp_dir)

    db_marker = Path(str(db_prefix) + ".dbtype")
    if not skip_existing or not db_marker.exists():
        run_command([executable, "createdb", str(structure_dir), str(db_prefix)], log_path=db_prefix.parent / "createdb.log", dry_run=dry_run)

    per_query_dir = raw_tsv.parent / "per_query"
    ensure_dir(per_query_dir)

    tsv_paths: List[Path] = []
    for r in records:
        qdir = per_query_dir / r.protein_id
        ensure_dir(qdir)
        q_struct_dir = qdir / "query_struct"
        ensure_dir(q_struct_dir)
        q_copy = q_struct_dir / r.pdb_path.name
        if not q_copy.exists():
            shutil.copy2(r.pdb_path, q_copy)
        q_db = qdir / "query_db"
        q_result = qdir / "result_db"
        q_tsv = qdir / f"{r.protein_id}.tsv"
        tsv_paths.append(q_tsv)

        if not skip_existing or not Path(str(q_db) + ".dbtype").exists():
            run_command([executable, "createdb", str(q_struct_dir), str(q_db)], log_path=qdir / "createdb.log", dry_run=dry_run)
        if not skip_existing or not Path(str(q_result) + ".dbtype").exists():
            cmd = [executable, "search", str(q_db), str(db_prefix), str(q_result), str(tmp_dir / r.protein_id), "--threads", str(threads)]
            cmd.extend(_build_search_extra_args(foldseek_cfg))
            run_command(cmd, log_path=qdir / "search.log", dry_run=dry_run)
        if not skip_existing or not q_tsv.exists():
            run_command([executable, "convertalis", str(q_db), str(db_prefix), str(q_result), str(q_tsv), "--format-output", fmt], log_path=qdir / "convertalis.log", dry_run=dry_run)

    if not skip_existing or not raw_tsv.exists():
        with open(raw_tsv, "w", encoding="utf-8") as out:
            for p in tsv_paths:
                if not p.exists():
                    continue
                with open(p, "r", encoding="utf-8") as f:
                    for line in f:
                        if line.strip():
                            out.write(line)
    return raw_tsv


def run_foldseek(config: Dict[str, Any], records: List[ProteinRecord]) -> Path:
    mode = config.get("foldseek", {}).get("mode", "all_vs_all")
    if mode == "all_vs_all":
        return run_foldseek_all_vs_all(config, records)
    if mode == "per_query":
        return run_foldseek_per_query(config, records)
    raise ValueError(f"Unsupported foldseek.mode: {mode}")
