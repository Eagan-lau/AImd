from __future__ import annotations

import csv
import json
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

from tapocket.core.schema import FoldseekHit, TemplateRecord, MCSARecord
from tapocket.utils.fs import clean_dir, ensure_dir, symlink_or_copy


def normalize_foldseek_target_id(target: str) -> str:
    """Normalize Foldseek target names to manifest IDs.

    Foldseek may report a full path, a filename, or a database entry name.
    Staging DBs created by this package use {template_id}.pdb and {mcsa_id}.pdb,
    but this function also handles older names such as *_protein.pdb and
    *_reference_clean.pdb.
    """
    value = str(target).strip()
    if not value:
        return value
    value = Path(value).name
    for suffix in [".pdb", ".cif", ".mmcif", ".ent"]:
        if value.lower().endswith(suffix):
            value = value[: -len(suffix)]
            break
    for suffix in ["_protein", "_reference_clean"]:
        if value.endswith(suffix):
            value = value[: -len(suffix)]
    return value


def _maybe_number(value: str) -> Any:
    value = value.strip()
    if value == "":
        return value
    try:
        if any(ch in value.lower() for ch in [".", "e"]):
            return float(value)
        return int(value)
    except ValueError:
        return value


class FoldseekRunner:
    def __init__(
        self,
        binary: str = "foldseek",
        tmp_dir: str | Path = "work/tmp/foldseek",
        format_fields: list[str] | None = None,
        extra_args: list[str] | None = None,
    ):
        self.binary = binary
        self.tmp_dir = Path(tmp_dir)
        self.format_fields = format_fields or ["query", "target", "qtmscore", "ttmscore", "alntmscore"]
        self.extra_args = extra_args or []

    def search(
        self,
        query_pdb: str | Path,
        db_path: str | Path,
        output_tsv: str | Path,
        max_seqs: int = 100,
        log_json: str | Path | None = None,
    ) -> Path:
        query_pdb = Path(query_pdb).resolve()
        db_path = Path(db_path).resolve()
        output_tsv = Path(output_tsv).resolve()
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        run_tmp = self.tmp_dir / f"search_{int(time.time() * 1000)}"
        run_tmp.mkdir(parents=True, exist_ok=True)

        command = [
            self.binary,
            "easy-search",
            str(query_pdb),
            str(db_path),
            str(output_tsv),
            str(run_tmp),
            "--format-output",
            ",".join(self.format_fields),
            "--max-seqs",
            str(max_seqs),
        ] + list(self.extra_args)

        started = time.time()
        completed = subprocess.run(command, capture_output=True, text=True)
        elapsed = time.time() - started
        if log_json:
            Path(log_json).parent.mkdir(parents=True, exist_ok=True)
            Path(log_json).write_text(
                json.dumps(
                    {
                        "command": command,
                        "returncode": completed.returncode,
                        "stdout": completed.stdout,
                        "stderr": completed.stderr,
                        "elapsed_seconds": elapsed,
                    },
                    indent=2,
                    ensure_ascii=False,
                ),
                encoding="utf-8",
            )
        if completed.returncode != 0:
            raise RuntimeError(
                "Foldseek search failed.\n"
                f"Command: {' '.join(command)}\n"
                f"STDERR:\n{completed.stderr}"
            )
        return output_tsv

    def parse_hits(self, tsv_path: str | Path) -> list[FoldseekHit]:
        return parse_foldseek_hits(tsv_path, self.format_fields)


def parse_foldseek_hits(tsv_path: str | Path, fields: list[str]) -> list[FoldseekHit]:
    path = Path(tsv_path)
    hits: list[FoldseekHit] = []
    if not path.exists() or path.stat().st_size == 0:
        return hits
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for idx, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            raw: dict[str, Any] = {}
            for i, field in enumerate(fields):
                raw[field] = _maybe_number(parts[i]) if i < len(parts) else ""
            query = str(raw.get("query", ""))
            target = str(raw.get("target", ""))
            hit = FoldseekHit(query=query, target=target, raw=raw, rank=idx)
            hit.normalized_target_id = normalize_foldseek_target_id(target)
            hit.compute_score()
            hits.append(hit)
    return hits


def select_hits(
    hits: list[FoldseekHit],
    min_qtmscore: float = 0.45,
    min_ttmscore: float = 0.45,
    min_alntmscore: float = 0.45,
    keep_top_k: int = 10,
    deduplicate_target: bool = True,
) -> list[FoldseekHit]:
    selected: list[FoldseekHit] = []
    for hit in hits:
        hit.compute_score()
        if (
            hit.qtmscore >= min_qtmscore
            or hit.ttmscore >= min_ttmscore
            or hit.alntmscore >= min_alntmscore
        ):
            selected.append(hit)
    selected.sort(key=lambda h: (h.score, h.bits, -h.evalue), reverse=True)

    if deduplicate_target:
        best_by_target: dict[str, FoldseekHit] = {}
        for hit in selected:
            key = hit.normalized_target_id or hit.target
            if key not in best_by_target:
                best_by_target[key] = hit
        selected = list(best_by_target.values())

    selected = selected[:keep_top_k]
    for rank, hit in enumerate(selected, start=1):
        hit.rank = rank
    return selected

def write_hits_tsv(hits: list[FoldseekHit], output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["rank", "query", "target", "normalized_target_id", "score", "qtmscore", "ttmscore", "alntmscore", "evalue", "bits"]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for hit in hits:
            writer.writerow(
                {
                    "rank": hit.rank,
                    "query": hit.query,
                    "target": hit.target,
                    "normalized_target_id": hit.normalized_target_id,
                    "score": hit.score,
                    "qtmscore": hit.qtmscore,
                    "ttmscore": hit.ttmscore,
                    "alntmscore": hit.alntmscore,
                    "evalue": hit.evalue,
                    "bits": hit.bits,
                }
            )


def _run_command(command: list[str], log_json: str | Path | None = None) -> None:
    started = time.time()
    completed = subprocess.run(command, capture_output=True, text=True)
    elapsed = time.time() - started
    if log_json:
        Path(log_json).parent.mkdir(parents=True, exist_ok=True)
        Path(log_json).write_text(
            json.dumps(
                {
                    "command": command,
                    "returncode": completed.returncode,
                    "stdout": completed.stdout,
                    "stderr": completed.stderr,
                    "elapsed_seconds": elapsed,
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(command)}\nSTDERR:\n{completed.stderr}")


def build_foldseek_db_from_template_records(
    records: list[TemplateRecord],
    root: str | Path,
    staging_dir: str | Path,
    output_db: str | Path,
    binary: str = "foldseek",
    create_index: bool = False,
    force: bool = False,
) -> dict[str, Any]:
    root = Path(root).resolve()
    staging_dir = Path(staging_dir).resolve()
    output_db = Path(output_db).resolve()
    if force:
        clean_dir(staging_dir)
    else:
        ensure_dir(staging_dir)
    output_db.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    for record in records:
        src = root / record.protein_path if not Path(record.protein_path).is_absolute() else Path(record.protein_path)
        if not src.exists():
            continue
        dst = staging_dir / f"{record.template_id}.pdb"
        symlink_or_copy(src, dst, overwrite=True)
        count += 1

    _run_command([binary, "createdb", str(staging_dir), str(output_db)], log_json=output_db.parent / "createdb_template_log.json")
    if create_index:
        tmp_dir = output_db.parent / "tmp_createindex_template"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        _run_command([binary, "createindex", str(output_db), str(tmp_dir)], log_json=output_db.parent / "createindex_template_log.json")
    return {"db": str(output_db), "staging_dir": str(staging_dir), "entries": count, "create_index": create_index}


def build_foldseek_db_from_mcsa_records(
    records: list[MCSARecord],
    root: str | Path,
    staging_dir: str | Path,
    output_db: str | Path,
    binary: str = "foldseek",
    create_index: bool = False,
    force: bool = False,
) -> dict[str, Any]:
    root = Path(root).resolve()
    staging_dir = Path(staging_dir).resolve()
    output_db = Path(output_db).resolve()
    if force:
        clean_dir(staging_dir)
    else:
        ensure_dir(staging_dir)
    output_db.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    for record in records:
        src = root / record.reference_protein_path if not Path(record.reference_protein_path).is_absolute() else Path(record.reference_protein_path)
        if not src.exists():
            continue
        dst = staging_dir / f"{record.mcsa_id}.pdb"
        symlink_or_copy(src, dst, overwrite=True)
        count += 1

    _run_command([binary, "createdb", str(staging_dir), str(output_db)], log_json=output_db.parent / "createdb_mcsa_log.json")
    if create_index:
        tmp_dir = output_db.parent / "tmp_createindex_mcsa"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        _run_command([binary, "createindex", str(output_db), str(tmp_dir)], log_json=output_db.parent / "createindex_mcsa_log.json")
    return {"db": str(output_db), "staging_dir": str(staging_dir), "entries": count, "create_index": create_index}
