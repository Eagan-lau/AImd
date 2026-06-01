from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from tapocket.core.schema import CandidatePocket, FoldseekHit, PipelineRunSummary
from tapocket.utils.fs import resolve_path
from tapocket.utils.pdb_utils import write_multimodel_pdb


def write_json(path: str | Path, data: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def export_hits_json(path: str | Path, hits: list[FoldseekHit]) -> None:
    write_json(path, [hit.to_dict() for hit in hits])


def export_candidates_json(path: str | Path, candidates: list[CandidatePocket]) -> None:
    write_json(path, [candidate.to_dict() for candidate in candidates])


def export_query_pocket_residues(path: str | Path, candidates: list[CandidatePocket]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pocket_id",
        "query_id",
        "source",
        "template_id",
        "template_pocket_id",
        "chain",
        "resi",
        "icode",
        "resn",
        "min_distance_to_mapped_pocket",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for candidate in candidates:
            for residue in candidate.query_pocket_residues:
                row = residue.to_dict()
                row.update({"pocket_id": candidate.pocket_id, "source": candidate.source})
                writer.writerow(row)


def export_boxes_tsv(path: str | Path, candidates: list[CandidatePocket]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pocket_id", "source", "center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
        "padding_angstrom", "final_score", "mcsa_matched", "ai_confidence",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for c in candidates:
            center = c.box.get("center", ["", "", ""])
            size = c.box.get("size", ["", "", ""])
            mcsa = c.mcsa_support or {}
            ai = c.ai_support or {}
            writer.writerow({
                "pocket_id": c.pocket_id,
                "source": c.source,
                "center_x": center[0] if len(center) > 0 else "",
                "center_y": center[1] if len(center) > 1 else "",
                "center_z": center[2] if len(center) > 2 else "",
                "size_x": size[0] if len(size) > 0 else "",
                "size_y": size[1] if len(size) > 1 else "",
                "size_z": size[2] if len(size) > 2 else "",
                "padding_angstrom": c.box.get("padding_angstrom", ""),
                "final_score": c.final_score,
                "mcsa_matched": mcsa.get("matched", ""),
                "ai_confidence": ai.get("confidence", ""),
            })


def export_summary_tsv(path: str | Path, candidates: list[CandidatePocket]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pocket_id",
        "query_id",
        "source",
        "final_score",
        "query_residue_count",
        "mapping_coverage",
        "mapping_method",
        "supporting_pocket_count",
        "supporting_template_count",
        "mcsa_matched",
        "mcsa_match_basis",
        "mcsa_atoms_in_box",
        "mcsa_fraction_in_box",
        "mcsa_min_distance_angstrom",
        "box_center",
        "box_size",
        "mapped_pocket_path",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for c in candidates:
            mcsa = c.mcsa_support or {}
            writer.writerow(
                {
                    "pocket_id": c.pocket_id,
                    "query_id": c.query_id,
                    "source": c.source,
                    "final_score": c.final_score,
                    "query_residue_count": len(c.query_pocket_residues),
                    "mapping_coverage": c.mapping_quality.get("mapping_coverage", ""),
                    "mapping_method": c.mapping_quality.get("method", ""),
                    "supporting_pocket_count": c.mapping_quality.get("supporting_pocket_count", 1),
                    "supporting_template_count": c.mapping_quality.get("supporting_template_count", 1),
                    "mcsa_matched": mcsa.get("matched", ""),
                    "mcsa_match_basis": mcsa.get("match_basis", ""),
                    "mcsa_atoms_in_box": mcsa.get("active_site_atoms_in_box", ""),
                    "mcsa_fraction_in_box": mcsa.get("active_site_fraction_in_box", ""),
                    "mcsa_min_distance_angstrom": mcsa.get("min_distance_angstrom", ""),
                    "box_center": c.box.get("center", ""),
                    "box_size": c.box.get("size", ""),
                    "mapped_pocket_path": c.mapped_pocket_path,
                }
            )


def export_final_pdb(path: str | Path, candidates: list[CandidatePocket], run_root: str | Path) -> None:
    path = Path(path)
    mapped_paths = []
    for c in candidates:
        if not c.mapped_pocket_path:
            continue
        p = resolve_path(c.mapped_pocket_path, run_root)
        if p.exists():
            mapped_paths.append(p)
    if mapped_paths:
        write_multimodel_pdb(mapped_paths, path)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("END\n", encoding="utf-8")


def export_run_summary(path: str | Path, summary: PipelineRunSummary) -> None:
    write_json(path, summary.to_dict())
