from __future__ import annotations

from pathlib import Path
from typing import Any

from tapocket.core.schema import CandidatePocket, FoldseekHit, TemplateRecord
from tapocket.databases.template_db import TemplateDB
from tapocket.mapping.aligner import get_aligner
from tapocket.utils.fs import rel_or_abs
from tapocket.utils.pdb_utils import (
    atom_box,
    box_center,
    box_size,
    extract_query_residues_by_ca_distance,
    extract_query_residues_by_heavy_atom_distance,
    extract_query_residues_near_pocket,
    parse_pdb_atoms,
)


def _template_pocket_id(template_id: str, pocket_path: Path) -> str:
    stem = pocket_path.stem
    prefix = f"{template_id}_"
    if stem.startswith(prefix):
        return stem[len(prefix):]
    return stem


def _box_dict(mapped_pocket_path: str | Path, padding: float) -> dict[str, Any]:
    atoms = parse_pdb_atoms(mapped_pocket_path, include_hydrogen=False, include_hetatm=True)
    if not atoms:
        return {}
    box = atom_box(atoms, padding=padding)
    center = box_center(box)
    size = box_size(box)
    return {
        "method": "mapped_pocket_atom_box",
        "padding_angstrom": padding,
        "center": [round(v, 4) for v in center],
        "size": [round(v, 4) for v in size],
        "xmin": round(box["xmin"], 4),
        "xmax": round(box["xmax"], 4),
        "ymin": round(box["ymin"], 4),
        "ymax": round(box["ymax"], 4),
        "zmin": round(box["zmin"], 4),
        "zmax": round(box["zmax"], 4),
    }


class PocketMapper:
    def __init__(self, config: Any, template_db: TemplateDB, query_id: str, run_root: Path):
        self.config = config
        self.template_db = template_db
        self.query_id = query_id
        self.run_root = Path(run_root).resolve()
        aligner_name = config.get("pocket_mapping", "aligner", default="pymol")
        save_sessions = bool(config.get("pocket_mapping", "save_alignment_sessions", default=False))
        self.aligner = get_aligner(aligner_name, save_sessions=save_sessions)
        self.extraction_method = str(config.get("pocket_mapping", "residue_extraction", "method", default="ca_distance"))
        self.ca_cutoff = float(config.get("pocket_mapping", "residue_extraction", "ca_distance_cutoff_angstrom", default=2.0))
        self.cutoff = float(config.get("pocket_mapping", "residue_extraction", "distance_cutoff_angstrom", default=4.5))
        self.include_h = bool(config.get("pocket_mapping", "residue_extraction", "include_hydrogen", default=False))
        self.qc_enabled = bool(config.get("pocket_mapping", "quality_control", "enabled", default=True))
        self.min_mapping_coverage = float(config.get("pocket_mapping", "quality_control", "min_mapping_coverage", default=0.40))
        self.min_query_residue_count = int(config.get("pocket_mapping", "quality_control", "min_query_residue_count", default=5))
        self.low_quality_action = str(config.get("pocket_mapping", "quality_control", "low_quality_action", default="drop"))
        self.box_padding = float(config.get("docking_box", "padding_angstrom", default=2.0))

    def map_template_hit(
        self,
        query_pdb: str | Path,
        hit: FoldseekHit,
        record: TemplateRecord,
        output_dir: str | Path,
    ) -> list[CandidatePocket]:
        template_protein = self.template_db.protein_path(record)
        pocket_paths = self.template_db.pocket_paths(record)
        if not pocket_paths:
            return []

        output_dir = Path(output_dir).resolve()
        mapped_dir = output_dir / "mapped_pockets"
        sessions_dir = output_dir / "sessions"
        prefix = f"{self.query_id}_{record.template_id}"
        session_path = sessions_dir / f"{prefix}.pse"

        mapping_result = self.aligner.map_objects_to_query(
            query_pdb=query_pdb,
            reference_pdb=template_protein,
            object_paths=pocket_paths,
            output_dir=mapped_dir,
            output_prefix=prefix,
            session_path=session_path,
        )

        candidates: list[CandidatePocket] = []
        for pocket_path, mapped_path in zip(pocket_paths, mapping_result["mapped_outputs"]):
            template_pocket_id = _template_pocket_id(record.template_id, Path(pocket_path))
            if self.extraction_method == "ca_distance":
                residues, mapping_quality = extract_query_residues_by_ca_distance(
                    query_pdb=query_pdb,
                    mapped_pocket_pdb=mapped_path,
                    query_id=self.query_id,
                    template_id=record.template_id,
                    template_pocket_id=template_pocket_id,
                    cutoff_angstrom=self.ca_cutoff,
                )
            elif self.extraction_method == "heavy_atom_distance":
                residues, mapping_quality = extract_query_residues_by_heavy_atom_distance(
                    query_pdb=query_pdb,
                    mapped_pocket_pdb=mapped_path,
                    query_id=self.query_id,
                    template_id=record.template_id,
                    template_pocket_id=template_pocket_id,
                    cutoff_angstrom=self.cutoff,
                    include_hydrogen=self.include_h,
                )
            else:
                residues = extract_query_residues_near_pocket(
                    query_pdb=query_pdb,
                    mapped_pocket_pdb=mapped_path,
                    query_id=self.query_id,
                    template_id=record.template_id,
                    template_pocket_id=template_pocket_id,
                    cutoff_angstrom=self.cutoff,
                    include_hydrogen=self.include_h,
                )
                mapping_quality = {
                    "method": self.extraction_method,
                    "query_residue_count": len(residues),
                    "mapping_coverage": 1.0 if residues else 0.0,
                }

            mapping_coverage = float(mapping_quality.get("mapping_coverage", 0.0) or 0.0)
            query_residue_count = int(mapping_quality.get("query_residue_count", len(residues)) or 0)
            mapping_pass = True
            drop_reason = None
            if self.qc_enabled:
                mapping_pass = (
                    mapping_coverage >= self.min_mapping_coverage
                    and query_residue_count >= self.min_query_residue_count
                )
                if not mapping_pass:
                    drop_reason = "low_mapping_quality"
            mapping_quality["quality_control_enabled"] = self.qc_enabled
            mapping_quality["min_mapping_coverage"] = self.min_mapping_coverage
            mapping_quality["min_query_residue_count"] = self.min_query_residue_count
            mapping_quality["mapping_pass"] = mapping_pass
            mapping_quality["low_quality_action"] = self.low_quality_action

            cand = CandidatePocket(
                query_id=self.query_id,
                template_id=record.template_id,
                template_pocket_id=template_pocket_id,
                source="foldseek_template",
                mapped_pocket_path=rel_or_abs(mapped_path, self.run_root),
                template_protein_path=rel_or_abs(template_protein, self.config.root),
                template_pocket_path=rel_or_abs(pocket_path, self.config.root),
                foldseek_hit=hit.to_dict(),
                query_pocket_residues=residues,
                mapping_quality=mapping_quality,
                box=_box_dict(mapped_path, self.box_padding),
                final_score=hit.score,
                kept=mapping_pass or self.low_quality_action != "drop",
                drop_reason=drop_reason,
            )
            candidates.append(cand)
        return candidates
