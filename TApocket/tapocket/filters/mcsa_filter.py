from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from tapocket.core.schema import CandidatePocket, FoldseekHit
from tapocket.databases.mcsa_db import MCSADB
from tapocket.mapping.aligner import get_aligner
from tapocket.utils.fs import rel_or_abs, resolve_path
from tapocket.utils.pdb_utils import (
    parse_pdb_atoms,
    min_distance_between_atom_sets,
    count_atom_pairs_within_cutoff,
    atom_box,
    count_atoms_in_box,
    box_center,
    box_size,
)


class MCSAFilter:
    def __init__(self, config: Any, mcsa_db: MCSADB, query_id: str, run_root: Path):
        self.config = config
        self.mcsa_db = mcsa_db
        self.query_id = query_id
        self.run_root = Path(run_root).resolve()
        aligner_name = config.get("pocket_mapping", "aligner", default="pymol")
        save_sessions = bool(config.get("pocket_mapping", "save_alignment_sessions", default=False))
        self.aligner = get_aligner(aligner_name, save_sessions=save_sessions)
        self.match_method = str(config.get("mcsa", "active_site_match", "method", default="box_overlap"))
        self.cutoff = float(config.get("mcsa", "active_site_match", "distance_cutoff_angstrom", default=4.5))
        self.min_matched_atoms = int(config.get("mcsa", "active_site_match", "min_matched_atoms", default=5))
        self.box_padding = float(config.get("mcsa", "active_site_match", "box_padding_angstrom", default=2.0))
        self.min_active_site_atoms_in_box = int(config.get("mcsa", "active_site_match", "min_active_site_atoms_in_box", default=1))
        self.min_active_site_fraction_in_box = float(config.get("mcsa", "active_site_match", "min_active_site_fraction_in_box", default=0.10))

    def filter_candidates(
        self,
        query_pdb: str | Path,
        candidates: list[CandidatePocket],
        mcsa_hits: list[FoldseekHit],
        output_dir: str | Path,
    ) -> tuple[list[CandidatePocket], dict[str, Any]]:
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        mapped_active_dir = output_dir / "mapped_active_sites"
        mapped_active_dir.mkdir(parents=True, exist_ok=True)

        report: dict[str, Any] = {
            "mcsa_enabled": True,
            "mcsa_hit_count": len(mcsa_hits),
            "mapped_active_sites": [],
            "matches": [],
            "active_site_in_candidate_pocket": False,
            "filter_action": "keep_all_pockets",
        }

        if not mcsa_hits:
            report["filter_action"] = self.config.get("mcsa", "filtering_policy", "if_no_mcsa_hit", default="keep_all_pockets")
            self._write_outputs(output_dir, report)
            return candidates, report

        mapped_active_paths: list[tuple[FoldseekHit, str]] = []
        for hit in mcsa_hits:
            record = self.mcsa_db.get(hit.normalized_target_id or hit.target)
            if not record:
                report["matches"].append({"mcsa_target": hit.target, "error": "not_found_in_mcsa_manifest"})
                continue
            reference_path = self.mcsa_db.reference_path(record)
            active_site_path = self.mcsa_db.active_site_path(record)
            prefix = f"{self.query_id}_{record.mcsa_id}"
            try:
                mapping = self.aligner.map_objects_to_query(
                    query_pdb=query_pdb,
                    reference_pdb=reference_path,
                    object_paths=[active_site_path],
                    output_dir=mapped_active_dir,
                    output_prefix=prefix,
                    session_path=output_dir / "sessions" / f"{prefix}.pse",
                )
            except Exception as exc:
                report["matches"].append({"mcsa_id": record.mcsa_id, "error": str(exc)})
                continue
            if not mapping["mapped_outputs"]:
                continue
            mapped_active = mapping["mapped_outputs"][0]
            mapped_active_paths.append((hit, mapped_active))
            report["mapped_active_sites"].append(
                {
                    "mcsa_id": record.mcsa_id,
                    "mcsa_hit": hit.to_dict(),
                    "reference_protein_path": rel_or_abs(reference_path, self.config.root),
                    "active_site_path": rel_or_abs(active_site_path, self.config.root),
                    "mapped_active_site_path": rel_or_abs(mapped_active, self.run_root),
                }
            )

        matched_pocket_ids: set[str] = set()
        for candidate in candidates:
            pocket_path = resolve_path(candidate.mapped_pocket_path, self.run_root)
            pocket_atoms = parse_pdb_atoms(pocket_path, include_hydrogen=False, include_hetatm=True)
            best_distance = float("inf")
            best_mcsa_id = None
            total_matched_atoms = 0
            best_atoms_in_box = 0
            best_fraction_in_box = 0.0
            best_box = None
            if pocket_atoms:
                pocket_box = atom_box(pocket_atoms, padding=self.box_padding)
            else:
                pocket_box = None

            for hit, mapped_active in mapped_active_paths:
                active_atoms = parse_pdb_atoms(mapped_active, include_hydrogen=False, include_hetatm=True)
                min_dist = min_distance_between_atom_sets(active_atoms, pocket_atoms)
                matched_atoms = count_atom_pairs_within_cutoff(active_atoms, pocket_atoms, self.cutoff)
                atoms_in_box = count_atoms_in_box(active_atoms, pocket_box) if pocket_box else 0
                fraction_in_box = (atoms_in_box / len(active_atoms)) if active_atoms else 0.0

                if min_dist < best_distance:
                    best_distance = min_dist
                    best_mcsa_id = hit.normalized_target_id or hit.target
                if atoms_in_box > best_atoms_in_box or fraction_in_box > best_fraction_in_box:
                    best_atoms_in_box = atoms_in_box
                    best_fraction_in_box = fraction_in_box
                    best_box = pocket_box
                total_matched_atoms += matched_atoms

            if self.match_method == "box_overlap":
                matched = (
                    best_atoms_in_box >= self.min_active_site_atoms_in_box
                    and best_fraction_in_box >= self.min_active_site_fraction_in_box
                )
                match_basis = "box_overlap"
            elif self.match_method == "box_or_distance":
                box_matched = (
                    best_atoms_in_box >= self.min_active_site_atoms_in_box
                    and best_fraction_in_box >= self.min_active_site_fraction_in_box
                )
                distance_matched = best_distance <= self.cutoff and total_matched_atoms >= self.min_matched_atoms
                matched = box_matched or distance_matched
                match_basis = "box_overlap" if box_matched else ("atom_distance" if distance_matched else "none")
            else:
                matched = best_distance <= self.cutoff and total_matched_atoms >= self.min_matched_atoms
                match_basis = "atom_distance"

            box_center_value = None
            box_size_value = None
            if best_box:
                box_center_value = [round(v, 4) for v in box_center(best_box)]
                box_size_value = [round(v, 4) for v in box_size(best_box)]

            match_record = {
                "query_id": self.query_id,
                "pocket_id": candidate.pocket_id,
                "template_id": candidate.template_id,
                "template_pocket_id": candidate.template_pocket_id,
                "matched": matched,
                "match_method": self.match_method,
                "match_basis": match_basis,
                "best_mcsa_id": best_mcsa_id,
                "min_distance_angstrom": None if best_distance == float("inf") else round(best_distance, 4),
                "matched_atom_pairs": total_matched_atoms,
                "active_site_atoms_in_box": best_atoms_in_box,
                "active_site_fraction_in_box": round(best_fraction_in_box, 6),
                "box_padding_angstrom": self.box_padding,
                "pocket_box_center": box_center_value,
                "pocket_box_size": box_size_value,
            }
            report["matches"].append(match_record)
            candidate.mcsa_support = match_record
            if matched:
                matched_pocket_ids.add(candidate.pocket_id)

        if matched_pocket_ids:
            report["active_site_in_candidate_pocket"] = True
            action = self.config.get("mcsa", "filtering_policy", "if_active_site_in_pocket", default="keep_matching_pockets_only")
            report["filter_action"] = action
            if action == "keep_matching_pockets_only":
                filtered: list[CandidatePocket] = []
                for candidate in candidates:
                    if candidate.pocket_id in matched_pocket_ids:
                        candidate.kept = True
                        candidate.drop_reason = None
                        filtered.append(candidate)
                    else:
                        candidate.kept = False
                        candidate.drop_reason = "mcsa_active_site_not_in_pocket"
                self._write_outputs(output_dir, report)
                return filtered, report
        else:
            report["active_site_in_candidate_pocket"] = False
            report["filter_action"] = self.config.get("mcsa", "filtering_policy", "if_no_active_site_in_pocket", default="keep_all_pockets")

        self._write_outputs(output_dir, report)
        return candidates, report

    def _write_outputs(self, output_dir: Path, report: dict[str, Any]) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "mcsa_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        matches = report.get("matches", [])
        if matches:
            all_fields = sorted({key for row in matches if isinstance(row, dict) for key in row.keys()})
            with (output_dir / "mcsa_pocket_matches.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=all_fields, delimiter="\t")
                writer.writeheader()
                for row in matches:
                    if isinstance(row, dict):
                        writer.writerow(row)
