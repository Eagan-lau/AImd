from __future__ import annotations

from pathlib import Path
from typing import Any

from tapocket.core.config import TapocketConfig
from tapocket.core.schema import TemplateRecord, MCSARecord
from tapocket.utils.fs import rel_or_abs
from tapocket.utils.jsonl import write_jsonl


def _first_or_none(paths: list[Path]) -> Path | None:
    return paths[0] if paths else None


def build_template_manifest(config: TapocketConfig) -> list[TemplateRecord]:
    root = config.root
    template_root = config.path("paths", "template_root")
    group_pattern = config.get("template_db", "group_pattern", default="protein_group_*")
    protein_suffix = config.get("template_db", "protein_suffix", default="_protein.pdb")
    pocket_pattern = config.get("template_db", "pocket_pattern", default="*_pocket-*.pdb")
    union_pocket_pattern = config.get("template_db", "union_pocket_pattern", default="*_pocket_union.pdb")
    use_union = bool(config.get("template_db", "use_union_pocket", default=False))
    union_as_fallback = bool(config.get("template_db", "union_pocket_as_fallback", default=True))

    records: list[TemplateRecord] = []
    if not template_root.exists():
        raise FileNotFoundError(f"Template root does not exist: {template_root}")

    for group_dir in sorted(p for p in template_root.glob(group_pattern) if p.is_dir()):
        template_id = group_dir.name
        protein_candidates = sorted(group_dir.glob(f"*{protein_suffix}"))
        if not protein_candidates:
            # Some datasets may use exact template_id_protein.pdb but suffix matching should catch it.
            continue
        protein_path = protein_candidates[0]
        pocket_paths = sorted(group_dir.glob(pocket_pattern))
        pocket_paths = [p for p in pocket_paths if "pocket_union" not in p.name]
        union_pocket = _first_or_none(sorted(group_dir.glob(union_pocket_pattern)))
        if use_union and union_pocket:
            pocket_paths = [union_pocket]
        elif not pocket_paths and union_as_fallback and union_pocket:
            pocket_paths = [union_pocket]

        records.append(
            TemplateRecord(
                template_id=template_id,
                group_dir=rel_or_abs(group_dir, root),
                protein_path=rel_or_abs(protein_path, root),
                pocket_paths=[rel_or_abs(p, root) for p in pocket_paths],
                union_pocket_path=rel_or_abs(union_pocket, root) if union_pocket else None,
                source_sample_ids_path=rel_or_abs(group_dir / "source_sample_ids.txt", root)
                if (group_dir / "source_sample_ids.txt").exists()
                else None,
                source_samples_csv_path=rel_or_abs(group_dir / "source_samples.csv", root)
                if (group_dir / "source_samples.csv").exists()
                else None,
            )
        )
    return records


def build_mcsa_manifest(config: TapocketConfig) -> list[MCSARecord]:
    root = config.root
    mcsa_root = config.path("paths", "mcsa_root")
    group_pattern = config.get("mcsa_db", "group_pattern", default="*")
    reference_pattern = config.get("mcsa_db", "reference_pattern", default="*_reference_clean.pdb")
    active_site_pattern = config.get("mcsa_db", "active_site_pattern", default="*_mcsa_active_site_all.pdb")

    records: list[MCSARecord] = []
    if not mcsa_root.exists():
        raise FileNotFoundError(f"M-CSA root does not exist: {mcsa_root}")

    for group_dir in sorted(p for p in mcsa_root.glob(group_pattern) if p.is_dir()):
        reference = _first_or_none(sorted(group_dir.glob(reference_pattern)))
        active = _first_or_none(sorted(group_dir.glob(active_site_pattern)))
        if not reference or not active:
            continue
        mcsa_id = reference.name.replace("_reference_clean.pdb", "")
        records.append(
            MCSARecord(
                mcsa_id=mcsa_id,
                group_dir=rel_or_abs(group_dir, root),
                reference_protein_path=rel_or_abs(reference, root),
                active_site_path=rel_or_abs(active, root),
            )
        )
    return records


def write_manifests(config: TapocketConfig) -> dict[str, Any]:
    template_records = build_template_manifest(config)
    mcsa_records = build_mcsa_manifest(config)
    template_manifest = config.path("paths", "template_manifest")
    mcsa_manifest = config.path("paths", "mcsa_manifest")

    write_jsonl(template_manifest, [r.to_dict() for r in template_records])
    write_jsonl(mcsa_manifest, [r.to_dict() for r in mcsa_records])

    return {
        "template_manifest": str(template_manifest),
        "template_count": len(template_records),
        "template_with_pockets": sum(1 for r in template_records if r.pocket_paths),
        "mcsa_manifest": str(mcsa_manifest),
        "mcsa_count": len(mcsa_records),
    }


def check_database_layout(config: TapocketConfig) -> dict[str, Any]:
    template_root = config.path("paths", "template_root")
    mcsa_root = config.path("paths", "mcsa_root")

    result: dict[str, Any] = {
        "project_root": str(config.root),
        "template_root": str(template_root),
        "mcsa_root": str(mcsa_root),
        "template_root_exists": template_root.exists(),
        "mcsa_root_exists": mcsa_root.exists(),
    }

    if template_root.exists():
        template_records = build_template_manifest(config)
        result.update(
            {
                "template_records_found": len(template_records),
                "template_records_with_pockets": sum(1 for r in template_records if r.pocket_paths),
                "template_records_without_pockets": sum(1 for r in template_records if not r.pocket_paths),
            }
        )
    if mcsa_root.exists():
        mcsa_records = build_mcsa_manifest(config)
        result.update({"mcsa_records_found": len(mcsa_records)})
    return result
