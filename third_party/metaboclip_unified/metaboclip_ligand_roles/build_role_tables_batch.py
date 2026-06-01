#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from metaboclip_ligand_roles.annotator_core import annotate_ligand, write_annotation
from metaboclip_ligand_roles.atom_map import recover_atom_map_by_coordinates, write_atom_map
from metaboclip_ligand_roles.molio import iter_ligand_sources
from metaboclip_ligand_roles.role_table import annotation_to_role_rows, write_role_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Build ligand role tables for a directory of ligands.")
    parser.add_argument("--ligand-source-dir", required=True)
    parser.add_argument("--prepared-pdbqt-dir", required=True)
    parser.add_argument("--rules", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--source-glob", default="*")
    parser.add_argument("--pdbqt-pattern", default="{ligand_id}.pdbqt")
    parser.add_argument("--groups", nargs="*", default=None)
    parser.add_argument("--max-dist", type=float, default=0.5)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    ann_dir = out_dir / "annotations"
    map_dir = out_dir / "atom_maps"
    table_dir = out_dir / "role_tables"
    report_rows = []

    for source_path in iter_ligand_sources(args.ligand_source_dir, args.source_glob):
        ligand_id = source_path.stem
        prepared_pdbqt = Path(args.prepared_pdbqt_dir) / args.pdbqt_pattern.format(ligand_id=ligand_id)
        status = "ok"
        group_count = 0
        role_rows = 0
        matched = 0
        error = ""
        try:
            annotation = annotate_ligand(source_path, ligand_id, args.rules, groups=args.groups)
            atom_map = recover_atom_map_by_coordinates(source_path, prepared_pdbqt, ligand_id, max_dist=args.max_dist)
            rows = annotation_to_role_rows(annotation, atom_map=atom_map)
            write_annotation(annotation, ann_dir / f"{ligand_id}.annotation.json")
            write_atom_map(atom_map, map_dir / f"{ligand_id}.atom_map.json")
            write_role_table(rows, table_dir / f"{ligand_id}.role_table.csv")
            group_count = len(annotation.functional_groups)
            role_rows = len(rows)
            matched = int(atom_map["matched_atoms"])
        except Exception as exc:
            status = "error"
            error = str(exc)
        report_rows.append(
            {
                "ligand_id": ligand_id,
                "status": status,
                "groups": group_count,
                "role_rows": role_rows,
                "matched_heavy_atoms": matched,
                "error": error,
            }
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "batch_report.csv"
    with open(report_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["ligand_id", "status", "groups", "role_rows", "matched_heavy_atoms", "error"])
        writer.writeheader()
        writer.writerows(report_rows)

    ok_count = sum(1 for row in report_rows if row["status"] == "ok")
    print(f"ligands={len(report_rows)}")
    print(f"ok={ok_count}")
    print(f"report={report_path}")


if __name__ == "__main__":
    main()
