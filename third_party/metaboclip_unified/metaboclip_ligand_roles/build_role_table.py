#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from metaboclip_ligand_roles.annotator_core import annotate_ligand, write_annotation
from metaboclip_ligand_roles.atom_map import recover_atom_map_by_coordinates, write_atom_map
from metaboclip_ligand_roles.role_table import annotation_to_role_rows, write_role_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a concise heavy-atom ligand role table with stable atom labels.")
    parser.add_argument("--ligand-source", required=True, help="Original ligand file before PDBQT conversion.")
    parser.add_argument("--prepared-pdbqt", required=True, help="Prepared ligand PDBQT before docking.")
    parser.add_argument("--ligand-id", required=True, help="Ligand identifier used in output tables.")
    parser.add_argument("--rules", required=True, help="Functional group rule YAML file.")
    parser.add_argument("--out-role-table", required=True, help="Output CSV with stable ligand atom labels.")
    parser.add_argument("--out-annotation", default=None, help="Optional full annotation JSON output.")
    parser.add_argument("--out-atom-map", default=None, help="Optional recovered atom map JSON output.")
    parser.add_argument("--groups", nargs="*", default=None, help="Optional subset of functional groups to detect.")
    parser.add_argument("--max-dist", type=float, default=0.5, help="Maximum source-to-PDBQT heavy-atom match distance.")
    args = parser.parse_args()

    annotation = annotate_ligand(args.ligand_source, args.ligand_id, args.rules, groups=args.groups)
    atom_map = recover_atom_map_by_coordinates(args.ligand_source, args.prepared_pdbqt, args.ligand_id, max_dist=args.max_dist)
    rows = annotation_to_role_rows(annotation, atom_map=atom_map)

    write_role_table(rows, args.out_role_table)
    if args.out_annotation:
        write_annotation(annotation, args.out_annotation)
    if args.out_atom_map:
        write_atom_map(atom_map, args.out_atom_map)

    print(f"ligand_id={args.ligand_id}")
    print(f"groups={len(annotation.functional_groups)}")
    print(f"role_rows={len(rows)}")
    print(f"matched_heavy_atoms={atom_map['matched_atoms']}")
    print(f"role_table={args.out_role_table}")


if __name__ == "__main__":
    main()
