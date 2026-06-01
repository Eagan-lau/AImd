#!/usr/bin/env python3
from __future__ import annotations

import argparse

from metaboclip_ligand_roles.role_coords import extract_pose_role_coordinates, write_pose_role_coordinates
from metaboclip_ligand_roles.role_table import read_role_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract pose coordinates for stable ligand atom labels.")
    parser.add_argument("--role-table", required=True)
    parser.add_argument("--docked-pdbqt", required=True)
    parser.add_argument("--out-pose-role-coords", required=True)
    parser.add_argument("--atom-labels", nargs="*", default=None)
    parser.add_argument("--atom-classes", nargs="*", default=None)
    args = parser.parse_args()

    rows = read_role_table(args.role_table)
    coord_rows = extract_pose_role_coordinates(
        rows,
        args.docked_pdbqt,
        atom_labels=args.atom_labels,
        atom_classes=args.atom_classes,
    )
    write_pose_role_coordinates(coord_rows, args.out_pose_role_coords)
    print(f"pose_role_rows={len(coord_rows)}")
    print(f"pose_role_coords={args.out_pose_role_coords}")


if __name__ == "__main__":
    main()
