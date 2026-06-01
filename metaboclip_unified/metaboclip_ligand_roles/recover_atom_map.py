#!/usr/bin/env python3
from __future__ import annotations

import argparse

from metaboclip_ligand_roles.atom_map import recover_atom_map_by_coordinates, write_atom_map


def main() -> None:
    parser = argparse.ArgumentParser(description="Recover source heavy atom to PDBQT atom order mapping by coordinates.")
    parser.add_argument("--ligand-source", required=True)
    parser.add_argument("--prepared-pdbqt", required=True)
    parser.add_argument("--ligand-id", required=True)
    parser.add_argument("--out-atom-map", required=True)
    parser.add_argument("--max-dist", type=float, default=0.5)
    args = parser.parse_args()

    atom_map = recover_atom_map_by_coordinates(args.ligand_source, args.prepared_pdbqt, args.ligand_id, max_dist=args.max_dist)
    write_atom_map(atom_map, args.out_atom_map)
    print(f"ligand_id={args.ligand_id}")
    print(f"matched_heavy_atoms={atom_map['matched_atoms']}")
    print(f"rmsd={atom_map['rmsd']:.6f}")
    print(f"atom_map={args.out_atom_map}")


if __name__ == "__main__":
    main()
