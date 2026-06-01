#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys

from metaboclip_ligand_roles.rule_metadata import iter_atom_label_records
from metaboclip_ligand_roles.rules import load_rules


FIELDS = ["atom_label", "group_id", "atom_role", "atom_class", "element", "description"]


def main() -> None:
    parser = argparse.ArgumentParser(description="List stable atom labels available in a rule library.")
    parser.add_argument("--rules", required=True)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    records = list(iter_atom_label_records(load_rules(args.rules)))
    if args.out:
        with open(args.out, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            for record in records:
                writer.writerow({key: record.get(key, "") for key in FIELDS})
        print(f"atom_labels={len(records)}")
        print(f"out={args.out}")
    else:
        writer = csv.DictWriter(sys.stdout, fieldnames=FIELDS)
        writer.writeheader()
        for record in records:
            writer.writerow({key: record.get(key, "") for key in FIELDS})


if __name__ == "__main__":
    main()
