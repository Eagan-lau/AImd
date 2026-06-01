from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from rdkit import Chem

from .models import FunctionalGroupMatch, LigandAnnotation
from .molio import read_molecule


ROLE_TABLE_FIELDS = [
    "ligand_id",
    "group_id",
    "instance_id",
    "atom_label",
    "atom_class",
    "atom_role",
    "source_atom_index",
    "element",
    "pdbqt_order",
    "subtype",
    "confidence",
]


def _atom_element(mol: Chem.Mol, atom_index: int) -> str:
    return mol.GetAtomWithIdx(int(atom_index)).GetSymbol()


def _metadata_for_role(match: FunctionalGroupMatch, role_name: str) -> Dict[str, Any]:
    metadata = match.role_metadata.get(role_name, {}) if match.role_metadata else {}
    return {
        "atom_label": metadata.get("atom_label") or f"{match.group_id}.{role_name}",
        "atom_class": metadata.get("atom_class") or "",
    }


def annotation_to_role_rows(
    annotation: LigandAnnotation,
    rules_source_mol: Optional[Chem.Mol] = None,
    atom_map: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    mol = rules_source_mol or read_molecule(annotation.source_file)
    map_data = (atom_map or {}).get("map", {})
    rows: List[Dict[str, Any]] = []
    for match in annotation.functional_groups:
        for atom_role, atom_index in sorted(match.roles.items()):
            source_idx = int(atom_index)
            map_entry = map_data.get(str(source_idx), {})
            metadata = _metadata_for_role(match, atom_role)
            rows.append(
                {
                    "ligand_id": annotation.ligand_id,
                    "group_id": match.group_id,
                    "instance_id": match.instance_id,
                    "atom_label": metadata["atom_label"],
                    "atom_class": metadata["atom_class"],
                    "atom_role": atom_role,
                    "source_atom_index": source_idx,
                    "element": _atom_element(mol, source_idx),
                    "pdbqt_order": map_entry.get("pdbqt_order", ""),
                    "subtype": match.subtype or "",
                    "confidence": match.confidence,
                }
            )
    rows.sort(
        key=lambda row: (
            row["ligand_id"],
            row["group_id"],
            row["instance_id"],
            row["atom_label"],
            int(row["source_atom_index"]),
        )
    )
    return rows


def write_role_table(rows: List[Dict[str, Any]], out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROLE_TABLE_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in ROLE_TABLE_FIELDS})


def load_annotation(path: str | Path) -> LigandAnnotation:
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    groups = [
        FunctionalGroupMatch(
            group_id=item["group_id"],
            instance_id=item["instance_id"],
            roles={k: int(v) for k, v in item.get("roles", {}).items()},
            atoms=[int(v) for v in item.get("atoms", [])],
            role_metadata=item.get("role_metadata", {}),
            evidence=item.get("evidence", {}),
            subtype=item.get("subtype"),
            confidence=float(item.get("confidence", 1.0)),
            priority=int(item.get("priority", 0)),
        )
        for item in data.get("functional_groups", [])
    ]
    return LigandAnnotation(
        ligand_id=data["ligand_id"],
        source_file=data["source_file"],
        index_base=int(data.get("index_base", 0)),
        framework_policy=data.get("framework_policy", {}),
        functional_groups=groups,
    )


def read_role_table(path: str | Path) -> List[Dict[str, Any]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))
