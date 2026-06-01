from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set

from .pdbqt import read_pdbqt_poses


POSE_ROLE_COORD_FIELDS = [
    "ligand_id",
    "pose_id",
    "affinity_kcal",
    "group_id",
    "instance_id",
    "atom_label",
    "atom_class",
    "atom_role",
    "source_atom_index",
    "element",
    "pdbqt_order",
    "x",
    "y",
    "z",
]


def extract_pose_role_coordinates(
    role_rows: List[Dict[str, Any]],
    docked_pdbqt: str | Path,
    atom_labels: Optional[Sequence[str]] = None,
    atom_classes: Optional[Sequence[str]] = None,
) -> List[Dict[str, Any]]:
    allowed_labels: Optional[Set[str]] = set(atom_labels) if atom_labels else None
    allowed_classes: Optional[Set[str]] = set(atom_classes) if atom_classes else None
    poses = read_pdbqt_poses(docked_pdbqt)
    out: List[Dict[str, Any]] = []
    for pose in poses:
        coords = pose.coord_by_order
        for row in role_rows:
            if allowed_labels is not None and row.get("atom_label") not in allowed_labels:
                continue
            if allowed_classes is not None and row.get("atom_class") not in allowed_classes:
                continue
            order_text = str(row.get("pdbqt_order", "")).strip()
            if not order_text:
                continue
            order = int(order_text)
            if order not in coords:
                continue
            xyz = coords[order]
            out.append(
                {
                    "ligand_id": row.get("ligand_id", ""),
                    "pose_id": pose.pose_id,
                    "affinity_kcal": "" if pose.affinity_kcal is None else pose.affinity_kcal,
                    "group_id": row.get("group_id", ""),
                    "instance_id": row.get("instance_id", ""),
                    "atom_label": row.get("atom_label", ""),
                    "atom_class": row.get("atom_class", ""),
                    "atom_role": row.get("atom_role", ""),
                    "source_atom_index": row.get("source_atom_index", ""),
                    "element": row.get("element", ""),
                    "pdbqt_order": order,
                    "x": float(xyz[0]),
                    "y": float(xyz[1]),
                    "z": float(xyz[2]),
                }
            )
    return out


def write_pose_role_coordinates(rows: List[Dict[str, Any]], out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=POSE_ROLE_COORD_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in POSE_ROLE_COORD_FIELDS})
