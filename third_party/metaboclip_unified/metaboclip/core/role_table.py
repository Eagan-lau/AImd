from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
from typing import Any

from metaboclip.core.atoms import Atom


@dataclass(frozen=True)
class LigandRoleAtom:
    ligand_id: str
    group_id: str
    instance_id: str
    atom_label: str
    atom_class: str
    atom_role: str
    source_atom_index: str
    element: str
    pdbqt_order: int
    subtype: str = ""
    confidence: float = 1.0

    def to_atom(self, pose_atoms_by_order: dict[int, Atom], site_name: str) -> Atom | None:
        atom = pose_atoms_by_order.get(self.pdbqt_order)
        if atom is None:
            return None
        extra = dict(atom.extra or {})
        extra.update({
            "ligand_id": self.ligand_id,
            "group_id": self.group_id,
            "instance_id": self.instance_id,
            "atom_label": self.atom_label,
            "atom_class": self.atom_class,
            "atom_role": self.atom_role,
            "source_atom_index": self.source_atom_index,
            "site_name": site_name,
        })
        return Atom(
            serial=atom.serial,
            atom_name=atom.atom_name,
            resn=atom.resn,
            chain=atom.chain,
            resi=atom.resi,
            element=atom.element,
            x=atom.x,
            y=atom.y,
            z=atom.z,
            source="ligand",
            role=site_name,
            extra=extra,
        )


def _as_int(value: Any, default: int = -1) -> int:
    try:
        return int(float(str(value)))
    except Exception:
        return default


def _as_float(value: Any, default: float = 1.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def read_role_table(path: str | Path) -> list[LigandRoleAtom]:
    rows: list[LigandRoleAtom] = []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            order = _as_int(row.get("pdbqt_order"))
            if order < 0:
                continue
            rows.append(
                LigandRoleAtom(
                    ligand_id=str(row.get("ligand_id", "")),
                    group_id=str(row.get("group_id", "")),
                    instance_id=str(row.get("instance_id", "")),
                    atom_label=str(row.get("atom_label", "")),
                    atom_class=str(row.get("atom_class", "")),
                    atom_role=str(row.get("atom_role", "")),
                    source_atom_index=str(row.get("source_atom_index", "")),
                    element=str(row.get("element", "")),
                    pdbqt_order=order,
                    subtype=str(row.get("subtype", "")),
                    confidence=_as_float(row.get("confidence"), 1.0),
                )
            )
    return rows


def select_role_rows(rows: list[LigandRoleAtom], labels: list[str] | None, classes: list[str] | None) -> list[LigandRoleAtom]:
    label_set = set(labels or [])
    class_set = set(classes or [])
    selected: list[LigandRoleAtom] = []
    for row in rows:
        if label_set and row.atom_label in label_set:
            selected.append(row)
            continue
        if class_set and row.atom_class in class_set:
            selected.append(row)
            continue
    return selected
