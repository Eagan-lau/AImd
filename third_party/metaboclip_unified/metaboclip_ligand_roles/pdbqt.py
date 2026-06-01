from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np


@dataclass
class PdbqtAtom:
    order: int
    atom_name: str
    atom_type: str
    element: str
    coord: np.ndarray


@dataclass
class PdbqtPose:
    pose_id: int
    affinity_kcal: Optional[float]
    atoms: List[PdbqtAtom]

    @property
    def coord_by_order(self) -> Dict[int, np.ndarray]:
        return {atom.order: atom.coord for atom in self.atoms}


def infer_pdbqt_element(atom_name: str, atom_type: str) -> str:
    token = (atom_type or atom_name or "").strip()
    if token in {"A", "C"}:
        return "C"
    if token == "OA":
        return "O"
    if token == "NA":
        return "N"
    if token == "SA":
        return "S"
    if token == "HD":
        return "H"
    two = token[:2].capitalize()
    if two in {"Cl", "Br", "Fe", "Mg", "Zn", "Mn", "Cu", "Co", "Ni", "Ca", "Na", "Li", "Al", "Si"}:
        return two
    if token:
        return token[0].upper()
    name = atom_name.strip()
    if len(name) >= 2 and name[:2].capitalize() in {"Cl", "Br", "Fe", "Mg", "Zn", "Mn", "Cu", "Co", "Ni", "Ca", "Na"}:
        return name[:2].capitalize()
    return name[0].upper() if name else ""


def _parse_atom_line(line: str, order: int) -> PdbqtAtom:
    atom_name = line[12:16].strip()
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    parts = line.split()
    atom_type = parts[-1] if parts else atom_name
    element = infer_pdbqt_element(atom_name, atom_type)
    return PdbqtAtom(order=order, atom_name=atom_name, atom_type=atom_type, element=element, coord=np.array([x, y, z]))


def parse_vina_affinity(line: str) -> Optional[float]:
    if "VINA RESULT:" not in line:
        return None
    try:
        after = line.split("VINA RESULT:", 1)[1].strip()
        return float(after.split()[0])
    except Exception:
        return None


def read_pdbqt_atoms(path: str | Path) -> List[PdbqtAtom]:
    atoms: List[PdbqtAtom] = []
    order = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                order += 1
                atoms.append(_parse_atom_line(line, order))
    return atoms


def read_pdbqt_poses(path: str | Path) -> List[PdbqtPose]:
    poses: List[PdbqtPose] = []
    current_atoms: List[PdbqtAtom] = []
    affinity: Optional[float] = None
    in_model = False
    pose_id = 0
    order = 0

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                in_model = True
                current_atoms = []
                affinity = None
                order = 0
                continue
            parsed_affinity = parse_vina_affinity(line)
            if parsed_affinity is not None:
                affinity = parsed_affinity
            if line.startswith(("ATOM", "HETATM")):
                if not in_model and pose_id == 0 and not current_atoms:
                    in_model = True
                order += 1
                current_atoms.append(_parse_atom_line(line, order))
            if line.startswith("ENDMDL"):
                pose_id += 1
                poses.append(PdbqtPose(pose_id=pose_id, affinity_kcal=affinity, atoms=current_atoms))
                in_model = False
                current_atoms = []
                affinity = None
                order = 0

    if current_atoms:
        pose_id += 1
        poses.append(PdbqtPose(pose_id=pose_id, affinity_kcal=affinity, atoms=current_atoms))
    return poses
