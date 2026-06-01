from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Any
import math
import re

import numpy as np
try:
    from scipy.spatial import cKDTree
except Exception:  # pragma: no cover
    cKDTree = None


@dataclass(frozen=True)
class Atom:
    serial: int
    atom_name: str
    resn: str
    chain: str
    resi: str
    element: str
    x: float
    y: float
    z: float
    source: str = "protein"
    role: str | None = None
    extra: dict | None = None

    @property
    def coord(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)

    @property
    def residue_key(self) -> tuple[str, str, str]:
        return (self.chain, self.resi, self.resn)

    @property
    def atom_key(self) -> tuple[str, str, str, str, int]:
        return (self.chain, self.resi, self.resn, self.atom_name, self.serial)


def distance(a: Atom, b: Atom) -> float:
    return float(np.linalg.norm(a.coord - b.coord))


_PDBQT_TYPE_TO_ELEMENT = {
    "A": "C",
    "C": "C",
    "N": "N",
    "NA": "N",
    "OA": "O",
    "OS": "O",
    "O": "O",
    "S": "S",
    "SA": "S",
    "H": "H",
    "HD": "H",
    "HS": "H",
    "F": "F",
    "CL": "Cl",
    "BR": "Br",
    "I": "I",
    "P": "P",
    "FE": "Fe",
    "ZN": "Zn",
    "MG": "Mg",
    "MN": "Mn",
    "CU": "Cu",
    "CO": "Co",
    "NI": "Ni",
    "CA": "Ca",
    "NA+": "Na",
    "K": "K",
}

_TWO_LETTER_ELEMENTS = {"Cl", "Br", "Fe", "Mg", "Zn", "Mn", "Cu", "Co", "Ni", "Ca", "Na"}


def infer_element(atom_name: str, atom_type: str | None = None) -> str:
    atom_type_clean = (atom_type or "").strip().replace(".", "").upper()
    if atom_type_clean in _PDBQT_TYPE_TO_ELEMENT:
        return _PDBQT_TYPE_TO_ELEMENT[atom_type_clean]
    name = (atom_name or "").strip()
    name = re.sub(r"^[0-9]+", "", name).replace(".", "")
    if not name:
        return ""
    two = name[:2].capitalize()
    if two in _TWO_LETTER_ELEMENTS:
        return two
    return name[0].upper()


def _parse_file_serial(line: str) -> int | None:
    try:
        value = line[6:11].strip()
        return int(value) if value else None
    except Exception:
        return None


def _parse_coordinates(line: str) -> tuple[float, float, float]:
    if len(line) >= 54:
        try:
            return float(line[30:38]), float(line[38:46]), float(line[46:54])
        except Exception:
            pass
    parts = line.split()
    numeric: list[float] = []
    for token in parts:
        try:
            numeric.append(float(token))
        except Exception:
            continue
    if len(numeric) < 3:
        raise ValueError(f"Cannot parse coordinates from line: {line.rstrip()}")
    # For normal PDB/PDBQT records, x/y/z are the first three non-serial floating values.
    # The last floating fields are often charge or occupancy, so prefer fixed columns first.
    if len(parts) >= 8:
        for start in range(5, min(len(parts) - 2, 9)):
            try:
                return float(parts[start]), float(parts[start + 1]), float(parts[start + 2])
            except Exception:
                continue
    return numeric[-3], numeric[-2], numeric[-1]


def _parse_atom_line(line: str, internal_order: int, source: str, serial_policy: str = "file") -> Atom:
    parts = line.split()
    file_serial = _parse_file_serial(line)
    atom_name = line[12:16].strip() if len(line) >= 16 else (parts[2] if len(parts) > 2 else "")
    resn = line[17:20].strip() if len(line) >= 20 else (parts[3] if len(parts) > 3 else "")
    chain = line[21:22].strip() if len(line) >= 22 else ""
    resi = line[22:26].strip() if len(line) >= 26 else (parts[4] if len(parts) > 4 else "")
    x, y, z = _parse_coordinates(line)
    atom_type = parts[-1] if parts else atom_name
    element = infer_element(atom_name, atom_type)
    if serial_policy == "order":
        serial = internal_order
    else:
        serial = file_serial if file_serial is not None else internal_order
    extra = {
        "file_serial": file_serial if file_serial is not None else "",
        "internal_order": internal_order,
    }
    if serial_policy == "order":
        extra["pdbqt_order"] = internal_order
    return Atom(serial, atom_name, resn, chain, resi, element, x, y, z, source=source, extra=extra)


def read_structure_atoms(path: str | Path, source: str = "protein") -> list[Atom]:
    atoms: list[Atom] = []
    order = 0
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                order += 1
                atoms.append(_parse_atom_line(line, order, source, serial_policy="file"))
    return atoms


def parse_pdbqt_poses(path: str | Path) -> list[dict[str, Any]]:
    poses: list[dict[str, Any]] = []
    current_lines: list[str] = []
    current_affinity: float | None = None

    def flush() -> None:
        nonlocal current_lines, current_affinity
        atom_lines = [ln for ln in current_lines if ln.startswith(("ATOM", "HETATM"))]
        if not atom_lines:
            current_lines = []
            current_affinity = None
            return
        atoms: list[Atom] = []
        for order, ln in enumerate(atom_lines, start=1):
            atoms.append(_parse_atom_line(ln, order, source="ligand", serial_policy="order"))
        poses.append({"pose_id": len(poses) + 1, "atoms": atoms, "affinity_kcal": current_affinity})
        current_lines = []
        current_affinity = None

    with open(path, "r", encoding="utf-8") as handle:
        saw_model = False
        for line in handle:
            if line.startswith("MODEL"):
                saw_model = True
                if current_lines:
                    flush()
                continue
            if line.startswith("ENDMDL"):
                flush()
                continue
            if line.startswith("REMARK VINA RESULT"):
                match = re.search(r"RESULT:\s*([-+]?\d+(?:\.\d+)?)", line)
                if match:
                    current_affinity = float(match.group(1))
            if line.startswith(("ATOM", "HETATM", "REMARK")):
                current_lines.append(line)
        if current_lines:
            flush()
    return poses


def parse_vina_out(path: str | Path) -> dict[int, float]:
    results: dict[int, float] = {}
    if not Path(path).exists():
        return results
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            m = re.match(r"^\s*(\d+)\s+([-+]?\d+(?:\.\d+)?)", line)
            if m:
                results[int(m.group(1))] = float(m.group(2))
                continue
            if "REMARK VINA RESULT" in line:
                value = re.search(r"RESULT:\s*([-+]?\d+(?:\.\d+)?)", line)
                if value:
                    results[len(results) + 1] = float(value.group(1))
    return results


def atom_matches_residue_atoms(atom: Atom, residue_atoms: dict[str, list[str]]) -> bool:
    atom_resn = atom.resn.upper()
    atom_name = atom.atom_name.upper()
    for resn, names in residue_atoms.items():
        if atom_resn == str(resn).upper() and atom_name in {str(n).upper() for n in names}:
            return True
    return False


def atom_matches_selector(atom: Atom, selector: str) -> bool:
    text = selector.strip()
    clauses = [c.strip() for c in re.split(r"\band\b", text) if c.strip()]
    for clause in clauses:
        toks = clause.split()
        if len(toks) < 2:
            continue
        key = toks[0].lower()
        vals = {v.upper() for v in toks[1].replace("+", " ").split()}
        if key == "elem" and atom.element.upper() not in vals:
            return False
        if key == "name" and atom.atom_name.upper() not in vals:
            return False
        if key == "resn" and atom.resn.upper() not in vals:
            return False
        if key == "resi" and atom.resi.upper() not in vals:
            return False
        if key == "chain" and atom.chain.upper() not in vals:
            return False
    return True


def select_atoms(atoms: Iterable[Atom], residue_atoms: dict | None = None, atom_selectors: list[str] | None = None) -> list[Atom]:
    selected: list[Atom] = []
    residue_atoms = residue_atoms or {}
    atom_selectors = atom_selectors or []
    for atom in atoms:
        ok = False
        if residue_atoms and atom_matches_residue_atoms(atom, residue_atoms):
            ok = True
        if atom_selectors and any(atom_matches_selector(atom, sel) for sel in atom_selectors):
            ok = True
        if ok:
            selected.append(atom)
    return selected


def nearest_atoms(atoms: Iterable[Atom], anchor: Atom, radius: float) -> list[tuple[float, Atom]]:
    atom_list = list(atoms)
    if not atom_list:
        return []
    if cKDTree is not None and len(atom_list) >= 200:
        coords = np.array([a.coord for a in atom_list], dtype=float)
        tree = cKDTree(coords)
        idxs = tree.query_ball_point(anchor.coord, radius)
        out = [(float(np.linalg.norm(coords[i] - anchor.coord)), atom_list[i]) for i in idxs]
    else:
        out = [(distance(atom, anchor), atom) for atom in atom_list if distance(atom, anchor) <= radius]
    out.sort(key=lambda x: x[0])
    return out


def angle_3pt(a: Atom, vertex: Atom, c: Atom) -> float:
    v1 = a.coord - vertex.coord
    v2 = c.coord - vertex.coord
    denom = float(np.linalg.norm(v1) * np.linalg.norm(v2))
    if denom <= 1e-12:
        raise ValueError("Cannot compute angle with zero-length vector")
    cosv = float(np.dot(v1, v2) / denom)
    cosv = max(-1.0, min(1.0, cosv))
    return float(math.degrees(math.acos(cosv)))
