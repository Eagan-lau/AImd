from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from rdkit import Chem
from scipy.optimize import linear_sum_assignment

from .molio import read_molecule, require_3d_coordinates
from .pdbqt import PdbqtAtom, read_pdbqt_atoms


def _source_heavy_atoms(mol: Chem.Mol) -> List[Dict[str, Any]]:
    conf = mol.GetConformer()
    out: List[Dict[str, Any]] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        out.append({"source_index": idx, "element": atom.GetSymbol(), "coord": np.array([pos.x, pos.y, pos.z])})
    return out


def _pdbqt_heavy_atoms(pdbqt_atoms: List[PdbqtAtom]) -> List[Dict[str, Any]]:
    return [
        {"pdbqt_order": atom.order, "element": atom.element, "coord": atom.coord, "atom_name": atom.atom_name, "atom_type": atom.atom_type}
        for atom in pdbqt_atoms
        if atom.element != "H"
    ]


def _match_element_block(src: List[Dict[str, Any]], dst: List[Dict[str, Any]], max_dist: float) -> Tuple[Dict[str, Any], List[float]]:
    if len(src) > len(dst):
        element = src[0]["element"] if src else "unknown"
        raise ValueError(f"Not enough PDBQT heavy atoms for element {element}: source={len(src)}, pdbqt={len(dst)}")
    cost = np.zeros((len(src), len(dst)), dtype=float)
    for i, s in enumerate(src):
        for j, d in enumerate(dst):
            cost[i, j] = float(np.linalg.norm(s["coord"] - d["coord"]))
    row_ind, col_ind = linear_sum_assignment(cost)
    mapping: Dict[str, Any] = {}
    distances: List[float] = []
    for i, j in zip(row_ind, col_ind):
        dist = float(cost[i, j])
        if dist > max_dist:
            raise ValueError(
                "Coordinate match exceeds threshold: "
                f"element={src[i]['element']}, source_index={src[i]['source_index']}, "
                f"pdbqt_order={dst[j]['pdbqt_order']}, distance={dist:.6f}, max_dist={max_dist}"
            )
        mapping[str(src[i]["source_index"])] = {
            "pdbqt_order": int(dst[j]["pdbqt_order"]),
            "element": src[i]["element"],
            "distance": dist,
            "atom_name": dst[j]["atom_name"],
            "atom_type": dst[j]["atom_type"],
        }
        distances.append(dist)
    return mapping, distances


def recover_atom_map_by_coordinates(
    ligand_source: str | Path,
    prepared_pdbqt: str | Path,
    ligand_id: str,
    max_dist: float = 0.5,
) -> Dict[str, Any]:
    mol = read_molecule(ligand_source)
    require_3d_coordinates(mol, ligand_source)
    source_atoms = _source_heavy_atoms(mol)
    pdbqt_atoms = _pdbqt_heavy_atoms(read_pdbqt_atoms(prepared_pdbqt))

    mapping: Dict[str, Any] = {}
    distances: List[float] = []
    elements = sorted(set(a["element"] for a in source_atoms))
    for element in elements:
        src = [a for a in source_atoms if a["element"] == element]
        dst = [a for a in pdbqt_atoms if a["element"] == element]
        block_map, block_dist = _match_element_block(src, dst, max_dist=max_dist)
        mapping.update(block_map)
        distances.extend(block_dist)

    rmsd = math.sqrt(sum(d * d for d in distances) / max(len(distances), 1))
    return {
        "ligand_id": ligand_id,
        "source_file": str(ligand_source),
        "pdbqt_file": str(prepared_pdbqt),
        "index_base": 0,
        "pdbqt_order_base": 1,
        "map_scope": "heavy_atoms_only",
        "match_method": "element_coordinate_assignment",
        "matched_atoms": len(mapping),
        "rmsd": rmsd,
        "max_distance": max(distances) if distances else None,
        "map": mapping,
    }


def write_atom_map(atom_map: Dict[str, Any], out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(atom_map, handle, indent=2, sort_keys=True)


def load_atom_map(path: str | Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)
