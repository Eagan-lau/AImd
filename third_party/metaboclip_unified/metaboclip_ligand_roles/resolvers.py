from __future__ import annotations

from collections import defaultdict, deque
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set

from rdkit import Chem

from .chem_utils import classify_ch_site, total_h_count
from .models import FunctionalGroupMatch
from .rule_metadata import build_role_metadata


def _ring_neighbors_in_ring(mol: Chem.Mol, atom_idx: int, ring_set: Set[int], element: Optional[str] = None) -> List[int]:
    atom = mol.GetAtomWithIdx(atom_idx)
    out = []
    for nbr in atom.GetNeighbors():
        idx = nbr.GetIdx()
        if idx in ring_set and (element is None or nbr.GetSymbol() == element):
            out.append(idx)
    return out


def _bfs_path_within_ring(mol: Chem.Mol, start: int, end: int, ring_set: Set[int], forbidden: Set[int]) -> Optional[List[int]]:
    queue = deque([[start]])
    visited = {start}
    while queue:
        path = queue.popleft()
        node = path[-1]
        if node == end:
            return path
        for nbr in mol.GetAtomWithIdx(node).GetNeighbors():
            idx = nbr.GetIdx()
            if idx not in ring_set or idx in forbidden or idx in visited:
                continue
            visited.add(idx)
            queue.append(path + [idx])
    return None


def _exocyclic_hetero_neighbors(mol: Chem.Mol, atom_idx: int, ring_set: Set[int], allowed: Set[str]) -> List[int]:
    out = []
    for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
        idx = nbr.GetIdx()
        if idx in ring_set:
            continue
        if nbr.GetSymbol() in allowed and nbr.GetAtomicNum() > 1:
            out.append(idx)
    return out


def _best_glycosidic_atom(mol: Chem.Mol, c1: int, candidates: Sequence[int]) -> Optional[int]:
    if not candidates:
        return None
    preferred = []
    for idx in candidates:
        atom = mol.GetAtomWithIdx(idx)
        heavy_neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() > 1 and n.GetIdx() != c1]
        if heavy_neighbors:
            preferred.append(idx)
    return int(preferred[0] if preferred else candidates[0])


def _oxygen_neighbor_outside_ring(mol: Chem.Mol, atom_idx: int, ring_set: Set[int]) -> Optional[int]:
    for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
        idx = nbr.GetIdx()
        if idx in ring_set:
            continue
        if nbr.GetSymbol() == "O":
            return idx
    return None


def _carbon_neighbor_outside_ring(mol: Chem.Mol, atom_idx: int, ring_set: Set[int]) -> Optional[int]:
    for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
        idx = nbr.GetIdx()
        if idx in ring_set:
            continue
        if nbr.GetSymbol() == "C":
            return idx
    return None


def resolve_glycosyl_pyranose(mol: Chem.Mol, group_id: str, rule: Dict[str, Any]) -> List[FunctionalGroupMatch]:
    params = (rule.get("detector") or {}).get("params", {})
    allowed = set(params.get("allowed_glycosidic_atoms", ["O", "N", "S"]))
    require_exocyclic = bool(params.get("require_exocyclic_glycosidic_atom", True))
    out: List[FunctionalGroupMatch] = []
    count = 0
    seen_c1: Set[int] = set()

    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6:
            continue
        ring_set = set(ring)
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "O"]
        ring_carbons = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "C"]
        if len(ring_oxygens) != 1 or len(ring_carbons) != 5:
            continue
        ring_o = ring_oxygens[0]
        o_adjacent_carbons = _ring_neighbors_in_ring(mol, ring_o, ring_set, element="C")
        if len(o_adjacent_carbons) != 2:
            continue

        for candidate_c1 in o_adjacent_carbons:
            glyco_candidates = _exocyclic_hetero_neighbors(mol, candidate_c1, ring_set, allowed)
            if require_exocyclic and not glyco_candidates:
                continue
            glyco_atom = _best_glycosidic_atom(mol, candidate_c1, glyco_candidates)
            if glyco_atom is None:
                continue
            c5 = [c for c in o_adjacent_carbons if c != candidate_c1][0]
            path = _bfs_path_within_ring(mol, candidate_c1, c5, ring_set, forbidden={ring_o})
            if path is None or len(path) != 5:
                continue
            c1, c2, c3, c4, c5_final = path
            if c1 in seen_c1:
                continue
            seen_c1.add(c1)
            c6 = _carbon_neighbor_outside_ring(mol, c5_final, ring_set)
            roles: Dict[str, int] = {
                "c1": int(c1),
                "glycosidic_atom": int(glyco_atom),
                "ring_o": int(ring_o),
                "c2": int(c2),
                "c3": int(c3),
                "c4": int(c4),
                "c5": int(c5_final),
            }
            if mol.GetAtomWithIdx(glyco_atom).GetSymbol() == "O":
                roles["glycosidic_o"] = int(glyco_atom)
            if c6 is not None:
                roles["c6"] = int(c6)
            for c_role in ["c2", "c3", "c4", "c6"]:
                if c_role in roles:
                    o_idx = _oxygen_neighbor_outside_ring(mol, roles[c_role], ring_set)
                    if o_idx is not None:
                        roles["o" + c_role[1:]] = int(o_idx)
            count += 1
            atoms = sorted(set(roles.values()))
            evidence = {
                "ring_size": 6,
                "ring_oxygen_count": 1,
                "ring_carbon_count": 5,
                "uses_heavy_atom_geometry": True,
            }
            out.append(
                FunctionalGroupMatch(
                    group_id=group_id,
                    instance_id=f"{group_id}_{count}",
                    roles=roles,
                    atoms=atoms,
                    role_metadata=build_role_metadata(group_id, rule, roles),
                    evidence=evidence,
                    subtype="hexopyranosyl_like",
                    confidence=float(rule.get("confidence", 1.0)),
                    priority=int(rule.get("priority", 0)),
                )
            )
    return out


def resolve_ch_sites(mol: Chem.Mol, group_id: str, rule: Dict[str, Any]) -> List[FunctionalGroupMatch]:
    out: List[FunctionalGroupMatch] = []
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        h_count = total_h_count(atom)
        if h_count <= 0:
            continue
        count += 1
        idx = atom.GetIdx()
        subtype = classify_ch_site(mol, idx)
        roles = {"c": int(idx)}
        evidence = {
            "total_h_count": int(h_count),
            "hydrogen_used_for_recognition_only": True,
            "uses_heavy_atom_geometry": True,
        }
        out.append(
            FunctionalGroupMatch(
                group_id=group_id,
                instance_id=f"{group_id}_{count}",
                roles=roles,
                atoms=[int(idx)],
                role_metadata=build_role_metadata(group_id, rule, roles),
                evidence=evidence,
                subtype=subtype,
                confidence=float(rule.get("confidence", 1.0)),
                priority=int(rule.get("priority", 0)),
            )
        )
    return out


def run_graph_resolvers(mol: Chem.Mol, rules: Dict[str, Any]) -> List[FunctionalGroupMatch]:
    out: List[FunctionalGroupMatch] = []
    for group_id, rule in rules.get("functional_groups", {}).items():
        detector = rule.get("detector", {})
        if detector.get("type") != "graph_resolver":
            continue
        resolver = detector.get("resolver")
        if resolver == "glycosyl_pyranose_resolver":
            out.extend(resolve_glycosyl_pyranose(mol, group_id, rule))
        elif resolver == "ch_site_resolver":
            out.extend(resolve_ch_sites(mol, group_id, rule))
        else:
            raise ValueError(f"Unknown graph resolver: {resolver}")
    return out
