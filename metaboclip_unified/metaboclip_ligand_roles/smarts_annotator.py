from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from rdkit import Chem

from .chem_utils import atom_symbol, classify_hydroxyl, first_heavy_neighbor, total_h_count
from .models import FunctionalGroupMatch
from .rule_metadata import build_role_metadata


def _resolve_role_rule(mol: Chem.Mol, roles: Dict[str, int], role_name: str, spec: Dict[str, Any]) -> Optional[int]:
    rule = spec.get("rule")
    if rule == "heavy_neighbor_of":
        source_role = spec["role"]
        exclude_roles = spec.get("exclude_roles", [])
        exclude_indices = [roles[r] for r in exclude_roles if r in roles]
        return first_heavy_neighbor(mol, roles[source_role], exclude=exclude_indices)
    return None


def _resolve_roles(mol: Chem.Mol, match: Tuple[int, ...], role_specs: Dict[str, Any]) -> Optional[Dict[str, int]]:
    roles: Dict[str, int] = {}
    pending: Dict[str, Any] = {}
    for role_name, spec in role_specs.items():
        if "match_atom" in spec:
            atom_index = match[int(spec["match_atom"])]
            if mol.GetAtomWithIdx(atom_index).GetAtomicNum() == 1:
                return None
            roles[role_name] = int(atom_index)
        else:
            pending[role_name] = spec

    for role_name, spec in pending.items():
        atom_index = _resolve_role_rule(mol, roles, role_name, spec)
        if atom_index is None and spec.get("required", True):
            return None
        if atom_index is not None:
            roles[role_name] = int(atom_index)
    return roles


def _build_evidence(mol: Chem.Mol, roles: Dict[str, int], rule: Dict[str, Any]) -> Dict[str, Any]:
    evidence: Dict[str, Any] = {}
    for name, spec in (rule.get("evidence") or {}).items():
        role = spec.get("role")
        if role not in roles:
            continue
        atom = mol.GetAtomWithIdx(roles[role])
        if "min_total_h" in spec:
            evidence[name] = total_h_count(atom) >= int(spec["min_total_h"])
        elif spec.get("type") == "total_h_count":
            evidence[name] = total_h_count(atom)
    return evidence


def _evidence_passes(evidence: Dict[str, Any], rule: Dict[str, Any]) -> bool:
    for name, spec in (rule.get("evidence") or {}).items():
        if spec.get("required", False) and evidence.get(name) is not True:
            return False
    return True


def _subtype_for_rule(mol: Chem.Mol, roles: Dict[str, int], group_id: str) -> Optional[str]:
    if group_id == "hydroxyl" and "o" in roles:
        return classify_hydroxyl(mol, roles["o"], roles.get("parent_atom"))
    return None


def annotate_smarts_groups(mol: Chem.Mol, rules: Dict[str, Any]) -> List[FunctionalGroupMatch]:
    matches_out: List[FunctionalGroupMatch] = []
    counts: Dict[str, int] = defaultdict(int)
    seen: Set[Tuple[str, Tuple[Tuple[str, int], ...]]] = set()

    for group_id, rule in rules.get("functional_groups", {}).items():
        detector = rule.get("detector", {})
        if detector.get("type") != "smarts":
            continue
        pattern = detector.get("pattern")
        if not pattern:
            continue
        query = Chem.MolFromSmarts(pattern)
        if query is None:
            raise ValueError(f"Invalid SMARTS for group {group_id}: {pattern}")
        for match in mol.GetSubstructMatches(query, uniquify=True):
            roles = _resolve_roles(mol, match, rule.get("roles", {}))
            if roles is None:
                continue
            evidence = _build_evidence(mol, roles, rule)
            if not _evidence_passes(evidence, rule):
                continue
            key = (group_id, tuple(sorted(roles.items())))
            if key in seen:
                continue
            seen.add(key)
            counts[group_id] += 1
            subtype = _subtype_for_rule(mol, roles, group_id)
            atoms = sorted(set(roles.values()))
            matches_out.append(
                FunctionalGroupMatch(
                    group_id=group_id,
                    instance_id=f"{group_id}_{counts[group_id]}",
                    roles=roles,
                    atoms=atoms,
                    role_metadata=build_role_metadata(group_id, rule, roles),
                    evidence=evidence,
                    subtype=subtype,
                    confidence=float(rule.get("confidence", 1.0)),
                    priority=int(rule.get("priority", 0)),
                )
            )
    return matches_out
