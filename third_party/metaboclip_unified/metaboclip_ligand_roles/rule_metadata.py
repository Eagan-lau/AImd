from __future__ import annotations

from typing import Any, Dict, Iterable


def build_role_metadata(group_id: str, rule: Dict[str, Any], roles: Dict[str, int]) -> Dict[str, Dict[str, Any]]:
    metadata: Dict[str, Dict[str, Any]] = {}
    role_specs = rule.get("roles", {}) or {}
    for role_name in roles:
        spec = role_specs.get(role_name, {}) or {}
        atom_label = spec.get("atom_label") or f"{group_id}.{role_name}"
        atom_class = spec.get("atom_class") or ""
        metadata[role_name] = {
            "atom_label": str(atom_label),
            "atom_class": str(atom_class),
            "description": str(spec.get("description", "")),
            "canonical_role": str(spec.get("canonical_role", role_name)),
        }
    return metadata


def iter_atom_label_records(rule_data: Dict[str, Any]) -> Iterable[Dict[str, Any]]:
    registry = rule_data.get("atom_labels", {}) or {}
    for atom_label, spec in sorted(registry.items()):
        yield {
            "atom_label": atom_label,
            "group_id": spec.get("group_id", ""),
            "atom_role": spec.get("atom_role", ""),
            "atom_class": spec.get("atom_class", ""),
            "element": spec.get("element", ""),
            "description": spec.get("description", ""),
        }
