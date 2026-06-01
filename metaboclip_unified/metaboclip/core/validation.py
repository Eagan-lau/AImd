from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from metaboclip.core.config import load_yaml


@dataclass
class ValidationReport:
    ok: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {"ok": self.ok, "errors": self.errors, "warnings": self.warnings}


def _feature_items(mechanism: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    cfg = mechanism.get("features") or mechanism.get("geometry_features", {}).get("features") or {}
    if isinstance(cfg, list):
        return [(str(f.get("name")), f) for f in cfg if isinstance(f, dict) and f.get("name")]
    if isinstance(cfg, dict):
        return [(str(k), v) for k, v in cfg.items() if isinstance(v, dict)]
    return []


def _geometry_ref_items(mechanism: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    cfg = mechanism.get("geometry_refs") or mechanism.get("geometry_references", {}).get("references") or {}
    if isinstance(cfg, list):
        return [(str(r.get("name")), r) for r in cfg if isinstance(r, dict) and r.get("name")]
    if isinstance(cfg, dict):
        return [(str(k), v) for k, v in cfg.items() if isinstance(v, dict)]
    return []


def _split_ref(ref: Any) -> tuple[str, str] | None:
    if not isinstance(ref, str) or "." not in ref:
        return None
    prefix, name = ref.split(".", 1)
    return prefix, name


def _check_atom_ref(ref: Any, ligand_sites: set[str], protein_roles: set[str], geometry_refs: set[str], errors: list[str], context: str) -> None:
    parsed = _split_ref(ref)
    if parsed is None:
        errors.append(f"{context}: invalid atom reference {ref!r}")
        return
    prefix, name = parsed
    if prefix == "ligand" and name not in ligand_sites:
        errors.append(f"{context}: unknown ligand site {name!r}")
    elif prefix == "protein" and name not in protein_roles:
        errors.append(f"{context}: unknown protein role {name!r}")
    elif prefix == "geometry" and name not in geometry_refs:
        errors.append(f"{context}: unknown geometry reference {name!r}")
    elif prefix not in {"ligand", "protein", "geometry"}:
        errors.append(f"{context}: unsupported reference prefix {prefix!r}")


def validate_mechanism_dict(mechanism: dict[str, Any], mechanism_path: Path | None = None) -> ValidationReport:
    errors: list[str] = []
    warnings: list[str] = []

    ligand_cfg = mechanism.get("ligand_sites") or {}
    if not isinstance(ligand_cfg, dict) or not ligand_cfg:
        errors.append("ligand_sites must be a non-empty mapping")
    ligand_sites = set(ligand_cfg.keys()) if isinstance(ligand_cfg, dict) else set()
    for site, spec in ligand_cfg.items() if isinstance(ligand_cfg, dict) else []:
        if not isinstance(spec, dict):
            errors.append(f"ligand_sites.{site} must be a mapping")
            continue
        if not spec.get("atom_labels") and not spec.get("atom_classes"):
            warnings.append(f"ligand_sites.{site} has neither atom_labels nor atom_classes")
        linked_to = spec.get("linked_to")
        if linked_to and linked_to not in ligand_sites:
            errors.append(f"ligand_sites.{site}.linked_to references unknown ligand site {linked_to!r}")

    templates = mechanism.get("protein_templates") or {}
    if isinstance(templates, list):
        template_keys = {str(t.get("key")) for t in templates if isinstance(t, dict) and t.get("key")}
    elif isinstance(templates, dict):
        template_keys = set(templates.keys())
    else:
        template_keys = set()
        if templates:
            errors.append("protein_templates must be a mapping or a list")

    roles_cfg = mechanism.get("protein_roles") or mechanism.get("protein_site_map", {}).get("role_chain") or []
    if not isinstance(roles_cfg, list):
        errors.append("protein_roles must be a list")
        roles_cfg = []
    protein_roles = {str(r.get("role")) for r in roles_cfg if isinstance(r, dict) and r.get("role")}
    seen_roles: set[str] = set()
    for idx, spec in enumerate(roles_cfg):
        if not isinstance(spec, dict):
            errors.append(f"protein_roles[{idx}] must be a mapping")
            continue
        role = spec.get("role")
        if not role:
            errors.append(f"protein_roles[{idx}] is missing role")
            continue
        from_ref = str(spec.get("from", ""))
        parsed = _split_ref(from_ref)
        if parsed is None:
            errors.append(f"protein_roles.{role}.from is invalid")
        else:
            prefix, name = parsed
            if prefix == "ligand" and name not in ligand_sites:
                errors.append(f"protein_roles.{role}.from references unknown ligand site {name!r}")
            elif prefix == "protein" and name not in seen_roles:
                errors.append(f"protein_roles.{role}.from references protein role {name!r} before it is resolved")
            elif prefix == "template" and name not in template_keys:
                errors.append(f"protein_roles.{role}.from references unknown template {name!r}")
            elif prefix not in {"ligand", "protein", "template"}:
                errors.append(f"protein_roles.{role}.from has unsupported prefix {prefix!r}")
        if not spec.get("residues") and not spec.get("residue_atoms") and not spec.get("atoms") and not spec.get("atom_selectors") and not spec.get("template_atom"):
            warnings.append(f"protein_roles.{role} has no candidate atom rule")
        if spec.get("collection") and spec.get("min_count") is None:
            warnings.append(f"protein_roles.{role} is a collection role without min_count")
        seen_roles.add(str(role))

    geometry_refs = {name for name, _ in _geometry_ref_items(mechanism)}
    for name, spec in _geometry_ref_items(mechanism):
        typ = spec.get("type")
        if typ == "best_fit_plane":
            ref = spec.get("atoms_from_collection") or spec.get("atoms_from")
            if ref:
                _check_atom_ref(ref, ligand_sites, protein_roles, geometry_refs, errors, f"geometry_refs.{name}")
        elif typ == "oriented_axis":
            _check_atom_ref(spec.get("origin"), ligand_sites, protein_roles, geometry_refs, errors, f"geometry_refs.{name}.origin")
            direction = spec.get("direction") or {}
            if isinstance(direction, dict) and direction.get("type") == "plane_normal":
                plane_ref = direction.get("plane")
                _check_atom_ref(plane_ref, ligand_sites, protein_roles, geometry_refs, errors, f"geometry_refs.{name}.direction.plane")
            orient = spec.get("orient") or {}
            if orient.get("reference"):
                _check_atom_ref(orient.get("reference"), ligand_sites, protein_roles, geometry_refs, errors, f"geometry_refs.{name}.orient.reference")

    for name, spec in _feature_items(mechanism):
        if spec.get("enabled", True) is False:
            continue
        typ = spec.get("type")
        if typ == "distance":
            atoms = spec.get("atoms") or []
            if len(atoms) != 2:
                errors.append(f"features.{name}.atoms must contain exactly two atom references")
            for ref in atoms:
                _check_atom_ref(ref, ligand_sites, protein_roles, geometry_refs, errors, f"features.{name}")
        elif typ == "angle_3pt":
            for key in ("a", "vertex", "c"):
                _check_atom_ref(spec.get(key), ligand_sites, protein_roles, geometry_refs, errors, f"features.{name}.{key}")
        elif typ == "axis_deviation":
            _check_atom_ref(spec.get("axis"), ligand_sites, protein_roles, geometry_refs, errors, f"features.{name}.axis")
            vector = spec.get("vector") or {}
            _check_atom_ref(vector.get("from"), ligand_sites, protein_roles, geometry_refs, errors, f"features.{name}.vector.from")
            _check_atom_ref(vector.get("to"), ligand_sites, protein_roles, geometry_refs, errors, f"features.{name}.vector.to")
        else:
            errors.append(f"features.{name}.type is unsupported or missing")
        score = spec.get("score") or {}
        if score and score.get("transform") in {"distance_piecewise", "angle_gaussian", "angle_deviation_linear"}:
            if score.get("weight") is None:
                warnings.append(f"features.{name}.score has no weight")

    return ValidationReport(ok=not errors, errors=errors, warnings=warnings)


def validate_mechanism_file(path: str | Path) -> ValidationReport:
    p = Path(path)
    mechanism = load_yaml(p)
    return validate_mechanism_dict(mechanism, p)
