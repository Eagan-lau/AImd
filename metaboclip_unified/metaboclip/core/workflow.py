from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import csv
import itertools
import re

from metaboclip.core.atoms import Atom, angle_3pt, distance, nearest_atoms, parse_pdbqt_poses, parse_vina_out, read_structure_atoms, select_atoms
from metaboclip.core.config import load_yaml
from metaboclip.core.geometry import best_fit_plane, oriented_axis, axis_deviation
from metaboclip.core.template_align import TemplateAlignmentError, load_template_atoms_with_alignment
from metaboclip.core.role_table import LigandRoleAtom, read_role_table, select_role_rows
from metaboclip.core.scoring import energy_linear, transform_score, weighted_mean, hierarchical_geometry_score


@dataclass
class CandidateState:
    ligand_sites: dict[str, Atom]
    protein_roles: dict[str, Any] = field(default_factory=dict)
    role_distances: dict[str, float] = field(default_factory=dict)
    features: dict[str, float] = field(default_factory=dict)
    feature_scores: dict[str, float] = field(default_factory=dict)
    level_scores: dict[int, float] = field(default_factory=dict)
    candidate_score: float = 0.0
    affinity_score: float = 0.0
    geometry_score: float = 0.0
    site_set_id: str = ""


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))


def _format_pattern(pattern: str, ligand_id: str, protein_name: str, ligand_extension: str = ".pdbqt") -> str:
    return pattern.format(ligand_id=ligand_id, protein_name=protein_name, ligand_extension=ligand_extension)


def split_protein_conformation_name(protein_name: str) -> tuple[str, str, int]:
    """Return protein_id, conformation_id, and conformation_index.

    Names ending with an underscore followed by digits are treated as
    alternative conformations. The unsuffixed structure is conformation 0.
    Examples:
      P001 -> (P001, 0, 0)
      P001_1 -> (P001, 1, 1)
      Tca01g00352.t1_pLDDT_91.1_pTM_0.91_5 ->
        (Tca01g00352.t1_pLDDT_91.1_pTM_0.91, 5, 5)
    """
    match = re.match(r"^(?P<base>.+)_(?P<idx>\d+)$", protein_name)
    if match:
        idx = int(match.group("idx"))
        return match.group("base"), str(idx), idx
    return protein_name, "0", 0


def _resolve_atom_ref(ref: str, state: CandidateState) -> Atom | list[Atom] | None:
    prefix, _, name = ref.partition(".")
    if prefix == "ligand":
        return state.ligand_sites.get(name)
    if prefix == "protein":
        return state.protein_roles.get(name)
    raise ValueError(f"Unsupported atom reference: {ref}")


def _first_atom(value: Any) -> Atom | None:
    if isinstance(value, Atom):
        return value
    if isinstance(value, list) and value:
        first = value[0]
        if isinstance(first, Atom):
            return first
    return None


def _role_values_as_atoms(value: Any) -> list[Atom]:
    if isinstance(value, Atom):
        return [value]
    if isinstance(value, list):
        return [v for v in value if isinstance(v, Atom)]
    return []


def _site_specs(mechanism: dict[str, Any]) -> dict[str, Any]:
    return mechanism.get("ligand_sites", {}) or {}


def _select_rows_for_site(spec: dict[str, Any], role_rows: list[LigandRoleAtom]) -> list[LigandRoleAtom]:
    labels = spec.get("atom_labels") or []
    classes = spec.get("atom_classes") or []
    return select_role_rows(role_rows, labels, classes)


def _same_group_instance(a: LigandRoleAtom, b: LigandRoleAtom) -> bool:
    return a.group_id == b.group_id and a.instance_id == b.instance_id


def _read_ligand_sites(mechanism: dict[str, Any], role_rows: list[LigandRoleAtom], pose_atoms: list[Atom]) -> list[dict[str, Atom]]:
    site_specs = _site_specs(mechanism)
    if not site_specs:
        return []
    pose_by_order = {a.serial: a for a in pose_atoms}
    resolved_rows: dict[str, list[LigandRoleAtom]] = {}
    for name, spec in site_specs.items():
        selected = _select_rows_for_site(spec, role_rows)
        if spec.get("required", True) and not selected:
            return []
        resolved_rows[name] = selected

    independent = [name for name, spec in site_specs.items() if not spec.get("linked_to")]
    if not independent:
        return []
    for name in independent:
        if site_specs[name].get("required", True) and not resolved_rows.get(name):
            return []

    independent_products = itertools.product(*[resolved_rows.get(name, []) for name in independent])
    states: list[dict[str, Atom]] = []
    for product in independent_products:
        row_by_site: dict[str, LigandRoleAtom] = dict(zip(independent, product))
        atom_by_site: dict[str, Atom] = {}
        ok = True
        for name, row in row_by_site.items():
            atom = row.to_atom(pose_by_order, name)
            if atom is None:
                ok = False
                break
            atom_by_site[name] = atom
        if not ok:
            continue

        pending = [name for name in site_specs if name not in row_by_site]
        partial: list[tuple[dict[str, LigandRoleAtom], dict[str, Atom]]] = [(dict(row_by_site), dict(atom_by_site))]
        for name in pending:
            spec = site_specs[name]
            linked_to = spec.get("linked_to")
            next_partial: list[tuple[dict[str, LigandRoleAtom], dict[str, Atom]]] = []
            for rows_map, atoms_map in partial:
                candidates = list(resolved_rows.get(name, []))
                if linked_to:
                    anchor_row = rows_map.get(str(linked_to))
                    if anchor_row is None:
                        if spec.get("required", False):
                            continue
                        next_partial.append((rows_map, atoms_map))
                        continue
                    candidates = [r for r in candidates if _same_group_instance(r, anchor_row)]
                if not candidates:
                    if spec.get("required", False):
                        continue
                    next_partial.append((rows_map, atoms_map))
                    continue
                for candidate_row in candidates:
                    atom = candidate_row.to_atom(pose_by_order, name)
                    if atom is None:
                        continue
                    nr = dict(rows_map)
                    na = dict(atoms_map)
                    nr[name] = candidate_row
                    na[name] = atom
                    next_partial.append((nr, na))
            partial = next_partial
        states.extend(atoms for _, atoms in partial)
    return states


def _candidate_atoms_for_role(protein_atoms: list[Atom], spec: dict[str, Any]) -> list[Atom]:
    residues = spec.get("residues") or spec.get("residue_atoms") or {}
    atoms = spec.get("atoms") or spec.get("atom_selectors") or []
    return select_atoms(protein_atoms, residues, atoms)


def _resolve_anchor(from_ref: str, state: CandidateState) -> list[Atom]:
    if from_ref.startswith("ligand.") or from_ref.startswith("protein."):
        value = _resolve_atom_ref(from_ref, state)
        return _role_values_as_atoms(value)
    return []


def _map_template_role(spec: dict[str, Any], protein_atoms: list[Atom], template_atoms: list[Atom]) -> tuple[Atom, float] | None:
    tmpl = spec.get("template_atom") or {}
    resn = str(tmpl.get("residue_name", "")).upper()
    resi = str(tmpl.get("residue_id", ""))
    atom_name = str(tmpl.get("atom_name", "")).upper()
    anchors = [a for a in template_atoms if (not resn or a.resn.upper() == resn) and (not resi or a.resi == resi) and (not atom_name or a.atom_name.upper() == atom_name)]
    if not anchors:
        return None
    anchor = anchors[0]
    search = spec.get("target_search") or {}
    radius = float(search.get("max_mapping_distance", search.get("search_radius", search.get("radius", 2.5))))
    candidates = select_atoms(
        protein_atoms,
        search.get("residues") or search.get("residue_atoms") or spec.get("residues") or {},
        search.get("atoms") or search.get("atom_selectors") or [],
    )
    if not candidates:
        return None
    hits = nearest_atoms(candidates, anchor, radius)
    if not hits:
        return None
    return hits[0][1], hits[0][0]


def _group_collection_hits(spec: dict[str, Any], hits: list[tuple[float, Atom]]) -> list[tuple[float, list[Atom]]]:
    min_count = int(spec.get("min_count", spec.get("min_atoms", 1)))
    same_residue = bool(spec.get("same_residue", False))
    preferred_count = spec.get("preferred_count")
    if same_residue:
        groups: dict[tuple[str, str, str], list[tuple[float, Atom]]] = {}
        for d, atom in hits:
            groups.setdefault(atom.residue_key, []).append((d, atom))
        out: list[tuple[float, list[Atom]]] = []
        for vals in groups.values():
            vals = sorted(vals, key=lambda x: x[0])
            if len(vals) < min_count:
                continue
            if preferred_count is not None:
                vals = vals[: int(preferred_count)]
            score = sum(d for d, _ in vals) / len(vals)
            out.append((score, [atom for _, atom in vals]))
        return sorted(out, key=lambda x: x[0])
    if len(hits) < min_count:
        return []
    if preferred_count is not None:
        hits = hits[: int(preferred_count)]
    return [(sum(d for d, _ in hits) / len(hits), [atom for _, atom in hits])]


def _expand_protein_roles(mechanism: dict[str, Any], protein_atoms: list[Atom], states: list[CandidateState], template_atoms_by_key: dict[str, list[Atom]]) -> list[CandidateState]:
    role_chain = mechanism.get("protein_roles") or mechanism.get("protein_site_map", {}).get("role_chain") or []
    expanded = states
    max_site_sets = int(mechanism.get("runtime", {}).get("max_site_sets_per_pose", 1000000))
    for spec in role_chain:
        role = spec["role"]
        required = bool(spec.get("required", True))
        radius = float(spec.get("radius", 0.0)) if spec.get("radius") is not None else None
        keep_top = int((spec.get("ranking") or {}).get("keep_top", spec.get("keep_top", 10)))
        collection = bool(spec.get("collection", spec.get("cardinality") == "many"))
        new_states: list[CandidateState] = []
        from_ref = str(spec.get("from", ""))

        if from_ref.startswith("template."):
            template_key = from_ref.split(".", 1)[1]
            template_atoms = template_atoms_by_key.get(template_key, [])
            for state in expanded:
                mapped = _map_template_role(spec, protein_atoms, template_atoms)
                if mapped is None:
                    if not required:
                        new_states.append(state)
                    continue
                atom, map_dist = mapped
                validate = spec.get("validate") or {}
                if validate:
                    anchors = _resolve_anchor(str(validate.get("from", "")), state)
                    if anchors:
                        maxd = float(validate.get("radius", 9999.0))
                        if min(distance(atom, a) for a in anchors) > maxd:
                            continue
                ns = CandidateState(dict(state.ligand_sites), dict(state.protein_roles), dict(state.role_distances))
                ns.protein_roles[role] = atom
                ns.role_distances[role] = map_dist
                new_states.append(ns)
            expanded = new_states[:max_site_sets]
            continue

        for state in expanded:
            anchors = _resolve_anchor(from_ref, state)
            if not anchors:
                if not required:
                    new_states.append(state)
                continue
            candidates = _candidate_atoms_for_role(protein_atoms, spec.get("candidates", spec))
            hits: list[tuple[float, Atom]] = []
            for anchor in anchors:
                if radius is None:
                    hits.extend((distance(atom, anchor), atom) for atom in candidates)
                else:
                    hits.extend(nearest_atoms(candidates, anchor, radius))
            dedup: dict[tuple, tuple[float, Atom]] = {}
            for d, atom in hits:
                key = atom.atom_key
                if key not in dedup or d < dedup[key][0]:
                    dedup[key] = (d, atom)
            hits = sorted(dedup.values(), key=lambda x: x[0])[:keep_top]
            if collection:
                groups = _group_collection_hits(spec, hits)
                if not groups:
                    if not required:
                        new_states.append(state)
                    continue
                for group_dist, group_atoms in groups[:keep_top]:
                    ns = CandidateState(dict(state.ligand_sites), dict(state.protein_roles), dict(state.role_distances))
                    ns.protein_roles[role] = group_atoms
                    ns.role_distances[role] = group_dist
                    new_states.append(ns)
            else:
                if not hits:
                    if not required:
                        new_states.append(state)
                    continue
                for d, atom in hits:
                    ns = CandidateState(dict(state.ligand_sites), dict(state.protein_roles), dict(state.role_distances))
                    ns.protein_roles[role] = atom
                    ns.role_distances[role] = d
                    new_states.append(ns)
        expanded = new_states[:max_site_sets]
    return expanded


def _load_template_atoms(mechanism: dict[str, Any], mechanism_path: Path | None, protein_file: Path, protein_atoms: list[Atom]) -> dict[str, list[Atom]]:
    templates = mechanism.get("protein_templates") or mechanism.get("protein_site_map", {}).get("templates") or {}
    if isinstance(templates, list):
        items = [(t.get("key", f"template_{i}"), t) for i, t in enumerate(templates) if isinstance(t, dict)]
    else:
        items = list(templates.items())
    out: dict[str, list[Atom]] = {}
    base_dir = mechanism_path.parent if mechanism_path else Path.cwd()
    for key, spec in items:
        path = Path(str(spec.get("path", "")))
        if not path.is_absolute():
            path = base_dir / path
        if not path.exists():
            if bool(spec.get("required", False)):
                raise FileNotFoundError(f"Required template not found: {path}")
            continue
        try:
            atoms, _quality = load_template_atoms_with_alignment(path, protein_file, protein_atoms, spec)
        except TemplateAlignmentError:
            if bool(spec.get("required", False)):
                raise
            continue
        out[str(key)] = atoms
    return out


def _build_geometry_refs(mechanism: dict[str, Any], state: CandidateState) -> dict[str, Any]:
    refs: dict[str, Any] = {}
    raw_refs = mechanism.get("geometry_refs") or mechanism.get("geometry_references", {}).get("references") or {}
    if isinstance(raw_refs, list):
        iterable = [(r.get("name"), r) for r in raw_refs if isinstance(r, dict) and r.get("name")]
    else:
        iterable = raw_refs.items()
    for name, spec in iterable:
        typ = spec.get("type")
        if typ == "best_fit_plane":
            atoms_ref = spec.get("atoms_from_collection") or spec.get("atoms_from")
            atoms = _role_values_as_atoms(_resolve_atom_ref(atoms_ref, state)) if atoms_ref else []
            if len(atoms) >= int(spec.get("min_atoms", 3)):
                refs[name] = best_fit_plane(atoms)
        elif typ == "oriented_axis":
            origin = _first_atom(_resolve_atom_ref(spec["origin"], state))
            if origin is None:
                continue
            direction_spec = spec.get("direction") or {}
            if isinstance(direction_spec, dict) and direction_spec.get("type") == "plane_normal":
                plane_name = str(direction_spec.get("plane", "")).replace("geometry.", "")
                plane = refs.get(plane_name)
                if plane is None:
                    continue
                direction = plane.normal
            else:
                continue
            orient = spec.get("orient") or {}
            away_atom = None
            if orient.get("method") == "away_from" and orient.get("reference"):
                away_atom = _first_atom(_resolve_atom_ref(orient["reference"], state))
            refs[name] = oriented_axis(origin, direction, away_atom)
    return refs


def _features_iter(mechanism: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    features_cfg = mechanism.get("features") or mechanism.get("geometry_features", {}).get("features") or {}
    if isinstance(features_cfg, list):
        return [(str(f.get("name")), f) for f in features_cfg if isinstance(f, dict) and f.get("name")]
    return [(str(k), v) for k, v in features_cfg.items()]


def _compute_features(mechanism: dict[str, Any], state: CandidateState) -> dict[str, float | None]:
    refs = _build_geometry_refs(mechanism, state)
    values: dict[str, float | None] = {}
    for name, spec in _features_iter(mechanism):
        if spec.get("enabled", True) is False:
            continue
        typ = spec.get("type")
        try:
            if typ == "distance":
                a_ref, b_ref = spec["atoms"]
                a = _first_atom(_resolve_atom_ref(a_ref, state))
                b = _first_atom(_resolve_atom_ref(b_ref, state))
                values[name] = distance(a, b) if a is not None and b is not None else None
            elif typ == "angle_3pt":
                a = _first_atom(_resolve_atom_ref(spec["a"], state))
                v = _first_atom(_resolve_atom_ref(spec["vertex"], state))
                c = _first_atom(_resolve_atom_ref(spec["c"], state))
                values[name] = angle_3pt(a, v, c) if a is not None and v is not None and c is not None else None
            elif typ == "axis_deviation":
                axis_name = str(spec["axis"]).replace("geometry.", "")
                axis = refs.get(axis_name)
                start = _first_atom(_resolve_atom_ref(spec["vector"]["from"], state))
                end = _first_atom(_resolve_atom_ref(spec["vector"]["to"], state))
                values[name] = axis_deviation(axis, start, end) if axis is not None and start is not None and end is not None else None
            else:
                values[name] = None
        except Exception:
            values[name] = None
    state.features = {k: v for k, v in values.items() if v is not None}
    return values


def _passes_gate(mechanism: dict[str, Any], features: dict[str, float | None]) -> bool:
    for name, spec in _features_iter(mechanism):
        if spec.get("enabled", True) is False:
            continue
        gate = spec.get("gate") or {}
        if not gate:
            continue
        required = bool(gate.get("required", spec.get("required", False)))
        value = features.get(name)
        if value is None:
            if required:
                return False
            continue
        if gate.get("min") is not None and value < float(gate["min"]):
            if required:
                return False
        if gate.get("max") is not None and value > float(gate["max"]):
            if required:
                return False
    return True


def _score_candidate(mechanism: dict[str, Any], profile: dict[str, Any], state: CandidateState, affinity_kcal: float | None) -> None:
    scoring_profile = profile.get("scoring") or profile
    aff_cfg = scoring_profile.get("affinity") or profile.get("affinity", {}).get("scoring") or {}
    affinity_score = energy_linear(affinity_kcal, float(aff_cfg.get("full", -7.0)), float(aff_cfg.get("zero", -3.0)))
    state.affinity_score = affinity_score

    level_terms: dict[int, list[tuple[float, float]]] = {}
    feature_scores: dict[str, float] = {}
    for name, spec in _features_iter(mechanism):
        if spec.get("enabled", True) is False:
            continue
        score_cfg = spec.get("score") or {}
        if not score_cfg or score_cfg.get("enabled", True) is False:
            continue
        if name not in state.features:
            continue
        cfg = dict(score_cfg)
        cfg["source"] = name
        val = transform_score(cfg, state.features)
        feature_scores[f"{name}_score"] = val
        level = int(score_cfg.get("level", spec.get("level", 1)))
        weight = float(score_cfg.get("weight", 1.0))
        level_terms.setdefault(level, []).append((val, weight))
    state.feature_scores = feature_scores
    state.level_scores = {level: weighted_mean(vals) or 0.0 for level, vals in level_terms.items()}
    geom_profile = scoring_profile.get("geometry") or profile.get("geometry_scoring", {}) or {}
    level_weights = {int(k): float(v) for k, v in (geom_profile.get("level_weights") or {1: 0.60, 2: 0.25, 3: 0.15}).items()}
    state.geometry_score = hierarchical_geometry_score(state.level_scores, level_weights)
    aff_weight = float((scoring_profile.get("affinity") or {}).get("weight", 0.30))
    geom_weight = float((scoring_profile.get("geometry") or {}).get("weight", 0.70))
    total_weight = aff_weight + geom_weight
    if total_weight <= 0:
        total_weight = 1.0
    state.candidate_score = 100.0 * ((aff_weight * affinity_score + geom_weight * state.geometry_score) / total_weight)


def _state_row(state: CandidateState, context: dict[str, Any]) -> dict[str, Any]:
    row = dict(context)
    row.update({
        "site_set_id": state.site_set_id,
        "candidate_score": round(state.candidate_score, 6),
        "affinity_score": round(state.affinity_score, 6),
        "geometry_score": round(state.geometry_score, 6),
    })
    for site, atom in state.ligand_sites.items():
        extra = atom.extra or {}
        row[f"ligand_{site}_label"] = extra.get("atom_label", "")
        row[f"ligand_{site}_class"] = extra.get("atom_class", "")
        row[f"ligand_{site}_instance"] = extra.get("instance_id", "")
        row[f"ligand_{site}_pdbqt_order"] = atom.serial
    for role, val in state.protein_roles.items():
        atom = _first_atom(val)
        if atom is not None:
            row[f"protein_{role}_residue"] = f"{atom.chain}:{atom.resn}{atom.resi}:{atom.atom_name}"
            row[f"protein_{role}_element"] = atom.element
            if isinstance(val, list):
                row[f"protein_{role}_count"] = len(val)
    for name, value in state.features.items():
        row[name] = round(value, 6)
    for name, value in state.feature_scores.items():
        row[name] = round(value, 6)
    for level, value in state.level_scores.items():
        row[f"level_{level}_score"] = round(value, 6)
    return row


def _ligand_detail_rows(state: CandidateState, context: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for site, atom in state.ligand_sites.items():
        extra = atom.extra or {}
        row = dict(context)
        row.update({
            "ligand_site": site,
            "group_id": extra.get("group_id", ""),
            "instance_id": extra.get("instance_id", ""),
            "atom_label": extra.get("atom_label", ""),
            "atom_class": extra.get("atom_class", ""),
            "atom_role": extra.get("atom_role", ""),
            "pdbqt_order": atom.serial,
            "file_serial": extra.get("file_serial", ""),
            "element": atom.element,
            "x": round(atom.x, 6),
            "y": round(atom.y, 6),
            "z": round(atom.z, 6),
        })
        rows.append(row)
    return rows


def _protein_detail_rows(state: CandidateState, context: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for role, value in state.protein_roles.items():
        atoms = _role_values_as_atoms(value)
        for i, atom in enumerate(atoms, start=1):
            row = dict(context)
            row.update({
                "protein_role": role,
                "role_atom_index": i,
                "chain": atom.chain,
                "resi": atom.resi,
                "resn": atom.resn,
                "atom_name": atom.atom_name,
                "element": atom.element,
                "x": round(atom.x, 6),
                "y": round(atom.y, 6),
                "z": round(atom.z, 6),
                "distance_to_anchor": round(state.role_distances.get(role, 0.0), 6) if role in state.role_distances else "",
            })
            rows.append(row)
    return rows


def _feature_meta(mechanism: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {name: spec for name, spec in _features_iter(mechanism)}


def _geometry_detail_rows(state: CandidateState, context: dict[str, Any], mechanism: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    meta = _feature_meta(mechanism)
    for feature_name, value in state.features.items():
        spec = meta.get(feature_name, {})
        row = dict(context)
        row.update({
            "feature_name": feature_name,
            "feature_type": spec.get("type", ""),
            "value": round(value, 6),
            "unit": spec.get("unit", "angstrom" if spec.get("type") == "distance" else "degree"),
            "level": int((spec.get("score") or {}).get("level", spec.get("level", 1))),
            "score_value": round(state.feature_scores.get(f"{feature_name}_score", 0.0), 6) if f"{feature_name}_score" in state.feature_scores else "",
        })
        rows.append(row)
    return rows


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = list(fieldnames or [])
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    if not fields:
        path.write_text("", encoding="utf-8")
        return
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


CANDIDATE_FIELDS = [
    "protein_id", "conformation_id", "conformation_name", "conformation_index",
    "ligand_id", "pose_id", "site_set_id", "affinity_kcal",
    "candidate_score", "affinity_score", "geometry_score",
]
LIGAND_DETAIL_FIELDS = [
    "protein_id", "conformation_id", "conformation_name", "conformation_index",
    "ligand_id", "pose_id", "site_set_id", "affinity_kcal",
    "ligand_site", "group_id", "instance_id", "atom_label", "atom_class",
    "atom_role", "pdbqt_order", "file_serial", "element", "x", "y", "z",
]
PROTEIN_DETAIL_FIELDS = [
    "protein_id", "conformation_id", "conformation_name", "conformation_index",
    "ligand_id", "pose_id", "site_set_id", "affinity_kcal",
    "protein_role", "role_atom_index", "chain", "resi", "resn", "atom_name",
    "element", "x", "y", "z", "distance_to_anchor",
]
GEOMETRY_DETAIL_FIELDS = [
    "protein_id", "conformation_id", "conformation_name", "conformation_index",
    "ligand_id", "pose_id", "site_set_id", "affinity_kcal",
    "feature_name", "feature_type", "value", "unit", "level", "score_value",
]
POSE_SCORE_FIELDS = [
    "protein_id", "conformation_id", "conformation_name", "conformation_index",
    "ligand_id", "pose_id", "pose_score", "best_site_set_id",
]
CONFORMATION_SCORE_FIELDS = [
    "protein_id", "conformation_id", "conformation_name", "conformation_index",
    "ligand_id", "conformation_score", "best_pose_id",
]
PROTEIN_SCORE_FIELDS = [
    "protein_id", "ligand_id", "protein_score", "quality_score", "coverage",
    "active_conformations", "total_conformations", "observed_conformations",
]


def run_single_pair(mechanism: dict[str, Any], profile: dict[str, Any], protein_file: Path, docked_pdbqt: Path, role_table: Path, out_dir: Path, mechanism_path: Path | None = None, vina_out: Path | None = None) -> dict[str, Any]:
    protein_atoms = read_structure_atoms(protein_file, source="protein")
    role_rows = read_role_table(role_table)
    poses = parse_pdbqt_poses(docked_pdbqt)
    affinities = parse_vina_out(vina_out) if vina_out else {}
    template_atoms_by_key = _load_template_atoms(mechanism, mechanism_path, protein_file, protein_atoms)
    all_rows: list[dict[str, Any]] = []
    ligand_detail_rows: list[dict[str, Any]] = []
    protein_detail_rows: list[dict[str, Any]] = []
    geometry_detail_rows: list[dict[str, Any]] = []
    protein_name = protein_file.stem
    protein_id, conformation_id, conformation_index = split_protein_conformation_name(protein_name)
    ligand_id = role_rows[0].ligand_id if role_rows else docked_pdbqt.stem.split("@")[0]
    for pose in poses:
        pose_id = pose["pose_id"]
        affinity = affinities.get(pose_id, pose.get("affinity_kcal"))
        ligand_site_maps = _read_ligand_sites(mechanism, role_rows, pose["atoms"])
        initial_states = [CandidateState(sites) for sites in ligand_site_maps]
        states = _expand_protein_roles(mechanism, protein_atoms, initial_states, template_atoms_by_key)
        candidate_idx = 0
        for state in states:
            features = _compute_features(mechanism, state)
            if not _passes_gate(mechanism, features):
                continue
            _score_candidate(mechanism, profile, state, affinity)
            candidate_idx += 1
            state.site_set_id = f"pose_{pose_id}_set_{candidate_idx}"
            context = {
                "protein_id": protein_id,
                "conformation_id": conformation_id,
                "conformation_name": protein_name,
                "conformation_index": conformation_index,
                "ligand_id": ligand_id,
                "pose_id": pose_id,
                "site_set_id": state.site_set_id,
                "affinity_kcal": affinity if affinity is not None else "",
            }
            row = _state_row(state, context)
            all_rows.append(row)
            ligand_detail_rows.extend(_ligand_detail_rows(state, context))
            protein_detail_rows.extend(_protein_detail_rows(state, context))
            geometry_detail_rows.extend(_geometry_detail_rows(state, context, mechanism))
    _write_csv(out_dir / "resolved_ligand_sites.csv", ligand_detail_rows, LIGAND_DETAIL_FIELDS)
    _write_csv(out_dir / "resolved_protein_roles.csv", protein_detail_rows, PROTEIN_DETAIL_FIELDS)
    _write_csv(out_dir / "geometry_features.csv", geometry_detail_rows, GEOMETRY_DETAIL_FIELDS)
    _write_csv(out_dir / "passing_candidates.csv", all_rows, CANDIDATE_FIELDS)
    _write_csv(out_dir / "candidate_scores.csv", all_rows, CANDIDATE_FIELDS)

    pose_rows: list[dict[str, Any]] = []
    by_pose: dict[tuple, list[dict[str, Any]]] = {}
    for row in all_rows:
        key = (row["protein_id"], row["conformation_id"], row["ligand_id"], row["pose_id"])
        by_pose.setdefault(key, []).append(row)
    for key, rows in by_pose.items():
        best = max(rows, key=lambda r: float(r["candidate_score"]))
        pose_rows.append({"protein_id": key[0], "conformation_id": key[1], "ligand_id": key[2], "pose_id": key[3], "pose_score": best["candidate_score"], "best_site_set_id": best["site_set_id"]})
    _write_csv(out_dir / "pose_scores.csv", pose_rows, POSE_SCORE_FIELDS)

    conf_rows: list[dict[str, Any]] = []
    by_conf: dict[tuple, list[dict[str, Any]]] = {}
    for row in pose_rows:
        key = (row["protein_id"], row["conformation_id"], row["ligand_id"])
        by_conf.setdefault(key, []).append(row)
    for key, rows in by_conf.items():
        best = max(rows, key=lambda r: float(r["pose_score"]))
        conf_rows.append({
            "protein_id": key[0],
            "conformation_id": key[1],
            "conformation_name": best.get("conformation_name", ""),
            "conformation_index": best.get("conformation_index", ""),
            "ligand_id": key[2],
            "conformation_score": best["pose_score"],
            "best_pose_id": best["pose_id"],
        })
    _write_csv(out_dir / "conformation_scores.csv", conf_rows, CONFORMATION_SCORE_FIELDS)

    return {
        "candidates": len(all_rows),
        "poses": len(pose_rows),
        "conformations": len(conf_rows),
        "proteins": 0,
        "conformation_rows": conf_rows,
    }


def run_directory(mechanism_path: Path, profile_path: Path, protein_dir: Path, docking_dir: Path, role_table_dir: Path, out_dir: Path, ligand_id: str | None = None) -> dict[str, Any]:
    mechanism = load_yaml(mechanism_path)
    profile = load_yaml(profile_path)
    ligand_ext = profile.get("input", {}).get("ligand_extension", ".pdbqt")
    ligand_pattern = profile.get("input", {}).get("ligand_pose_pattern", "{ligand_id}@{protein_name}{ligand_extension}")
    out_pattern = profile.get("input", {}).get("docking_output_pattern", "{ligand_id}_{protein_name}_cavity_1.out")
    role_pattern = profile.get("input", {}).get("ligand_role_table_pattern", "{ligand_id}.role_table.csv")
    protein_files = sorted(list(protein_dir.glob("*.pdbqt")) + list(protein_dir.glob("*.pdb")))
    if ligand_id:
        ligand_ids = [ligand_id]
    else:
        ligand_ids = sorted(p.name.replace(".role_table.csv", "") for p in role_table_dir.glob("*.role_table.csv"))
    summary = {"pairs": 0, "candidates": 0, "proteins": 0}
    all_conformation_rows: list[dict[str, Any]] = []
    for pid in ligand_ids:
        role_table = role_table_dir / role_pattern.format(ligand_id=pid)
        if not role_table.exists():
            continue
        for protein_file in protein_files:
            protein_name = protein_file.stem
            docked = docking_dir / _format_pattern(ligand_pattern, pid, protein_name, ligand_ext)
            if not docked.exists():
                continue
            vina_out = docking_dir / out_pattern.format(ligand_id=pid, protein_name=protein_name)
            pair_out = out_dir / _safe_name(pid) / _safe_name(protein_name)
            result = run_single_pair(mechanism, profile, protein_file, docked, role_table, pair_out, mechanism_path, vina_out if vina_out.exists() else None)
            summary["pairs"] += 1
            summary["candidates"] += int(result.get("candidates", 0))
            all_conformation_rows.extend(result.get("conformation_rows", []))

    _write_csv(out_dir / "merged_conformation_scores.csv", all_conformation_rows, CONFORMATION_SCORE_FIELDS)

    protein_rows: list[dict[str, Any]] = []
    by_protein: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in all_conformation_rows:
        key = (str(row["protein_id"]), str(row["ligand_id"]))
        by_protein.setdefault(key, []).append(row)

    agg_cfg = profile.get("aggregation", {}).get("protein", {}) or {}
    total_confs = int(agg_cfg.get("total_conformations", profile.get("runtime", {}).get("total_conformations", 1)))
    threshold = float(agg_cfg.get("coverage_threshold", 70.0))
    target = float(agg_cfg.get("coverage_target", 0.70))
    alpha = float(agg_cfg.get("alpha", 0.30))
    for key, rows in by_protein.items():
        scores = [float(r["conformation_score"]) for r in rows]
        quality = max(scores) if scores else 0.0
        active = sum(1 for s in scores if s >= threshold)
        observed = len({str(r.get("conformation_id", "")) for r in rows})
        coverage = active / float(total_confs) if total_confs else 0.0
        coverage_factor = (1.0 - alpha) + alpha * min(coverage / target, 1.0) if target > 0 else 1.0
        protein_rows.append({
            "protein_id": key[0],
            "ligand_id": key[1],
            "protein_score": round(quality * coverage_factor, 6),
            "quality_score": round(quality, 6),
            "coverage": round(coverage, 6),
            "active_conformations": active,
            "total_conformations": total_confs,
            "observed_conformations": observed,
        })
    _write_csv(out_dir / "protein_scores.csv", protein_rows, PROTEIN_SCORE_FIELDS)
    summary["proteins"] = len(protein_rows)
    return summary
