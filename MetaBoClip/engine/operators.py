from __future__ import annotations

import math
from typing import Any

from core.pymol_utils import dist


OBJECT_OPS = {
    "select_first_atom",
    "select_first_atom_any",
    "select_atoms",
    "select_atoms_within",
    "select_nearest_residue",
    "select_atom_from_residue",
    "select_nearest_from_collection",
    "select_atom_triplet",
    "item_from_collection",
    "plane_normal",
    "oriented_normal",
    "pose_clash",
    "collect_hydroxyl_sites",
    "collect_protic_nucleophiles",
    "collect_ch_pairs",
    "collect_carbonyl_sites",
}

FEATURE_OPS = {
    "distance",
    "angle",
    "axis_deviation",
    "alias",
}


def is_ref_object(value: Any) -> bool:
    return isinstance(value, dict) and "__kind__" in value


def is_missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, list) and len(value) == 0:
        return True
    if isinstance(value, dict) and not value:
        return True
    return False


def format_selector(template: str, selector_vars: dict[str, Any]) -> str:
    if not isinstance(template, str):
        return template
    safe_vars = {k: (v if isinstance(v, (str, int, float)) else str(v)) for k, v in selector_vars.items()}
    return template.format(**safe_vars)


def atom_ref_from_pymol(atom, object_name: str | None = None) -> dict[str, Any]:
    symbol = getattr(atom, "symbol", None) or getattr(atom, "elem", None) or ""
    atom_id = int(getattr(atom, "id", 0))
    # PyMOL keeps two useful numeric identifiers:
    #   - index: session/object atom index, used by the original ACH script and PSE selections
    #   - id:    PDB/PDBQT atom serial, stable across ligand conformations when present
    pymol_index = int(getattr(atom, "index", atom_id))
    return {
        "__kind__": "atom",
        "object": object_name,
        "atom_id": atom_id,
        "pymol_index": pymol_index,
        "index": pymol_index,
        "name": str(atom.name),
        "symbol": str(symbol),
        "resn": str(getattr(atom, "resn", "")),
        "resi": str(getattr(atom, "resi", "")),
        "chain": str(getattr(atom, "chain", "")),
        "coord": [float(atom.coord[0]), float(atom.coord[1]), float(atom.coord[2])],
    }


def residue_ref(object_name: str, chain: str, resi: str, resn: str) -> dict[str, Any]:
    chain_part = f"chain {chain} and " if chain and str(chain).strip() else ""
    selector = f"({object_name} and {chain_part}resi {resi} and resn {resn})"
    return {
        "__kind__": "residue",
        "object": object_name,
        "chain": str(chain),
        "resi": str(resi),
        "resn": str(resn),
        "selector": selector,
    }


def vector_ref(values: list[float]) -> dict[str, Any]:
    return {"__kind__": "vector", "value": [float(v) for v in values]}


def coord_of(obj: Any) -> list[float]:
    if isinstance(obj, dict):
        kind = obj.get("__kind__")
        if kind in {"atom", "point"}:
            return [float(x) for x in obj["coord"]]
    raise RuntimeError(f"Object does not have coordinates: {obj}")


def vector_of(obj: Any) -> list[float]:
    if isinstance(obj, dict) and obj.get("__kind__") == "vector":
        return [float(x) for x in obj["value"]]
    raise RuntimeError(f"Object is not a vector: {obj}")


def attr_of(obj: Any, attr: str) -> Any:
    if is_ref_object(obj):
        if attr in obj:
            return obj[attr]
        if attr == "id":
            return obj.get("atom_id")
        if attr == "element":
            return obj.get("symbol")
    if isinstance(obj, dict) and attr in obj:
        return obj[attr]
    raise RuntimeError(f"Attribute '{attr}' not found on object {obj}")


def _selected_atoms(cmd, selector: str) -> list[dict[str, Any]]:
    try:
        if cmd.count_atoms(selector) <= 0:
            return []
        model = cmd.get_model(selector, 1)
    except Exception:
        return []
    return [atom_ref_from_pymol(a, object_name=str(getattr(a, "model", ""))) for a in model.atom]


def _first_atom_from_selector(cmd, selector: str, nth: int = 0) -> dict[str, Any] | None:
    atoms = _selected_atoms(cmd, selector)
    if len(atoms) <= nth:
        return None
    return atoms[nth]


def _nearest(items: list[Any], reference: Any, max_distance: float | None = None) -> Any:
    ref_xyz = coord_of(reference)
    best = None
    for item in items:
        d = dist(coord_of(item), ref_xyz)
        if max_distance is not None and d > float(max_distance):
            continue
        if best is None or d < best[0]:
            best = (d, item)
    return best[1] if best is not None else None


def _plane_normal(a_xyz, b_xyz, c_xyz) -> list[float]:
    v1 = [b_xyz[i] - a_xyz[i] for i in range(3)]
    v2 = [c_xyz[i] - a_xyz[i] for i in range(3)]
    n = [
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0],
    ]
    nn = math.sqrt(sum(x * x for x in n)) + 1e-12
    return [x / nn for x in n]


def _axis_deviation(axis, vec) -> float:
    nv = math.sqrt(sum(x * x for x in vec)) + 1e-12
    cosv = max(-1.0, min(1.0, sum(axis[i] * vec[i] for i in range(3)) / nv))
    angle = math.degrees(math.acos(cosv))
    return min(angle, 180.0 - angle)


def resolve_name(name: str, named: dict[str, Any]) -> Any:
    if name not in named:
        raise RuntimeError(f"Unknown reference '{name}'. Available: {sorted(named)}")
    return named[name]


def execute_object_op(cmd, op: str, spec: dict[str, Any], named: dict[str, Any], selector_vars: dict[str, Any]) -> Any:
    if op == "select_first_atom":
        selector = format_selector(spec["selector"], selector_vars)
        nth = int(spec.get("nth", 0))
        return _first_atom_from_selector(cmd, selector, nth=nth)

    if op == "select_first_atom_any":
        selectors = [format_selector(s, selector_vars) for s in spec.get("selectors", [])]
        if not selectors:
            selectors = [format_selector(s, selector_vars) for s in spec.get("selector_any", [])]
        nth = int(spec.get("nth", 0))
        for selector in selectors:
            atom = _first_atom_from_selector(cmd, selector, nth=nth)
            if atom is not None:
                return atom
        return None

    if op == "select_atoms":
        selector = format_selector(spec["selector"], selector_vars)
        return _selected_atoms(cmd, selector)

    if op == "select_atoms_within":
        selector = format_selector(spec["candidate_selector"], selector_vars)
        reference = resolve_name(spec["reference"], named)
        radius = float(spec["radius"])
        ref_xyz = coord_of(reference)
        return [a for a in _selected_atoms(cmd, selector) if dist(a["coord"], ref_xyz) <= radius]

    if op == "select_nearest_residue":
        selector = format_selector(spec["candidate_selector"], selector_vars)
        atom_name = spec.get("residue_atom_name", "CA")
        if f"name {atom_name}" not in selector:
            selector = f"({selector}) and name {atom_name}"
        reference = resolve_name(spec["reference"], named)
        ref_xyz = coord_of(reference)
        atoms = _selected_atoms(cmd, selector)
        best = None
        for atom in atoms:
            d = dist(atom["coord"], ref_xyz)
            if best is None or d < best[0]:
                best = (d, atom)
        if best is None:
            return None
        atom = best[1]
        return residue_ref(atom["object"], atom["chain"], atom["resi"], atom["resn"])

    if op == "select_atom_from_residue":
        residue = resolve_name(spec["residue"], named)
        atom_name = spec["atom_name"]
        selector = f"({residue['selector']}) and name {atom_name}"
        return _first_atom_from_selector(cmd, selector)

    if op == "select_nearest_from_collection":
        collection = resolve_name(spec["collection"], named)
        reference = resolve_name(spec["reference"], named)
        max_distance = spec.get("max_distance")
        return _nearest(collection, reference, max_distance=max_distance)

    if op == "select_atom_triplet":
        selectors = [format_selector(s, selector_vars) for s in spec.get("preferred_selectors", [])]
        atoms = []
        seen_ids = set()
        for selector in selectors:
            atom = _first_atom_from_selector(cmd, selector)
            if atom is None:
                continue
            if atom["atom_id"] in seen_ids:
                continue
            atoms.append(atom)
            seen_ids.add(atom["atom_id"])
        if len(atoms) < 3 and spec.get("fallback_selector"):
            fallback = format_selector(spec["fallback_selector"], selector_vars)
            fallback_atoms = _selected_atoms(cmd, fallback)
            for atom in fallback_atoms:
                if atom["atom_id"] in seen_ids:
                    continue
                atoms.append(atom)
                seen_ids.add(atom["atom_id"])
                if len(atoms) >= 3:
                    break
        return atoms[:3]

    if op == "item_from_collection":
        collection = resolve_name(spec["collection"], named)
        idx = int(spec["index"])
        if len(collection) <= idx:
            return None
        return collection[idx]

    if op == "plane_normal":
        a = resolve_name(spec["a"], named)
        b = resolve_name(spec["b"], named)
        c = resolve_name(spec["c"], named)
        return vector_ref(_plane_normal(coord_of(a), coord_of(b), coord_of(c)))

    if op == "oriented_normal":
        normal = resolve_name(spec["normal"], named)
        origin = resolve_name(spec["origin"], named)
        reference = resolve_name(spec["reference"], named)
        axis = vector_of(normal)
        origin_xyz = coord_of(origin)
        ref_xyz = coord_of(reference)
        ref_vec = [ref_xyz[i] - origin_xyz[i] for i in range(3)]
        dot = sum(axis[i] * ref_vec[i] for i in range(3))
        if dot > 0:
            axis = [-x for x in axis]
        return vector_ref(axis)

    if op == "pose_clash":
        pose_obj = selector_vars["pose_obj"]
        protein_obj = selector_vars["protein_obj"]
        cutoff = float(spec.get("cutoff", 2.0))
        clash_sel = f"({protein_obj} and not hydro) within {cutoff} of ({pose_obj} and not hydro)"
        try:
            return bool(cmd.count_atoms(clash_sel) > 0)
        except Exception:
            return False

    if op == "collect_hydroxyl_sites":
        pose_obj = selector_vars["pose_obj"]
        selector = spec.get("selector") or f"({pose_obj} and name O) and (neighbor ({pose_obj} and name H))"
        selector = format_selector(selector, selector_vars)
        return _selected_atoms(cmd, selector)

    if op == "collect_protic_nucleophiles":
        pose_obj = selector_vars["pose_obj"]
        elements = spec.get("elements", ["O", "N", "S"])
        atoms = []
        for elem in elements:
            selector = f"({pose_obj} and elem {elem}) and (neighbor ({pose_obj} and elem H))"
            atoms.extend(_selected_atoms(cmd, selector))
        uniq = {}
        for atom in atoms:
            uniq[atom["atom_id"]] = atom
        return list(uniq.values())

    if op == "collect_ch_pairs":
        pose_obj = selector_vars["pose_obj"]
        cutoff = float(spec.get("h_cutoff", 1.25))
        model = cmd.get_model(pose_obj)
        carbons = [a for a in model.atom if (getattr(a, "symbol", None) or getattr(a, "elem", None)) == "C"]
        hydrogens = [a for a in model.atom if (getattr(a, "symbol", None) or getattr(a, "elem", None)) == "H"]
        pairs = []
        for h in hydrogens:
            h_ref = atom_ref_from_pymol(h, object_name=pose_obj)
            best = None
            for c in carbons:
                c_ref = atom_ref_from_pymol(c, object_name=pose_obj)
                d = dist(h_ref["coord"], c_ref["coord"])
                if d <= cutoff and (best is None or d < best[0]):
                    best = (d, c_ref)
            if best is not None:
                pairs.append({"lig_h_atom": h_ref, "parent_atom": best[1]})
        return pairs

    if op == "collect_carbonyl_sites":
        pose_obj = selector_vars["pose_obj"]
        cutoff = float(spec.get("neighbor_cutoff", 2.0))
        required_oxygens = int(spec.get("exact_oxygens", spec.get("required_oxygens", spec.get("min_oxygens", 2))))
        oxygen_count_mode = str(spec.get("oxygen_count_mode", "exact")).lower()
        require_o_bond_neighbor = bool(spec.get("require_o_bond_neighbor", True))

        # Match the original ACH gating script:
        #   1) start from carbons that PyMOL sees as bonded neighbors of oxygen
        #   2) within NEAR_O_RADIUS around each carbon, require exactly two O atoms
        #   3) use the nearest O for the Bürgi-Dunitz angle
        carbon_selector = f"({pose_obj} and elem C)"
        if require_o_bond_neighbor:
            carbon_selector = f"({pose_obj} and elem C and neighbor ({pose_obj} and elem O))"

        try:
            carbon_model = cmd.get_model(carbon_selector)
            oxygen_model = cmd.get_model(f"{pose_obj} and elem O")
            ligand_model = cmd.get_model(pose_obj)
        except Exception:
            return []

        carbons = [atom_ref_from_pymol(a, object_name=pose_obj) for a in carbon_model.atom]
        oxygens = [atom_ref_from_pymol(a, object_name=pose_obj) for a in oxygen_model.atom]
        ligand_atom_ids = ",".join(str(int(getattr(a, "index", getattr(a, "id", 0)))) for a in ligand_model.atom)
        sites = []
        for idx, carbon in enumerate(carbons, start=1):
            o_neighbors = [o for o in oxygens if dist(carbon["coord"], o["coord"]) <= cutoff]
            if oxygen_count_mode in {"exact", "eq", "=="}:
                if len(o_neighbors) != required_oxygens:
                    continue
            else:
                if len(o_neighbors) < required_oxygens:
                    continue

            chosen = min(o_neighbors, key=lambda o: dist(carbon["coord"], o["coord"]))
            sites.append({
                "carbonyl_c": carbon,
                "carbonyl_o": chosen,
                "carbon_index": carbon.get("pymol_index", carbon.get("atom_id")),
                "carbon_id": carbon.get("atom_id"),
                "selection_name": f"{pose_obj}_Csel_{idx}",
                "ligand_atom_ids": ligand_atom_ids,
            })
        return sites

    raise RuntimeError(f"Unknown object operator: {op}")


def compute_feature(op: str, spec: dict[str, Any], named: dict[str, Any]) -> Any:
    if op == "distance":
        a = resolve_name(spec["a"], named)
        b = resolve_name(spec["b"], named)
        return float(dist(coord_of(a), coord_of(b)))

    if op == "angle":
        a = resolve_name(spec["a"], named)
        b = resolve_name(spec["b"], named)
        c = resolve_name(spec["c"], named)
        a_xyz = coord_of(a)
        b_xyz = coord_of(b)
        c_xyz = coord_of(c)
        bax = a_xyz[0] - b_xyz[0]
        bay = a_xyz[1] - b_xyz[1]
        baz = a_xyz[2] - b_xyz[2]
        bcx = c_xyz[0] - b_xyz[0]
        bcy = c_xyz[1] - b_xyz[1]
        bcz = c_xyz[2] - b_xyz[2]
        dot = bax * bcx + bay * bcy + baz * bcz
        na = math.sqrt(bax * bax + bay * bay + baz * baz) + 1e-12
        nb = math.sqrt(bcx * bcx + bcy * bcy + bcz * bcz) + 1e-12
        cosv = max(-1.0, min(1.0, dot / (na * nb)))
        return float(math.degrees(math.acos(cosv)))

    if op == "axis_deviation":
        axis = vector_of(resolve_name(spec["axis"], named))
        origin = resolve_name(spec["origin"], named)
        target = resolve_name(spec["target"], named)
        vec = [coord_of(target)[i] - coord_of(origin)[i] for i in range(3)]
        return float(_axis_deviation(axis, vec))

    if op == "alias":
        return resolve_name(spec["source"], named)

    raise RuntimeError(f"Unknown feature operator: {op}")
