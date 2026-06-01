from __future__ import annotations

from pathlib import Path
import csv
import re
from typing import Any

try:
    import yaml
except Exception:  # pragma: no cover
    yaml = None


PML_HEADER = """reinitialize
set retain_order, 1
set auto_zoom, off
set sphere_scale, 0.35
set stick_radius, 0.16
set dash_width, 2.0
set dash_gap, 0.20
set label_size, 16
"""


def _read_rows(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _safe_name(value: str | int | float | None) -> str:
    text = re.sub(r"[^A-Za-z0-9_]+", "_", str(value or "").strip())
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "item"


def _abs_path(path: Path | str | None, base_dir: Path | None = None) -> Path | None:
    if path is None:
        return None
    p = Path(path)
    if not p.is_absolute() and base_dir is not None:
        p = base_dir / p
    return p.expanduser().resolve()


def _require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    if path.is_dir():
        raise IsADirectoryError(f"{label} is a directory: {path}")


def _quote_path(path: Path | str) -> str:
    text = str(path).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{text}"'


def _quote_text(value: str) -> str:
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def _parse_pdbqt_atom_line(line: str, serial: int) -> str:
    atom_name = line[12:16].strip() if len(line) >= 16 else "C"
    resn = line[17:20].strip() if len(line) >= 20 else "LIG"
    chain = line[21:22].strip() if len(line) >= 22 else "L"
    resi = line[22:26].strip() if len(line) >= 26 else "1"
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except Exception:
        parts = line.split()
        x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
    element = ""
    if len(line) >= 78:
        element = line[76:78].strip()
    if not element:
        parts = line.split()
        raw = parts[-1] if parts else atom_name
        raw = raw.replace(".", "")
        upper = raw.upper()
        if upper in {"OA", "OS"}:
            element = "O"
        elif upper in {"NA", "NS"}:
            element = "N"
        elif upper == "A":
            element = "C"
        elif len(raw) >= 2 and raw[:2].capitalize() in {"Cl", "Br", "Fe", "Mg", "Zn", "Mn", "Cu", "Co", "Ni", "Ca", "Na", "K"}:
            element = raw[:2].capitalize()
        elif raw:
            element = raw[0].upper()
    return f"HETATM{serial:5d} {atom_name:<4s} {resn:>3s} {chain or 'L':1s}{str(resi):>4s}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"


def _pdbqt_to_pdb(in_file: Path, out_pdb: Path) -> None:
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    serial = 0
    with open(in_file, "r", encoding="utf-8", errors="ignore") as inp, open(out_pdb, "w", encoding="utf-8") as out:
        for line in inp:
            if line.startswith(("ATOM", "HETATM")):
                serial += 1
                out.write(_parse_pdbqt_atom_line(line, serial))
        out.write("END\n")


def _extract_pose_to_pdb(docked_pdbqt: Path, pose_id: int, out_pdb: Path) -> None:
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    current_pose = 0
    saw_model = False
    selected_lines: list[str] = []
    one_pose_lines: list[str] = []
    with open(docked_pdbqt, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                saw_model = True
                current_pose += 1
                one_pose_lines = []
                continue
            if line.startswith("ENDMDL"):
                if current_pose == pose_id:
                    selected_lines = one_pose_lines
                    break
                one_pose_lines = []
                continue
            if line.startswith(("ATOM", "HETATM")):
                if saw_model:
                    one_pose_lines.append(line)
                else:
                    selected_lines.append(line)
        if not saw_model and pose_id != 1:
            selected_lines = []
    if not selected_lines:
        raise ValueError(f"Pose {pose_id} was not found in {docked_pdbqt}")
    with open(out_pdb, "w", encoding="utf-8") as out:
        for serial, line in enumerate(selected_lines, start=1):
            out.write(_parse_pdbqt_atom_line(line, serial))
        out.write("END\n")


def _choose_candidate(
    candidate_rows: list[dict[str, str]],
    site_set_id: str | None = None,
    pose_id: str | None = None,
    mode: str = "best",
) -> dict[str, str] | None:
    rows = list(candidate_rows)
    if site_set_id:
        rows = [r for r in rows if str(r.get("site_set_id", "")) == str(site_set_id)]
    if pose_id:
        rows = [r for r in rows if str(r.get("pose_id", "")) == str(pose_id)]
    if not rows:
        return None
    if mode == "first":
        return rows[0]
    def score(row: dict[str, str]) -> float:
        try:
            return float(row.get("candidate_score", "nan"))
        except Exception:
            return float("nan")
    scored = [(score(r), i, r) for i, r in enumerate(rows)]
    scored = [(s, i, r) for s, i, r in scored if s == s]
    if not scored:
        return rows[0]
    return max(scored, key=lambda x: (x[0], -x[1]))[2]


def _filter_context_rows(rows: list[dict[str, str]], candidate: dict[str, str] | None) -> list[dict[str, str]]:
    if candidate is None:
        return rows
    pose = str(candidate.get("pose_id", ""))
    site = str(candidate.get("site_set_id", ""))
    out = []
    for row in rows:
        if pose and str(row.get("pose_id", "")) != pose:
            continue
        if site and str(row.get("site_set_id", "")) != site:
            continue
        out.append(row)
    return out


def _ligand_selection_name(row: dict[str, str], index: int) -> str:
    site = _safe_name(row.get("ligand_site") or "ligand")
    label = _safe_name(row.get("atom_label") or row.get("atom_role") or str(index))
    inst = _safe_name(row.get("instance_id") or str(index))
    return f"LIG_{site}_{label}_{inst}_{index}"


def _protein_selection_name(row: dict[str, str], index: int) -> str:
    role = _safe_name(row.get("protein_role") or row.get("role") or f"role_{index}")
    atom = _safe_name(row.get("atom_name") or str(index))
    resn = _safe_name(row.get("resn") or "res")
    resi = _safe_name(row.get("resi") or "")
    return f"PROT_{role}_{resn}{resi}_{atom}_{index}"


def _protein_selection_expression(row: dict[str, str]) -> str:
    parts = ["protein"]
    chain = str(row.get("chain", "")).strip()
    resi = str(row.get("resi", "")).strip()
    atom = str(row.get("atom_name", "")).strip()
    resn = str(row.get("resn", "")).strip()
    if chain:
        parts.append(f"chain {chain}")
    if resi:
        parts.append(f"resi {resi}")
    if resn:
        parts.append(f"resn {resn}")
    if atom:
        parts.append(f"name {atom}")
    return " and ".join(parts)


def _load_mechanism(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    _require_file(path, "mechanism")
    if yaml is None:
        raise RuntimeError("PyYAML is required to use --mechanism in export-pymol")
    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return data if isinstance(data, dict) else {}


def _feature_items(mechanism: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    features = mechanism.get("features") or mechanism.get("geometry_features", {}).get("features") or {}
    if isinstance(features, list):
        return [(str(f.get("name")), f) for f in features if isinstance(f, dict) and f.get("name")]
    if isinstance(features, dict):
        return [(str(k), v) for k, v in features.items() if isinstance(v, dict)]
    return []


def _ref_to_selection(ref: str, ligand_site_selections: dict[str, list[str]], protein_role_selections: dict[str, list[str]]) -> str | None:
    prefix, _, name = str(ref).partition(".")
    if prefix == "ligand":
        sels = ligand_site_selections.get(name) or []
    elif prefix == "protein":
        sels = protein_role_selections.get(name) or []
    else:
        sels = []
    if not sels:
        return None
    return sels[0] if len(sels) == 1 else "(" + " or ".join(sels) + ")"


def _pml_select_ligand(name: str, order: str, label: str, coord: tuple[float, float, float] | None) -> list[str]:
    lines = []
    lines.append(f"select {name}, ligand and (id {order} or index {order})")
    lines.append(f"show spheres, {name}")
    lines.append(f"color red, {name}")
    lines.append(f"label {name}, \"{_quote_text(label)}\"")
    if coord is not None:
        marker = f"MARKER_{name}"
        x, y, z = coord
        lines.append("python")
        lines.append("from pymol import cmd")
        lines.append(f"if cmd.count_atoms('{name}') == 0:")
        lines.append(f"    cmd.pseudoatom('{marker}', pos=[{x:.3f}, {y:.3f}, {z:.3f}], label='{_quote_text(label)}')")
        lines.append(f"    cmd.show('spheres', '{marker}')")
        lines.append(f"    cmd.color('red', '{marker}')")
        lines.append("python end")
    return lines


def _pml_select_protein(name: str, expr: str, label: str, coord: tuple[float, float, float] | None) -> list[str]:
    lines = []
    lines.append(f"select {name}, {expr}")
    lines.append(f"show sticks, {name}")
    lines.append(f"show spheres, {name}")
    lines.append(f"color cyan, {name}")
    lines.append(f"label {name}, \"{_quote_text(label)}\"")
    if coord is not None:
        marker = f"MARKER_{name}"
        x, y, z = coord
        lines.append("python")
        lines.append("from pymol import cmd")
        lines.append(f"if cmd.count_atoms('{name}') == 0:")
        lines.append(f"    cmd.pseudoatom('{marker}', pos=[{x:.3f}, {y:.3f}, {z:.3f}], label='{_quote_text(label)}')")
        lines.append(f"    cmd.show('spheres', '{marker}')")
        lines.append(f"    cmd.color('cyan', '{marker}')")
        lines.append("python end")
    return lines


def _float_coord(row: dict[str, str]) -> tuple[float, float, float] | None:
    try:
        return float(row["x"]), float(row["y"]), float(row["z"])
    except Exception:
        return None




def _candidate_sort_score(row: dict[str, str]) -> float:
    for key in ("candidate_score", "pose_score", "conformation_score", "score"):
        try:
            value = float(row.get(key, "nan"))
            if value == value:
                return value
        except Exception:
            pass
    return float("-inf")


def _ranked_candidate_rows(rows: list[dict[str, str]], mode: str = "pose", top_n: int | None = None) -> list[tuple[int, dict[str, str]]]:
    valid = [r for r in rows if str(r.get("pose_id", "")).strip()]
    if mode == "pose":
        best_by_pose: dict[str, dict[str, str]] = {}
        for row in valid:
            pose = str(row.get("pose_id", "")).strip()
            current = best_by_pose.get(pose)
            if current is None or _candidate_sort_score(row) > _candidate_sort_score(current):
                best_by_pose[pose] = row
        valid = list(best_by_pose.values())
    elif mode != "candidate":
        raise ValueError(f"Unsupported export mode: {mode}")
    valid.sort(key=_candidate_sort_score, reverse=True)
    if top_n and top_n > 0:
        valid = valid[:top_n]
    return [(i, row) for i, row in enumerate(valid, start=1)]


def _candidate_file_stem(rank: int, row: dict[str, str]) -> str:
    pose = _safe_name(row.get("pose_id") or "pose")
    site = _safe_name(row.get("site_set_id") or "site")
    score = row.get("candidate_score") or row.get("score") or ""
    try:
        score_text = f"score_{float(score):.2f}"
    except Exception:
        score_text = "score_na"
    return f"top_{rank:03d}_pose_{pose}_site_{site}_{score_text}"


def export_pymol_script(
    protein_file: Path,
    docked_ligand_file: Path,
    resolved_ligand_sites: Path,
    resolved_protein_roles: Path,
    out_pml: Path,
    candidate_scores: Path | None = None,
    passing_candidates: Path | None = None,
    geometry_features: Path | None = None,
    mechanism_file: Path | None = None,
    candidate_site_set_id: str | None = None,
    pose_id: str | None = None,
    select_mode: str = "best",
    save_pse_name: str | None = None,
) -> None:
    out_pml = _abs_path(out_pml) or Path("candidate_view.pml").resolve()
    out_dir = out_pml.parent
    protein_file = _abs_path(protein_file)
    docked_ligand_file = _abs_path(docked_ligand_file)
    resolved_ligand_sites = _abs_path(resolved_ligand_sites)
    resolved_protein_roles = _abs_path(resolved_protein_roles)
    candidate_scores = _abs_path(candidate_scores) if candidate_scores else None
    passing_candidates = _abs_path(passing_candidates) if passing_candidates else None
    geometry_features = _abs_path(geometry_features) if geometry_features else None
    mechanism_file = _abs_path(mechanism_file) if mechanism_file else None
    save_pse = _abs_path(save_pse_name, out_dir) if save_pse_name else None

    assert protein_file is not None and docked_ligand_file is not None and resolved_ligand_sites is not None and resolved_protein_roles is not None
    _require_file(protein_file, "protein")
    _require_file(docked_ligand_file, "docked PDBQT")
    _require_file(resolved_ligand_sites, "resolved ligand sites")
    _require_file(resolved_protein_roles, "resolved protein roles")

    ligand_rows_all = _read_rows(resolved_ligand_sites)
    protein_rows_all = _read_rows(resolved_protein_roles)
    candidate_rows = _read_rows(candidate_scores) or _read_rows(passing_candidates)
    candidate = _choose_candidate(candidate_rows, candidate_site_set_id, pose_id, select_mode) if candidate_rows else None
    if candidate is None and (candidate_site_set_id or pose_id):
        candidate = {"site_set_id": candidate_site_set_id or "", "pose_id": pose_id or ""}
    if candidate is None:
        raise ValueError("No candidate was selected. Provide candidate_scores.csv, passing_candidates.csv, --site-set-id, or --pose-id.")

    ligand_rows = _filter_context_rows(ligand_rows_all, candidate)
    protein_rows = _filter_context_rows(protein_rows_all, candidate)
    if not ligand_rows:
        raise ValueError("No ligand rows matched the selected candidate")
    if not protein_rows:
        raise ValueError("No protein rows matched the selected candidate")

    selected_pose_id = pose_id or candidate.get("pose_id") or ligand_rows[0].get("pose_id") or "1"
    try:
        selected_pose_int = int(str(selected_pose_id))
    except Exception:
        selected_pose_int = 1

    selected_site_set_id = str(candidate.get("site_set_id", ""))
    out_dir.mkdir(parents=True, exist_ok=True)
    protein_pdb = out_dir / f"{out_pml.stem}_protein_view.pdb"
    ligand_pdb = out_dir / f"{out_pml.stem}_ligand_pose_{selected_pose_int}.pdb"
    _pdbqt_to_pdb(protein_file, protein_pdb)
    _extract_pose_to_pdb(docked_ligand_file, selected_pose_int, ligand_pdb)

    mechanism = _load_mechanism(mechanism_file)
    geom_rows = _filter_context_rows(_read_rows(geometry_features), candidate) if geometry_features else []
    geom_value_by_name = {r.get("feature_name", ""): r for r in geom_rows}

    ligand_site_selections: dict[str, list[str]] = {}
    ligand_label_selections: dict[str, list[str]] = {}
    protein_role_selections: dict[str, list[str]] = {}
    pml_lines: list[str] = []
    pml_lines.extend(PML_HEADER.rstrip().splitlines())
    pml_lines.append(f"load {_quote_path(protein_pdb)}, protein")
    pml_lines.append(f"load {_quote_path(ligand_pdb)}, ligand")
    pml_lines.append("select PROTEIN_STRUCTURE, protein")
    pml_lines.append("select LIGAND_POSE_SELECTED, ligand")
    pml_lines.append("hide everything")
    pml_lines.append("show lines, protein")
    pml_lines.append("show sticks, ligand")
    pml_lines.append("color gray70, protein")
    pml_lines.append("color orange, ligand")
    pml_lines.append("show sticks, LIGAND_POSE_SELECTED")
    pml_lines.append(f"# selected_pose_id = {selected_pose_int}")
    pml_lines.append(f"# selected_site_set_id = {selected_site_set_id}")

    individual_ligand_sels: list[str] = []
    for i, row in enumerate(ligand_rows, start=1):
        order = str(row.get("pdbqt_order", "")).strip()
        if not order:
            continue
        sel_name = _ligand_selection_name(row, i)
        site = row.get("ligand_site", "")
        label = row.get("atom_label") or row.get("atom_role") or site or sel_name
        if site:
            ligand_site_selections.setdefault(site, []).append(sel_name)
        if row.get("atom_label"):
            ligand_label_selections.setdefault(row.get("atom_label", ""), []).append(sel_name)
        individual_ligand_sels.append(sel_name)
        pml_lines.extend(_pml_select_ligand(sel_name, order, f"ligand.{site}:{label}", _float_coord(row)))

    individual_protein_sels: list[str] = []
    for i, row in enumerate(protein_rows, start=1):
        sel_name = _protein_selection_name(row, i)
        role = row.get("protein_role") or row.get("role") or ""
        if role:
            protein_role_selections.setdefault(role, []).append(sel_name)
        individual_protein_sels.append(sel_name)
        expr = _protein_selection_expression(row)
        res_label = f"protein.{role}:{row.get('resn','')}{row.get('resi','')} {row.get('atom_name','')}"
        pml_lines.extend(_pml_select_protein(sel_name, expr, res_label, _float_coord(row)))

    for site, sels in ligand_site_selections.items():
        if sels:
            name = f"LIGAND_SITE_{_safe_name(site)}"
            pml_lines.append(f"select {name}, " + " or ".join(sels))
            pml_lines.append(f"show spheres, {name}")
            pml_lines.append(f"color red, {name}")
    for atom_label, sels in ligand_label_selections.items():
        if sels:
            name = f"LIGAND_LABEL_{_safe_name(atom_label)}"
            pml_lines.append(f"select {name}, " + " or ".join(sels))
    for role, sels in protein_role_selections.items():
        if sels:
            name = f"PROTEIN_ROLE_{_safe_name(role)}"
            pml_lines.append(f"select {name}, " + " or ".join(sels))
            pml_lines.append(f"show sticks, {name}")
            pml_lines.append(f"show spheres, {name}")
            pml_lines.append(f"color cyan, {name}")

    if individual_ligand_sels:
        pml_lines.append("select CANDIDATE_LIGAND_ATOMS, " + " or ".join(individual_ligand_sels))
        pml_lines.append("show spheres, CANDIDATE_LIGAND_ATOMS")
        pml_lines.append("color red, CANDIDATE_LIGAND_ATOMS")
        pml_lines.append("create CANDIDATE_LIGAND_ATOMS_OBJECT, CANDIDATE_LIGAND_ATOMS")
        pml_lines.append("show spheres, CANDIDATE_LIGAND_ATOMS_OBJECT")
    if individual_protein_sels:
        pml_lines.append("select CANDIDATE_PROTEIN_ATOMS, " + " or ".join(individual_protein_sels))
        pml_lines.append("show sticks, CANDIDATE_PROTEIN_ATOMS")
        pml_lines.append("show spheres, CANDIDATE_PROTEIN_ATOMS")
        pml_lines.append("color cyan, CANDIDATE_PROTEIN_ATOMS")
        pml_lines.append("create CANDIDATE_PROTEIN_ATOMS_OBJECT, CANDIDATE_PROTEIN_ATOMS")
        pml_lines.append("show spheres, CANDIDATE_PROTEIN_ATOMS_OBJECT")
    if individual_ligand_sels or individual_protein_sels:
        pml_lines.append("select CANDIDATE_ALL_ATOMS, " + " or ".join(individual_ligand_sels + individual_protein_sels))

    for feature_name, spec in _feature_items(mechanism):
        if spec.get("enabled", True) is False:
            continue
        ftype = spec.get("type")
        value_row = geom_value_by_name.get(feature_name, {})
        value_text = value_row.get("value", "")
        if ftype == "distance":
            atoms = spec.get("atoms") or []
            if len(atoms) != 2:
                continue
            sel_a = _ref_to_selection(atoms[0], ligand_site_selections, protein_role_selections)
            sel_b = _ref_to_selection(atoms[1], ligand_site_selections, protein_role_selections)
            if sel_a and sel_b:
                name = f"DIST_{_safe_name(feature_name)}"
                pml_lines.append(f"distance {name}, {sel_a}, {sel_b}")
                pml_lines.append(f"hide labels, {name}")
                if value_text:
                    pml_lines.append(f"# {feature_name} = {value_text}")
        elif ftype == "angle_3pt":
            sel_a = _ref_to_selection(spec.get("a", ""), ligand_site_selections, protein_role_selections)
            sel_v = _ref_to_selection(spec.get("vertex", ""), ligand_site_selections, protein_role_selections)
            sel_c = _ref_to_selection(spec.get("c", ""), ligand_site_selections, protein_role_selections)
            if sel_a and sel_v and sel_c:
                name = f"ANGLE_{_safe_name(feature_name)}"
                pml_lines.append(f"angle {name}, {sel_a}, {sel_v}, {sel_c}")
                if value_text:
                    pml_lines.append(f"# {feature_name} = {value_text}")
        elif ftype == "axis_deviation":
            if value_text:
                pml_lines.append(f"# {feature_name} axis_deviation = {value_text}")

    pml_lines.append("python")
    pml_lines.append("from pymol import cmd")
    pml_lines.append("for _sel in ['CANDIDATE_LIGAND_ATOMS','CANDIDATE_PROTEIN_ATOMS','CANDIDATE_ALL_ATOMS']:")
    pml_lines.append("    if cmd.count_atoms(_sel) > 0:")
    pml_lines.append("        print('selection_count', _sel, cmd.count_atoms(_sel))")
    pml_lines.append("print('protein_atoms', cmd.count_atoms('protein'))")
    pml_lines.append("print('ligand_atoms', cmd.count_atoms('ligand'))")
    pml_lines.append("python end")
    focus = "CANDIDATE_ALL_ATOMS" if individual_ligand_sels or individual_protein_sels else "all"
    pml_lines.append(f"zoom {focus}, 10")
    pml_lines.append("deselect")
    if save_pse:
        save_pse.parent.mkdir(parents=True, exist_ok=True)
        pml_lines.append("python")
        pml_lines.append("from pathlib import Path")
        pml_lines.append("from pymol import cmd")
        pml_lines.append(f"_pse_path = r'{str(save_pse)}'")
        pml_lines.append("Path(_pse_path).parent.mkdir(parents=True, exist_ok=True)")
        pml_lines.append("cmd.save(_pse_path)")
        pml_lines.append("python end")
    pml_lines.append("# Run this script with: pymol -cq <script.pml>")

    with open(out_pml, "w", encoding="utf-8") as handle:
        handle.write("\n".join(pml_lines) + "\n")



def export_pymol_batch_scripts(
    protein_file: Path,
    docked_ligand_file: Path,
    resolved_ligand_sites: Path,
    resolved_protein_roles: Path,
    candidate_scores: Path,
    out_dir: Path,
    passing_candidates: Path | None = None,
    geometry_features: Path | None = None,
    mechanism_file: Path | None = None,
    top_n: int | None = None,
    mode: str = "pose",
    make_pse: bool = True,
    pml_suffix: str = ".pml",
    pse_suffix: str = ".pse",
) -> list[dict[str, str]]:
    out_dir = _abs_path(out_dir) or Path("pymol_candidates").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    candidate_scores = _abs_path(candidate_scores)
    assert candidate_scores is not None
    _require_file(candidate_scores, "candidate scores")
    rows = _read_rows(candidate_scores)
    if not rows and passing_candidates is not None:
        rows = _read_rows(_abs_path(passing_candidates))
    ranked = _ranked_candidate_rows(rows, mode=mode, top_n=top_n)
    if not ranked:
        raise ValueError("No candidate rows were available for batch PyMOL export")

    report: list[dict[str, str]] = []
    for rank, row in ranked:
        stem = _candidate_file_stem(rank, row)
        out_pml = out_dir / f"{stem}{pml_suffix}"
        out_pse = out_dir / f"{stem}{pse_suffix}" if make_pse else None
        export_pymol_script(
            protein_file=protein_file,
            docked_ligand_file=docked_ligand_file,
            resolved_ligand_sites=resolved_ligand_sites,
            resolved_protein_roles=resolved_protein_roles,
            out_pml=out_pml,
            candidate_scores=candidate_scores,
            passing_candidates=passing_candidates,
            geometry_features=geometry_features,
            mechanism_file=mechanism_file,
            candidate_site_set_id=row.get("site_set_id") or None,
            pose_id=row.get("pose_id") or None,
            select_mode="first",
            save_pse_name=str(out_pse) if out_pse else None,
        )
        report.append({
            "rank": str(rank),
            "pose_id": str(row.get("pose_id", "")),
            "site_set_id": str(row.get("site_set_id", "")),
            "candidate_score": str(row.get("candidate_score", "")),
            "pml": str(out_pml),
            "pse": str(out_pse) if out_pse else "",
        })
    report_path = out_dir / "pymol_export_report.csv"
    with open(report_path, "w", encoding="utf-8", newline="") as handle:
        fieldnames = ["rank", "pose_id", "site_set_id", "candidate_score", "pml", "pse"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(report)
    return report
