from __future__ import annotations

import base64
import csv
import json
import math
import os
import random
from collections import Counter
from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple


def clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "null"}:
        return ""
    return text


def safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(clean(value)))
    except Exception:
        return default


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(clean(value))
    except Exception:
        return default


def natural_key(text: Any) -> Tuple[Any, ...]:
    value = clean(text)
    parts: List[Any] = []
    current = ""
    numeric = False
    for ch in value:
        if ch.isdigit():
            if current and not numeric:
                parts.append(current)
                current = ""
            numeric = True
            current += ch
        else:
            if current and numeric:
                parts.append(int(current))
                current = ""
            numeric = False
            current += ch
    if current:
        parts.append(int(current) if numeric else current)
    return tuple(parts)


def html_escape(text: Any) -> str:
    return escape(clean(text), quote=True)


def read_csv(path: Path) -> List[Dict[str, str]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, rows: List[Dict[str, Any]], skip_keys: Optional[Set[str]] = None) -> None:
    skip_keys = skip_keys or set()
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in skip_keys and key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            out: Dict[str, Any] = {}
            for key in keys:
                value = row.get(key, "")
                out[key] = json.dumps(value, ensure_ascii=True) if isinstance(value, (dict, list)) else value
            writer.writerow(out)


def first_existing(network_dir: Path, names: Sequence[str]) -> Path:
    for name in names:
        path = Path(network_dir) / name
        if path.exists():
            return path
    raise FileNotFoundError(f"None of these files were found in {network_dir}: {', '.join(names)}")


def image_to_data_uri(path: Optional[Path]) -> str:
    if not path:
        return ""
    path = Path(path)
    if not path.exists() or not path.is_file():
        return ""
    suffix = path.suffix.lower()
    if suffix == ".png":
        mime = "image/png"
    elif suffix in {".jpg", ".jpeg"}:
        mime = "image/jpeg"
    elif suffix == ".svg":
        mime = "image/svg+xml"
    else:
        mime = "application/octet-stream"
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{data}"


def make_mol_svg(smiles: str, width: int = 260, height: int = 180) -> str:
    smiles = clean(smiles)
    if not smiles:
        return ""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svg = svg.replace("<?xml version='1.0' encoding='iso-8859-1'?>", "")
        svg = svg.replace("svg:", "")
        return svg.strip()
    except Exception:
        return ""


def read_scaffold_json(path: Path, json_index_base: int) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError("The scaffold JSON must be a list of molecule records.")
    rows: List[Dict[str, Any]] = []
    smiles_index: Dict[str, Dict[str, Any]] = {}
    for idx, rec in enumerate(data):
        core = rec.get("core", {}) or {}
        chain = rec.get("chain", {}) or {}
        group = clean(core.get("abstract_core", "NA")) or "NA"
        row = {
            "json_index0": idx,
            "json_index1": idx + 1,
            "fallback_id": str(idx + json_index_base),
            "scaffold_group_raw": group,
            "scaffold_group": group,
            "scaffold_label": f"Scaffold_{group}",
            "core_smarts": clean(core.get("smarts", "")),
            "smiles_json": clean(rec.get("smiles", "")),
            "chain": {str(k): clean(v) for k, v in chain.items()},
        }
        rows.append(row)
        if row["smiles_json"] and row["smiles_json"] not in smiles_index:
            smiles_index[row["smiles_json"]] = row
    return rows, smiles_index


def load_valid_molecule_map(valid_molecules: Optional[Path], network_dir: Path) -> Dict[str, str]:
    candidates: List[Path] = []
    if valid_molecules:
        candidates.append(Path(valid_molecules))
    candidates.extend([
        network_dir.parent / "compute" / "transformnet.valid_molecules.csv",
        network_dir.parent / "taxane_template_linkage.valid_molecules.csv",
        network_dir / "taxane_template_linkage.valid_molecules.csv",
    ])
    for path in candidates:
        if not path.exists():
            continue
        out: Dict[str, str] = {}
        for row in read_csv(path):
            matrix = clean(row.get("matrix_index") or row.get("valid_index") or row.get("row_index"))
            source = clean(row.get("source_id") or row.get("source_no") or row.get("molecule_id"))
            if matrix and source:
                out[matrix] = source
        if out:
            return out
    return {}


def infer_node_id(row: Dict[str, str]) -> str:
    for key in ("source_id", "molecule_id", "node_id", "id", "source_no"):
        value = clean(row.get(key))
        if value:
            return value
    value = clean(row.get("row_index"))
    if value:
        return value
    raise ValueError("Node table must contain source_id, molecule_id, node_id, id, source_no, or row_index.")


def infer_edge_ids(row: Dict[str, str]) -> Tuple[str, str, str]:
    pairs = [
        ("node1_id", "node2_id", "node_id"),
        ("source_id", "target_id", "source_target"),
        ("substrate_id", "product_id", "substrate_product"),
        ("source", "target", "source_target"),
        ("node1", "node2", "node"),
    ]
    for a, b, mode in pairs:
        av = clean(row.get(a))
        bv = clean(row.get(b))
        if av and bv:
            return av, bv, mode
    raise ValueError("Edge table must contain supported endpoint columns.")


def choose_name(row: Dict[str, str]) -> str:
    for key in ("molecule_name", "compound", "name", "input_molecule_name", "node_name"):
        value = clean(row.get(key))
        if value:
            return value
    return ""


def choose_smiles(row: Dict[str, str]) -> str:
    for key in ("smiles", "smiles_raw", "standardized_smiles", "canonical_smiles"):
        value = clean(row.get(key))
        if value:
            return value
    return ""


def attach_scaffold(node_row: Dict[str, str], scaffold_rows: List[Dict[str, Any]], smiles_index: Dict[str, Dict[str, Any]], json_index_base: int) -> Tuple[Dict[str, Any], Dict[str, str]]:
    existing_chain = {k: clean(v) for k, v in node_row.items() if k.startswith("R") and k[1:].isdigit() and clean(v)}
    if clean(node_row.get("scaffold_group")) and existing_chain:
        group = clean(node_row.get("scaffold_group"))
        info = {
            "scaffold_group_raw": group,
            "scaffold_group": group,
            "scaffold_label": clean(node_row.get("scaffold_label") or f"Scaffold_{group}"),
            "core_smarts": clean(node_row.get("core_smarts")),
            "smiles_json": clean(node_row.get("scaffold_smiles") or choose_smiles(node_row)),
            "chain": existing_chain,
            "json_index0": clean(node_row.get("scaffold_json_index")),
            "json_index1": "",
        }
        return info, {"mapping_status": "node_table", "mapping_reason": "node_table_has_scaffold_and_chain"}
    smi = choose_smiles(node_row)
    if smi in smiles_index:
        return smiles_index[smi], {"mapping_status": "exact_smiles", "mapping_reason": "node_smiles_matches_json_smiles"}
    for key in ("source_id", "source_no", "molecule_id"):
        value = clean(node_row.get(key))
        if not value:
            continue
        idx = safe_int(value, -10**9) - json_index_base
        if 0 <= idx < len(scaffold_rows):
            return scaffold_rows[idx], {"mapping_status": "id_fallback", "mapping_reason": f"{key}_with_base_{json_index_base}"}
    for key in ("row_index", "matrix_index"):
        value = clean(node_row.get(key))
        if not value:
            continue
        idx = safe_int(value, -10**9)
        if 0 <= idx < len(scaffold_rows):
            return scaffold_rows[idx], {"mapping_status": "row_index_fallback", "mapping_reason": key}
    empty = {"scaffold_group_raw": "NA", "scaffold_group": "NA", "scaffold_label": "Scaffold_NA", "core_smarts": "", "smiles_json": "", "chain": {}, "json_index0": "", "json_index1": ""}
    return empty, {"mapping_status": "unmapped", "mapping_reason": "no_supported_mapping"}


def compare_chains(source: Dict[str, Any], target: Dict[str, Any]) -> Tuple[str, int, str, str, List[Dict[str, str]]]:
    s_chain = source.get("chain", {}) or {}
    t_chain = target.get("chain", {}) or {}
    positions = sorted(set(s_chain) | set(t_chain), key=natural_key)
    changed: List[str] = []
    detail: Dict[str, Dict[str, str]] = {}
    text: List[str] = []
    long_rows: List[Dict[str, str]] = []
    for pos in positions:
        sv = clean(s_chain.get(pos))
        tv = clean(t_chain.get(pos))
        if sv != tv:
            changed.append(pos)
            detail[pos] = {"source": sv, "target": tv}
            text.append(f"{pos}:{sv}->{tv}")
            long_rows.append({"position": pos, "source_sidechain": sv, "target_sidechain": tv})
    return ";".join(changed), len(changed), json.dumps(detail, ensure_ascii=True), " | ".join(text), long_rows


def build_node_maps(raw_nodes: List[Dict[str, str]], label_column: str) -> Tuple[Dict[str, Dict[str, str]], Dict[str, str], Dict[str, str], Dict[str, str], Dict[str, str]]:
    node_id_to_raw: Dict[str, Dict[str, str]] = {}
    row_index_to_node_id: Dict[str, str] = {}
    source_to_node_id: Dict[str, str] = {}
    label_to_node_id: Dict[str, str] = {}
    matrix_to_node_id: Dict[str, str] = {}
    for idx, row in enumerate(raw_nodes):
        nid = infer_node_id(row)
        row["__node_id"] = nid
        node_id_to_raw[nid] = row
        row_index_to_node_id[str(idx)] = nid
        row_index_to_node_id[str(idx + 1)] = nid
        for key in ("row_index", "matrix_index"):
            value = clean(row.get(key))
            if value:
                row_index_to_node_id[value] = nid
                if key == "matrix_index":
                    matrix_to_node_id[value] = nid
        for key in ("source_id", "source_no", "molecule_id"):
            value = clean(row.get(key))
            if value and value not in source_to_node_id:
                source_to_node_id[value] = nid
        for key in (label_column, "label", "node_label", "name"):
            value = clean(row.get(key))
            if value and value not in label_to_node_id:
                label_to_node_id[value] = nid
    return node_id_to_raw, row_index_to_node_id, source_to_node_id, label_to_node_id, matrix_to_node_id


def resolve_value(value: str, mode: str, node_ids: Set[str], row_index_to_node_id: Dict[str, str], matrix_to_source: Dict[str, str], source_to_node_id: Dict[str, str], label_to_node_id: Dict[str, str], matrix_to_node_id: Dict[str, str]) -> str:
    value = clean(value)
    if not value:
        return ""
    if mode == "node_id":
        return value if value in node_ids else ""
    if mode == "source_id":
        return source_to_node_id.get(value, value if value in node_ids else "")
    if mode == "row_index":
        return row_index_to_node_id.get(value, "")
    if mode == "label":
        return label_to_node_id.get(value, "")
    if mode == "matrix_index":
        if value in matrix_to_node_id:
            return matrix_to_node_id[value]
        mapped_source = matrix_to_source.get(value)
        if mapped_source:
            return source_to_node_id.get(mapped_source, mapped_source if mapped_source in node_ids else "")
        return row_index_to_node_id.get(value, "")
    return ""


def choose_edge_id_space(raw_edges: List[Dict[str, str]], requested: str, node_ids: Set[str], row_index_to_node_id: Dict[str, str], matrix_to_source: Dict[str, str], source_to_node_id: Dict[str, str], label_to_node_id: Dict[str, str], matrix_to_node_id: Dict[str, str]) -> Tuple[str, Dict[str, int]]:
    modes = ["source_id", "node_id", "matrix_index", "row_index", "label"] if requested == "auto" else [requested]
    scores: Dict[str, int] = {}
    sample = raw_edges[: min(1000, len(raw_edges))]
    for mode in modes:
        ok = 0
        for row in sample:
            try:
                a, b, _ = infer_edge_ids(row)
            except Exception:
                continue
            ra = resolve_value(a, mode, node_ids, row_index_to_node_id, matrix_to_source, source_to_node_id, label_to_node_id, matrix_to_node_id)
            rb = resolve_value(b, mode, node_ids, row_index_to_node_id, matrix_to_source, source_to_node_id, label_to_node_id, matrix_to_node_id)
            if ra in node_ids and rb in node_ids:
                ok += 1
        scores[mode] = ok
    if requested != "auto":
        return requested, scores
    best = max(scores.items(), key=lambda x: (x[1], x[0]))[0] if scores else "source_id"
    return best, scores


def get_template_count(row: Dict[str, str]) -> int:
    for key in ("template_count_total", "template_count", "event_count_total", "edge_weight", "weight"):
        value = clean(row.get(key))
        if value:
            return max(1, safe_int(value, 1))
    return 1


def get_event_count(row: Dict[str, str], default: int) -> int:
    for key in ("event_count_total", "event_count", "template_count_total", "template_count"):
        value = clean(row.get(key))
        if value:
            return max(1, safe_int(value, default))
    return default


def layer_match_factor(layer: str) -> float:
    if layer == "strict":
        return 1.0
    if layer == "relaxed":
        return 0.6
    if layer == "original":
        return 0.5
    return 1.0


def high_contrast_palette(index: int) -> str:
    colors = [
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
        "#FFFF33", "#A65628", "#F781BF", "#1B9E77", "#D95F02",
        "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
        "#666666", "#8DD3C7", "#FB8072", "#80B1D3", "#B3DE69",
    ]
    return colors[index % len(colors)]


def muted_palette(index: int) -> str:
    colors = ["#9ECAE1", "#A1D99B", "#BCBDDC", "#FDD0A2", "#C7E9C0", "#FCAE91", "#C6DBEF", "#D9D9D9", "#FDD49E", "#BDBDBD"]
    return colors[index % len(colors)]


def build_colors_by_group(group_counts: Counter, central_groups: int) -> Tuple[Dict[str, str], Set[str]]:
    ordered = [g for g, _ in sorted(group_counts.items(), key=lambda x: (-x[1], natural_key(x[0])))]
    central = set(ordered[: max(0, central_groups)])
    colors: Dict[str, str] = {}
    major_i = 0
    minor_i = 0
    for g in ordered:
        if g in central:
            colors[g] = high_contrast_palette(major_i)
            major_i += 1
        else:
            colors[g] = muted_palette(minor_i)
            minor_i += 1
    colors.setdefault("NA", "#999999")
    return colors, central


def compute_circular_group_layout(nodes: List[Dict[str, Any]], edges: List[Dict[str, Any]], layout_seed: int, central_groups: int) -> None:
    rng = random.Random(layout_seed)
    width, height = 1450.0, 980.0
    cx, cy = width / 2.0, height / 2.0
    group_counts = Counter(clean(n.get("scaffold_group") or "NA") for n in nodes)
    ordered = [g for g, _ in sorted(group_counts.items(), key=lambda x: (-x[1], natural_key(x[0])))]
    central = ordered[: max(0, central_groups)]
    outer = ordered[max(0, central_groups):]
    centers: Dict[str, Tuple[float, float]] = {}
    if central:
        if len(central) == 1:
            centers[central[0]] = (cx, cy)
        else:
            inner_r = min(230.0, 105.0 + 18.0 * len(central))
            for i, g in enumerate(central):
                angle = -math.pi / 2 + 2 * math.pi * i / len(central)
                centers[g] = (cx + inner_r * math.cos(angle), cy + inner_r * math.sin(angle))
    if outer:
        outer_r = 410.0
        for i, g in enumerate(outer):
            angle = -math.pi / 2 + 2 * math.pi * i / len(outer)
            centers[g] = (cx + outer_r * math.cos(angle), cy + outer_r * math.sin(angle))
    for g in ordered:
        centers.setdefault(g, (cx, cy))
    group_nodes: Dict[str, List[Dict[str, Any]]] = {}
    for n in nodes:
        group_nodes.setdefault(clean(n.get("scaffold_group") or "NA"), []).append(n)
    central_set = set(central)
    for g, ns in group_nodes.items():
        gcx, gcy = centers.get(g, (cx, cy))
        count = len(ns)
        local_r = min(250.0 if g in central_set else 160.0, 34.0 + 11.5 * math.sqrt(max(count, 1)))
        ordered_nodes = sorted(ns, key=lambda x: natural_key(x.get("label") or x.get("id")))
        for i, n in enumerate(ordered_nodes):
            ring = 0 if count <= 18 else int(i / max(1, math.ceil(math.sqrt(count))))
            ring_factor = 0.42 + 0.22 * ring
            rr = min(local_r, local_r * ring_factor + rng.uniform(-12, 12))
            angle = 2 * math.pi * i / max(1, count) + rng.uniform(-0.20, 0.20)
            n["x"] = gcx + rr * math.cos(angle)
            n["y"] = gcy + rr * math.sin(angle)
            n["layout_group_center_x"] = gcx
            n["layout_group_center_y"] = gcy
            n["major_scaffold_group"] = g in central_set
    id_to_node = {str(n.get("id")): n for n in nodes}
    edge_pairs = [(id_to_node.get(str(e.get("source"))), id_to_node.get(str(e.get("target")))) for e in edges[:8000]]
    edge_pairs = [(a, b) for a, b in edge_pairs if a is not None and b is not None]
    for _ in range(95):
        disp: Dict[str, List[float]] = {str(n.get("id")): [0.0, 0.0] for n in nodes}
        for i in range(len(nodes)):
            a = nodes[i]
            aid = str(a.get("id"))
            for j in range(i + 1, len(nodes)):
                b = nodes[j]
                bid = str(b.get("id"))
                vx = float(a["x"]) - float(b["x"])
                vy = float(a["y"]) - float(b["y"])
                d2 = vx * vx + vy * vy + 130.0
                d = math.sqrt(d2)
                f = 1450.0 / d2
                fx = vx / d * f
                fy = vy / d * f
                disp[aid][0] += fx
                disp[aid][1] += fy
                disp[bid][0] -= fx
                disp[bid][1] -= fy
        for a, b in edge_pairs:
            aid = str(a.get("id")); bid = str(b.get("id"))
            vx = float(b["x"]) - float(a["x"])
            vy = float(b["y"]) - float(a["y"])
            k = 0.008
            disp[aid][0] += vx * k
            disp[aid][1] += vy * k
            disp[bid][0] -= vx * k
            disp[bid][1] -= vy * k
        for n in nodes:
            nid = str(n.get("id"))
            gx = float(n.get("layout_group_center_x", cx)) - float(n["x"])
            gy = float(n.get("layout_group_center_y", cy)) - float(n["y"])
            d = disp[nid]
            d[0] += gx * 0.0045
            d[1] += gy * 0.0045
            n["x"] += max(-8.0, min(8.0, d[0]))
            n["y"] += max(-8.0, min(8.0, d[1]))


def edge_table_path(network_dir: Path) -> Path:
    return first_existing(network_dir, ["transformnet_network.collapsed_edges.csv", "transformnet_network.edges.csv"])


def nodes_table_path(network_dir: Path) -> Path:
    return first_existing(network_dir, ["transformnet_network.nodes.csv"])


def build_visual_edges(
    raw_edges: List[Dict[str, str]],
    layer_name: str,
    node_id_to_raw: Dict[str, Dict[str, str]],
    row_index_to_node_id: Dict[str, str],
    source_to_node_id: Dict[str, str],
    label_to_node_id: Dict[str, str],
    matrix_to_node_id: Dict[str, str],
    matrix_to_source: Dict[str, str],
    selected_edge_id_space: str,
    included_ids: Set[str],
    node_info: Dict[str, Dict[str, Any]],
    max_edges: int,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    node_ids = set(node_id_to_raw)
    visual_edges: List[Dict[str, Any]] = []
    sidechain_rows: List[Dict[str, Any]] = []
    unresolved_edges: List[Dict[str, Any]] = []
    for raw in raw_edges:
        raw_a, raw_b, _ = infer_edge_ids(raw)
        source_id = resolve_value(raw_a, selected_edge_id_space, node_ids, row_index_to_node_id, matrix_to_source, source_to_node_id, label_to_node_id, matrix_to_node_id)
        target_id = resolve_value(raw_b, selected_edge_id_space, node_ids, row_index_to_node_id, matrix_to_source, source_to_node_id, label_to_node_id, matrix_to_node_id)
        if source_id not in included_ids or target_id not in included_ids:
            unresolved_edges.append({**raw, "raw_edge_node1": raw_a, "raw_edge_node2": raw_b, "resolved_source": source_id, "resolved_target": target_id, "edge_id_space": selected_edge_id_space, "edge_layer": layer_name})
            continue
        if len(visual_edges) >= max_edges:
            break
        s = node_info[source_id]
        t = node_info[target_id]
        changed, changed_count, detail, text, longs = compare_chains(s, t)
        template_count = get_template_count(raw)
        event_count = get_event_count(raw, template_count)
        same_group = clean(s.get("scaffold_group")) == clean(t.get("scaffold_group"))
        match_factor = layer_match_factor(layer_name)
        structural_simplicity = 1.0 / max(1, changed_count)
        tss_raw = math.log1p(template_count) * match_factor * structural_simplicity
        edge = {
            "source": source_id,
            "target": target_id,
            "source_label": clean(s.get("label") or s.get("source_id") or source_id),
            "target_label": clean(t.get("label") or t.get("source_id") or target_id),
            "node1": source_id,
            "node2": target_id,
            "raw_edge_node1": raw_a,
            "raw_edge_node2": raw_b,
            "directionality": clean(raw.get("directionality") or f"{source_id}->{target_id}"),
            "template_count_total": template_count,
            "event_count_total": event_count,
            "source_scaffold_group": s.get("scaffold_group"),
            "target_scaffold_group": t.get("scaffold_group"),
            "same_group": same_group,
            "changed_positions": changed,
            "changed_position_count": changed_count,
            "sidechain_change_detail": detail,
            "sidechain_change_text": text,
            "edge_id_mode": selected_edge_id_space,
            "edge_layer": layer_name,
            "match_factor": round(match_factor, 3),
            "structural_simplicity": round(structural_simplicity, 3),
            "transformation_support_score_raw": round(tss_raw, 6),
        }
        visual_edges.append(edge)
        for lr in longs:
            sidechain_rows.append({"source_id": edge["source_label"], "target_id": edge["target_label"], "position": lr["position"], "source_sidechain": lr["source_sidechain"], "target_sidechain": lr["target_sidechain"], "source_scaffold_group": s.get("scaffold_group"), "target_scaffold_group": t.get("scaffold_group"), "same_group": same_group, "template_count_total": template_count, "edge_layer": layer_name})
    return visual_edges, sidechain_rows, unresolved_edges


def normalize_support_scores(edge_layers: Dict[str, List[Dict[str, Any]]]) -> None:
    all_edges = [e for layer in edge_layers.values() for e in layer]
    max_raw = max([safe_float(e.get("transformation_support_score_raw"), 0.0) for e in all_edges] or [0.0])
    for edge in all_edges:
        raw = safe_float(edge.get("transformation_support_score_raw"), 0.0)
        edge["transformation_support_score"] = round(raw / max_raw, 4) if max_raw > 0 else 0.0
        if edge["transformation_support_score"] >= 0.8:
            edge["support_bin"] = "high"
        elif edge["transformation_support_score"] >= 0.5:
            edge["support_bin"] = "medium"
        else:
            edge["support_bin"] = "low"


def build_html(nodes: List[Dict[str, Any]], edge_layers: Dict[str, List[Dict[str, Any]]], title: str, summary: Dict[str, Any], edge_scores: Dict[str, int], scaffold_image_uri: str = "", label_size: int = 10) -> str:
    payload = json.dumps({"nodes": nodes, "edge_layers": edge_layers, "summary": summary, "edge_scores": edge_scores, "scaffold_image_uri": scaffold_image_uri}, ensure_ascii=True).replace("</", "<\\/")
    safe_title = html_escape(title)
    html = r'''<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
html, body { margin: 0; height: 100%; font-family: Arial, sans-serif; }
#app { display: grid; grid-template-columns: minmax(620px, 1fr) 440px; height: 100vh; }
#stage { position: relative; background: #fafafa; overflow: hidden; display: grid; grid-template-rows: auto 1fr; }
#stageHeader { display: flex; align-items: flex-start; gap: 12px; padding: 8px 10px; border-bottom: 1px solid #ddd; background: #fff; min-height: 56px; max-height: 165px; overflow: auto; }
#scaffoldImageWrap { flex: 0 0 auto; }
#scaffoldImageWrap img { max-height: 148px; max-width: 520px; display: block; }
#topLegend { display: flex; flex-wrap: wrap; gap: 4px 10px; align-items: flex-start; font-size: 11px; line-height: 1.35; padding-top: 2px; }
#net { width: 100%; height: 100%; display: block; background: #fafafa; }
aside { border-left: 1px solid #ddd; padding: 12px; overflow: auto; background: #fff; }
h2 { margin: 0 0 8px 0; font-size: 18px; }
h3 { margin: 14px 0 6px 0; font-size: 14px; }
.small { color: #555; font-size: 12px; line-height: 1.4; }
button { margin: 4px 4px 4px 0; padding: 5px 8px; font-family: Arial, sans-serif; }
input, select { width: 100%; box-sizing: border-box; margin-top: 4px; padding: 5px; font-family: Arial, sans-serif; }
pre { white-space: pre-wrap; word-break: break-word; background: #f6f6f6; border: 1px solid #ddd; padding: 8px; border-radius: 4px; font-size: 11px; }
.result-item { padding: 4px; border-bottom: 1px solid #eee; cursor: pointer; font-size: 12px; }
.result-item:hover { background: #f2f6ff; }
.swatch { display: inline-block; width: 12px; height: 12px; border-radius: 3px; margin-right: 5px; vertical-align: middle; }
.row { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; }
.edge { stroke-linecap: round; cursor: pointer; }
.node { cursor: move; }
.node.selected { }
.node.path { }
.edge.path { opacity: 0.95 !important; }
.edge.highlight { opacity: 0.95 !important; }
.node.highlight { }
.label { font-size: __LABEL_SIZE__px; fill: #111; pointer-events: none; text-anchor: middle; dominant-baseline: central; font-family: Arial, sans-serif; paint-order: stroke; stroke: #fff; stroke-width: 2px; stroke-linejoin: round; }
.detail-card { border: 1px solid #ddd; border-radius: 6px; padding: 8px; margin: 8px 0; }
.detail-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; }
.mol-svg svg { width: 100%; height: auto; max-height: 190px; }
.kv { font-size: 12px; line-height: 1.5; }
.kv b { display: inline-block; min-width: 112px; }
</style>
</head>
<body>
<div id="app">
<div id="stage"><div id="stageHeader"><div id="scaffoldImageWrap"></div><div id="topLegend"></div></div><svg id="net" xmlns="http://www.w3.org/2000/svg"></svg></div>
<aside>
<h2>__TITLE__</h2>
<h3>Edge layer</h3>
<select id="edgeLayer"></select>
<h3>Explorer controls</h3>
<label>Explorer mode</label><select id="explorerMode"><option value="standard">Standard filter</option><option value="transform">Transformation Explorer</option></select>
<div class="row"><div><label>Direction</label><select id="directionMode"><option value="downstream">downstream</option><option value="upstream">upstream</option><option value="both">both</option></select></div><div><label>Hop depth</label><select id="hopDepth"><option value="1">1-hop</option><option value="2">2-hop</option><option value="3">3-hop</option></select></div></div>
<label>Display mode</label><select id="displayMode"><option value="filter">filter</option><option value="highlight">highlight</option></select>
<div class="row"><div><label>Scaffold group</label><select id="scaffold"><option value="all">all</option></select></div><div><label>Changed position</label><select id="position"><option value="all">all</option></select></div></div>
<label>Starting or matched side-chain</label><select id="sidechain"><option value="all">all</option></select>
<label><input id="cross" type="checkbox"> Show cross-scaffold edges only</label>
<h3>Search nodes</h3>
<input id="search" placeholder="Multiple exact IDs, separated by comma, space, or newline">
<label><input id="fuzzySearch" type="checkbox"> Enable fuzzy search when no exact ID is found</label>
<button id="applySearch">Search</button><button id="clearSearch">Clear</button><button id="fitSelected">Fit selected</button>
<h3>Focus mode</h3>
<select id="focusMode"><option value="all">All graph</option><option value="selected">Selected nodes only</option><option value="onehop">Selected plus one-hop neighbors</option><option value="paths">Shortest paths among selected nodes</option></select>
<label><input id="showLabels" type="checkbox" checked> Show node labels</label>
<button id="relax">Relax visible layout</button><button id="fit">Fit graph</button><button id="reset">Reset filters</button>
<h3>Search results</h3>
<div id="results"></div>
<h3>Details</h3>
<div id="details" class="small">Click a node or edge.</div>
</aside>
</div>
<script>
const GRAPH = __DATA__;
const NODE_LABEL_SIZE = __LABEL_SIZE__;
const svgNS = 'http://www.w3.org/2000/svg';
const svg = document.getElementById('net');
const details = document.getElementById('details');
let nodes = (GRAPH.nodes || []).map(x => Object.assign({}, x));
let edgeLayers = GRAPH.edge_layers || {};
let currentLayer = Object.keys(edgeLayers)[0] || 'current';
let edges = [];
let byId = new Map(nodes.map(n => [String(n.id), n]));
let selectedIds = new Set();
let pathNodeIds = new Set();
let pathEdgeKeys = new Set();
let highlightedEdgeKeys = new Set();
let highlightedNodeIds = new Set();
let currentViewBox = null;
let isPanning = false;
let panStart = null;
let draggedNode = null;
let edgeEls = [];
let nodeEls = [];
function esc(x) { return String(x == null ? '' : x); }
function html(x) { return esc(x).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c])); }
function edgeKey(a,b) { return [String(a), String(b)].sort().join('||'); }
function directedEdgeKey(e) { return String(e.source) + '->' + String(e.target) + '::' + String(e.edge_layer || currentLayer); }
function sourceLabel(n) { return String(n.label || n.source_id || n.id); }
function parseChanges(e) {
  if (Array.isArray(e.parsed_changes)) return e.parsed_changes;
  const out = [];
  const txt = String(e.sidechain_change_text || '');
  txt.split('|').forEach(part => {
    const p = part.trim();
    if (!p) return;
    const idx = p.indexOf(':');
    const arrow = p.indexOf('->');
    if (idx > 0 && arrow > idx) out.push({position: p.slice(0, idx).trim(), from: p.slice(idx + 1, arrow).trim(), to: p.slice(arrow + 2).trim()});
  });
  e.parsed_changes = out;
  return out;
}
function allLayerNames() { return Object.keys(edgeLayers); }
function mergeOverlay() {
  const names = allLayerNames();
  if (names.length <= 1) return (edgeLayers[names[0]] || []).map(e => Object.assign({}, e));
  const map = new Map();
  names.forEach(name => {
    (edgeLayers[name] || []).forEach(e => {
      const key = String(e.source) + '||' + String(e.target);
      if (!map.has(key)) {
        const x = Object.assign({}, e);
        x.edge_layer = name;
        x.edge_layer_members = [name];
        map.set(key, x);
      } else {
        const x = map.get(key);
        if (!x.edge_layer_members.includes(name)) x.edge_layer_members.push(name);
        x.edge_layer = x.edge_layer_members.length > 1 ? 'both' : x.edge_layer_members[0];
        x.template_count_total = Math.max(Number(x.template_count_total || 0), Number(e.template_count_total || 0));
        x.transformation_support_score = Math.max(Number(x.transformation_support_score || 0), Number(e.transformation_support_score || 0));
      }
    });
  });
  return Array.from(map.values());
}
function setCurrentEdges() {
  const selected = document.getElementById('edgeLayer') ? document.getElementById('edgeLayer').value : currentLayer;
  currentLayer = selected || currentLayer;
  edges = currentLayer === 'overlay' ? mergeOverlay() : (edgeLayers[currentLayer] || []).map(e => Object.assign({}, e));
  edges.forEach(parseChanges);
  updatePositionOptions();
  updateSidechainOptions();
}
function exactFields(n) { return [n.source_id,n.label,n.id,n.molecule_id].map(x => esc(x).toLowerCase()).filter(Boolean); }
function fuzzyFields(n) { return [n.id,n.label,n.name,n.scaffold_group,n.scaffold_label,n.source_id,n.molecule_id,n.smiles].map(x => esc(x).toLowerCase()).join(' '); }
function queryTokens() { return document.getElementById('search').value.split(/[\n,;\s]+/).map(x => x.trim().toLowerCase()).filter(Boolean); }
function searchNodes() {
  const toks = queryTokens();
  if (!toks.length) return [];
  const exact = [];
  const seen = new Set();
  for (const t of toks) {
    for (const n of nodes) {
      if (exactFields(n).includes(t) && !seen.has(String(n.id))) { exact.push(n); seen.add(String(n.id)); }
    }
  }
  if (exact.length || !document.getElementById('fuzzySearch').checked) return exact;
  const hits = [];
  for (const n of nodes) {
    const text = fuzzyFields(n);
    if (toks.some(t => text.includes(t))) hits.push(n);
  }
  return hits;
}
function changeMatches(e, mode, pos, chain, direction) {
  const changes = parseChanges(e);
  if (!changes.length && (pos !== 'all' || chain !== 'all')) return false;
  return changes.some(c => {
    if (pos !== 'all' && c.position !== pos) return false;
    if (chain === 'all') return true;
    if (mode === 'transform') {
      if (direction === 'downstream') return c.from === chain;
      if (direction === 'upstream') return c.to === chain;
      return c.from === chain || c.to === chain;
    }
    return c.from === chain || c.to === chain;
  });
}
function scaffoldPass(e) {
  const sg = document.getElementById('scaffold').value;
  if (sg === 'all') return true;
  const s = byId.get(String(e.source));
  const t = byId.get(String(e.target));
  return (s && String(s.scaffold_group) === sg) || (t && String(t.scaffold_group) === sg);
}
function baseEdgePass(e, ignoreChange=false) {
  if (!scaffoldPass(e)) return false;
  if (document.getElementById('cross').checked && e.same_group) return false;
  if (ignoreChange) return true;
  const mode = document.getElementById('explorerMode').value;
  const pos = document.getElementById('position').value;
  const chain = document.getElementById('sidechain').value;
  const direction = document.getElementById('directionMode').value;
  if (pos === 'all' && chain === 'all') return true;
  return changeMatches(e, mode, pos, chain, direction);
}
function seedEdgesForExplorer() {
  const pos = document.getElementById('position').value;
  const chain = document.getElementById('sidechain').value;
  const direction = document.getElementById('directionMode').value;
  return edges.filter(e => scaffoldPass(e) && !((document.getElementById('cross').checked) && e.same_group) && changeMatches(e, 'transform', pos, chain, direction));
}
function explorerEdgeKeys() {
  if (document.getElementById('explorerMode').value !== 'transform') return new Set(edges.filter(e => baseEdgePass(e)).map(directedEdgeKey));
  const direction = document.getElementById('directionMode').value;
  const depth = Number(document.getElementById('hopDepth').value || 1);
  const result = new Set();
  const visitedNodes = new Set();
  let frontier = new Set();
  const seeds = seedEdgesForExplorer();
  seeds.forEach(e => {
    result.add(directedEdgeKey(e));
    visitedNodes.add(String(e.source)); visitedNodes.add(String(e.target));
    if (direction === 'downstream') frontier.add(String(e.target));
    else if (direction === 'upstream') frontier.add(String(e.source));
    else { frontier.add(String(e.source)); frontier.add(String(e.target)); }
  });
  for (let level=1; level<depth; level++) {
    const next = new Set();
    edges.forEach(e => {
      if (!baseEdgePass(e, true)) return;
      const s=String(e.source), t=String(e.target);
      let take=false, addNode='';
      if (direction === 'downstream' && frontier.has(s)) { take=true; addNode=t; }
      else if (direction === 'upstream' && frontier.has(t)) { take=true; addNode=s; }
      else if (direction === 'both' && (frontier.has(s) || frontier.has(t))) { take=true; addNode=frontier.has(s) ? t : s; }
      if (take) { result.add(directedEdgeKey(e)); if (!visitedNodes.has(addNode)) { visitedNodes.add(addNode); next.add(addNode); } }
    });
    frontier = next;
    if (!frontier.size) break;
  }
  return result;
}
function matchedEdgeKeys() { return explorerEdgeKeys(); }
function baseFilteredEdges(applyFocus, matchedOnly) {
  let es = edges.filter(e => baseEdgePass(e, document.getElementById('explorerMode').value === 'transform'));
  const display = document.getElementById('displayMode').value;
  if (matchedOnly || display === 'filter') {
    const keys = matchedEdgeKeys();
    es = es.filter(e => keys.has(directedEdgeKey(e)));
  }
  if (!applyFocus) return es;
  const mode = document.getElementById('focusMode').value;
  if (!selectedIds.size || mode === 'all') return es;
  if (mode === 'selected') es = es.filter(e => selectedIds.has(String(e.source)) && selectedIds.has(String(e.target)));
  else if (mode === 'onehop') es = es.filter(e => selectedIds.has(String(e.source)) || selectedIds.has(String(e.target)));
  else if (mode === 'paths') es = es.filter(e => pathEdgeKeys.has(edgeKey(e.source, e.target)));
  return es;
}
function visible() {
  const display = document.getElementById('displayMode').value;
  const matched = matchedEdgeKeys();
  highlightedEdgeKeys = matched;
  highlightedNodeIds = new Set();
  edges.forEach(e => { if (matched.has(directedEdgeKey(e))) { highlightedNodeIds.add(String(e.source)); highlightedNodeIds.add(String(e.target)); } });
  const es = baseFilteredEdges(true, display === 'filter');
  const used = new Set();
  es.forEach(e => { used.add(String(e.source)); used.add(String(e.target)); });
  const sg = document.getElementById('scaffold').value;
  let ns = nodes.filter(n => sg === 'all' || String(n.scaffold_group) === sg || (display === 'highlight'));
  if (display === 'filter') ns = ns.filter(n => used.has(String(n.id)) || (es.length === 0 && (sg === 'all' || String(n.scaffold_group) === sg)));
  const mode = document.getElementById('focusMode').value;
  if (selectedIds.size && mode === 'selected') ns = ns.filter(n => selectedIds.has(String(n.id)));
  else if (selectedIds.size && mode === 'onehop') ns = ns.filter(n => selectedIds.has(String(n.id)) || used.has(String(n.id)));
  else if (selectedIds.size && mode === 'paths') ns = ns.filter(n => pathNodeIds.has(String(n.id)));
  return {nodes: ns, edges: es};
}
function buildAdj(activeEdges) { const adj = new Map(); function add(a,b){ if(!adj.has(a)) adj.set(a,[]); adj.get(a).push(b); } activeEdges.forEach(e=>{const a=String(e.source),b=String(e.target);add(a,b);add(b,a);}); return adj; }
function bfsPath(start, goal, adj) { start=String(start); goal=String(goal); if(start===goal)return[start]; const q=[start], prev=new Map([[start,null]]); for(let qi=0; qi<q.length; qi++){const u=q[qi]; for(const v of(adj.get(u)||[])){ if(prev.has(v))continue; prev.set(v,u); if(v===goal){const path=[v]; let cur=u; while(cur!==null){path.push(cur); cur=prev.get(cur);} return path.reverse();} q.push(v);}} return []; }
function updatePaths(){ pathNodeIds.clear(); pathEdgeKeys.clear(); const ids=Array.from(selectedIds); if(ids.length<2)return; const base=baseFilteredEdges(false, false); const adj=buildAdj(base); for(let i=0;i<ids.length;i++)for(let j=i+1;j<ids.length;j++){const p=bfsPath(ids[i],ids[j],adj); for(let k=0;k<p.length;k++)pathNodeIds.add(p[k]); for(let k=0;k<p.length-1;k++)pathEdgeKeys.add(edgeKey(p[k],p[k+1]));}}
function boundsFor(ns){ if(!ns.length)return null; const xs=ns.map(n=>Number(n.x||0)); const ys=ns.map(n=>Number(n.y||0)); return [Math.min(...xs),Math.max(...xs),Math.min(...ys),Math.max(...ys)]; }
function setViewBoxFromBounds(b){ if(!b)return; const [minx,maxx,miny,maxy]=b; const pad=45; const w=Math.max(10,maxx-minx+pad*2); const h=Math.max(10,maxy-miny+pad*2); currentViewBox=[minx-pad,miny-pad,w,h]; svg.setAttribute('viewBox', currentViewBox.join(' ')); }
function fitGraph(){ setViewBoxFromBounds(boundsFor(visible().nodes)); }
function fitSelected(){ const ns=nodes.filter(n=>selectedIds.has(String(n.id))||pathNodeIds.has(String(n.id))||highlightedNodeIds.has(String(n.id))); setViewBoxFromBounds(boundsFor(ns.length?ns:visible().nodes)); }
function centerOnSelectedNode(){
  const ns=Array.from(selectedIds).map(id=>byId.get(String(id))).filter(Boolean);
  if(!ns.length){ fitGraph(); return; }
  if(ns.length>5){ fitSelected(); return; }
  const xs=ns.map(n=>Number(n.x||0));
  const ys=ns.map(n=>Number(n.y||0));
  const cx=xs.reduce((a,b)=>a+b,0)/xs.length;
  const cy=ys.reduce((a,b)=>a+b,0)/ys.length;
  const minx=Math.min(...xs), maxx=Math.max(...xs), miny=Math.min(...ys), maxy=Math.max(...ys);
  const spreadX=Math.max(110,maxx-minx);
  const spreadY=Math.max(95,maxy-miny);
  const w=ns.length===1?280:Math.max(320,spreadX*2.8);
  const h=ns.length===1?230:Math.max(260,spreadY*2.8);
  currentViewBox=[cx-w/2,cy-h/2,w,h];
  svg.setAttribute('viewBox',currentViewBox.join(' '));
}
function svgPoint(evt){ if(!currentViewBox)fitGraph(); const rect=svg.getBoundingClientRect(); const vb=currentViewBox; return [vb[0]+(evt.clientX-rect.left)/rect.width*vb[2], vb[1]+(evt.clientY-rect.top)/rect.height*vb[3]]; }
function edgeStrokeWidth(e){ return 0.7 + Math.min(5.0, Math.log1p(Number(e.template_count_total||1))*0.9); }
function getNodeScaffold(nodeId){
  const n = byId.get(String(nodeId));
  if (!n) return '';
  return String(n.scaffold_group == null ? '' : n.scaffold_group).trim();
}
function isSameScaffoldEdge(e){
  let sg = String(e.source_scaffold_group == null ? '' : e.source_scaffold_group).trim();
  let tg = String(e.target_scaffold_group == null ? '' : e.target_scaffold_group).trim();
  if (!sg || sg === 'NA') sg = getNodeScaffold(e.source);
  if (!tg || tg === 'NA') tg = getNodeScaffold(e.target);
  if (sg && tg && sg !== 'NA' && tg !== 'NA') return sg === tg;
  const flag = String(e.same_group == null ? '' : e.same_group).toLowerCase();
  return flag === 'true' || flag === '1' || flag === 'yes';
}
function edgeColor(e){ return isSameScaffoldEdge(e) ? '#8f8f8f' : '#111111'; }
function edgeDimColor(e){ return isSameScaffoldEdge(e) ? '#d0d0d0' : '#a8a8a8'; }
function edgeDashArray(e){
  const layer = String(e.edge_layer || e.match_mode || e.profile || '').toLowerCase();
  if (layer === 'relaxed' || layer === 'loose') return '5,4';
  return '';
}
function nodeBaseColor(n){ return n.color || '#4E79A7'; }
function addLine(g,e,s,t){
  const line=document.createElementNS(svgNS,'line');
  line.setAttribute('x1',s.x); line.setAttribute('y1',s.y); line.setAttribute('x2',t.x); line.setAttribute('y2',t.y);
  const key=directedEdgeKey(e);
  const isPath=pathEdgeKeys.has(edgeKey(e.source,e.target));
  const isHi=highlightedEdgeKeys.has(key);
  const highlightMode=document.getElementById('displayMode').value==='highlight';
  const active=!highlightMode || isHi || isPath;
  line.setAttribute('class','edge'+(isPath?' path':'')+(isHi?' highlight':''));
  line.setAttribute('stroke',isPath?'#1f77b4':(active?edgeColor(e):edgeDimColor(e)));
  line.setAttribute('stroke-width',String(isHi?Math.max(2.5,edgeStrokeWidth(e)):edgeStrokeWidth(e)));
  line.setAttribute('opacity',active?'0.82':'0.16');
  line.setAttribute('fill','none');
  const dash=edgeDashArray(e);
  if(dash) line.setAttribute('stroke-dasharray',dash); else line.removeAttribute('stroke-dasharray');
  line.removeAttribute('marker-start');
  line.removeAttribute('marker-end');
  line.addEventListener('click',evt=>{evt.stopPropagation();showEdge(e);});
  g.appendChild(line); edgeEls.push({el:line,edge:e});
}
function addNode(g,n,showLabels){
  const c=document.createElementNS(svgNS,'circle');
  const sid=String(n.id);
  const isSelected=selectedIds.has(sid);
  const isPath=pathNodeIds.has(sid);
  const isHi=highlightedNodeIds.has(sid);
  const highlightMode=document.getElementById('displayMode').value==='highlight';
  const active=!highlightMode || isHi || isSelected || isPath;
  const baseRadius=Number(n.radius||6);
  const drawRadius=baseRadius;
  const fillColor=nodeBaseColor(n);
  c.setAttribute('cx',n.x); c.setAttribute('cy',n.y); c.setAttribute('r',drawRadius);
  c.setAttribute('fill',fillColor);
  c.setAttribute('fill-opacity',active?'1':'0.24');
  c.setAttribute('stroke',fillColor);
  c.setAttribute('stroke-opacity',active?'1':'0.32');
  c.setAttribute('stroke-width',isSelected?'3.2':(isPath?'2.4':(isHi?'2.0':'1.2')));
  c.style.stroke = fillColor;
  c.style.strokeWidth = (isSelected?'3.2':(isPath?'2.4':(isHi?'2.0':'1.2')));
  c.setAttribute('class','node'+(isSelected?' selected':'')+(isPath?' path':'')+(isHi?' highlight':''));
  c.addEventListener('mousedown',evt=>{evt.stopPropagation();draggedNode=n;});
  c.addEventListener('click',evt=>{evt.stopPropagation();showNode(n);});
  g.appendChild(c);
  let label=null;
  if(showLabels){
    label=document.createElementNS(svgNS,'text');
    label.setAttribute('x',Number(n.x));
    label.setAttribute('y',Number(n.y));
    label.setAttribute('class','label');
    label.setAttribute('opacity',active?'1':'0.35');
    label.setAttribute('style','font-size:'+NODE_LABEL_SIZE+'px;font-weight:normal;');
    label.textContent=String(n.label||n.source_id||n.id);
    g.appendChild(label);
  }
  nodeEls.push({circle:c,label:label,node:n});
}
function updatePositions(){ edgeEls.forEach(obj=>{const e=obj.edge,s=byId.get(String(e.source)),t=byId.get(String(e.target)); if(!s||!t)return; obj.el.setAttribute('x1',s.x); obj.el.setAttribute('y1',s.y); obj.el.setAttribute('x2',t.x); obj.el.setAttribute('y2',t.y);}); nodeEls.forEach(obj=>{const n=obj.node; obj.circle.setAttribute('cx',n.x); obj.circle.setAttribute('cy',n.y); if(obj.label){obj.label.setAttribute('x',Number(n.x)); obj.label.setAttribute('y',Number(n.y));}}); }
function render(){ const data=visible(); while(svg.firstChild)svg.removeChild(svg.firstChild); edgeEls=[]; nodeEls=[]; const g=document.createElementNS(svgNS,'g'); svg.appendChild(g); data.edges.forEach(e=>{const s=byId.get(String(e.source)),t=byId.get(String(e.target)); if(s&&t)addLine(g,e,s,t);}); const showLabels=document.getElementById('showLabels').checked; data.nodes.forEach(n=>addNode(g,n,showLabels)); }
function relaxVisibleLayout(steps=90){ const data=visible(); const ns=data.nodes, es=data.edges; const idset=new Set(ns.map(n=>String(n.id))); for(let step=0; step<steps; step++){const dx=new Map(ns.map(n=>[String(n.id),0])); const dy=new Map(ns.map(n=>[String(n.id),0])); for(let i=0;i<ns.length;i++)for(let j=i+1;j<ns.length;j++){const a=ns[i],b=ns[j]; let vx=Number(a.x)-Number(b.x),vy=Number(a.y)-Number(b.y); const d2=vx*vx+vy*vy+40,d=Math.sqrt(d2),f=1100/d2; vx=vx/d*f; vy=vy/d*f; dx.set(String(a.id),dx.get(String(a.id))+vx); dy.set(String(a.id),dy.get(String(a.id))+vy); dx.set(String(b.id),dx.get(String(b.id))-vx); dy.set(String(b.id),dy.get(String(b.id))-vy);} es.forEach(e=>{const a=byId.get(String(e.source)),b=byId.get(String(e.target)); if(!a||!b||!idset.has(String(a.id))||!idset.has(String(b.id)))return; const vx=Number(b.x)-Number(a.x),vy=Number(b.y)-Number(a.y),k=0.010; dx.set(String(a.id),dx.get(String(a.id))+vx*k); dy.set(String(a.id),dy.get(String(a.id))+vy*k); dx.set(String(b.id),dx.get(String(b.id))-vx*k); dy.set(String(b.id),dy.get(String(b.id))-vy*k);}); ns.forEach(n=>{const id=String(n.id); n.x+=Math.max(-10,Math.min(10,dx.get(id))); n.y+=Math.max(-10,Math.min(10,dy.get(id)));}); } updatePositions(); }
function directEdgesAmongSelected(){ const ids=Array.from(selectedIds); if(ids.length<2)return[]; const idset=new Set(ids); return edges.filter(e=>idset.has(String(e.source))&&idset.has(String(e.target))); }
function selectedSummaryHtml(){ const ns=Array.from(selectedIds).map(id=>byId.get(String(id))).filter(Boolean); const direct=directEdgesAmongSelected(); let out='<div class="kv"><b>Selected:</b> '+ns.map(n=>html(sourceLabel(n))).join(', ')+'</div>'; out+='<div class="kv"><b>Direct edges:</b> '+direct.length+'</div>'; direct.slice(0,20).forEach(e=>{const s=byId.get(String(e.source)),t=byId.get(String(e.target)); out+='<div class="detail-card"><div class="kv"><b>source_id:</b> '+html(e.source_label)+'</div><div class="kv"><b>target_id:</b> '+html(e.target_label)+'</div><div class="kv"><b>source scaffold:</b> '+html(s?s.scaffold_group:e.source_scaffold_group)+'</div><div class="kv"><b>target scaffold:</b> '+html(t?t.scaffold_group:e.target_scaffold_group)+'</div><div class="kv"><b>template support:</b> '+html(e.template_count_total||'')+'</div><div class="kv"><b>TSS:</b> '+html(e.transformation_support_score||'')+'</div><div class="kv"><b>R positions:</b> '+html(e.changed_positions)+'</div><div class="kv"><b>Sidechains:</b> '+html(e.sidechain_change_text)+'</div><div class="kv"><b>source SMILES:</b> '+html(s?(s.smiles||s.smiles_json):'')+'</div><div class="kv"><b>target SMILES:</b> '+html(t?(t.smiles||t.smiles_json):'')+'</div></div>';}); details.innerHTML=out; }
function showNode(n){ let out='<div class="detail-card"><div class="kv"><b>source_id:</b> '+html(sourceLabel(n))+'</div><div class="kv"><b>scaffold:</b> '+html(n.scaffold_group)+'</div><div class="kv"><b>degree:</b> '+html(n.degree_undirected)+'</div><div class="kv"><b>in_degree:</b> '+html(n.in_degree)+'</div><div class="kv"><b>out_degree:</b> '+html(n.out_degree)+'</div><div class="kv"><b>name:</b> '+html(n.name)+'</div><div class="kv"><b>SMILES:</b> '+html(n.smiles||n.smiles_json)+'</div></div>'; out+='<div class="detail-card"><h3>Molecule</h3>'+(n.mol_svg?'<div class="mol-svg">'+n.mol_svg+'</div>':'<pre>'+html(n.smiles||n.smiles_json)+'</pre>')+'</div>'; details.innerHTML=out; }
function showEdge(e){ const s=byId.get(String(e.source)); const t=byId.get(String(e.target)); const sLabel=e.source_label||(s?sourceLabel(s):e.source); const tLabel=e.target_label||(t?sourceLabel(t):e.target); const sSmiles=s?(s.smiles||s.smiles_json||''):''; const tSmiles=t?(t.smiles||t.smiles_json||''):''; let out='<div class="detail-card"><div class="kv"><b>source_id:</b> '+html(sLabel)+'</div><div class="kv"><b>target_id:</b> '+html(tLabel)+'</div><div class="kv"><b>source scaffold:</b> '+html(s?s.scaffold_group:e.source_scaffold_group)+'</div><div class="kv"><b>target scaffold:</b> '+html(t?t.scaffold_group:e.target_scaffold_group)+'</div><div class="kv"><b>edge layer:</b> '+html(e.edge_layer)+'</div><div class="kv"><b>template support:</b> '+html(e.template_count_total||'')+'</div><div class="kv"><b>TSS:</b> '+html(e.transformation_support_score||'')+'</div><div class="kv"><b>changed count:</b> '+html(e.changed_position_count||'')+'</div><div class="kv"><b>R positions:</b> '+html(e.changed_positions)+'</div><div class="kv"><b>Sidechains:</b> '+html(e.sidechain_change_text)+'</div><div class="kv"><b>source SMILES:</b> '+html(sSmiles)+'</div><div class="kv"><b>target SMILES:</b> '+html(tSmiles)+'</div></div>'; out+='<div class="detail-grid"><div class="detail-card"><h3>Source molecule</h3>'+(s&&s.mol_svg?'<div class="mol-svg">'+s.mol_svg+'</div>':'<pre>'+html(sSmiles)+'</pre>')+'</div><div class="detail-card"><h3>Target molecule</h3>'+(t&&t.mol_svg?'<div class="mol-svg">'+t.mol_svg+'</div>':'<pre>'+html(tSmiles)+'</pre>')+'</div></div>'; details.innerHTML=out; }
function renderSearchResults(hits){ const box=document.getElementById('results'); box.innerHTML=''; hits.slice(0,80).forEach(n=>{const d=document.createElement('div'); d.className='result-item'; d.textContent=String(sourceLabel(n))+' | scaffold '+String(n.scaffold_group); d.onclick=()=>{selectedIds.add(String(n.id)); updatePaths(); selectedSummaryHtml(); render();}; box.appendChild(d);}); if(hits.length>80){const d=document.createElement('div'); d.className='small'; d.textContent='Showing first 80 of '+hits.length+' matches.'; box.appendChild(d);} }
function updatePositionOptions(){ const current=document.getElementById('position').value; const set=new Set(); edges.forEach(e=>parseChanges(e).forEach(c=>set.add(c.position))); const sel=document.getElementById('position'); sel.innerHTML='<option value="all">all</option>'; Array.from(set).sort((a,b)=>a.localeCompare(b,undefined,{numeric:true})).forEach(v=>{const o=document.createElement('option'); o.value=v; o.textContent=v; sel.appendChild(o);}); if(Array.from(set).includes(current))sel.value=current; }
function updateSidechainOptions(){ const current=document.getElementById('sidechain').value; const pos=document.getElementById('position').value; const direction=document.getElementById('directionMode').value; const set=new Set(); edges.forEach(e=>parseChanges(e).forEach(c=>{if(pos!=='all'&&c.position!==pos)return; if(direction==='upstream')set.add(c.to); else if(direction==='downstream')set.add(c.from); else {set.add(c.from); set.add(c.to);}})); const sel=document.getElementById('sidechain'); sel.innerHTML='<option value="all">all</option>'; Array.from(set).sort().forEach(v=>{const o=document.createElement('option'); o.value=v; o.textContent=v; sel.appendChild(o);}); if(Array.from(set).includes(current))sel.value=current; }
function initSelectors(){ const layerSel=document.getElementById('edgeLayer'); const names=allLayerNames(); names.forEach(n=>{const o=document.createElement('option'); o.value=n; o.textContent=n; layerSel.appendChild(o);}); if(names.length>1){const o=document.createElement('option'); o.value='overlay'; o.textContent='overlay'; layerSel.appendChild(o);} currentLayer=names.includes('strict')?'strict':(names[0]||'current'); layerSel.value=currentLayer; setCurrentEdges(); const groups=Array.from(new Set(nodes.map(n=>String(n.scaffold_group)))).sort((a,b)=>a.localeCompare(b,undefined,{numeric:true})); const sg=document.getElementById('scaffold'); groups.forEach(g=>{const o=document.createElement('option'); o.value=g; o.textContent=g; sg.appendChild(o);}); }
function initLegend(){ const box=document.getElementById('topLegend'); const groups=Array.from(new Set(nodes.map(n=>String(n.scaffold_group)))).sort((a,b)=>a.localeCompare(b,undefined,{numeric:true})); groups.forEach(g=>{const n=nodes.find(x=>String(x.scaffold_group)===g); const div=document.createElement('div'); div.className='small'; div.innerHTML='<span class="swatch" style="background:'+(n?n.color:'#999')+'"></span>Scaffold '+html(g); box.appendChild(div);}); }
function refresh(){ updateSidechainOptions(); updatePaths(); render(); fitGraph(); }
function setupEvents(){ document.getElementById('edgeLayer').addEventListener('change',()=>{setCurrentEdges(); render(); fitGraph();}); document.getElementById('applySearch').onclick=()=>{selectedIds.clear(); const hits=searchNodes(); hits.forEach(n=>selectedIds.add(String(n.id))); renderSearchResults(hits); updatePaths(); selectedSummaryHtml(); render(); centerOnSelectedNode();}; document.getElementById('clearSearch').onclick=()=>{selectedIds.clear(); pathNodeIds.clear(); pathEdgeKeys.clear(); document.getElementById('search').value=''; document.getElementById('results').innerHTML=''; details.innerHTML='Click a node or edge.'; render(); fitGraph();}; document.getElementById('fitSelected').onclick=fitSelected; document.getElementById('fit').onclick=fitGraph; document.getElementById('relax').onclick=()=>relaxVisibleLayout(90); document.getElementById('reset').onclick=()=>{['scaffold','position','sidechain'].forEach(id=>document.getElementById(id).value='all'); document.getElementById('cross').checked=false; document.getElementById('focusMode').value='all'; document.getElementById('explorerMode').value='standard'; document.getElementById('directionMode').value='downstream'; document.getElementById('hopDepth').value='1'; document.getElementById('displayMode').value='filter'; selectedIds.clear(); pathNodeIds.clear(); pathEdgeKeys.clear(); updateSidechainOptions(); render(); fitGraph();}; ['focusMode','scaffold','position','sidechain','cross','showLabels','explorerMode','directionMode','hopDepth','displayMode'].forEach(id=>document.getElementById(id).addEventListener('change',refresh)); svg.addEventListener('click',()=>{details.innerHTML='Click a node or edge.';}); svg.addEventListener('wheel',evt=>{evt.preventDefault(); if(!currentViewBox)fitGraph(); const vb=currentViewBox.slice(); const factor=evt.deltaY<0?0.88:1.14; const rect=svg.getBoundingClientRect(); const mx=(evt.clientX-rect.left)/rect.width; const my=(evt.clientY-rect.top)/rect.height; const wx=vb[0]+mx*vb[2]; const wy=vb[1]+my*vb[3]; vb[2]*=factor; vb[3]*=factor; vb[0]=wx-mx*vb[2]; vb[1]=wy-my*vb[3]; currentViewBox=vb; svg.setAttribute('viewBox',vb.join(' '));}); svg.addEventListener('mousedown',evt=>{if(draggedNode)return; isPanning=true; panStart=[evt.clientX,evt.clientY,currentViewBox?currentViewBox.slice():null];}); window.addEventListener('mouseup',()=>{isPanning=false; draggedNode=null;}); window.addEventListener('mousemove',evt=>{if(draggedNode){const p=svgPoint(evt); draggedNode.x=p[0]; draggedNode.y=p[1]; updatePositions(); return;} if(!isPanning||!panStart||!panStart[2])return; const rect=svg.getBoundingClientRect(); const vb=panStart[2].slice(); const dx=(evt.clientX-panStart[0])/rect.width*vb[2]; const dy=(evt.clientY-panStart[1])/rect.height*vb[3]; vb[0]-=dx; vb[1]-=dy; currentViewBox=vb; svg.setAttribute('viewBox',vb.join(' '));}); }
function initScaffoldImage(){ const wrap=document.getElementById('scaffoldImageWrap'); if(!wrap||!GRAPH.scaffold_image_uri)return; const img=document.createElement('img'); img.src=GRAPH.scaffold_image_uri; img.alt='Scaffold reference'; wrap.appendChild(img); }
function startup(){ initScaffoldImage(); initSelectors(); initLegend(); setupEvents(); render(); fitGraph(); }
startup();
</script>
</body>
</html>'''
    return html.replace("__TITLE__", safe_title).replace("__DATA__", payload).replace("__LABEL_SIZE__", str(int(label_size)))


def visualize_existing_network(
    network_dir: Path,
    scaffold_json: Path,
    output_html: Path,
    valid_molecules: Optional[Path] = None,
    json_index_base: int = 1,
    edge_id_space: str = "auto",
    max_nodes: int = 2500,
    max_edges: int = 12000,
    layout_seed: int = 17,
    central_groups: int = 5,
    label_column: str = "source_id",
    title: str = "MolSpaceHub TransformNet interactive view",
    make_svg: bool = True,
    scaffold_image: Optional[Path] = None,
    strict_network_dir: Optional[Path] = None,
    relaxed_network_dir: Optional[Path] = None,
    node_size_scale: float = 1.35,
    label_size: int = 10,
) -> Dict[str, str]:
    network_dir = Path(network_dir)
    node_size_scale = safe_float(os.environ.get("TRANSFORMNET_NODE_SIZE_SCALE", node_size_scale), node_size_scale)
    label_size = safe_int(os.environ.get("TRANSFORMNET_LABEL_SIZE", label_size), label_size)
    nodes_dir = Path(strict_network_dir) if strict_network_dir else network_dir
    raw_nodes = read_csv(nodes_table_path(nodes_dir))
    scaffold_rows, smiles_index = read_scaffold_json(scaffold_json, json_index_base=json_index_base)
    matrix_to_source = load_valid_molecule_map(valid_molecules, nodes_dir)
    node_id_to_raw, row_index_to_node_id, source_to_node_id, label_to_node_id, matrix_to_node_id = build_node_maps(raw_nodes, label_column)

    visual_nodes: List[Dict[str, Any]] = []
    mapping_qc: List[Dict[str, Any]] = []
    temp_nodes: List[Tuple[int, Dict[str, str], str, Dict[str, Any], Dict[str, str], str]] = []
    group_counts: Counter = Counter()
    for idx, row in enumerate(raw_nodes[:max_nodes]):
        nid = row["__node_id"]
        scaf, qc = attach_scaffold(row, scaffold_rows, smiles_index, json_index_base=json_index_base)
        group = clean(scaf.get("scaffold_group") or "NA")
        group_counts[group] += 1
        temp_nodes.append((idx, row, nid, scaf, qc, group))
    color_map, major_groups = build_colors_by_group(group_counts, central_groups)

    for idx, row, nid, scaf, qc, group in temp_nodes:
        degree = safe_int(row.get("degree_undirected") or row.get("degree") or row.get("in_degree_directed"), 0)
        base_radius = 5.0 + min(14.0, 1.7 * math.sqrt(max(0, degree)))
        display_radius = round(base_radius * max(0.2, node_size_scale), 3)
        label = clean(row.get(label_column)) if label_column in row else nid
        if not label:
            label = clean(row.get("source_id") or row.get("molecule_id") or nid)
        smiles = choose_smiles(row) or clean(scaf.get("smiles_json"))
        node = {
            "id": nid,
            "label": label,
            "row_index": idx,
            "source_id": clean(row.get("source_id")) or label,
            "molecule_id": clean(row.get("molecule_id")),
            "matrix_index": clean(row.get("matrix_index")),
            "name": choose_name(row),
            "smiles": smiles,
            "scaffold_group": group,
            "scaffold_label": clean(scaf.get("scaffold_label") or f"Scaffold_{group}"),
            "scaffold_group_raw": clean(scaf.get("scaffold_group_raw")),
            "degree_undirected": degree,
            "in_degree": safe_int(row.get("in_degree_directed") or row.get("in_degree"), 0),
            "out_degree": safe_int(row.get("out_degree_directed") or row.get("out_degree"), 0),
            "component_id": clean(row.get("component_id")),
            "core_smarts": clean(scaf.get("core_smarts")),
            "smiles_json": clean(scaf.get("smiles_json")),
            "json_index0": clean(scaf.get("json_index0")),
            "json_index1": clean(scaf.get("json_index1")),
            "mapping_status": qc["mapping_status"],
            "mapping_reason": qc["mapping_reason"],
            "color": color_map.get(group, "#999999"),
            "major_scaffold_group": group in major_groups,
            "radius": display_radius,
            "chain": scaf.get("chain", {}),
            "mol_svg": make_mol_svg(smiles) if make_svg else "",
        }
        visual_nodes.append(node)
        mapping_qc.append({"node_id": nid, "label": label, "scaffold_group": group, "json_index0": node["json_index0"], "json_index1": node["json_index1"], "mapping_status": qc["mapping_status"], "mapping_reason": qc["mapping_reason"]})

    included_ids = {clean(n["id"]) for n in visual_nodes}
    node_info = {clean(n["id"]): n for n in visual_nodes}

    if strict_network_dir or relaxed_network_dir:
        edge_inputs: Dict[str, Path] = {}
        if strict_network_dir:
            edge_inputs["strict"] = edge_table_path(Path(strict_network_dir))
        if relaxed_network_dir:
            edge_inputs["relaxed"] = edge_table_path(Path(relaxed_network_dir))
    else:
        edge_inputs = {"current": edge_table_path(network_dir)}

    edge_layers: Dict[str, List[Dict[str, Any]]] = {}
    all_sidechain_rows: List[Dict[str, Any]] = []
    all_unresolved_edges: List[Dict[str, Any]] = []
    selected_edge_id_space = edge_id_space
    edge_scores: Dict[str, int] = {}
    for layer, path in edge_inputs.items():
        raw_edges = read_csv(path)
        selected_edge_id_space, layer_scores = choose_edge_id_space(raw_edges, edge_id_space, set(node_id_to_raw), row_index_to_node_id, matrix_to_source, source_to_node_id, label_to_node_id, matrix_to_node_id)
        edge_scores[layer] = layer_scores.get(selected_edge_id_space, 0)
        visual_edges, sidechain_rows, unresolved_edges = build_visual_edges(raw_edges, layer, node_id_to_raw, row_index_to_node_id, source_to_node_id, label_to_node_id, matrix_to_node_id, matrix_to_source, selected_edge_id_space, included_ids, node_info, max_edges)
        edge_layers[layer] = visual_edges
        all_sidechain_rows.extend(sidechain_rows)
        all_unresolved_edges.extend(unresolved_edges)

    normalize_support_scores(edge_layers)
    layout_edges = edge_layers.get("strict") or next(iter(edge_layers.values()), [])
    compute_circular_group_layout(visual_nodes, layout_edges, layout_seed, central_groups)

    summary = {
        "nodes": len(visual_nodes),
        "edge_layers": {k: len(v) for k, v in edge_layers.items()},
        "edges_total": sum(len(v) for v in edge_layers.values()),
        "scaffold_groups": len({clean(n.get("scaffold_group")) for n in visual_nodes}),
        "mapped_nodes": sum(1 for n in visual_nodes if n.get("mapping_status") != "unmapped"),
        "major_scaffold_groups": sorted(list(major_groups), key=natural_key),
        "selected_edge_id_space": selected_edge_id_space,
        "edge_id_scores": edge_scores,
        "unresolved_edges": len(all_unresolved_edges),
        "html_renderer": "self_contained_transformation_explorer_svg",
        "molecule_svg_enabled": bool(make_svg),
        "support_score_definition": "normalized log1p(template_support) * match_factor * 1/max(1, changed_position_count)",
    }
    image_path = scaffold_image if scaffold_image is not None else (network_dir / "taxane_scaffold.png")
    scaffold_image_uri = image_to_data_uri(image_path)
    summary["scaffold_image_embedded"] = bool(scaffold_image_uri)
    html = build_html(visual_nodes, edge_layers, title=title, summary=summary, edge_scores=edge_scores, scaffold_image_uri=scaffold_image_uri)
    output_html = Path(output_html)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text(html, encoding="utf-8")
    prefix = output_html.parent / "transformnet_network"
    write_csv(prefix.parent / "transformnet_network.visual_nodes.csv", visual_nodes, skip_keys={"chain", "mol_svg"})
    flat_edges = []
    for layer, es in edge_layers.items():
        flat_edges.extend(es)
    write_csv(prefix.parent / "transformnet_network.visual_edges.csv", flat_edges)
    write_csv(prefix.parent / "transformnet_network.sidechain_changes.long.csv", all_sidechain_rows)
    write_csv(prefix.parent / "transformnet_network.scaffold_mapping_qc.csv", mapping_qc)
    write_csv(prefix.parent / "transformnet_network.unresolved_edges.csv", all_unresolved_edges)
    (prefix.parent / "transformnet_network.visual_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=True), encoding="utf-8")
    return {"html": str(output_html), "visual_nodes": str(prefix.parent / "transformnet_network.visual_nodes.csv"), "visual_edges": str(prefix.parent / "transformnet_network.visual_edges.csv"), "sidechain_changes": str(prefix.parent / "transformnet_network.sidechain_changes.long.csv"), "mapping_qc": str(prefix.parent / "transformnet_network.scaffold_mapping_qc.csv"), "unresolved_edges": str(prefix.parent / "transformnet_network.unresolved_edges.csv"), "summary": str(prefix.parent / "transformnet_network.visual_summary.json")}
