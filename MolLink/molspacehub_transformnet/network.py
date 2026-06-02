from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import networkx as nx
import pandas as pd

from .io import load_annotation_table
from .utils import ensure_dir, semicolon_join, write_json


@dataclass
class NetworkBuildResult:
    paths: dict[str, Path]
    qc: dict[str, Any]


def _detect_join_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    normalized = {str(c).strip().lower(): c for c in df.columns}
    for cand in candidates:
        key = cand.lower()
        if key in normalized:
            return normalized[key]
    return None


def _load_computed_dir(computed_dir: str | Path, prefix: str = "transformnet") -> tuple[pd.DataFrame, pd.DataFrame]:
    computed_dir = Path(computed_dir)
    valid_path = computed_dir / f"{prefix}.valid_molecules.csv"
    edges_path = computed_dir / f"{prefix}.directed_edges.csv"
    if not valid_path.exists():
        candidates = sorted(computed_dir.glob("*.valid_molecules.csv"))
        if candidates:
            valid_path = candidates[0]
    if not edges_path.exists():
        candidates = sorted(computed_dir.glob("*.directed_edges.csv"))
        if candidates:
            edges_path = candidates[0]
    if not valid_path.exists() or not edges_path.exists():
        raise FileNotFoundError(f"Could not find computed valid_molecules/direct_edges in {computed_dir}")
    return pd.read_csv(valid_path), pd.read_csv(edges_path)


def _merge_annotation(
    nodes: pd.DataFrame,
    annotation_table: str | Path | None = None,
    annotation_sheet: str | int | None = None,
    annotation_id_column: str | None = None,
    molecule_id_column: str = "molecule_id",
) -> pd.DataFrame:
    if not annotation_table:
        return nodes
    ann = load_annotation_table(annotation_table, sheet_name=annotation_sheet)
    if ann.empty:
        return nodes
    ann = ann.copy()

    annotation_id_column = annotation_id_column or _detect_join_column(
        ann,
        ["molecule_id", "source_id", "source_no", "library_index", "id", "ID", "compound_id"],
    )
    if not annotation_id_column:
        raise ValueError(f"Could not detect annotation ID column. Available columns: {ann.columns.tolist()}")

    left = nodes.copy()
    left["_join_id"] = left[molecule_id_column].astype(str)
    ann["_join_id"] = ann[annotation_id_column].astype(str)

    # Keep annotation columns under original names unless collision occurs.
    rename = {}
    for col in ann.columns:
        if col == "_join_id":
            continue
        if col in left.columns:
            rename[col] = f"annotation_{col}"
    ann = ann.rename(columns=rename)
    return left.merge(ann, on="_join_id", how="left").drop(columns=["_join_id"])


def _infer_group_column(nodes: pd.DataFrame, group_column: str | None = None) -> str | None:
    if group_column:
        return group_column if group_column in nodes.columns else None
    candidates = ["scaffold_group", "annotation_scaffold_group", "group", "cluster", "class", "family"]
    return _detect_join_column(nodes, candidates)


def _infer_position_columns(nodes: pd.DataFrame, position_columns: list[str] | None = None) -> list[str]:
    if position_columns:
        return [c for c in position_columns if c in nodes.columns]
    out = []
    for c in nodes.columns:
        text = str(c)
        if re.fullmatch(r"R\d+", text, flags=re.IGNORECASE) or re.fullmatch(r"C\d+", text, flags=re.IGNORECASE):
            out.append(c)
    return sorted(out, key=lambda x: (re.sub(r"\d+", "", str(x)), int(re.findall(r"\d+", str(x))[0]) if re.findall(r"\d+", str(x)) else 0))


def _node_metrics(valid: pd.DataFrame, directed_edges: pd.DataFrame) -> pd.DataFrame:
    nodes = valid.copy()
    nodes["in_degree_directed"] = 0
    nodes["out_degree_directed"] = 0
    nodes["degree_undirected"] = 0
    if directed_edges.empty:
        return nodes

    out_counts = directed_edges.groupby("substrate_matrix_index").size()
    in_counts = directed_edges.groupby("product_matrix_index").size()
    edge_pairs = set()
    for _, e in directed_edges.iterrows():
        s = int(e["substrate_matrix_index"])
        p = int(e["product_matrix_index"])
        if s == p:
            edge_pairs.add((s, p))
        else:
            edge_pairs.add(tuple(sorted((s, p))))
    deg = {}
    for a, b in edge_pairs:
        deg[a] = deg.get(a, 0) + 1
        if b != a:
            deg[b] = deg.get(b, 0) + 1

    nodes["in_degree_directed"] = nodes["matrix_index"].map(in_counts).fillna(0).astype(int)
    nodes["out_degree_directed"] = nodes["matrix_index"].map(out_counts).fillna(0).astype(int)
    nodes["degree_undirected"] = nodes["matrix_index"].map(deg).fillna(0).astype(int)
    return nodes


def _component_summary(nodes: pd.DataFrame, directed_edges: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    G = nx.Graph()
    for _, n in nodes.iterrows():
        G.add_node(int(n["matrix_index"]))
    for _, e in directed_edges.iterrows():
        G.add_edge(int(e["substrate_matrix_index"]), int(e["product_matrix_index"]))
    comp_map = {}
    records = []
    for comp_id, members in enumerate(nx.connected_components(G), start=1):
        members = sorted(members)
        for m in members:
            comp_map[m] = comp_id
        records.append({
            "component_id": comp_id,
            "size": len(members),
            "members": semicolon_join(members, max_items=100),
        })
    nodes = nodes.copy()
    nodes["component_id"] = nodes["matrix_index"].map(comp_map).fillna(-1).astype(int)
    comps = pd.DataFrame(records).sort_values(["size", "component_id"], ascending=[False, True]) if records else pd.DataFrame(columns=["component_id", "size", "members"])
    return nodes, comps


def _collapse_edges(
    nodes: pd.DataFrame,
    directed_edges: pd.DataFrame,
    group_column: str | None,
    position_columns: list[str],
) -> pd.DataFrame:
    if directed_edges.empty:
        return pd.DataFrame()

    node_by_idx = nodes.set_index("matrix_index", drop=False)
    records = []
    for (a, b), g in directed_edges.assign(
        _a=directed_edges[["substrate_matrix_index", "product_matrix_index"]].min(axis=1),
        _b=directed_edges[["substrate_matrix_index", "product_matrix_index"]].max(axis=1),
    ).groupby(["_a", "_b"], sort=True):
        a, b = int(a), int(b)
        n1 = node_by_idx.loc[a]
        n2 = node_by_idx.loc[b]
        forward = ((g["substrate_matrix_index"] == a) & (g["product_matrix_index"] == b)).any()
        reverse = ((g["substrate_matrix_index"] == b) & (g["product_matrix_index"] == a)).any()
        if forward and reverse:
            directionality = "bidirectional"
        elif forward:
            directionality = f"{n1['molecule_id']}->{n2['molecule_id']}"
        else:
            directionality = f"{n2['molecule_id']}->{n1['molecule_id']}"

        changed = []
        change_detail = []
        for col in position_columns:
            v1 = str(n1.get(col, ""))
            v2 = str(n2.get(col, ""))
            if v1 != v2:
                changed.append(col)
                change_detail.append(f"{col}:{v1}->{v2}")

        group1 = n1.get(group_column, "") if group_column else ""
        group2 = n2.get(group_column, "") if group_column else ""
        records.append({
            "node1_matrix_index": a,
            "node2_matrix_index": b,
            "node1_id": n1["molecule_id"],
            "node2_id": n2["molecule_id"],
            "node1_name": n1.get("molecule_name", n1["molecule_id"]),
            "node2_name": n2.get("molecule_name", n2["molecule_id"]),
            "directionality": directionality,
            "has_bidirectional_support": bool(forward and reverse),
            "directed_edge_count": int(len(g)),
            "event_count_total": int(g.get("event_count", pd.Series([len(g)])).sum()),
            "template_count_total": int(g.get("template_count", pd.Series([0])).sum()),
            "template_uids": semicolon_join(";".join(map(str, g.get("template_uids", []))).split(";"), max_items=80),
            "product_keys": semicolon_join(";".join(map(str, g.get("product_keys", []))).split(";"), max_items=80),
            "group_column": group_column or "",
            "node1_group": group1,
            "node2_group": group2,
            "same_group": bool(group_column and str(group1) == str(group2)),
            "changed_positions": semicolon_join(changed),
            "changed_position_count": len(changed),
            "change_detail": semicolon_join(change_detail, max_items=80),
        })
    return pd.DataFrame(records)


def _group_edge_summary(collapsed_edges: pd.DataFrame, group_column: str | None) -> pd.DataFrame:
    if not group_column or collapsed_edges.empty:
        return pd.DataFrame()
    df = collapsed_edges.copy()
    df["group_pair"] = df.apply(lambda r: " -- ".join(sorted([str(r["node1_group"]), str(r["node2_group"])])), axis=1)
    out = (
        df.groupby(["node1_group", "node2_group"], dropna=False)
        .agg(
            edge_count=("node1_id", "size"),
            bidirectional_edges=("has_bidirectional_support", "sum"),
            mean_changed_positions=("changed_position_count", "mean"),
            total_event_count=("event_count_total", "sum"),
            total_template_count=("template_count_total", "sum"),
        )
        .reset_index()
    )
    return out


def _write_graphs(nodes: pd.DataFrame, edges: pd.DataFrame, output_dir: Path, prefix: str) -> dict[str, Path]:
    paths = {}
    G = nx.Graph()
    for _, n in nodes.iterrows():
        attrs = {k: ("" if pd.isna(v) else v) for k, v in n.to_dict().items()}
        G.add_node(str(n["molecule_id"]), **attrs)
    for _, e in edges.iterrows():
        attrs = {k: ("" if pd.isna(v) else v) for k, v in e.to_dict().items() if k not in {"node1_id", "node2_id"}}
        G.add_edge(str(e["node1_id"]), str(e["node2_id"]), **attrs)
    paths["graphml"] = output_dir / f"{prefix}.graphml"
    paths["gexf"] = output_dir / f"{prefix}.gexf"
    try:
        nx.write_graphml(G, paths["graphml"])
    except Exception:
        paths.pop("graphml", None)
    try:
        nx.write_gexf(G, paths["gexf"])
    except Exception:
        paths.pop("gexf", None)
    return paths


def _write_html_preview(nodes: pd.DataFrame, edges: pd.DataFrame, path: Path, max_edges: int = 3000) -> Path:
    node_cols = [c for c in ["molecule_id", "molecule_name", "degree_undirected", "component_id"] if c in nodes.columns]
    node_records = nodes[node_cols].head(5000).to_dict(orient="records") if node_cols else []
    edge_cols = ["node1_id", "node2_id", "directionality", "template_count_total", "changed_positions"]
    if edges.empty:
        edge_records = []
    else:
        edge_records = edges.reindex(columns=edge_cols, fill_value="").head(max_edges).to_dict(orient="records")
    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>MolSpaceHub TransformNet preview</title>
<style>body{{font-family:Arial,sans-serif;margin:24px}} table{{border-collapse:collapse;font-size:12px}} td,th{{border:1px solid #ccc;padding:4px 6px}} th{{background:#f3f3f3}}</style></head>
<body>
<h1>MolSpaceHub TransformNet preview</h1>
<p>Nodes: {len(nodes)}; collapsed edges: {len(edges)}. This lightweight preview lists the first {min(max_edges, len(edges))} edges.</p>
<h2>Edges</h2>
<table><thead><tr><th>node1</th><th>node2</th><th>directionality</th><th>template_count_total</th><th>changed_positions</th></tr></thead><tbody>
{''.join(f"<tr><td>{e['node1_id']}</td><td>{e['node2_id']}</td><td>{e['directionality']}</td><td>{e['template_count_total']}</td><td>{e['changed_positions']}</td></tr>" for e in edge_records)}
</tbody></table>
</body></html>"""
    path.write_text(html, encoding="utf-8")
    return path


def _write_excel(path: Path, sheets: dict[str, pd.DataFrame]) -> Path:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, df in sheets.items():
            safe = name[:31]
            df.to_excel(writer, sheet_name=safe, index=False)
    return path


def build_network_from_computed(
    computed_dir: str | Path | None = None,
    valid_molecules: str | Path | None = None,
    directed_edges: str | Path | None = None,
    output_dir: str | Path = "outputs/network",
    annotation_table: str | Path | None = None,
    annotation_sheet: str | int | None = None,
    annotation_id_column: str | None = None,
    group_column: str | None = None,
    position_columns: list[str] | None = None,
    prefix: str = "transformnet_network",
    write_excel: bool = True,
    write_graphs: bool = True,
    write_html: bool = True,
) -> NetworkBuildResult:
    output_dir = ensure_dir(output_dir)

    if computed_dir:
        valid, directed = _load_computed_dir(computed_dir)
    else:
        if not valid_molecules or not directed_edges:
            raise ValueError("Provide either computed_dir or both valid_molecules and directed_edges.")
        valid = pd.read_csv(valid_molecules)
        directed = pd.read_csv(directed_edges)

    if "matrix_index" not in valid.columns:
        valid["matrix_index"] = range(len(valid))
    if directed.empty:
        directed = pd.DataFrame(columns=["substrate_matrix_index", "product_matrix_index"])

    nodes = _node_metrics(valid, directed)
    nodes = _merge_annotation(
        nodes,
        annotation_table=annotation_table,
        annotation_sheet=annotation_sheet,
        annotation_id_column=annotation_id_column,
    )
    group_column_eff = _infer_group_column(nodes, group_column=group_column)
    position_columns_eff = _infer_position_columns(nodes, position_columns=position_columns)
    nodes, components = _component_summary(nodes, directed)
    collapsed = _collapse_edges(nodes, directed, group_column_eff, position_columns_eff)
    group_edges = _group_edge_summary(collapsed, group_column_eff)

    paths: dict[str, Path] = {
        "nodes": output_dir / f"{prefix}.nodes.csv",
        "directed_edges": output_dir / f"{prefix}.directed_edges.csv",
        "collapsed_edges": output_dir / f"{prefix}.collapsed_edges.csv",
        "components": output_dir / f"{prefix}.components.csv",
        "group_edges": output_dir / f"{prefix}.group_edges.csv",
        "summary": output_dir / f"{prefix}.summary.json",
    }
    nodes.to_csv(paths["nodes"], index=False)
    directed.to_csv(paths["directed_edges"], index=False)
    collapsed.to_csv(paths["collapsed_edges"], index=False)
    components.to_csv(paths["components"], index=False)
    group_edges.to_csv(paths["group_edges"], index=False)

    if write_excel:
        paths["excel"] = output_dir / f"{prefix}.xlsx"
        _write_excel(paths["excel"], {
            "Nodes": nodes,
            "Collapsed_Edges": collapsed,
            "Directed_Edges": directed,
            "Components": components,
            "Group_Edges": group_edges,
            "QC": pd.DataFrame([{
                "nodes": len(nodes),
                "directed_edges": len(directed),
                "collapsed_edges": len(collapsed),
                "components": len(components),
                "group_column": group_column_eff or "",
                "position_columns": ";".join(position_columns_eff),
            }]),
        })

    if write_graphs:
        paths.update(_write_graphs(nodes, collapsed, output_dir, prefix))

    if write_html:
        paths["html"] = _write_html_preview(nodes, collapsed, output_dir / f"{prefix}.html")

    qc = {
        "module_version": "0.6.1",
        "mode": "network_build",
        "nodes": int(len(nodes)),
        "directed_edges": int(len(directed)),
        "collapsed_edges": int(len(collapsed)),
        "components": int(len(components)),
        "largest_component_size": int(components["size"].max()) if not components.empty else 0,
        "isolated_nodes": int((nodes["degree_undirected"] == 0).sum()) if "degree_undirected" in nodes else 0,
        "group_column": group_column_eff,
        "position_columns": position_columns_eff,
        "same_group_edges": int(collapsed["same_group"].sum()) if "same_group" in collapsed else 0,
        "cross_group_edges": int((~collapsed["same_group"]).sum()) if "same_group" in collapsed and not collapsed.empty else 0,
        "outputs": {k: str(v) for k, v in paths.items()},
    }
    write_json(qc, paths["summary"])
    return NetworkBuildResult(paths=paths, qc=qc)
