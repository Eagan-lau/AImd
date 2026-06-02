from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .utils import ensure_dir, write_json


@dataclass
class HotspotResult:
    paths: dict[str, Path]
    qc: dict[str, Any]


def _clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "null"}:
        return ""
    return text


def _first_existing(base: Path, names: list[str]) -> Path:
    for name in names:
        path = base / name
        if path.exists():
            return path
    raise FileNotFoundError(f"None of these files were found in {base}: {', '.join(names)}")


def _prepare_group_rows(sidechain: pd.DataFrame, include_cross_edges: bool) -> pd.DataFrame:
    if sidechain.empty:
        return pd.DataFrame(columns=[
            "group", "position", "source_sidechain", "target_sidechain", "source_id", "target_id",
            "same_group", "template_count_total", "group_side"
        ])
    rows = []
    for _, r in sidechain.iterrows():
        source_group = _clean(r.get("source_scaffold_group")) or "NA"
        target_group = _clean(r.get("target_scaffold_group")) or "NA"
        same = str(_clean(r.get("same_group"))).lower() in {"true", "1", "yes"} or source_group == target_group
        base = {
            "position": _clean(r.get("position")),
            "source_sidechain": _clean(r.get("source_sidechain")),
            "target_sidechain": _clean(r.get("target_sidechain")),
            "source_id": _clean(r.get("source_id")),
            "target_id": _clean(r.get("target_id")),
            "same_group": same,
            "template_count_total": float(r.get("template_count_total", 1) or 1),
        }
        if same:
            rows.append({**base, "group": source_group, "group_side": "same_group"})
        elif include_cross_edges:
            rows.append({**base, "group": source_group, "group_side": "source_group"})
            rows.append({**base, "group": target_group, "group_side": "target_group"})
    return pd.DataFrame(rows)


def compute_group_position_hotspots(
    network_dir: str | Path,
    output_dir: str | Path | None = None,
    prefix: str = "transformnet_network",
    include_cross_edges: bool = True,
) -> HotspotResult:
    """Compute group-level R-position hotspot coefficients from visualized side-chain changes.

    The primary count is the number of side-chain change records assigned to a scaffold group
    for each R position. The coefficient is normalized within each scaffold group:

        hotspot_coefficient = group_position_edge_count / total_group_position_edge_count

    Cross-scaffold edges can be assigned to both endpoint groups. Same-scaffold edges are
    assigned once to that shared group.
    """
    network_dir = Path(network_dir)
    output_dir = ensure_dir(output_dir or network_dir)
    changes_path = _first_existing(network_dir, [
        f"{prefix}.sidechain_changes.long.csv",
        "transformnet_network.sidechain_changes.long.csv",
    ])
    try:
        sidechain = pd.read_csv(changes_path)
    except pd.errors.EmptyDataError:
        sidechain = pd.DataFrame(columns=[
            "source_id", "target_id", "position", "source_sidechain", "target_sidechain",
            "source_scaffold_group", "target_scaffold_group", "same_group", "template_count_total"
        ])
    group_rows = _prepare_group_rows(sidechain, include_cross_edges=include_cross_edges)

    paths: dict[str, Path] = {
        "group_position_rows": output_dir / f"{prefix}.group_position_change_rows.csv",
        "group_position_hotspots": output_dir / f"{prefix}.group_position_hotspots.csv",
        "group_position_transitions": output_dir / f"{prefix}.group_position_transitions.csv",
        "summary": output_dir / f"{prefix}.hotspot_summary.json",
    }
    group_rows.to_csv(paths["group_position_rows"], index=False)

    if group_rows.empty:
        hotspots = pd.DataFrame(columns=[
            "scaffold_group", "position", "group_position_edge_count", "template_support_sum",
            "unique_sidechain_transitions", "total_group_position_edge_count", "hotspot_coefficient",
            "rank_within_group"
        ])
        transitions = pd.DataFrame(columns=[
            "scaffold_group", "position", "source_sidechain", "target_sidechain", "edge_count", "template_support_sum"
        ])
    else:
        group_rows["transition"] = group_rows["source_sidechain"].astype(str) + "->" + group_rows["target_sidechain"].astype(str)
        gp = group_rows.groupby(["group", "position"], dropna=False)
        hotspots = gp.agg(
            group_position_edge_count=("position", "size"),
            template_support_sum=("template_count_total", "sum"),
            unique_sidechain_transitions=("transition", "nunique"),
            unique_source_molecules=("source_id", "nunique"),
            unique_target_molecules=("target_id", "nunique"),
            same_group_records=("same_group", "sum"),
        ).reset_index().rename(columns={"group": "scaffold_group"})
        total = hotspots.groupby("scaffold_group")["group_position_edge_count"].sum().rename("total_group_position_edge_count")
        hotspots = hotspots.merge(total, on="scaffold_group", how="left")
        hotspots["hotspot_coefficient"] = hotspots["group_position_edge_count"] / hotspots["total_group_position_edge_count"].replace(0, pd.NA)
        hotspots["rank_within_group"] = hotspots.groupby("scaffold_group")["group_position_edge_count"].rank(method="dense", ascending=False).astype(int)
        hotspots = hotspots.sort_values(["scaffold_group", "rank_within_group", "position"])

        transitions = group_rows.groupby(["group", "position", "source_sidechain", "target_sidechain"], dropna=False).agg(
            edge_count=("position", "size"),
            template_support_sum=("template_count_total", "sum"),
            unique_source_molecules=("source_id", "nunique"),
            unique_target_molecules=("target_id", "nunique"),
        ).reset_index().rename(columns={"group": "scaffold_group"}).sort_values(["scaffold_group", "position", "edge_count"], ascending=[True, True, False])

    hotspots.to_csv(paths["group_position_hotspots"], index=False)
    transitions.to_csv(paths["group_position_transitions"], index=False)
    qc = {
        "network_dir": str(network_dir),
        "changes_path": str(changes_path),
        "sidechain_change_rows": int(len(sidechain)),
        "group_position_rows": int(len(group_rows)),
        "hotspot_rows": int(len(hotspots)),
        "transition_rows": int(len(transitions)),
        "include_cross_edges": bool(include_cross_edges),
        "coefficient_definition": "group_position_edge_count / total_group_position_edge_count",
        "outputs": {k: str(v) for k, v in paths.items()},
    }
    write_json(qc, paths["summary"])
    return HotspotResult(paths=paths, qc=qc)
