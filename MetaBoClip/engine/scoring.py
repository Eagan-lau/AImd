from __future__ import annotations

from typing import Any

import pandas as pd

from core.scoring import energy_linear, distance_piecewise, gaussian_score, coupled_mean_positive, aggregate_standard


SCORING_OPS = {
    "energy_linear",
    "distance_piecewise",
    "gaussian",
    "piecewise_linear",
    "coupled_mean_positive",
    "mean_positive",
    "weighted_sum",
}


def _piecewise_scalar(x: Any, points: list[list[float]]) -> float:
    if pd.isna(x):
        return 0.0
    xv = float(x)
    pts = [(float(a), float(b)) for a, b in points]
    pts.sort(key=lambda t: t[0])
    if xv <= pts[0][0]:
        return float(pts[0][1])
    if xv >= pts[-1][0]:
        return float(pts[-1][1])
    for i in range(len(pts) - 1):
        x1, y1 = pts[i]
        x2, y2 = pts[i + 1]
        if x1 <= xv <= x2:
            if x2 == x1:
                return float(y2)
            t = (xv - x1) / (x2 - x1)
            return float(y1 + t * (y2 - y1))
    return float(pts[-1][1])


def evaluate_scoring_op(df: pd.DataFrame, spec: dict[str, Any]) -> pd.Series:
    op = spec.get("op") or spec.get("type") or spec.get("expression")

    if op == "energy_linear":
        src = spec["source"]
        full = float(spec.get("full", spec.get("min", -7.0)))
        zero = float(spec.get("zero", spec.get("max", -3.0)))
        return df[src].apply(lambda x: energy_linear(x, full, zero))

    if op == "distance_piecewise":
        src = spec["source"]
        best = float(spec["best"])
        cutoff = float(spec["cutoff"])
        return df[src].apply(lambda x: distance_piecewise(x, best, cutoff))

    if op == "gaussian":
        src = spec["source"]
        target = float(spec["target"])
        sigma = float(spec["sigma"])
        flat = float(spec.get("flat", 0.0))
        fold180 = bool(spec.get("fold180", False))
        return df[src].apply(lambda x: gaussian_score(x, target, sigma, flat=flat, fold180=fold180))

    if op == "piecewise_linear":
        src = spec["source"]
        points = spec["points"]
        return df[src].apply(lambda x: _piecewise_scalar(x, points))

    if op in {"coupled_mean_positive", "mean_positive"}:
        sources = spec.get("sources") or spec.get("inputs")
        if not sources or len(sources) != 2:
            raise RuntimeError(f"{op} currently expects exactly two sources")
        a, b = sources
        return pd.Series([coupled_mean_positive(x, y) for x, y in zip(df[a], df[b])], index=df.index)

    if op == "weighted_sum":
        terms = spec.get("terms", [])
        scale = float(spec.get("scale", 1.0))
        total = pd.Series([0.0] * len(df), index=df.index, dtype=float)
        for term in terms:
            src = term["source"]
            weight = float(term.get("weight", 1.0))
            total = total + weight * df[src].astype(float)
        return total * scale

    raise RuntimeError(f"Unknown scoring operator: {op}")


def score_rows(rows: list[dict[str, Any]], cfg: dict[str, Any]) -> dict[str, pd.DataFrame]:
    df = pd.DataFrame(rows)
    if df.empty:
        return {
            "pose_scores": pd.DataFrame(),
            "conformation_scores": pd.DataFrame(),
            "protein_scores": pd.DataFrame(),
        }

    scoring_cfg = cfg.get("scoring", {})
    for block in scoring_cfg.get("transforms", []):
        df[block["name"]] = evaluate_scoring_op(df, block)

    for block in scoring_cfg.get("intermediates", []):
        df[block["name"]] = evaluate_scoring_op(df, block)

    pose_spec = scoring_cfg.get("pose_score", {"op": "weighted_sum", "terms": [], "scale": 100.0})
    if "scale" not in pose_spec:
        pose_spec = dict(pose_spec)
        pose_spec["scale"] = 100.0
    df["s_pose"] = evaluate_scoring_op(df, pose_spec)

    agg = cfg.get("aggregation", {})
    conf_df, prot_df = aggregate_standard(
        df,
        total_confs=int(agg.get("total_confs", 6)),
        cover_t=float(agg.get("cover_t", 70.0)),
        alpha=float(agg.get("alpha", 0.30)),
    )
    return {
        "pose_scores": df,
        "conformation_scores": conf_df,
        "protein_scores": prot_df,
    }
