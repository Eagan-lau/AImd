from __future__ import annotations

import math
from typing import Any


def clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def energy_linear(value: float | None, full: float = -7.0, zero: float = -3.0) -> float:
    if value is None:
        return 0.0
    x = float(value)
    if x <= full:
        return 1.0
    if x >= zero:
        return 0.0
    return clamp((zero - x) / (zero - full))


def distance_piecewise(value: float | None, best: float, cutoff: float) -> float:
    if value is None:
        return 0.0
    x = float(value)
    if x <= best:
        return 1.0
    if x >= cutoff:
        return 0.0
    return clamp((cutoff - x) / (cutoff - best))


def angle_gaussian(value: float | None, target: float, sigma: float, flat: float = 0.0, fold180: bool = False) -> float:
    if value is None:
        return 0.0
    x = float(value)
    d = abs(x - target)
    if fold180:
        d = min(d, 180.0 - d)
    if d <= flat:
        return 1.0
    d = d - flat
    return float(math.exp(-0.5 * (d / sigma) ** 2))


def angle_deviation_linear(value: float | None, cutoff: float) -> float:
    if value is None:
        return 0.0
    x = abs(float(value))
    if x >= cutoff:
        return 0.0
    return clamp((cutoff - x) / cutoff)


def transform_score(transform: dict[str, Any], features: dict[str, float]) -> float:
    source = transform.get("source")
    value = features.get(source) if source else None
    typ = str(transform.get("transform") or transform.get("type") or "distance_piecewise")
    if typ == "energy_linear":
        return energy_linear(value, float(transform.get("full", -7.0)), float(transform.get("zero", -3.0)))
    if typ == "distance_piecewise":
        return distance_piecewise(value, float(transform.get("best", 3.2)), float(transform.get("cutoff", 6.0)))
    if typ == "angle_gaussian":
        return angle_gaussian(value, float(transform.get("target", 107.0)), float(transform.get("sigma", 20.0)), float(transform.get("flat", 0.0)), bool(transform.get("fold180", False)))
    if typ == "angle_deviation_linear":
        return angle_deviation_linear(value, float(transform.get("cutoff", 45.0)))
    raise ValueError(f"Unknown scoring transform: {typ}")


def weighted_mean(values: list[tuple[float, float]]) -> float | None:
    filtered = [(v, w) for v, w in values if v is not None and w > 0]
    if not filtered:
        return None
    total = sum(w for _, w in filtered)
    if total <= 0:
        return None
    return sum(v * w for v, w in filtered) / total


def hierarchical_geometry_score(level_scores: dict[int, float], level_weights: dict[int, float]) -> float:
    items = [(score, float(level_weights.get(level, 0.0))) for level, score in level_scores.items() if score is not None and float(level_weights.get(level, 0.0)) > 0]
    if not items:
        return 0.0
    total = sum(w for _, w in items)
    return sum(score * w for score, w in items) / total
