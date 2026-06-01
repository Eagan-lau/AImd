from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np

from metaboclip.core.atoms import Atom


@dataclass(frozen=True)
class Plane:
    point: np.ndarray
    normal: np.ndarray


@dataclass(frozen=True)
class Axis:
    origin: np.ndarray
    direction: np.ndarray


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n <= 1e-12:
        raise ValueError("Zero-length vector")
    return v / n


def best_fit_plane(atoms: list[Atom]) -> Plane:
    if len(atoms) < 3:
        raise ValueError("At least three atoms are required for a best-fit plane")
    coords = np.array([a.coord for a in atoms], dtype=float)
    centroid = coords.mean(axis=0)
    centered = coords - centroid
    _, _, vh = np.linalg.svd(centered)
    normal = _unit(vh[-1])
    return Plane(point=centroid, normal=normal)


def oriented_axis(origin: Atom, direction: np.ndarray, away_from: Atom | None = None) -> Axis:
    vec = _unit(direction)
    if away_from is not None:
        away_vec = away_from.coord - origin.coord
        if float(np.dot(vec, away_vec)) > 0.0:
            vec = -vec
    return Axis(origin=origin.coord, direction=vec)


def axis_deviation(axis: Axis, start: Atom, end: Atom, fold180: bool = True) -> float:
    vec = _unit(end.coord - start.coord)
    dot = max(-1.0, min(1.0, float(np.dot(axis.direction, vec))))
    theta = float(math.degrees(math.acos(dot)))
    if fold180:
        theta = min(theta, 180.0 - theta)
    return theta
