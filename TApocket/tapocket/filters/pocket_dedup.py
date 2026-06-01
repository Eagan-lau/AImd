from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from tapocket.core.schema import CandidatePocket


def residue_set(candidate: CandidatePocket) -> set[tuple[str, str, str, str]]:
    return {r.key() for r in candidate.query_pocket_residues}


def residue_overlap(a: CandidatePocket, b: CandidatePocket) -> float:
    ra = residue_set(a)
    rb = residue_set(b)
    if not ra or not rb:
        return 0.0
    return len(ra & rb) / min(len(ra), len(rb))


def _candidate_sort_key(c: CandidatePocket) -> tuple[float, float, int]:
    coverage = float(c.mapping_quality.get("mapping_coverage", 0.0) or 0.0)
    return (float(c.final_score or 0.0), coverage, len(c.query_pocket_residues))


def deduplicate_pockets(
    candidates: list[CandidatePocket],
    overlap_threshold: float = 0.70,
    representative: str = "best_score",
) -> tuple[list[CandidatePocket], dict[str, Any]]:
    """Greedy pocket-level deduplication using query residue overlap.

    Two pockets are considered redundant if
    |residues_A ∩ residues_B| / min(|A|, |B|) >= overlap_threshold.
    The representative keeps the supporting evidence from all cluster members.
    """
    if not candidates:
        return [], {"enabled": True, "cluster_count": 0, "clusters": []}

    sorted_candidates = sorted(candidates, key=_candidate_sort_key, reverse=True)
    used: set[int] = set()
    representatives: list[CandidatePocket] = []
    clusters: list[dict[str, Any]] = []

    for i, seed in enumerate(sorted_candidates):
        if i in used:
            continue
        used.add(i)
        members: list[CandidatePocket] = [seed]
        for j in range(i + 1, len(sorted_candidates)):
            if j in used:
                continue
            other = sorted_candidates[j]
            ov = residue_overlap(seed, other)
            if ov >= overlap_threshold:
                used.add(j)
                members.append(other)

        if representative == "largest_residue_count":
            rep = max(members, key=lambda c: (len(c.query_pocket_residues),) + _candidate_sort_key(c))
        else:
            rep = max(members, key=_candidate_sort_key)

        cluster_id = f"pocket_cluster_{len(representatives) + 1:03d}"
        rep.cluster_id = cluster_id
        rep.custom_pocket_id = cluster_id
        rep.supporting_pockets = [
            {
                "pocket_id": m.pocket_id,
                "template_id": m.template_id,
                "template_pocket_id": m.template_pocket_id,
                "source": m.source,
                "score": m.final_score,
                "mapping_coverage": m.mapping_quality.get("mapping_coverage"),
                "query_residue_count": len(m.query_pocket_residues),
            }
            for m in members
        ]
        rep.mapping_quality["deduplicated"] = True
        rep.mapping_quality["supporting_pocket_count"] = len(members)
        rep.mapping_quality["supporting_template_count"] = len({m.template_id for m in members})
        representatives.append(rep)
        clusters.append(
            {
                "cluster_id": cluster_id,
                "representative_original_pocket_id": f"{rep.template_id}:{rep.template_pocket_id}",
                "supporting_pocket_count": len(members),
                "supporting_template_count": len({m.template_id for m in members}),
                "members": [f"{m.template_id}:{m.template_pocket_id}" for m in members],
            }
        )

    report = {
        "enabled": True,
        "method": "query_residue_overlap",
        "overlap_threshold": overlap_threshold,
        "input_pocket_count": len(candidates),
        "cluster_count": len(representatives),
        "clusters": clusters,
    }
    return representatives, report
