from __future__ import annotations

from typing import Any

from tapocket.core.schema import CandidatePocket


def mean_mapping_coverage(candidates: list[CandidatePocket]) -> float:
    values = [float(c.mapping_quality.get("mapping_coverage", 0.0) or 0.0) for c in candidates]
    return sum(values) / len(values) if values else 0.0


def should_run_ai(config: Any, final_candidates: list[CandidatePocket], selected_template_count: int) -> tuple[bool, str]:
    if not bool(config.get("ai_models", "enabled", default=False)):
        return False, "ai_disabled"
    mode = str(config.get("ai_models", "mode", default="fallback_only"))
    if mode in {"always", "rerank", "merge"}:
        return True, mode

    trigger = config.get("ai_models", "trigger", default={}) or {}
    if mode == "fallback_only":
        if bool(trigger.get("no_template_hit", True)) and selected_template_count == 0:
            return True, "no_template_hit"
        if bool(trigger.get("no_final_pocket", True)) and len(final_candidates) == 0:
            return True, "no_final_pocket"
        min_count = int(trigger.get("min_final_pocket_count", 1))
        if len(final_candidates) < min_count:
            return True, "below_min_final_pocket_count"
        if bool(trigger.get("low_mapping_quality", False)):
            threshold = float(trigger.get("min_mean_mapping_coverage", 0.40))
            if mean_mapping_coverage(final_candidates) < threshold:
                return True, "low_mapping_quality"
    return False, "template_result_available"
