from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, Optional

import yaml


def load_rules(path: str | Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if "functional_groups" not in data:
        raise ValueError("Rules file must contain a functional_groups section")
    return data


def filter_rules(rule_data: Dict[str, Any], groups: Optional[Iterable[str]] = None) -> Dict[str, Any]:
    if not groups:
        return rule_data
    wanted = set(groups)
    out = dict(rule_data)
    out["functional_groups"] = {
        key: value for key, value in rule_data.get("functional_groups", {}).items() if key in wanted
    }
    return out
