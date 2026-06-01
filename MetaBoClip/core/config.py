from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if cfg is None:
        return {}
    if not isinstance(cfg, dict):
        raise RuntimeError("Config root must be a mapping")
    return cfg


def dump_config(cfg: dict, path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, sort_keys=False, allow_unicode=True)
    return out


def deep_merge(base: Any, override: Any) -> Any:
    if isinstance(base, dict) and isinstance(override, dict):
        merged = {k: deep_merge(base[k], override[k]) if k in base else override[k] for k in override}
        for k, v in base.items():
            if k not in merged:
                merged[k] = v
        return merged
    if isinstance(base, list) and isinstance(override, list):
        return override
    return override if override is not None else base


def normalize_family_key(value: str | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    return text.lower().replace(" ", "_").replace("-", "_")


def normalize_family_section(cfg: dict, fallback_key: str | None = None, fallback_name: str | None = None) -> dict:
    family = dict(cfg.get("family") or {})
    key = (
        family.get("plugin")
        or family.get("key")
        or family.get("runner")
        or family.get("name")
        or fallback_key
        or fallback_name
    )
    norm = normalize_family_key(key)
    if norm:
        family.setdefault("plugin", norm)
        family.setdefault("key", norm)
        family.setdefault("runner", norm)
    if fallback_name and "name" not in family:
        family["name"] = fallback_name
    cfg["family"] = family
    return cfg
