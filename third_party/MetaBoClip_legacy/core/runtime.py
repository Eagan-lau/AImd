from __future__ import annotations

from pathlib import Path

from core.config import deep_merge, normalize_family_key, normalize_family_section
from core.plugin_manager import PluginManager
from core.registry import build_default_registry
from engine.operators import OBJECT_OPS, FEATURE_OPS
from engine.scoring import SCORING_OPS


def _family_from_cfg(cfg: dict | None) -> str | None:
    family = dict((cfg or {}).get("family") or {})
    return normalize_family_key(
        family.get("plugin") or family.get("key") or family.get("runner") or family.get("name")
    )


def resolve_plugin_and_config(raw_cfg: dict | None, plugin_dirs=None, explicit_family: str | None = None):
    pm = PluginManager(plugin_dirs=plugin_dirs)
    family_key = normalize_family_key(explicit_family) or _family_from_cfg(raw_cfg)
    if not family_key:
        raise RuntimeError("Could not determine family plugin. Provide --family or set family.plugin/key/runner/name in the config.")
    plugin = pm.get(family_key)
    default_cfg = plugin.load_default_config()
    merged_cfg = deep_merge(default_cfg, raw_cfg or {})
    merged_cfg = normalize_family_section(merged_cfg, fallback_key=plugin.key, fallback_name=plugin.display_name)
    reg = build_default_registry()
    return pm, plugin, merged_cfg, reg


def show_primitives(plugin_dirs=None):
    print(f"object_ops: {sorted(OBJECT_OPS)}")
    print(f"feature_ops: {sorted(FEATURE_OPS)}")
    print(f"scoring_ops: {sorted(SCORING_OPS)}")
    pm = PluginManager(plugin_dirs=plugin_dirs)
    print(f"family_specs: {sorted(pm.plugins.keys())}")


def run_family(cfg, protein_dir, docking_dir, out_dir, plugin_dirs=None, explicit_family: str | None = None):
    _, plugin, merged_cfg, reg = resolve_plugin_and_config(cfg, plugin_dirs=plugin_dirs, explicit_family=explicit_family)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    return plugin.execute(
        cfg=merged_cfg,
        protein_dir=Path(protein_dir),
        docking_dir=Path(docking_dir),
        out_dir=out_dir,
        registry=reg,
    )
