#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from core.config import dump_config, load_config
from core.plugin_manager import PluginManager
from core.runtime import run_family, show_primitives, resolve_plugin_and_config


def _print_families(pm: PluginManager):
    for info in pm.describe_all():
        print(f"{info['key']:<12} {info['display_name']:<12} v{info['version']}  {info['description']}")


def main():
    ap = argparse.ArgumentParser(description="YAML-driven catalytic family runner")
    ap.add_argument("--instance", help="Path to family YAML config overlay")
    ap.add_argument("--family", help="Family plugin key, e.g. ugt / act / cyp450")
    ap.add_argument("--protein-dir", help="Directory containing protein files")
    ap.add_argument("--docking-dir", help="Single file_N directory")
    ap.add_argument("--out-dir", help="Output directory")
    ap.add_argument("--plugin-dir", action="append", default=[], help="Extra family-spec directory. Can be provided multiple times.")
    ap.add_argument("--list-families", action="store_true", help="List discovered family plugins")
    ap.add_argument("--show-family", help="Show resolved metadata for one family spec")
    ap.add_argument("--write-config", help="Write the default config for --family to this path")
    ap.add_argument("--dump-effective-config", help="Write the merged effective config used for this run")
    ap.add_argument("--show-primitives", action="store_true", help="Show registered primitives and family specs")
    args = ap.parse_args()

    pm = PluginManager(plugin_dirs=args.plugin_dir)

    if args.list_families:
        _print_families(pm)
        return

    if args.show_family:
        plugin = pm.get(args.show_family)
        print(json.dumps(plugin.describe(), indent=2, ensure_ascii=False))
        return

    if args.write_config:
        if not args.family:
            raise RuntimeError("--write-config requires --family")
        plugin = pm.get(args.family)
        cfg = plugin.load_default_config()
        out_path = dump_config(cfg, args.write_config)
        print(f"Default config written to: {out_path}")
        return

    if args.show_primitives:
        show_primitives(plugin_dirs=args.plugin_dir)
        return

    raw_cfg = load_config(args.instance) if args.instance else {}
    _, plugin, merged_cfg, _ = resolve_plugin_and_config(raw_cfg, plugin_dirs=args.plugin_dir, explicit_family=args.family)

    if args.dump_effective_config:
        out_path = dump_config(merged_cfg, args.dump_effective_config)
        print(f"Effective config written to: {out_path}")

    if not args.protein_dir or not args.docking_dir or not args.out_dir:
        raise RuntimeError("--protein-dir, --docking-dir and --out-dir are required for execution")

    result = run_family(
        cfg=merged_cfg,
        protein_dir=args.protein_dir,
        docking_dir=args.docking_dir,
        out_dir=args.out_dir,
        plugin_dirs=args.plugin_dir,
        explicit_family=plugin.key,
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
