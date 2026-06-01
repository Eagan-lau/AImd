#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
from pathlib import Path

from .bridge import run_metaboclip_bridge_core
from .utils import load_yaml, write_json


def run_metaboclip_bridge(config_path: str | Path) -> Path:
    config_path = Path(config_path).resolve()
    config = load_yaml(config_path)
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    print(f"[MetaBoClipBridge] AImd root: {root}")
    final_path = run_metaboclip_bridge_core(config)
    write_json(final_path.parent / "run_config_snapshot.json", config)
    return final_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run AImd MetaBoClipBridge with the unified MetaboClip backend")
    parser.add_argument("--config", required=True, help="Path to configs/MetaBoClip/metaboclip_bridge.yaml")
    args = parser.parse_args()
    run_metaboclip_bridge(args.config)


if __name__ == "__main__":
    main()
