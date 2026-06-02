#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from .runner import run_mollink_from_config


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run MolLink as the AImd ligand input transformation step")
    parser.add_argument("--config", default="configs/MolLink/mollink.yaml", help="Path to the MolLink bridge config")
    parser.add_argument("--mode", choices=["auto", "csv_only", "template"], default=None, help="Override configured run mode")
    parser.add_argument("--max-molecules", type=int, default=None, help="Limit molecule rows for a small smoke test")
    parser.add_argument("--max-templates", type=int, default=None, help="Limit reaction templates for a small smoke test")
    args = parser.parse_args(argv)

    result = run_mollink_from_config(
        Path(args.config),
        mode=args.mode,
        max_molecules=args.max_molecules,
        max_templates=args.max_templates,
    )
    print("[MolLink] run completed")
    print(f"[MolLink] mode: {result.summary.get('mode', '')}")
    print(f"[MolLink] molecules: {result.summary.get('valid_molecules', 0)} valid / {result.summary.get('invalid_molecules', 0)} invalid")
    print(f"[MolLink] directed edges: {result.summary.get('directed_edges', 0)}")
    print(f"[MolLink] ligand source manifest: {result.paths.get('ligand_source_manifest')}")
    print(f"[MolLink] summary: {result.paths.get('run_summary')}")


if __name__ == "__main__":
    main()
