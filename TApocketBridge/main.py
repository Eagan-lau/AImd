#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
from pathlib import Path

from .runner import run_tapocket_batch


def main() -> None:
    parser = argparse.ArgumentParser(description="AImd/TApocketBridge: batch TApocket pocket prediction from RGPC protein_manifest.csv")
    parser.add_argument("--config", required=True, help="Path to configs/TApocket/tapocket_batch.yaml")
    args = parser.parse_args()
    run_tapocket_batch(Path(args.config))


if __name__ == "__main__":
    main()
