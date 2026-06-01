#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import argparse
import os
from pathlib import Path

def main() -> None:
    parser = argparse.ArgumentParser(description="Create or update a symlink in AImd/third_party/bin.")
    parser.add_argument("tool", help="Tool name, for example hipmcl")
    parser.add_argument("target", help="Absolute path to the executable")
    parser.add_argument("--root", default=".", help="AImd root directory")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    target = Path(args.target).expanduser().resolve()
    if not target.exists():
        raise FileNotFoundError(f"Target does not exist: {target}")
    bin_dir = root / "third_party" / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    link = bin_dir / args.tool
    if link.exists() or link.is_symlink():
        link.unlink()
    os.symlink(target, link)
    print(f"Created symlink: {link} -> {target}")
if __name__ == "__main__":
    main()
