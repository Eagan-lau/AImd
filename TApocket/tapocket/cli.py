from __future__ import annotations

import argparse
import json
import os
import subprocess
import shutil
from pathlib import Path

from tapocket.core.config import load_config
from tapocket.core.pipeline_template import TApocketPipeline
from tapocket.databases.manifest import check_database_layout, write_manifests
from tapocket.databases.template_db import TemplateDB
from tapocket.databases.mcsa_db import MCSADB
from tapocket.retrievers.foldseek import (
    build_foldseek_db_from_template_records,
    build_foldseek_db_from_mcsa_records,
)


DEFAULT_GVP_REPO_URL = "https://github.com/Mickdub/gvp.git"
DEFAULT_GVP_GH_REPO = "Mickdub/gvp"


def cmd_check_layout(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    result = check_database_layout(config)
    print(json.dumps(result, indent=2, ensure_ascii=False))


def cmd_build_manifest(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    result = write_manifests(config)
    print(json.dumps(result, indent=2, ensure_ascii=False))


def cmd_build_index(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    binary = config.get("foldseek", "binary", default="foldseek")
    results = []

    if args.db in {"template", "all"}:
        template_manifest = config.path("paths", "template_manifest")
        if not template_manifest.exists():
            print("Template manifest not found. Building manifests first...")
            write_manifests(config)
        template_db = TemplateDB.from_manifest(template_manifest, config.root)
        result = build_foldseek_db_from_template_records(
            records=template_db.records,
            root=config.root,
            staging_dir=config.path("paths", "staging_template_inputs"),
            output_db=config.path("paths", "foldseek_template_db"),
            binary=binary,
            create_index=args.create_index,
            force=args.force,
        )
        result["kind"] = "template"
        results.append(result)

    if args.db in {"mcsa", "all"}:
        mcsa_manifest = config.path("paths", "mcsa_manifest")
        if not mcsa_manifest.exists():
            print("M-CSA manifest not found. Building manifests first...")
            write_manifests(config)
        mcsa_db = MCSADB.from_manifest(mcsa_manifest, config.root)
        result = build_foldseek_db_from_mcsa_records(
            records=mcsa_db.records,
            root=config.root,
            staging_dir=config.path("paths", "staging_mcsa_inputs"),
            output_db=config.path("paths", "foldseek_mcsa_db"),
            binary=binary,
            create_index=args.create_index,
            force=args.force,
        )
        result["kind"] = "mcsa"
        results.append(result)

    print(json.dumps(results, indent=2, ensure_ascii=False))


def _command_exists(name: str) -> bool:
    try:
        subprocess.run(["bash", "-lc", f"command -v {name}"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except Exception:
        return False


def cmd_setup_pocketminer(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    dest_rel = Path(args.dest)
    dest = root / dest_rel

    default_repo_url = "https://github.com/Mickdub/gvp.git"
    default_gh_repo = "Mickdub/gvp"
    repo_url = args.repo_url or os.environ.get("GVP_POCKET_PRED_REPO_URL", default_repo_url)
    gh_repo = args.gh_repo or os.environ.get("GVP_POCKET_PRED_GH_REPO", default_gh_repo)
    method = args.method

    dest.parent.mkdir(parents=True, exist_ok=True)
    if (dest / ".git").exists():
        print(f"gvp-pocket_pred already exists: {dest}")
    elif dest.exists() and any(dest.iterdir()):
        raise SystemExit(f"Destination exists and is not empty: {dest}")
    else:
        if dest.exists():
            dest.rmdir()
        if method == "gh":
            if shutil.which("gh") is None:
                raise SystemExit("GitHub CLI 'gh' not found. Install gh or use --method git.")
            print(f"Cloning PocketMiner with GitHub CLI: gh repo clone {gh_repo} {dest}")
            subprocess.run(["gh", "repo", "clone", gh_repo, str(dest)], check=True, cwd=str(root))
        else:
            if shutil.which("git") is None:
                raise SystemExit("git not found. Install git first.")
            print(f"Cloning PocketMiner with git: git clone {repo_url} {dest}")
            subprocess.run(["git", "clone", repo_url, str(dest)], check=True, cwd=str(root))

    (dest / "models" / "pocketminer").mkdir(parents=True, exist_ok=True)
    result = {
        "status": "ok",
        "method": method,
        "repo_url": repo_url if method == "git" else None,
        "gh_repo": gh_repo if method == "gh" else None,
        "destination": str(dest),
        "pocketminer_root": str(dest_rel / "src"),
        "pocketminer_checkpoint": str(dest_rel / "models" / "pocketminer"),
        "examples": [
            "tapocket setup-pocketminer",
            "tapocket setup-pocketminer --method gh",
            "gh repo clone Mickdub/gvp third_party/gvp-pocket_pred",
        ],
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))

def cmd_run(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    pipeline = TApocketPipeline(config)
    summary = pipeline.run(query=args.query, run_id=args.run_id)
    print(json.dumps(summary.to_dict(), indent=2, ensure_ascii=False))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="TApocket template-first pocket pipeline with optional AI adapters")
    parser.add_argument("--version", action="version", version="TApocket full v4 0.6.0")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("check-layout", help="Check TemplateDB and M-CSA directory layout")
    p.add_argument("--config", default="configs/tapocket_template_v1.yaml")
    p.set_defaults(func=cmd_check_layout)

    p = sub.add_parser("build-manifest", help="Scan local databases and write manifest JSONL files")
    p.add_argument("--config", default="configs/tapocket_template_v1.yaml")
    p.set_defaults(func=cmd_build_manifest)

    p = sub.add_parser("build-index", help="Build Foldseek DBs from manifest records")
    p.add_argument("--config", default="configs/tapocket_template_v1.yaml")
    p.add_argument("--db", choices=["template", "mcsa", "all"], default="all")
    p.add_argument("--create-index", action="store_true", help="Run foldseek createindex after createdb")
    p.add_argument("--force", action="store_true", help="Clean staging folders before creating symlinks")
    p.set_defaults(func=cmd_build_index)

    p = sub.add_parser("setup-pocketminer", help="Clone PocketMiner/gvp into third_party/gvp-pocket_pred")
    p.add_argument("--repo-url", default=None, help="Git URL. Default: https://github.com/Mickdub/gvp.git; can also use GVP_POCKET_PRED_REPO_URL.")
    p.add_argument("--gh-repo", default=None, help="GitHub CLI repo name. Default: Mickdub/gvp; can also use GVP_POCKET_PRED_GH_REPO.")
    p.add_argument("--method", choices=["git", "gh"], default="git", help="Clone method. Default: git. Use gh for `gh repo clone Mickdub/gvp`.")
    p.add_argument("--root", default=".", help="TApocket project root; default is current directory")
    p.add_argument("--dest", default="third_party/gvp-pocket_pred", help="Destination relative to project root")
    p.set_defaults(func=cmd_setup_pocketminer)

    p = sub.add_parser("run", help="Run TApocket for one query PDB")
    p.add_argument("--config", default="configs/tapocket_template_v1.yaml")
    p.add_argument("--query", required=True, help="Input query protein PDB")
    p.add_argument("--run-id", default=None)
    p.set_defaults(func=cmd_run)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
