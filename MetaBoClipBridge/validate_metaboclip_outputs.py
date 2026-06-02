#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate AImd MetaBoClipBridge outputs")
    ap.add_argument("--out-dir", default="data/data_output/metaboclip/results")
    ap.add_argument("--require-scored", action="store_true", help="Require scored unified aggregate outputs, not only a dry-run report")
    args = ap.parse_args()
    out = Path(args.out_dir)
    required = [
        out / "metaboclip_run_manifest.csv",
        out / "metaboclip_final_ranking.csv",
        out / "metaboclip_report.json",
    ]
    if args.require_scored:
        required.extend([
            out / "metaboclip_protein_scores_all.csv",
            out / "metaboclip_conformation_scores_all.csv",
            out / "metaboclip_pose_scores_all.csv",
            out / "metaboclip_candidate_scores_all.csv",
            out / "metaboclip_geometry_features_all.csv",
            out / "metaboclip_resolved_ligand_sites_all.csv",
            out / "metaboclip_resolved_protein_roles_all.csv",
        ])
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing MetaBoClipBridge outputs:\n" + "\n".join(missing))
    print(f"MetaBoClipBridge outputs look complete: {out}")

if __name__ == "__main__":
    main()
