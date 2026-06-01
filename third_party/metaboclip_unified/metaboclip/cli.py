from __future__ import annotations

import argparse
import json
from pathlib import Path

from metaboclip.core.config import package_path, load_yaml
from metaboclip.core.workflow import run_directory, run_single_pair
from metaboclip.paper_locked.runner import PaperLockedRunner, prepare_original_scripts, verify_original_scripts, write_json_report


CURATED_FAMILIES = {"ugt", "act", "cyp450", "fe2og", "ach"}


def _default_profile() -> Path:
    return package_path("config", "profiles", "default_profile.yaml")


def _curated_mechanism(family: str) -> Path:
    fam = family.lower()
    if fam not in CURATED_FAMILIES:
        raise ValueError(f"Unknown curated family: {family}")
    return package_path("config", "families", fam, "mechanism.yaml")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="metaboclip", description="MetaboClip ligand-role and catalytic geometry scoring toolkit")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="Run generic YAML mechanism mode")
    run.add_argument("--mechanism", required=True, type=Path)
    run.add_argument("--profile", default=_default_profile(), type=Path)
    run.add_argument("--protein-dir", required=True, type=Path)
    run.add_argument("--docking-dir", required=True, type=Path)
    run.add_argument("--role-table-dir", required=True, type=Path)
    run.add_argument("--out-dir", required=True, type=Path)
    run.add_argument("--ligand-id", default=None)

    validate = sub.add_parser("validate-mechanism", help="Validate a mechanism YAML file before running")
    validate.add_argument("--mechanism", required=True, type=Path)
    validate.add_argument("--report", default=None, type=Path)

    curated = sub.add_parser("run-curated", help="Run built-in curated family mode")
    curated.add_argument("--family", required=True, choices=sorted(CURATED_FAMILIES))
    curated.add_argument("--profile", default=_default_profile(), type=Path)
    curated.add_argument("--protein-dir", required=True, type=Path)
    curated.add_argument("--docking-dir", required=True, type=Path)
    curated.add_argument("--role-table-dir", required=True, type=Path)
    curated.add_argument("--out-dir", required=True, type=Path)
    curated.add_argument("--ligand-id", default=None)

    pl_prepare = sub.add_parser("prepare-paper-locked", help="Prepare external original scripts from an archive")
    pl_prepare.add_argument("--archive", required=True, type=Path)
    pl_prepare.add_argument("--dest", required=True, type=Path)
    pl_prepare.add_argument("--overwrite", action="store_true")
    pl_prepare.add_argument("--report", default=None, type=Path)

    pl_verify = sub.add_parser("verify-paper-locked", help="Verify external original scripts against paper-locked checksums")
    pl_verify.add_argument("--original-root", required=True, type=Path)
    pl_verify.add_argument("--family", default=None)
    pl_verify.add_argument("--report", default=None, type=Path)

    pl_run = sub.add_parser("run-paper-locked", help="Run byte-for-byte external original scripts")
    pl_run.add_argument("--family", required=True, choices=sorted(CURATED_FAMILIES | {"p450", "2-odd", "2odd", "dac"}))
    pl_run.add_argument("--original-root", required=True, type=Path)
    pl_run.add_argument("--stage", choices=["gate", "score", "both"], default="both")
    pl_run.add_argument("--file-range", default=None)
    pl_run.add_argument("--docking-base-dir", default=None, type=Path)
    pl_run.add_argument("--out-dir", default=None, type=Path)
    pl_run.add_argument("--python", default=None)
    pl_run.add_argument("--skip-checksums", action="store_true")
    pl_run.add_argument("--dry-run", action="store_true")
    pl_run.add_argument("--report", default=None, type=Path)
    pl_run.add_argument("--extra-gate-arg", action="append", default=[])
    pl_run.add_argument("--extra-score-arg", action="append", default=[])

    export_pml = sub.add_parser("export-pymol", help="Export a PyMOL script for passing candidates")
    export_pml.add_argument("--protein", required=True, type=Path)
    export_pml.add_argument("--docked-pdbqt", required=True, type=Path)
    export_pml.add_argument("--resolved-ligand-sites", required=True, type=Path)
    export_pml.add_argument("--resolved-protein-roles", required=True, type=Path)
    export_pml.add_argument("--candidate-scores", default=None, type=Path)
    export_pml.add_argument("--passing-candidates", default=None, type=Path)
    export_pml.add_argument("--geometry-features", default=None, type=Path)
    export_pml.add_argument("--mechanism", default=None, type=Path)
    export_pml.add_argument("--out-pml", required=True, type=Path)
    export_pml.add_argument("--site-set-id", default=None)
    export_pml.add_argument("--pose-id", default=None)
    export_pml.add_argument("--select", choices=["best", "first"], default="best")
    export_pml.add_argument("--save-pse", default=None)


    export_batch = sub.add_parser("export-pymol-batch", help="Export one PyMOL script per passing pose or candidate")
    export_batch.add_argument("--protein", required=True, type=Path)
    export_batch.add_argument("--docked-pdbqt", required=True, type=Path)
    export_batch.add_argument("--resolved-ligand-sites", required=True, type=Path)
    export_batch.add_argument("--resolved-protein-roles", required=True, type=Path)
    export_batch.add_argument("--candidate-scores", required=True, type=Path)
    export_batch.add_argument("--passing-candidates", default=None, type=Path)
    export_batch.add_argument("--geometry-features", default=None, type=Path)
    export_batch.add_argument("--mechanism", default=None, type=Path)
    export_batch.add_argument("--out-dir", required=True, type=Path)
    export_batch.add_argument("--top-n", type=int, default=0)
    export_batch.add_argument("--mode", choices=["pose", "candidate"], default="pose")
    export_batch.add_argument("--no-pse", action="store_true")

    single = sub.add_parser("run-single", help="Run one protein-ligand pair")
    single.add_argument("--mechanism", required=True, type=Path)
    single.add_argument("--profile", default=_default_profile(), type=Path)
    single.add_argument("--protein", required=True, type=Path)
    single.add_argument("--docked-pdbqt", required=True, type=Path)
    single.add_argument("--role-table", required=True, type=Path)
    single.add_argument("--out-dir", required=True, type=Path)
    single.add_argument("--vina-out", default=None, type=Path)

    args = parser.parse_args(argv)
    if args.command == "validate-mechanism":
        from metaboclip.core.validation import validate_mechanism_file
        report = validate_mechanism_file(args.mechanism).to_dict()
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"Mechanism validation ok={report['ok']} errors={len(report['errors'])} warnings={len(report['warnings'])}")
        return 0 if report["ok"] else 2
    if args.command == "run":
        result = run_directory(args.mechanism, args.profile, args.protein_dir, args.docking_dir, args.role_table_dir, args.out_dir, args.ligand_id)
        print(f"Completed generic run: pairs={result['pairs']} candidates={result['candidates']}")
        return 0
    if args.command == "run-curated":
        mechanism = _curated_mechanism(args.family)
        result = run_directory(mechanism, args.profile, args.protein_dir, args.docking_dir, args.role_table_dir, args.out_dir, args.ligand_id)
        print(f"Completed curated run: family={args.family} pairs={result['pairs']} candidates={result['candidates']}")
        return 0
    if args.command == "prepare-paper-locked":
        report = prepare_original_scripts(args.archive, args.dest, overwrite=args.overwrite)
        write_json_report(report, args.report)
        print(f"Prepared paper-locked scripts under: {args.dest}")
        return 0
    if args.command == "verify-paper-locked":
        report = verify_original_scripts(args.original_root, args.family)
        write_json_report(report, args.report)
        print(f"Paper-locked verification ok={report['ok']}")
        return 0
    if args.command == "run-paper-locked":
        runner = PaperLockedRunner(args.original_root, args.family, python_executable=args.python, skip_checksums=args.skip_checksums)
        report = runner.run(stage=args.stage, file_range=args.file_range, docking_base_dir=args.docking_base_dir, out_dir=args.out_dir, dry_run=args.dry_run, extra_gate_args=args.extra_gate_arg, extra_score_args=args.extra_score_arg)
        write_json_report(report, args.report)
        print(f"Completed paper-locked run: family={report['family']} stage={args.stage} dry_run={args.dry_run}")
        return 0
    if args.command == "export-pymol":
        from metaboclip.pymol_exporter import export_pymol_script
        export_pymol_script(
            protein_file=args.protein,
            docked_ligand_file=args.docked_pdbqt,
            resolved_ligand_sites=args.resolved_ligand_sites,
            resolved_protein_roles=args.resolved_protein_roles,
            candidate_scores=args.candidate_scores,
            passing_candidates=args.passing_candidates,
            geometry_features=args.geometry_features,
            mechanism_file=args.mechanism,
            out_pml=args.out_pml,
            candidate_site_set_id=args.site_set_id,
            pose_id=args.pose_id,
            select_mode=args.select,
            save_pse_name=args.save_pse,
        )
        print(f"Wrote PyMOL script: {args.out_pml}")
        return 0
    if args.command == "export-pymol-batch":
        from metaboclip.pymol_exporter import export_pymol_batch_scripts
        report = export_pymol_batch_scripts(
            protein_file=args.protein,
            docked_ligand_file=args.docked_pdbqt,
            resolved_ligand_sites=args.resolved_ligand_sites,
            resolved_protein_roles=args.resolved_protein_roles,
            candidate_scores=args.candidate_scores,
            passing_candidates=args.passing_candidates,
            geometry_features=args.geometry_features,
            mechanism_file=args.mechanism,
            out_dir=args.out_dir,
            top_n=args.top_n if args.top_n > 0 else None,
            mode=args.mode,
            make_pse=not args.no_pse,
        )
        print(f"Wrote PyMOL exports: {len(report)} files under {args.out_dir}")
        return 0
    if args.command == "run-single":
        mechanism = load_yaml(args.mechanism)
        profile = load_yaml(args.profile)
        result = run_single_pair(mechanism, profile, args.protein, args.docked_pdbqt, args.role_table, args.out_dir, args.mechanism, args.vina_out)
        print(f"Completed single run: candidates={result['candidates']}")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
