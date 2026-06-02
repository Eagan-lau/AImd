from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import yaml

from .compute import compute_transformations
from .hotspot import compute_group_position_hotspots
from .io import inspect_molecule_library
from .network import build_network_from_computed
from .templates import validate_templates
from .visualize import visualize_existing_network


def _load_config(path: str | None) -> dict:
    if not path:
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _cfg_get(cfg: dict, *keys, default=None):
    cur: Any = cfg
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def _profile_defaults(profile: str) -> dict[str, str]:
    if profile == "strict":
        return {"match_method": "full_inchikey", "self_loop_policy": "drop"}
    if profile == "relaxed":
        return {"match_method": "connectivity_inchikey", "self_loop_policy": "drop"}
    if profile == "original":
        return {"match_method": "connectivity_inchikey", "self_loop_policy": "keep"}
    raise ValueError(f"Unsupported profile: {profile}")




def _print_invalid_smiles_records(result, max_records: int = 20) -> None:
    path = result.paths.get("invalid_molecules") if hasattr(result, "paths") else None
    count = int(result.qc.get("invalid_molecules", 0)) if hasattr(result, "qc") else 0
    if not path or count <= 0:
        return
    print(f"Invalid SMILES records: {count}")
    print(f"invalid_molecules: {path}")
    try:
        with open(path, "r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            for i, row in enumerate(reader):
                if i >= max_records:
                    remaining = count - max_records
                    if remaining > 0:
                        print(f"  ... {remaining} additional invalid records omitted from console output")
                    break
                source_id = row.get("source_id") or row.get("molecule_id") or row.get("row_index") or ""
                reason = row.get("drop_reason") or row.get("status") or "invalid"
                smiles = row.get("smiles_raw") or row.get("smiles") or ""
                error = row.get("error") or ""
                print(f"  source_id={source_id}; reason={reason}; smiles={smiles}; error={error}")
    except Exception as exc:
        print(f"  could not print invalid records: {exc}")

def _add_molecule_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--molecule-table", default=None, help="Molecule library table: CSV/TSV/XLSX/SDF/SMI.")
    parser.add_argument("--smiles-column", default=None, help="SMILES column; auto-detected when omitted.")
    parser.add_argument("--id-column", default=None, help="Molecule ID column; auto-detected when omitted.")
    parser.add_argument("--name-column", default=None, help="Molecule name column; auto-detected when omitted.")
    parser.add_argument("--sheet-name", default=None, help="Excel sheet name or index.")


def _add_compute_args(parser: argparse.ArgumentParser, output_required: bool = True) -> None:
    _add_molecule_args(parser)
    parser.add_argument("--templates", default=None, help="Reaction SMARTS template file: TXT/CSV/TSV/XLSX.")
    parser.add_argument("--template-column", default=None, help="Reaction SMARTS column for tabular template files.")
    parser.add_argument("--output-dir", required=output_required)
    parser.add_argument("--prefix", default="transformnet")
    parser.add_argument("--profile", default=None, choices=["strict", "relaxed", "original"], help="Preset matching profile.")
    parser.add_argument("--match-method", default=None, choices=["full_inchikey", "connectivity_inchikey", "canonical_smiles", "canonical_smiles_no_stereo"])
    parser.add_argument("--duplicate-key-policy", default=None, choices=["all", "first", "drop", "error"])
    parser.add_argument("--self-loop-policy", default=None, choices=["drop", "keep"])
    parser.add_argument("--max-molecules", type=int, default=None)
    parser.add_argument("--max-templates", type=int, default=None)
    parser.add_argument("--allow-multireactant", action="store_true", help="Attempt templates with more than one reactant.")
    parser.add_argument("--no-prefilter", action="store_true", help="Disable reactant-template substructure prefilter.")
    parser.add_argument("--largest-fragment", action="store_true")
    parser.add_argument("--uncharge", action="store_true")
    parser.add_argument("--remove-hs", action="store_true")
    parser.add_argument("--no-cleanup", action="store_true")
    parser.add_argument("--show-rdkit-logs", action="store_true")
    parser.add_argument("--write-unmatched-products", action="store_true")
    parser.add_argument("--no-product-failures", action="store_true")


def _add_build_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", default=None)
    parser.add_argument("--computed-dir", default=None)
    parser.add_argument("--valid-molecules", default=None)
    parser.add_argument("--directed-edges", default=None)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--annotation-table", default=None)
    parser.add_argument("--annotation-sheet", default=None)
    parser.add_argument("--annotation-id-column", default=None)
    parser.add_argument("--group-column", default=None)
    parser.add_argument("--position-columns", default=None, help="Comma-separated position columns to compare, e.g. R1,R2,R3.")
    parser.add_argument("--prefix", default="transformnet_network")
    parser.add_argument("--no-excel", action="store_true")
    parser.add_argument("--no-graphs", action="store_true")
    parser.add_argument("--no-html", action="store_true")


def _add_visualize_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--network-dir", required=True, type=Path, help="Directory containing TransformNet network CSV files.")
    parser.add_argument("--strict-network-dir", default=None, type=Path, help="Optional strict network directory for a dual-layer HTML view.")
    parser.add_argument("--relaxed-network-dir", default=None, type=Path, help="Optional relaxed network directory for a dual-layer HTML view.")
    parser.add_argument("--scaffold-json", required=True, type=Path, help="Path to taxane_lib_v2.json or a compatible scaffold JSON file.")
    parser.add_argument("--output-html", required=True, type=Path, help="Output interactive HTML path.")
    parser.add_argument("--valid-molecules", default=None, type=Path, help="Optional matrix-index to source-id mapping CSV.")
    parser.add_argument("--json-index-base", default=1, type=int, choices=[0, 1], help="Fallback ID base for scaffold JSON row mapping.")
    parser.add_argument("--edge-id-space", default="auto", choices=["auto", "source_id", "node_id", "matrix_index", "row_index", "label"], help="How edge endpoint values should be resolved to node IDs.")
    parser.add_argument("--max-nodes", default=2500, type=int, help="Maximum nodes included in the HTML view.")
    parser.add_argument("--max-edges", default=12000, type=int, help="Maximum edges included in the HTML view.")
    parser.add_argument("--layout-seed", default=17, type=int, help="Deterministic layout seed.")
    parser.add_argument("--central-groups", default=5, type=int, help="Number of largest scaffold groups placed near the center.")
    parser.add_argument("--label-column", default="source_id", help="Node column used for the visible node label.")
    parser.add_argument("--title", default="MolSpaceHub TransformNet interactive view", help="HTML page title.")
    parser.add_argument("--scaffold-image", default=None, type=Path, help="Optional scaffold image path. If omitted, network-dir/taxane_scaffold.png is used when present.")
    parser.add_argument("--no-molecule-svg", action="store_true", help="Disable optional RDKit molecule SVG depictions.")
    parser.add_argument("--node-size-scale", default=1.35, type=float, help="Scale factor for node radius in the interactive HTML view.")
    parser.add_argument("--label-size", default=10, type=int, help="Node label font size in pixels for the interactive HTML view.")


def _position_columns(value: str | None) -> list[str] | None:
    if not value:
        return None
    return [x.strip() for x in value.split(",") if x.strip()]


def _compute_from_args(args: argparse.Namespace, cfg: dict, output_dir: str | Path | None = None, profile: str | None = None):
    molecule_table = args.molecule_table or _cfg_get(cfg, "inputs", "molecule_table")
    templates = args.templates or _cfg_get(cfg, "inputs", "template_file")
    if not molecule_table:
        raise ValueError("--molecule-table is required unless provided in config.")
    if not templates:
        raise ValueError("--templates is required unless provided in config.")
    cols = cfg.get("columns", {})
    compute_cfg = cfg.get("compute", {})
    standardization = cfg.get("standardization", {})
    profile_name = profile or args.profile or compute_cfg.get("profile") or "strict"
    defaults = _profile_defaults(profile_name)
    return compute_transformations(
        molecule_table=molecule_table,
        template_path=templates,
        output_dir=output_dir or args.output_dir,
        smiles_column=args.smiles_column or cols.get("smiles"),
        id_column=args.id_column or cols.get("id"),
        name_column=args.name_column or cols.get("name"),
        sheet_name=args.sheet_name or cols.get("sheet_name"),
        template_column=args.template_column or cols.get("template"),
        match_method=args.match_method or compute_cfg.get("match_method") or defaults["match_method"],
        duplicate_key_policy=args.duplicate_key_policy or compute_cfg.get("duplicate_key_policy", "all"),
        self_loop_policy=args.self_loop_policy or compute_cfg.get("self_loop_policy") or defaults["self_loop_policy"],
        single_reactant_only=not args.allow_multireactant and compute_cfg.get("single_reactant_only", True),
        prefilter_substructure=not args.no_prefilter and compute_cfg.get("prefilter_substructure", True),
        max_molecules=args.max_molecules if args.max_molecules is not None else compute_cfg.get("max_molecules"),
        max_templates=args.max_templates if args.max_templates is not None else compute_cfg.get("max_templates"),
        standardize_cleanup=not args.no_cleanup and standardization.get("cleanup", True),
        largest_fragment=args.largest_fragment or standardization.get("largest_fragment", False),
        uncharge=args.uncharge or standardization.get("uncharge", False),
        remove_hs=args.remove_hs or standardization.get("remove_hs", False),
        quiet_rdkit=not args.show_rdkit_logs,
        write_raw_product_failures=not args.no_product_failures,
        write_unmatched_products=args.write_unmatched_products or compute_cfg.get("write_unmatched_products", False),
        prefix=args.prefix,
    )


def _build_network(args: argparse.Namespace, cfg: dict, computed_dir: str | Path | None = None, output_dir: str | Path | None = None):
    network_cfg = cfg.get("network", {})
    return build_network_from_computed(
        computed_dir=computed_dir or args.computed_dir or _cfg_get(cfg, "outputs", "computed_dir"),
        valid_molecules=args.valid_molecules,
        directed_edges=args.directed_edges,
        output_dir=output_dir or args.output_dir,
        annotation_table=args.annotation_table or _cfg_get(cfg, "inputs", "annotation_table"),
        annotation_sheet=args.annotation_sheet or _cfg_get(cfg, "columns", "annotation_sheet"),
        annotation_id_column=args.annotation_id_column or _cfg_get(cfg, "columns", "annotation_id"),
        group_column=args.group_column or network_cfg.get("group_column"),
        position_columns=_position_columns(args.position_columns) or network_cfg.get("position_columns"),
        prefix=args.prefix,
        write_excel=not args.no_excel,
        write_graphs=not args.no_graphs,
        write_html=not args.no_html,
    )


def _visualize(args: argparse.Namespace):
    return visualize_existing_network(
        network_dir=args.network_dir,
        scaffold_json=args.scaffold_json,
        output_html=args.output_html,
        valid_molecules=args.valid_molecules,
        json_index_base=args.json_index_base,
        edge_id_space=args.edge_id_space,
        max_nodes=args.max_nodes,
        max_edges=args.max_edges,
        layout_seed=args.layout_seed,
        central_groups=args.central_groups,
        label_column=args.label_column,
        title=args.title,
        make_svg=not args.no_molecule_svg,
        scaffold_image=args.scaffold_image,
        strict_network_dir=args.strict_network_dir,
        relaxed_network_dir=args.relaxed_network_dir,
        node_size_scale=args.node_size_scale,
        label_size=args.label_size,
    )


def _run_one_profile(args: argparse.Namespace, cfg: dict, profile: str, base_output: Path) -> dict[str, Any]:
    profile_dir = base_output if args.run_profile != "both" else base_output / profile
    compute_dir = profile_dir / "compute"
    network_dir = profile_dir / "network"
    compute_dir.mkdir(parents=True, exist_ok=True)
    network_dir.mkdir(parents=True, exist_ok=True)
    comp_res = _compute_from_args(args, cfg, output_dir=compute_dir, profile=profile)
    _print_invalid_smiles_records(comp_res)
    net_args = argparse.Namespace(
        computed_dir=str(compute_dir),
        valid_molecules=None,
        directed_edges=None,
        output_dir=str(network_dir),
        annotation_table=None,
        annotation_sheet=None,
        annotation_id_column=None,
        group_column=None,
        position_columns=None,
        prefix="transformnet_network",
        no_excel=False,
        no_graphs=False,
        no_html=False,
    )
    net_res = _build_network(net_args, cfg, computed_dir=compute_dir, output_dir=network_dir)
    html_path = network_dir / "transformnet_network.interactive.html"
    viz_res = visualize_existing_network(
        network_dir=network_dir,
        scaffold_json=args.scaffold_json,
        output_html=html_path,
        valid_molecules=compute_dir / "transformnet.valid_molecules.csv",
        json_index_base=args.json_index_base,
        edge_id_space=args.edge_id_space,
        max_nodes=args.max_nodes,
        max_edges=args.max_edges,
        layout_seed=args.layout_seed,
        central_groups=args.central_groups,
        label_column=args.label_column,
        title=args.title,
        make_svg=not args.no_molecule_svg,
        scaffold_image=args.scaffold_image,
        node_size_scale=args.node_size_scale,
        label_size=args.label_size,
    )
    hs_res = compute_group_position_hotspots(network_dir=network_dir, output_dir=network_dir, prefix="transformnet_network", include_cross_edges=not args.exclude_cross_edges_from_hotspots)
    return {
        "profile": profile,
        "compute_summary": str(comp_res.paths["summary"]),
        "network_summary": str(net_res.paths["summary"]),
        "visual_summary": viz_res.get("summary", ""),
        "hotspot_summary": str(hs_res.paths["summary"]),
        "interactive_html": str(html_path),
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="molspace-transformnet",
        description="Integrated molecular transformation-network construction, scaffold-aware visualization, and R-position hotspot analysis module.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    inspect = sub.add_parser("inspect-library", help="Inspect molecule table column detection and basic row counts.")
    inspect.add_argument("--config", default=None)
    _add_molecule_args(inspect)
    inspect.add_argument("--output", default=None, help="Optional JSON output path.")

    val = sub.add_parser("validate-templates", help="Validate reaction SMARTS templates.")
    val.add_argument("--templates", required=True)
    val.add_argument("--template-column", default=None)
    val.add_argument("--max-templates", type=int, default=1000)
    val.add_argument("--single-reactant-only", action="store_true")
    val.add_argument("--output", required=True)

    comp = sub.add_parser("compute", help="Compute directed transformation edges from molecule library and reaction templates.")
    comp.add_argument("--config", default=None)
    _add_compute_args(comp)

    build = sub.add_parser("build", help="Build network tables from computed edge files.")
    _add_build_args(build)

    viz = sub.add_parser("visualize", help="Create interactive scaffold-aware network HTML from existing network outputs.")
    _add_visualize_args(viz)

    hot = sub.add_parser("hotspot", help="Compute scaffold-group R-position hotspot coefficients from visualized side-chain changes.")
    hot.add_argument("--network-dir", required=True, type=Path)
    hot.add_argument("--output-dir", default=None, type=Path)
    hot.add_argument("--prefix", default="transformnet_network")
    hot.add_argument("--exclude-cross-edges", action="store_true")

    run_all = sub.add_parser("run-all", help="End-to-end compute, build, visualize, and hotspot analysis.")
    run_all.add_argument("--config", default=None)
    _add_compute_args(run_all, output_required=False)
    run_all.add_argument("--output-root", required=True, type=Path, help="Root output directory for all stages.")
    run_all.add_argument("--scaffold-json", required=True, type=Path)
    run_all.add_argument("--json-index-base", default=1, type=int, choices=[0, 1])
    run_all.add_argument("--edge-id-space", default="auto", choices=["auto", "source_id", "node_id", "matrix_index", "row_index", "label"])
    run_all.add_argument("--max-nodes", default=2500, type=int)
    run_all.add_argument("--max-edges", default=12000, type=int)
    run_all.add_argument("--layout-seed", default=17, type=int)
    run_all.add_argument("--central-groups", default=5, type=int)
    run_all.add_argument("--label-column", default="source_id")
    run_all.add_argument("--title", default="MolSpaceHub TransformNet interactive view")
    run_all.add_argument("--scaffold-image", default=None, type=Path)
    run_all.add_argument("--no-molecule-svg", action="store_true")
    run_all.add_argument("--node-size-scale", default=1.35, type=float, help="Scale factor for node radius in the interactive HTML view.")
    run_all.add_argument("--label-size", default=10, type=int, help="Node label font size in pixels for the interactive HTML view.")
    run_all.add_argument("--exclude-cross-edges-from-hotspots", action="store_true")
    run_all.add_argument("--run-profile", default="strict", choices=["strict", "relaxed", "original", "both"], help="End-to-end profile. Use both to run strict and relaxed profiles.")

    compat = sub.add_parser("compute-linkage", help="Alias for compute; writes linkage_new.csv and from_to_list.csv as compatibility outputs.")
    compat.add_argument("--config", default=None)
    _add_compute_args(compat)

    args = parser.parse_args(argv)

    if args.command == "inspect-library":
        cfg = _load_config(args.config)
        cols = cfg.get("columns", {})
        molecule_table = args.molecule_table or _cfg_get(cfg, "inputs", "molecule_table")
        if not molecule_table:
            raise ValueError("--molecule-table is required unless provided in config.")
        report = inspect_molecule_library(
            molecule_table,
            smiles_column=args.smiles_column or cols.get("smiles"),
            id_column=args.id_column or cols.get("id"),
            name_column=args.name_column or cols.get("name"),
            sheet_name=args.sheet_name or cols.get("sheet_name"),
        )
        if args.output:
            Path(args.output).parent.mkdir(parents=True, exist_ok=True)
            Path(args.output).write_text(json.dumps(report, ensure_ascii=True, indent=2), encoding="utf-8")
        else:
            print(json.dumps(report, ensure_ascii=True, indent=2))
        return

    if args.command == "validate-templates":
        df = validate_templates(args.templates, template_column=args.template_column, max_templates=args.max_templates, single_reactant_only=args.single_reactant_only)
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.suffix.lower() in {".xlsx", ".xls"}:
            df.to_excel(out, index=False)
        else:
            df.to_csv(out, index=False)
        print(f"Wrote template validation table to {out}")
        return

    if args.command in {"compute", "compute-linkage"}:
        cfg = _load_config(args.config)
        res = _compute_from_args(args, cfg)
        print("Transformation computation complete.")
        _print_invalid_smiles_records(res)
        for k, v in res.paths.items():
            print(f"{k}: {v}")
        print("QC:")
        for k, v in res.qc.items():
            if k != "outputs":
                print(f"  {k}: {v}")
        return

    if args.command == "build":
        cfg = _load_config(args.config)
        res = _build_network(args, cfg)
        print("Network build complete.")
        for k, v in res.paths.items():
            print(f"{k}: {v}")
        print("QC:")
        for k, v in res.qc.items():
            if k != "outputs":
                print(f"  {k}: {v}")
        return

    if args.command == "visualize":
        res = _visualize(args)
        print("Interactive TransformNet visualization complete.")
        for k, v in res.items():
            print(f"{k}: {v}")
        return

    if args.command == "hotspot":
        res = compute_group_position_hotspots(
            network_dir=args.network_dir,
            output_dir=args.output_dir,
            prefix=args.prefix,
            include_cross_edges=not args.exclude_cross_edges,
        )
        print("Group-position hotspot analysis complete.")
        for k, v in res.paths.items():
            print(f"{k}: {v}")
        return

    if args.command == "run-all":
        cfg = _load_config(args.config)
        profiles = ["strict", "relaxed"] if args.run_profile == "both" else [args.run_profile]
        summaries = []
        for profile in profiles:
            summaries.append(_run_one_profile(args, cfg, profile, args.output_root))
        compare_summary = {}
        if args.run_profile == "both":
            compare_dir = args.output_root / "compare"
            compare_dir.mkdir(parents=True, exist_ok=True)
            compare_html = compare_dir / "transformnet_network.interactive.html"
            compare_viz = visualize_existing_network(
                network_dir=args.output_root / "strict" / "network",
                strict_network_dir=args.output_root / "strict" / "network",
                relaxed_network_dir=args.output_root / "relaxed" / "network",
                scaffold_json=args.scaffold_json,
                output_html=compare_html,
                valid_molecules=args.output_root / "strict" / "compute" / "transformnet.valid_molecules.csv",
                json_index_base=args.json_index_base,
                edge_id_space=args.edge_id_space,
                max_nodes=args.max_nodes,
                max_edges=args.max_edges,
                layout_seed=args.layout_seed,
                central_groups=args.central_groups,
                label_column=args.label_column,
                title=args.title + " strict/relaxed comparison",
                make_svg=not args.no_molecule_svg,
                scaffold_image=args.scaffold_image,
                node_size_scale=args.node_size_scale,
                label_size=args.label_size,
            )
            compare_summary = {"interactive_html": str(compare_html), "visual_summary": compare_viz.get("summary", "")}
        args.output_root.mkdir(parents=True, exist_ok=True)
        summary_path = args.output_root / "run_all_summary.json"
        summary_path.write_text(json.dumps({"runs": summaries, "compare": compare_summary}, indent=2, ensure_ascii=True), encoding="utf-8")
        print("End-to-end TransformNet run complete.")
        print(f"summary: {summary_path}")
        for item in summaries:
            print(f"profile: {item['profile']}")
            print(f"  interactive_html: {item['interactive_html']}")
            print(f"  hotspot_summary: {item['hotspot_summary']}")
        if compare_summary:
            print(f"compare_interactive_html: {compare_summary['interactive_html']}")
        return


if __name__ == "__main__":
    main()
