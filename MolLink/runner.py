from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .molspacehub_transformnet import __version__ as MOLLINK_VERSION
from .molspacehub_transformnet.compute import compute_transformations
from .molspacehub_transformnet.io import load_molecule_table
from .molspacehub_transformnet.network import build_network_from_computed


@dataclass
class MolLinkRunResult:
    paths: dict[str, Path]
    summary: dict[str, Any]


def _find_project_root(start: Path) -> Path:
    for candidate in [start, *start.parents]:
        if (candidate / "AGENTS.md").is_file() and (candidate / "data").is_dir():
            return candidate.resolve()
    return Path.cwd().resolve()


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _resolve_path(value: str | Path | None, root: Path) -> Path | None:
    if value is None or str(value).strip() == "":
        return None
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def _relative(path: Path | None, root: Path) -> str:
    if path is None:
        return ""
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)


def _write_json(payload: dict[str, Any], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=True, indent=2)
    return path


def _profile_defaults(profile: str) -> dict[str, str]:
    if profile == "strict":
        return {"match_method": "full_inchikey", "self_loop_policy": "drop"}
    if profile == "relaxed":
        return {"match_method": "connectivity_inchikey", "self_loop_policy": "drop"}
    if profile == "original":
        return {"match_method": "connectivity_inchikey", "self_loop_policy": "keep"}
    raise ValueError(f"Unsupported MolLink profile: {profile}")


def _empty_directed_edges() -> pd.DataFrame:
    return pd.DataFrame(columns=[
        "substrate_matrix_index",
        "product_matrix_index",
        "substrate_id",
        "product_id",
        "event_count",
        "template_count",
        "template_uids",
        "template_classes",
        "product_key_count",
        "product_keys",
        "product_smiles_examples",
        "has_ambiguous_product_key",
        "max_target_ambiguity",
    ])


def _make_ligand_source_manifest(valid: pd.DataFrame, molecule_table: Path, root: Path) -> pd.DataFrame:
    manifest = pd.DataFrame()
    manifest["ligand_id"] = valid["source_id"].astype(str)
    manifest["source_id"] = valid["source_id"].astype(str)
    manifest["molecule_id"] = valid["molecule_id"].astype(str)
    manifest["molecule_name"] = valid["molecule_name"].astype(str)
    manifest["smiles"] = valid["smiles_raw"].astype(str)
    manifest["ligand_source_type"] = "smiles_csv"
    manifest["source_table"] = _relative(molecule_table, root)
    manifest["row_index"] = valid["row_index"].astype(int)
    manifest["transform_status"] = "loaded"
    return manifest


def _write_ligand_source_manifest(valid: pd.DataFrame, molecule_table: Path, root: Path, path: Path) -> Path:
    manifest = _make_ligand_source_manifest(valid, molecule_table, root)
    path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(path, index=False)
    return path


def _write_csv_only_outputs(
    molecule_table: Path,
    output_root: Path,
    root: Path,
    columns: dict[str, Any],
    prefix: str,
    max_molecules: int | None,
    template_status: str,
) -> tuple[dict[str, Path], dict[str, Any]]:
    compute_dir = output_root / "compute"
    network_dir = output_root / "network"
    compute_dir.mkdir(parents=True, exist_ok=True)
    network_dir.mkdir(parents=True, exist_ok=True)

    molecules = load_molecule_table(
        molecule_table,
        smiles_column=columns.get("smiles"),
        id_column=columns.get("id"),
        name_column=columns.get("name"),
        sheet_name=columns.get("sheet_name"),
    )
    if max_molecules is not None:
        molecules = molecules.iloc[:max_molecules].copy()

    empty_smiles = molecules["smiles_raw"].astype(str).str.strip().eq("")
    valid = molecules.loc[~empty_smiles].copy()
    invalid = molecules.loc[empty_smiles].copy()
    valid["matrix_index"] = range(len(valid))
    valid["standardized_smiles"] = valid["smiles_raw"].astype(str)
    valid["transform_status"] = "csv_loaded"
    if invalid.empty:
        invalid = pd.DataFrame(columns=list(molecules.columns) + ["drop_reason"])
    else:
        invalid["drop_reason"] = "empty_smiles"

    duplicate_keys = pd.DataFrame(columns=["molecule_id", "duplicate_group", "drop_reason"])
    key_index = valid[["matrix_index", "source_id", "molecule_id"]].copy()
    key_index["match_key"] = ""
    template_qc = pd.DataFrame(columns=["template_uid", "template_index", "reaction_smarts", "status", "error"])
    edge_events = pd.DataFrame(columns=[
        "substrate_matrix_index",
        "product_matrix_index",
        "substrate_id",
        "product_id",
        "template_uid",
        "reaction_class",
        "reaction_name",
        "product_key",
        "product_smiles",
        "match_method",
    ])
    directed_edges = _empty_directed_edges()
    from_to = pd.DataFrame(columns=["product_matrix_index", "substrate_matrix_index"])

    n = int(len(valid))
    zeros = np.zeros((n, n), dtype=np.int16)
    zero_float = np.zeros((n, n), dtype=np.float64)
    source_ids = valid["source_id"].astype(str).tolist()

    paths = {
        "compute_dir": compute_dir,
        "network_dir": network_dir,
        "valid_molecules": compute_dir / f"{prefix}.valid_molecules.csv",
        "invalid_molecules": compute_dir / f"{prefix}.invalid_molecules.csv",
        "duplicate_keys": compute_dir / f"{prefix}.duplicate_keys.csv",
        "key_index": compute_dir / f"{prefix}.key_index.csv",
        "template_qc": compute_dir / f"{prefix}.template_qc.csv",
        "edge_events": compute_dir / f"{prefix}.edge_events.csv",
        "directed_edges": compute_dir / f"{prefix}.directed_edges.csv",
        "from_to_list": compute_dir / "from_to_list.csv",
        "adjacency_directed_binary": compute_dir / f"{prefix}.adjacency_directed_binary.csv",
        "adjacency_directed_weighted": compute_dir / f"{prefix}.adjacency_directed_weighted.csv",
        "adjacency_undirected_binary": compute_dir / f"{prefix}.adjacency_undirected_binary.csv",
        "adjacency_undirected_weighted": compute_dir / f"{prefix}.adjacency_undirected_weighted.csv",
        "linkage_new": compute_dir / "linkage_new.csv",
        "linkage_new_matrix_index": compute_dir / "linkage_new.matrix_index.csv",
        "summary": compute_dir / f"{prefix}.summary.json",
        "ligand_source_manifest": output_root / "ligand_source_manifest.csv",
        "run_summary": output_root / "mollink_run_summary.json",
    }

    valid.to_csv(paths["valid_molecules"], index=False)
    invalid.to_csv(paths["invalid_molecules"], index=False)
    duplicate_keys.to_csv(paths["duplicate_keys"], index=False)
    key_index.to_csv(paths["key_index"], index=False)
    template_qc.to_csv(paths["template_qc"], index=False)
    edge_events.to_csv(paths["edge_events"], index=False)
    directed_edges.to_csv(paths["directed_edges"], index=False)
    from_to.to_csv(paths["from_to_list"], index=False)
    pd.DataFrame(zeros).to_csv(paths["adjacency_directed_binary"])
    pd.DataFrame(zero_float).to_csv(paths["adjacency_directed_weighted"])
    pd.DataFrame(zeros).to_csv(paths["adjacency_undirected_binary"])
    pd.DataFrame(zero_float).to_csv(paths["adjacency_undirected_weighted"])
    pd.DataFrame(zeros, index=source_ids, columns=source_ids).to_csv(paths["linkage_new"])
    pd.DataFrame(zeros).to_csv(paths["linkage_new_matrix_index"])
    _write_ligand_source_manifest(valid, molecule_table, root, paths["ligand_source_manifest"])

    compute_summary = {
        "module": "MolLink",
        "module_version": MOLLINK_VERSION,
        "mode": "csv_only_no_template",
        "molecule_table": _relative(molecule_table, root),
        "template_status": template_status,
        "input_molecule_rows": int(len(molecules)),
        "valid_molecules": int(len(valid)),
        "invalid_molecules": int(len(invalid)),
        "directed_edges": 0,
        "edge_events": 0,
        "message": "No reaction template file was used; transformation edges were not inferred.",
    }
    _write_json(compute_summary, paths["summary"])

    return paths, compute_summary


def _build_network(paths: dict[str, Path], cfg: dict[str, Any], root: Path) -> dict[str, Any]:
    network_cfg = cfg.get("network", {}) or {}
    try:
        result = build_network_from_computed(
            computed_dir=paths["compute_dir"],
            output_dir=paths["network_dir"],
            annotation_table=_resolve_path(cfg.get("inputs", {}).get("annotation_table"), root),
            annotation_sheet=cfg.get("columns", {}).get("annotation_sheet"),
            annotation_id_column=cfg.get("columns", {}).get("annotation_id"),
            group_column=network_cfg.get("group_column"),
            position_columns=network_cfg.get("position_columns"),
            prefix=network_cfg.get("prefix", "transformnet_network"),
            write_excel=bool(network_cfg.get("write_excel", False)),
            write_graphs=bool(network_cfg.get("write_graphs", False)),
            write_html=bool(network_cfg.get("write_html", True)),
        )
    except ModuleNotFoundError as exc:
        if "openpyxl" in str(exc):
            raise RuntimeError("MolLink network Excel output requires openpyxl. Install openpyxl or set network.write_excel: false.") from exc
        raise
    for key, value in result.paths.items():
        paths[f"network_{key}"] = value
    return result.qc


def _run_template_mode(
    molecule_table: Path,
    template_file: Path,
    output_root: Path,
    root: Path,
    cfg: dict[str, Any],
    max_molecules: int | None,
    max_templates: int | None,
) -> tuple[dict[str, Path], dict[str, Any]]:
    compute_cfg = cfg.get("compute", {}) or {}
    columns = cfg.get("columns", {}) or {}
    standardization = cfg.get("standardization", {}) or {}
    prefix = compute_cfg.get("prefix", "transformnet")
    profile_name = compute_cfg.get("profile", "strict")
    defaults = _profile_defaults(profile_name)
    compute_dir = output_root / "compute"
    network_dir = output_root / "network"
    compute_dir.mkdir(parents=True, exist_ok=True)
    network_dir.mkdir(parents=True, exist_ok=True)

    result = compute_transformations(
        molecule_table=molecule_table,
        template_path=template_file,
        output_dir=compute_dir,
        smiles_column=columns.get("smiles"),
        id_column=columns.get("id"),
        name_column=columns.get("name"),
        sheet_name=columns.get("sheet_name"),
        template_column=columns.get("template"),
        match_method=compute_cfg.get("match_method") or defaults["match_method"],
        duplicate_key_policy=compute_cfg.get("duplicate_key_policy", "all"),
        self_loop_policy=compute_cfg.get("self_loop_policy") or defaults["self_loop_policy"],
        single_reactant_only=bool(compute_cfg.get("single_reactant_only", True)),
        prefilter_substructure=bool(compute_cfg.get("prefilter_substructure", True)),
        max_molecules=max_molecules if max_molecules is not None else compute_cfg.get("max_molecules"),
        max_templates=max_templates if max_templates is not None else compute_cfg.get("max_templates"),
        standardize_cleanup=bool(standardization.get("cleanup", True)),
        largest_fragment=bool(standardization.get("largest_fragment", False)),
        uncharge=bool(standardization.get("uncharge", False)),
        remove_hs=bool(standardization.get("remove_hs", False)),
        quiet_rdkit=bool(compute_cfg.get("quiet_rdkit", True)),
        write_raw_product_failures=bool(compute_cfg.get("write_product_failures", True)),
        write_unmatched_products=bool(compute_cfg.get("write_unmatched_products", False)),
        prefix=prefix,
    )
    paths = dict(result.paths)
    paths["compute_dir"] = compute_dir
    paths["network_dir"] = network_dir
    paths["ligand_source_manifest"] = output_root / "ligand_source_manifest.csv"
    paths["run_summary"] = output_root / "mollink_run_summary.json"

    valid = pd.read_csv(paths["valid_molecules"])
    _write_ligand_source_manifest(valid, molecule_table, root, paths["ligand_source_manifest"])
    summary = dict(result.qc)
    summary.update({
        "module": "MolLink",
        "module_version": MOLLINK_VERSION,
        "mode": "template_based_transformnet_compute",
        "molecule_table": _relative(molecule_table, root),
        "template_file": _relative(template_file, root),
    })
    return paths, summary


def run_mollink(
    config: dict[str, Any],
    config_path: str | Path | None = None,
    mode: str | None = None,
    max_molecules: int | None = None,
    max_templates: int | None = None,
) -> MolLinkRunResult:
    start = Path(config_path).resolve().parent if config_path else Path.cwd()
    detected_root = _find_project_root(start)
    project_cfg = config.get("project", {}) or {}
    root_value = project_cfg.get("root")
    root = detected_root if not root_value else _resolve_path(root_value, detected_root)
    if root is None:
        root = detected_root

    inputs = config.get("inputs", {}) or {}
    outputs = config.get("outputs", {}) or {}
    compute_cfg = config.get("compute", {}) or {}
    columns = config.get("columns", {}) or {}

    molecule_table = _resolve_path(inputs.get("molecule_table"), root)
    if molecule_table is None or not molecule_table.is_file():
        raise FileNotFoundError(f"MolLink molecule_table is missing: {inputs.get('molecule_table')}")

    output_root = _resolve_path(outputs.get("output_root", "data/data_output/ligand_transformation"), root)
    if output_root is None:
        raise ValueError("MolLink output_root could not be resolved.")
    output_root.mkdir(parents=True, exist_ok=True)

    configured_mode = mode or compute_cfg.get("mode", "auto")
    if configured_mode not in {"auto", "csv_only", "template"}:
        raise ValueError("MolLink compute.mode must be one of: auto, csv_only, template")
    template_file = _resolve_path(inputs.get("template_file"), root)
    require_templates = bool(compute_cfg.get("require_templates", False))
    prefix = compute_cfg.get("prefix", "transformnet")

    template_exists = bool(template_file and template_file.is_file())
    if configured_mode == "template" and not template_exists:
        raise FileNotFoundError(f"MolLink template_file is required for template mode but was not found: {inputs.get('template_file')}")
    if require_templates and not template_exists:
        raise FileNotFoundError(f"MolLink require_templates is true but template_file was not found: {inputs.get('template_file')}")

    if configured_mode == "csv_only" or not template_exists:
        template_status = "disabled" if not template_file else "missing"
        paths, summary = _write_csv_only_outputs(
            molecule_table=molecule_table,
            output_root=output_root,
            root=root,
            columns=columns,
            prefix=prefix,
            max_molecules=max_molecules if max_molecules is not None else compute_cfg.get("max_molecules"),
            template_status=template_status,
        )
    else:
        paths, summary = _run_template_mode(
            molecule_table=molecule_table,
            template_file=template_file,
            output_root=output_root,
            root=root,
            cfg=config,
            max_molecules=max_molecules,
            max_templates=max_templates,
        )

    network_qc = _build_network(paths, config, root)
    summary.update({
        "output_root": _relative(output_root, root),
        "ligand_source_manifest": _relative(paths["ligand_source_manifest"], root),
        "network_summary": network_qc,
    })
    _write_json(summary, paths["run_summary"])
    return MolLinkRunResult(paths=paths, summary=summary)


def run_mollink_from_config(
    config_path: str | Path,
    mode: str | None = None,
    max_molecules: int | None = None,
    max_templates: int | None = None,
) -> MolLinkRunResult:
    path = Path(config_path).resolve()
    config = _load_yaml(path)
    return run_mollink(
        config,
        config_path=path,
        mode=mode,
        max_molecules=max_molecules,
        max_templates=max_templates,
    )
