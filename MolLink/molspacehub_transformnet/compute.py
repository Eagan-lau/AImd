from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tqdm import tqdm

from .chem import StandardizationOptions, make_key_index, mol_descriptors, mol_from_smiles, mol_to_key, standardize_library, standardize_mol
from .io import load_molecule_table
from .logging import set_rdkit_log_state
from .templates import compile_templates
from .utils import ensure_dir, semicolon_join, write_json


@dataclass
class ComputeResult:
    paths: dict[str, Path]
    qc: dict[str, Any]


def _require_rdkit():
    try:
        from rdkit import Chem
    except Exception as exc:
        raise RuntimeError("RDKit is required for transformation prediction.") from exc
    return Chem


def _safe_product_record(
    product_mol,
    match_method: str,
    standardization: StandardizationOptions,
) -> tuple[dict[str, Any] | None, str, str]:
    """Return (record, failure_type, failure_error)."""
    Chem = _require_rdkit()
    try:
        product_mol = standardize_mol(product_mol, standardization)
    except Exception as exc:
        return None, "product_standardization_failed", str(exc)
    try:
        product_smiles = Chem.MolToSmiles(product_mol, isomericSmiles=True, canonical=True)
    except Exception as exc:
        return None, "product_smiles_failed", str(exc)
    key = mol_to_key(product_mol, method=match_method)
    if not key:
        return None, f"product_{match_method}_failed", ""
    desc = mol_descriptors(product_mol)
    return {
        "product_smiles": product_smiles,
        "product_key": key,
        "product_formula": desc.get("formula", ""),
        "product_exact_mass": desc.get("exact_mass", None),
    }, "", ""


def _aggregate_events(edge_events: pd.DataFrame) -> pd.DataFrame:
    if edge_events.empty:
        return pd.DataFrame(columns=[
            "substrate_matrix_index", "product_matrix_index", "substrate_id", "product_id",
            "event_count", "template_count", "template_uids", "template_classes",
            "product_key_count", "product_keys", "product_smiles_examples",
            "has_ambiguous_product_key", "max_target_ambiguity"
        ])

    grouped = []
    for (sidx, pidx), g in edge_events.groupby(["substrate_matrix_index", "product_matrix_index"], sort=True):
        grouped.append({
            "substrate_matrix_index": int(sidx),
            "product_matrix_index": int(pidx),
            "substrate_id": g["substrate_id"].iloc[0],
            "product_id": g["product_id"].iloc[0],
            "substrate_name": g.get("substrate_name", pd.Series([""])).iloc[0],
            "product_name": g.get("product_name", pd.Series([""])).iloc[0],
            "event_count": int(len(g)),
            "template_count": int(g["template_uid"].nunique()),
            "template_uids": semicolon_join(g["template_uid"], max_items=50),
            "template_classes": semicolon_join(g.get("reaction_class", []), max_items=20),
            "template_names": semicolon_join(g.get("reaction_name", []), max_items=20),
            "product_key_count": int(g["product_key"].nunique()),
            "product_keys": semicolon_join(g["product_key"], max_items=50),
            "product_smiles_examples": semicolon_join(g["product_smiles"], max_items=10),
            "has_ambiguous_product_key": bool((g["target_ambiguity_count"] > 1).any()),
            "max_target_ambiguity": int(g["target_ambiguity_count"].max()),
            "match_method": g["match_method"].iloc[0],
        })
    return pd.DataFrame(grouped)


def _matrix_from_edges(valid: pd.DataFrame, directed_edges: pd.DataFrame, weight_col: str = "event_count"):
    n = int(len(valid))
    binary_directed = np.zeros((n, n), dtype=np.int16)
    weighted_directed = np.zeros((n, n), dtype=np.float64)
    if not directed_edges.empty:
        for _, edge in directed_edges.iterrows():
            s = int(edge["substrate_matrix_index"])
            p = int(edge["product_matrix_index"])
            if s < 0 or p < 0 or s >= n or p >= n:
                continue
            binary_directed[s, p] = 1
            weighted_directed[s, p] += float(edge.get(weight_col, 1))
    binary_undirected = np.maximum(binary_directed, binary_directed.T)
    weighted_undirected = weighted_directed + weighted_directed.T
    return binary_directed, weighted_directed, binary_undirected, weighted_undirected


def compute_transformations(
    molecule_table: str | Path,
    template_path: str | Path,
    output_dir: str | Path,
    smiles_column: str | None = None,
    id_column: str | None = None,
    name_column: str | None = None,
    sheet_name: str | int | None = None,
    template_column: str | None = None,
    match_method: str = "full_inchikey",
    duplicate_key_policy: str = "all",
    self_loop_policy: str = "drop",
    single_reactant_only: bool = True,
    prefilter_substructure: bool = True,
    max_molecules: int | None = None,
    max_templates: int | None = None,
    standardize_cleanup: bool = True,
    largest_fragment: bool = False,
    uncharge: bool = False,
    remove_hs: bool = False,
    quiet_rdkit: bool = True,
    write_raw_product_failures: bool = True,
    write_unmatched_products: bool = False,
    prefix: str = "transformnet",
) -> ComputeResult:
    """Compute a strict molecular transformation network from a molecule table and reaction SMARTS.

    Strict default:
    - `match_method='full_inchikey'` preserves stereochemical distinctions.
    - `self_loop_policy='drop'` removes substrate->itself edges.
    - duplicate target keys are retained but explicitly flagged unless duplicate_key_policy is changed.
    """
    start = time.time()
    if self_loop_policy not in {"drop", "keep"}:
        raise ValueError("self_loop_policy must be 'drop' or 'keep'")
    set_rdkit_log_state(quiet=quiet_rdkit)
    Chem = _require_rdkit()
    output_dir = ensure_dir(output_dir)

    std = StandardizationOptions(
        sanitize=True,
        cleanup=standardize_cleanup,
        largest_fragment=largest_fragment,
        uncharge=uncharge,
        remove_hs=remove_hs,
    )

    molecules_raw = load_molecule_table(
        molecule_table,
        smiles_column=smiles_column,
        id_column=id_column,
        name_column=name_column,
        sheet_name=sheet_name,
    )
    valid, invalid, duplicate_keys = standardize_library(molecules_raw, match_method=match_method, standardization=std)
    if max_molecules is not None:
        valid = valid.iloc[:max_molecules].copy()
        valid["matrix_index"] = range(len(valid))

    key_to_targets, key_index_records = make_key_index(valid, duplicate_policy=duplicate_key_policy)

    templates, template_qc = compile_templates(
        template_path=template_path,
        template_column=template_column,
        max_templates=max_templates,
        single_reactant_only=single_reactant_only,
        require_product=True,
    )

    edge_event_records: list[dict[str, Any]] = []
    product_failure_records: list[dict[str, Any]] = []
    unmatched_product_records: list[dict[str, Any]] = []

    reaction_run_errors = 0
    product_failure_counts: dict[str, int] = {}
    substructure_prefilter_skips = 0
    product_tuples_seen = 0
    products_seen = 0
    unique_predicted_product_keys: set[str] = set()
    seen_events: set[tuple[int, int, str, str]] = set()

    iterator = valid.iterrows()
    for _, row in tqdm(iterator, total=len(valid), desc="predicting transformations"):
        substrate_idx = int(row["matrix_index"])
        substrate_id = str(row["molecule_id"])
        substrate_name = str(row.get("molecule_name", substrate_id))
        substrate_mol = mol_from_smiles(str(row["standardized_smiles"]), sanitize=True)
        if substrate_mol is None:
            continue

        for tmpl in templates:
            rxn = tmpl["rxn"]
            template_uid = str(tmpl.get("template_uid", tmpl.get("template_index", "")))
            reactant_template = tmpl.get("reactant_template")
            if prefilter_substructure and reactant_template is not None:
                try:
                    if not substrate_mol.HasSubstructMatch(reactant_template):
                        substructure_prefilter_skips += 1
                        continue
                except Exception:
                    pass

            try:
                product_tuples = rxn.RunReactants((substrate_mol,))
            except Exception as exc:
                reaction_run_errors += 1
                product_failure_counts["reaction_run_failed"] = product_failure_counts.get("reaction_run_failed", 0) + 1
                if write_raw_product_failures:
                    product_failure_records.append({
                        "substrate_matrix_index": substrate_idx,
                        "substrate_id": substrate_id,
                        "template_uid": template_uid,
                        "failure_type": "reaction_run_failed",
                        "error": str(exc),
                    })
                continue

            if not product_tuples:
                continue

            for tuple_index, product_tuple in enumerate(product_tuples):
                product_tuples_seen += 1
                for product_position, product_mol in enumerate(product_tuple):
                    products_seen += 1
                    product_rec, failure_type, error = _safe_product_record(product_mol, match_method=match_method, standardization=std)
                    if product_rec is None:
                        product_failure_counts[failure_type] = product_failure_counts.get(failure_type, 0) + 1
                        if write_raw_product_failures:
                            product_failure_records.append({
                                "substrate_matrix_index": substrate_idx,
                                "substrate_id": substrate_id,
                                "template_uid": template_uid,
                                "failure_type": failure_type,
                                "error": error,
                                "tuple_index": tuple_index,
                                "product_position": product_position,
                            })
                        continue

                    product_key = str(product_rec["product_key"])
                    unique_predicted_product_keys.add(product_key)
                    target_indices = key_to_targets.get(product_key, [])
                    if not target_indices:
                        if write_unmatched_products:
                            unmatched_product_records.append({
                                "substrate_matrix_index": substrate_idx,
                                "substrate_id": substrate_id,
                                "template_uid": template_uid,
                                "product_key": product_key,
                                "product_smiles": product_rec["product_smiles"],
                                "product_formula": product_rec.get("product_formula", ""),
                                "product_exact_mass": product_rec.get("product_exact_mass", None),
                            })
                        continue

                    for target_idx in target_indices:
                        if target_idx == substrate_idx and self_loop_policy == "drop":
                            continue
                        target_row = valid.iloc[int(target_idx)]
                        event_key = (substrate_idx, int(target_idx), template_uid, product_key)
                        if event_key in seen_events:
                            continue
                        seen_events.add(event_key)
                        edge_event_records.append({
                            "substrate_matrix_index": substrate_idx,
                            "product_matrix_index": int(target_idx),
                            "substrate_id": substrate_id,
                            "product_id": str(target_row["molecule_id"]),
                            "substrate_name": substrate_name,
                            "product_name": str(target_row.get("molecule_name", target_row["molecule_id"])),
                            "template_uid": template_uid,
                            "template_index": int(tmpl.get("template_index", -1)),
                            "source_line": tmpl.get("source_line", ""),
                            "reaction_class": tmpl.get("reaction_class", tmpl.get("class", "")),
                            "reaction_name": tmpl.get("reaction_name", tmpl.get("name", "")),
                            "enzyme_family": tmpl.get("enzyme_family", ""),
                            "product_key": product_key,
                            "product_smiles": product_rec["product_smiles"],
                            "product_formula": product_rec.get("product_formula", ""),
                            "product_exact_mass": product_rec.get("product_exact_mass", None),
                            "target_ambiguity_count": int(len(target_indices)),
                            "is_ambiguous_product_key": bool(len(target_indices) > 1),
                            "match_method": match_method,
                            "product_tuple_size": int(len(product_tuple)),
                            "product_position": int(product_position),
                        })

    edge_events = pd.DataFrame(edge_event_records)
    directed_edges = _aggregate_events(edge_events)
    binary_directed, weighted_directed, binary_undirected, weighted_undirected = _matrix_from_edges(valid, directed_edges)

    # Original TaxLink-compatible files use product,substrate ordering for from_to.
    from_to = edge_events[["product_matrix_index", "substrate_matrix_index"]].drop_duplicates().copy() if not edge_events.empty else pd.DataFrame(columns=["product_matrix_index", "substrate_matrix_index"])

    paths: dict[str, Path] = {
        "valid_molecules": output_dir / f"{prefix}.valid_molecules.csv",
        "invalid_molecules": output_dir / f"{prefix}.invalid_molecules.csv",
        "duplicate_keys": output_dir / f"{prefix}.duplicate_keys.csv",
        "key_index": output_dir / f"{prefix}.key_index.csv",
        "template_qc": output_dir / f"{prefix}.template_qc.csv",
        "edge_events": output_dir / f"{prefix}.edge_events.csv",
        "directed_edges": output_dir / f"{prefix}.directed_edges.csv",
        "from_to_list": output_dir / "from_to_list.csv",
        "adjacency_directed_binary": output_dir / f"{prefix}.adjacency_directed_binary.csv",
        "adjacency_directed_weighted": output_dir / f"{prefix}.adjacency_directed_weighted.csv",
        "adjacency_undirected_binary": output_dir / f"{prefix}.adjacency_undirected_binary.csv",
        "adjacency_undirected_weighted": output_dir / f"{prefix}.adjacency_undirected_weighted.csv",
        "linkage_new": output_dir / "linkage_new.csv",
        "linkage_new_matrix_index": output_dir / "linkage_new.matrix_index.csv",
        "summary": output_dir / f"{prefix}.summary.json",
    }

    valid.to_csv(paths["valid_molecules"], index=False)
    invalid.to_csv(paths["invalid_molecules"], index=False)
    duplicate_keys.to_csv(paths["duplicate_keys"], index=False)
    key_index_records.to_csv(paths["key_index"], index=False)
    template_qc.to_csv(paths["template_qc"], index=False)
    edge_events.to_csv(paths["edge_events"], index=False)
    directed_edges.to_csv(paths["directed_edges"], index=False)
    from_to.to_csv(paths["from_to_list"], index=False)
    pd.DataFrame(binary_directed).to_csv(paths["adjacency_directed_binary"])
    pd.DataFrame(weighted_directed).to_csv(paths["adjacency_directed_weighted"])
    pd.DataFrame(binary_undirected).to_csv(paths["adjacency_undirected_binary"])
    pd.DataFrame(weighted_undirected).to_csv(paths["adjacency_undirected_weighted"])
    source_ids = [str(x) for x in valid.get("source_id", valid.get("molecule_id")).tolist()] if not valid.empty else []
    linkage_by_source = pd.DataFrame(binary_undirected, index=source_ids, columns=source_ids)
    linkage_by_source.index.name = "source_id"
    linkage_by_source.to_csv(paths["linkage_new"])
    pd.DataFrame(binary_undirected).to_csv(paths["linkage_new_matrix_index"])

    if write_raw_product_failures:
        product_failures = pd.DataFrame(product_failure_records)
        paths["product_failures"] = output_dir / f"{prefix}.product_failures.csv"
        product_failures.to_csv(paths["product_failures"], index=False)
    else:
        product_failures = pd.DataFrame()

    if write_unmatched_products:
        unmatched_products = pd.DataFrame(unmatched_product_records)
        paths["unmatched_products"] = output_dir / f"{prefix}.unmatched_products.csv"
        unmatched_products.to_csv(paths["unmatched_products"], index=False)
    else:
        unmatched_products = pd.DataFrame()

    qc = {
        "module_version": "0.6.1",
        "input_table_encoding": str(molecules_raw.attrs.get("source_encoding", "")),
        "mode": "template_based_transformnet_compute",
        "molecule_table": str(molecule_table),
        "template_path": str(template_path),
        "input_molecule_rows": int(len(molecules_raw)),
        "valid_molecules": int(len(valid)),
        "invalid_molecules": int(len(invalid)),
        "duplicate_match_keys": int(len(duplicate_keys)),
        "duplicate_key_policy": duplicate_key_policy,
        "match_method": match_method,
        "self_loop_policy": self_loop_policy,
        "single_reactant_only": bool(single_reactant_only),
        "prefilter_substructure": bool(prefilter_substructure),
        "standardization": {
            "cleanup": bool(standardize_cleanup),
            "largest_fragment": bool(largest_fragment),
            "uncharge": bool(uncharge),
            "remove_hs": bool(remove_hs),
        },
        "templates_seen": int(len(template_qc)),
        "compiled_templates": int(len(templates)),
        "template_status_counts": template_qc["status"].value_counts().to_dict() if not template_qc.empty else {},
        "substructure_prefilter_skips": int(substructure_prefilter_skips),
        "reaction_run_errors": int(reaction_run_errors),
        "product_tuples_seen": int(product_tuples_seen),
        "products_seen": int(products_seen),
        "unique_predicted_product_keys": int(len(unique_predicted_product_keys)),
        "product_failure_counts": {str(k): int(v) for k, v in product_failure_counts.items()},
        "product_failure_records": int(len(product_failures)),
        "unmatched_product_records": int(len(unmatched_products)),
        "edge_events": int(len(edge_events)),
        "directed_edges": int(len(directed_edges)),
        "from_to_records": int(len(from_to)),
        "matrix_shape": [int(len(valid)), int(len(valid))],
        "directed_binary_nonzero": int(binary_directed.sum()),
        "undirected_binary_nonzero": int(binary_undirected.sum()),
        "self_loops_directed": int(np.trace(binary_directed)),
        "runtime_seconds": round(time.time() - start, 3),
        "outputs": {k: str(v) for k, v in paths.items()},
    }
    write_json(qc, paths["summary"])
    return ComputeResult(paths=paths, qc=qc)
