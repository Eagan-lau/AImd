from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


def require_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors, inchi
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except Exception as exc:
        raise RuntimeError(
            "RDKit is required. Recommended installation: conda install -c conda-forge rdkit"
        ) from exc
    return Chem, Descriptors, rdMolDescriptors, inchi, rdMolStandardize


MATCH_METHODS = {
    "full_inchikey",
    "connectivity_inchikey",
    "canonical_smiles",
    "canonical_smiles_no_stereo",
}


@dataclass
class StandardizationOptions:
    sanitize: bool = True
    cleanup: bool = True
    largest_fragment: bool = False
    uncharge: bool = False
    remove_hs: bool = False


def standardize_mol(mol, options: StandardizationOptions | None = None):
    Chem, _Descriptors, _rdMolDescriptors, _inchi, rdMolStandardize = require_rdkit()
    options = options or StandardizationOptions()
    if mol is None:
        return None
    out = mol
    try:
        if options.sanitize:
            Chem.SanitizeMol(out)
        if options.cleanup:
            out = rdMolStandardize.Cleanup(out)
        if options.largest_fragment:
            out = rdMolStandardize.LargestFragmentChooser().choose(out)
        if options.uncharge:
            out = rdMolStandardize.Uncharger().uncharge(out)
        if options.remove_hs:
            out = Chem.RemoveHs(out)
        if options.sanitize:
            Chem.SanitizeMol(out)
    except Exception:
        raise
    return out


def mol_from_smiles(smiles: str, sanitize: bool = True):
    Chem, *_ = require_rdkit()
    if smiles is None:
        return None
    text = str(smiles).strip()
    if not text or text.lower() == "nan":
        return None
    try:
        return Chem.MolFromSmiles(text, sanitize=sanitize)
    except Exception:
        return None


def mol_to_key(mol, method: str = "full_inchikey") -> str | None:
    Chem, _Descriptors, _rdMolDescriptors, inchi, _std = require_rdkit()
    if method not in MATCH_METHODS:
        raise ValueError(f"Unsupported match method: {method}. Choose from {sorted(MATCH_METHODS)}")
    try:
        if method == "full_inchikey":
            return inchi.MolToInchiKey(mol) or None
        if method == "connectivity_inchikey":
            key = inchi.MolToInchiKey(mol)
            return key.split("-")[0] if key else None
        if method == "canonical_smiles":
            return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        if method == "canonical_smiles_no_stereo":
            return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
    except Exception:
        return None
    return None


def mol_descriptors(mol) -> dict[str, Any]:
    _Chem, Descriptors, rdMolDescriptors, _inchi, _std = require_rdkit()
    try:
        formula = rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        formula = ""
    try:
        mass = float(Descriptors.ExactMolWt(mol))
    except Exception:
        mass = None
    return {"formula": formula, "exact_mass": mass}


def standardize_library(
    molecules: pd.DataFrame,
    match_method: str = "full_inchikey",
    standardization: StandardizationOptions | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Parse, standardize, key, and describe a molecule table.

    Returns (valid, invalid, duplicate_keys).
    """
    Chem, *_ = require_rdkit()
    standardization = standardization or StandardizationOptions()

    valid_records: list[dict[str, Any]] = []
    invalid_records: list[dict[str, Any]] = []
    for _, row in molecules.iterrows():
        record = row.to_dict()
        smiles = str(row.get("smiles_raw", "")).strip()
        mol = mol_from_smiles(smiles, sanitize=False)
        if mol is None:
            record.update({"status": "invalid", "drop_reason": "smiles_parse_failed"})
            invalid_records.append(record)
            continue
        try:
            mol = standardize_mol(mol, standardization)
        except Exception as exc:
            record.update({"status": "invalid", "drop_reason": "standardization_failed", "error": str(exc)})
            invalid_records.append(record)
            continue
        try:
            canonical = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        except Exception as exc:
            record.update({"status": "invalid", "drop_reason": "canonical_smiles_failed", "error": str(exc)})
            invalid_records.append(record)
            continue
        key = mol_to_key(mol, method=match_method)
        if not key:
            record.update({"status": "invalid", "drop_reason": f"{match_method}_failed"})
            invalid_records.append(record)
            continue
        desc = mol_descriptors(mol)
        record.update({
            "matrix_index": len(valid_records),
            "status": "valid",
            "standardized_smiles": canonical,
            "match_key": key,
            "match_method": match_method,
            **desc,
        })
        valid_records.append(record)

    valid = pd.DataFrame(valid_records)
    invalid = pd.DataFrame(invalid_records)

    if valid.empty:
        dup = pd.DataFrame(columns=["match_key", "duplicate_count", "molecule_ids"])
    else:
        dup = (
            valid.groupby("match_key", dropna=False)
            .agg(
                duplicate_count=("molecule_id", "size"),
                molecule_ids=("molecule_id", lambda x: ";".join(map(str, x))),
                source_ids=("source_id", lambda x: ";".join(map(str, x))),
            )
            .reset_index()
        )
        dup = dup[dup["duplicate_count"] > 1].copy()

    return valid, invalid, dup


def make_key_index(valid: pd.DataFrame, duplicate_policy: str = "all") -> tuple[dict[str, list[int]], pd.DataFrame]:
    """Build match key -> matrix index map.

    duplicate_policy:
    - all: keep all duplicate-key targets and flag ambiguous matches
    - first: keep only first occurrence for each key
    - drop: remove all keys that are duplicated
    - error: raise if duplicate keys exist
    """
    if duplicate_policy not in {"all", "first", "drop", "error"}:
        raise ValueError("duplicate_policy must be one of: all, first, drop, error")
    if valid.empty:
        return {}, pd.DataFrame()
    counts = valid["match_key"].value_counts()
    duplicated_keys = set(counts[counts > 1].index)
    if duplicate_policy == "error" and duplicated_keys:
        raise ValueError(f"Duplicate match keys found under strict duplicate_policy=error: {len(duplicated_keys)} duplicated keys")
    mapping: dict[str, list[int]] = {}
    used_records = []
    for _, row in valid.iterrows():
        key = row["match_key"]
        idx = int(row["matrix_index"])
        if duplicate_policy == "drop" and key in duplicated_keys:
            continue
        if duplicate_policy == "first" and key in mapping:
            continue
        mapping.setdefault(key, []).append(idx)
        used_records.append({"match_key": key, "matrix_index": idx, "molecule_id": row["molecule_id"], "is_duplicated_key": key in duplicated_keys})
    return mapping, pd.DataFrame(used_records)
