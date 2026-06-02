from __future__ import annotations

from pathlib import Path
from typing import Iterable, Any

import numpy as np
import pandas as pd

SMILES_CANDIDATES = [
    "Isomeric SMILES",
    "isomeric_smiles",
    "SMILES",
    "Smiles",
    "smiles",
    "canonical_smiles",
    "Canonical SMILES",
]
ID_CANDIDATES = [
    "molecule_id",
    "Molecule ID",
    "source_no",
    "source_id",
    "library_index",
    "compound_id",
    "Compound ID",
    "Unnamed: 0",
    "No",
    "no",
    "id",
    "ID",
]
NAME_CANDIDATES = ["Compound", "compound", "name", "Name", "molecule_name", "Molecule Name"]


def detect_column(columns: Iterable[str], candidates: list[str], required: bool = True, label: str = "column") -> str | None:
    cols = list(columns)
    normalized = {str(c).strip().lower(): c for c in cols}
    for cand in candidates:
        key = cand.strip().lower()
        if key in normalized:
            return normalized[key]
    if required:
        raise ValueError(f"Could not detect {label}. Available columns: {cols}")
    return None


CSV_ENCODINGS = ["utf-8", "utf-8-sig", "gbk", "latin1"]


def read_csv_with_fallback(path: str | Path, sep: str = ",") -> pd.DataFrame:
    """Read a delimited text file with common encoding fallbacks.

    The selected encoding is stored in DataFrame.attrs["source_encoding"].
    """
    path = Path(path)
    errors: list[str] = []
    for encoding in CSV_ENCODINGS:
        try:
            df = pd.read_csv(path, sep=sep, encoding=encoding)
            df.attrs["source_encoding"] = encoding
            return df
        except UnicodeDecodeError as exc:
            errors.append(f"{encoding}: {exc}")
            continue
    raise UnicodeDecodeError("unknown", b"", 0, 1, f"Failed to read {path} with encodings: {'; '.join(errors)}")


def open_text_with_fallback(path: str | Path):
    path = Path(path)
    for encoding in CSV_ENCODINGS:
        try:
            handle = path.open("r", encoding=encoding)
            handle.read(1024)
            handle.seek(0)
            return handle
        except UnicodeDecodeError:
            try:
                handle.close()
            except Exception:
                pass
            continue
    return path.open("r", encoding="utf-8", errors="replace")


def read_table(path: str | Path, sheet_name: str | int | None = None) -> pd.DataFrame:
    """Read CSV/TSV/TXT/XLSX/SDF into a DataFrame.

    TXT is treated as CSV unless the suffix is .smi or .smiles, in which case each
    non-empty line is parsed as `SMILES [optional_id]`.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()

    if suffix in {".xlsx", ".xls"}:
        if isinstance(sheet_name, str) and sheet_name.isdigit():
            sheet_name = int(sheet_name)
        return pd.read_excel(path, sheet_name=0 if sheet_name is None else sheet_name)

    if suffix in {".csv", ".tsv"}:
        sep = "\t" if suffix == ".tsv" else ","
        return read_csv_with_fallback(path, sep=sep)

    if suffix in {".smi", ".smiles"}:
        rows = []
        with open_text_with_fallback(path) as handle:
            for i, line in enumerate(handle, start=1):
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                parts = text.split()
                rows.append({"molecule_id": parts[1] if len(parts) > 1 else i, "smiles": parts[0]})
        return pd.DataFrame(rows)

    if suffix == ".txt":
        # Try CSV first; fall back to line-by-line SMILES.
        try:
            df = read_csv_with_fallback(path)
            if len(df.columns) > 1:
                return df
        except Exception:
            pass
        rows = []
        with open_text_with_fallback(path) as handle:
            for i, line in enumerate(handle, start=1):
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                parts = text.split()
                rows.append({"molecule_id": parts[1] if len(parts) > 1 else i, "smiles": parts[0]})
        return pd.DataFrame(rows)

    if suffix == ".sdf":
        return read_sdf_table(path)

    raise ValueError(f"Unsupported table file type: {path}")


def read_sdf_table(path: str | Path) -> pd.DataFrame:
    try:
        from rdkit import Chem
    except Exception as exc:
        raise RuntimeError("RDKit is required to read SDF molecule libraries.") from exc
    path = Path(path)
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    rows: list[dict[str, Any]] = []
    for idx, mol in enumerate(supplier, start=1):
        if mol is None:
            rows.append({"molecule_id": idx, "smiles": "", "sdf_parse_error": True})
            continue
        props = {p: mol.GetProp(p) for p in mol.GetPropNames()}
        smi = Chem.MolToSmiles(mol, isomericSmiles=True)
        props.setdefault("molecule_id", props.get("_Name", idx))
        props["smiles"] = smi
        rows.append(props)
    return pd.DataFrame(rows)


def _make_unique_ids(raw_ids: pd.Series, row_indices: pd.Series) -> pd.Series:
    ids = raw_ids.astype(str).replace({"nan": "", "None": ""})
    ids = ids.where(ids.str.len() > 0, row_indices.add(1).astype(str))
    duplicated = ids.duplicated(keep=False)
    if duplicated.any():
        ids = ids.where(~duplicated, ids + "__row" + row_indices.add(1).astype(str))
    return ids


def load_molecule_table(
    path: str | Path,
    smiles_column: str | None = None,
    id_column: str | None = None,
    name_column: str | None = None,
    sheet_name: str | int | None = None,
    keep_all_input_columns: bool = True,
) -> pd.DataFrame:
    raw = read_table(path, sheet_name=sheet_name)
    source_encoding = raw.attrs.get("source_encoding", "")
    if raw.empty:
        raise ValueError(f"Molecule table is empty: {path}")

    smiles_column = smiles_column or detect_column(raw.columns, SMILES_CANDIDATES, label="SMILES column")
    id_column = id_column or detect_column(raw.columns, ID_CANDIDATES, required=False, label="molecule ID column")
    if not name_column:
        columns = list(raw.columns)
        if len(columns) >= 3:
            name_column = columns[2]
        elif columns:
            name_column = columns[0]

    out = pd.DataFrame()
    out["row_index"] = np.arange(len(raw), dtype=int)
    if id_column and id_column in raw.columns:
        out["source_id"] = raw[id_column].astype(str)
    else:
        out["source_id"] = out["row_index"].add(1).astype(str)
    out["molecule_id"] = _make_unique_ids(out["source_id"], out["row_index"])
    if name_column and name_column in raw.columns:
        out["molecule_name"] = raw[name_column].astype(str)
    else:
        out["molecule_name"] = out["molecule_id"]
    out["smiles_raw"] = raw[smiles_column].astype(str).replace({"nan": ""})
    out.attrs["source_encoding"] = source_encoding

    if keep_all_input_columns:
        for col in raw.columns:
            safe_col = str(col)
            if safe_col in out.columns:
                safe_col = f"input_{safe_col}"
            out[safe_col] = raw[col]
    return out


def inspect_molecule_library(
    molecule_table: str | Path,
    smiles_column: str | None = None,
    id_column: str | None = None,
    name_column: str | None = None,
    sheet_name: str | int | None = None,
) -> dict[str, Any]:
    df = load_molecule_table(
        molecule_table,
        smiles_column=smiles_column,
        id_column=id_column,
        name_column=name_column,
        sheet_name=sheet_name,
    )
    return {
        "path": str(molecule_table),
        "rows": int(len(df)),
        "detected_columns": {
            "molecule_id": "molecule_id",
            "molecule_name": "molecule_name",
            "smiles_raw": "smiles_raw",
        },
        "empty_smiles_rows": int((df["smiles_raw"].str.strip() == "").sum()),
        "duplicate_source_ids": int(df["source_id"].duplicated(keep=False).sum()),
        "duplicate_molecule_ids_after_disambiguation": int(df["molecule_id"].duplicated(keep=False).sum()),
        "first_rows": df[["molecule_id", "molecule_name", "smiles_raw"]].head(5).to_dict(orient="records"),
    }


def load_annotation_table(path: str | Path | None, sheet_name: str | int | None = None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    df = read_table(path, sheet_name=sheet_name)
    return df


def write_dataframe_csv(df: pd.DataFrame, path: str | Path, index: bool = False) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=index)
    return p
