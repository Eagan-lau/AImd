from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from tqdm import tqdm

from .io import read_table


TEMPLATE_COLUMN_CANDIDATES = [
    "template",
    "Template",
    "reaction_smarts",
    "Reaction SMARTS",
    "smarts",
    "SMARTS",
    "rxn_smarts",
]


def iter_template_lines(path: str | Path, max_templates: int | None = None):
    path = Path(path)
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            text = line.strip()
            if not text or text.startswith("#") or text.lower() in {"template", "smarts", "reaction_smarts"}:
                continue
            yield line_no, text
            count += 1
            if max_templates is not None and count >= max_templates:
                break


def detect_template_column(columns, template_column: str | None = None) -> str:
    if template_column:
        if template_column not in columns:
            raise ValueError(f"Template column '{template_column}' not found. Available columns: {list(columns)}")
        return template_column
    normalized = {str(c).strip().lower(): c for c in columns}
    for cand in TEMPLATE_COLUMN_CANDIDATES:
        key = cand.strip().lower()
        if key in normalized:
            return normalized[key]
    raise ValueError(f"Could not detect reaction SMARTS/template column. Available columns: {list(columns)}")


def load_template_table(
    template_path: str | Path,
    template_column: str | None = None,
    max_templates: int | None = None,
) -> pd.DataFrame:
    """Load reaction SMARTS templates from txt/csv/tsv/xlsx.

    TXT files are treated as one reaction SMARTS per line. Tabular files may include
    metadata columns such as template_id, reaction_class, reaction_name, enzyme_family.
    """
    path = Path(template_path)
    suffix = path.suffix.lower()

    if suffix in {".txt", ".smarts", ".rxnsmarts"}:
        records = []
        for idx, (line_no, smarts) in enumerate(iter_template_lines(path, max_templates=max_templates)):
            # Allow tab-delimited "id<TAB>smarts" if the second token looks like a reaction.
            parts = smarts.split("\t")
            if len(parts) >= 2 and ">>" in parts[-1]:
                template_id = parts[0]
                smarts_value = parts[-1]
            else:
                template_id = f"T{idx:06d}"
                smarts_value = smarts
            records.append({
                "template_uid": template_id,
                "template_index": idx,
                "source_line": line_no,
                "reaction_smarts": smarts_value,
            })
        return pd.DataFrame(records)

    df = read_table(path)
    if max_templates is not None:
        df = df.head(max_templates).copy()
    col = detect_template_column(df.columns, template_column=template_column)
    df = df.copy()
    if "template_uid" not in df.columns:
        if "template_id" in df.columns:
            df["template_uid"] = df["template_id"].astype(str)
        elif "id" in df.columns:
            df["template_uid"] = df["id"].astype(str)
        else:
            df["template_uid"] = [f"T{i:06d}" for i in range(len(df))]
    df["template_index"] = range(len(df))
    df["reaction_smarts"] = df[col].astype(str)
    if "source_line" not in df.columns:
        df["source_line"] = df["template_index"] + 1
    return df


def compile_templates(
    template_path: str | Path,
    template_column: str | None = None,
    max_templates: int | None = None,
    single_reactant_only: bool = True,
    require_product: bool = True,
) -> tuple[list[dict[str, Any]], pd.DataFrame]:
    """Compile and QC reaction SMARTS templates."""
    from rdkit.Chem import AllChem

    templates = load_template_table(template_path, template_column=template_column, max_templates=max_templates)
    compiled: list[dict[str, Any]] = []
    qc_records: list[dict[str, Any]] = []

    for _, row in tqdm(templates.iterrows(), total=len(templates), desc="compiling templates"):
        smarts = str(row["reaction_smarts"]).strip()
        status = "ok"
        error = ""
        n_reactants = None
        n_products = None
        rxn = None
        reactant_template = None
        try:
            if ">>" not in smarts:
                status = "missing_reaction_arrow"
            else:
                rxn = AllChem.ReactionFromSmarts(smarts)
                if rxn is None:
                    status = "invalid_reaction"
                else:
                    try:
                        rxn.Initialize()
                    except Exception:
                        # RDKit may still allow RunReactants after a noisy Initialize warning.
                        pass
                    n_reactants = int(rxn.GetNumReactantTemplates())
                    n_products = int(rxn.GetNumProductTemplates())
                    if single_reactant_only and n_reactants != 1:
                        status = "skipped_multi_reactant"
                    elif require_product and n_products < 1:
                        status = "skipped_no_product"
                    else:
                        try:
                            reactant_template = rxn.GetReactantTemplate(0) if n_reactants >= 1 else None
                        except Exception:
                            reactant_template = None
                        rec = row.to_dict()
                        rec.update({
                            "rxn": rxn,
                            "reactant_template": reactant_template,
                            "n_reactants": n_reactants,
                            "n_products": n_products,
                        })
                        compiled.append(rec)
        except Exception as exc:
            status = "parse_error"
            error = str(exc)

        q = {k: v for k, v in row.to_dict().items() if k not in {"rxn", "reactant_template"}}
        q.update({
            "status": status,
            "n_reactants": n_reactants,
            "n_products": n_products,
            "error": error,
        })
        qc_records.append(q)

    return compiled, pd.DataFrame(qc_records)


def validate_templates(
    template_path: str | Path,
    template_column: str | None = None,
    max_templates: int | None = 1000,
    single_reactant_only: bool = False,
) -> pd.DataFrame:
    compiled, qc = compile_templates(
        template_path=template_path,
        template_column=template_column,
        max_templates=max_templates,
        single_reactant_only=single_reactant_only,
        require_product=True,
    )
    qc["would_be_compiled"] = qc["status"].eq("ok")
    qc.attrs["compiled_count"] = len(compiled)
    return qc
