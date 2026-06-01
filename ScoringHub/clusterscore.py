#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import math
import re
import zipfile
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

from .utils import ensure_dir, numeric, parse_vina_log, read_csv, resolve_path, safe_eval_formula, write_csv, write_json


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return str(value).strip() in {"", "NA", "NaN", "nan", "None"}


def _to_float(value: Any) -> float | None:
    return numeric(value, default=None)


def _read_protein_manifest(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None or not path.exists():
        return {}
    rows = read_csv(path)
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        pid = row.get("protein_id", "")
        if pid:
            out[pid] = row
    return out


def _rows_from_manifest(config: dict[str, Any]) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    docking_manifest = resolve_path(config.get("paths", {}).get("docking_result_manifest", "data/docking_out/docking_result_manifest.csv"), root)
    if docking_manifest is None or not docking_manifest.exists():
        return []
    rows = read_csv(docking_manifest)
    if not rows:
        return []
    num_aff = int(config.get("affinity_extraction", {}).get("num_affinities", 9) or 9)
    out: list[dict[str, Any]] = []
    for row in rows:
        best = _to_float(row.get("best_affinity"))
        parsed: dict[str, Any] = {}
        if best is None and row.get("log_path"):
            parsed = parse_vina_log(row["log_path"], num_affinities=num_aff)
            affs = parsed.get("affinities", [])
            best = affs[0] if affs else None
            if affs:
                row["affinities"] = ",".join(map(str, affs))
        if best is None:
            continue
        clean = dict(row)
        clean["best_affinity"] = best
        for key in ["grid_size", "grid_space", "exhaustiveness", "random_seed"]:
            if _is_missing(clean.get(key)) and parsed:
                clean[key] = parsed.get(key, "")
        out.append(clean)
    return out


def _parse_job_id_from_log_name(log_path: Path) -> dict[str, str]:
    stem = log_path.stem
    parts = stem.split("@")
    row = {"job_id": stem, "log_path": str(log_path), "batch_id": log_path.parent.name}
    if len(parts) >= 4:
        row.update({"ligand_id": parts[0], "protein_id": parts[1], "pocket_id": parts[2], "conformer_id": parts[3]})
    elif len(parts) >= 2:
        row.update({"ligand_id": parts[0], "protein_id": parts[1]})
    else:
        legacy = stem.split("_", 1)
        if len(legacy) == 2:
            row.update({"ligand_id": legacy[0], "protein_id": legacy[1]})
        else:
            row.update({"ligand_id": "", "protein_id": stem})
    return row


def _rows_from_logs(config: dict[str, Any]) -> list[dict[str, Any]]:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    docking_out = resolve_path(config.get("paths", {}).get("docking_out_dir", "data/docking_out"), root)
    if docking_out is None or not docking_out.exists():
        return []
    num_aff = int(config.get("affinity_extraction", {}).get("num_affinities", 9) or 9)
    rows: list[dict[str, Any]] = []
    for fp in sorted(docking_out.glob("file_*/*.out")):
        parsed = parse_vina_log(fp, num_affinities=num_aff)
        affs = parsed.get("affinities", [])
        if not affs:
            continue
        base = _parse_job_id_from_log_name(fp)
        rows.append({
            **base,
            "best_affinity": affs[0],
            "affinities": ",".join(map(str, affs)),
            "n_affinities": len(affs),
            "grid_size": parsed.get("grid_size", ""),
            "grid_space": parsed.get("grid_space", ""),
            "exhaustiveness": parsed.get("exhaustiveness", ""),
            "random_seed": parsed.get("random_seed", ""),
            "status": "success",
        })
    return rows


def load_affinity_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    rows = _rows_from_manifest(config)
    if not rows and bool(config.get("affinity_extraction", {}).get("scan_logs_when_manifest_missing", True)):
        rows = _rows_from_logs(config)
    if not rows:
        raise RuntimeError("No docking affinity results found from docking_result_manifest.csv or data/docking_out/file_*/*.out")

    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    protein_manifest = resolve_path(config.get("paths", {}).get("protein_manifest", "data/protein/protein_manifest.csv"), root)
    protein_map = _read_protein_manifest(protein_manifest)
    out: list[dict[str, Any]] = []
    for row in rows:
        pid = row.get("protein_id", "")
        if not row.get("cluster_id") and pid in protein_map:
            row["cluster_id"] = protein_map[pid].get("cluster_id", "")
        row["best_affinity"] = _to_float(row.get("best_affinity"))
        if row["best_affinity"] is not None:
            out.append(row)
    if not out:
        raise RuntimeError("Docking files were found, but no best affinity values could be parsed")
    return out


def reduce_to_best_pair_rows(config: dict[str, Any], rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    reducer = config.get("affinity_extraction", {}).get("duplicate_policy", "best_min_affinity")
    if reducer == "keep_all":
        return rows
    best: dict[tuple[str, str], dict[str, Any]] = {}
    for row in rows:
        key = (str(row.get("ligand_id", "")), str(row.get("protein_id", "")))
        aff = _to_float(row.get("best_affinity"))
        if aff is None:
            continue
        if key not in best or aff < float(best[key]["best_affinity"]):
            best[key] = dict(row)
            best[key]["best_affinity"] = aff
    return sorted(best.values(), key=lambda r: (str(r.get("protein_id", "")), str(r.get("ligand_id", ""))))


def make_affinity_matrix(best_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[str]]:
    ligands = sorted({str(r.get("ligand_id", "")) for r in best_rows if str(r.get("ligand_id", ""))})
    proteins = sorted({str(r.get("protein_id", "")) for r in best_rows if str(r.get("protein_id", ""))})
    meta: dict[str, dict[str, Any]] = {p: {"name": p, "cluster": ""} for p in proteins}
    values: dict[tuple[str, str], float] = {}
    for r in best_rows:
        pid = str(r.get("protein_id", ""))
        lig = str(r.get("ligand_id", ""))
        if not pid or not lig:
            continue
        meta.setdefault(pid, {"name": pid, "cluster": ""})
        if not meta[pid].get("cluster"):
            meta[pid]["cluster"] = r.get("cluster_id", "")
        values[(pid, lig)] = float(r["best_affinity"])
    matrix: list[dict[str, Any]] = []
    for pid in proteins:
        row = dict(meta[pid])
        for ligand in ligands:
            row[ligand] = values.get((pid, ligand), "")
        matrix.append(row)
    return matrix, ligands


def compute_protein_binding_classes(config: dict[str, Any], matrix: list[dict[str, Any]], ligands: list[str]) -> list[dict[str, Any]]:
    scoring = config.get("clusterscore", {})
    strong_threshold = float(scoring.get("strong_threshold", -7.0))
    weak_threshold = float(scoring.get("weak_threshold", -3.0))
    missing_policy = scoring.get("missing_policy", "ignore")
    rows: list[dict[str, Any]] = []
    for row in matrix:
        values: list[float] = []
        for ligand in ligands:
            val = _to_float(row.get(ligand))
            if val is None:
                if missing_policy == "ignore":
                    continue
                val = float(scoring.get("missing_fill_value", 0.0))
            values.append(val)
        total_count = len(values)
        strong_count = sum(1 for v in values if v <= strong_threshold)
        weak_count = sum(1 for v in values if strong_threshold < v <= weak_threshold)
        none_count = sum(1 for v in values if v > weak_threshold)
        rows.append({
            "protein_id": row.get("name", ""),
            "cluster_id": row.get("cluster", ""),
            "total_count": total_count,
            "strong_count": strong_count,
            "weak_count": weak_count,
            "none_count": none_count,
            "strong_ratio": strong_count / total_count if total_count else 0.0,
            "weak_ratio": weak_count / total_count if total_count else 0.0,
            "none_ratio": none_count / total_count if total_count else 0.0,
            "best_affinity_min": min(values) if values else "",
            "mean_affinity": sum(values) / total_count if total_count else "",
        })
    return rows


def compute_cluster_score(config: dict[str, Any], protein_stats: list[dict[str, Any]]) -> list[dict[str, Any]]:
    scoring = config.get("clusterscore", {})
    alpha = float(scoring.get("alpha", 0.2))
    beta = float(scoring.get("beta", 0.2))
    expression = str(scoring.get("formula", "total_strong_count + alpha * total_weak_count + n_proteins * (avg_strong_ratio + beta * avg_weak_ratio)"))
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in protein_stats:
        groups.setdefault(str(row.get("cluster_id", "")), []).append(row)
    rows: list[dict[str, Any]] = []
    for cluster_id, group in groups.items():
        n = len(group)
        strong_ratios = [float(r.get("strong_ratio", 0) or 0) for r in group]
        weak_ratios = [float(r.get("weak_ratio", 0) or 0) for r in group]
        none_ratios = [float(r.get("none_ratio", 0) or 0) for r in group]
        best_vals = [_to_float(r.get("best_affinity_min")) for r in group]
        best_vals = [v for v in best_vals if v is not None]
        mean_vals = [_to_float(r.get("mean_affinity")) for r in group]
        mean_vals = [v for v in mean_vals if v is not None]
        variables = {
            "alpha": alpha,
            "beta": beta,
            "n_proteins": n,
            "avg_strong_ratio": sum(strong_ratios) / n if n else 0.0,
            "avg_weak_ratio": sum(weak_ratios) / n if n else 0.0,
            "avg_none_ratio": sum(none_ratios) / n if n else 0.0,
            "total_strong_count": sum(int(r.get("strong_count", 0) or 0) for r in group),
            "total_weak_count": sum(int(r.get("weak_count", 0) or 0) for r in group),
            "total_none_count": sum(int(r.get("none_count", 0) or 0) for r in group),
            "total_binding_count": sum(int(r.get("total_count", 0) or 0) for r in group),
            "best_affinity_min": min(best_vals) if best_vals else 0.0,
            "mean_affinity": sum(mean_vals) / len(mean_vals) if mean_vals else 0.0,
        }
        composite_score = safe_eval_formula(expression, variables)
        rows.append({"cluster": cluster_id, **variables, "composite_score": composite_score, "formula": expression})
    sort_by = scoring.get("sort_by", "composite_score")
    ascending = bool(scoring.get("sort_ascending", False))
    rows.sort(key=lambda r: r.get(sort_by, 0), reverse=not ascending)
    for i, row in enumerate(rows, start=1):
        row["rank"] = i
    # Put rank first.
    return [{"rank": r.pop("rank"), **r} for r in rows]


def _excel_col(n: int) -> str:
    s = ""
    while n:
        n, rem = divmod(n - 1, 26)
        s = chr(65 + rem) + s
    return s


def _sheet_xml(rows: list[dict[str, Any]], fieldnames: list[str]) -> str:
    xml_rows = []
    header = {name: name for name in fieldnames}
    all_rows = [header] + rows
    for r_idx, row in enumerate(all_rows, start=1):
        cells = []
        for c_idx, key in enumerate(fieldnames, start=1):
            value = row.get(key, "")
            ref = f"{_excel_col(c_idx)}{r_idx}"
            num = _to_float(value)
            if value != "" and num is not None and not isinstance(value, bool) and r_idx != 1:
                cells.append(f'<c r="{ref}"><v>{num}</v></c>')
            else:
                text = escape("" if value is None else str(value))
                cells.append(f'<c r="{ref}" t="inlineStr"><is><t>{text}</t></is></c>')
        xml_rows.append(f'<row r="{r_idx}">' + "".join(cells) + "</row>")
    return '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n' + \
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">' + \
        '<sheetData>' + "".join(xml_rows) + '</sheetData></worksheet>'


def _safe_sheet_name(name: str, used: set[str]) -> str:
    cleaned = re.sub(r"[\\/*?:\[\]]", "_", name)[:31] or "Sheet"
    base = cleaned
    i = 1
    while cleaned in used:
        suffix = f"_{i}"
        cleaned = (base[:31 - len(suffix)] + suffix)[:31]
        i += 1
    used.add(cleaned)
    return cleaned


def _fieldnames(rows: list[dict[str, Any]]) -> list[str]:
    keys: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in keys:
                keys.append(key)
    return keys


def write_xlsx(path: Path, sheets: dict[str, list[dict[str, Any]]]) -> Path:
    ensure_dir(path.parent)
    sheet_names: list[str] = []
    used: set[str] = set()
    normalized: list[tuple[str, list[dict[str, Any]], list[str]]] = []
    for name, rows in sheets.items():
        sheet_name = _safe_sheet_name(name, used)
        sheet_names.append(sheet_name)
        normalized.append((sheet_name, rows, _fieldnames(rows)))
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
''' + "".join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>\n' for i in range(1, len(normalized)+1)) + "</Types>")
        zf.writestr("_rels/.rels", '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>''')
        sheets_xml = "".join(f'<sheet name="{escape(name)}" sheetId="{i}" r:id="rId{i}"/>' for i, (name, _, _) in enumerate(normalized, start=1))
        zf.writestr("xl/workbook.xml", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n' +
                    '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>' + sheets_xml + '</sheets></workbook>')
        rels_xml = "".join(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1, len(normalized)+1))
        zf.writestr("xl/_rels/workbook.xml.rels", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n' +
                    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' + rels_xml + '</Relationships>')
        for i, (_name, rows, fields) in enumerate(normalized, start=1):
            zf.writestr(f"xl/worksheets/sheet{i}.xml", _sheet_xml(rows, fields))
    return path


def write_optional_plots(config: dict[str, Any], out_dir: Path, cluster_rows: list[dict[str, Any]]) -> list[str]:
    plot_cfg = config.get("plots", {})
    if not bool(plot_cfg.get("enabled", False)) or not cluster_rows:
        return []
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as exc:
        write_json(out_dir / "plot_warning.json", {"warning": f"matplotlib not available: {exc}"})
        return []
    top_n = int(plot_cfg.get("top_n", 20) or 20)
    top = cluster_rows[:top_n]
    fig_dir = ensure_dir(out_dir / "figures")
    plt.figure(figsize=(max(8, min(18, 0.5 * len(top))), 6), dpi=int(plot_cfg.get("dpi", 300) or 300))
    plt.bar([str(r.get("cluster", "")) for r in top], [float(r.get("composite_score", 0) or 0) for r in top])
    plt.xlabel("Cluster")
    plt.ylabel("ClusterScore")
    plt.title(f"Top {len(top)} Clusters by ClusterScore")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    pdf_path = fig_dir / f"top{top_n}_clusterscore.pdf"
    png_path = fig_dir / f"top{top_n}_clusterscore.png"
    plt.savefig(pdf_path)
    plt.savefig(png_path)
    plt.close()
    return [str(pdf_path), str(png_path)]


def run_clusterscore_core(config: dict[str, Any]) -> Path:
    root = Path(config.get("paths", {}).get("aimd_root", ".")).resolve()
    out_dir = resolve_path(config.get("paths", {}).get("clusterscore_out_dir", "data/scoring/ClusterScore"), root)
    assert out_dir is not None
    ensure_dir(out_dir)

    raw_rows = load_affinity_rows(config)
    best_rows = reduce_to_best_pair_rows(config, raw_rows)
    matrix, ligands = make_affinity_matrix(best_rows)
    protein_stats = compute_protein_binding_classes(config, matrix, ligands)
    cluster_rows = compute_cluster_score(config, protein_stats)
    top_n = int(config.get("clusterscore", {}).get("top_n", 20) or 20)
    top_rows = cluster_rows[:top_n]

    write_csv(out_dir / "best_affinity_long.csv", best_rows)
    write_csv(out_dir / "protein_ligand_matrix.csv", matrix, fieldnames=["name", "cluster"] + ligands)
    write_csv(out_dir / "protein_binding_statistics.csv", protein_stats)
    write_csv(out_dir / "cluster_binding_statistics.csv", cluster_rows)
    write_csv(out_dir / f"top{top_n}_clusters.csv", top_rows)

    excel_path = out_dir / config.get("output", {}).get("excel_name", "clusterscore_results.xlsx")
    write_xlsx(excel_path, {
        "best_affinity_long": best_rows,
        "protein_ligand_matrix": matrix,
        "protein_binding_stats": protein_stats,
        "cluster_score_all": cluster_rows,
        f"top{top_n}_clusters": top_rows,
    })
    plot_paths = write_optional_plots(config, out_dir, cluster_rows)

    report = {
        "n_raw_docking_rows": len(raw_rows),
        "n_best_ligand_protein_pairs": len(best_rows),
        "n_proteins": len(matrix),
        "n_ligands": len(ligands),
        "n_clusters": len(cluster_rows),
        "excel": str(excel_path),
        "top_clusters_csv": str(out_dir / f"top{top_n}_clusters.csv"),
        "all_clusters_csv": str(out_dir / "cluster_binding_statistics.csv"),
        "plots": plot_paths,
    }
    write_json(out_dir / "clusterscore_report.json", report)
    return excel_path
