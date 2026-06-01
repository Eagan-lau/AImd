from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from tapocket.core.schema import CandidatePocket, QueryPocketResidue
from tapocket.utils.fs import rel_or_abs


def _residue_from_dict(query_id: str, model_name: str, pocket_id: str, row: dict[str, Any]) -> QueryPocketResidue:
    return QueryPocketResidue(
        query_id=query_id,
        template_id=model_name,
        template_pocket_id=pocket_id,
        chain=str(row.get("chain", "_") or "_"),
        resi=str(row.get("resi", row.get("residue_number", ""))),
        icode=str(row.get("icode", "") or ""),
        resn=str(row.get("resn", row.get("residue_name", "")) or ""),
        min_distance_to_mapped_pocket=float(row.get("distance", row.get("min_distance_to_mapped_pocket", 0.0)) or 0.0),
    )


def parse_tapocket_json(
    output_json: str | Path,
    query_id: str,
    model_name: str,
    run_root: str | Path,
    output_dir: str | Path,
) -> list[CandidatePocket]:
    """Parse the standard AI output JSON.

    Expected schema:
      {"pockets": [{"pocket_id": "ai_pocket_001", "confidence": 0.9,
                    "residues": [{"chain":"A","resi":"57","resn":"HIS"}],
                    "box": {"center": [...], "size": [...]},
                    "pocket_pdb": "optional/path.pdb"}]}
    """
    output_json = Path(output_json)
    if not output_json.exists():
        raise FileNotFoundError(output_json)
    data = json.loads(output_json.read_text(encoding="utf-8"))
    pockets = data.get("pockets", []) if isinstance(data, dict) else data
    out: list[CandidatePocket] = []
    for idx, item in enumerate(pockets, start=1):
        pocket_id = str(item.get("pocket_id") or f"ai_pocket_{idx:03d}")
        residues = [_residue_from_dict(query_id, model_name, pocket_id, r) for r in item.get("residues", [])]
        pocket_pdb = item.get("pocket_pdb") or item.get("pocket_path") or ""
        if pocket_pdb:
            p = Path(str(pocket_pdb))
            if not p.is_absolute():
                p = Path(output_dir) / p
            pocket_pdb = rel_or_abs(p, run_root)
        confidence = float(item.get("confidence", item.get("score", 0.0)) or 0.0)
        cand = CandidatePocket(
            query_id=query_id,
            template_id=model_name,
            template_pocket_id=pocket_id,
            custom_pocket_id=pocket_id,
            source=model_name,
            mapped_pocket_path=str(pocket_pdb),
            query_pocket_residues=residues,
            mapping_quality={"method": "ai_model", "query_residue_count": len(residues)},
            ai_support={"model_name": model_name, "confidence": confidence, "raw": item},
            box=item.get("box", {}) or {},
            final_score=confidence,
            extra={"ai_raw": item},
        )
        out.append(cand)
    return out


def locate_ai_output(output_dir: str | Path, configured: str | None = None) -> Path:
    output_dir = Path(output_dir)
    if configured:
        candidate = Path(configured)
        if not candidate.is_absolute():
            candidate = output_dir / candidate
        if candidate.exists():
            return candidate
        raise FileNotFoundError(candidate)
    for name in ["tapocket_ai_pockets.json", "final_pockets.json", "pockets.json", "predicted_pockets.json"]:
        candidate = output_dir / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"No AI output JSON found in {output_dir}")
