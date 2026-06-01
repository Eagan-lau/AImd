from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Optional

from .models import LigandAnnotation
from .molio import read_molecule
from .resolvers import run_graph_resolvers
from .rules import filter_rules, load_rules
from .smarts_annotator import annotate_smarts_groups


DEFAULT_FRAMEWORK_POLICY = {
    "geometry_atom_policy": "heavy_atoms_only",
    "hydrogen_policy": "hydrogen_is_used_for_recognition_only",
    "pdbqt_policy": "pdbqt_is_used_as_pose_coordinate_container",
}


def annotate_ligand(
    ligand_path: str | Path,
    ligand_id: str,
    rules_path: str | Path,
    groups: Optional[Iterable[str]] = None,
) -> LigandAnnotation:
    mol = read_molecule(ligand_path)
    rules = filter_rules(load_rules(rules_path), groups=groups)
    matches = []
    matches.extend(annotate_smarts_groups(mol, rules))
    matches.extend(run_graph_resolvers(mol, rules))
    matches.sort(key=lambda m: (-m.priority, m.group_id, m.instance_id))
    return LigandAnnotation(
        ligand_id=ligand_id,
        source_file=str(ligand_path),
        index_base=0,
        framework_policy=DEFAULT_FRAMEWORK_POLICY,
        functional_groups=matches,
    )


def write_annotation(annotation: LigandAnnotation, out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(annotation.to_dict(), handle, indent=2, sort_keys=True)
