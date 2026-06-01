from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any
import os
import re
import shutil
import subprocess
import tempfile

import numpy as np

from metaboclip.core.atoms import Atom, read_structure_atoms


class TemplateAlignmentError(RuntimeError):
    pass


def write_atoms_as_pdb(atoms: list[Atom], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for i, atom in enumerate(atoms, start=1):
            name = atom.atom_name[:4].rjust(4)
            resn = atom.resn[:3].rjust(3)
            chain = (atom.chain or "A")[:1]
            try:
                resi_int = int(re.sub(r"[^0-9-]", "", atom.resi) or i)
            except Exception:
                resi_int = i
            element = (atom.element or atom.atom_name[:1]).rjust(2)
            handle.write(
                f"ATOM  {i:5d} {name} {resn} {chain}{resi_int:4d}    "
                f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}  1.00  0.00          {element}\n"
            )
        handle.write("END\n")


def parse_tmalign_matrix(matrix_file: Path) -> tuple[np.ndarray, np.ndarray]:
    rows: list[tuple[float, list[float]]] = []
    with open(matrix_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) >= 5 and parts[0] in {"1", "2", "3"}:
                try:
                    t = float(parts[1])
                    u = [float(parts[2]), float(parts[3]), float(parts[4])]
                    rows.append((t, u))
                except Exception:
                    continue
    if len(rows) != 3:
        raise TemplateAlignmentError(f"Failed to parse TM-align matrix: {matrix_file}")
    translation = np.array([r[0] for r in rows], dtype=float)
    rotation = np.array([r[1] for r in rows], dtype=float)
    return rotation, translation


def parse_tmalign_quality(output_text: str) -> dict[str, Any]:
    out: dict[str, Any] = {}
    aligned = re.search(r"Aligned length=\s*(\d+),\s*RMSD=\s*([-+]?\d+(?:\.\d+)?)", output_text)
    if aligned:
        out["aligned_length"] = int(aligned.group(1))
        out["rmsd"] = float(aligned.group(2))
    scores = re.findall(r"TM-score=\s*([-+]?\d+(?:\.\d+)?)", output_text)
    if scores:
        out["tm_score_1"] = float(scores[0])
        if len(scores) > 1:
            out["tm_score_2"] = float(scores[1])
        out["tm_score"] = max(float(v) for v in scores)
    identity = re.search(r"Seq_ID=n_identical/n_aligned=\s*([-+]?\d+(?:\.\d+)?)", output_text)
    if identity:
        out["seq_identity_aligned"] = float(identity.group(1))
    return out


def transform_atoms(atoms: list[Atom], rotation: np.ndarray, translation: np.ndarray, source: str = "template") -> list[Atom]:
    transformed: list[Atom] = []
    for atom in atoms:
        coord = rotation @ atom.coord + translation
        extra = dict(atom.extra or {})
        extra["template_original_x"] = atom.x
        extra["template_original_y"] = atom.y
        extra["template_original_z"] = atom.z
        transformed.append(
            replace(
                atom,
                x=float(coord[0]),
                y=float(coord[1]),
                z=float(coord[2]),
                source=source,
                extra=extra,
            )
        )
    return transformed


def run_tmalign(template_atoms: list[Atom], target_atoms: list[Atom], executable: str = "TMalign") -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    exe = shutil.which(executable)
    if exe is None:
        raise TemplateAlignmentError(f"TM-align executable not found: {executable}")
    with tempfile.TemporaryDirectory(prefix="metaboclip_tmalign_") as tmp:
        tmpdir = Path(tmp)
        template_pdb = tmpdir / "template.pdb"
        target_pdb = tmpdir / "target.pdb"
        matrix_file = tmpdir / "matrix.txt"
        write_atoms_as_pdb(template_atoms, template_pdb)
        write_atoms_as_pdb(target_atoms, target_pdb)
        cmd = [exe, str(template_pdb), str(target_pdb), "-m", str(matrix_file)]
        proc = subprocess.run(cmd, cwd=tmpdir, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise TemplateAlignmentError(proc.stderr.strip() or proc.stdout.strip() or "TM-align failed")
        rotation, translation = parse_tmalign_matrix(matrix_file)
        quality = parse_tmalign_quality(proc.stdout + "\n" + proc.stderr)
        quality["backend"] = "tmalign"
        return rotation, translation, quality


def quality_passes(quality: dict[str, Any], spec: dict[str, Any]) -> bool:
    cfg = spec.get("quality") or {}
    if cfg.get("min_tm_score") is not None and quality.get("tm_score") is not None:
        if float(quality["tm_score"]) < float(cfg["min_tm_score"]):
            return False
    if cfg.get("max_rmsd") is not None and quality.get("rmsd") is not None:
        if float(quality["rmsd"]) > float(cfg["max_rmsd"]):
            return False
    if cfg.get("min_aligned_length") is not None and quality.get("aligned_length") is not None:
        if int(quality["aligned_length"]) < int(cfg["min_aligned_length"]):
            return False
    return True


def load_template_atoms_with_alignment(
    template_path: Path,
    target_protein_file: Path,
    target_atoms: list[Atom],
    spec: dict[str, Any],
) -> tuple[list[Atom], dict[str, Any]]:
    template_atoms = read_structure_atoms(template_path, source="template")
    align_to_protein = bool(spec.get("align_to_protein", False))
    method = str(spec.get("align_method", "none")).lower()
    backend = str(spec.get("align_backend", spec.get("backend", ""))).lower()
    if not align_to_protein or method in {"none", "identity", "false"}:
        return template_atoms, {"backend": "identity", "aligned": False}
    if method in {"structure_global", "tmalign"} or backend == "tmalign":
        executable = str(spec.get("executable", "TMalign"))
        rotation, translation, quality = run_tmalign(template_atoms, target_atoms, executable=executable)
        if not quality_passes(quality, spec):
            raise TemplateAlignmentError(f"Template alignment quality failed for {template_path}")
        return transform_atoms(template_atoms, rotation, translation, source="template_aligned"), quality
    raise TemplateAlignmentError(f"Unsupported template alignment method: {method}")
