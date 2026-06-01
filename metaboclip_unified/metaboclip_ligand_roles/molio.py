from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

from rdkit import Chem

SUPPORTED_SOURCE_EXTENSIONS = {".sdf", ".sd", ".mol", ".mol2", ".smi", ".smiles"}


def read_molecule(path: str | Path, sanitize: bool = True) -> Chem.Mol:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".sdf", ".sd"}:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=sanitize)
        mol = next((m for m in supplier if m is not None), None)
    elif suffix == ".mol":
        mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=sanitize)
    elif suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(path), removeHs=False, sanitize=sanitize)
    elif suffix in {".smi", ".smiles"}:
        text = path.read_text(encoding="utf-8").strip().splitlines()
        if not text:
            raise ValueError(f"Empty SMILES file: {path}")
        first = text[0].strip().split()[0]
        mol = Chem.MolFromSmiles(first, sanitize=sanitize)
    else:
        raise ValueError(f"Unsupported ligand source format: {suffix}")

    if mol is None:
        raise ValueError(f"Failed to read ligand source: {path}")
    return mol


def iter_ligand_sources(source_dir: str | Path, source_glob: str = "*") -> Iterable[Path]:
    source_dir = Path(source_dir)
    for path in sorted(source_dir.glob(source_glob)):
        if path.is_file() and path.suffix.lower() in SUPPORTED_SOURCE_EXTENSIONS:
            yield path


def heavy_atom_indices(mol: Chem.Mol) -> List[int]:
    return [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]


def require_3d_coordinates(mol: Chem.Mol, source_file: str | Path) -> None:
    if mol.GetNumConformers() == 0:
        raise ValueError(f"Ligand source has no 3D coordinates: {source_file}")
