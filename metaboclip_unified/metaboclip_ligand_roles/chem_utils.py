from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Set

from rdkit import Chem


def atom_symbol(mol: Chem.Mol, atom_index: int) -> str:
    return mol.GetAtomWithIdx(int(atom_index)).GetSymbol()


def total_h_count(atom: Chem.Atom) -> int:
    return int(atom.GetTotalNumHs(includeNeighbors=True))


def heavy_neighbors(mol: Chem.Mol, atom_index: int) -> List[int]:
    atom = mol.GetAtomWithIdx(int(atom_index))
    return [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]


def first_heavy_neighbor(mol: Chem.Mol, atom_index: int, exclude: Optional[Iterable[int]] = None) -> Optional[int]:
    excluded: Set[int] = set(exclude or [])
    for idx in heavy_neighbors(mol, atom_index):
        if idx not in excluded:
            return idx
    return None


def double_bonded_oxygen_of(mol: Chem.Mol, atom_index: int) -> Optional[int]:
    atom = mol.GetAtomWithIdx(int(atom_index))
    for nbr in atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
        if nbr.GetSymbol() == "O" and bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
            return nbr.GetIdx()
    return None


def classify_hydroxyl(mol: Chem.Mol, oxygen_idx: int, parent_idx: Optional[int]) -> str:
    if parent_idx is None:
        return "unknown_hydroxyl"
    parent = mol.GetAtomWithIdx(parent_idx)
    if parent.GetIsAromatic():
        return "phenolic_hydroxyl"
    if parent.GetSymbol() == "C":
        if double_bonded_oxygen_of(mol, parent_idx) is not None:
            return "carboxylic_acid_hydroxyl"
        return "aliphatic_hydroxyl"
    return "heteroatom_bound_hydroxyl"


def is_carbonyl_carbon(mol: Chem.Mol, atom_index: int) -> bool:
    return double_bonded_oxygen_of(mol, atom_index) is not None


def classify_ch_site(mol: Chem.Mol, carbon_idx: int) -> str:
    atom = mol.GetAtomWithIdx(carbon_idx)
    h_count = total_h_count(atom)
    if atom.GetIsAromatic():
        return "aromatic_ch"
    if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
        return "alkenyl_ch"
    neighbor_indices = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
    if any(mol.GetAtomWithIdx(n).GetIsAromatic() for n in neighbor_indices):
        return "benzylic_ch"
    if any(mol.GetAtomWithIdx(n).GetSymbol() in {"O", "N", "S"} for n in neighbor_indices):
        return "alpha_heteroatom_ch"
    if any(is_carbonyl_carbon(mol, n) for n in neighbor_indices):
        return "alpha_carbonyl_ch"
    if h_count >= 3:
        return "methyl_ch"
    if h_count == 2:
        return "methylene_ch2"
    if h_count == 1:
        return "methine_ch"
    return "carbon_without_h"
