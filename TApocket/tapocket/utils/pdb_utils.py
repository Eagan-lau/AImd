from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import math

from tapocket.core.schema import QueryPocketResidue


@dataclass(frozen=True)
class Atom:
    record: str
    atom_name: str
    resn: str
    chain: str
    resi: str
    icode: str
    x: float
    y: float
    z: float
    element: str

    @property
    def residue_key(self) -> tuple[str, str, str, str]:
        return (self.chain, self.resi, self.icode, self.resn)

    @property
    def coord(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


def _safe_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return float("nan")


def parse_pdb_atoms(path: str | Path, include_hydrogen: bool = False, include_hetatm: bool = True) -> list[Atom]:
    atoms: list[Atom] = []
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not (line.startswith("ATOM  ") or (include_hetatm and line.startswith("HETATM"))):
                continue
            atom_name = line[12:16].strip()
            resn = line[17:20].strip()
            chain = line[21:22].strip() or "_"
            resi = line[22:26].strip()
            icode = line[26:27].strip()
            x = _safe_float(line[30:38].strip())
            y = _safe_float(line[38:46].strip())
            z = _safe_float(line[46:54].strip())
            if math.isnan(x) or math.isnan(y) or math.isnan(z):
                continue
            element = line[76:78].strip() if len(line) >= 78 else ""
            if not element:
                # Good enough for filtering hydrogens in most PDB files.
                element = "".join(ch for ch in atom_name if ch.isalpha())[:1].upper()
            if not include_hydrogen and (element.upper() == "H" or atom_name.upper().startswith("H")):
                continue
            atoms.append(Atom(line[0:6].strip(), atom_name, resn, chain, resi, icode, x, y, z, element))
    return atoms


def min_distance_between_atom_sets(atoms_a: list[Atom], atoms_b: list[Atom]) -> float:
    if not atoms_a or not atoms_b:
        return float("inf")
    best2 = float("inf")
    for a in atoms_a:
        ax, ay, az = a.coord
        for b in atoms_b:
            dx = ax - b.x
            dy = ay - b.y
            dz = az - b.z
            d2 = dx * dx + dy * dy + dz * dz
            if d2 < best2:
                best2 = d2
    return math.sqrt(best2)


def count_atom_pairs_within_cutoff(atoms_a: list[Atom], atoms_b: list[Atom], cutoff: float) -> int:
    if not atoms_a or not atoms_b:
        return 0
    cutoff2 = cutoff * cutoff
    count = 0
    for a in atoms_a:
        ax, ay, az = a.coord
        for b in atoms_b:
            dx = ax - b.x
            dy = ay - b.y
            dz = az - b.z
            if dx * dx + dy * dy + dz * dz <= cutoff2:
                count += 1
    return count


def group_atoms_by_residue(atoms: Iterable[Atom]) -> dict[tuple[str, str, str, str], list[Atom]]:
    grouped: dict[tuple[str, str, str, str], list[Atom]] = {}
    for atom in atoms:
        grouped.setdefault(atom.residue_key, []).append(atom)
    return grouped



def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def count_unique_residues(atoms: Iterable[Atom]) -> int:
    return len({atom.residue_key for atom in atoms})


def extract_ca_atoms(atoms: Iterable[Atom]) -> list[Atom]:
    return [atom for atom in atoms if atom.atom_name.strip().upper() == "CA"]


def extract_query_residues_by_ca_distance(
    query_pdb: str | Path,
    mapped_pocket_pdb: str | Path,
    query_id: str,
    template_id: str,
    template_pocket_id: str,
    cutoff_angstrom: float = 3.0,
) -> tuple[list[QueryPocketResidue], dict[str, float | int | str | bool]]:
    """Find query pocket residues using the original CA-to-CA logic.

    For each mapped template-pocket CA atom, this finds the nearest query CA atom.
    If the nearest query CA is within cutoff_angstrom, the corresponding query residue
    is considered a spatially mapped pocket residue. Mapping coverage is calculated as
    matched template-pocket CA atoms divided by all template-pocket CA atoms.
    """
    query_atoms = parse_pdb_atoms(query_pdb, include_hydrogen=False, include_hetatm=False)
    pocket_atoms = parse_pdb_atoms(mapped_pocket_pdb, include_hydrogen=False, include_hetatm=True)
    query_ca_atoms = extract_ca_atoms(query_atoms)
    pocket_ca_atoms = extract_ca_atoms(pocket_atoms)

    best_by_query_residue: dict[tuple[str, str, str, str], float] = {}
    matched_template_ca_count = 0
    distance_values: list[float] = []

    for pocket_ca in pocket_ca_atoms:
        best_atom = None
        best_dist = float("inf")
        for query_ca in query_ca_atoms:
            dist = _distance(pocket_ca.coord, query_ca.coord)
            if dist < best_dist:
                best_dist = dist
                best_atom = query_ca
        if best_atom is not None and best_dist <= cutoff_angstrom:
            matched_template_ca_count += 1
            distance_values.append(best_dist)
            key = best_atom.residue_key
            if key not in best_by_query_residue or best_dist < best_by_query_residue[key]:
                best_by_query_residue[key] = best_dist

    residues: list[QueryPocketResidue] = []
    for (chain, resi, icode, resn), dist in best_by_query_residue.items():
        residues.append(
            QueryPocketResidue(
                query_id=query_id,
                template_id=template_id,
                template_pocket_id=template_pocket_id,
                chain=chain,
                resi=resi,
                icode=icode,
                resn=resn,
                min_distance_to_mapped_pocket=round(dist, 4),
            )
        )
    residues.sort(key=lambda r: (r.chain, int(r.resi) if str(r.resi).lstrip('-').isdigit() else 10**9, r.icode, r.resn))

    template_ca_count = len(pocket_ca_atoms)
    coverage = (matched_template_ca_count / template_ca_count) if template_ca_count else 0.0
    mean_distance = (sum(distance_values) / len(distance_values)) if distance_values else float("inf")
    sorted_distances = sorted(distance_values)
    if sorted_distances:
        mid = len(sorted_distances) // 2
        if len(sorted_distances) % 2:
            median_distance = sorted_distances[mid]
        else:
            median_distance = (sorted_distances[mid - 1] + sorted_distances[mid]) / 2.0
    else:
        median_distance = float("inf")

    quality = {
        "method": "ca_distance",
        "ca_distance_cutoff_angstrom": cutoff_angstrom,
        "template_pocket_ca_count": template_ca_count,
        "matched_template_ca_count": matched_template_ca_count,
        "query_residue_count": len(residues),
        "mapping_coverage": round(coverage, 6),
        "mean_ca_distance_angstrom": None if mean_distance == float("inf") else round(mean_distance, 4),
        "median_ca_distance_angstrom": None if median_distance == float("inf") else round(median_distance, 4),
    }
    return residues, quality


def extract_query_residues_by_heavy_atom_distance(
    query_pdb: str | Path,
    mapped_pocket_pdb: str | Path,
    query_id: str,
    template_id: str,
    template_pocket_id: str,
    cutoff_angstrom: float = 4.5,
    include_hydrogen: bool = False,
) -> tuple[list[QueryPocketResidue], dict[str, float | int | str | bool]]:
    """Find query pocket residues by heavy-atom distance and return mapping quality."""
    query_atoms = parse_pdb_atoms(query_pdb, include_hydrogen=include_hydrogen, include_hetatm=False)
    pocket_atoms = parse_pdb_atoms(mapped_pocket_pdb, include_hydrogen=include_hydrogen, include_hetatm=True)
    residues_by_key = group_atoms_by_residue(query_atoms)

    out: list[QueryPocketResidue] = []
    distance_values: list[float] = []
    for (chain, resi, icode, resn), atoms in residues_by_key.items():
        min_dist = min_distance_between_atom_sets(atoms, pocket_atoms)
        if min_dist <= cutoff_angstrom:
            distance_values.append(min_dist)
            out.append(
                QueryPocketResidue(
                    query_id=query_id,
                    template_id=template_id,
                    template_pocket_id=template_pocket_id,
                    chain=chain,
                    resi=resi,
                    icode=icode,
                    resn=resn,
                    min_distance_to_mapped_pocket=round(min_dist, 4),
                )
            )
    out.sort(key=lambda r: (r.chain, int(r.resi) if str(r.resi).lstrip('-').isdigit() else 10**9, r.icode, r.resn))

    template_residue_count = count_unique_residues(pocket_atoms)
    matched_count = len(out)
    coverage = (matched_count / template_residue_count) if template_residue_count else 0.0
    sorted_distances = sorted(distance_values)
    median_distance = None
    if sorted_distances:
        mid = len(sorted_distances) // 2
        median_distance = sorted_distances[mid] if len(sorted_distances) % 2 else (sorted_distances[mid - 1] + sorted_distances[mid]) / 2.0

    quality = {
        "method": "heavy_atom_distance",
        "distance_cutoff_angstrom": cutoff_angstrom,
        "template_pocket_residue_count": template_residue_count,
        "matched_template_residue_count": matched_count,
        "query_residue_count": len(out),
        "mapping_coverage": round(coverage, 6),
        "median_min_distance_angstrom": None if median_distance is None else round(median_distance, 4),
    }
    return out, quality


def atom_box(atoms: list[Atom], padding: float = 0.0) -> dict[str, float]:
    if not atoms:
        return {
            "xmin": float("inf"), "xmax": float("-inf"),
            "ymin": float("inf"), "ymax": float("-inf"),
            "zmin": float("inf"), "zmax": float("-inf"),
            "padding": padding,
        }
    xs = [a.x for a in atoms]
    ys = [a.y for a in atoms]
    zs = [a.z for a in atoms]
    return {
        "xmin": min(xs) - padding,
        "xmax": max(xs) + padding,
        "ymin": min(ys) - padding,
        "ymax": max(ys) + padding,
        "zmin": min(zs) - padding,
        "zmax": max(zs) + padding,
        "padding": padding,
    }


def atom_in_box(atom: Atom, box: dict[str, float]) -> bool:
    return (
        box["xmin"] <= atom.x <= box["xmax"]
        and box["ymin"] <= atom.y <= box["ymax"]
        and box["zmin"] <= atom.z <= box["zmax"]
    )


def count_atoms_in_box(atoms: list[Atom], box: dict[str, float]) -> int:
    return sum(1 for atom in atoms if atom_in_box(atom, box))


def box_center(box: dict[str, float]) -> tuple[float, float, float]:
    return ((box["xmin"] + box["xmax"]) / 2.0, (box["ymin"] + box["ymax"]) / 2.0, (box["zmin"] + box["zmax"]) / 2.0)


def box_size(box: dict[str, float]) -> tuple[float, float, float]:
    return (box["xmax"] - box["xmin"], box["ymax"] - box["ymin"], box["zmax"] - box["zmin"])

def extract_query_residues_near_pocket(
    query_pdb: str | Path,
    mapped_pocket_pdb: str | Path,
    query_id: str,
    template_id: str,
    template_pocket_id: str,
    cutoff_angstrom: float = 4.5,
    include_hydrogen: bool = False,
) -> list[QueryPocketResidue]:
    query_atoms = parse_pdb_atoms(query_pdb, include_hydrogen=include_hydrogen, include_hetatm=False)
    pocket_atoms = parse_pdb_atoms(mapped_pocket_pdb, include_hydrogen=include_hydrogen, include_hetatm=True)
    residues = group_atoms_by_residue(query_atoms)

    out: list[QueryPocketResidue] = []
    for (chain, resi, icode, resn), atoms in residues.items():
        min_dist = min_distance_between_atom_sets(atoms, pocket_atoms)
        if min_dist <= cutoff_angstrom:
            out.append(
                QueryPocketResidue(
                    query_id=query_id,
                    template_id=template_id,
                    template_pocket_id=template_pocket_id,
                    chain=chain,
                    resi=resi,
                    icode=icode,
                    resn=resn,
                    min_distance_to_mapped_pocket=round(min_dist, 4),
                )
            )
    out.sort(key=lambda r: (r.chain, int(r.resi) if str(r.resi).lstrip('-').isdigit() else 10**9, r.icode, r.resn))
    return out


def write_multimodel_pdb(input_paths: list[str | Path], output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as out:
        for idx, pdb_path in enumerate(input_paths, start=1):
            out.write(f"MODEL     {idx}\n")
            with Path(pdb_path).open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if line.startswith(("ATOM  ", "HETATM", "TER")):
                        out.write(line)
            out.write("ENDMDL\n")
        out.write("END\n")
