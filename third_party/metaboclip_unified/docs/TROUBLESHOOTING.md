# Troubleshooting

## Completed generic run but candidates=0

Possible causes:

```text
1. role_table.csv does not contain required atom_label or atom_class.
2. role_table.csv is an old format without atom_label and atom_class.
3. mechanism.yaml ligand_sites do not match actual role table labels.
4. required protein role was not resolved.
5. required cofactor or donor atoms are missing from protein PDBQT.
6. gate thresholds are too strict.
```

Check role table:

```bash
head -n 1 /path/to/ligand_roles/1.role_table.csv
cut -d, -f4,5 /path/to/ligand_roles/1.role_table.csv | sort | uniq
```

Check mechanism ligand sites:

```bash
grep -n -A40 "ligand_sites" /path/to/mechanism.yaml
```

Check protein residue names:

```bash
python - <<'PY'
from pathlib import Path
from collections import Counter
protein_dir = Path('/path/to/protein_pdbqt')
counter = Counter()
for path in protein_dir.glob('*.pdbqt'):
    for line in open(path, errors='ignore'):
        if line.startswith(('ATOM', 'HETATM')):
            counter[line[17:20].strip()] += 1
for resn, count in counter.most_common(100):
    print(resn, count)
PY
```

## ACT candidates=0

ACT often requires ligand nucleophile labels and acyl donor/cofactor atoms.

Check ligand role table for:

```text
hydroxyl.o
protic_nucleophile.atom
hydroxyl_oxygen
protic_heteroatom
```

Check protein PDBQT for donor/cofactor residues such as:

```text
COA
ACO
ACT
A2C
H5L
MLC
```

Example:

```bash
grep -E "COA|ACO|ACT|A2C|H5L|MLC" /path/to/protein_pdbqt/*.pdbqt | head
```

If no donor atoms exist, strict ACT YAML may produce no candidates.

## PSE saves but opens empty

Check PSE content:

```bash
pymol -cq -d "load /path/to/candidate_view.pse; print('all', cmd.count_atoms('all')); print(cmd.get_names()); quit"
```

If output is:

```text
all 0
[]
```

Regenerate PML with the latest exporter. The PML should save with:

```python
cmd.save(_pse_path)
```

## PSE looks empty but has atoms

Run inside PyMOL:

```pml
zoom CANDIDATE_ALL_ATOMS, 10
show spheres, CANDIDATE_ALL_ATOMS
show sticks, ligand
show sticks, protein
```

PDBQT protein structures may not display cartoon well.

## File paths do not match conformation

Make sure the protein, docked PDBQT, and result CSV folder refer to the same conformation.

Example:

```text
protein:        ProteinA_2.pdbqt
docked PDBQT:   ligand@ProteinA_2.pdbqt
result folder:  results/.../ProteinA_2/
```

Do not mix conformation 2 structures with conformation 0 result tables.

## Editable install fails

See:

```text
docs/PACKAGE_DISCOVERY_FIX.md
```
