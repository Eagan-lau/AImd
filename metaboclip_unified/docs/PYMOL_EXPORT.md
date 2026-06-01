# PyMOL Single-Candidate Export

PyMOL export is a post-processing step. It does not participate in MetaboClip core calculation.

## Purpose

`export-pymol` exports one selected passing candidate to PML/PSE.

It uses calculation outputs to save:

```text
selected ligand pose
ligand atom_label selections
protein role selections
candidate atom selections
distance objects
angle objects
```

## Command

```bash
metaboclip export-pymol \
  --protein /path/to/protein_conf.pdbqt \
  --docked-pdbqt /path/to/ligand@protein_conf.pdbqt \
  --candidate-scores /path/to/candidate_scores.csv \
  --geometry-features /path/to/geometry_features.csv \
  --mechanism /path/to/mechanism.yaml \
  --resolved-ligand-sites /path/to/resolved_ligand_sites.csv \
  --resolved-protein-roles /path/to/resolved_protein_roles.csv \
  --out-pml /path/to/candidate_view.pml \
  --save-pse /path/to/candidate_view.pse
```

Then run PyMOL:

```bash
pymol -cq /path/to/candidate_view.pml
```

## What is exported

Default selection:

```text
best candidate by candidate_score
```

Objects and selections saved in the session include:

```text
protein
ligand
LIGAND_SITE_<site>
LIGAND_LABEL_<atom_label>
PROTEIN_ROLE_<role>
CANDIDATE_LIGAND_ATOMS
CANDIDATE_PROTEIN_ATOMS
CANDIDATE_ALL_ATOMS
CANDIDATE_LIGAND_ATOMS_OBJECT
CANDIDATE_PROTEIN_ATOMS_OBJECT
DIST_<feature_name>
ANGLE_<feature_name>
```

## Validate PSE content

```bash
pymol -cq -d "load /path/to/candidate_view.pse; print('all', cmd.count_atoms('all')); print(cmd.get_names()); quit"
```

Expected output:

```text
all > 0
['protein', 'ligand', ...]
```

If `all 0` appears, regenerate the PML with the latest exporter and make sure the PML uses:

```python
cmd.save(_pse_path)
```

not command-line `save path, all, 0, pse`.

## Common issue: PDBQT display

If cartoon is not visible, use sticks:

```pml
hide everything
show sticks, protein
show sticks, ligand
show spheres, CANDIDATE_ALL_ATOMS
zoom CANDIDATE_ALL_ATOMS, 10
```

PDBQT can be loaded as PDB-like input for visualization, but the exporter may internally convert selected PDBQT poses to temporary PDB files for PyMOL.
