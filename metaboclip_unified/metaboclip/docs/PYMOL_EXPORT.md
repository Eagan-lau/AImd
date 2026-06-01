# PyMOL Export

`metaboclip export-pymol` is a visualization-only step. It does not calculate catalytic geometry. It reads MetaboClip result tables, selects one passing candidate, extracts the selected ligand pose, writes persistent PyMOL selections, draws distance and angle objects, and writes a PML script.

The core calculation remains PyMOL-free.

## Required inputs

- `--protein`: protein PDBQT or PDB file for the selected conformation.
- `--docked-pdbqt`: multi-pose docked ligand PDBQT file.
- `--candidate-scores`: candidate score table from MetaboClip.
- `--resolved-ligand-sites`: ligand atom-label coordinates from MetaboClip.
- `--resolved-protein-roles`: protein role atoms from MetaboClip.
- `--geometry-features`: geometry feature table from MetaboClip.
- `--mechanism`: mechanism YAML used for the run.
- `--out-pml`: output PyMOL script.
- `--save-pse`: optional PSE path written into the PML script.

## Example

```bash
metaboclip export-pymol \
  --protein /path/to/protein_conf_2.pdbqt \
  --docked-pdbqt /path/to/ligand@protein_conf_2.pdbqt \
  --candidate-scores /path/to/candidate_scores.csv \
  --geometry-features /path/to/geometry_features.csv \
  --mechanism /path/to/mechanism.yaml \
  --resolved-ligand-sites /path/to/resolved_ligand_sites.csv \
  --resolved-protein-roles /path/to/resolved_protein_roles.csv \
  --out-pml /path/to/candidate_view.pml \
  --save-pse /path/to/candidate_view.pse
```

Then run:

```bash
pymol -cq /path/to/candidate_view.pml
```

## What is saved

The PML creates persistent selections such as:

- `LIGAND_POSE_SELECTED`
- `PROTEIN_STRUCTURE`
- `LIGAND_SITE_<site>`
- `LIGAND_LABEL_<atom_label>`
- `PROTEIN_ROLE_<role>`
- `CANDIDATE_LIGAND_ATOMS`
- `CANDIDATE_PROTEIN_ATOMS`
- `CANDIDATE_ALL_ATOMS`

It also creates object copies:

- `CANDIDATE_LIGAND_ATOMS_OBJECT`
- `CANDIDATE_PROTEIN_ATOMS_OBJECT`

These object copies are useful when the PyMOL selection panel is hidden or when an empty selection is difficult to diagnose.

## Important behavior

The exporter selects the best candidate by `candidate_score` unless `--site-set-id` or `--pose-id` is provided. It extracts only the selected ligand pose from the multi-pose PDBQT file. It does not load all ligand poses.

All paths written to PML are absolute paths.
