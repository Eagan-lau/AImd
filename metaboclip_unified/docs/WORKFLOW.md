# Workflow

## Overview

MetaboClip separates chemistry, protein role resolution, geometry gating, scoring, and visualization.

```text
Original ligand files
  -> LigandRoleMap
  -> role_table.csv

Docked ligand PDBQT files
  -> LigandPoseParser
  -> pose-specific ligand atom coordinates

Protein PDBQT files
  -> ProteinSiteMap
  -> resolved protein roles

mechanism.yaml
  -> ligand site selection
  -> protein role chain
  -> geometry features
  -> gate and score settings

CatalyticGate
  -> passing_candidates.csv

CatalyticScorer
  -> candidate_scores.csv
  -> pose_scores.csv
  -> conformation_scores.csv
  -> protein_scores.csv

PyMOLExporter
  -> optional PML/PSE files for passing candidates
```

## Step 1: Generate ligand role tables

Single ligand:

```bash
metaboclip-build-role-table \
  --ligand-source /path/to/original_ligands/1.sdf \
  --prepared-pdbqt /path/to/prepared_ligands/1.pdbqt \
  --ligand-id 1 \
  --rules rules/functional_groups.yaml \
  --out-role-table /path/to/ligand_roles/1.role_table.csv \
  --out-annotation /path/to/ligand_roles/1.annotation.json \
  --out-atom-map /path/to/ligand_roles/1.atom_map.json
```

Batch:

```bash
metaboclip-build-role-tables-batch \
  --ligand-source-dir /path/to/original_ligands \
  --prepared-pdbqt-dir /path/to/prepared_ligands \
  --rules rules/functional_groups.yaml \
  --out-dir /path/to/ligand_roles
```

The main file used by MetaboClip is:

```text
*.role_table.csv
```

Required columns:

```text
ligand_id
group_id
instance_id
atom_label
atom_class
atom_role
source_atom_index
element
pdbqt_order
subtype
confidence
```

## Step 2: Inspect ligand atom labels

Inspect role tables before writing mechanism YAML:

```bash
head -n 5 /path/to/ligand_roles/1.role_table.csv
cut -d, -f4,5 /path/to/ligand_roles/1.role_table.csv | sort | uniq
```

Common labels and classes:

```text
hydroxyl.o                      hydroxyl_oxygen
protic_nucleophile.atom         protic_heteroatom
acetyl.acyl_c                   acyl_carbonyl_c
ester.carbonyl_c                ester_carbonyl_c
glycosyl_pyranose.c1            glycosyl_anomeric_c
c_h_site.c                      hydrogen_bearing_carbon
```

## Step 3: Edit mechanism.yaml

Edit only mechanism-relevant sections:

```text
ligand_sites
protein_templates
protein_roles
geometry_refs
features
```

Do not put runtime output paths or family metadata in this file unless needed.

## Step 4: Validate mechanism.yaml

```bash
metaboclip validate-mechanism \
  --mechanism /path/to/mechanism.yaml
```

## Step 5: Run MetaboClip

Generic mode:

```bash
metaboclip run \
  --mechanism /path/to/mechanism.yaml \
  --profile metaboclip/config/profiles/default_profile.yaml \
  --protein-dir /path/to/protein_pdbqt \
  --docking-dir /path/to/docking_results/file_1 \
  --role-table-dir /path/to/ligand_roles \
  --out-dir /path/to/results
```

Curated built-in mode:

```bash
metaboclip run-curated \
  --family act \
  --profile metaboclip/config/profiles/default_profile.yaml \
  --protein-dir /path/to/protein_pdbqt \
  --docking-dir /path/to/docking_results/file_1 \
  --role-table-dir /path/to/ligand_roles \
  --out-dir /path/to/results/generic_act
```

## Step 6: Check outputs

Conformation-level output folders contain:

```text
resolved_ligand_sites.csv
resolved_protein_roles.csv
geometry_features.csv
passing_candidates.csv
candidate_scores.csv
pose_scores.csv
conformation_scores.csv
```

Top-level output directory contains:

```text
merged_conformation_scores.csv
protein_scores.csv
summary.json
```

## Step 7: Export PyMOL PSE files

Single best candidate:

```bash
metaboclip export-pymol ...
pymol -cq candidate_view.pml
```

All passing poses, one best candidate per pose:

```bash
metaboclip export-pymol-batch --mode pose ...
for pml in /path/to/pymol_candidates/*.pml; do pymol -cq "$pml"; done
```
