# PyMOL Batch Export

Use `export-pymol-batch` to export all passing poses or candidates as separate PML/PSE files.

## Recommended mode

Default recommended mode:

```text
--mode pose
```

This exports one best candidate per passing pose.

If a pose has multiple passing site sets, only the best candidate for that pose is exported.

## Command

```bash
metaboclip export-pymol-batch \
  --protein /path/to/protein_conf.pdbqt \
  --docked-pdbqt /path/to/ligand@protein_conf.pdbqt \
  --candidate-scores /path/to/candidate_scores.csv \
  --geometry-features /path/to/geometry_features.csv \
  --mechanism /path/to/mechanism.yaml \
  --resolved-ligand-sites /path/to/resolved_ligand_sites.csv \
  --resolved-protein-roles /path/to/resolved_protein_roles.csv \
  --out-dir /path/to/pymol_candidates \
  --mode pose
```

Limit output count:

```bash
--top-n 20
```

Export every passing candidate instead of one best per pose:

```bash
--mode candidate
```

## Generate PSE files

The batch command writes multiple PML files. Run PyMOL on them:

```bash
for pml in /path/to/pymol_candidates/*.pml; do
  pymol -cq "$pml"
done
```

## Output naming

Files are named with rank, pose ID, site set ID, and score:

```text
top_001_pose_3_site_pose_3_set_2_score_93.72.pml
top_001_pose_3_site_pose_3_set_2_score_93.72.pse

top_002_pose_1_site_pose_1_set_5_score_88.14.pml
top_002_pose_1_site_pose_1_set_5_score_88.14.pse
```

Meaning:

```text
top_001      candidate rank among exported poses
pose_3       pose ID in the docking PDBQT
site_...     protein role path / site set
score_93.72  CandidateScore
```

## Export report

The output directory includes:

```text
pymol_export_report.csv
```

Columns include:

```text
rank
pose_id
site_set_id
candidate_score
pml
pse
```

## Validate all PSE files

```bash
for pse in /path/to/pymol_candidates/*.pse; do
  pymol -cq -d "load $pse; print('$pse', cmd.count_atoms('all')); quit"
done
```

Each valid PSE should have atom count greater than 0.
