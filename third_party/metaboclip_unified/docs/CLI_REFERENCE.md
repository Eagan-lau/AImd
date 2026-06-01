# CLI Reference

## Validate mechanism YAML

```bash
metaboclip validate-mechanism \
  --mechanism /path/to/mechanism.yaml
```

Checks references in `ligand_sites`, `protein_roles`, `geometry_refs`, and `features`.

## Run generic YAML backend

```bash
metaboclip run \
  --mechanism /path/to/mechanism.yaml \
  --profile /path/to/default_profile.yaml \
  --protein-dir /path/to/protein_pdbqt \
  --docking-dir /path/to/docking_results/file_1 \
  --role-table-dir /path/to/ligand_roles \
  --out-dir /path/to/results
```

## Run built-in curated family backend

```bash
metaboclip run-curated \
  --family act \
  --profile /path/to/default_profile.yaml \
  --protein-dir /path/to/protein_pdbqt \
  --docking-dir /path/to/docking_results/file_1 \
  --role-table-dir /path/to/ligand_roles \
  --out-dir /path/to/results/generic_act
```

Supported built-in families:

```text
ugt
act
cyp450
fe2og
ach
```

## Run one protein-ligand pair

```bash
metaboclip run-single \
  --mechanism /path/to/mechanism.yaml \
  --profile /path/to/default_profile.yaml \
  --protein /path/to/protein_conf.pdbqt \
  --docked-pdbqt /path/to/ligand@protein_conf.pdbqt \
  --role-table /path/to/ligand.role_table.csv \
  --out-dir /path/to/single_result
```

Use this for debugging.

## Export one best candidate to PyMOL

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

pymol -cq /path/to/candidate_view.pml
```

## Export passing poses to separate PSE files

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

Then execute generated scripts:

```bash
for pml in /path/to/pymol_candidates/*.pml; do
  pymol -cq "$pml"
done
```

Modes:

```text
--mode pose       One best candidate per pose.
--mode candidate  Every passing candidate.
```

Limit export count:

```bash
--top-n 20
```

## Paper-locked backend

Prepare original scripts:

```bash
metaboclip prepare-paper-locked \
  --archive /path/to/original_scoring_modules.zip \
  --dest /path/to/paper_locked_originals \
  --overwrite \
  --report paper_locked_prepare_report.json
```

Verify checksums:

```bash
metaboclip verify-paper-locked \
  --original-root /path/to/paper_locked_originals \
  --report verify_all.json
```

Run external original scripts:

```bash
metaboclip run-paper-locked \
  --family act \
  --original-root /path/to/paper_locked_originals \
  --stage both \
  --file-range 1-100 \
  --out-dir /path/to/paper_locked_logs/act \
  --report paper_locked_act_run.json
```
