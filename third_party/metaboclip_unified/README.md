# MetaboClip Unified

MetaboClip is a catalytic geometry framework for linking ligand atom labels, protein catalytic roles, docking poses, hard geometry gates, and hierarchical scores.

The framework has three operating modes:

1. Generic YAML mode for new enzyme families.
2. Built-in curated YAML mode for the five current families.
3. Paper-locked external-runner mode for exact manuscript-result reproduction.

## Core modules

```text
LigandRoleMap
  Original ligand chemistry -> stable ligand atom labels and PDBQT atom order mapping.

LigandPoseParser
  Multi-pose ligand PDBQT -> pose-specific atom coordinates by pdbqt_order.

ProteinSiteMap
  ligand.<site>, protein.<role>, or template.<template> -> stepwise protein role resolution.

TemplateAligner
  Structure-first template alignment with TMalign or US-align, followed by local target atom search.

GeometryEngine
  distance, angle_3pt, best_fit_plane, oriented_axis, and axis_deviation.

CatalyticGate
  Hard filtering using required feature thresholds.

CatalyticScorer
  CandidateScore = 100 * (0.30 * AffinityScore + 0.70 * HierarchicalGeometryScore).

PyMOLExporter
  Optional visualization only. It exports selected passing candidates as PML/PSE.

PaperLockedBackend
  External original-script runner for exact reproduction.
```

## Recommended workflow

```text
1. Generate ligand role tables with LigandRoleMap.
2. Inspect key ligand atom labels when needed.
3. Write or edit mechanism.yaml.
4. Run MetaboClip generic or curated YAML mode.
5. Inspect candidate_scores.csv and protein_scores.csv.
6. Export top passing poses as PyMOL PSE files when needed.
7. Use paper-locked mode only for exact original-script reproduction.
```

## Install

### Option A: Conda environment

RDKit is most stable in a conda environment.

```bash
conda create -n metaboclip python=3.10 -y
conda activate metaboclip
conda install -c conda-forge rdkit -y
pip install -e .
```

### Option B: Existing conda environment

If you already have a conda environment, install RDKit and then install MetaboClip in the same environment.

```bash
conda activate aimd
conda install -c conda-forge rdkit -y
pip install -e .
```

### Option C: Python venv

Use this only if RDKit is available inside the venv.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

If editable installation fails with a setuptools flat-layout package discovery error, see:

```text
docs/PACKAGE_DISCOVERY_FIX.md
```

### Optional tools

TMalign or US-align is only required for template-based protein role mapping.

PyMOL is only required for PSE visualization export. The core MetaboClip scoring workflow does not require PyMOL.

## Quick start with example data

Download and unpack the example data package before running the commands below.

Example data archive:
https://pan.baidu.com/s/1dEnCIN6fr6d0am7CyzCk9w?pwd=ah1w
Access code: ah1w

Expected directory layout:

```text
metaboclip-unified/
  data_input/
    ligand_raw/
      1.mol2
    ligand_pdbqt/
      1.pdbqt
    protein_pdbqt/
      Tca01g00352.t1_pLDDT=91.1_pTM=0.91_2.pdbqt
    docking_results/
      file_1/
        1@Tca01g00352.t1_pLDDT=91.1_pTM=0.91_2.pdbqt
  data_output/
```

### 1. Build ligand role table

```bash
metaboclip-build-role-table \
  --ligand-source data_input/ligand_raw/1.mol2 \
  --prepared-pdbqt data_input/ligand_pdbqt/1.pdbqt \
  --ligand-id 1 \
  --rules rules/functional_groups.yaml \
  --out-role-table data_output/ligand_roles/1.role_table.csv \
  --out-annotation data_output/ligand_roles/1.annotation.json \
  --out-atom-map data_output/ligand_roles/1.atom_map.json
```

### 2. Run catalytic geometry scoring

```bash
metaboclip run \
  --mechanism metaboclip/config/families/act/mechanism.yaml \
  --profile metaboclip/config/profiles/default_profile.yaml \
  --protein-dir data_input/protein_pdbqt \
  --docking-dir data_input/docking_results/file_1 \
  --role-table-dir data_output/ligand_roles \
  --out-dir data_output/results/generic_act
```

### 3. Export passing poses to PyMOL sessions

```bash
metaboclip export-pymol-batch \
  --protein data_input/protein_pdbqt/Tca01g00352.t1_pLDDT=91.1_pTM=0.91_2.pdbqt \
  --docked-pdbqt data_input/docking_results/file_1/1@Tca01g00352.t1_pLDDT=91.1_pTM=0.91_2.pdbqt \
  --candidate-scores data_output/results/generic_act/1/Tca01g00352.t1_pLDDT_91.1_pTM_0.91_2/candidate_scores.csv \
  --geometry-features data_output/results/generic_act/1/Tca01g00352.t1_pLDDT_91.1_pTM_0.91_2/geometry_features.csv \
  --mechanism metaboclip/config/families/act/mechanism.yaml \
  --resolved-ligand-sites data_output/results/generic_act/1/Tca01g00352.t1_pLDDT_91.1_pTM_0.91_2/resolved_ligand_sites.csv \
  --resolved-protein-roles data_output/results/generic_act/1/Tca01g00352.t1_pLDDT_91.1_pTM_0.91_2/resolved_protein_roles.csv \
  --out-dir data_output/results/generic_act/1/Tca01g00352.t1_pLDDT_91.1_pTM_0.91_2/pymol_candidates \
  --mode pose
```

Run PyMOL to generate PSE files:

```bash
for pml in data_output/results/generic_act/1/Tca01g00352.t1_pLDDT_91.1_pTM_0.91_2/pymol_candidates/*.pml; do
  pymol -cq "$pml"
done
```

The batch exporter writes one PML/PSE pair per selected passing pose. File names include top rank, pose_id, site_set_id, and candidate_score.

## Run a generic mechanism

```bash
metaboclip run \
  --mechanism examples/new_family/mechanism.yaml \
  --profile metaboclip/config/profiles/default_profile.yaml \
  --protein-dir /path/to/protein_pdbqt \
  --docking-dir /path/to/docking_results/file_1 \
  --role-table-dir /path/to/ligand_role_tables \
  --out-dir /path/to/results/new_family
```

## Run a built-in curated family

```bash
metaboclip run-curated \
  --family act \
  --profile metaboclip/config/profiles/default_profile.yaml \
  --protein-dir /path/to/protein_pdbqt \
  --docking-dir /path/to/docking_results/file_1 \
  --role-table-dir /path/to/ligand_role_tables \
  --out-dir /path/to/results/generic_act
```

## Validate a mechanism file

```bash
metaboclip validate-mechanism \
  --mechanism metaboclip/config/families/act/mechanism.yaml
```

## Export passing poses to PyMOL PSE files

### Export one best candidate

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

### Export all passing poses, one PSE per pose

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

for pml in /path/to/pymol_candidates/*.pml; do
  pymol -cq "$pml"
done
```

Use `--mode pose` to export the best candidate per pose. Use `--mode candidate` to export every passing candidate.

## Paper-locked reproduction mode

Paper-locked mode runs external original scripts without modifying them. Use this mode only when exact manuscript-result reproduction is required.

```bash
metaboclip prepare-paper-locked \
  --archive /path/to/original_scoring_modules.zip \
  --dest paper_locked_originals \
  --overwrite

metaboclip verify-paper-locked \
  --original-root paper_locked_originals

metaboclip run-paper-locked \
  --family act \
  --original-root paper_locked_originals \
  --stage both \
  --file-range 1-100 \
  --out-dir paper_locked_logs/act
```

## Important output files

```text
resolved_ligand_sites.csv      Ligand atom-label coordinates for each pose.
resolved_protein_roles.csv     Protein role atoms resolved by ProteinSiteMap.
geometry_features.csv          Distances, angles, axis deviations, and scores.
passing_candidates.csv         Candidates that passed CatalyticGate.
candidate_scores.csv           Candidate-level scores.
pose_scores.csv                Pose-level scores.
conformation_scores.csv        Conformation-level scores.
merged_conformation_scores.csv Merged conformation table across protein conformations.
protein_scores.csv             Protein-level scores across conformations.
```

## Documentation map

```text
docs/INSTALLATION.md              Installation and environment setup.
docs/WORKFLOW.md                  End-to-end workflow.
docs/CLI_REFERENCE.md             Command reference.
docs/LIGAND_ROLE_MAP.md           Ligand atom-label generation.
docs/PROTEIN_SITE_MAP.md          Protein role resolution.
docs/MECHANISM_YAML_MANUAL.md     mechanism.yaml authoring guide.
docs/SCORING.md                   Gate and scoring logic.
docs/CONFORMATION_AGGREGATION.md  Conformation and protein-level scores.
docs/TEMPLATE_ALIGNMENT.md        TMalign / US-align template mapping.
docs/PYMOL_EXPORT.md              Single-candidate PSE export.
docs/PYMOL_BATCH_EXPORT.md        Batch PSE export for passing poses.
docs/PAPER_LOCKED_BACKEND.md      Exact reproduction using original scripts.
docs/TROUBLESHOOTING.md           Common runtime issues.
docs/PACKAGE_DISCOVERY_FIX.md     Editable install flat-layout fix.
```
