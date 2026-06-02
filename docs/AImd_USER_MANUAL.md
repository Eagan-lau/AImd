# AImd User Manual

## Overview

AImd connects protein clustering, pocket prediction, docking, cluster scoring, refinement, and catalytic scoring through manifest files.

Current workflow:

```text
RGPC -> TApocketBridge -> DockingHub -> ClusterScore -> RefinementHub -> refined DockingHub -> MetaBoClipBridge
```

MetaboClip is a core scientific component of AImd. `MetaBoClipBridge` is the active adapter for the unified MetaboClip backend under `metaboclip_unified`; the old implementation is not part of the clean deliverable package.

## Data Layout

AImd separates starting inputs from generated workflow outputs:

```text
data/
  data_input/
    protein/
    ligand/
    cofactor/
    workflow/
  data_output/
    cluster/
    protein_batches/
    ligand_transformation/
    pocket/
    ensemble/
    cofactor_mapped/
    receptor/
    docking_configs/
    docking_tasks/
    docking_out/
    scoring/
    refinement/
    refined/
    metaboclip/
```

`data_input` is reserved for user-provided starting files and manual workflow control tables. `data_output` is reserved for module-generated files.

Prepare or migrate the local data layout before running modules:

```bash
python scripts/migrate_data_layout.py --root .
```

The helper copies legacy local input folders into the canonical layout and writes project-root-relative protein and ligand input manifests. It does not delete legacy folders.

## Environment

Install Python dependencies:

```bash
pip install -r requirements.txt
```

Core unified MetaboClip scoring requires Python imports for PyYAML, NumPy, SciPy, and RDKit. PyMOL is not required for core scoring; it is optional for visualization/export features or other modules that explicitly call it.

Validate the repository layout:

```bash
python validate_aimd_layout.py --root .
```

## Main Commands

Run modules one step at a time:

```bash
python run_rgpc.py --config configs/RGPC/rgpc.yaml
python run_tapocket_batch.py --config configs/TApocket/tapocket_batch.yaml
python run_docking.py --config configs/Docking/docking.yaml
python run_clusterscore.py --config configs/Scoring/cluster_score.yaml
python run_refinement.py --config configs/Refinement/refine_from_clusterscore.yaml
python run_metaboclip_bridge.py --config configs/MetaBoClip/metaboclip_bridge.yaml
```

Run the full workflow only after all inputs are prepared:

```bash
python run_full_iterative_metaboclip.py --config configs/workflows/full_iterative_metaboclip.yaml
```

Do not use full docking or production workflows as smoke tests.

## Required Data Flow

The refined docking manifest consumed by `MetaBoClipBridge` is:

```text
data/data_output/refined/docking_out/docking_result_manifest.csv
```

Important columns include:

```csv
job_id,ligand_id,protein_id,cluster_id,batch_id,conformer_id,pocket_id,pocket_rank,receptor_pdbqt_path,ligand_pdbqt_path,config_path,out_pose_path,log_path,center_x,center_y,center_z,size_x,size_y,size_z,status,return_code,message,best_affinity,affinities,n_affinities,grid_size,grid_space,exhaustiveness,random_seed,pose_exists
```

The stable final ranking output is:

```text
data/data_output/metaboclip/results/metaboclip_final_ranking.csv
```

## Unified MetaboClip Configuration

The active bridge config is:

```text
configs/MetaBoClip/metaboclip_bridge.yaml
```

Core fields:

```yaml
backend: unified
paths:
  metaboclip_project_dir: metaboclip_unified
  metaboclip_profile: metaboclip_unified/metaboclip/config/profiles/default_profile.yaml
  unified_output_dir: data/data_output/metaboclip/unified_runs
  role_table_dir: data/data_output/metaboclip/ligand_roles/role_tables
  annotation_dir: data/data_output/metaboclip/ligand_roles/annotations
  atom_map_dir: data/data_output/metaboclip/ligand_roles/atom_maps
```

Family mechanisms are configured under `mechanisms`. The bridge calls the Python API `metaboclip.core.workflow.run_directory` for scoring. `run_single_pair` is detected and available, but directory execution is preferred because it preserves unified backend protein aggregation.

## Role Tables

Role-table behavior is controlled by:

```yaml
role_tables:
  mode: existing
```

Allowed modes:

```text
existing
generate
auto
```

`existing` requires prebuilt role tables under `role_table_dir`. `generate` requires an original ligand source file and prepared ligand PDBQT file. `auto` reuses existing role tables and generates missing ones when the required source inputs exist.

Generated or reused role assets are recorded as:

```text
role_table_path
annotation_json_path
atom_map_json_path
```

The atom map and downstream catalytic geometry use heavy atoms only.

## Outputs

Unified backend artifacts:

```text
data/data_output/metaboclip/unified_runs/
```

AImd-compatible aggregate outputs:

```text
data/data_output/metaboclip/results/metaboclip_protein_scores_all.csv
data/data_output/metaboclip/results/metaboclip_conformation_scores_all.csv
data/data_output/metaboclip/results/metaboclip_pose_scores_all.csv
data/data_output/metaboclip/results/metaboclip_candidate_scores_all.csv
data/data_output/metaboclip/results/metaboclip_passing_candidates_all.csv
data/data_output/metaboclip/results/metaboclip_geometry_features_all.csv
data/data_output/metaboclip/results/metaboclip_resolved_ligand_sites_all.csv
data/data_output/metaboclip/results/metaboclip_resolved_protein_roles_all.csv
data/data_output/metaboclip/results/metaboclip_run_manifest.csv
data/data_output/metaboclip/results/metaboclip_report.json
```

The final ranking preserves AImd metadata and uses real unified score fields when available:

```csv
protein_score,quality_score,coverage
```

The bridge does not fabricate unavailable legacy scientific columns:

```csv
protein_score_norm,max_s_r
```

## Minimal Smoke Test Guidance

The repository includes a deterministic MetaboClip bridge smoke test that uses tiny unified-backend fixtures and writes only to pytest temporary directories:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:metaboclip_unified \
  python -m pytest tests/test_metaboclip_bridge_smoke.py -q
```

A smoke test should use tiny inputs only. It should verify:

1. `MetaBoClipBridge` imports.
2. `metaboclip_unified` can be located.
3. mechanism YAML and profile YAML paths resolve.
4. role-table paths, annotation JSON paths, and atom-map JSON paths are handled.
5. the bridge can write an AImd-compatible final ranking or dry-run report.

Do not run full docking or a full production AImd workflow for smoke testing.
