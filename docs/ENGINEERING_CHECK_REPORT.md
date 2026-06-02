# AImd Engineering Check Report

## Current Integration State

The active MetaboClip integration uses `MetaBoClipBridge` as the adapter layer and `metaboclip_unified` as the backend. The old implementation has been removed from the clean deliverable package.

The bridge preserves the refined docking manifest interface and writes AImd-compatible outputs under:

```text
data/data_output/metaboclip/results/
```

The stable final ranking file remains:

```text
data/data_output/metaboclip/results/metaboclip_final_ranking.csv
```

## Active Backend Checks

Expected active backend:

```text
metaboclip_unified
```

Expected Python APIs:

```text
metaboclip.core.workflow.run_directory
metaboclip.core.workflow.run_single_pair
```

Expected role-table APIs:

```text
metaboclip_ligand_roles.annotator_core.annotate_ligand
metaboclip_ligand_roles.atom_map.recover_atom_map_by_coordinates
metaboclip_ligand_roles.role_table.annotation_to_role_rows
```

## Manifest Compatibility

`MetaBoClipBridge` consumes:

```text
data/data_output/refined/docking_out/docking_result_manifest.csv
```

Important fields preserved when present:

```csv
job_id,ligand_id,protein_id,cluster_id,batch_id,conformer_id,pocket_id,pocket_rank,receptor_pdbqt_path,ligand_pdbqt_path,config_path,out_pose_path,log_path,center_x,center_y,center_z,size_x,size_y,size_z,status,return_code,message,best_affinity,affinities,n_affinities,grid_size,grid_space,exhaustiveness,random_seed,pose_exists
```

Original docking `status` and `message` are retained as `source_status` and `source_message`; bridge execution uses its own `status` and `message`.

## Output Compatibility

The bridge aggregates unified outputs into:

```text
metaboclip_protein_scores_all.csv
metaboclip_conformation_scores_all.csv
metaboclip_pose_scores_all.csv
metaboclip_candidate_scores_all.csv
metaboclip_passing_candidates_all.csv
metaboclip_geometry_features_all.csv
metaboclip_resolved_ligand_sites_all.csv
metaboclip_resolved_protein_roles_all.csv
metaboclip_run_manifest.csv
metaboclip_report.json
metaboclip_final_ranking.csv
```

The final ranking uses real unified score columns when present:

```csv
protein_score,quality_score,coverage
```

Unavailable legacy score columns are not generated:

```csv
protein_score_norm,max_s_r
```

## Scientific Constraint

Unified role-table generation may use ligand hydrogens for functional-group recognition or atom eligibility checks. Downstream catalytic geometry must use heavy atoms only. The generated atom map is heavy-atom only.

## Recommended Safe Checks

```bash
python validate_aimd_layout.py --root .
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:metaboclip_unified python3 -c "from MetaBoClipBridge.main import run_metaboclip_bridge; from metaboclip.core.workflow import run_directory, run_single_pair; print('imports ok')"
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:metaboclip_unified python3 -c "from MetaBoClipBridge.utils import load_yaml; from MetaBoClipBridge.bridge import validate_unified_bridge_config; cfg = load_yaml('configs/MetaBoClip/metaboclip_bridge.yaml'); print(validate_unified_bridge_config(cfg)['backend'])"
```

Use tiny fixtures or dry-run mode for smoke tests. Do not run full docking or full production workflows as validation checks.

Deterministic MetaboClip bridge smoke test:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:metaboclip_unified python -m pytest tests/test_metaboclip_bridge_smoke.py -q
```
