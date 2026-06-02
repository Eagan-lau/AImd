# AImd Module Interface Specification

This document defines the stable manifest interfaces between AImd modules. Internal tool implementations may change, but these interfaces should remain compatible unless a future migration explicitly changes them.

## Workflow

```text
RGPC
  -> TApocketBridge
  -> DockingHub
  -> ClusterScore
  -> RefinementHub
  -> refined DockingHub
  -> MetaBoClipBridge
```

MetaboClip is a core scientific component of the workflow. The active backend is `metaboclip_unified`, called through `MetaBoClipBridge`; the old implementation is not part of the clean deliverable package.

## protein_manifest.csv

Source: `RGPC`

Downstream modules: `TApocketBridge`, `DockingHub`, `RefinementHub`

Recommended columns:

```csv
protein_id,cluster_id,batch_id,source_pdb,protein_path,is_representative,file_action,status
```

## pocket_manifest.csv

Source: `TApocketBridge`

Downstream module: `DockingHub`

Key columns:

```csv
protein_id,cluster_id,batch_id,pocket_id,pocket_rank,center_x,center_y,center_z,size_x,size_y,size_z,final_score,protein_path,pocket_pdb_path,pocket_json_path,box_yaml_path,status
```

## ligand_manifest.csv

Source: molecule preparation module or manual preparation

Downstream module: `DockingHub`; optional lookup source for `MetaBoClipBridge`

Recommended columns:

```csv
ligand_id,batch_id,ligand_path,smiles,name,status
```

For unified MetaboClip role-table generation, the ligand manifest or ligand source manifest may also provide:

```csv
ligand_source_path,source_ligand_path,source_path,sdf_path,mol_path,mol2_path,prepared_ligand_pdbqt_path,pdbqt_path
```

Existing column names must not be renamed.

## docking_result_manifest.csv

Source: `DockingHub`

Downstream modules: `ClusterScore`, `MetaBoClipBridge`

Key columns:

```csv
job_id,ligand_id,protein_id,cluster_id,batch_id,conformer_id,pocket_id,pocket_rank,receptor_pdbqt_path,ligand_pdbqt_path,config_path,out_pose_path,log_path,center_x,center_y,center_z,size_x,size_y,size_z,status,return_code,message,best_affinity,affinities,n_affinities,grid_size,grid_space,exhaustiveness,random_seed,pose_exists
```

`MetaBoClipBridge` preserves these fields when present and adds bridge/unified backend status fields in downstream outputs.

## MetaBoClipBridge Inputs

Default refined docking input:

```text
data/data_output/refined/docking_out/docking_result_manifest.csv
```

Required for scoring rows:

```csv
ligand_id,protein_id,receptor_pdbqt_path,out_pose_path,log_path,best_affinity
```

Recommended metadata:

```csv
job_id,cluster_id,batch_id,conformer_id,pocket_id,pocket_rank,ligand_pdbqt_path,config_path,status,return_code,message,affinities,n_affinities,grid_size,grid_space,exhaustiveness,random_seed,pose_exists
```

## Unified MetaboClip Inputs

Configured in `configs/MetaBoClip/metaboclip_bridge.yaml`:

```yaml
backend: unified
paths:
  metaboclip_project_dir: metaboclip_unified
  metaboclip_profile: metaboclip_unified/metaboclip/config/profiles/default_profile.yaml
  unified_output_dir: data/data_output/metaboclip/unified_runs
  ligand_manifest: data/data_input/ligand/ligand_manifest.csv
  ligand_source_manifest: data/data_input/ligand/ligand_source_manifest.csv
  role_table_dir: data/data_output/metaboclip/ligand_roles/role_tables
  annotation_dir: data/data_output/metaboclip/ligand_roles/annotations
  atom_map_dir: data/data_output/metaboclip/ligand_roles/atom_maps
mechanisms:
  cyp450: metaboclip_unified/metaboclip/config/families/cyp450/mechanism.yaml
```

Role generation requires an original ligand source file and a prepared ligand PDBQT file. Role generation writes role tables, annotation JSON, and atom-map JSON. The atom map is heavy-atom only.

## MetaBoClipBridge Outputs

Stable final output:

```text
data/data_output/metaboclip/results/metaboclip_final_ranking.csv
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

Unified backend artifacts are written under:

```text
data/data_output/metaboclip/unified_runs/
```

The final ranking uses real unified score columns when available:

```csv
protein_score,quality_score,coverage
```

Legacy scientific columns that cannot be safely derived are not generated:

```csv
protein_score_norm,max_s_r
```

## Replacement Principles

1. Preserve manifest readability.
2. Preserve existing input column names.
3. Use project-root-relative or explicitly configured paths.
4. Keep failed rows visible with clear `status` and `message` fields.
5. Do not restore the legacy MetaboClip backend as the active workflow.
6. Do not use ligand hydrogen coordinates for downstream catalytic geometry.
