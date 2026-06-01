# AImd Manifest Interface Specification

This file records the stable AImd manifest interfaces used by the integrated workflow. The new unified MetaboClip backend is integrated internally through `MetaBoClipBridge`; upstream and downstream manifest paths remain stable.

## Overall Data Flow

```text
RGPC
  outputs:
    data/protein/protein_manifest.csv
    data/cluster/RGPC/clusters.tsv
    data/cluster/RGPC/representatives.tsv

TApocketBridge
  inputs:
    data/protein/protein_manifest.csv
  outputs:
    data/pocket/pocket_manifest.csv
    data/pocket/tapocket_run_manifest.csv

DockingHub broad
  inputs:
    data/protein/protein_manifest.csv
    data/pocket/pocket_manifest.csv
    data/ligand/ligand_manifest.csv
  outputs:
    data/docking_out/docking_result_manifest.csv

ClusterScore
  inputs:
    data/docking_out/docking_result_manifest.csv
  outputs:
    data/scoring/ClusterScore/clusterscore_results.xlsx
    data/scoring/ClusterScore/top10_clusters.csv

RefinementHub
  inputs:
    data/scoring/ClusterScore/top10_clusters.csv
    data/protein/protein_manifest.csv
  outputs:
    data/refinement/selected_protein_manifest.csv
    data/refinement/refined_docking.generated.yaml
    data/refined/docking_out/docking_result_manifest.csv

MetaBoClipBridge
  inputs:
    data/refined/docking_out/docking_result_manifest.csv
  outputs:
    data/metaboclip/results/metaboclip_final_ranking.csv
```

## protein_manifest.csv

Recommended columns:

```csv
protein_id,cluster_id,batch_id,source_pdb,protein_path,is_representative,file_action,status
```

## pocket_manifest.csv

Key columns:

```csv
protein_id,cluster_id,batch_id,pocket_id,pocket_rank,center_x,center_y,center_z,size_x,size_y,size_z,final_score,protein_path,pocket_pdb_path,pocket_json_path,box_yaml_path,status
```

## ligand_manifest.csv

Recommended columns:

```csv
ligand_id,batch_id,ligand_path,smiles,name,status
```

Optional columns used by `MetaBoClipBridge` for role-table generation:

```csv
ligand_source_path,source_ligand_path,source_path,sdf_path,mol_path,mol2_path,prepared_ligand_pdbqt_path,pdbqt_path
```

## docking_result_manifest.csv

Source: `DockingHub`

Downstream modules: `ClusterScore`, `MetaBoClipBridge`

Key columns:

```csv
job_id,ligand_id,protein_id,cluster_id,batch_id,conformer_id,pocket_id,pocket_rank,receptor_pdbqt_path,ligand_pdbqt_path,config_path,out_pose_path,log_path,center_x,center_y,center_z,size_x,size_y,size_z,status,return_code,message,best_affinity,affinities,n_affinities,grid_size,grid_space,exhaustiveness,random_seed,pose_exists
```

`MetaBoClipBridge` preserves these fields when present. Bridge-level status is reported in `status` and `message`; original docking status and message are also retained as `source_status` and `source_message`.

## MetaBoClipBridge Configuration

Active backend:

```yaml
backend: unified
```

Required unified backend paths:

```yaml
paths:
  metaboclip_project_dir: metaboclip_unified
  metaboclip_profile: metaboclip_unified/metaboclip/config/profiles/default_profile.yaml
  unified_output_dir: data/metaboclip/unified_runs
  ligand_manifest: data/ligand/ligand_manifest.csv
  ligand_source_manifest: data/ligand/ligand_source_manifest.csv
  role_table_dir: data/metaboclip/ligand_roles/role_tables
  annotation_dir: data/metaboclip/ligand_roles/annotations
  atom_map_dir: data/metaboclip/ligand_roles/atom_maps
```

Required mechanism mapping:

```yaml
mechanisms:
  cyp450: metaboclip_unified/metaboclip/config/families/cyp450/mechanism.yaml
  fe2og: metaboclip_unified/metaboclip/config/families/fe2og/mechanism.yaml
  ugt: metaboclip_unified/metaboclip/config/families/ugt/mechanism.yaml
  act: metaboclip_unified/metaboclip/config/families/act/mechanism.yaml
  ach: metaboclip_unified/metaboclip/config/families/ach/mechanism.yaml
```

Role-table modes:

```yaml
role_tables:
  mode: existing
```

Allowed values are `existing`, `generate`, and `auto`.

## MetaBoClipBridge Outputs

Stable final result:

```text
data/metaboclip/results/metaboclip_final_ranking.csv
```

Aggregate outputs:

```text
data/metaboclip/results/metaboclip_protein_scores_all.csv
data/metaboclip/results/metaboclip_conformation_scores_all.csv
data/metaboclip/results/metaboclip_pose_scores_all.csv
data/metaboclip/results/metaboclip_candidate_scores_all.csv
data/metaboclip/results/metaboclip_passing_candidates_all.csv
data/metaboclip/results/metaboclip_geometry_features_all.csv
data/metaboclip/results/metaboclip_resolved_ligand_sites_all.csv
data/metaboclip/results/metaboclip_resolved_protein_roles_all.csv
data/metaboclip/results/metaboclip_run_manifest.csv
data/metaboclip/results/metaboclip_report.json
```

Final ranking is sorted by real unified score fields when available:

```csv
protein_score,quality_score,coverage,family,ligand_id,protein_id
```

Legacy score columns that are not safely derivable from the unified backend are not fabricated:

```csv
protein_score_norm,max_s_r
```

## Scientific Constraint

Catalytic geometry must use heavy atoms only. Ligand hydrogens may be used by functional-group detection or eligibility checks, but downstream catalytic geometry must not use ligand hydrogen coordinates.
