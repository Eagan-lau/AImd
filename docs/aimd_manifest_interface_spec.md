# AImd Manifest Interface Specification

This file records the stable AImd manifest interfaces used by the integrated workflow. The new unified MetaboClip backend is integrated internally through `MetaBoClipBridge`; upstream and downstream manifest paths remain stable.

## Overall Data Flow

```text
MolLink
  inputs:
    data/data_input/ligand/taxane_molecules.csv
  outputs:
    data/data_output/ligand_transformation/ligand_source_manifest.csv

RGPC
  outputs:
    data/data_output/protein_batches/protein_manifest.csv
    data/data_output/cluster/RGPC/clusters.tsv
    data/data_output/cluster/RGPC/representatives.tsv

TApocketBridge
  inputs:
    data/data_output/protein_batches/protein_manifest.csv
  outputs:
    data/data_output/pocket/pocket_manifest.csv
    data/data_output/pocket/tapocket_run_manifest.csv

DockingHub broad
  inputs:
    data/data_output/protein_batches/protein_manifest.csv
    data/data_output/pocket/pocket_manifest.csv
    data/data_output/ligand_preparation/ligand_manifest.csv
  outputs:
    data/data_output/ensemble/conformer_manifest.csv
    data/data_output/docking_out/docking_result_manifest.csv

ClusterScore
  inputs:
    data/data_output/docking_out/docking_result_manifest.csv
  outputs:
    data/data_output/scoring/ClusterScore/clusterscore_results.xlsx
    data/data_output/scoring/ClusterScore/top10_clusters.csv

RefinementHub
  inputs:
    data/data_output/scoring/ClusterScore/top10_clusters.csv
    data/data_output/protein_batches/protein_manifest.csv
  outputs:
    data/data_output/refinement/selected_protein_manifest.csv
    data/data_output/refinement/refined_docking.generated.yaml
    data/data_output/refined/docking_out/docking_result_manifest.csv

MetaBoClipBridge
  inputs:
    data/data_output/refined/docking_out/docking_result_manifest.csv
  outputs:
    data/data_output/metaboclip/results/metaboclip_final_ranking.csv
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

## ligand_source_manifest.csv

Source: `MolLink`

Downstream modules: ligand preparation; optional lookup source for `MetaBoClipBridge`

Recommended columns:

```csv
ligand_id,source_id,molecule_id,molecule_name,smiles,ligand_source_type,source_table,row_index,transform_status
```

## ligand_manifest.csv

Recommended columns:

```csv
ligand_id,batch_id,ligand_path,pdbqt_path,source_csv,source_row,smiles,name,sdf_path,pdb_path,embedding_status,optimization_status,prepare_ligand_return_code,status,message
```

Optional columns used by `MetaBoClipBridge` for role-table generation:

```csv
ligand_source_path,source_ligand_path,source_path,mol_path,mol2_path,prepared_ligand_pdbqt_path
```

## conformer_manifest.csv

Source: `DockingHub.ensemble`

Downstream modules: `DockingHub.receptor`, `DockingHub.tasks`

Recommended columns:

```csv
protein_id,cluster_id,batch_id,conformer_id,structure_path,original_structure_path,unaligned_structure_path,coordinate_frame,alignment_status,alignment_message,source,status,message
```

Optional AlphaFlow adapter columns:

```csv
alphaflow_name,alphaflow_project_dir,alphaflow_msa_dir,alphaflow_templates_dir
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
  unified_output_dir: data/data_output/metaboclip/unified_runs
  ligand_manifest: data/data_input/ligand/ligand_manifest.csv
  ligand_source_manifest: data/data_input/ligand/ligand_source_manifest.csv
  role_table_dir: data/data_output/metaboclip/ligand_roles/role_tables
  annotation_dir: data/data_output/metaboclip/ligand_roles/annotations
  atom_map_dir: data/data_output/metaboclip/ligand_roles/atom_maps
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
data/data_output/metaboclip/results/metaboclip_final_ranking.csv
```

Aggregate outputs:

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
