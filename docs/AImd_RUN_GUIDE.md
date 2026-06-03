# AImd Run Guide

This guide describes what each AImd module does, which files it consumes, which files it writes, and how to run the workflow in a controlled order. AImd is intended to work as a plug-and-play declarative package: provide protein and ligand inputs, select tool plugins and module switches in YAML, and run the full workflow command when the environment is ready.

## Workflow Summary

```text
MolLink
  -> RGPC
  -> TApocketBridge
  -> DockingHub broad docking
  -> ClusterScore
  -> RefinementHub and refined DockingHub
  -> MetaBoClipBridge
  -> MetaBoClipHub
```

The workflow separates user inputs from generated outputs:

```text
data/data_input/
data/data_output/
```

Place only starting files and manual control tables under `data/data_input`. Keep all generated module outputs under `data/data_output`.

## Environment

Use the AImd environment for normal Python modules:

```bash
conda deactivate
conda activate aimd
```

Install or verify Python dependencies:

```bash
pip install -r requirements.txt
python validate_aimd_layout.py --root .
python third_party/check_tools.py --config third_party/tools.yaml --root .
```

Some external tools may use their own environments. AlphaFlow is normally run from an AlphaFlow-specific environment. AImd records the required commands and paths through configuration instead of vendoring the third-party model runtime.

## Input Layout

Protein starting structures:

```text
data/data_input/protein/file_1/*.pdb
```

Ligand starting table:

```text
data/data_input/ligand/taxane_molecules.csv
```

Cofactor templates:

```text
data/data_input/cofactor/file_1/*.pdb
```

Optional manual workflow tables:

```text
data/data_input/workflow/
```

## 1. MolLink

Purpose:

MolLink is the ligand transformation module. It reads the starting molecule table, records ligand source metadata, and builds a transformation network when reaction templates are available.

Main input:

```text
data/data_input/ligand/taxane_molecules.csv
```

Main command:

```bash
python run_mollink.py --config configs/MolLink/mollink.yaml
```

Main outputs:

```text
data/data_output/ligand_transformation/ligand_source_manifest.csv
data/data_output/ligand_transformation/transformation_edges.csv
data/data_output/ligand_transformation/transformation_nodes.csv
```

If no reaction template is configured, MolLink still completes in `csv_only_no_template` mode. It records molecules from the CSV and does not invent transformation edges.

## 2. RGPC

Purpose:

RGPC batches protein structures, computes structural similarity, and clusters proteins. Foldseek and HipMCL are the main third-party tools for large-scale protein similarity and graph clustering.

Main input:

```text
data/data_input/protein/file_1/*.pdb
```

Main command:

```bash
python run_rgpc.py --config configs/RGPC/rgpc.yaml
```

Main outputs:

```text
data/data_output/protein_batches/protein_manifest.csv
data/data_output/cluster/
```

The protein manifest is consumed by TApocketBridge and DockingHub.

## 3. TApocketBridge

Purpose:

TApocketBridge predicts or maps catalytic pockets for each protein. It writes pocket boxes used by DockingHub.

Main command:

```bash
python run_tapocket_batch.py --config configs/TApocket/tapocket_batch.yaml
```

Main output:

```text
data/data_output/pocket/pocket_manifest.csv
```

The pocket manifest must contain pocket centers and box sizes:

```csv
protein_id,pocket_id,pocket_rank,center_x,center_y,center_z,size_x,size_y,size_z,status
```

## 4. DockingHub Ligand Preparation

Purpose:

DockingHub converts ligand SMILES from the molecule CSV into 3D structures, minimizes them with RDKit, and converts them to PDBQT for docking.

Main config:

```text
configs/Docking/docking.yaml
```

Main inputs:

```text
data/data_input/ligand/taxane_molecules.csv
```

Main outputs:

```text
data/data_output/ligand_preparation/sdf/
data/data_output/ligand_preparation/pdb/
data/data_output/ligand_preparation/pdbqt/
data/data_output/ligand_preparation/ligand_manifest.csv
```

The ligand preparation manifest is also used by MetaBoClipBridge to generate ligand role tables.

## 5. DockingHub Broad Docking

Purpose:

DockingHub prepares receptors, generates docking task files, runs AutoDock Vina, and writes the docking result manifest.

Main command:

```bash
python run_docking.py --config configs/Docking/docking.yaml
```

Main outputs:

```text
data/data_output/receptor/receptor_manifest.csv
data/data_output/docking_tasks/docking_task_manifest.csv
data/data_output/docking_out/docking_result_manifest.csv
```

The docking result manifest preserves AImd metadata such as:

```csv
job_id,ligand_id,protein_id,cluster_id,batch_id,conformer_id,pocket_id,pocket_rank,receptor_pdbqt_path,ligand_pdbqt_path,config_path,out_pose_path,log_path,best_affinity,status
```

## 6. Cofactor-Aware Receptor Handling

Purpose:

DockingHub can transfer cofactors from template structures into target receptors. The transfer is gated by local cofactor-pocket geometry after global structural alignment.

Main config block:

```yaml
cofactor:
  enabled: true
  alignment:
    pocket_validation:
      enabled: true
      pocket_radius: 5.0
      rmsd_metric: ca
      local_rmsd_cutoff: 10.0
      min_mapped_residues: 5
      min_pocket_coverage: 0.70
      require_pass_before_transfer: true
```

Decision logic:

1. Select a cofactor template by Foldseek or the first available template.
2. Align the template protein to the target protein with PyMOL.
3. Identify template pocket residues with heavy atoms within `pocket_radius` of cofactor heavy atoms.
4. Map template pocket residues to target residues through the global alignment.
5. Compute local CA RMSD using mapped pocket residues.
6. Transfer the cofactor only if local RMSD, mapped residue count, and pocket coverage pass the configured gates.

Main output:

```text
data/data_output/refined/cofactor_mapped/cofactor_manifest.csv
```

Important cofactor validation fields:

```csv
cofactor_pocket_radius,cofactor_pocket_residue_count,mapped_pocket_residue_count,cofactor_pocket_coverage,cofactor_site_ca_rmsd,cofactor_local_rmsd_cutoff,cofactor_transfer_pass,cofactor_validation_status,cofactor_validation_message,transfer_mode
```

If validation fails and `continue_without_cofactor_on_error` is true, DockingHub continues with the original receptor and marks the row as `success_no_cofactor`.

## 7. ClusterScore

Purpose:

ClusterScore summarizes broad docking results and selects top protein clusters or families for refined docking.

Main command:

```bash
python run_clusterscore.py --config configs/Scoring/cluster_score.yaml
```

Main outputs:

```text
data/data_output/scoring/ClusterScore/top10_clusters.csv
data/data_output/scoring/ClusterScore/cluster_binding_statistics.csv
```

## 8. RefinementHub and Refined Docking

Purpose:

RefinementHub selects proteins or clusters for refined docking and generates a refined DockingHub configuration. Refined DockingHub can use AlphaFlow conformers, cofactor validation, receptor preparation, and Vina docking.

Main command:

```bash
python run_refinement.py --config configs/Refinement/refine_from_clusterscore.yaml
```

Main refined outputs:

```text
data/data_output/refinement/selected_protein_manifest.csv
data/data_output/refined/cofactor_mapped/cofactor_manifest.csv
data/data_output/refined/receptor/receptor_manifest.csv
data/data_output/refined/docking_tasks/docking_task_manifest.csv
data/data_output/refined/docking_out/docking_result_manifest.csv
```

If AlphaFlow conformers already exist, a controlled engineering run can reuse `data/data_output/ensemble/conformer_manifest.csv` and write refined docking outputs under `data/data_output/refined/`.

## 9. MetaBoClipBridge

Purpose:

MetaBoClipBridge stages refined docking outputs, generates or reuses ligand role assets, calls the MetaBoClipHub scoring layer, and writes AImd-compatible catalytic scoring outputs.

Main command:

```bash
PYTHONPATH=.:metaboclip_unified \
python run_metaboclip_bridge.py --config configs/MetaBoClip/metaboclip_bridge.yaml
```

Main input:

```text
data/data_output/refined/docking_out/docking_result_manifest.csv
```

Role assets:

```text
data/data_output/metaboclip/ligand_roles/role_tables/
data/data_output/metaboclip/ligand_roles/annotations/
data/data_output/metaboclip/ligand_roles/atom_maps/
```

Unified backend outputs:

```text
data/data_output/metaboclip/unified_runs/
```

AImd-facing outputs:

```text
data/data_output/metaboclip/results/metaboclip_run_manifest.csv
data/data_output/metaboclip/results/metaboclip_final_ranking.csv
data/data_output/metaboclip/results/metaboclip_report.json
```

Aggregate score tables:

```text
data/data_output/metaboclip/results/metaboclip_protein_scores_all.csv
data/data_output/metaboclip/results/metaboclip_conformation_scores_all.csv
data/data_output/metaboclip/results/metaboclip_pose_scores_all.csv
data/data_output/metaboclip/results/metaboclip_candidate_scores_all.csv
data/data_output/metaboclip/results/metaboclip_geometry_features_all.csv
data/data_output/metaboclip/results/metaboclip_resolved_ligand_sites_all.csv
data/data_output/metaboclip/results/metaboclip_resolved_protein_roles_all.csv
```

MetaBoClipBridge does not fabricate unavailable compatibility score columns such as:

```csv
protein_score_norm,max_s_r
```

It preserves AImd metadata and uses real MetaBoClipHub score fields such as:

```csv
protein_score,quality_score,coverage
```

## 10. MetaBoClipHub

Purpose:

MetaBoClipHub is the AImd-facing result hub for MetaboClip outputs. It organizes role assets, runtime run folders, aggregate score tables, reports, and final rankings under:

```text
data/data_output/metaboclip/
```

The MetaBoClipHub artifacts are produced by the `run_metaboclip_bridge.py` command. No separate scoring command is required for this output layer.

## 11. Validation Commands

Compile Python modules with the AImd Python environment:

```bash
/home/yugengliu/anaconda3/envs/aimd/bin/python -m compileall DockingHub MetaBoClipBridge tests
```

Run repository-owned smoke tests:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:metaboclip_unified \
python3 -m pytest -q tests
```

Validate MetaBoClipBridge outputs:

```bash
python MetaBoClipBridge/validate_metaboclip_outputs.py --out-dir data/data_output/metaboclip/results
```

Full-repository `pytest` may collect third-party model tests under vendored tool directories. Use the repository-owned `tests/` directory for AImd smoke testing unless the third-party test environments are explicitly configured.

## 12. Practical Run Order

For a controlled engineering run from starting inputs:

```bash
conda deactivate
conda activate aimd
python validate_aimd_layout.py --root .
python third_party/check_tools.py --config third_party/tools.yaml --root .
python run_mollink.py --config configs/MolLink/mollink.yaml
python run_rgpc.py --config configs/RGPC/rgpc.yaml
python run_tapocket_batch.py --config configs/TApocket/tapocket_batch.yaml
python run_docking.py --config configs/Docking/docking.yaml
python run_clusterscore.py --config configs/Scoring/cluster_score.yaml
python run_refinement.py --config configs/Refinement/refine_from_clusterscore.yaml
PYTHONPATH=.:metaboclip_unified python run_metaboclip_bridge.py --config configs/MetaBoClip/metaboclip_bridge.yaml
```

For an engineering rerun that already has AlphaFlow conformers and ligand preparation outputs, reuse:

```text
data/data_output/ensemble/conformer_manifest.csv
data/data_output/ligand_preparation/ligand_manifest.csv
data/data_output/pocket/pocket_manifest.csv
data/data_input/cofactor/file_1/*.pdb
```

Then run refined cofactor-aware ensemble docking and MetaBoClipBridge against `data/data_output/refined/`.
