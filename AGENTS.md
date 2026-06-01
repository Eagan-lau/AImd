# AImd Integration Instructions for Codex

## Project context

This repository is the AImd engineering package.

AImd treats MetaboClip as a core scientific component of the workflow.

The active unified MetaboClip core logic is located at:

third_party/metaboclip_unified

The old MetaboClip implementation has been removed from the clean deliverable package.

The task is to keep the unified MetaboClip logic integrated into AImd through MetaBoClipBridge.

## Source-of-truth priority

When project documents, code comments, legacy documentation, or README sections conflict, use this AGENTS.md file as the highest-priority project instruction.

Old README sections, old MetaBoClip documentation, and legacy comments may describe outdated behavior. Do not restore old behavior unless explicitly instructed.

## AImd interface specification

Before modifying MetaBoClipBridge, inspect the current AImd manifest interfaces and preserve upstream and downstream compatibility.

If the file exists, read:

docs/aimd_manifest_interface_spec.md

This file defines the current AImd module interfaces and manifest schemas. The integration of the new MetaboClip backend must preserve these interfaces unless explicitly instructed.

Internal backend replacement is allowed. Upstream and downstream manifest compatibility must be preserved.

## Current AImd module flow

The expected high-level AImd data flow is:

1. RGPC generates protein and cluster manifests.
2. TApocketBridge consumes the protein manifest and generates the pocket manifest.
3. DockingHub consumes protein, pocket, and ligand manifests and generates docking results.
4. ClusterScore consumes the broad docking results and selects candidate clusters.
5. RefinementHub generates selected protein manifests and refined docking outputs.
6. MetaBoClipBridge consumes refined docking results and generates MetaboClip-compatible ranking outputs.

The current MetaBoClipBridge input is expected to be:

data/refined/docking_out/docking_result_manifest.csv

The current MetaBoClipBridge output must remain compatible with:

data/metaboclip/results/metaboclip_final_ranking.csv

## Important AImd manifest fields

The refined docking manifest may include the following important fields:

job_id
ligand_id
protein_id
cluster_id
batch_id
conformer_id
pocket_id
pocket_rank
receptor_pdbqt_path
ligand_pdbqt_path
config_path
out_pose_path
log_path
center_x
center_y
center_z
size_x
size_y
size_z
status
return_code
message
best_affinity
affinities
n_affinities
grid_size
grid_space
exhaustiveness
random_seed
pose_exists

MetaBoClipBridge must preserve as much AImd metadata as possible when adapting the refined docking manifest to the new MetaboClip backend.

## Integration principle

Do not blindly overwrite source files.

Use MetaBoClipBridge as the adapter layer between AImd and the unified MetaboClip core.

Keep third_party/metaboclip_unified as a standalone core module.

Do not modify third_party/metaboclip_unified unless a minimal compatibility wrapper, import fix, or packaging fix is strictly required.

Do not restore or reintroduce the old MetaBoClip logic.

## New MetaboClip logic

The new MetaboClip logic includes:

1. Ligand functional-group detection.
2. Ligand role-table generation.
3. Ligand annotation JSON output.
4. Ligand atom-map JSON output.
5. Family-specific YAML mechanism parsing.
6. Protein-side catalytic residue extraction.
7. Ligand-side reactive atom assignment.
8. Pose-level catalytic geometry filtering.
9. Final scoring and ranking.

## Critical scientific constraint

Catalytic geometry calculations must use heavy atoms only.

Ligand hydrogen atoms may be used only during functional-group detection to identify atom eligibility and avoid misclassification, such as distinguishing hydroxyl oxygen from carbonyl oxygen or identifying carbon atoms bearing hydrogen for CH sites.

After the functional group or reactive site is confirmed, downstream catalytic geometry must be computed using heavy atoms, not ligand hydrogen atoms.

Do not introduce geometry calculations that depend on ligand hydrogen coordinates.

## Integration requirements

MetaBoClipBridge must preserve and correctly pass:

1. protein_id
2. ligand_id
3. pose_id
4. pocket_id when available
5. docking score or best_affinity
6. protein pose path or receptor path
7. ligand pose path
8. ligand source path when available
9. prepared ligand PDBQT path when available
10. role-table path
11. annotation JSON path
12. atom-map JSON path
13. catalytic geometry outputs
14. final ranking outputs

The integrated AImd workflow must keep stable downstream output conventions unless explicitly instructed.

## New backend input requirements

The new MetaboClip backend may require information that is not fully present in the old AImd refined docking manifest.

If required inputs are missing, do not hard-code local paths. Add project-root-relative or config-driven fields.

Potential required inputs include:

1. ligand_source_dir
2. prepared_ligand_pdbqt_dir
3. role_table_dir
4. annotation_dir
5. atom_map_dir
6. mechanism_config_root
7. metabo_clip_profile
8. staged_protein_dir
9. staged_docking_dir
10. unified_output_dir

When adding new config fields, update the bridge config and documentation together.

## Output compatibility rules

The new MetaboClip backend may produce output files such as:

resolved_ligand_sites.csv
resolved_protein_roles.csv
geometry_features.csv
candidate_scores.csv
pose_scores.csv
merged_conformation_scores.csv
protein_scores.csv

MetaBoClipBridge must translate the new backend outputs into stable AImd-compatible outputs under:

data/metaboclip/results/

The final AImd-facing ranking file should remain:

data/metaboclip/results/metaboclip_final_ranking.csv

Do not invent legacy score columns if they cannot be mapped safely.

If legacy columns such as protein_score_norm or max_s_r are required by downstream modules, explicitly document how they are derived from new backend outputs. If no valid derivation exists, preserve the raw new backend score columns and mark the legacy columns as deprecated or unavailable.

## Path and configuration rules

Do not introduce hard-coded local paths into reusable source code.

Do not hard-code paths such as:

/media/yugengliu
DATA3
AImd_integrated_clean_manual_v4_third_party

Use one of the following:

1. project-root-relative paths
2. configuration files
3. command-line arguments
4. environment variables only when documented

If an absolute path is needed for a temporary local test, keep it only in a temporary note or one-off command, not in reusable source code, configs, examples, tests, or documentation.

## Safety constraints

Do not change the core scientific scoring logic unless explicitly instructed.

Do not delete, rewrite, or batch-format large data directories.

Do not scan or modify large docking output directories unless explicitly instructed.

Do not run full-scale docking, full-scale screening, or full AImd production workflows unless explicitly instructed.

Use minimal smoke tests only.

Do not run commands that may modify user data unless the user explicitly approves them.

Do not recreate the removed old MetaBoClip implementation unless explicitly instructed.

All generated source code, comments, README text, configuration examples, YAML files, JSON files, CSV headers, plot labels, legends, titles, and documentation must be in English only.

## Required workflow for integration work

For integration work:

1. Inspect the current AImd MetaBoClip-related code first.
2. Inspect third_party/metaboclip_unified first.
3. Do not modify files during the first analysis step.
4. Report the current AImd MetaBoClip interface.
5. Report the new MetaboClip CLI or Python API.
6. Compare current AImd expectations with the new MetaboClip input and output schema.
7. Report old MetaBoClip assumptions that remain in AImd.
8. Propose the smallest integration plan.
9. List files that need modification before editing.
10. Wait for explicit confirmation before editing files.
11. After editing, run a minimal smoke test.
12. Summarize the exact commands used.
13. Summarize git diff.

## Required workflow for adapter design

Before editing files, produce a concrete adapter design that includes:

1. Current AImd call chain into MetaBoClipBridge.
2. Current required input manifest columns for MetaBoClipBridge.
3. Required new MetaboClip inputs.
4. Which required new inputs are already available from AImd manifests.
5. Which required new inputs need new config fields.
6. Mapping from AImd refined docking manifest to new MetaboClip inputs.
7. Mapping from new MetaboClip outputs to AImd-compatible outputs.
8. Which legacy result columns can be safely preserved.
9. Which legacy result columns must not be invented.
10. Exact files that need modification.
11. Minimal smoke test plan.

Do not edit files during this adapter-design step.

## Bug-fixing workflow

For bug fixing:

1. Reproduce the bug using the exact command provided by the user.
2. Identify the smallest failing module.
3. Report the root cause before editing files.
4. Propose the smallest safe fix.
5. Modify only the necessary files.
6. Preserve the new MetaboClip logic.
7. Preserve the output schema unless explicitly instructed.
8. Run the smallest possible smoke test.
9. Summarize git diff.

## Documentation requirements

Documentation must describe the current new MetaboClip integration, not the old legacy implementation.

README and example commands must explain:

1. How AImd calls MetaBoClipBridge.
2. How MetaBoClipBridge calls third_party/metaboclip_unified.
3. Required input files.
4. Required configuration files.
5. Expected output files.
6. Minimal smoke-test commands.
7. Troubleshooting for common path, manifest, dependency, and missing role-table errors.

Documentation must not contain Chinese characters.

## Dependency rules

Do not assume PyMOL is required for core MetaboClip scoring.

PyMOL-related commands should be treated as export or visualization functionality only unless the new backend explicitly requires them for a specific operation.

If the new backend requires RDKit, SciPy, NumPy, pandas, PyYAML, or other dependencies, update the environment documentation or requirements file only after verifying the import requirements.

Do not install dependencies automatically unless explicitly instructed.

## Git workflow

Before editing files, check:

git status --short

After editing files, summarize:

1. files changed
2. purpose of each change
3. commands run
4. smoke-test result
5. git diff summary

Do not commit changes unless explicitly instructed.

## GitHub delivery workflow

For this project, the user has explicitly requested GitHub synchronization after completed engineering updates.

After a controlled update is implemented and the relevant safe checks pass:

1. Commit the completed change set.
2. Push it to `origin/main`.
3. Report the commit hash and push result.

Do not push half-finished work, failed validation states, generated cache files, large data/output directories, model weights, or temporary files. If the user explicitly says not to commit or not to push for a task, follow that newer instruction.

## Minimal smoke-test policy

A minimal smoke test should use tiny example inputs only.

It must not run full docking, full screening, or large-scale AImd workflows.

The smoke test should verify only the integration boundary:

1. AImd can call MetaBoClipBridge.
2. MetaBoClipBridge can locate the new backend.
3. Required role-table, annotation, and atom-map paths are handled.
4. The new backend can produce or expose expected output files.
5. MetaBoClipBridge can produce an AImd-compatible final ranking file or a clearly documented dry-run report.

If no suitable smoke test exists, propose one before creating files.
