# AImd Integration Instructions for Codex

## Project context

This repository is the AImd engineering package.

The new MetaboClip core logic is located at:

third_party/metaboclip_unified

The old MetaboClip implementation is preserved only as legacy code at:

third_party/MetaBoClip_legacy

The task is to integrate the new MetaboClip logic into AImd through MetaBoClipBridge.

## Integration principle

Do not blindly overwrite source files.

Use MetaBoClipBridge as the adapter layer between AImd and the new MetaboClip core.

Do not restore or reintroduce the old MetaboClip logic as the active workflow.

## New MetaboClip logic

The new MetaboClip logic includes:

1. Ligand functional-group detection.
2. Ligand role-table generation.
3. Ligand annotation JSON output.
4. Ligand atom-map JSON output.
5. Family-specific YAML rule parsing.
6. Protein-side catalytic residue extraction.
7. Ligand-side reactive atom assignment.
8. Pose-level catalytic geometry filtering.
9. Final scoring and ranking.

## Critical scientific constraint

Catalytic geometry calculations must use heavy atoms only.

Ligand hydrogen atoms may be used only during functional-group detection to identify atom eligibility and avoid misclassification.

After the functional group or reactive site is confirmed, downstream catalytic geometry must be computed using heavy atoms, not ligand hydrogen atoms.

## Integration requirements

MetaBoClipBridge must preserve and correctly pass:

1. protein_id
2. ligand_id
3. pose_id
4. pocket_id when available
5. docking score
6. protein pose path
7. ligand pose path
8. role-table path
9. annotation JSON path
10. atom-map JSON path
11. catalytic geometry outputs
12. final ranking outputs

## Safety constraints

Do not change the core scientific scoring logic unless explicitly instructed.

Do not delete, rewrite, or batch-format large data directories.

Do not introduce hard-coded local paths into reusable source code.

All generated source code, comments, README text, configuration examples, YAML files, JSON files, CSV headers, plot labels, legends, titles, and documentation must be in English only.

## Required workflow

For integration work:

1. Inspect the current AImd MetaBoClip-related code first.
2. Inspect third_party/metaboclip_unified first.
3. Do not modify files during the first analysis step.
4. Report the current AImd MetaBoClip interface.
5. Report the new MetaboClip CLI or Python API.
6. Propose the smallest integration plan.
7. List files that need modification before editing.
8. After editing, run a minimal smoke test.
9. Summarize git diff.
