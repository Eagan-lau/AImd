# AImd Integration Instructions for Codex

## Project context

This repository is the AImd engineering package.

The new MetaboClip core logic is located at:

third_party/metaboclip_unified

The old MetaboClip implementation is preserved only as legacy code at:

third_party/MetaBoClip_legacy

The task is to integrate the new MetaboClip logic into AImd through MetaBoClipBridge.

## Source-of-truth priority

When project documents or code comments conflict, use this AGENTS.md file as the highest-priority project instruction.

The old README, old MetaBoClip documentation, and legacy comments may describe outdated behavior. Do not restore old behavior unless explicitly instructed.

## Integration principle

Do not blindly overwrite source files.

Use MetaBoClipBridge as the adapter layer between AImd and the new MetaboClip core.

Keep third_party/metaboclip_unified as a standalone core module.

Do not modify third_party/metaboclip_unified unless a minimal compatibility wrapper or import fix is strictly required.

Do not restore or reintroduce the old MetaboClip logic as the active workflow.

Do not delete third_party/MetaBoClip_legacy unless explicitly instructed.

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

Ligand hydrogen atoms may be used only during functional-group detection to identify atom eligibility and avoid misclassification, such as distinguishing hydroxyl oxygen from carbonyl oxygen or identifying carbon atoms bearing hydrogen for CH sites.

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

The integrated AImd workflow must keep stable downstream output conventions unless explicitly instructed.

## Path and configuration rules

Do not introduce hard-coded local paths into reusable source code.

Do not hard-code paths such as:

/media/yugengliu
DATA3
AImd_integrated_clean_manual_v4_third_party

Use project-root-relative paths, configuration files, or command-line arguments.

If an absolute path is needed for a local test, keep it only in temporary notes or user-specific examples, not in reusable source code.

## Safety constraints

Do not change the core scientific scoring logic unless explicitly instructed.

Do not delete, rewrite, or batch-format large data directories.

Do not scan or modify large docking output directories unless explicitly instructed.

Do not run full-scale docking, full-scale screening, or full AImd production workflows unless explicitly instructed.

Use minimal smoke tests only.

All generated source code, comments, README text, configuration examples, YAML files, JSON files, CSV headers, plot labels, legends, titles, and documentation must be in English only.

## Required workflow

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
7. Troubleshooting for common path, manifest, and dependency errors.
