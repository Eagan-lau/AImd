# AImd

AImd is a modular engineering package for protein clustering, pocket prediction, docking, cluster scoring, refined docking, and catalytic scoring.

The current workflow is:

```text
RGPC -> TApocketBridge -> DockingHub -> ClusterScore -> RefinementHub -> refined DockingHub -> MetaBoClipBridge
```

MetaboClip is a core scientific component of the AImd workflow. `MetaBoClipBridge` calls the unified MetaboClip core under `metaboclip_unified`; the old implementation is not part of the clean deliverable package.

## Quick Start

```bash
cd AImd
pip install -r requirements.txt
python validate_aimd_layout.py --root .
```

Run the deterministic MetaboClip bridge smoke test:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:metaboclip_unified \
  python -m pytest tests/test_metaboclip_bridge_smoke.py -q
```

Run the full workflow only after all required inputs and external tools are available:

```bash
python run_full_iterative_metaboclip.py \
  --config configs/workflows/full_iterative_metaboclip.yaml
```

Run modules individually:

```bash
python run_rgpc.py --config configs/RGPC/rgpc.yaml
python run_tapocket_batch.py --config configs/TApocket/tapocket_batch.yaml
python run_docking.py --config configs/Docking/docking.yaml
python run_clusterscore.py --config configs/Scoring/cluster_score.yaml
python run_refinement.py --config configs/Refinement/refine_from_clusterscore.yaml
python run_metaboclip_bridge.py --config configs/MetaBoClip/metaboclip_bridge.yaml
```

## MetaboClip Integration

The active MetaboClip configuration is:

```text
configs/MetaBoClip/metaboclip_bridge.yaml
```

The bridge consumes the refined docking manifest:

```text
data/refined/docking_out/docking_result_manifest.csv
```

It preserves AImd metadata, prepares unified backend inputs, handles ligand role tables, calls `metaboclip.core.workflow.run_directory`, and writes AImd-compatible outputs under:

```text
data/metaboclip/results/
```

The stable final ranking path is:

```text
data/metaboclip/results/metaboclip_final_ranking.csv
```

Unified backend run artifacts are written under:

```text
data/metaboclip/unified_runs/
```

## Role Tables

`role_tables.mode` supports:

```text
existing
generate
auto
```

`existing` reuses prebuilt role tables. `generate` builds role tables from original ligand source files and prepared ligand PDBQT files. `auto` reuses existing role tables and generates only missing ones when source inputs are available.

Role generation uses the unified ligand role APIs and records:

```text
role_table_path
annotation_json_path
atom_map_json_path
```

Downstream catalytic geometry uses heavy atoms only. Ligand hydrogen atoms may be used by the unified role detector for functional-group recognition, but hydrogen coordinates are not used for downstream catalytic geometry.

## Documentation

Primary documents:

```text
docs/AImd_USER_MANUAL.md
docs/MODULE_INTERFACE_SPEC.md
docs/aimd_manifest_interface_spec.md
docs/ENGINEERING_CHECK_REPORT.md
```

## Third-Party Tools

AImd keeps external tool metadata under `third_party/`. To check command availability:

```bash
python third_party/check_tools.py --config third_party/tools.yaml --root .
```

PyMOL is not required for core unified MetaboClip scoring. It is optional for visualization/export paths or other modules that explicitly use it.
