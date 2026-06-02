#!/usr/bin/env bash
set -euo pipefail

# Run from AImd root:
#   cd AImd

python scripts/migrate_data_layout.py --root .
python validate_aimd_layout.py --root .

# Deterministic MetaboClip bridge smoke test
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:metaboclip_unified \
  python -m pytest tests/test_metaboclip_bridge_smoke.py -q

# 1. Protein structure clustering
python run_rgpc.py --config configs/RGPC/rgpc.yaml

# 2. Pocket prediction on original protein structures
python run_tapocket_batch.py --config configs/TApocket/tapocket_batch.yaml

# 3. Broad docking
python run_docking.py --config configs/Docking/docking.yaml

# 4. ClusterScore
python run_clusterscore.py --config configs/Scoring/cluster_score.yaml

# 5. Refined top-cluster docking with ensemble/cofactor settings
python run_refinement.py --config configs/Refinement/refine_from_clusterscore.yaml

# 6. Unified MetaboClip role handling, catalytic geometry scoring, and AImd result translation
python run_metaboclip_bridge.py --config configs/MetaBoClip/metaboclip_bridge.yaml

# Or run the full workflow:
python run_full_iterative_metaboclip.py --config configs/workflows/full_iterative_metaboclip.yaml
