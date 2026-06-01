#!/usr/bin/env bash
set -euo pipefail

# Run from the TApocket project root.
tapocket check-layout --config configs/tapocket_template_v1.yaml
tapocket build-manifest --config configs/tapocket_template_v1.yaml
tapocket build-index --config configs/tapocket_template_v1.yaml --db all --force

tapocket run \
  --config configs/tapocket_template_v1.yaml \
  --query examples/query/example.pdb \
  --run-id example_001
