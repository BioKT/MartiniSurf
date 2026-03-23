#!/usr/bin/env bash
set -euo pipefail

EXAMPLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${EXAMPLE_DIR}/../.." && pwd)"
cd "${EXAMPLE_DIR}"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

python -m martinisurf \
  --complex-config input/complex_config.yaml \
  --surface-mode 4-1 \
  --surface-bead P4 \
  --dx 0.47 \
  --lx 15 \
  --ly 15 \
  --dist 0.8 \
  --substrate input/ETO.gro \
  --substrate-count 10 \
  --solvate \
  --ionize \
  --salt-conc 0.15 \
  --outdir Simulation_Files
