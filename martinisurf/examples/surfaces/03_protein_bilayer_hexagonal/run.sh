#!/usr/bin/env bash
set -euo pipefail

EXAMPLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${EXAMPLE_DIR}/../../.." && pwd)"
cd "${EXAMPLE_DIR}"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

python -m martinisurf \
  --pdb 2A3D \
  --dssp \
  --moltype Protein \
  --surface-mode 2-1 \
  --lx 15 \
  --ly 15 \
  --surface-layers 2 \
  --surface-dist-z 0.382 \
  --surface-bead P4 C1 \
  --anchor A 1 \
  --dist 1.0 \
  --merge A \
  --outdir Simulation_Files
