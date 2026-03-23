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
  --surface-mode cnt \
  --lx 1 \
  --ly 1 \
  --cnt-numrings 24 \
  --anchor A 1 \
  --dist 1.0 \
  --merge A \
  --outdir Simulation_Files
