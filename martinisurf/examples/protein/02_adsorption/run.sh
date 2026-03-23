#!/usr/bin/env bash
set -euo pipefail

EXAMPLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${EXAMPLE_DIR}/../.." && pwd)"
cd "${EXAMPLE_DIR}"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

python -m martinisurf \
  --pdb inputs/1RJW.pdb \
  --dssp \
  --moltype Protein \
  --surface-mode 4-1 \
  --surface-bead P4 \
  --dx 0.47 \
  --lx 15 \
  --ly 15 \
  --anchor A 8 10 11 \
  --anchor D 8 10 11 \
  --ads-mode \
  --dist 1.0 \
  --merge A,B,C,D
