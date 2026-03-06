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
  --lx 15 \
  --ly 15 \
  --dx 0.47 \
  --surface-bead P4 \
  --linker inputs/EPOXY.gro \
  --linker-group A 8 10 11 \
  --linker-group D 8 10 11 \
  --surface-linkers 12 \
  --merge A,B,C,D
