#!/usr/bin/env bash
set -euo pipefail

EXAMPLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${EXAMPLE_DIR}/../.." && pwd)"
cd "${EXAMPLE_DIR}"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

python -m martinisurf \
  --dna \
  --pdb inputs/4C64.pdb \
  --surface inputs/surface.gro \
  --linker inputs/ALK.gro \
  --linker-group A 1 \
  --merge A,B
