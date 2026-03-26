#!/usr/bin/env bash
set -euo pipefail

EXAMPLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${EXAMPLE_DIR}"

martinisurf \
  --dna \
  --dnatype ds-stiff \
  --pdb inputs/4C64.pdb \
  --surface-mode 4-1 \
  --surface-layers 1 \
  --lx 10 \
  --ly 10 \
  --dx 0.27 \
  --surface-bead C1 \
  --linker inputs/ALK.gro \
  --linker-group A 1 \
  --solvate \
  --ionize \
  --salt-conc 0.15 \
  --freeze-water-fraction 0.2 \
  --freeze-water-seed 42 \
  --merge A,B         
