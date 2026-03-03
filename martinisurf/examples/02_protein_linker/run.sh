#!/usr/bin/env bash
set -euo pipefail

martinisurf \
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
