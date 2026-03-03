#!/usr/bin/env bash
set -euo pipefail

martinisurf \
  --pdb inputs/1RJW.pdb \
  --dssp \
  --go  \
  --moltype Protein \
  --surface-mode 4-1 \
  --surface-bead P4 \
  --dx 0.47 \
  --lx 15 \
  --ly 15 \
  --anchor A 8 10 11 \
  --anchor D 8 10 11 \
  --dist 10 \
  --solvate \
  --ionize \
  --salt-conc 0.15 \
  --merge A,B,C,D
