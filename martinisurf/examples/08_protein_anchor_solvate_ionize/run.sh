#!/usr/bin/env bash
set -euo pipefail

martinisurf \
  --pdb inputs/1RJW.pdb \
  --dssp \
  --moltype Protein \
  --surface inputs/surface.gro \
  --anchor 1 8 10 11 \
  --anchor 2 1025 1027 1028 \
  --dist 10 \
  --solvate \
  --ionize \
  --salt-conc 0.15 \
  --merge A,B,C,D
