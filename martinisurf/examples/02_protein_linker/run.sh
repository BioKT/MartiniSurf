#!/usr/bin/env bash
set -euo pipefail

martinisurf \
  --pdb inputs/1RJW.pdb \
  --go \
  --moltype Protein \
  --surface inputs/surface.gro \
  --linker inputs/linker.gro \
  --linker-group 1 8 10 11 \
  --linker-group 2 1025 1027 1028 \
  --merge A,B,C,D
