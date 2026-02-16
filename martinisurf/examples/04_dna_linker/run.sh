#!/usr/bin/env bash
set -euo pipefail

martinisurf \
  --dna \
  --pdb inputs/4C64.pdb \
  --surface inputs/surface.gro \
  --linker inputs/linker.gro \
  --linker-group 1 1 \
  --merge A,B
