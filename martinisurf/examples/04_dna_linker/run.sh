#!/usr/bin/env bash
set -euo pipefail

martinisurf \
  --dna \
  --pdb inputs/4C64.pdb \
  --surface inputs/surface.gro \
  --linker inputs/ALK.gro \
  --linker-group A 1 \
  --merge A,B
