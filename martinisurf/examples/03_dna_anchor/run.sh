#!/usr/bin/env bash
set -euo pipefail

martinisurf \
  --dna \
  --pdb inputs/4C64.pdb \
  --surface inputs/surface.gro \
  --anchor 1 1 \
  --anchor 2 24 \
  --dist 10 \
  --merge A,B
