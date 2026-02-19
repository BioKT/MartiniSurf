#!/usr/bin/env bash
set -euo pipefail

martinisurf \
  --complex-config input/complex_config.yaml \
  --surface-mode 4-1 \
  --surface-bead P4 \
  --dx 0.47 \
  --lx 15 \
  --ly 15 \
  --dist 10 \
  --substrate input/ETO.gro \
  --substrate-count 10 \
  --solvate \
  --ionize \
  --salt-conc 0.15 \
  --outdir Simulation_Files
