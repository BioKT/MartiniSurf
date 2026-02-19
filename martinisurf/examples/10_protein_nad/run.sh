#!/usr/bin/env bash
set -euo pipefail

martinisurf \
  --complex-config input/complex_config.yaml \
  --lx 20 \
  --ly 20 \
  --dist 10 \
  --substrate input/ETO.gro \
  --substrate-count 10 \
  --solvate \
  --ionize \
  --salt-conc 0.15 \
  --outdir Simulation_Files
