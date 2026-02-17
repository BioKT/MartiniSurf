#!/usr/bin/env bash
set -euo pipefail

mkdir -p outputs

# 2-1 mode (hex mapping)
python -m martinisurf.surface_builder \
  --mode 2-1 \
  --lx 8 \
  --ly 8 \
  --dx 0.47 \
  --bead C1 \
  --resname SRF \
  --charge 0.0 \
  --output outputs/surface_2to1

# 4-1 mode (honeycomb carbon)
python -m martinisurf.surface_builder \
  --mode 4-1 \
  --lx 8 \
  --ly 8 \
  --dx 0.142 \
  --layers 1 \
  --dist-z 0.335 \
  --bead C1 \
  --resname SRF \
  --charge 0.0 \
  --output outputs/surface_4to1

echo
ls -lh outputs/*.gro outputs/*.itp
