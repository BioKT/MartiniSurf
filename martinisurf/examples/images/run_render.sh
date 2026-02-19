#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXAMPLES_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="$SCRIPT_DIR"

if ! command -v vmd >/dev/null 2>&1; then
  echo "Error: vmd not found in PATH"
  exit 1
fi

vmd -dispdev text \
  -e "$SCRIPT_DIR/render_examples_vmd.tcl" \
  -args "$EXAMPLES_DIR" "$OUT_DIR"

echo "Render completed. Images in: $OUT_DIR"
