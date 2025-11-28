#!/usr/bin/env bash

set -e

echo "===================================="
echo "  🧼 Running full code formatter..."
echo "===================================="

TARGET="surfmartini"

# 1. Remove unused imports & unused variables
echo "🔹 Running autoflake..."
autoflake --remove-all-unused-imports \
          --remove-unused-variables \
          --in-place -r $TARGET

# 2. Fix indentation, whitespace, line breaks, basic PEP8
echo "🔹 Running autopep8..."
autopep8 --max-line-length 88 \
         --aggressive --aggressive \
         --in-place -r $TARGET

# 3. Sort imports (PEP8 + logical groups)
echo "🔹 Running isort..."
isort $TARGET

# 4. Apply Black for consistent styling
echo "🔹 Running black..."
black --line-length 88 $TARGET

echo "===================================="
echo "  ✔ All formatting complete!"
echo "===================================="

