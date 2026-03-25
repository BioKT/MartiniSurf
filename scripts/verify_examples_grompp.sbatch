#!/bin/bash
#SBATCH --job-name=martinisurf_grompp
#SBATCH --qos=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CHECK_ROOT="${CHECK_ROOT:-${REPO_ROOT}/grompp_checks}"

if command -v module >/dev/null 2>&1; then
  module load GROMACS/2024-foss-2023a || true
fi

export GMX_BIN="${GMX_BIN:-gmx_mpi}"

if ! command -v "${GMX_BIN}" >/dev/null 2>&1; then
  echo "ERROR: ${GMX_BIN} not found. Load your GROMACS 2024 module or set GMX_BIN." >&2
  exit 1
fi

mkdir -p "${CHECK_ROOT}"

echo "============================================"
echo "MartiniSurf grompp verification job"
echo "Repo root : ${REPO_ROOT}"
echo "GROMACS   : ${GMX_BIN}"
echo "Outputs   : ${CHECK_ROOT}"
echo "Host      : $(hostname)"
echo "PWD       : ${PWD}"
echo "============================================"

cd "${REPO_ROOT}"
bash "${REPO_ROOT}/scripts/verify_examples_grompp.sh" "${CHECK_ROOT}"
