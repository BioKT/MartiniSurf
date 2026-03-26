#!/bin/bash
#SBATCH --job-name=martinisurf_grompp
#SBATCH --qos=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-${PWD}}"
REPO_ROOT="${REPO_ROOT_OVERRIDE:-${SUBMIT_DIR}}"
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
echo "Submit dir: ${SUBMIT_DIR}"
echo "GROMACS   : ${GMX_BIN}"
echo "Outputs   : ${CHECK_ROOT}"
echo "Host      : $(hostname)"
echo "PWD       : ${PWD}"
echo "============================================"

cd "${REPO_ROOT}"
bash "${REPO_ROOT}/verify_examples_grompp.sh" "${CHECK_ROOT}"
