#!/bin/bash
#SBATCH --qos=regular
#SBATCH --job-name=05_example
#SBATCH --mem=200gb
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=ALL

set -euo pipefail

if command -v module >/dev/null 2>&1; then
  module load GROMACS/2024-foss-2023a || true
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}}"
cd "${BASE_DIR}"

SYSTEM_TAG="ADH"
REPLICA=1

if [[ -d "Simulation_Files/0_topology" && -d "Simulation_Files/1_mdp" && -d "Simulation_Files/2_system" ]]; then
  TOP_DIR="Simulation_Files/0_topology"
  MDP_DIR="Simulation_Files/1_mdp"
  SYS_DIR="Simulation_Files/2_system"
elif [[ -d "Simulation/0_topology" && -d "Simulation/1_mdp" && -d "Simulation/2_system" ]]; then
  TOP_DIR="Simulation/0_topology"
  MDP_DIR="Simulation/1_mdp"
  SYS_DIR="Simulation/2_system"
elif [[ -d "0_topology" && -d "1_mdp" && -d "2_system" ]]; then
  TOP_DIR="0_topology"
  MDP_DIR="1_mdp"
  SYS_DIR="2_system"
else
  echo "ERROR: no simulation folder set found from BASE_DIR=${BASE_DIR}" >&2
  echo "Expected one of:" >&2
  echo "  - ${BASE_DIR}/Simulation_Files/{0_topology,1_mdp,2_system}" >&2
  echo "  - ${BASE_DIR}/Simulation/{0_topology,1_mdp,2_system}" >&2
  echo "  - ${BASE_DIR}/{0_topology,1_mdp,2_system}" >&2
  exit 1
fi

pick_first_existing() {
  for candidate in "$@"; do
    if [[ -f "${candidate}" ]]; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

must_find_file() {
  local label="$1"
  shift
  local found
  found="$(pick_first_existing "$@" || true)"
  if [[ -z "${found}" ]]; then
    echo "ERROR: could not find ${label}. Tried:" >&2
    for candidate in "$@"; do
      echo "  - ${candidate}" >&2
    done
    exit 1
  fi
  echo "${found}"
}

TOPOLOGY="$(must_find_file "topology" \
  "${TOP_DIR}/system_final_res.top" \
  "${TOP_DIR}/system_final.top" \
  "${TOP_DIR}/system_res.top" \
  "${TOP_DIR}/system_anchor.top" \
  "${TOP_DIR}/system.top")"

INPUT_STRUCTURE="$(must_find_file "input structure" \
  "${SYS_DIR}/system_final.gro" \
  "${SYS_DIR}/final_system.gro" \
  "${SYS_DIR}/system.gro" \
  "${SYS_DIR}/immobilized_system.gro")"

INDEX_FILE="${TOP_DIR}/index.ndx"
INDEX_ARGS=()
if [[ -f "${INDEX_FILE}" ]]; then
  INDEX_ARGS=(-n "${INDEX_FILE}")
fi

pick_mdp() {
  local stage="$1"
  must_find_file "${stage} mdp" \
    "${MDP_DIR}/${stage}.mdp" \
    "${MDP_DIR}/${stage}_dna.mdp"
}

MIN_MDP="$(pick_mdp minimization)"
NVT_MDP="$(pick_mdp nvt)"
NPT_MDP="$(pick_mdp npt)"
DEP_MDP="$(pick_mdp deposition)"
PROD_MDP="$(pick_mdp production)"

if command -v gmx_mpi >/dev/null 2>&1; then
  GMX_BIN="gmx_mpi"
elif command -v gmx >/dev/null 2>&1; then
  GMX_BIN="gmx"
else
  echo "Error: no GROMACS binary found (gmx_mpi/gmx)." >&2
  exit 1
fi

REPLICA_PADDED="$(printf "%02d" "${REPLICA}")"
MIN_NAME="${SYSTEM_TAG}_min_r${REPLICA_PADDED}"
NVT_NAME="${SYSTEM_TAG}_nvt_r${REPLICA_PADDED}"
NPT_NAME="${SYSTEM_TAG}_npt_r${REPLICA_PADDED}"
DEP_NAME="${SYSTEM_TAG}_dep_r${REPLICA_PADDED}"
PROD_NAME="${SYSTEM_TAG}_prod_r${REPLICA_PADDED}"

NTOMP="${SLURM_CPUS_PER_TASK:-1}"

run_mdrun() {
  if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v srun >/dev/null 2>&1; then
    srun "${GMX_BIN}" mdrun -ntomp "${NTOMP}" -deffnm "$1" -v
  else
    "${GMX_BIN}" mdrun -ntomp "${NTOMP}" -deffnm "$1" -v
  fi
}

echo "===================================================="
echo "MartiniSurf Job"
echo "System   : ${SYSTEM_TAG}"
echo "Replica  : r${REPLICA_PADDED}"
echo "Base dir : ${BASE_DIR}"
echo "Topology : ${TOPOLOGY}"
echo "Input GRO: ${INPUT_STRUCTURE}"
echo "GMX      : ${GMX_BIN}"
echo "===================================================="

"${GMX_BIN}" grompp \
  -p "${TOPOLOGY}" \
  -f "${MIN_MDP}" \
  -c "${INPUT_STRUCTURE}" \
  -o "${MIN_NAME}.tpr" \
  -maxwarn 3
run_mdrun "${MIN_NAME}"

"${GMX_BIN}" grompp \
  -p "${TOPOLOGY}" \
  -f "${NVT_MDP}" \
  -c "${MIN_NAME}.gro" \
  -r "${MIN_NAME}.gro" \
  -o "${NVT_NAME}.tpr" \
  "${INDEX_ARGS[@]}" \
  -maxwarn 3
run_mdrun "${NVT_NAME}"

"${GMX_BIN}" grompp \
  -p "${TOPOLOGY}" \
  -f "${NPT_MDP}" \
  -c "${NVT_NAME}.gro" \
  -r "${NVT_NAME}.gro" \
  -t "${NVT_NAME}.cpt" \
  -o "${NPT_NAME}.tpr" \
  "${INDEX_ARGS[@]}" \
  -maxwarn 3
run_mdrun "${NPT_NAME}"

"${GMX_BIN}" grompp \
  -p "${TOPOLOGY}" \
  -f "${DEP_MDP}" \
  -c "${NPT_NAME}.gro" \
  -r "${NPT_NAME}.gro" \
  -t "${NPT_NAME}.cpt" \
  -o "${DEP_NAME}.tpr" \
  "${INDEX_ARGS[@]}" \
  -maxwarn 3
run_mdrun "${DEP_NAME}"

"${GMX_BIN}" grompp \
  -p "${TOPOLOGY}" \
  -f "${PROD_MDP}" \
  -c "${DEP_NAME}.gro" \
  -r "${DEP_NAME}.gro" \
  -t "${DEP_NAME}.cpt" \
  -o "${PROD_NAME}.tpr" \
  "${INDEX_ARGS[@]}" \
  -maxwarn 3
run_mdrun "${PROD_NAME}"

echo "===================================================="
echo "Simulation completed successfully."
echo "Final output prefix: ${PROD_NAME}"
echo "===================================================="
