#!/bin/bash
#SBATCH --qos=regular
#SBATCH --job-name=12_soft
#SBATCH --mem=200gb
#SBATCH --cpus-per-task=4
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

SYSTEM_TAG="DNA"
REPLICA=1

# Lista de cargas corregida (sin comas erróneas)
CHARGES=(0 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1)

if [[ -d "Simulation_Files/0_topology" && -d "Simulation_Files/1_mdp" && -d "Simulation_Files/2_system" ]]; then
  SRC_TOP_DIR="Simulation_Files/0_topology"
  SRC_MDP_DIR="Simulation_Files/1_mdp"
  SRC_SYS_DIR="Simulation_Files/2_system"
elif [[ -d "Simulation/0_topology" && -d "Simulation/1_mdp" && -d "Simulation/2_system" ]]; then
  SRC_TOP_DIR="Simulation/0_topology"
  SRC_MDP_DIR="Simulation/1_mdp"
  SRC_SYS_DIR="Simulation/2_system"
elif [[ -d "0_topology" && -d "1_mdp" && -d "2_system" ]]; then
  SRC_TOP_DIR="0_topology"
  SRC_MDP_DIR="1_mdp"
  SRC_SYS_DIR="2_system"
else
  echo "ERROR: no simulation folder set found from BASE_DIR=${BASE_DIR}" >&2
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

first_layer_atom_count() {
  local gro_path="$1"
  python - "$gro_path" <<'PY'
from collections import Counter
from pathlib import Path
import sys

gro_path = Path(sys.argv[1])
lines = gro_path.read_text().splitlines()[2:-1]
if not lines:
    raise SystemExit("surface.gro has no atoms")

z_values = [round(float(line[36:44]), 3) for line in lines]
counts = Counter(z_values)
if len(counts) == 1:
    layer_z, atom_count = next(iter(counts.items()))
    print(
        f"WARN: surface.gro contains a single layer at z={layer_z:.3f}; "
        "the charge scan will apply the selected charge to the whole surface.",
        file=sys.stderr,
    )
    print(atom_count)
    raise SystemExit(0)
if len(counts) < 1:
    raise SystemExit("surface.gro does not contain any detectable layer")

layer_z = min(counts)
print(counts[layer_z])
PY
}

patch_surface_charge() {
  local itp_path="$1"
  local first_layer_atoms="$2"
  local charge_value="$3"
  local tmp_path
  tmp_path="$(mktemp)"

  awk -v n_layer="${first_layer_atoms}" -v c="${charge_value}" '
    BEGIN { in_atoms = 0 }
    {
      if ($0 ~ /^[[:space:]]*\[/) {
        section = tolower($0)
        gsub(/[[:space:]]/, "", section)
        in_atoms = (section == "[atoms]")
        print
        next
      }
      if (in_atoms && $1 ~ /^[0-9]+$/) {
        charge = $7 + 0.0
        if ($1 <= n_layer) {
          charge = c
        }
        printf "%6d %6s %6d %-6s %6s %6d %8.4f\n", $1, $2, $3, $4, $5, $6, charge
        next
      }
      print
    }
  ' "${itp_path}" > "${tmp_path}"

  mv "${tmp_path}" "${itp_path}"
}

if command -v gmx_mpi >/dev/null 2>&1; then
  GMX_BIN="gmx_mpi"
elif command -v gmx >/dev/null 2>&1; then
  GMX_BIN="gmx"
else
  echo "Error: no GROMACS binary found (gmx_mpi/gmx)." >&2
  exit 1
fi

NTOMP="${SLURM_CPUS_PER_TASK:-1}"

run_mdrun() {
  local deffnm="$1"
  if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v srun >/dev/null 2>&1; then
    srun "${GMX_BIN}" mdrun -ntomp "${NTOMP}" -deffnm "${deffnm}" -v
  else
    "${GMX_BIN}" mdrun -ntomp "${NTOMP}" -deffnm "${deffnm}" -v
  fi
}

REPLICA_PADDED="$(printf "%02d" "${REPLICA}")"
SCAN_ROOT="${BASE_DIR}/charge_scan_runs"
mkdir -p "${SCAN_ROOT}"

for CHARGE in "${CHARGES[@]}"; do
  # Mantener el valor real sin redondearlo a 2 decimales
  CHARGE_LABEL="${CHARGE}"
  CHARGE_TAG="q${CHARGE_LABEL//./p}"
  CASE_DIR="${SCAN_ROOT}/${CHARGE_TAG}"

  echo "===================================================="
  echo "Preparing case ${CHARGE_TAG} (first-layer surface charge = ${CHARGE_LABEL})"
  echo "Case dir: ${CASE_DIR}"
  echo "===================================================="

  mkdir -p "${CASE_DIR}"
  rm -rf "${CASE_DIR}/0_topology" "${CASE_DIR}/1_mdp" "${CASE_DIR}/2_system"
  cp -a "${SRC_TOP_DIR}" "${CASE_DIR}/0_topology"
  cp -a "${SRC_MDP_DIR}" "${CASE_DIR}/1_mdp"
  cp -a "${SRC_SYS_DIR}" "${CASE_DIR}/2_system"

  TOP_DIR="${CASE_DIR}/0_topology"
  MDP_DIR="${CASE_DIR}/1_mdp"
  SYS_DIR="${CASE_DIR}/2_system"

  SURFACE_ITP="$(must_find_file "surface.itp" \
    "${TOP_DIR}/system_itp/surface.itp" \
    "${TOP_DIR}/surface.itp")"
  SURFACE_GRO="$(must_find_file "surface.gro" \
    "${SYS_DIR}/surface.gro")"
  FIRST_LAYER_ATOMS="$(first_layer_atom_count "${SURFACE_GRO}")"
  patch_surface_charge "${SURFACE_ITP}" "${FIRST_LAYER_ATOMS}" "${CHARGE_LABEL}"

  EQUIL_TOPOLOGY="$(must_find_file "equilibration topology" \
    "${TOP_DIR}/system_final.top" \
    "${TOP_DIR}/system.top")"

  PROD_TOPOLOGY="$(pick_first_existing \
    "${TOP_DIR}/system_final_res.top" \
    "${TOP_DIR}/system_res.top" \
    "${TOP_DIR}/system_anchor.top" || true)"
  if [[ -z "${PROD_TOPOLOGY}" ]]; then
    PROD_TOPOLOGY="${EQUIL_TOPOLOGY}"
    echo "WARN: production restricted topology not found; using equilibration topology for production."
  fi

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

  MIN_MDP="$(must_find_file "minimization mdp" \
    "${MDP_DIR}/minimization_dna.mdp" \
    "${MDP_DIR}/minimization.mdp")"

  NVT_MDP="$(must_find_file "nvt mdp" \
    "${MDP_DIR}/nvt_dna.mdp" \
    "${MDP_DIR}/nvt.mdp")"

  DEP_MDP="$(must_find_file "deposition mdp" \
    "${MDP_DIR}/deposition_dna.mdp" \
    "${MDP_DIR}/deposition.mdp")"

  PROD_MDP="$(must_find_file "production mdp" \
    "${MDP_DIR}/production_dna.mdp" \
    "${MDP_DIR}/production.mdp")"

  MIN_NAME="${SYSTEM_TAG}_${CHARGE_TAG}_min_r${REPLICA_PADDED}"
  NVT_NAME="${SYSTEM_TAG}_${CHARGE_TAG}_nvt_r${REPLICA_PADDED}"
  DEP_NAME="${SYSTEM_TAG}_${CHARGE_TAG}_dep_r${REPLICA_PADDED}"
  PROD_NAME="${SYSTEM_TAG}_${CHARGE_TAG}_prod_r${REPLICA_PADDED}"

  pushd "${CASE_DIR}" >/dev/null

  "${GMX_BIN}" grompp \
    -p "${EQUIL_TOPOLOGY}" \
    -f "${MIN_MDP}" \
    -c "${INPUT_STRUCTURE}" \
    -r "${INPUT_STRUCTURE}" \
    -o "${MIN_NAME}.tpr" \
    -maxwarn 3
  run_mdrun "${MIN_NAME}"

  "${GMX_BIN}" grompp \
    -p "${EQUIL_TOPOLOGY}" \
    -f "${NVT_MDP}" \
    -c "${MIN_NAME}.gro" \
    -r "${MIN_NAME}.gro" \
    -o "${NVT_NAME}.tpr" \
    "${INDEX_ARGS[@]}" \
    -maxwarn 3
  run_mdrun "${NVT_NAME}"

  "${GMX_BIN}" grompp \
    -p "${EQUIL_TOPOLOGY}" \
    -f "${DEP_MDP}" \
    -c "${NVT_NAME}.gro" \
    -r "${NVT_NAME}.gro" \
    -t "${NVT_NAME}.cpt" \
    -o "${DEP_NAME}.tpr" \
    "${INDEX_ARGS[@]}" \
    -maxwarn 3
  run_mdrun "${DEP_NAME}"

  "${GMX_BIN}" grompp \
    -p "${PROD_TOPOLOGY}" \
    -f "${PROD_MDP}" \
    -c "${DEP_NAME}.gro" \
    -r "${DEP_NAME}.gro" \
    -t "${DEP_NAME}.cpt" \
    -o "${PROD_NAME}.tpr" \
    "${INDEX_ARGS[@]}" \
    -maxwarn 3
  run_mdrun "${PROD_NAME}"

  popd >/dev/null

  echo "Completed ${CHARGE_TAG}. Final prefix: ${CASE_DIR}/${PROD_NAME}"
done

echo "===================================================="
echo "Charge scan completed successfully."
echo "Results directory: ${SCAN_ROOT}"
echo "===================================================="
