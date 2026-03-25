#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'EOF'
Usage:
  bash scripts/verify_examples_grompp.sh [CHECK_ROOT]

Description:
  Runs a GROMACS grompp consistency check for every generated
  martinisurf/examples/**/Simulation_Files directory in the repository.

Arguments:
  CHECK_ROOT   Optional output directory for logs and generated .tpr files.
               Default: ./grompp_checks

Environment:
  GMX_BIN      Optional GROMACS binary to use, for example:
               GMX_BIN=gmx_mpi

Notes:
  - Run this on the HPC or environment where your real GROMACS version is available.
  - The script uses the generated GRO file as both -c and -r so restrained topologies
    can be checked without requiring an extra reference file.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  show_help
  exit 0
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CHECK_ROOT="${1:-${REPO_ROOT}/grompp_checks}"

if [[ -z "${GMX_BIN:-}" ]]; then
  if command -v gmx_mpi >/dev/null 2>&1; then
    GMX_BIN="gmx_mpi"
  elif command -v gmx >/dev/null 2>&1; then
    GMX_BIN="gmx"
  else
    echo "ERROR: no GROMACS binary found. Set GMX_BIN or load GROMACS first." >&2
    exit 1
  fi
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

mkdir -p "${CHECK_ROOT}"

total=0
passed=0
failed=0

while IFS= read -r simdir; do
  total=$((total + 1))

  top_dir="${simdir}/0_topology"
  mdp_dir="${simdir}/1_mdp"
  sys_dir="${simdir}/2_system"

  top_path="$(pick_first_existing \
    "${top_dir}/system_final_res.top" \
    "${top_dir}/system_final.top" \
    "${top_dir}/system_res.top" \
    "${top_dir}/system_anchor.top" \
    "${top_dir}/system.top" || true)"

  gro_path="$(pick_first_existing \
    "${sys_dir}/system_final.gro" \
    "${sys_dir}/final_system.gro" \
    "${sys_dir}/immobilized_system.gro" \
    "${sys_dir}/system.gro" || true)"

  mdp_path="$(pick_first_existing \
    "${mdp_dir}/minimization_dna.mdp" \
    "${mdp_dir}/minimization.mdp" \
    "${mdp_dir}/nvt_dna.mdp" \
    "${mdp_dir}/nvt.mdp" || true)"

  rel_name="${simdir#${REPO_ROOT}/}"
  safe_name="${rel_name//\//__}"
  log_path="${CHECK_ROOT}/${safe_name}.log"
  out_tpr="${CHECK_ROOT}/${safe_name}.tpr"
  out_mdp="${CHECK_ROOT}/${safe_name}_mdout.mdp"

  if [[ -z "${top_path}" || -z "${gro_path}" || -z "${mdp_path}" ]]; then
    echo "[FAIL] ${rel_name}" | tee "${log_path}"
    {
      echo "Missing required input."
      echo "TOP: ${top_path:-missing}"
      echo "GRO: ${gro_path:-missing}"
      echo "MDP: ${mdp_path:-missing}"
    } >> "${log_path}"
    failed=$((failed + 1))
    continue
  fi

  index_args=()
  if [[ -f "${top_dir}/index.ndx" ]]; then
    index_args=(-n "${top_dir}/index.ndx")
  fi

  echo "[CHECK] ${rel_name}"

  if "${GMX_BIN}" grompp \
    -p "${top_path}" \
    -f "${mdp_path}" \
    -c "${gro_path}" \
    -r "${gro_path}" \
    -o "${out_tpr}" \
    -po "${out_mdp}" \
    "${index_args[@]}" \
    -maxwarn 3 > "${log_path}" 2>&1; then
    echo "[OK] ${rel_name}"
    passed=$((passed + 1))
  else
    echo "[FAIL] ${rel_name}"
    failed=$((failed + 1))
  fi
done < <(find "${REPO_ROOT}/martinisurf/examples" -type d -name Simulation_Files | sort)

echo
echo "GROMPP check summary"
echo "  GROMACS : ${GMX_BIN}"
echo "  Outputs : ${CHECK_ROOT}"
echo "  Total   : ${total}"
echo "  Passed  : ${passed}"
echo "  Failed  : ${failed}"

if [[ "${failed}" -ne 0 ]]; then
  exit 1
fi
