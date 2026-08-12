from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


DEFAULT_OUTPUT_TAG = "Protein MD"
DEFAULT_XTC_WRITE_EVERY_PS = 10.0
DEFAULT_GROMPP_MAXWARN = 3
DEFAULT_STAGE_SETTINGS = {
    "nvt": {"enabled": True, "dt_ps": 0.001, "time_ns": 0.02},
    "npt": {"enabled": True, "dt_ps": 0.005, "time_ns": 0.02},
    "deposition": {"enabled": True, "dt_ps": 0.01, "time_ns": 0.05},
    "production": {"enabled": True, "dt_ps": 0.01, "time_ns": 0.1},
}
STAGE_ORDER = ["nvt", "npt", "deposition", "production"]


@dataclass
class ShortMDStage:
    name: str
    enabled: bool
    dt_ps: float
    time_ns: float


@dataclass
class ShortMDConfig:
    output_tag: str = DEFAULT_OUTPUT_TAG
    xtc_write_every_ps: float = DEFAULT_XTC_WRITE_EVERY_PS
    grompp_maxwarn: int = DEFAULT_GROMPP_MAXWARN
    stages: tuple[ShortMDStage, ...] = ()


@dataclass
class ShortMDStageResult:
    name: str
    grompp_elapsed_s: float
    mdrun_elapsed_s: float
    ns_day: float | None = None
    hour_ns: float | None = None
    dt_ps: float | None = None
    time_ns: float | None = None
    nsteps: int | None = None
    saved_frame_dt_ps: float | None = None


@dataclass
class ShortMDResult:
    returncode: int
    stdout: str
    stderr: str
    work_dir: Path
    last_gro: Path | None
    last_tpr: Path | None
    last_xtc: Path | None
    stage_results: list[ShortMDStageResult]
    command_log: list[str]


def detect_gmx(extra_tool_dirs: list[Path] | None = None) -> str | None:
    path = _path_with_tools(extra_tool_dirs)
    return shutil.which("gmx", path=path) or shutil.which("gmx_mpi", path=path)


def selected_stages(config: ShortMDConfig) -> list[ShortMDStage]:
    return [stage for stage in config.stages if stage.enabled]


def validate_stage_order(stages: list[ShortMDStage]) -> list[str]:
    names = [stage.name for stage in stages]
    errors: list[str] = []
    if "npt" in names and "nvt" not in names:
        errors.append("NPT requires NVT to be enabled.")
    if ("deposition" in names or "production" in names) and "npt" not in names:
        errors.append("Deposition and Production require NPT to be enabled.")
    return errors


def preview_rows(config: ShortMDConfig) -> list[dict[str, object]]:
    rows = []
    for stage in selected_stages(config):
        nsteps = safe_nsteps(stage.time_ns, stage.dt_ps)
        xtc_stride = max(1, int(round(config.xtc_write_every_ps / stage.dt_ps)))
        rows.append(
            {
                "Stage": stage.name.upper(),
                "dt (ps)": stage.dt_ps,
                "time (ns)": stage.time_ns,
                "nsteps": nsteps,
                "saved frame dt (ps)": round(xtc_stride * stage.dt_ps, 6),
            }
        )
    return rows


def run_short_md(sim_root: Path, config: ShortMDConfig, repo_root: Path, extra_tool_dirs: list[Path] | None = None) -> ShortMDResult:
    top_dir = sim_root / "0_topology"
    mdp_dir = sim_root / "1_mdp"
    system_dir = sim_root / "2_system"
    work_dir = sim_root / "3_colab_md"
    work_dir.mkdir(parents=True, exist_ok=True)

    stdout_parts: list[str] = []
    stderr_parts: list[str] = []
    command_log: list[str] = []
    stage_results: list[ShortMDStageResult] = []

    def emit(text: str) -> None:
        stdout_parts.append(text.rstrip() + "\n")

    if not sim_root.exists():
        raise FileNotFoundError(f"Missing folder: {sim_root}. Generate the system first.")

    stages = selected_stages(config)
    errors = validate_stage_order(stages)
    if errors:
        raise ValueError(" ".join(errors))
    if not stages:
        raise ValueError("No protocol stage selected.")

    gmx = detect_gmx(extra_tool_dirs)
    if not gmx:
        raise FileNotFoundError("GROMACS not found. Install gromacs or ensure gmx/gmx_mpi is in PATH.")

    protocol_suffix = "_dna" if sorted(mdp_dir.glob("*_dna.mdp")) else ""
    initial_gro = _pick_existing(
        [
            system_dir / "final_system.gro",
            system_dir / "system_final.gro",
            system_dir / "immobilized_system.gro",
            system_dir / "system.gro",
            system_dir / "solvated_system.gro",
        ]
    )
    if initial_gro is None:
        raise FileNotFoundError("Could not find an initial GRO file in Simulation_Files/2_system.")

    equilibrium_top = _pick_existing([top_dir / "system_final.top", top_dir / "system.top"])
    if equilibrium_top is None:
        raise FileNotFoundError("Missing required equilibrium topology file (system_final.top or system.top).")

    production_top = _pick_existing([top_dir / "system_final_res.top", top_dir / "system_res.top"])
    if production_top is None:
        production_top = equilibrium_top
        emit("Warning: restrained production topology not found; production will use equilibrium topology.")

    index_ndx = top_dir / "index.ndx"
    if not index_ndx.exists():
        raise FileNotFoundError(f"Missing required index file: {index_ndx}")

    emit("Selected stages: " + " -> ".join(["minimization"] + [stage.name for stage in stages]))
    emit(f"Initial GRO    : {initial_gro}")
    emit(f"Topology(eq)   : {equilibrium_top}")
    emit(f"Topology(prod) : {production_top}")
    emit(f"Index          : {index_ndx}")
    emit("Execution mode : CPU")

    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root) + os.pathsep + env.get("PYTHONPATH", "")
    env["PATH"] = _path_with_tools(extra_tool_dirs)

    prev_gro = initial_gro
    prev_cpt: Path | None = None
    total_stages = 1 + len(stages)

    try:
        emit(f"\n[Stage 1/{total_stages}] MINIMIZATION")
        em_src = _mdp_path_for(mdp_dir, "minimization", protocol_suffix)
        em_dst = work_dir / f"{_safe_tag(config.output_tag)}_minimization.mdp"
        em_dst.write_text(em_src.read_text())
        em_base = work_dir / f"{_safe_tag(config.output_tag)}_minimization"
        em_cmd = [
            gmx,
            "grompp",
            "-f",
            str(em_dst),
            "-c",
            str(prev_gro),
            "-r",
            str(prev_gro),
            "-p",
            str(equilibrium_top),
            "-n",
            str(index_ndx),
            "-o",
            str(em_base.with_suffix(".tpr")),
            "-maxwarn",
            str(int(config.grompp_maxwarn)),
        ]
        _, grompp_elapsed = _run_cmd(em_cmd, top_dir, env, stdout_parts, stderr_parts, command_log)
        em_res, mdrun_elapsed = _run_cmd([gmx, "mdrun", "-deffnm", str(em_base)], work_dir, env, stdout_parts, stderr_parts, command_log)
        ns_day, hour_ns = extract_mdrun_performance(em_res.stdout + "\n" + em_res.stderr)
        stage_results.append(ShortMDStageResult("minimization", grompp_elapsed, mdrun_elapsed, ns_day, hour_ns))
        prev_gro = em_base.with_suffix(".gro")
        prev_cpt = em_base.with_suffix(".cpt")

        for index, stage in enumerate(stages, start=2):
            emit(f"\n[Stage {index}/{total_stages}] {stage.name.upper()} | dt={stage.dt_ps} ps | time={stage.time_ns} ns")
            nsteps = safe_nsteps(stage.time_ns, stage.dt_ps)
            xtc_stride = max(1, int(round(config.xtc_write_every_ps / stage.dt_ps)))
            mdp_src = _mdp_path_for(mdp_dir, stage.name, protocol_suffix)
            mdp_dst = work_dir / f"{_safe_tag(config.output_tag)}_{stage.name}.mdp"
            mdp_dst.write_text(patch_mdp_text(mdp_src.read_text(), stage.dt_ps, nsteps, xtc_stride))
            base = work_dir / f"{_safe_tag(config.output_tag)}_{stage.name}"
            stage_top = production_top if stage.name == "production" else equilibrium_top
            grompp_cmd = [
                gmx,
                "grompp",
                "-f",
                str(mdp_dst),
                "-c",
                str(prev_gro),
                "-r",
                str(prev_gro),
                "-p",
                str(stage_top),
                "-n",
                str(index_ndx),
                "-o",
                str(base.with_suffix(".tpr")),
                "-maxwarn",
                str(int(config.grompp_maxwarn)),
            ]
            if prev_cpt and prev_cpt.exists():
                grompp_cmd.extend(["-t", str(prev_cpt)])
            _, grompp_elapsed = _run_cmd(grompp_cmd, top_dir, env, stdout_parts, stderr_parts, command_log)
            mdrun_res, mdrun_elapsed = _run_cmd([gmx, "mdrun", "-deffnm", str(base)], work_dir, env, stdout_parts, stderr_parts, command_log)
            ns_day, hour_ns = extract_mdrun_performance(mdrun_res.stdout + "\n" + mdrun_res.stderr)
            prev_gro = base.with_suffix(".gro")
            prev_cpt = base.with_suffix(".cpt")
            stage_results.append(
                ShortMDStageResult(
                    stage.name,
                    grompp_elapsed,
                    mdrun_elapsed,
                    ns_day,
                    hour_ns,
                    stage.dt_ps,
                    stage.time_ns,
                    nsteps,
                    xtc_stride * stage.dt_ps,
                )
            )
    except subprocess.CalledProcessError as exc:
        return ShortMDResult(
            exc.returncode,
            "".join(stdout_parts),
            "".join(stderr_parts),
            work_dir,
            prev_gro if "prev_gro" in locals() else None,
            None,
            None,
            stage_results,
            command_log,
        )

    last_stage = stages[-1].name
    last_base = work_dir / f"{_safe_tag(config.output_tag)}_{last_stage}"
    emit("Short MD finished successfully.")
    emit(f"Last stage GRO: {last_base.with_suffix('.gro')}")
    return ShortMDResult(
        0,
        "".join(stdout_parts),
        "".join(stderr_parts),
        work_dir,
        last_base.with_suffix(".gro"),
        last_base.with_suffix(".tpr"),
        last_base.with_suffix(".xtc"),
        stage_results,
        command_log,
    )


def patch_mdp_text(text: str, dt_ps: float, nsteps: int, nstxoutc: int) -> str:
    updates = {
        "dt": f"{dt_ps:.6f}",
        "nsteps": str(int(nsteps)),
        "nstxout-compressed": str(int(nstxoutc)),
    }
    seen: set[str] = set()
    out_lines: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(";") or "=" not in stripped:
            out_lines.append(line)
            continue
        key = stripped.split("=", 1)[0].strip()
        if key in updates:
            out_lines.append(f"{key:<24} = {updates[key]}")
            seen.add(key)
        else:
            out_lines.append(line)
    for key, value in updates.items():
        if key not in seen:
            out_lines.append(f"{key:<24} = {value}")
    return "\n".join(out_lines) + "\n"


def safe_nsteps(time_ns: float, dt_ps: float) -> int:
    if dt_ps <= 0:
        raise ValueError("Time step must be > 0 ps.")
    if time_ns <= 0:
        raise ValueError("Time must be > 0 ns.")
    return int(round((time_ns * 1000.0) / dt_ps))


def extract_mdrun_performance(text: str) -> tuple[float | None, float | None]:
    match = re.search(r"Performance:\s+([0-9.+-Ee]+)\s+ns/day,\s+([0-9.+-Ee]+)\s+hours/ns", text)
    if not match:
        return None, None
    return float(match.group(1)), float(match.group(2))


def format_elapsed(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f} s"
    minutes, sec = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)} min {sec:.0f} s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)} h {int(minutes)} min {sec:.0f} s"


def _run_cmd(
    cmd: list[str],
    cwd: Path,
    env: dict[str, str],
    stdout_parts: list[str],
    stderr_parts: list[str],
    command_log: list[str],
) -> tuple[subprocess.CompletedProcess[str], float]:
    command_log.append(" ".join(cmd))
    stdout_parts.append("+ " + " ".join(cmd) + "\n")
    start = time.time()
    result = subprocess.run(cmd, cwd=cwd, env=env, text=True, capture_output=True, check=False)
    elapsed = time.time() - start
    if result.stdout:
        stdout_parts.append(result.stdout[-4000:] + "\n")
    if result.stderr:
        stderr_parts.append(result.stderr[-4000:] + "\n")
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
    return result, elapsed


def _path_with_tools(extra_tool_dirs: list[Path] | None) -> str:
    tool_dirs = [Path(sys.executable).resolve().parent]
    if extra_tool_dirs:
        tool_dirs.extend(extra_tool_dirs)
    return os.pathsep.join(str(path) for path in tool_dirs) + os.pathsep + os.environ.get("PATH", "")


def _pick_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def _mdp_path_for(mdp_dir: Path, stage_name: str, suffix: str) -> Path:
    candidate = mdp_dir / f"{stage_name}{suffix}.mdp"
    if not candidate.exists():
        if stage_name == "deposition":
            raise FileNotFoundError(
                f"Missing MDP file: {candidate}. Deposition is not generated in Adsorption mode."
            )
        raise FileNotFoundError(f"Missing MDP file: {candidate}")
    return candidate


def _safe_tag(value: str) -> str:
    tag = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return tag or DEFAULT_OUTPUT_TAG
