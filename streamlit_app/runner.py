from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class RunResult:
    returncode: int
    stdout: str
    stderr: str
    command: list[str]


def run_martinisurf(args: list[str], workdir: Path, repo_root: Path) -> RunResult:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root) + os.pathsep + env.get("PYTHONPATH", "")
    tool_paths = [
        str(Path(sys.executable).resolve().parent),
        str(repo_root / ".venv" / "bin"),
    ]
    env["PATH"] = os.pathsep.join(tool_paths + [env.get("PATH", "")])
    command = [sys.executable, "-m", "martinisurf", *args]
    process = subprocess.run(
        command,
        cwd=workdir,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    return RunResult(
        returncode=process.returncode,
        stdout=process.stdout,
        stderr=process.stderr,
        command=command,
    )
