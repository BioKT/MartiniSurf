import os
import shutil
import subprocess


def _is_python2(exe: str) -> bool:
    try:
        res = subprocess.run(
            [exe, "-c", "import sys; print(sys.version_info[0])"],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return False
    return res.returncode == 0 and res.stdout.strip() == "2"


def find_python2(candidates=None):
    env_exe = os.environ.get("MARTINISURF_PYTHON2", "").strip()
    default_candidates = ["python2.7", "python2", "/usr/bin/python2.7", "/usr/bin/python2"]
    ordered_candidates = ([env_exe] if env_exe else []) + list(candidates or default_candidates)

    incompatible = []
    for exe in ordered_candidates:
        resolved = shutil.which(exe) or exe
        if not shutil.which(resolved) and not os.path.exists(resolved):
            continue
        if _is_python2(resolved):
            return resolved
        incompatible.append(exe)

    incompatible_msg = ""
    if incompatible:
        incompatible_msg = (
            "\nDetected non-Python2 executables for candidates: "
            + ", ".join(incompatible)
            + "."
        )
    raise RuntimeError(
        "❌ DNA mode requires python2.7, but no python2 interpreter was found.\n"
        "You can force a path with MARTINISURF_PYTHON2=/path/to/python2.7\n"
        + incompatible_msg
        + "\n"
        "Install it with:\n"
        "    conda install -n martinisurf-dna python=2.7\n"
    )
