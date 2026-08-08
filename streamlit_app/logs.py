from __future__ import annotations


def summarize_log(stdout: str, stderr: str, returncode: int | None) -> list[tuple[str, str]]:
    text = (stdout + "\n" + stderr).strip()
    rows: list[tuple[str, str]] = []

    if returncode is not None:
        rows.append(("Status", "Completed" if returncode == 0 else "Failed"))

    for marker, label in [
        ("Downloaded", "Input"),
        ("Cleaned PDB", "Cleaning"),
        ("Protein mode", "Mode"),
        ("Surface", "Surface"),
        ("Topology", "Topology"),
        ("Generated", "Output"),
    ]:
        match = next((line.strip() for line in text.splitlines() if marker in line), "")
        if match:
            rows.append((label, match))

    if not rows and text:
        rows.append(("Log", text.splitlines()[-1]))

    return rows
