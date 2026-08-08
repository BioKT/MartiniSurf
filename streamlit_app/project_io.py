from __future__ import annotations

import json
from pathlib import Path
from typing import Any


PROJECT_KEYS = [
    "project_name",
    "preset_name",
    "remote_id",
    "moltype",
    "ff",
    "merge_text",
    "dssp",
    "go",
    "elastic",
    "position_restraints",
    "surface_mode",
    "surface_geometry",
    "surface_beads",
    "lx",
    "ly",
    "dx",
    "charge",
    "surface_layers",
    "surface_stacking",
    "surface_dist_z",
    "graphite_layers",
    "cnt_numrings",
    "cnt_ringsize",
    "orientation_mode",
    "anchors_text",
    "dist",
    "balance_low_z",
    "balance_low_z_fraction",
    "histag",
    "linker_path_text",
    "linker_groups_text",
    "linker_prot_dist",
    "linker_surf_dist",
    "invert_linker",
    "surface_linkers",
    "solvate",
    "ionize",
    "salt_conc",
    "water_mix",
    "substrate_count",
    "martinize_extra",
]


def export_state(session_state: Any) -> str:
    payload = {key: session_state[key] for key in PROJECT_KEYS if key in session_state}
    payload["schema"] = "martinisurf-studio-v1"
    return json.dumps(payload, indent=2, sort_keys=True)


def import_state(uploaded_file: Any) -> dict[str, Any]:
    data = json.loads(uploaded_file.getvalue().decode("utf-8"))
    if data.get("schema") != "martinisurf-studio-v1":
        raise ValueError("Unsupported MartiniSurf Studio project file.")
    return {key: data[key] for key in PROJECT_KEYS if key in data}
