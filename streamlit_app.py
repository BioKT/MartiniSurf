from __future__ import annotations

import base64
import html
import shutil
import sys
import tempfile
from dataclasses import fields
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from streamlit_app.archive import make_zip
from streamlit_app.command_builder import BuildConfig, build_args, shell_command
from streamlit_app.logs import summarize_log
from streamlit_app.linker_generator import (
    generate_linker_from_smiles,
    parse_linker_gro,
    parse_martini_mapper_report,
    safe_molecule_name,
    smiles_to_svg,
)
from streamlit_app.molecular_viewer import (
    file_caption,
    find_viewable_structures,
    render_linker_mapping,
    render_martini_step4_viewer,
    render_molecule,
    render_remote_molecule,
    render_structure_preview,
)
from streamlit_app.project_io import export_state, import_state
from streamlit_app.runner import run_martinisurf
from streamlit_app.structure import StructureSummary, summarize_structure
from streamlit_app.theme import apply_theme
from streamlit_app.validators import check_external_tools, validate_config


REPO_ROOT = Path(__file__).resolve().parent
RUNS_DIR = REPO_ROOT / "streamlit_runs"
APP_DEPLOY_REVISION = "2026-08-10-local-full"
TOOL_DIRS = [Path(sys.executable).resolve().parent, REPO_ROOT / ".venv" / "bin"]

STEPS = ["Structure", "Model", "Surface", "Orientation", "Environment", "Review & Build"]

CHAIN_LABELS = list("ABCDEFGHI")
PROTEIN_COLAB_DEFAULTS_VERSION = 3
SURFACE_DEFAULTS_VERSION = 5
PERSISTED_STATE_KEYS = {
    "project_name",
    "active_step",
    "remote_id",
    "moltype",
    "ff",
    "merge_text",
    "merge_chain_count",
    "dssp",
    "go",
    "go_eps",
    "elastic",
    "maxwarn",
    "position_restraints",
    "surface_workflow",
    "carbon_surface_kind",
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
    "anchor_group_count",
    "dist",
    "balance_low_z",
    "balance_low_z_fraction",
    "histag",
    "linker_path_text",
    "linker_groups_text",
    "linker_group_count",
    "linker_prot_dist",
    "linker_surf_dist",
    "invert_linker",
    "surface_linkers",
    "linker_generator_smiles",
    "linker_generator_molname",
    "generated_linker_path",
    "generated_linker_itp_path",
    "generated_linker_beads",
    "generated_linker_mapping",
    "linker_surface_bead",
    "linker_protein_bead",
    "linker_generator_log",
    "linker_generator_message",
    "solvate",
    "ionize",
    "salt_conc",
    "water_mix",
    "enable_water_mix",
    "sw_water_percent",
    "tw_water_percent",
    "water_mix_seed",
    "substrate_count",
    "martinize_extra",
    "_structure_path",
    "_surface_path",
    "_linker_path",
    "_linker_itp_path",
    "_substrate_path",
    "_substrate_itp_path",
}

HEXAGONAL_SURFACE_DEFAULTS = {
    "surface_mode": "4-1",
    "surface_geometry": "planar",
    "surface_beads": "P4 P4",
    "lx": 6.0,
    "ly": 6.0,
    "dx": 0.5,
    "charge": 0.0,
    "surface_layers": 2,
    "surface_stacking": "hcp",
}

CARBON_SURFACE_DEFAULTS = {
    "Graphene": {
        "surface_mode": "graphene",
        "surface_geometry": "planar",
        "surface_beads": "C1",
        "lx": 6.0,
        "ly": 6.0,
        "graphite_layers": 5,
        "cnt_numrings": 24,
        "cnt_ringsize": 9,
    },
    "Graphite": {
        "surface_mode": "graphite",
        "surface_geometry": "planar",
        "surface_beads": "C1",
        "lx": 6.0,
        "ly": 6.0,
        "graphite_layers": 5,
        "cnt_numrings": 24,
        "cnt_ringsize": 9,
    },
    "Nanotubes": {
        "surface_mode": "cnt-martini3",
        "surface_geometry": "3d",
        "surface_beads": "C1",
        "lx": 6.0,
        "ly": 6.0,
        "graphite_layers": 5,
        "cnt_numrings": 24,
        "cnt_ringsize": 9,
    },
}


def _init_state() -> None:
    defaults = {
        "project_name": "MartiniSurf Protein - 1UBQ",
        "active_step": "Structure",
        "remote_id": "1UBQ",
        "moltype": "Protein",
        "ff": "martini3001",
        "merge_text": "A",
        "merge_chain_count": 1,
        "dssp": True,
        "go": True,
        "go_eps": 9.415,
        "elastic": False,
        "maxwarn": 1,
        "position_restraints": "backbone",
        "surface_workflow": "Hexagonal Lattice",
        "carbon_surface_kind": "Graphene",
        "surface_mode": "4-1",
        "surface_geometry": "planar",
        "surface_beads": "P4 P4",
        "lx": 6.0,
        "ly": 6.0,
        "dx": 0.5,
        "charge": 0.0,
        "surface_layers": 2,
        "surface_stacking": "hcp",
        "surface_dist_z": 0.0,
        "graphite_layers": 5,
        "cnt_numrings": 24,
        "cnt_ringsize": 9,
        "orientation_mode": "Anchor",
        "anchors_text": "A 73",
        "anchor_group_count": 1,
        "anchor_chain_0": "A",
        "anchor_residues_0": "73",
        "anchor_chain_1": "A",
        "anchor_residues_1": "",
        "dist": 1.0,
        "balance_low_z": False,
        "balance_low_z_fraction": 0.2,
        "histag": False,
        "linker_path_text": "martinisurf/examples/protein/03_linker_surface_decoration/inputs/EPOXY.gro",
        "linker_groups_text": "A 73",
        "linker_group_count": 1,
        "linker_group_chain_0": "A",
        "linker_group_residues_0": "73",
        "linker_prot_dist": 0.0,
        "linker_surf_dist": 0.0,
        "invert_linker": False,
        "surface_linkers": 0,
        "linker_generator_smiles": "",
        "linker_generator_molname": "LINKER",
        "generated_linker_path": "",
        "generated_linker_itp_path": "",
        "generated_linker_beads": [],
        "generated_linker_mapping": [],
        "linker_surface_bead": "",
        "linker_protein_bead": "",
        "linker_generator_log": "",
        "linker_generator_message": "",
        "solvate": True,
        "ionize": True,
        "salt_conc": 0.15,
        "water_mix": "",
        "enable_water_mix": False,
        "sw_water_percent": 0,
        "tw_water_percent": 0,
        "water_mix_seed": 42,
        "substrate_count": 0,
        "martinize_extra": "",
    }
    for key, value in defaults.items():
        st.session_state.setdefault(key, value)
    if not str(st.session_state.get("moltype", "")).strip():
        st.session_state.moltype = "Protein"
    if float(st.session_state.get("go_eps") or 0) <= 0:
        st.session_state.go_eps = 9.415
    if st.session_state.get("_protein_colab_defaults_version") != PROTEIN_COLAB_DEFAULTS_VERSION:
        for key in [
            "project_name",
            "remote_id",
            "moltype",
            "merge_text",
            "merge_chain_count",
            "go",
            "go_eps",
            "maxwarn",
            "orientation_mode",
            "anchors_text",
            "anchor_group_count",
            "anchor_chain_0",
            "anchor_residues_0",
            "anchor_chain_1",
            "anchor_residues_1",
            "dist",
            "linker_groups_text",
            "linker_group_count",
            "surface_workflow",
            "carbon_surface_kind",
            "surface_linkers",
            "linker_generator_smiles",
            "linker_generator_molname",
            "generated_linker_path",
            "generated_linker_itp_path",
            "generated_linker_beads",
            "generated_linker_mapping",
            "linker_surface_bead",
            "linker_protein_bead",
            "linker_generator_log",
            "linker_generator_message",
            "solvate",
            "ionize",
            "salt_conc",
            "enable_water_mix",
            "sw_water_percent",
            "tw_water_percent",
            "water_mix_seed",
            "substrate_count",
        ]:
            st.session_state[key] = defaults[key]
        st.session_state["_protein_colab_defaults_version"] = PROTEIN_COLAB_DEFAULTS_VERSION
    if st.session_state.get("_surface_defaults_version") != SURFACE_DEFAULTS_VERSION:
        _apply_surface_workflow_defaults()
        st.session_state["_surface_defaults_version"] = SURFACE_DEFAULTS_VERSION


def _lines(value: str) -> list[str]:
    return [line.strip() for line in value.splitlines() if line.strip()]


def _choice(value: object, default: str) -> str:
    if isinstance(value, str) and value:
        return value
    if isinstance(value, (list, tuple)) and value:
        return str(value[0])
    return default


def _versioned_key(name: str) -> str:
    return f"{name}_ui_colab_v{PROTEIN_COLAB_DEFAULTS_VERSION}_{SURFACE_DEFAULTS_VERSION}"


def _commit_transient_widget_state() -> None:
    aliases = {
        _versioned_key("moltype"): "moltype",
        _versioned_key("merge_chain_count"): "merge_chain_count",
        _versioned_key("dssp"): "dssp",
        _versioned_key("go"): "go",
        _versioned_key("go_eps"): "go_eps",
        _versioned_key("position_restraints"): "position_restraints",
        _versioned_key("surface_workflow"): "surface_workflow",
        _versioned_key("surface_mode_hex"): "surface_mode",
        _versioned_key("dx_hex"): "dx",
        _versioned_key("surface_beads_hex"): "surface_beads",
        _versioned_key("lx_hex"): "lx",
        _versioned_key("ly_hex"): "ly",
        _versioned_key("surface_layers_hex"): "surface_layers",
        _versioned_key("surface_stacking_hex"): "surface_stacking",
        _versioned_key("charge_hex"): "charge",
        _versioned_key("carbon_surface_kind"): "carbon_surface_kind",
        _versioned_key("surface_linkers"): "surface_linkers",
        _versioned_key("surface_geometry"): "surface_geometry",
        _versioned_key("solvate"): "solvate",
        _versioned_key("ionize"): "ionize",
        _versioned_key("salt_conc"): "salt_conc",
        _versioned_key("enable_water_mix"): "enable_water_mix",
        _versioned_key("sw_water_percent"): "sw_water_percent",
        _versioned_key("tw_water_percent"): "tw_water_percent",
        _versioned_key("water_mix_seed"): "water_mix_seed",
        _versioned_key("substrate_count"): "substrate_count",
    }
    for kind in ["Graphene", "Graphite", "Nanotubes"]:
        aliases[_versioned_key(f"lx_carbon_{kind}")] = "lx"
        aliases[_versioned_key(f"ly_carbon_{kind}")] = "ly"
    aliases[_versioned_key("graphite_layers")] = "graphite_layers"
    aliases[_versioned_key("cnt_numrings")] = "cnt_numrings"
    aliases[_versioned_key("cnt_ringsize")] = "cnt_ringsize"
    for widget_key, state_key in aliases.items():
        if widget_key in st.session_state:
            st.session_state[state_key] = st.session_state[widget_key]
    if (
        st.session_state.get("_anchor_rows_loaded_from_text") == st.session_state.get("anchors_text")
        and any(key.startswith("anchor_chain_") or key.startswith("anchor_residues_") for key in st.session_state)
    ):
        _sync_anchor_text()
    if (
        st.session_state.get("_linker_group_rows_loaded_from_text") == st.session_state.get("linker_groups_text")
        and any(key.startswith("linker_group_chain_") or key.startswith("linker_group_residues_") for key in st.session_state)
    ):
        _sync_linker_group_text()


def _preserve_app_state() -> None:
    keys = set(PERSISTED_STATE_KEYS)
    for key in list(st.session_state.keys()):
        if key.startswith("anchor_chain_") or key.startswith("anchor_residues_"):
            keys.add(key)
        if key.startswith("linker_group_chain_") or key.startswith("linker_group_residues_"):
            keys.add(key)
    for key in keys:
        if key in st.session_state:
            st.session_state[key] = st.session_state[key]


def _image_data_uri(path: Path) -> str:
    if not path.exists():
        return ""
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _save_upload(uploaded_file, target_dir: Path, state_key: str | None = None) -> Path | None:
    if uploaded_file is None:
        remembered = st.session_state.get(state_key) if state_key else None
        if remembered:
            path = Path(str(remembered))
            if path.exists():
                return path
        return None
    path = target_dir / uploaded_file.name
    path.write_bytes(uploaded_file.getbuffer())
    if state_key:
        st.session_state[state_key] = str(path)
    return path


def _remember_path(key: str, path: Path | None) -> None:
    if path:
        st.session_state[key] = str(path)


def _chain_count(summary: StructureSummary) -> int | None:
    if summary.chains:
        return max(1, min(len(summary.chains), len(CHAIN_LABELS)))
    if st.session_state.remote_id.strip().upper() == "1UBQ":
        return 1
    return None


def _merge_group_from_count(count: int) -> str:
    return ",".join(CHAIN_LABELS[: max(1, min(count, len(CHAIN_LABELS)))])


def _sync_anchor_text() -> None:
    count = int(st.session_state.get("anchor_group_count", 1))
    groups = []
    for index in range(count):
        chain = str(st.session_state.get(f"anchor_chain_{index}", "A")).strip() or "A"
        residues = str(st.session_state.get(f"anchor_residues_{index}", "")).strip()
        if residues:
            groups.append(f"{chain} {residues}")
    st.session_state.anchors_text = "\n".join(groups)
    st.session_state["_anchor_rows_loaded_from_text"] = st.session_state.anchors_text


def _remove_anchor_row(remove_index: int) -> None:
    count = int(st.session_state.get("anchor_group_count", 1))
    if count <= 1:
        return
    for index in range(remove_index, count - 1):
        st.session_state[f"anchor_chain_{index}"] = st.session_state.get(f"anchor_chain_{index + 1}", "A")
        st.session_state[f"anchor_residues_{index}"] = st.session_state.get(f"anchor_residues_{index + 1}", "")
    st.session_state.pop(f"anchor_chain_{count - 1}", None)
    st.session_state.pop(f"anchor_residues_{count - 1}", None)
    st.session_state.anchor_group_count = count - 1
    _sync_anchor_text()


def _add_anchor_row() -> None:
    count = int(st.session_state.get("anchor_group_count", 1))
    if count >= 20:
        return
    st.session_state.anchor_group_count = count + 1
    st.session_state[f"anchor_chain_{count}"] = "A"
    st.session_state[f"anchor_residues_{count}"] = ""
    _sync_anchor_text()


def _load_anchor_rows_from_text() -> None:
    raw_text = st.session_state.get("anchors_text", "")
    expected_count = int(st.session_state.get("anchor_group_count", 1))
    keys_present = all(
        f"anchor_chain_{index}" in st.session_state and f"anchor_residues_{index}" in st.session_state
        for index in range(expected_count)
    )
    if keys_present and st.session_state.get("_anchor_rows_loaded_from_text") == raw_text:
        return
    groups = _lines(raw_text) or ["A 73"]
    st.session_state.anchor_group_count = len(groups)
    for index, group in enumerate(groups):
        parts = group.split(maxsplit=1)
        st.session_state[f"anchor_chain_{index}"] = parts[0] if parts else "A"
        st.session_state[f"anchor_residues_{index}"] = parts[1] if len(parts) > 1 else ""
    st.session_state["_anchor_rows_loaded_from_text"] = raw_text


def _sync_linker_group_text() -> None:
    count = int(st.session_state.get("linker_group_count", 1))
    rows = []
    for index in range(count):
        chain = str(st.session_state.get(f"linker_group_chain_{index}", "")).strip().upper()
        residues = str(st.session_state.get(f"linker_group_residues_{index}", "")).strip()
        if chain and residues:
            rows.append(f"{chain} {residues}")
    st.session_state.linker_groups_text = "\n".join(rows)
    st.session_state["_linker_group_rows_loaded_from_text"] = st.session_state.linker_groups_text


def _remove_linker_group_row(index_to_remove: int) -> None:
    count = int(st.session_state.get("linker_group_count", 1))
    if count <= 1:
        st.session_state[f"linker_group_chain_{index_to_remove}"] = "A"
        st.session_state[f"linker_group_residues_{index_to_remove}"] = ""
        _sync_linker_group_text()
        return
    for index in range(index_to_remove, count - 1):
        st.session_state[f"linker_group_chain_{index}"] = st.session_state.get(f"linker_group_chain_{index + 1}", "A")
        st.session_state[f"linker_group_residues_{index}"] = st.session_state.get(f"linker_group_residues_{index + 1}", "")
    st.session_state.pop(f"linker_group_chain_{count - 1}", None)
    st.session_state.pop(f"linker_group_residues_{count - 1}", None)
    st.session_state.linker_group_count = count - 1
    _sync_linker_group_text()


def _add_linker_group_row() -> None:
    count = int(st.session_state.get("linker_group_count", 1))
    if count >= 20:
        return
    st.session_state.linker_group_count = count + 1
    st.session_state[f"linker_group_chain_{count}"] = "A"
    st.session_state[f"linker_group_residues_{count}"] = ""
    _sync_linker_group_text()


def _load_linker_group_rows_from_text() -> None:
    raw_text = st.session_state.get("linker_groups_text", "")
    expected_count = int(st.session_state.get("linker_group_count", 1))
    keys_present = all(
        f"linker_group_chain_{index}" in st.session_state and f"linker_group_residues_{index}" in st.session_state
        for index in range(expected_count)
    )
    if keys_present and st.session_state.get("_linker_group_rows_loaded_from_text") == raw_text:
        return
    groups = _lines(raw_text) or ["A 73"]
    st.session_state.linker_group_count = len(groups)
    for index, group in enumerate(groups):
        parts = group.split(maxsplit=1)
        st.session_state[f"linker_group_chain_{index}"] = parts[0] if parts else "A"
        st.session_state[f"linker_group_residues_{index}"] = parts[1] if len(parts) > 1 else ""
    st.session_state["_linker_group_rows_loaded_from_text"] = raw_text


def _apply_values(values: dict[str, object]) -> None:
    for key, value in values.items():
        st.session_state[key] = value


def _apply_carbon_defaults() -> None:
    kind = st.session_state.get("carbon_surface_kind", "Graphene")
    if kind == "Nanotube":
        kind = "Nanotubes"
        st.session_state.carbon_surface_kind = kind
    _apply_values(CARBON_SURFACE_DEFAULTS.get(kind, CARBON_SURFACE_DEFAULTS["Graphene"]))


def _apply_surface_workflow_defaults() -> None:
    workflow = st.session_state.get("surface_workflow", "Hexagonal Lattice")
    if workflow == "Hexagonal Lattice":
        _apply_values(HEXAGONAL_SURFACE_DEFAULTS)
    elif workflow == "Carbon Nanomaterial":
        _apply_carbon_defaults()


def _resolve_optional_path(value: str, workdir: Path) -> Path | None:
    clean = value.strip()
    if not clean:
        return None
    path = Path(clean)
    if not path.is_absolute():
        path = REPO_ROOT / path
    if path.exists():
        return path
    candidate = workdir / clean
    return candidate if candidate.exists() else path


def _bead_dicts(beads: list[object]) -> list[dict[str, object]]:
    payload = []
    for bead in beads:
        payload.append(
            {
                "index": int(getattr(bead, "index")),
                "name": str(getattr(bead, "name")),
                "resid": int(getattr(bead, "resid")),
                "resname": str(getattr(bead, "resname")),
                "x": float(getattr(bead, "x")),
                "y": float(getattr(bead, "y")),
                "z": float(getattr(bead, "z")),
                "label": str(getattr(bead, "label")),
            }
        )
    return payload


def _mapping_dicts(rows: list[object]) -> list[dict[str, object]]:
    payload = []
    for row in rows:
        payload.append(
            {
                "Bead": str(getattr(row, "bead_label")),
                "Martini type": str(getattr(row, "martini_type")),
                "Atom indices": ", ".join(str(index) for index in getattr(row, "atom_indices_1based")),
            }
        )
    return payload


def _linker_bead_options() -> list[str]:
    beads = st.session_state.get("generated_linker_beads") or []
    return [str(bead.get("label") or f"{bead.get('index')}: {bead.get('name')}") for bead in beads if isinstance(bead, dict)]


def _refresh_generated_linker_from_path() -> None:
    path_text = st.session_state.get("generated_linker_path") or st.session_state.get("linker_path_text", "")
    if not path_text:
        return
    path = _resolve_optional_path(str(path_text), Path(st.session_state.run_root))
    if not path or not path.exists() or path.suffix.lower() != ".gro":
        return
    try:
        parsed_beads = parse_linker_gro(path)
        beads = _bead_dicts(parsed_beads)
    except OSError:
        return
    if beads:
        st.session_state.generated_linker_path = str(path)
        st.session_state.generated_linker_beads = beads
        report_path = path.with_suffix(".txt")
        mapping = _mapping_dicts(parse_martini_mapper_report(report_path, parsed_beads))
        if mapping:
            st.session_state.generated_linker_mapping = mapping
        options = _linker_bead_options()
        if options:
            st.session_state.setdefault("linker_protein_bead", options[0])
            st.session_state.setdefault("linker_surface_bead", options[-1])


def _sync_linker_orientation_from_selected_beads() -> None:
    options = _linker_bead_options()
    if len(options) < 2:
        return
    st.session_state.invert_linker = False
    first, last = options[0], options[-1]
    surface_bead = st.session_state.get("linker_surface_bead", last)
    protein_bead = st.session_state.get("linker_protein_bead", first)
    if protein_bead == first and surface_bead == last:
        st.session_state.invert_linker = False
    elif protein_bead == last and surface_bead == first:
        st.session_state.invert_linker = True


def _render_linker_generator() -> None:
    with st.expander("Linker Generator", expanded=False):
        st.caption("Generate a Martini 3 linker from SMILES and select the terminal beads used by MartiniSurf.")
        smiles_col, name_col = st.columns([2.2, 0.8], gap="large")
        smiles_col.text_input("SMILES", key="linker_generator_smiles", placeholder="C1CO1")
        name_col.text_input("Linker name", key="linker_generator_molname")
        if st.button("Generate linker", type="primary", width="stretch"):
            molname = safe_molecule_name(st.session_state.linker_generator_molname)
            output_dir = Path(st.session_state.run_root) / "inputs" / "linker_generator" / molname
            result = generate_linker_from_smiles(st.session_state.linker_generator_smiles, molname, output_dir)
            st.session_state.linker_generator_message = result.message
            st.session_state.linker_generator_log = result.log
            if result.ok and result.gro_path and result.itp_path:
                st.session_state.generated_linker_path = str(result.gro_path)
                st.session_state.generated_linker_itp_path = str(result.itp_path)
                st.session_state.generated_linker_beads = _bead_dicts(result.beads)
                st.session_state.generated_linker_mapping = _mapping_dicts(
                    parse_martini_mapper_report(result.gro_path.with_suffix(".txt"), result.beads)
                )
                st.session_state.linker_path_text = str(result.gro_path)
                st.session_state._linker_path = str(result.gro_path)
                st.session_state._linker_itp_path = str(result.itp_path)
                options = _linker_bead_options()
                if options:
                    st.session_state.linker_protein_bead = options[0]
                    st.session_state.linker_surface_bead = options[-1]
                _sync_linker_orientation_from_selected_beads()
                st.success(f"Linker generated with {result.generator}.")
            else:
                st.error(result.message)

        _refresh_generated_linker_from_path()
        options = _linker_bead_options()
        generated_path = Path(st.session_state.generated_linker_path) if st.session_state.get("generated_linker_path") else None
        if generated_path and generated_path.exists() and options:
            preview_col, beads_col = st.columns([1.35, 1], gap="large")
            with preview_col:
                svg = smiles_to_svg(st.session_state.linker_generator_smiles)
                if svg:
                    components.html(
                        f"<div style='background:#f7fbff;border-radius:16px;padding:12px;margin-bottom:12px'>{svg}</div>",
                        height=350,
                    )
                    st.caption("Atom indices in the 2D structure and mapping table start at 1.")
                render_linker_mapping(generated_path, st.session_state.generated_linker_beads, height=360)
            with beads_col:
                mapping_rows = st.session_state.get("generated_linker_mapping") or []
                if mapping_rows:
                    st.markdown("##### Bead mapping")
                    st.dataframe(mapping_rows, hide_index=True, width="stretch")
                else:
                    st.dataframe(
                        [
                            {"Bead": bead["label"], "Residue": bead["resname"], "x": bead["x"], "y": bead["y"], "z": bead["z"]}
                            for bead in st.session_state.generated_linker_beads
                        ],
                        hide_index=True,
                        width="stretch",
                    )
                    st.caption("This linker has no Martini Mapper chemistry report beside the `.gro` file.")
                if st.session_state.linker_surface_bead not in options:
                    st.session_state.linker_surface_bead = options[-1]
                if st.session_state.linker_protein_bead not in options:
                    st.session_state.linker_protein_bead = options[0]
                st.selectbox("Bead oriented toward surface", options, key="linker_surface_bead")
                st.selectbox("Bead connected to protein", options, key="linker_protein_bead")
                _sync_linker_orientation_from_selected_beads()
                st.caption(f"Linker file: `{generated_path.name}`")
        elif st.session_state.get("linker_generator_message"):
            st.info(st.session_state.linker_generator_message)
        if st.session_state.get("linker_generator_log"):
            with st.expander("Generator log"):
                st.code(st.session_state.linker_generator_log, language="text")


def _step_state(step: str, errors: list[str], tool_warnings: list[str], has_output: bool) -> str:
    active_index = STEPS.index(st.session_state.active_step)
    index = STEPS.index(step)
    if has_output:
        return "completed"
    if step == "Review & Build" and errors:
        return "error"
    if step == "Review & Build" and tool_warnings:
        return "warning"
    if index < active_index:
        return "completed"
    if index == active_index:
        return "active"
    return "pending"


def _render_sidebar(errors: list[str], tool_warnings: list[str], has_output: bool) -> None:
    with st.sidebar:
        logo = REPO_ROOT / "logo.png"
        st.image(str(logo), width="stretch")
        selected = st.radio("Workflow step", STEPS, key="active_step", label_visibility="collapsed")
        completed = sum(1 for step in STEPS if _step_state(step, errors, tool_warnings, has_output) == "completed")
        st.progress(completed / len(STEPS), text=f"{completed} of {len(STEPS)}")


def _render_topbar(config: BuildConfig, errors: list[str], tool_warnings: list[str], has_output: bool) -> None:
    if st.session_state.get("building"):
        status = "Building"
    elif st.session_state.get("last_returncode") not in (None, 0):
        status = "Error"
    elif has_output:
        status = "Completed"
    elif errors:
        status = "Not configured"
    elif tool_warnings:
        status = "Warning"
    else:
        status = "Ready"

    st.markdown(
        f"""
        <!-- MartiniSurf deploy revision: {APP_DEPLOY_REVISION} -->
        <div class="ms-topbar">
          <div class="ms-top-brand">
            <div>
              <div class="ms-top-title">MartiniSurf</div>
              <div class="ms-top-subtitle">{html.escape(st.session_state.project_name)}</div>
            </div>
          </div>
          <div class="ms-ready {status.lower().replace(' ', '-')}"><span></span>{html.escape(status)}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    tools = st.columns(3)
    tools[0].link_button("Documentation", "https://biokt.github.io/MartiniSurf/", width="stretch")
    with tools[1].popover("Import config", use_container_width=True):
        import_file = st.file_uploader("Project JSON", type=["json"])
        if import_file is not None:
            try:
                for key, value in import_state(import_file).items():
                    st.session_state[key] = value
                st.success("Project configuration imported.")
                st.rerun()
            except Exception as exc:
                st.error(f"Could not import project: {exc}")
    tools[2].download_button(
        "Export project",
        data=export_state(st.session_state),
        file_name="martinisurf-project.json",
        mime="application/json",
        width="stretch",
    )


def _render_structure_step(upload_dir: Path) -> tuple[Path | None, Path | None, Path | None, Path | None, Path | None, Path | None, StructureSummary]:
    left, right = st.columns([0.9, 1.1], gap="large")
    with left:
        st.markdown('<div class="ms-panel-title">Structure input</div>', unsafe_allow_html=True)
        st.text_input("Project name", key="project_name")
        uploaded_structure = st.file_uploader("Protein structure", type=["pdb", "cif", "mmcif"], key="structure_upload")
        st.text_input("PDB ID or UniProt ID", key="remote_id")
        with st.expander("Optional auxiliary files"):
            uploaded_surface = st.file_uploader("Existing surface .gro", type=["gro"], key="surface_upload")
            uploaded_linker = st.file_uploader("Linker .gro", type=["gro"], key="linker_upload")
            uploaded_linker_itp = st.file_uploader("Linker .itp", type=["itp"], key="linker_itp_upload")
            uploaded_substrate = st.file_uploader("Substrate .gro", type=["gro"], key="substrate_upload")
            uploaded_substrate_itp = st.file_uploader("Substrate .itp", type=["itp"], key="substrate_itp_upload")

    structure_path = _save_upload(uploaded_structure, upload_dir, "_structure_path")
    surface_path = _save_upload(uploaded_surface, upload_dir, "_surface_path")
    linker_path = _save_upload(uploaded_linker, upload_dir, "_linker_path")
    linker_itp_path = _save_upload(uploaded_linker_itp, upload_dir, "_linker_itp_path")
    if linker_path and linker_itp_path:
        matching_itp = linker_path.with_suffix(".itp")
        if linker_itp_path != matching_itp:
            shutil.copyfile(linker_itp_path, matching_itp)
            linker_itp_path = matching_itp
    substrate_path = _save_upload(uploaded_substrate, upload_dir, "_substrate_path")
    substrate_itp_path = _save_upload(uploaded_substrate_itp, upload_dir, "_substrate_itp_path")
    summary = summarize_structure(structure_path, st.session_state.remote_id, upload_dir / "preview")

    with right:
        _render_structure_summary(summary, structure_path)

    return structure_path, surface_path, linker_path, linker_itp_path, substrate_path, substrate_itp_path, summary


def _render_structure_summary(summary: StructureSummary, structure_path: Path | None) -> None:
    st.markdown('<div class="ms-panel-title">Structure preview</div>', unsafe_allow_html=True)
    preview_metrics = st.columns(2)
    preview_metrics[0].metric("Source", summary.source)
    preview_metrics[1].metric("Chains", len(summary.chains) if summary.chains else "-")

    if structure_path and structure_path.suffix.lower() in {".pdb", ".gro"}:
        render_structure_preview(structure_path, height=430)
    elif summary.format == "remote":
        render_remote_molecule(summary.source, height=430)
    else:
        st.markdown(
            """
            <div class="ms-empty-preview">
              <strong>Preview pending</strong>
              <span>Local PDB/GRO structures can be inspected here before build. Remote IDs are fetched by MartiniSurf during execution.</span>
            </div>
            """,
            unsafe_allow_html=True,
        )

    amino_acid_rows = [
        {"Chain": item.chain, "Amino acids": item.residues}
        for item in summary.chains
        if item.residues
    ]
    if amino_acid_rows:
        st.dataframe(amino_acid_rows, hide_index=True, width="stretch")

    if summary.ligands:
        st.info("Ligands/cofactors detected: " + ", ".join(summary.ligands[:12]))
    for note in summary.notes:
        st.warning(note)


def _render_model_step(summary: StructureSummary) -> None:
    st.markdown('<div class="ms-panel-title">Model decisions</div>', unsafe_allow_html=True)
    available_count = _chain_count(summary)
    max_count = available_count or len(CHAIN_LABELS)
    if int(st.session_state.merge_chain_count) > max_count:
        st.session_state.merge_chain_count = max_count
    st.session_state.moltype = st.text_input(
        "Molecule name",
        value=st.session_state.moltype,
        key=_versioned_key("moltype"),
    )
    chain_options = list(range(1, max_count + 1))
    merge_count = st.selectbox(
        "Merge chains",
        chain_options,
        index=chain_options.index(int(st.session_state.merge_chain_count)),
        key=_versioned_key("merge_chain_count"),
        format_func=lambda value: f"{value} chain" if value == 1 else f"{value} chains",
    )
    st.session_state.merge_chain_count = int(merge_count)
    st.session_state.merge_text = _merge_group_from_count(int(st.session_state.merge_chain_count))
    st.markdown('<div class="ms-model-toggles"></div>', unsafe_allow_html=True)
    c1, c2 = st.columns(2)
    st.session_state.dssp = c1.toggle(
        "DSSP",
        value=bool(st.session_state.dssp),
        key=_versioned_key("dssp"),
        help="Passes DSSP handling to martinize2.",
    )
    st.session_state.go = c2.toggle(
        "GōMartini",
        value=bool(st.session_state.go),
        key=_versioned_key("go"),
        help="Enables the Gō model support exposed by MartiniSurf/martinize2.",
    )
    st.session_state.elastic = False
    st.session_state.go_eps = st.number_input(
        "Go_Epsilon",
        min_value=0.0,
        value=float(st.session_state.go_eps),
        key=_versioned_key("go_eps"),
        step=0.001,
        format="%.3f",
    )
    position_options = ["backbone", "all", "none"]
    st.session_state.position_restraints = st.segmented_control(
        "Position restraints",
        position_options,
        default=_choice(st.session_state.position_restraints, "backbone"),
        key=_versioned_key("position_restraints"),
    )
    with st.expander("Advanced model options"):
        st.text_area("Extra martinize2 args", key="martinize_extra", height=80)


def _render_surface_step() -> None:
    st.markdown('<div class="ms-panel-title">Surface setup</div>', unsafe_allow_html=True)
    workflow_options = ["Hexagonal Lattice", "Carbon Nanomaterial", "Upload Surface File"]
    selected_workflow = st.selectbox(
        "Surface workflow",
        workflow_options,
        index=workflow_options.index(st.session_state.surface_workflow) if st.session_state.surface_workflow in workflow_options else 0,
        key=_versioned_key("surface_workflow"),
    )
    if selected_workflow != st.session_state.surface_workflow:
        st.session_state.surface_workflow = selected_workflow
        _apply_surface_workflow_defaults()
    st.caption("All generated surface workflows use Martini 3 compatible bead definitions.")
    if st.session_state.surface_workflow == "Hexagonal Lattice":
        st.session_state.surface_geometry = "planar"
        st.markdown("#### Hexagonal Lattice")
        a, b, c = st.columns(3)
        topology_options = ["4-1", "2-1"]
        st.session_state.surface_mode = a.selectbox(
            "Lattice topology",
            topology_options,
            index=topology_options.index(st.session_state.surface_mode) if st.session_state.surface_mode in topology_options else 0,
            key=_versioned_key("surface_mode_hex"),
        )
        st.session_state.dx = b.number_input(
            "Lattice bond length (nm)",
            min_value=0.01,
            value=float(st.session_state.dx),
            key=_versioned_key("dx_hex"),
            step=0.01,
        )
        st.session_state.surface_beads = c.text_input("Surface bead type", value=st.session_state.surface_beads, key=_versioned_key("surface_beads_hex"))
        d, e, f = st.columns(3)
        st.session_state.lx = d.number_input("Surface size X (nm)", min_value=0.1, value=float(st.session_state.lx), key=_versioned_key("lx_hex"), step=0.5)
        st.session_state.ly = e.number_input("Surface size Y (nm)", min_value=0.1, value=float(st.session_state.ly), key=_versioned_key("ly_hex"), step=0.5)
        st.session_state.surface_layers = f.number_input("Lattice layers", min_value=1, value=int(st.session_state.surface_layers), key=_versioned_key("surface_layers_hex"), step=1)
        g, h = st.columns(2)
        stacking_options = ["hcp", "fcc"]
        st.session_state.surface_stacking = g.selectbox(
            "Lattice stacking",
            stacking_options,
            index=stacking_options.index(st.session_state.surface_stacking) if st.session_state.surface_stacking in stacking_options else 0,
            key=_versioned_key("surface_stacking_hex"),
            format_func=lambda value: "HCP / ABAB" if value == "hcp" else "FCC / ABC",
        )
        st.session_state.charge = h.number_input("Charge", value=float(st.session_state.charge), key=_versioned_key("charge_hex"), step=0.1)
    elif st.session_state.surface_workflow == "Carbon Nanomaterial":
        st.markdown("#### Carbon Nanomaterial")
        carbon_options = ["Graphene", "Graphite", "Nanotubes"]
        kind = st.selectbox(
            "Carbon surface kind",
            carbon_options,
            index=carbon_options.index(st.session_state.carbon_surface_kind) if st.session_state.carbon_surface_kind in carbon_options else 0,
            key=_versioned_key("carbon_surface_kind"),
        )
        if kind != st.session_state.carbon_surface_kind:
            st.session_state.carbon_surface_kind = kind
            _apply_carbon_defaults()
        st.session_state.surface_beads = "C1"
        st.session_state.surface_geometry = "3d" if kind == "Nanotubes" else "planar"
        st.session_state.surface_mode = {"Graphene": "graphene", "Graphite": "graphite", "Nanotubes": "cnt-martini3"}[kind]
        a, b = st.columns(2)
        st.session_state.lx = a.number_input("Carbon edge length X (nm)", min_value=0.1, value=float(st.session_state.lx), key=_versioned_key(f"lx_carbon_{kind}"), step=0.5)
        st.session_state.ly = b.number_input("Carbon edge length Y (nm)", min_value=0.1, value=float(st.session_state.ly), key=_versioned_key(f"ly_carbon_{kind}"), step=0.5)
        if kind == "Graphite":
            st.session_state.graphite_layers = st.number_input("Graphite layers", min_value=1, value=int(st.session_state.graphite_layers), key=_versioned_key("graphite_layers"), step=1)
        if kind == "Nanotubes":
            c, d = st.columns(2)
            st.session_state.cnt_numrings = c.number_input("CNT rings", min_value=1, value=int(st.session_state.cnt_numrings), key=_versioned_key("cnt_numrings"), step=1)
            st.session_state.cnt_ringsize = d.number_input("CNT ring size", min_value=1, value=int(st.session_state.cnt_ringsize), key=_versioned_key("cnt_ringsize"), step=1)
    else:
        st.markdown("#### Upload Surface File")
        st.info("Upload `surface.gro` in Structure > Optional auxiliary files. MartiniSurf will use that file instead of generating a local lattice.")
    st.session_state.surface_linkers = st.number_input(
        "Linkers decorating the surface",
        min_value=0,
        value=int(st.session_state.surface_linkers),
        key=_versioned_key("surface_linkers"),
        step=1,
        help="Number of linker molecules placed as surface decoration. Use 0 when the linker is only for protein immobilization.",
    )
    _render_linker_generator()
    with st.expander("Advanced surface options"):
        geometry_options = ["planar", "3d"]
        st.session_state.surface_geometry = st.segmented_control(
            "Surface geometry",
            geometry_options,
            default=_choice(st.session_state.surface_geometry, "planar"),
            key=_versioned_key("surface_geometry"),
        )
    _render_surface_preview()


def _render_surface_preview() -> None:
    st.markdown(
        """
        <div class="ms-surface-preview">
          <span>Lightweight geometric preview</span>
          <div></div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _render_orientation_step(config: BuildConfig) -> None:
    st.markdown('<div class="ms-panel-title">Orientation</div>', unsafe_allow_html=True)
    mode_col, dist_col, tag_col = st.columns([1.4, 1, 1], gap="large")
    mode_col.segmented_control("Mode", ["Anchor", "Linker", "Adsorption"], key="orientation_mode")
    dist_col.number_input("Target distance (nm)", min_value=0.0, key="dist", step=0.1)
    tag_col.toggle("His-tag orientation", key="histag")
    if st.session_state.orientation_mode == "Linker":
        _render_linker_group_editor()
        linker_path = _resolve_optional_path(st.session_state.get("linker_path_text", ""), Path(st.session_state.run_root))
        if linker_path and linker_path.exists():
            st.caption(f"Active linker: `{linker_path.name}`")
        else:
            st.warning("Generate a linker in Surface > Linker Generator or upload a linker file in Structure before building.")
        with st.expander("Advanced linker parameters"):
            c, d = st.columns(2)
            c.number_input("Linker-protein distance (nm)", min_value=0.0, key="linker_prot_dist", step=0.1)
            d.number_input("Linker-surface distance (nm)", min_value=0.0, key="linker_surf_dist", step=0.1)
    else:
        _render_anchor_group_editor()
        if st.session_state.orientation_mode == "Adsorption":
            st.info("Adsorption uses anchor-based orientation but skips pull/restraint topology generation.")


def _render_linker_group_editor() -> None:
    _load_linker_group_rows_from_text()
    count = int(st.session_state.get("linker_group_count", 1))
    for index in range(count):
        st.session_state.setdefault(f"linker_group_chain_{index}", "A")
        st.session_state.setdefault(f"linker_group_residues_{index}", "73" if index == 0 else "")
        chain_col, residue_col, remove_col, add_col = st.columns([0.45, 1.45, 0.22, 0.22], gap="small")
        chain_col.text_input("Chain", key=f"linker_group_chain_{index}", max_chars=2, on_change=_sync_linker_group_text)
        residue_col.text_input(
            "Linker residues",
            key=f"linker_group_residues_{index}",
            placeholder="73 or 8 10 11",
            on_change=_sync_linker_group_text,
        )
        remove_col.button(
            "-",
            key=f"remove_linker_group_{index}",
            on_click=_remove_linker_group_row,
            args=(index,),
            width="stretch",
        )
        add_col.button(
            "+",
            key=f"add_linker_group_{index}",
            on_click=_add_linker_group_row,
            width="stretch",
        )
    _sync_linker_group_text()


def _render_anchor_group_editor() -> None:
    _load_anchor_rows_from_text()
    count = int(st.session_state.get("anchor_group_count", 1))
    for index in range(count):
        st.session_state.setdefault(f"anchor_chain_{index}", "A")
        st.session_state.setdefault(f"anchor_residues_{index}", "73" if index == 0 else "")
        chain_col, residue_col, remove_col, add_col = st.columns([0.45, 1.45, 0.22, 0.22], gap="small")
        chain_col.text_input("Chain", key=f"anchor_chain_{index}", max_chars=2, on_change=_sync_anchor_text)
        residue_col.text_input(
            "Anchor residues",
            key=f"anchor_residues_{index}",
            placeholder="73 or 8 10 11",
            on_change=_sync_anchor_text,
        )
        remove_col.button(
            "-",
            key=f"remove_anchor_{index}",
            help="Remove this anchor group",
            disabled=count <= 1,
            on_click=_remove_anchor_row,
            args=(index,),
        )
        if index == count - 1:
            add_col.button("+", key=f"add_anchor_{index}", help="Add another anchor group", on_click=_add_anchor_row)
    _sync_anchor_text()


def _render_workspace_preview(config: BuildConfig) -> None:
    orientation_mode = _choice(config.orientation_mode, "Anchor")
    anchor_count = len(config.anchors) if orientation_mode != "Linker" else len(config.linker_groups)
    st.markdown(
        f"""
        <div class="ms-canvas">
          <div class="ms-legend">
            <span><b class="cyan"></b>Protein</span>
            <span><b class="amber"></b>Selected groups</span>
            <span><b class="silver"></b>Surface</span>
          </div>
          <div class="ms-protein-shape"><i></i><i></i><i></i><i></i><i></i><i></i></div>
          <div class="ms-anchor-callout left">Groups: {anchor_count}</div>
          <div class="ms-anchor-callout right">{html.escape(orientation_mode)}</div>
          <div class="ms-surface-grid"></div>
          <div class="ms-axis">Z</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _render_environment_step() -> None:
    st.markdown('<div class="ms-panel-title">Environment</div>', unsafe_allow_html=True)
    a, b = st.columns(2)
    st.session_state.solvate = a.toggle("Solvate", value=bool(st.session_state.solvate), key=_versioned_key("solvate"))
    st.session_state.ionize = b.toggle("Ionize", value=bool(st.session_state.ionize), key=_versioned_key("ionize"))
    st.session_state.salt_conc = st.number_input(
        "Salt concentration (M)",
        min_value=0.0,
        value=float(st.session_state.salt_conc),
        key=_versioned_key("salt_conc"),
        step=0.01,
    )
    st.session_state.enable_water_mix = st.toggle(
        "Enable Martini 3 water mixing",
        value=bool(st.session_state.enable_water_mix),
        key=_versioned_key("enable_water_mix"),
    )
    mix_a, mix_b = st.columns(2)
    st.session_state.sw_water_percent = mix_a.slider(
        "SW water (%)",
        min_value=0,
        max_value=100,
        value=int(st.session_state.sw_water_percent),
        key=_versioned_key("sw_water_percent"),
    )
    max_tw = max(0, 100 - int(st.session_state.sw_water_percent))
    if int(st.session_state.tw_water_percent) > max_tw:
        st.session_state.tw_water_percent = max_tw
    st.session_state.tw_water_percent = mix_b.slider(
        "TW water (%)",
        min_value=0,
        max_value=max_tw,
        value=int(st.session_state.tw_water_percent),
        key=_versioned_key("tw_water_percent"),
    )
    st.session_state.water_mix_seed = st.number_input("Water mix seed", min_value=0, value=int(st.session_state.water_mix_seed), key=_versioned_key("water_mix_seed"), step=1)
    st.session_state.substrate_count = st.number_input("Substrate count", min_value=0, value=int(st.session_state.substrate_count), key=_versioned_key("substrate_count"), step=1)


def _render_review_step(config: BuildConfig, errors: list[str], tool_warnings: list[str], outdir: Path) -> None:
    left, right = st.columns([1.2, 0.9], gap="large")
    with left:
        _render_results(outdir, config)
        _render_command(config)
    with right:
        for error in errors:
            st.error(error)
        for warning in tool_warnings:
            st.warning(warning)
        run = st.button("Build system", type="primary", width="stretch", disabled=bool(errors or tool_warnings))
        reset = st.button("Reset working folder", width="stretch")
        if reset:
            shutil.rmtree(Path(st.session_state.run_root), ignore_errors=True)
            st.session_state.pop("run_root", None)
            st.rerun()
        if run:
            st.session_state["building"] = True
            with st.status("Running MartiniSurf", expanded=True) as status:
                result = run_martinisurf(build_args(config), Path(st.session_state.run_root), REPO_ROOT)
                st.session_state["last_stdout"] = result.stdout
                st.session_state["last_stderr"] = result.stderr
                st.session_state["last_returncode"] = result.returncode
                st.session_state["building"] = False
                if result.returncode == 0:
                    status.update(label="System generated", state="complete")
                    st.session_state["zip_path"] = str(make_zip(outdir, Path(st.session_state.run_root) / "Simulation_Files.zip"))
                    st.rerun()
                else:
                    status.update(label="Run failed", state="error")
        _render_log()


def _render_results(outdir: Path, config: BuildConfig) -> None:
    st.markdown('<div class="ms-panel-title">Generated structure</div>', unsafe_allow_html=True)
    structures = find_viewable_structures(outdir)
    if structures:
        structure_path = structures[0]
        st.caption(file_caption(structure_path))
        if structure_path.suffix.lower() == ".gro":
            with st.expander("Viewer Options", expanded=False):
                show_connectivity = st.toggle(
                    "Show protein topology connectivity",
                    value=True,
                    key="viewer_show_connectivity",
                    help="Draws short protein bonds as thin green cylinders. Surface lattice bonds are not drawn.",
                )
                a, b = st.columns(2)
                bead_radius = a.number_input("Bead radius", min_value=0.05, max_value=3.0, value=0.85, step=0.05, key="viewer_bead_radius")
                bond_radius = b.number_input("Bond radius", min_value=0.01, max_value=1.0, value=0.2, step=0.01, key="viewer_bond_radius")
                c, d = st.columns(2)
                highlight_radius = c.number_input("Highlight radius", min_value=0.05, max_value=3.0, value=1.1, step=0.05, key="viewer_highlight_radius")
                topology_bond_max_nm = d.number_input("Topology bond max (nm)", min_value=0.05, max_value=2.0, value=0.75, step=0.05, key="viewer_topology_bond_max_nm")
            selection_text = "\n".join(config.linker_groups if config.orientation_mode == "Linker" else config.anchors)
            stats = render_martini_step4_viewer(
                structure_path,
                outdir,
                selection_text,
                height=800,
                show_connectivity=show_connectivity,
                bead_radius=bead_radius,
                bond_radius=bond_radius,
                highlight_radius=highlight_radius,
                topology_bond_max_nm=topology_bond_max_nm,
            )
            st.caption(
                f"Viewer: {stats['bonds']} protein bonds drawn, "
                f"{stats['skipped_long']} long topology contacts skipped, "
                f"{stats['highlighted']} immobilization residues highlighted."
            )
        else:
            render_molecule(structure_path, height=520)
    else:
        st.info("No generated files yet.")
    _render_simulation_download(outdir)


def _render_simulation_download(outdir: Path) -> None:
    if not outdir.exists():
        return
    zip_path = Path(st.session_state["zip_path"]) if st.session_state.get("zip_path") else Path(st.session_state.run_root) / "Simulation_Files.zip"
    if not zip_path.exists():
        zip_path = make_zip(outdir, Path(st.session_state.run_root) / "Simulation_Files.zip")
        st.session_state["zip_path"] = str(zip_path)
    st.download_button(
        "Download MartiniSurf simulation files",
        zip_path.read_bytes(),
        "Simulation_Files.zip",
        "application/zip",
        width="stretch",
    )


def _render_command(config: BuildConfig) -> None:
    with st.expander("Reproducible command", expanded=True):
        st.code(shell_command(build_args(config)), language="bash")


def _render_log() -> None:
    returncode = st.session_state.get("last_returncode")
    if returncode is None:
        return
    stdout = st.session_state.get("last_stdout", "")
    stderr = st.session_state.get("last_stderr", "")
    st.markdown('<div class="ms-panel-title">Run summary</div>', unsafe_allow_html=True)
    for label, value in summarize_log(stdout, stderr, returncode):
        st.markdown(
            f'<div class="ms-log-row"><div class="ms-log-label">{html.escape(label)}</div><div class="ms-log-value">{html.escape(value)}</div></div>',
            unsafe_allow_html=True,
        )
    with st.expander("Raw execution details"):
        st.code((stdout + "\n" + stderr).strip(), language="text")


def _build_config(
    input_path: Path | str,
    surface_path: Path | None,
    linker_upload_path: Path | None,
    substrate_path: Path | None,
    substrate_itp_path: Path | None,
) -> BuildConfig:
    orientation_mode = _choice(st.session_state.orientation_mode, "Anchor")
    surface_workflow = st.session_state.surface_workflow
    carbon_kind = st.session_state.carbon_surface_kind
    water_mix = ""
    if st.session_state.enable_water_mix:
        sw_fraction = int(st.session_state.sw_water_percent) / 100.0
        tw_fraction = int(st.session_state.tw_water_percent) / 100.0
        w_fraction = max(0.0, 1.0 - sw_fraction - tw_fraction)
        water_mix = f"W:{w_fraction:.4f},SW:{sw_fraction:.4f},TW:{tw_fraction:.4f}"
    values = {
        "input_path": input_path,
        "workdir": Path(st.session_state.run_root),
        "moltype": st.session_state.moltype,
        "ff": st.session_state.ff,
        "dssp": st.session_state.dssp,
        "go": st.session_state.go,
        "go_eps": st.session_state.go_eps,
        "elastic": False,
        "maxwarn": st.session_state.maxwarn,
        "position_restraints": _choice(st.session_state.position_restraints, "backbone"),
        "merge_groups": _lines(st.session_state.merge_text),
        "surface_workflow": surface_workflow,
        "surface_mode": st.session_state.surface_mode,
        "surface_geometry": _choice(st.session_state.surface_geometry, "planar"),
        "surface_path": surface_path if surface_workflow == "Upload Surface File" else None,
        "lx": st.session_state.lx,
        "ly": st.session_state.ly,
        "dx": st.session_state.dx,
        "surface_beads": st.session_state.surface_beads.split(),
        "charge": st.session_state.charge,
        "surface_layers": st.session_state.surface_layers if surface_workflow == "Hexagonal Lattice" else None,
        "surface_stacking": st.session_state.surface_stacking,
        "surface_dist_z": st.session_state.surface_dist_z or None,
        "graphite_layers": st.session_state.graphite_layers if surface_workflow == "Carbon Nanomaterial" and carbon_kind == "Graphite" else None,
        "cnt_numrings": st.session_state.cnt_numrings if surface_workflow == "Carbon Nanomaterial" and carbon_kind == "Nanotubes" else None,
        "cnt_ringsize": st.session_state.cnt_ringsize if surface_workflow == "Carbon Nanomaterial" and carbon_kind == "Nanotubes" else None,
        "orientation_mode": orientation_mode,
        "anchors": _lines(st.session_state.anchors_text),
        "dist": st.session_state.dist,
        "ads_mode": orientation_mode == "Adsorption",
        "balance_low_z": False,
        "balance_low_z_fraction": None,
        "histag": st.session_state.histag,
        "linker_path": linker_upload_path or _resolve_optional_path(st.session_state.linker_path_text, Path(st.session_state.run_root)),
        "linker_groups": _lines(st.session_state.linker_groups_text),
        "linker_prot_dist": st.session_state.linker_prot_dist or None,
        "linker_surf_dist": st.session_state.linker_surf_dist or None,
        "linker_protein_bead": st.session_state.linker_protein_bead or None,
        "linker_surface_bead": st.session_state.linker_surface_bead or None,
        "invert_linker": st.session_state.invert_linker,
        "surface_linkers": st.session_state.surface_linkers,
        "substrate_path": substrate_path,
        "substrate_itp_path": substrate_itp_path,
        "substrate_count": st.session_state.substrate_count,
        "solvate": st.session_state.solvate,
        "ionize": st.session_state.ionize,
        "salt_conc": st.session_state.salt_conc,
        "water_mix": water_mix,
        "water_mix_seed": st.session_state.water_mix_seed,
        "martinize_extra_args": _lines(st.session_state.martinize_extra),
    }
    known_fields = {field.name for field in fields(BuildConfig)}
    return BuildConfig(**{key: value for key, value in values.items() if key in known_fields})


def main() -> None:
    st.set_page_config(page_title="MartiniSurf", page_icon="MS", layout="wide", initial_sidebar_state="expanded")
    apply_theme()
    _init_state()
    _commit_transient_widget_state()
    _preserve_app_state()
    RUNS_DIR.mkdir(exist_ok=True)
    st.session_state.setdefault("run_root", tempfile.mkdtemp(prefix="martinisurf_", dir=RUNS_DIR))
    upload_dir = Path(st.session_state.run_root) / "inputs"
    upload_dir.mkdir(parents=True, exist_ok=True)

    structure_path = st.session_state.get("_structure_path")
    surface_path = st.session_state.get("_surface_path")
    linker_upload_path = st.session_state.get("_linker_path")
    linker_itp_upload_path = st.session_state.get("_linker_itp_path")
    substrate_path = st.session_state.get("_substrate_path")
    substrate_itp_path = st.session_state.get("_substrate_itp_path")

    input_path: Path | str = Path(structure_path) if structure_path else st.session_state.remote_id
    summary = summarize_structure(
        Path(structure_path) if structure_path else None,
        st.session_state.remote_id,
        upload_dir / "preview",
    )
    config = _build_config(
        input_path,
        Path(surface_path) if surface_path else None,
        Path(linker_upload_path) if linker_upload_path else None,
        Path(substrate_path) if substrate_path else None,
        Path(substrate_itp_path) if substrate_itp_path else None,
    )
    errors = validate_config(config)
    tool_warnings = check_external_tools(config, TOOL_DIRS)
    outdir = Path(st.session_state.run_root) / config.outdir
    has_output = outdir.exists()

    _render_sidebar(errors, tool_warnings, has_output)
    _render_topbar(config, errors, tool_warnings, has_output)

    if st.session_state.active_step == "Structure":
        paths = _render_structure_step(upload_dir)
        structure_path, surface_path, linker_upload_path, linker_itp_upload_path, substrate_path, substrate_itp_path, summary = paths
        _remember_path("_structure_path", structure_path)
        _remember_path("_surface_path", surface_path)
        _remember_path("_linker_path", linker_upload_path)
        _remember_path("_linker_itp_path", linker_itp_upload_path)
        _remember_path("_substrate_path", substrate_path)
        _remember_path("_substrate_itp_path", substrate_itp_path)
    elif st.session_state.active_step == "Model":
        _render_model_step(summary)
    elif st.session_state.active_step == "Surface":
        _render_surface_step()
    elif st.session_state.active_step == "Orientation":
        _render_orientation_step(config)
    elif st.session_state.active_step == "Environment":
        _render_environment_step()
    else:
        _render_review_step(config, errors, tool_warnings, outdir)


if __name__ == "__main__":
    main()
