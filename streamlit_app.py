from __future__ import annotations

import base64
import html
import shutil
import sys
import tempfile
from dataclasses import fields
from pathlib import Path

import streamlit as st

from streamlit_app.archive import make_zip
from streamlit_app.command_builder import BuildConfig, build_args, shell_command
from streamlit_app.logs import summarize_log
from streamlit_app.molecular_viewer import file_caption, find_viewable_structures, render_molecule
from streamlit_app.project_io import export_state, import_state
from streamlit_app.runner import run_martinisurf
from streamlit_app.structure import StructureSummary, summarize_structure
from streamlit_app.theme import apply_theme
from streamlit_app.validators import check_external_tools, validate_config


REPO_ROOT = Path(__file__).resolve().parent
RUNS_DIR = REPO_ROOT / "streamlit_runs"
TOOL_DIRS = [Path(sys.executable).resolve().parent, REPO_ROOT / ".venv" / "bin"]

STEPS = ["Structure", "Model", "Surface", "Orientation", "Environment", "Review & Build"]

PRESETS = {
    "Agarose-like": {"surface_mode": "4-1", "surface_beads": "P4", "dx": 0.47, "lx": 15.0, "ly": 15.0},
    "Graphene sheet": {"surface_mode": "graphene", "surface_beads": "C1", "dx": 0.47, "lx": 15.0, "ly": 15.0},
    "Graphite slab": {"surface_mode": "graphite", "surface_beads": "C1", "dx": 0.47, "lx": 15.0, "ly": 15.0},
    "Custom": {"surface_mode": "4-1", "surface_beads": "P4", "dx": 0.47, "lx": 15.0, "ly": 15.0},
}


def _init_state() -> None:
    defaults = {
        "project_name": "BsADH - 1RJW",
        "active_step": "Structure",
        "preset_name": "Agarose-like",
        "remote_id": "1RJW",
        "moltype": "Protein",
        "ff": "martini3001",
        "merge_text": "A,B,C,D",
        "dssp": True,
        "go": False,
        "elastic": False,
        "position_restraints": "backbone",
        "surface_mode": "4-1",
        "surface_geometry": "planar",
        "surface_beads": "P4",
        "lx": 15.0,
        "ly": 15.0,
        "dx": 0.47,
        "charge": 0.0,
        "surface_layers": 0,
        "surface_stacking": "hcp",
        "surface_dist_z": 0.0,
        "graphite_layers": 0,
        "cnt_numrings": 0,
        "cnt_ringsize": 0,
        "orientation_mode": "Anchor",
        "anchors_text": "A 8 10 11\nD 8 10 11",
        "dist": 1.0,
        "balance_low_z": False,
        "balance_low_z_fraction": 0.2,
        "histag": False,
        "linker_path_text": "martinisurf/examples/protein/03_linker_surface_decoration/inputs/EPOXY.gro",
        "linker_groups_text": "A 8 10 11\nD 8 10 11",
        "linker_prot_dist": 0.0,
        "linker_surf_dist": 0.0,
        "invert_linker": False,
        "surface_linkers": 0,
        "solvate": False,
        "ionize": False,
        "salt_conc": 0.15,
        "water_mix": "",
        "substrate_count": 0,
        "martinize_extra": "",
    }
    for key, value in defaults.items():
        st.session_state.setdefault(key, value)


def _lines(value: str) -> list[str]:
    return [line.strip() for line in value.splitlines() if line.strip()]


def _image_data_uri(path: Path) -> str:
    if not path.exists():
        return ""
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _save_upload(uploaded_file, target_dir: Path) -> Path | None:
    if uploaded_file is None:
        return None
    path = target_dir / uploaded_file.name
    path.write_bytes(uploaded_file.getbuffer())
    return path


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
        st.markdown('<div class="ms-sidebar-title">MartiniSurf Studio</div>', unsafe_allow_html=True)
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

    badges = ["Protein", "Martini 3", "Gō model" if config.go else "Standard model"]

    badge_html = "".join(f"<span>{html.escape(badge)}</span>" for badge in badges)
    st.markdown(
        f"""
        <div class="ms-topbar">
          <div class="ms-top-brand">
            <img src="{_image_data_uri(REPO_ROOT / 'logo.png')}" />
            <div>
              <div class="ms-top-title">MartiniSurf Studio</div>
              <div class="ms-top-subtitle">{html.escape(st.session_state.project_name)}</div>
            </div>
          </div>
          <div class="ms-top-pills">{badge_html}</div>
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
        file_name="martinisurf-studio-project.json",
        mime="application/json",
        width="stretch",
    )


def _render_structure_step(upload_dir: Path) -> tuple[Path | None, Path | None, Path | None, Path | None, Path | None, StructureSummary]:
    left, right = st.columns([0.9, 1.1], gap="large")
    with left:
        st.markdown('<div class="ms-panel-title">Structure input</div>', unsafe_allow_html=True)
        st.text_input("Project name", key="project_name")
        uploaded_structure = st.file_uploader("Protein structure", type=["pdb", "cif", "mmcif"])
        st.text_input("PDB ID or UniProt ID", key="remote_id")
        with st.expander("Optional auxiliary files"):
            uploaded_surface = st.file_uploader("Existing surface .gro", type=["gro"])
            uploaded_linker = st.file_uploader("Linker .gro", type=["gro"])
            uploaded_substrate = st.file_uploader("Substrate .gro", type=["gro"])
            uploaded_substrate_itp = st.file_uploader("Substrate .itp", type=["itp"])

    structure_path = _save_upload(uploaded_structure, upload_dir)
    surface_path = _save_upload(uploaded_surface, upload_dir)
    linker_path = _save_upload(uploaded_linker, upload_dir)
    substrate_path = _save_upload(uploaded_substrate, upload_dir)
    substrate_itp_path = _save_upload(uploaded_substrate_itp, upload_dir)
    summary = summarize_structure(structure_path, st.session_state.remote_id)

    with right:
        _render_structure_summary(summary, structure_path)

    return structure_path, surface_path, linker_path, substrate_path, substrate_itp_path, summary


def _render_structure_summary(summary: StructureSummary, structure_path: Path | None) -> None:
    st.markdown('<div class="ms-panel-title">Structure preview</div>', unsafe_allow_html=True)
    top_metrics = st.columns(2)
    bottom_metrics = st.columns(2)
    top_metrics[0].metric("Source", summary.source)
    top_metrics[1].metric("Type", summary.molecule_type)
    bottom_metrics[0].metric("Chains", len(summary.chains) if summary.chains else "-")
    bottom_metrics[1].metric("Residues", summary.residue_count or "-")

    if structure_path and structure_path.suffix.lower() in {".pdb", ".gro"}:
        render_molecule(structure_path, height=430)
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

    if summary.chains:
        st.dataframe(
            [{"Chain": item.chain, "Residues": item.residues, "Atoms": item.atoms or "-"} for item in summary.chains],
            hide_index=True,
            width="stretch",
        )
    if summary.ligands:
        st.info("Ligands/cofactors detected: " + ", ".join(summary.ligands[:12]))
    for note in summary.notes:
        st.warning(note)


def _render_model_step() -> None:
    st.markdown('<div class="ms-panel-title">Model decisions</div>', unsafe_allow_html=True)
    a, b = st.columns(2)
    a.text_input("Molecule name", key="moltype")
    b.text_input("Force field", key="ff")
    st.text_area("Merge chains", key="merge_text", height=80)
    c1, c2, c3 = st.columns(3)
    c1.toggle("DSSP", key="dssp", help="Passes DSSP handling to martinize2.")
    c2.toggle("GōMartini", key="go", help="Enables the Gō model support exposed by MartiniSurf/martinize2.")
    c3.toggle("Elastic network", key="elastic", help="Adds elastic-network options supported by martinize2.")
    st.segmented_control("Position restraints", ["backbone", "all", "none"], key="position_restraints")
    with st.expander("Advanced model options"):
        st.text_area("Extra martinize2 args", key="martinize_extra", height=80)


def _apply_surface_preset() -> None:
    preset = PRESETS[st.session_state.preset_name]
    for key, value in preset.items():
        st.session_state[key] = value


def _render_surface_step() -> None:
    st.markdown('<div class="ms-panel-title">Surface setup</div>', unsafe_allow_html=True)
    st.selectbox("Surface preset", list(PRESETS), key="preset_name", on_change=_apply_surface_preset)
    st.selectbox(
        "Surface mode",
        ["2-1", "4-1", "graphene", "graphene-periodic", "graphene-finite", "graphite", "cnt", "cnt-m2", "cnt-martini2", "cnt-m3", "cnt-martini3"],
        key="surface_mode",
    )
    st.segmented_control("Surface geometry", ["planar", "3d"], key="surface_geometry")
    a, b, c = st.columns(3)
    a.text_input("Bead type", key="surface_beads")
    b.number_input("Charge", key="charge", step=0.1)
    c.number_input("dx", min_value=0.01, key="dx", step=0.01)
    d, e = st.columns(2)
    d.number_input("lx", min_value=0.1, key="lx", step=0.5)
    e.number_input("ly", min_value=0.1, key="ly", step=0.5)
    with st.expander("Advanced surface options"):
        st.number_input("Surface layers", min_value=0, key="surface_layers", step=1)
        st.selectbox("Stacking", ["hcp", "fcc"], key="surface_stacking")
        st.number_input("Surface layer distance", min_value=0.0, key="surface_dist_z", step=0.01)
        st.number_input("Graphite layers", min_value=0, key="graphite_layers", step=1)
        st.number_input("CNT rings", min_value=0, key="cnt_numrings", step=1)
        st.number_input("CNT ring size", min_value=0, key="cnt_ringsize", step=1)
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
    left, right = st.columns([1.15, 0.85], gap="large")
    with left:
        _render_workspace_preview(config)
    with right:
        st.markdown('<div class="ms-panel-title">Orientation</div>', unsafe_allow_html=True)
        st.segmented_control("Mode", ["Anchor", "Linker", "Adsorption"], key="orientation_mode")
        st.number_input("Target distance (nm)", min_value=0.0, key="dist", step=0.1)
        if st.session_state.orientation_mode == "Linker":
            st.text_input("Linker path", key="linker_path_text")
            st.text_area("Linker groups", key="linker_groups_text", height=96)
            st.toggle("Invert linker", key="invert_linker")
            st.number_input("Linker-protein distance", min_value=0.0, key="linker_prot_dist", step=0.1)
            st.number_input("Linker-surface distance", min_value=0.0, key="linker_surf_dist", step=0.1)
            st.number_input("Additional surface linkers", min_value=0, key="surface_linkers", step=1)
        else:
            st.text_area("Anchor groups", key="anchors_text", height=96)
            st.toggle("Balance low-Z", key="balance_low_z")
            if st.session_state.balance_low_z:
                st.number_input("Balance low-Z fraction", min_value=0.01, max_value=1.0, key="balance_low_z_fraction", step=0.01)
            st.toggle("His-tag orientation", key="histag")
            if st.session_state.orientation_mode == "Adsorption":
                st.info("Adsorption uses anchor-based orientation but skips pull/restraint topology generation.")
        _render_orientation_quality(config)


def _render_workspace_preview(config: BuildConfig) -> None:
    anchor_count = len(config.anchors) if config.orientation_mode != "Linker" else len(config.linker_groups)
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
          <div class="ms-anchor-callout right">{html.escape(config.orientation_mode)}</div>
          <div class="ms-surface-grid"></div>
          <div class="ms-axis">Z</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _render_orientation_quality(config: BuildConfig) -> None:
    groups = config.linker_groups if config.orientation_mode == "Linker" else config.anchors
    checks = [
        ("Orientation mode selected", bool(config.orientation_mode)),
        ("Groups provided", bool(groups)),
        ("Target distance set", config.dist >= 0),
        ("Orientation generated", False),
    ]
    st.markdown('<div class="ms-quality"><strong>Orientation quality</strong>', unsafe_allow_html=True)
    for label, ok in checks:
        state = "ok" if ok else "pending"
        st.markdown(f'<div class="ms-quality-row {state}"><span></span>{html.escape(label)}</div>', unsafe_allow_html=True)
    st.caption("These are configuration checks only. Scientific validation requires inspecting generated coordinates.")
    st.markdown("</div>", unsafe_allow_html=True)


def _render_environment_step() -> None:
    st.markdown('<div class="ms-panel-title">Environment</div>', unsafe_allow_html=True)
    a, b = st.columns(2)
    a.toggle("Solvate", key="solvate")
    b.toggle("Ionize", key="ionize")
    st.number_input("Salt concentration (M)", min_value=0.0, key="salt_conc", step=0.01)
    st.text_input("Water mix", key="water_mix", placeholder="SW:0.10,TW:0.10")
    st.number_input("Substrate count", min_value=0, key="substrate_count", step=1)


def _render_review_step(config: BuildConfig, errors: list[str], tool_warnings: list[str], outdir: Path) -> None:
    left, right = st.columns([1.2, 0.9], gap="large")
    with left:
        _render_results(outdir)
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
                else:
                    status.update(label="Run failed", state="error")
        _render_log()


def _render_results(outdir: Path) -> None:
    st.markdown('<div class="ms-panel-title">Generated structure</div>', unsafe_allow_html=True)
    structures = find_viewable_structures(outdir)
    if structures:
        labels = [file_caption(path) for path in structures]
        selected = st.selectbox("Structure", labels)
        render_molecule(structures[labels.index(selected)], height=520)
    else:
        st.info("No generated files yet.")
    zip_path = Path(st.session_state["zip_path"]) if st.session_state.get("zip_path") else None
    if zip_path and zip_path.exists():
        st.download_button("Download Simulation_Files.zip", zip_path.read_bytes(), "Simulation_Files.zip", "application/zip", width="stretch")


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
    values = {
        "input_path": input_path,
        "workdir": Path(st.session_state.run_root),
        "moltype": st.session_state.moltype,
        "ff": st.session_state.ff,
        "dssp": st.session_state.dssp,
        "go": st.session_state.go,
        "elastic": st.session_state.elastic,
        "position_restraints": st.session_state.position_restraints,
        "merge_groups": _lines(st.session_state.merge_text),
        "surface_mode": st.session_state.surface_mode,
        "surface_geometry": st.session_state.surface_geometry,
        "surface_path": surface_path,
        "lx": st.session_state.lx,
        "ly": st.session_state.ly,
        "dx": st.session_state.dx,
        "surface_beads": st.session_state.surface_beads.split(),
        "charge": st.session_state.charge,
        "surface_layers": st.session_state.surface_layers or None,
        "surface_stacking": st.session_state.surface_stacking,
        "surface_dist_z": st.session_state.surface_dist_z or None,
        "graphite_layers": st.session_state.graphite_layers or None,
        "cnt_numrings": st.session_state.cnt_numrings or None,
        "cnt_ringsize": st.session_state.cnt_ringsize or None,
        "orientation_mode": st.session_state.orientation_mode,
        "anchors": _lines(st.session_state.anchors_text),
        "dist": st.session_state.dist,
        "ads_mode": st.session_state.orientation_mode == "Adsorption",
        "balance_low_z": st.session_state.balance_low_z,
        "balance_low_z_fraction": st.session_state.balance_low_z_fraction if st.session_state.balance_low_z else None,
        "histag": st.session_state.histag,
        "linker_path": linker_upload_path or _resolve_optional_path(st.session_state.linker_path_text, Path(st.session_state.run_root)),
        "linker_groups": _lines(st.session_state.linker_groups_text),
        "linker_prot_dist": st.session_state.linker_prot_dist or None,
        "linker_surf_dist": st.session_state.linker_surf_dist or None,
        "invert_linker": st.session_state.invert_linker,
        "surface_linkers": st.session_state.surface_linkers,
        "substrate_path": substrate_path,
        "substrate_itp_path": substrate_itp_path,
        "substrate_count": st.session_state.substrate_count,
        "solvate": st.session_state.solvate,
        "ionize": st.session_state.ionize,
        "salt_conc": st.session_state.salt_conc,
        "water_mix": st.session_state.water_mix,
        "martinize_extra_args": _lines(st.session_state.martinize_extra),
    }
    known_fields = {field.name for field in fields(BuildConfig)}
    return BuildConfig(**{key: value for key, value in values.items() if key in known_fields})


def main() -> None:
    st.set_page_config(page_title="MartiniSurf Studio", page_icon="MS", layout="wide", initial_sidebar_state="expanded")
    apply_theme()
    _init_state()
    RUNS_DIR.mkdir(exist_ok=True)
    st.session_state.setdefault("run_root", tempfile.mkdtemp(prefix="martinisurf_", dir=RUNS_DIR))
    upload_dir = Path(st.session_state.run_root) / "inputs"
    upload_dir.mkdir(parents=True, exist_ok=True)

    structure_path = st.session_state.get("_structure_path")
    surface_path = st.session_state.get("_surface_path")
    linker_upload_path = st.session_state.get("_linker_path")
    substrate_path = st.session_state.get("_substrate_path")
    substrate_itp_path = st.session_state.get("_substrate_itp_path")

    input_path: Path | str = Path(structure_path) if structure_path else st.session_state.remote_id
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
        structure_path, surface_path, linker_upload_path, substrate_path, substrate_itp_path, _ = paths
        st.session_state["_structure_path"] = str(structure_path) if structure_path else ""
        st.session_state["_surface_path"] = str(surface_path) if surface_path else ""
        st.session_state["_linker_path"] = str(linker_upload_path) if linker_upload_path else ""
        st.session_state["_substrate_path"] = str(substrate_path) if substrate_path else ""
        st.session_state["_substrate_itp_path"] = str(substrate_itp_path) if substrate_itp_path else ""
    elif st.session_state.active_step == "Model":
        _render_model_step()
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
