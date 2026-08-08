from __future__ import annotations

import base64
import html
import shutil
import sys
import tempfile
from pathlib import Path

import streamlit as st

from streamlit_app.archive import make_zip
from streamlit_app.command_builder import BuildConfig, build_args, shell_command
from streamlit_app.logs import summarize_log
from streamlit_app.molecular_viewer import file_caption, find_viewable_structures, render_molecule
from streamlit_app.runner import run_martinisurf
from streamlit_app.theme import apply_theme
from streamlit_app.validators import check_external_tools, validate_config


REPO_ROOT = Path(__file__).resolve().parent
RUNS_DIR = REPO_ROOT / "streamlit_runs"
TOOL_DIRS = [Path(sys.executable).resolve().parent, REPO_ROOT / ".venv" / "bin"]


PRESETS = {
    "Anchor": {
        "anchors": "A 8 10 11\nD 8 10 11",
        "merge": "A,B,C,D",
        "ads": False,
        "linker": "",
        "linker_groups": "",
        "solvate": False,
        "ionize": False,
        "go": False,
    },
    "Adsorption": {
        "anchors": "A 8 10 11\nD 8 10 11",
        "merge": "A,B,C,D",
        "ads": True,
        "linker": "",
        "linker_groups": "",
        "solvate": False,
        "ionize": False,
        "go": False,
    },
    "Linker": {
        "anchors": "",
        "merge": "A,B,C,D",
        "ads": False,
        "linker": "martinisurf/examples/protein/03_linker_surface_decoration/inputs/EPOXY.gro",
        "linker_groups": "A 8 10 11\nD 8 10 11",
        "solvate": False,
        "ionize": False,
        "go": False,
    },
    "Full Prep": {
        "anchors": "A 8 10 11\nD 8 10 11",
        "merge": "A,B,C,D",
        "ads": False,
        "linker": "",
        "linker_groups": "",
        "solvate": True,
        "ionize": True,
        "go": True,
    },
}


def _lines(value: str) -> list[str]:
    return [line.strip() for line in value.splitlines() if line.strip()]


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


def _image_data_uri(path: Path) -> str:
    if not path.exists():
        return ""
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _status_dot(label: str, active: bool = False, done: bool = False) -> str:
    class_name = "done" if done else ("active" if active else "")
    return f'<div class="ms-step {class_name}"><span></span><p>{html.escape(label)}</p></div>'


def _show_sidebar_progress(active_step: str, has_output: bool) -> None:
    steps = ["Structure", "Model", "Surface", "Orientation", "Environment", "Review & Build"]
    done_steps = {"Structure", "Model", "Surface"} if active_step != "Structure" else set()
    if has_output:
        done_steps = set(steps)
    markup = "\n".join(
        _status_dot(step, active=step == active_step and not has_output, done=step in done_steps)
        for step in steps
    )
    st.markdown(f'<div class="ms-steps">{markup}</div>', unsafe_allow_html=True)
    completed = len(done_steps)
    st.progress(completed / len(steps), text=f"{completed} of {len(steps)}")


def _show_topbar(logo: Path, preset_name: str, input_label: str, go: bool, has_output: bool) -> None:
    st.markdown(
        f"""
        <div class="ms-topbar">
          <div class="ms-top-brand">
            <img src="{_image_data_uri(logo)}" />
            <div>
              <div class="ms-top-title">MartiniSurf Studio</div>
              <div class="ms-top-subtitle">{html.escape(input_label or "Protein design")}</div>
            </div>
          </div>
          <div class="ms-top-pills">
            <span>{html.escape(preset_name)}</span>
            <span>Martini 3</span>
            <span>{'Go model' if go else 'Standard model'}</span>
          </div>
          <div class="ms-ready {'on' if has_output else ''}"><span></span>{'Generated' if has_output else 'Ready'}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _show_design_canvas(outdir: Path, config: BuildConfig) -> None:
    structures = find_viewable_structures(outdir)
    if structures:
        st.markdown('<div class="ms-panel-title">Molecular workspace</div>', unsafe_allow_html=True)
        _show_viewer(outdir)
        return

    anchor_count = len(config.anchors) if config.orientation_mode != "Linker" else len(config.linker_groups)
    st.markdown(
        f"""
        <div class="ms-canvas">
          <div class="ms-legend">
            <span><b class="cyan"></b>Protein</span>
            <span><b class="amber"></b>Anchors</span>
            <span><b class="silver"></b>Surface</span>
          </div>
          <div class="ms-protein-shape">
            <i></i><i></i><i></i><i></i><i></i><i></i>
          </div>
          <div class="ms-anchor-callout left">Groups: {anchor_count}</div>
          <div class="ms-anchor-callout right">{html.escape(config.orientation_mode)}</div>
          <div class="ms-surface-grid"></div>
          <div class="ms-axis">Z</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _show_spec_cards(config: BuildConfig) -> None:
    solvent = "Water + ions" if config.ionize else ("Water" if config.solvate else "Dry setup")
    st.markdown(
        f"""
        <div class="ms-spec-grid">
          <div class="ms-spec"><span>Surface</span><strong>{html.escape(config.surface_mode)} · {' '.join(config.surface_beads)}</strong></div>
          <div class="ms-spec"><span>Box</span><strong>{config.lx:g} × {config.ly:g} nm</strong></div>
          <div class="ms-spec"><span>Solvent</span><strong>{html.escape(solvent)}</strong></div>
          <div class="ms-spec"><span>Output</span><strong>GROMACS-ready</strong></div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _show_log(stdout: str, stderr: str, returncode: int | None) -> None:
    if returncode is None:
        return

    st.subheader("Run summary")
    for label, value in summarize_log(stdout, stderr, returncode):
        st.markdown(
            f"""
            <div class="ms-log-row">
              <div class="ms-log-label">{html.escape(label)}</div>
              <div class="ms-log-value">{html.escape(value)}</div>
            </div>
            """,
            unsafe_allow_html=True,
        )

    with st.expander("Raw execution details"):
        st.code((stdout + "\n" + stderr).strip() or f"Return code: {returncode}", language="text")


def _show_viewer(outdir: Path) -> None:
    structures = find_viewable_structures(outdir)
    if not structures:
        st.info("Generate a system to activate the molecular viewer.")
        return

    labels = [file_caption(path) for path in structures]
    selected_label = st.selectbox("Structure", labels, index=0)
    selected_path = structures[labels.index(selected_label)]
    render_molecule(selected_path)


def _show_results(outdir: Path, zip_path: Path | None) -> None:
    st.subheader("Output")
    if not outdir.exists():
        st.info("No output folder was created yet.")
        return

    files = sorted(path.relative_to(outdir.parent) for path in outdir.rglob("*") if path.is_file())
    cols = st.columns(4)
    cols[0].metric("Files", len(files))
    cols[1].metric("Topology", "yes" if (outdir / "0_topology" / "system.top").exists() else "no")
    cols[2].metric("System", "yes" if (outdir / "2_system").exists() else "no")
    cols[3].metric("MDP", "yes" if (outdir / "1_mdp").exists() else "no")

    with st.expander("File browser", expanded=False):
        top_files = [str(path) for path in files[:120]]
        st.code("\n".join(top_files) or "No files found.", language="text")

    if zip_path and zip_path.exists():
        st.download_button(
            "Download Simulation_Files.zip",
            data=zip_path.read_bytes(),
            file_name="Simulation_Files.zip",
            mime="application/zip",
            width="stretch",
        )


def main() -> None:
    st.set_page_config(
        page_title="MartiniSurf Protein Designer",
        page_icon="MS",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    apply_theme()
    RUNS_DIR.mkdir(exist_ok=True)

    with st.sidebar:
        logo = REPO_ROOT / "logo.png"
        if logo.exists():
            st.image(str(logo), width="stretch")
        st.markdown('<div class="ms-sidebar-title">MartiniSurf</div>', unsafe_allow_html=True)
        preset_name = st.selectbox("Preset", list(PRESETS), index=0)
        preset = PRESETS[preset_name]
        st.caption("Design files only. No molecular dynamics run.")

    run_root = Path(st.session_state.get("run_root", tempfile.mkdtemp(prefix="martinisurf_", dir=RUNS_DIR)))
    st.session_state["run_root"] = str(run_root)
    uploads_dir = run_root / "inputs"
    uploads_dir.mkdir(parents=True, exist_ok=True)

    top_slot = st.container()
    workspace_col, controls_col = st.columns([1.55, 0.85], gap="large")
    workspace_slot = workspace_col.container()

    with controls_col:
        st.markdown('<div class="ms-control-title">Build controls</div>', unsafe_allow_html=True)
        input_tab, model_tab, surface_tab, orientation_tab, optional_tab = st.tabs(
            ["Structure", "Model", "Surface", "Orient", "Env"]
        )

        with input_tab:
            uploaded_pdb = st.file_uploader("Protein structure", type=["pdb", "cif", "mmcif"])
            remote_id = st.text_input("RCSB or UniProt ID", value="1RJW" if uploaded_pdb is None else "")
            uploaded_surface = st.file_uploader("Existing surface .gro", type=["gro"])
            uploaded_linker = st.file_uploader("Linker .gro", type=["gro"])
            uploaded_substrate = st.file_uploader("Substrate .gro", type=["gro"])
            uploaded_substrate_itp = st.file_uploader("Substrate .itp", type=["itp"])

        with model_tab:
            moltype = st.text_input("Molecule name", value="Protein")
            ff = st.text_input("Force field", value="martini3001")
            merge_text = st.text_area("Merge groups", value=preset["merge"], height=72)
            martinize_cols = st.columns(3)
            dssp = martinize_cols[0].toggle("DSSP", value=True)
            go = martinize_cols[1].toggle("GoMartini", value=bool(preset["go"]))
            elastic = martinize_cols[2].toggle("Elastic", value=False)
            position_restraints = st.segmented_control(
                "Position restraints",
                ["backbone", "all", "none"],
                default="backbone",
            )

        with surface_tab:
            surface_mode = st.selectbox(
                "Mode",
                ["4-1", "2-1", "graphene", "graphene-periodic", "graphene-finite", "graphite", "cnt"],
                index=0,
            )
            surface_geometry = st.segmented_control("Geometry", ["planar", "3d"], default="planar")
            surface_beads = st.text_input("Beads", value="P4")
            d1, d2 = st.columns(2)
            lx = d1.number_input("X size nm", min_value=0.1, value=15.0, step=0.5)
            ly = d2.number_input("Y size nm", min_value=0.1, value=15.0, step=0.5)
            d3, d4 = st.columns(2)
            dx = d3.number_input("Spacing nm", min_value=0.01, value=0.47, step=0.01)
            charge = d4.number_input("Charge", value=0.0, step=0.1)
            with st.expander("Layer, CNT, graphite"):
                surface_layers = st.number_input("Layers", min_value=0, value=0, step=1)
                surface_stacking = st.selectbox("Stacking", ["hcp", "fcc"], index=0)
                graphite_layers = st.number_input("Graphite layers", min_value=0, value=0, step=1)
                cnt_numrings = st.number_input("CNT rings", min_value=0, value=0, step=1)
                cnt_ringsize = st.number_input("CNT ring size", min_value=0, value=0, step=1)

        with orientation_tab:
            orientation_mode = st.segmented_control(
                "Mode",
                ["Anchor", "Linker", "Adsorption"],
                default="Linker" if preset_name == "Linker" else ("Adsorption" if preset["ads"] else "Anchor"),
            )
            dist = st.number_input("Target distance nm", min_value=0.0, value=1.0, step=0.1)
            if orientation_mode == "Linker":
                default_linker = preset["linker"]
                linker_path_text = st.text_input("Linker path", value=default_linker if uploaded_linker is None else "")
                linker_groups_text = st.text_area("Linker groups", value=preset["linker_groups"], height=104)
                l1, l2 = st.columns(2)
                linker_prot_dist = l1.number_input("Protein-linker dist", min_value=0.0, value=0.0, step=0.1)
                linker_surf_dist = l2.number_input("Surface-linker dist", min_value=0.0, value=0.0, step=0.1)
                surface_linkers = st.number_input("Surface linkers", min_value=0, value=12 if preset_name == "Linker" else 0, step=1)
                invert_linker = st.toggle("Invert linker", value=False)
                anchors_text = ""
                ads_mode = False
                balance_low_z = False
                histag = False
            else:
                anchors_text = st.text_area("Anchor groups", value=preset["anchors"], height=104)
                ads_mode = orientation_mode == "Adsorption"
                balance_low_z = st.toggle("Balance low-Z region", value=False)
                histag = st.toggle("His-tag orientation", value=False)
                linker_path_text = ""
                linker_groups_text = ""
                linker_prot_dist = None
                linker_surf_dist = None
                invert_linker = False
                surface_linkers = 0

        with optional_tab:
            solvate = st.toggle("Solvate", value=bool(preset["solvate"]))
            ionize = st.toggle("Ionize", value=bool(preset["ionize"]))
            salt_conc = st.number_input("Salt M", min_value=0.0, value=0.15, step=0.01)
            water_mix = st.text_input("Water mix", placeholder="SW:0.10,TW:0.10")
            substrate_count = st.number_input("Substrate count", min_value=0, value=0, step=1)
            martinize_extra = st.text_area("Advanced martinize2 args", height=72)

    pdb_path = _save_upload(uploaded_pdb, uploads_dir)
    surface_path = _save_upload(uploaded_surface, uploads_dir)
    linker_upload_path = _save_upload(uploaded_linker, uploads_dir)
    substrate_path = _save_upload(uploaded_substrate, uploads_dir)
    substrate_itp_path = _save_upload(uploaded_substrate_itp, uploads_dir)

    input_path = pdb_path if pdb_path else remote_id.strip()
    config = BuildConfig(
        input_path=input_path,
        workdir=run_root,
        moltype=moltype,
        ff=ff,
        dssp=dssp,
        go=go,
        elastic=elastic,
        position_restraints=position_restraints or "backbone",
        merge_groups=_lines(merge_text),
        surface_mode=surface_mode,
        surface_geometry=surface_geometry,
        surface_path=surface_path,
        lx=lx,
        ly=ly,
        dx=dx,
        surface_beads=surface_beads.split(),
        charge=charge,
        surface_layers=surface_layers or None,
        surface_stacking=surface_stacking,
        graphite_layers=graphite_layers or None,
        cnt_numrings=cnt_numrings or None,
        cnt_ringsize=cnt_ringsize or None,
        orientation_mode=orientation_mode,
        anchors=_lines(anchors_text),
        dist=dist,
        ads_mode=ads_mode,
        balance_low_z=balance_low_z,
        histag=histag,
        linker_path=linker_upload_path or _resolve_optional_path(linker_path_text, run_root),
        linker_groups=_lines(linker_groups_text),
        linker_prot_dist=linker_prot_dist or None,
        linker_surf_dist=linker_surf_dist or None,
        invert_linker=invert_linker,
        surface_linkers=int(surface_linkers),
        substrate_path=substrate_path,
        substrate_itp_path=substrate_itp_path,
        substrate_count=int(substrate_count),
        solvate=solvate,
        ionize=ionize,
        salt_conc=salt_conc,
        water_mix=water_mix,
        martinize_extra_args=_lines(martinize_extra),
    )
    args = build_args(config)
    errors = validate_config(config)
    tool_warnings = check_external_tools(config, TOOL_DIRS)
    outdir = run_root / config.outdir
    has_output = outdir.exists()

    with st.sidebar:
        _show_sidebar_progress("Review & Build" if has_output else "Orientation", has_output)

    with top_slot:
        _show_topbar(logo, preset_name, str(input_path), go, has_output)

    with workspace_slot:
        _show_design_canvas(outdir, config)
        _show_spec_cards(config)

    build_col, output_col = st.columns([1.35, 0.95], gap="large")
    with build_col:
        st.markdown('<div class="ms-command-bar">', unsafe_allow_html=True)
        with st.expander("Reproducible command", expanded=False):
            st.code(shell_command(args), language="bash")
        st.markdown("</div>", unsafe_allow_html=True)

        if errors:
            for error in errors:
                st.error(error)
        if tool_warnings:
            for warning in tool_warnings:
                st.warning(warning)

    with output_col:
        action_a, action_b = st.columns([1.6, 1])
        run = action_a.button(
            "Build system",
            type="primary",
            width="stretch",
            disabled=bool(errors or tool_warnings),
        )
        reset = action_b.button("Reset", width="stretch")
        if reset:
            shutil.rmtree(run_root, ignore_errors=True)
            st.session_state.pop("run_root", None)
            st.rerun()

        if run:
            with st.status("Running MartiniSurf", expanded=True) as status:
                result = run_martinisurf(args, run_root, REPO_ROOT)
                st.session_state["last_stdout"] = result.stdout
                st.session_state["last_stderr"] = result.stderr
                st.session_state["last_returncode"] = result.returncode
                if result.returncode == 0:
                    status.update(label="System generated", state="complete")
                    zip_path = make_zip(run_root / config.outdir, run_root / "Simulation_Files.zip")
                    st.session_state["zip_path"] = str(zip_path)
                else:
                    status.update(label="Run failed", state="error")

        zip_path = Path(st.session_state["zip_path"]) if st.session_state.get("zip_path") else None
        stdout = st.session_state.get("last_stdout", "")
        stderr = st.session_state.get("last_stderr", "")
        returncode = st.session_state.get("last_returncode")
        _show_log(stdout, stderr, returncode)
        _show_results(outdir, zip_path)


if __name__ == "__main__":
    main()
