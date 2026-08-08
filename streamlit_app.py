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
        st.title("Protein Designer")
        preset_name = st.selectbox("Preset", list(PRESETS), index=0)
        preset = PRESETS[preset_name]
        st.caption("Generate MartiniSurf design files without running molecular dynamics.")

    st.markdown(
        f"""
        <div class="ms-header">
          <div>
            <div class="ms-kicker">Protein workflow</div>
            <div class="ms-brand">MartiniSurf</div>
            <p class="ms-soft">Build coarse-grained protein-surface systems, inspect the generated structure, and export reproducible setup files.</p>
            <div class="ms-chip-row">
              <span class="ms-chip">Protein</span>
              <span class="ms-chip">Martini 3</span>
              <span class="ms-chip">No mdrun</span>
              <span class="ms-chip">Local preview</span>
            </div>
          </div>
          <div class="ms-hero-logo"><img src="{_image_data_uri(logo)}" /></div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    run_root = Path(st.session_state.get("run_root", tempfile.mkdtemp(prefix="martinisurf_", dir=RUNS_DIR)))
    st.session_state["run_root"] = str(run_root)
    uploads_dir = run_root / "inputs"
    uploads_dir.mkdir(parents=True, exist_ok=True)

    input_col, setup_col = st.columns([1.05, 0.95], gap="large")

    with input_col:
        st.subheader("Input")
        uploaded_pdb = st.file_uploader("Protein structure", type=["pdb", "cif", "mmcif"])
        remote_id = st.text_input("RCSB or UniProt ID", value="1RJW" if uploaded_pdb is None else "")
        uploaded_surface = st.file_uploader("Existing surface .gro", type=["gro"])
        uploaded_linker = st.file_uploader("Linker .gro", type=["gro"])
        uploaded_substrate = st.file_uploader("Substrate .gro", type=["gro"])
        uploaded_substrate_itp = st.file_uploader("Substrate .itp", type=["itp"])

    with setup_col:
        st.subheader("Protein")
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

    surface_tab, orientation_tab, optional_tab, run_tab = st.tabs(
        ["Surface", "Orientation", "Optional Prep", "Run"]
    )

    with surface_tab:
        c1, c2, c3 = st.columns(3)
        surface_mode = c1.selectbox("Mode", ["4-1", "2-1", "graphene", "graphene-periodic", "graphene-finite", "graphite", "cnt"], index=0)
        surface_geometry = c2.selectbox("Geometry", ["planar", "3d"], index=0)
        surface_beads = c3.text_input("Beads", value="P4")
        d1, d2, d3 = st.columns(3)
        lx = d1.number_input("X size nm", min_value=0.1, value=15.0, step=0.5)
        ly = d2.number_input("Y size nm", min_value=0.1, value=15.0, step=0.5)
        dx = d3.number_input("Spacing nm", min_value=0.01, value=0.47, step=0.01)
        a1, a2, a3 = st.columns(3)
        charge = a1.number_input("Charge", value=0.0, step=0.1)
        surface_layers = a2.number_input("Layers", min_value=0, value=0, step=1)
        surface_stacking = a3.selectbox("Stacking", ["hcp", "fcc"], index=0)
        with st.expander("CNT and graphite"):
            graphite_layers = st.number_input("Graphite layers", min_value=0, value=0, step=1)
            cnt_numrings = st.number_input("CNT rings", min_value=0, value=0, step=1)
            cnt_ringsize = st.number_input("CNT ring size", min_value=0, value=0, step=1)

    with orientation_tab:
        orientation_mode = st.radio(
            "Mode",
            ["Anchor", "Adsorption", "Linker"],
            index=["Anchor", "Adsorption", "Linker"].index("Linker" if preset_name == "Linker" else ("Adsorption" if preset["ads"] else "Anchor")),
            horizontal=True,
        )
        dist = st.number_input("Surface distance nm", min_value=0.0, value=1.0 if preset_name != "Full Prep" else 1.0, step=0.1)
        if orientation_mode == "Linker":
            default_linker = preset["linker"]
            linker_path_text = st.text_input("Linker path", value=default_linker if uploaded_linker is None else "")
            linker_groups_text = st.text_area("Linker groups", value=preset["linker_groups"], height=104)
            l1, l2, l3 = st.columns(3)
            linker_prot_dist = l1.number_input("Protein-linker dist", min_value=0.0, value=0.0, step=0.1)
            linker_surf_dist = l2.number_input("Surface-linker dist", min_value=0.0, value=0.0, step=0.1)
            surface_linkers = l3.number_input("Surface linkers", min_value=0, value=12 if preset_name == "Linker" else 0, step=1)
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
        o1, o2, o3 = st.columns(3)
        solvate = o1.toggle("Solvate", value=bool(preset["solvate"]))
        ionize = o2.toggle("Ionize", value=bool(preset["ionize"]))
        salt_conc = o3.number_input("Salt M", min_value=0.0, value=0.15, step=0.01)
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

    with run_tab:
        st.subheader("Design snapshot")
        snap_cols = st.columns(5)
        snap_cols[0].metric("Preset", preset_name)
        snap_cols[1].metric("Surface", surface_mode)
        snap_cols[2].metric("Orientation", orientation_mode)
        snap_cols[3].metric("Input", str(input_path)[:12] or "none")
        snap_cols[4].metric("Solvent", "on" if solvate else "off")
        with st.expander("Reproducible command", expanded=False):
            st.code(shell_command(args), language="bash")
        if errors:
            for error in errors:
                st.error(error)
        if tool_warnings:
            for warning in tool_warnings:
                st.warning(warning)
        run = st.button(
            "Generate System",
            type="primary",
            width="stretch",
            disabled=bool(errors or tool_warnings),
        )
        if st.button("Reset working folder", width="stretch"):
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
        outdir = run_root / config.outdir
        viewer_col, output_col = st.columns([1.25, 0.75], gap="large")
        with viewer_col:
            st.subheader("Molecular viewer")
            _show_viewer(outdir)
        with output_col:
            stdout = st.session_state.get("last_stdout", "")
            stderr = st.session_state.get("last_stderr", "")
            returncode = st.session_state.get("last_returncode")
            _show_log(stdout, stderr, returncode)
            _show_results(outdir, zip_path)


if __name__ == "__main__":
    main()
