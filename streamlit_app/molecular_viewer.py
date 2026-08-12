from __future__ import annotations

import html
import io
import json
import subprocess
import tempfile
from pathlib import Path

import streamlit.components.v1 as components


VIEWABLE_SUFFIXES = {".pdb": "pdb", ".gro": "gro"}
MARTINISURF_PINK = "#FF4FA3"
MARTINISURF_ANCHOR_PINK = "#D91B72"
MAX_TRAJECTORY_GIF_FRAMES = 10


def find_viewable_structures(outdir: Path) -> list[Path]:
    system_dir = outdir / "2_system"
    if not system_dir.exists():
        return []

    preferred = [
        "immobilized_system.gro",
    ]
    files = [path for path in system_dir.iterdir() if path.suffix.lower() in VIEWABLE_SUFFIXES]
    by_name = {path.name: path for path in files}
    return [by_name[name] for name in preferred if name in by_name]


def render_molecule(path: Path, height: int = 520) -> None:
    mol_data = path.read_text(errors="replace")
    fmt = VIEWABLE_SUFFIXES.get(path.suffix.lower(), "pdb")
    element_id = "viewer_" + "".join(ch if ch.isalnum() else "_" for ch in path.name)

    script = f"""
    <div class="viewer-shell">
      <div id="{element_id}" class="viewer"></div>
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const element = document.getElementById({json.dumps(element_id)});
      const viewer = $3Dmol.createViewer(element, {{ backgroundColor: "#07131C" }});
      viewer.addModel({json.dumps(mol_data)}, {json.dumps(fmt)});
      viewer.setStyle({{}}, {{ sphere: {{ radius: 0.18, colorscheme: "Jmol" }} }});
      viewer.addStyle({{resn: "SRF"}}, {{ sphere: {{ radius: 0.22, color: "#ff4fa7" }} }});
      viewer.addStyle({{resn: "W"}}, {{ sphere: {{ radius: 0.08, color: "#62b6d8" }} }});
      viewer.zoomTo();
      viewer.render();
    </script>
    <style>
      html, body {{
        margin: 0;
        padding: 0;
        overflow: hidden;
        background: #07131C;
      }}
      .viewer-shell {{
        width: 100%;
        height: {height}px;
        box-sizing: border-box;
        border: 1px solid rgba(116, 152, 170, 0.28);
        border-radius: 16px;
        overflow: hidden;
        background: #07131C;
        box-shadow: inset 0 0 42px rgba(53, 201, 211, 0.06);
      }}
      .viewer {{
        width: 100%;
        height: {height}px;
        overflow: hidden;
      }}
    </style>
    """
    components.html(script, height=height + 2)


def render_structure_preview(path: Path, height: int = 430) -> None:
    mol_data = path.read_text(errors="replace")
    fmt = VIEWABLE_SUFFIXES.get(path.suffix.lower(), "pdb")
    element_id = "viewer_preview_" + "".join(ch if ch.isalnum() else "_" for ch in path.name)

    script = f"""
    <div class="viewer-shell">
      <div id="{element_id}" class="viewer"></div>
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const element = document.getElementById({json.dumps(element_id)});
      const viewer = $3Dmol.createViewer(element, {{ backgroundColor: "#07131C" }});
      viewer.addModel({json.dumps(mol_data)}, {json.dumps(fmt)});
      viewer.setStyle({{}}, {{ cartoon: {{ color: {json.dumps(MARTINISURF_PINK)} }} }});
      viewer.addStyle({{hetflag: true}}, {{ stick: {{ radius: 0.16, colorscheme: "Jmol" }} }});
      viewer.zoomTo();
      viewer.render();
    </script>
    <style>
      html, body {{
        margin: 0;
        padding: 0;
        overflow: hidden;
        background: #07131C;
      }}
      .viewer-shell {{
        width: 100%;
        height: {height}px;
        box-sizing: border-box;
        border: 1px solid rgba(116, 152, 170, 0.28);
        border-radius: 16px;
        overflow: hidden;
        background: #07131C;
        box-shadow: inset 0 0 42px rgba(53, 201, 211, 0.06);
      }}
      .viewer {{
        width: 100%;
        height: {height}px;
        overflow: hidden;
      }}
    </style>
    """
    components.html(script, height=height + 2)


def render_remote_molecule(identifier: str, height: int = 430, color: str = MARTINISURF_PINK) -> None:
    clean_id = "".join(ch for ch in identifier.strip().upper() if ch.isalnum())
    if not clean_id:
        return
    element_id = f"viewer_remote_{clean_id}"
    script = f"""
    <div class="viewer-shell">
      <div id="{element_id}" class="viewer"></div>
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const element = document.getElementById({json.dumps(element_id)});
      const viewer = $3Dmol.createViewer(element, {{ backgroundColor: "#07131C" }});
      $3Dmol.download("pdb:{clean_id}", viewer, {{}}, function() {{
        viewer.setStyle({{}}, {{ cartoon: {{ color: {json.dumps(color)} }} }});
        viewer.addStyle({{hetflag: true}}, {{ stick: {{ radius: 0.16, colorscheme: "Jmol" }} }});
        viewer.zoomTo();
        viewer.render();
      }});
    </script>
    <style>
      html, body {{
        margin: 0;
        padding: 0;
        overflow: hidden;
        background: #07131C;
      }}
      .viewer-shell {{
        width: 100%;
        height: {height}px;
        box-sizing: border-box;
        border: 1px solid rgba(116, 152, 170, 0.28);
        border-radius: 16px;
        overflow: hidden;
        background: #07131C;
        box-shadow: inset 0 0 42px rgba(53, 201, 211, 0.06);
      }}
      .viewer {{
        width: 100%;
        height: {height}px;
        overflow: hidden;
      }}
    </style>
    """
    components.html(script, height=height + 2)


def _parse_selection_groups(raw_text: str) -> list[list[str]]:
    groups = []
    for raw_line in (raw_text or "").replace(";", "\n").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = line.replace(",", " ").split()
        if len(parts) >= 2:
            groups.append(parts)
    return groups


def _build_chain_to_global_resid_map(cleaned_pdb_path: Path) -> dict[tuple[str, int], int]:
    mapping = {}
    seen = set()
    global_resid = 0
    for raw in cleaned_pdb_path.read_text(errors="replace").splitlines():
        if not raw.startswith(("ATOM", "HETATM")) or len(raw) < 26:
            continue
        chain = (raw[21].strip() or "_").upper()
        try:
            local_resid = int(raw[22:26].strip())
        except ValueError:
            continue
        key = (chain, local_resid)
        if key in seen:
            continue
        seen.add(key)
        global_resid += 1
        mapping[key] = global_resid
    return mapping


def immobilization_resids_from_selection(raw_selection: str, outdir: Path) -> list[int]:
    groups = _parse_selection_groups(raw_selection)
    if not groups:
        return []

    chain_mode_items = []
    numeric_resids = []
    for parts in groups:
        head = parts[0].upper()
        if head.lstrip("+-").isdigit():
            numeric_resids.extend(int(tok) for tok in parts[1:] if tok.lstrip("+-").isdigit())
        else:
            chain_mode_items.extend((head, int(tok)) for tok in parts[1:] if tok.lstrip("+-").isdigit())

    resolved = set(numeric_resids)
    cleaned_pdb = outdir / "2_system" / "cleaned_input.pdb"
    if chain_mode_items and cleaned_pdb.exists():
        chain_map = _build_chain_to_global_resid_map(cleaned_pdb)
        for chain, local_resid in chain_mode_items:
            global_resid = chain_map.get((chain.upper(), local_resid))
            if global_resid is not None:
                resolved.add(global_resid)

    return sorted(resid for resid in resolved if resid > 0)


def _parse_gro_atoms(gro_path: Path) -> list[dict[str, float | int | str]]:
    lines = gro_path.read_text(errors="replace").splitlines()
    atoms = []
    for raw in lines[2:-1]:
        if len(raw) < 44:
            continue
        try:
            atoms.append(
                {
                    "serial": int(raw[15:20]),
                    "resid": int(raw[0:5]),
                    "resn": raw[5:10].strip(),
                    "name": raw[10:15].strip(),
                    "x": float(raw[20:28]),
                    "y": float(raw[28:36]),
                    "z": float(raw[36:44]),
                }
            )
        except ValueError:
            continue
    return atoms


def _component_atom_serials(atoms: list[dict[str, float | int | str]]) -> dict[str, list[int] | bool]:
    surface_resn = {"SRF", "GRA", "CNT"}
    water_resn = {"W", "WF", "SW", "TW", "SOL"}
    ion_resn = {"NA", "CL", "ION", "K", "CA", "MG", "ZN", "LI", "RB", "CS", "BA", "SR", "F", "BR", "I"}
    by_residue: dict[tuple[int, str], set[str]] = {}
    for atom in atoms:
        resid = int(atom["resid"])
        resn = str(atom["resn"]).strip().upper()
        name = str(atom["name"]).strip().upper()
        by_residue.setdefault((resid, resn), set()).add(name)
    protein_residues = {
        residue
        for residue, names in by_residue.items()
        if "BB" in names or any(name.startswith("SC") for name in names) or any(name.startswith("BB") and name[2:].isdigit() for name in names)
    }

    surface: list[int] = []
    water: list[int] = []
    ions: list[int] = []
    linker: list[int] = []
    backbone: list[int] = []
    sidechain: list[int] = []

    serials = [int(atom["serial"]) for atom in atoms]
    for atom in atoms:
        serial = int(atom["serial"])
        resid = int(atom["resid"])
        resn = str(atom["resn"]).strip().upper()
        name = str(atom["name"]).strip().upper()
        if resn in surface_resn:
            surface.append(serial)
        elif resn in water_resn:
            water.append(serial)
        elif resn in ion_resn:
            ions.append(serial)
        elif (resid, resn) in protein_residues:
            if name == "BB" or (name.startswith("BB") and name[2:].isdigit()):
                backbone.append(serial)
            else:
                sidechain.append(serial)
        else:
            linker.append(serial)
    return {
        "surface": surface,
        "water": water,
        "ions": ions,
        "linker": linker,
        "backbone": backbone,
        "sidechain": sidechain,
        "contiguous": bool(serials) and serials == list(range(1, len(serials) + 1)),
    }


def _highlighted_protein_serials(atoms: list[dict[str, float | int | str]], highlighted_resids: list[int]) -> list[int]:
    if not highlighted_resids:
        return []
    highlighted = set(int(resid) for resid in highlighted_resids)
    component_serials = _component_atom_serials(atoms)
    protein_serials = set(component_serials["backbone"])  # type: ignore[arg-type]
    return [
        int(atom["serial"])
        for atom in atoms
        if int(atom["resid"]) in highlighted and int(atom["serial"]) in protein_serials
        and _is_backbone_bead_name(str(atom["name"]))
    ]


def _is_backbone_bead_name(name: str) -> bool:
    clean = name.strip().upper()
    return clean == "BB" or (clean.startswith("BB") and clean[2:].isdigit())


def _parse_itp_bonds(itp_path: Path) -> list[tuple[int, int]]:
    bonds = []
    section = None
    skip_elastic_block = False
    skip_go_block = False
    for raw in itp_path.read_text(errors="ignore").splitlines():
        raw_upper = raw.upper()
        stripped = raw.strip()
        if stripped.startswith("[") and "]" in stripped:
            section = stripped.strip("[]").strip().lower()
            skip_elastic_block = False
            skip_go_block = False
            continue

        directive = stripped.upper()
        if directive.startswith("#IFDEF GO_VIRT") or directive.startswith("#IFDEF RUBBER_BANDS"):
            skip_go_block = True
            continue
        if skip_go_block:
            if directive.startswith("#ENDIF"):
                skip_go_block = False
            continue

        if stripped.startswith(";"):
            if "LONG ELASTIC BONDS" in raw_upper or "SHORT ELASTIC BONDS" in raw_upper:
                skip_elastic_block = True
            elif "BACKBONE BONDS" in raw_upper or "SIDE CHAIN BONDS" in raw_upper:
                skip_elastic_block = False
            continue
        if skip_elastic_block:
            continue

        line = raw.split(";", 1)[0].strip()
        if not line or line.startswith("#") or section not in {"bonds", "constraints"}:
            continue
        if "RUBBER" in line.upper() or "GO" in line.upper():
            continue

        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            left, right = int(parts[0]), int(parts[1])
        except ValueError:
            continue
        bonds.append((left, right))
    return bonds


def _topology_bonds_for_system(outdir: Path) -> list[tuple[int, int]]:
    itp_dir = outdir / "0_topology" / "system_itp"
    if not itp_dir.exists():
        return []

    skip_words = ("surface", "linker", "posre", "martini_", "go_", "ions", "solvents")
    bonds = []
    for itp_path in sorted(itp_dir.glob("*.itp")):
        lower_name = itp_path.name.lower()
        if any(word in lower_name for word in skip_words):
            continue
        bonds.extend(_parse_itp_bonds(itp_path))
    return sorted(set(tuple(sorted(pair)) for pair in bonds if pair[0] != pair[1]))


def _distance2(atom_a: dict[str, float | int | str], atom_b: dict[str, float | int | str]) -> float:
    return (
        (float(atom_a["x"]) - float(atom_b["x"])) ** 2
        + (float(atom_a["y"]) - float(atom_b["y"])) ** 2
        + (float(atom_a["z"]) - float(atom_b["z"])) ** 2
    )


def _short_bond_cylinders(
    atoms_by_serial: dict[int, dict[str, float | int | str]],
    bonds: list[tuple[int, int]],
    max_distance_nm: float,
) -> tuple[list[dict[str, dict[str, float]]], int]:
    cylinders = []
    skipped_long = 0
    max_distance2 = float(max_distance_nm) ** 2
    for left, right in bonds:
        atom_a = atoms_by_serial.get(left)
        atom_b = atoms_by_serial.get(right)
        if atom_a is None or atom_b is None:
            continue
        if _distance2(atom_a, atom_b) > max_distance2:
            skipped_long += 1
            continue
        cylinders.append(
            {
                "start": {"x": 10.0 * float(atom_a["x"]), "y": 10.0 * float(atom_a["y"]), "z": 10.0 * float(atom_a["z"])},
                "end": {"x": 10.0 * float(atom_b["x"]), "y": 10.0 * float(atom_b["y"]), "z": 10.0 * float(atom_b["z"])},
            }
        )
    return cylinders, skipped_long


def render_martini_step4_viewer(
    path: Path,
    outdir: Path,
    selection_text: str,
    *,
    height: int = 800,
    show_connectivity: bool = True,
    bead_radius: float = 0.85,
    bond_radius: float = 0.2,
    highlight_radius: float = 1.1,
    topology_bond_max_nm: float = 0.75,
) -> dict[str, int]:
    mol_data = path.read_text(errors="replace")
    fmt = VIEWABLE_SUFFIXES.get(path.suffix.lower(), "gro")
    element_id = "viewer_" + "".join(ch if ch.isalnum() else "_" for ch in path.name)

    atoms = _parse_gro_atoms(path) if path.suffix.lower() == ".gro" else []
    atoms_by_serial = {int(atom["serial"]): atom for atom in atoms}
    cylinders, skipped_long = ([], 0)
    if show_connectivity and atoms_by_serial:
        cylinders, skipped_long = _short_bond_cylinders(
            atoms_by_serial,
            _topology_bonds_for_system(outdir),
            topology_bond_max_nm,
        )
    highlighted_resids = immobilization_resids_from_selection(selection_text, outdir)
    highlighted_serials = _highlighted_protein_serials(atoms, highlighted_resids)
    component_resnames = _short_md_component_resnames(path)

    script = f"""
    <div class="viewer-shell step4">
      <div id="{element_id}" class="viewer"></div>
      <div class="viewer-badge">Visual quality check</div>
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const element = document.getElementById({json.dumps(element_id)});
      const viewer = $3Dmol.createViewer(element, {{ backgroundColor: "#07131C" }});
      viewer.addModel({json.dumps(mol_data)}, {json.dumps(fmt)});
      viewer.setStyle({{}}, {{ sphere: {{ radius: {float(bead_radius):.4f} }} }});
      const components = {json.dumps(component_resnames)};
      if (components.surface.length) {{
        viewer.setStyle({{resn: components.surface}}, {{ sphere: {{ radius: {float(bead_radius):.4f}, color: "#42C7D5" }} }});
      }}
      if (components.water.length) {{
        viewer.setStyle({{resn: components.water}}, {{ sphere: {{ radius: {float(bead_radius):.4f}, color: "lightgray", opacity: 0.55 }} }});
      }}
      if (components.ions.length) {{
        viewer.setStyle({{resn: components.ions}}, {{ sphere: {{ radius: {float(bead_radius):.4f}, color: "limegreen" }} }});
      }}
      if (components.linker.length) {{
        viewer.setStyle({{resn: components.linker}}, {{ sphere: {{ radius: {float(bead_radius):.4f}, color: "yellow" }} }});
      }}
      if (components.protein.length) {{
        viewer.setStyle({{resn: components.protein}}, {{ sphere: {{ radius: {float(bead_radius):.4f}, color: "#FF9DCC" }} }});
      }}
      viewer.setStyle({{atom: "BB"}}, {{ sphere: {{ radius: {float(bead_radius):.4f}, color: "#FF4FA3" }} }});
      viewer.setStyle({{atom: /^BB\\d+$/}}, {{ sphere: {{ radius: {float(bead_radius):.4f}, color: "#FF4FA3" }} }});
      viewer.setStyle({{atom: /^SC/}}, {{ sphere: {{ radius: {float(bead_radius):.4f}, color: "#FF9DCC" }} }});
      for (const cylinder of {json.dumps(cylinders)}) {{
        viewer.addCylinder({{
          start: cylinder.start,
          end: cylinder.end,
          radius: {float(bond_radius):.4f},
          color: "#F4F5F7",
          fromCap: 1,
          toCap: 1
        }});
      }}
      const highlightedResids = {json.dumps(highlighted_resids)};
      if (highlightedResids.length) {{
        viewer.addStyle({{ resi: highlightedResids, atom: "BB" }}, {{ sphere: {{ color: "{MARTINISURF_ANCHOR_PINK}", radius: {float(highlight_radius):.4f} }} }});
        viewer.addStyle({{ resi: highlightedResids, atom: /^BB\\d+$/ }}, {{ sphere: {{ color: "{MARTINISURF_ANCHOR_PINK}", radius: {float(highlight_radius):.4f} }} }});
      }}
      viewer.zoomTo();
      viewer.render();
    </script>
    <style>
      html, body {{
        margin: 0;
        padding: 0;
        overflow: hidden;
        background: #07131C;
      }}
      .viewer-shell {{
        position: relative;
        width: 100%;
        height: {height}px;
        box-sizing: border-box;
        border: 1px solid rgba(116, 152, 170, 0.28);
        border-radius: 16px;
        overflow: hidden;
        background: #07131C;
        box-shadow: inset 0 0 46px rgba(53, 201, 211, 0.06);
      }}
      .viewer {{
        width: 100%;
        height: {height}px;
        overflow: hidden;
      }}
      .viewer-badge {{
        position: absolute;
        left: 14px;
        bottom: 14px;
        padding: 8px 10px;
        border: 1px solid rgba(53, 201, 211, 0.30);
        border-radius: 999px;
        background: rgba(7, 19, 28, 0.78);
        color: #F3F7FA;
        font: 700 12px/1.2 sans-serif;
        pointer-events: none;
      }}
    </style>
    """
    components.html(script, height=height + 2)
    return {"bonds": len(cylinders), "skipped_long": skipped_long, "highlighted": len(highlighted_serials)}


def render_short_md_trajectory(
    gro_path: Path | None,
    tpr_path: Path | None,
    xtc_path: Path | None,
    outdir: Path,
    selection_text: str,
    gmx_bin: str | None,
    *,
    frame_stride: int = 1,
    height: int = 700,
    show_protein: bool = True,
    show_surface: bool = True,
    show_linker: bool = False,
    show_water: bool = False,
    show_ions: bool = False,
) -> dict[str, int | str]:
    element_id = "viewer_short_md_" + "".join(ch if ch.isalnum() else "_" for ch in str(gro_path or tpr_path or "trajectory"))
    highlighted_resids = immobilization_resids_from_selection(selection_text, outdir)
    component_resnames = _short_md_component_resnames(gro_path)
    gro_atoms = _parse_gro_atoms(gro_path) if gro_path and gro_path.exists() else []
    highlighted_serials = _highlighted_protein_serials(gro_atoms, highlighted_resids)
    pdb_data = ""
    model_mode = "structure"
    frame_count = 0
    conversion_error = ""

    if tpr_path and xtc_path and tpr_path.exists() and xtc_path.exists() and gmx_bin:
        try:
            with tempfile.TemporaryDirectory() as tmp:
                pdb_path = Path(tmp) / "short_md_frames.pdb"
                cmd = [
                    gmx_bin,
                    "trjconv",
                    "-s",
                    str(tpr_path),
                    "-f",
                    str(xtc_path),
                    "-o",
                    str(pdb_path),
                    "-skip",
                    str(max(1, int(frame_stride))),
                ]
                res = subprocess.run(cmd, input="0\n", text=True, capture_output=True, check=False)
                if res.returncode == 0 and pdb_path.exists() and pdb_path.stat().st_size > 0:
                    pdb_data = pdb_path.read_text(errors="replace")
                    model_mode = "trajectory"
                    frame_count = pdb_data.count("MODEL")
                else:
                    conversion_error = (res.stderr or res.stdout or "Trajectory conversion failed.")[-600:]
        except OSError as exc:
            conversion_error = str(exc)

    if not pdb_data and gro_path and gro_path.exists():
        pdb_data = gro_path.read_text(errors="replace")

    if not pdb_data:
        return {"frames": 0, "mode": "missing", "error": conversion_error or "No viewable Short MD structure was found."}

    script = f"""
    <div class="viewer-shell short-md">
      <div id="{element_id}" class="viewer"></div>
      <div class="viewer-badge">{'Production trajectory' if model_mode == 'trajectory' else 'Final production structure'}</div>
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const element = document.getElementById({json.dumps(element_id)});
      const viewer = $3Dmol.createViewer(element, {{ backgroundColor: "#0E0D11" }});
      const data = {json.dumps(pdb_data)};
      if ({json.dumps(model_mode)} === "trajectory") {{
        viewer.addModelsAsFrames(data, "pdb");
      }} else {{
        viewer.addModel(data, {json.dumps(VIEWABLE_SUFFIXES.get((gro_path or Path('x.gro')).suffix.lower(), 'gro'))});
      }}
      viewer.setStyle({{}}, {{ sphere: {{ hidden: true }} }});
      const components = {json.dumps(component_resnames)};
      if ({json.dumps(bool(show_surface))} && components.surface.length) {{
        viewer.setStyle({{resn: components.surface}}, {{ sphere: {{ radius: 0.935, color: "#42C7D5" }} }});
      }}
      if ({json.dumps(bool(show_water))} && components.water.length) {{
        viewer.setStyle({{resn: components.water}}, {{ sphere: {{ radius: 0.462, color: "lightgray", opacity: 0.55 }} }});
      }}
      if ({json.dumps(bool(show_ions))} && components.ions.length) {{
        viewer.setStyle({{resn: components.ions}}, {{ sphere: {{ radius: 0.605, color: "limegreen" }} }});
      }}
      if ({json.dumps(bool(show_linker))} && components.linker.length) {{
        viewer.setStyle({{resn: components.linker}}, {{ sphere: {{ radius: 0.836, color: "yellow" }} }});
      }}
      if ({json.dumps(bool(show_protein))}) {{
        viewer.setStyle({{atom: "BB"}}, {{ sphere: {{ radius: 0.902, color: "#FF4FA3" }} }});
        viewer.setStyle({{atom: /^BB\\d+$/}}, {{ sphere: {{ radius: 0.902, color: "#FF4FA3" }} }});
        viewer.setStyle({{atom: /^SC/}}, {{ sphere: {{ radius: 0.902, color: "#FF9DCC" }} }});
      }}
      if ({json.dumps(bool(show_protein))}) {{
        const highlightedResids = {json.dumps(highlighted_resids)};
        if (highlightedResids.length) {{
          viewer.addStyle({{ resi: highlightedResids, atom: "BB" }}, {{ sphere: {{ color: "{MARTINISURF_ANCHOR_PINK}", radius: 1.155 }} }});
          viewer.addStyle({{ resi: highlightedResids, atom: /^BB\\d+$/ }}, {{ sphere: {{ color: "{MARTINISURF_ANCHOR_PINK}", radius: 1.155 }} }});
        }}
      }}
      viewer.zoomTo();
      if ({json.dumps(model_mode)} === "trajectory") {{
        viewer.animate({{loop: "forward", reps: 0}});
      }}
      viewer.render();
    </script>
    <style>
      html, body {{
        margin: 0;
        padding: 0;
        overflow: hidden;
        background: #0E0D11;
      }}
      .viewer-shell {{
        position: relative;
        width: 100%;
        height: {height}px;
        box-sizing: border-box;
        border: 1px solid rgba(255, 98, 180, 0.28);
        border-radius: 16px;
        overflow: hidden;
        background: #0E0D11;
        box-shadow: inset 0 0 46px rgba(255, 98, 180, 0.08);
      }}
      .viewer {{
        width: 100%;
        height: {height}px;
        overflow: hidden;
      }}
      .viewer-badge {{
        position: absolute;
        left: 14px;
        bottom: 14px;
        padding: 8px 10px;
        border: 1px solid rgba(255, 98, 180, 0.35);
        border-radius: 999px;
        background: rgba(14, 13, 17, 0.78);
        color: #FAF5F8;
        font: 700 12px/1.2 sans-serif;
        pointer-events: none;
      }}
    </style>
    """
    components.html(script, height=height + 2)
    return {"frames": frame_count, "mode": model_mode, "error": conversion_error}


def generate_short_md_trajectory_gif(
    gro_path: Path | None,
    tpr_path: Path | None,
    xtc_path: Path | None,
    outdir: Path,
    selection_text: str,
    gmx_bin: str | None,
    *,
    frame_count_hint: int | None = None,
    show_protein: bool = True,
    show_surface: bool = True,
    show_linker: bool = False,
    show_water: bool = False,
    show_ions: bool = False,
    width: int = 960,
    height: int = 700,
) -> tuple[bytes | None, str | None, int]:
    if not gro_path or not tpr_path or not xtc_path or not gmx_bin:
        return None, "A GRO, TPR, XTC, and GROMACS are required to generate the trajectory GIF.", 0
    if not gro_path.exists() or not tpr_path.exists() or not xtc_path.exists():
        return None, "The selected trajectory files are not available anymore.", 0

    try:
        from PIL import Image
        from playwright.sync_api import sync_playwright
    except ImportError as exc:
        return None, f"GIF generation needs Playwright and Pillow in the Streamlit environment: {exc}", 0

    highlighted_resids = immobilization_resids_from_selection(selection_text, outdir)
    gro_atoms = _parse_gro_atoms(gro_path)
    component_resnames = _short_md_component_resnames(gro_path)
    frame_count_source = _trajectory_frame_count(xtc_path, gmx_bin) or frame_count_hint or MAX_TRAJECTORY_GIF_FRAMES
    stride = max(1, int(frame_count_source) // MAX_TRAJECTORY_GIF_FRAMES)

    with tempfile.TemporaryDirectory() as tmp_name:
        tmp = Path(tmp_name)
        pdb_path = tmp / "trajectory_sample.pdb"
        cmd = [
            gmx_bin,
            "trjconv",
            "-s",
            str(tpr_path),
            "-f",
            str(xtc_path),
            "-o",
            str(pdb_path),
            "-skip",
            str(stride),
        ]
        res = subprocess.run(cmd, input="0\n", text=True, capture_output=True, check=False)
        if res.returncode != 0 or not pdb_path.exists() or pdb_path.stat().st_size == 0:
            detail = (res.stderr or res.stdout or "Trajectory conversion failed.")[-600:]
            return None, detail, 0

        pdb_data = _sample_pdb_models(pdb_path.read_text(errors="replace"), MAX_TRAJECTORY_GIF_FRAMES)
        frame_count = max(1, pdb_data.count("MODEL"))
        html_path = tmp / "trajectory_gif.html"
        html_path.write_text(
            _short_md_gif_html(
                pdb_data,
                component_resnames,
                highlighted_resids,
                frame_count,
                show_protein=show_protein,
                show_surface=show_surface,
                show_linker=show_linker,
                show_water=show_water,
                show_ions=show_ions,
                width=width,
                height=height,
            ),
            encoding="utf-8",
        )

        frames = []
        try:
            with sync_playwright() as play:
                browser = play.chromium.launch(
                    executable_path=_chrome_executable(),
                    args=["--no-sandbox"],
                )
                page = browser.new_page(viewport={"width": width, "height": height}, device_scale_factor=1)
                page.goto(html_path.as_uri())
                page.wait_for_function("window.__martinisurfGifReady === true", timeout=20000)
                target = page.locator("#viewer")
                for frame_index in range(frame_count):
                    page.evaluate("frameIndex => window.__setMartiniSurfGifFrame(frameIndex)", frame_index)
                    page.wait_for_timeout(160)
                    png = target.screenshot(type="png")
                    frames.append(Image.open(io.BytesIO(png)).convert("RGB"))
                browser.close()
        except Exception as exc:  # noqa: BLE001 - report browser rendering failures in Streamlit.
            return None, f"Trajectory GIF rendering failed: {exc}", 0

        if not frames:
            return None, "No frames were rendered for the trajectory GIF.", 0

        gif = io.BytesIO()
        frames[0].save(gif, format="GIF", save_all=True, append_images=frames[1:], duration=240, loop=0, optimize=True)
        return gif.getvalue(), None, len(frames)


def _trajectory_frame_count(xtc_path: Path, gmx_bin: str) -> int | None:
    res = subprocess.run([gmx_bin, "check", "-f", str(xtc_path)], text=True, capture_output=True, check=False)
    if res.returncode != 0:
        return None
    for line in (res.stdout + "\n" + res.stderr).splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "Coords":
            try:
                return int(parts[1])
            except ValueError:
                return None
    return None


def _chrome_executable() -> str | None:
    for candidate in ("/usr/bin/google-chrome", "/usr/bin/chromium", "/usr/bin/chromium-browser"):
        path = Path(candidate)
        if path.exists():
            return str(path)
    return None


def _sample_pdb_models(pdb_data: str, max_frames: int) -> str:
    frames = _split_pdb_models(pdb_data)
    if not frames:
        return pdb_data
    if len(frames) <= max_frames:
        return "".join(frames)
    sampled = [frames[index] for index in _sample_indices(len(frames), max_frames)]
    return "".join(sampled)


def _split_pdb_models(pdb_data: str) -> list[str]:
    frames: list[str] = []
    current: list[str] = []
    in_model = False
    for line in pdb_data.splitlines(keepends=True):
        if line.startswith("MODEL"):
            if current:
                frames.append("".join(current))
            current = [line]
            in_model = True
        elif in_model:
            current.append(line)
            if line.startswith("ENDMDL"):
                frames.append("".join(current))
                current = []
                in_model = False
    if current:
        frames.append("".join(current))
    return frames


def _sample_indices(total: int, target: int) -> list[int]:
    if total <= target:
        return list(range(total))
    if target <= 1:
        return [0]
    return sorted({round(index * (total - 1) / (target - 1)) for index in range(target)})


def _short_md_gif_html(
    pdb_data: str,
    component_resnames: dict[str, list[str]],
    highlighted_resids: list[int],
    frame_count: int,
    *,
    show_protein: bool,
    show_surface: bool,
    show_linker: bool,
    show_water: bool,
    show_ions: bool,
    width: int,
    height: int,
) -> str:
    return f"""
    <!doctype html>
    <html>
      <head>
        <meta charset="utf-8">
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <style>
          html, body {{
            margin: 0;
            padding: 0;
            width: {width}px;
            height: {height}px;
            overflow: hidden;
            background: #0E0D11;
          }}
          #viewer {{
            width: {width}px;
            height: {height}px;
            background: #0E0D11;
          }}
        </style>
      </head>
      <body>
        <div id="viewer"></div>
        <script>
          const viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{ backgroundColor: "#0E0D11" }});
          viewer.addModelsAsFrames({json.dumps(pdb_data)}, "pdb");
          viewer.setStyle({{}}, {{ sphere: {{ hidden: true }} }});
          const components = {json.dumps(component_resnames)};
          if ({json.dumps(bool(show_surface))} && components.surface.length) {{
            viewer.setStyle({{resn: components.surface}}, {{ sphere: {{ radius: 0.935, color: "#42C7D5" }} }});
          }}
          if ({json.dumps(bool(show_water))} && components.water.length) {{
            viewer.setStyle({{resn: components.water}}, {{ sphere: {{ radius: 0.462, color: "lightgray", opacity: 0.55 }} }});
          }}
          if ({json.dumps(bool(show_ions))} && components.ions.length) {{
            viewer.setStyle({{resn: components.ions}}, {{ sphere: {{ radius: 0.605, color: "limegreen" }} }});
          }}
          if ({json.dumps(bool(show_linker))} && components.linker.length) {{
            viewer.setStyle({{resn: components.linker}}, {{ sphere: {{ radius: 0.836, color: "yellow" }} }});
          }}
          if ({json.dumps(bool(show_protein))}) {{
            viewer.setStyle({{atom: "BB"}}, {{ sphere: {{ radius: 0.902, color: "#FF4FA3" }} }});
            viewer.setStyle({{atom: /^BB\\d+$/}}, {{ sphere: {{ radius: 0.902, color: "#FF4FA3" }} }});
            viewer.setStyle({{atom: /^SC/}}, {{ sphere: {{ radius: 0.902, color: "#FF9DCC" }} }});
          }}
          if ({json.dumps(bool(show_protein))}) {{
            const highlightedResids = {json.dumps(highlighted_resids)};
            if (highlightedResids.length) {{
              viewer.addStyle({{ resi: highlightedResids, atom: "BB" }}, {{ sphere: {{ color: "{MARTINISURF_ANCHOR_PINK}", radius: 1.155 }} }});
              viewer.addStyle({{ resi: highlightedResids, atom: /^BB\\d+$/ }}, {{ sphere: {{ color: "{MARTINISURF_ANCHOR_PINK}", radius: 1.155 }} }});
            }}
          }}
          viewer.zoomTo();
          window.__setMartiniSurfGifFrame = async function(frameIndex) {{
            await viewer.setFrame(frameIndex);
            viewer.render();
          }};
          viewer.setFrame(0).then(() => {{
            viewer.rotate(-90, "x");
            viewer.render();
            window.__martinisurfGifReady = true;
          }});
        </script>
      </body>
    </html>
    """


def _short_md_component_resnames(gro_path: Path | None) -> dict[str, list[str]]:
    surface_resn = {"SRF", "GRA", "CNT"}
    water_resn = {"W", "WF", "SW", "TW", "SOL"}
    ion_resn = {"NA", "CL", "ION", "K", "CA", "MG", "ZN", "LI", "RB", "CS", "BA", "SR", "F", "BR", "I"}
    groups = {"surface": set(), "water": set(), "ions": set(), "protein": set(), "linker": set()}
    if not gro_path or not gro_path.exists():
        return {key: sorted(value) for key, value in groups.items()}

    by_residue: dict[tuple[int, str], set[str]] = {}
    for atom in _parse_gro_atoms(gro_path):
        resid = int(atom["resid"])
        resn = str(atom["resn"]).strip()
        name = str(atom["name"]).strip().upper()
        by_residue.setdefault((resid, resn), set()).add(name)

    for (_resid, resn), names in by_residue.items():
        if resn in surface_resn:
            groups["surface"].add(resn)
        elif resn in water_resn:
            groups["water"].add(resn)
        elif resn in ion_resn:
            groups["ions"].add(resn)
        elif "BB" in names or any(name.startswith("SC") for name in names) or any(name.startswith("BB") and name[2:].isdigit() for name in names):
            groups["protein"].add(resn)
        else:
            groups["linker"].add(resn)
    return {key: sorted(value) for key, value in groups.items()}


def render_linker_mapping(path: Path, beads: list[object], height: int = 360) -> None:
    mol_data = path.read_text(errors="replace")
    fmt = VIEWABLE_SUFFIXES.get(path.suffix.lower(), "gro")
    element_id = "viewer_linker_" + "".join(ch if ch.isalnum() else "_" for ch in path.name)
    labels = []
    for bead in beads:
        try:
            if isinstance(bead, dict):
                label = str(bead.get("label") or f"{bead.get('index', '')}: {bead.get('name', '')}".strip())
                x = bead.get("x")
                y = bead.get("y")
                z = bead.get("z")
            else:
                label = getattr(bead, "label", str(getattr(bead, "name", "")))
                x = getattr(bead, "x")
                y = getattr(bead, "y")
                z = getattr(bead, "z")
            labels.append(
                {
                    "label": label,
                    "x": 10.0 * float(x),
                    "y": 10.0 * float(y),
                    "z": 10.0 * float(z),
                }
            )
        except (TypeError, ValueError):
            continue

    script = f"""
    <div class="viewer-shell linker">
      <div id="{element_id}" class="viewer"></div>
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const element = document.getElementById({json.dumps(element_id)});
      const viewer = $3Dmol.createViewer(element, {{ backgroundColor: "#07131C" }});
      viewer.addModel({json.dumps(mol_data)}, {json.dumps(fmt)});
      viewer.setStyle({{}}, {{ sphere: {{ radius: 0.55, color: "#35c9d3" }} }});
      const labels = {json.dumps(labels)};
      for (const bead of labels) {{
        viewer.addLabel(bead.label, {{
          position: {{x: bead.x, y: bead.y, z: bead.z}},
          fontColor: "#07131C",
          backgroundColor: "#F7FBFF",
          borderColor: "#35c9d3",
          borderThickness: 1,
          fontSize: 13,
          inFront: true
        }});
      }}
      viewer.zoomTo();
      viewer.render();
    </script>
    <style>
      html, body {{
        margin: 0;
        padding: 0;
        overflow: hidden;
        background: #07131C;
      }}
      .viewer-shell {{
        width: 100%;
        height: {height}px;
        box-sizing: border-box;
        border: 1px solid rgba(116, 152, 170, 0.28);
        border-radius: 16px;
        overflow: hidden;
        background: #07131C;
        box-shadow: inset 0 0 42px rgba(53, 201, 211, 0.06);
      }}
      .viewer {{
        width: 100%;
        height: {height}px;
        overflow: hidden;
      }}
    </style>
    """
    components.html(script, height=height + 2)


def file_caption(path: Path) -> str:
    size_kb = path.stat().st_size / 1024
    return f"{html.escape(path.name)} - {size_kb:.1f} KB"
