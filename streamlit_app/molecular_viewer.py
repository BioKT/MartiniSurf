from __future__ import annotations

import html
import json
from pathlib import Path

import streamlit.components.v1 as components


VIEWABLE_SUFFIXES = {".pdb": "pdb", ".gro": "gro"}


def find_viewable_structures(outdir: Path) -> list[Path]:
    system_dir = outdir / "2_system"
    if not system_dir.exists():
        return []

    preferred = [
        "immobilized_system.gro",
        "system_final.gro",
        "Protein_cg.pdb",
        "Protein_cg.gro",
        "surface.gro",
        "cleaned_input.pdb",
    ]
    files = [path for path in system_dir.iterdir() if path.suffix.lower() in VIEWABLE_SUFFIXES]
    by_name = {path.name: path for path in files}
    ordered = [by_name[name] for name in preferred if name in by_name]
    ordered.extend(sorted(path for path in files if path.name not in preferred))
    return ordered


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


def render_remote_molecule(identifier: str, height: int = 430) -> None:
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
        viewer.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }} }});
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
      for (const cylinder of {json.dumps(cylinders)}) {{
        viewer.addCylinder({{
          start: cylinder.start,
          end: cylinder.end,
          radius: {float(bond_radius):.4f},
          color: "green",
          fromCap: 1,
          toCap: 1
        }});
      }}
      const highlightedResids = {json.dumps(highlighted_resids)};
      if (highlightedResids.length) {{
        viewer.addStyle({{ resi: highlightedResids }}, {{ sphere: {{ color: "red", radius: {float(highlight_radius):.4f} }} }});
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
    return {"bonds": len(cylinders), "skipped_long": skipped_long, "highlighted": len(highlighted_resids)}


def file_caption(path: Path) -> str:
    size_kb = path.stat().st_size / 1024
    return f"{html.escape(path.name)} - {size_kb:.1f} KB"
