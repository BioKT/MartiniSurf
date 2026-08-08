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
        "system_final.gro",
        "immobilized_system.gro",
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
      const viewer = $3Dmol.createViewer(element, {{ backgroundColor: "#fbfdfb" }});
      viewer.addModel({json.dumps(mol_data)}, {json.dumps(fmt)});
      viewer.setStyle({{}}, {{ sphere: {{ radius: 0.18, colorscheme: "Jmol" }} }});
      viewer.addStyle({{resn: "SRF"}}, {{ sphere: {{ radius: 0.22, color: "#ff4fa7" }} }});
      viewer.addStyle({{resn: "W"}}, {{ sphere: {{ radius: 0.08, color: "#62b6d8" }} }});
      viewer.zoomTo();
      viewer.render();
    </script>
    <style>
      .viewer-shell {{
        width: 100%;
        height: {height}px;
        border: 1px solid #ecd4df;
        border-radius: 8px;
        overflow: hidden;
        background: #fbfdfb;
      }}
      .viewer {{
        width: 100%;
        height: {height}px;
      }}
    </style>
    """
    components.html(script, height=height + 2)


def file_caption(path: Path) -> str:
    size_kb = path.stat().st_size / 1024
    return f"{html.escape(path.name)} - {size_kb:.1f} KB"
