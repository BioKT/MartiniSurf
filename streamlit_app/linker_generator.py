from __future__ import annotations

import importlib.util
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class LinkerBead:
    index: int
    name: str
    resid: int
    resname: str
    x: float
    y: float
    z: float

    @property
    def label(self) -> str:
        return f"{self.index}: {self.name}"


@dataclass(frozen=True)
class LinkerGenerationResult:
    ok: bool
    gro_path: Path | None
    itp_path: Path | None
    beads: list[LinkerBead]
    message: str
    log: str
    generator: str | None = None


def safe_molecule_name(value: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9_]+", "_", value.strip())
    clean = clean.strip("_")
    return clean or "LINKER"


def parse_linker_gro(path: Path) -> list[LinkerBead]:
    lines = path.read_text(errors="replace").splitlines()
    beads: list[LinkerBead] = []
    for raw in lines[2:-1]:
        try:
            if len(raw) >= 44:
                resid = int(raw[0:5])
                resname = raw[5:10].strip()
                name = raw[10:15].strip()
                index = int(raw[15:20])
                x = float(raw[20:28])
                y = float(raw[28:36])
                z = float(raw[36:44])
            else:
                parts = raw.split()
                if len(parts) < 6:
                    continue
                resid_resname = parts[0]
                resid_text = "".join(ch for ch in resid_resname if ch.isdigit()) or "1"
                resname = resid_resname[len(resid_text):] or "LNK"
                resid = int(resid_text)
                name = parts[1]
                index = int(parts[2])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
            beads.append(LinkerBead(index=index, name=name, resid=resid, resname=resname, x=x, y=y, z=z))
        except ValueError:
            continue
    return beads


def smiles_to_svg(smiles: str) -> str:
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        return ""

    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return ""
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(520, 320)
    options = drawer.drawOptions()
    options.addAtomIndices = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _module_exists(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


def _generator_candidates(smiles: str, molname: str, output_dir: Path) -> list[tuple[str, list[str]]]:
    candidates: list[tuple[str, list[str]]] = []
    if _module_exists("martini_mapper"):
        candidates.append(
            (
                "Martini Mapper",
                [sys.executable, "-m", "martini_mapper", molname, smiles, "--no-xtb", "--out-dir", str(output_dir)],
            )
        )
    if _module_exists("Martini_mapper"):
        candidates.append(
            (
                "Martini Mapper",
                [sys.executable, "-m", "Martini_mapper", molname, smiles, "--no-xtb", "--out-dir", str(output_dir)],
            )
        )
    if _module_exists("auto_martiniM3"):
        candidates.append(
            (
                "AutoMartini M3",
                [sys.executable, "-m", "auto_martiniM3", "--mode", "run", "--smi", smiles, "--mol", molname],
            )
        )
    executable = shutil.which("martini_mapper")
    if executable:
        candidates.append(("Martini Mapper", [executable, molname, smiles, "--no-xtb", "--out-dir", str(output_dir)]))
    return candidates


def generate_linker_from_smiles(smiles: str, molname: str, output_dir: Path) -> LinkerGenerationResult:
    smiles = smiles.strip()
    molname = safe_molecule_name(molname)
    output_dir.mkdir(parents=True, exist_ok=True)
    if not smiles:
        return LinkerGenerationResult(False, None, None, [], "SMILES is empty.", "", None)

    candidates = _generator_candidates(smiles, molname, output_dir)
    if not candidates:
        message = (
            "No linker generator was found in this environment. Install Martini Mapper "
            "or AutoMartini M3 to generate linker .gro/.itp files from SMILES."
        )
        return LinkerGenerationResult(False, None, None, [], message, "", None)

    logs = []
    for generator, command in candidates:
        completed = subprocess.run(
            command,
            cwd=output_dir,
            text=True,
            capture_output=True,
            check=False,
            timeout=180,
        )
        log = "\n".join(part for part in [completed.stdout, completed.stderr] if part).strip()
        logs.append(f"$ {' '.join(command)}\n{log}".strip())
        gro_path = output_dir / f"{molname}.gro"
        itp_path = output_dir / f"{molname}.itp"
        if completed.returncode == 0 and gro_path.exists() and itp_path.exists():
            beads = parse_linker_gro(gro_path)
            if beads:
                return LinkerGenerationResult(True, gro_path, itp_path, beads, "Linker generated.", log, generator)
            return LinkerGenerationResult(False, gro_path, itp_path, [], "Generated linker has no readable beads.", log, generator)

    return LinkerGenerationResult(
        False,
        None,
        None,
        [],
        "The linker generator ran but did not produce the expected .gro and .itp files.",
        "\n\n".join(logs),
        candidates[-1][0],
    )
