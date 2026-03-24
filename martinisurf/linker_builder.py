from pathlib import Path
import shutil
import numpy as np

def is_smiles(value: str) -> bool:
    return not Path(value).exists()

def load_gro_coords(gro_file):
    coords = []
    atoms = []
    with open(gro_file, "r") as fh:
        lines = fh.readlines()[2:-1]

    for line in lines:
        resid = int(line[0:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomid = int(line[15:20])
        x = float(line[20:28]) * 10.0
        y = float(line[28:36]) * 10.0
        z = float(line[36:44]) * 10.0
        atoms.append((resid, resname, atomname, atomid))
        coords.append([x, y, z])

    return np.array(coords, float), atoms


def build_linker(linker_input: str,
                 molname: str,
                 is_dna: bool,
                 workdir: Path):
    # Reserved for future DNA-specific linker handling; keep in signature for API stability.
    _ = is_dna

    workdir = Path(workdir)
    linker_gro = workdir / f"{molname}.gro"

    if is_smiles(linker_input):
        raise NotImplementedError(
            "SMILES linker generation not implemented yet."
        )

    else:
        src = Path(linker_input)
        if not src.exists():
            raise FileNotFoundError(f"Linker file not found: {src}")

        shutil.copy(src, linker_gro)

    coords, atoms = load_gro_coords(linker_gro)

    return linker_gro, coords, atoms
