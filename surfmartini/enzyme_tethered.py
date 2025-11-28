#!/usr/bin/env python3
"""
MartiniSurf Orientation Assistant

Tools to:
- Load GROMACS .gro coordinates
- Convert PDB → GRO
- Compute anchor-residue centroids
- Automatically orient CG enzymes above flat surfaces
- Build the final enzyme + surface GRO file

All functions are typed and documented for HPC workflows.
"""

import argparse
from typing import List, Tuple, Sequence

import numpy as np


# ======================================================================
# GRO LOADING
# ======================================================================
def load_gro_coords(gro_file: str) -> Tuple[np.ndarray, List[Tuple[int, str, str, int]]]:
    """
    Load coordinates and atom records from a GROMACS .gro file.

    Parameters
    ----------
    gro_file : str
        Path to the input .gro file.

    Returns
    -------
    tuple
        (coords, atoms) where:
        - coords : np.ndarray (N, 3) in Å
        - atoms  : list of (resid, resname, atomname, atomid)

    Notes
    -----
    Coordinates in .gro are in nm; internally we use Å.
    """
    coords: List[List[float]] = []
    atoms: List[Tuple[int, str, str, int]] = []

    with open(gro_file, "r") as fh:
        lines = fh.readlines()[2:-1]  # skip title + count + last line (box)

    for line in lines:
        try:
            resid = int(line[0:5])
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            atomid = int(line[15:20])

            x = float(line[20:28]) * 10.0  # nm → Å
            y = float(line[28:36]) * 10.0
            z = float(line[36:44]) * 10.0

            atoms.append((resid, resname, atomname, atomid))
            coords.append([x, y, z])

        except ValueError:
            continue

    return np.array(coords, dtype=float), atoms


# ======================================================================
# PDB → GRO CONVERSION
# ======================================================================
def convert_pdb_to_gro(pdb_file: str, gro_file: str) -> str:
    """
    Convert a PDB file into a simplified .gro file.

    Parameters
    ----------
    pdb_file : str
        Path to the input .pdb file.
    gro_file : str
        Output path for the resulting .gro file.

    Returns
    -------
    str
        The path to the written GRO file.
    """
    atoms: List[Tuple[int, str, str, int]] = []
    coords: List[List[float]] = []

    with open(pdb_file, "r") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atomname = line[12:16].strip()
                resname = line[17:20].strip()
                resid = int(line[22:26])

                x = float(line[30:38]) / 10.0
                y = float(line[38:46]) / 10.0
                z = float(line[46:54]) / 10.0

                atoms.append((resid, resname, atomname, len(atoms) + 1))
                coords.append([x, y, z])

    coords_np = np.array(coords, dtype=float)

    with open(gro_file, "w") as fh:
        fh.write("Converted from PDB by MartiniSurf\n")
        fh.write(f"{len(coords_np):5d}\n")

        resid_map: dict[int, int] = {}
        new_resid = 0

        for (resid, resname, atom, atomid), (x, y, z) in zip(atoms, coords_np):
            if resid not in resid_map:
                new_resid += 1
                resid_map[resid] = new_resid

            fh.write(
                f"{resid_map[resid]:5d}{resname:<5}{atom:>5}{atomid:5d}"
                f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
            )

        fh.write("   10.00000   10.00000   10.00000\n")

    return gro_file


# ======================================================================
# SELECT RESIDUES / CENTROIDS
# ======================================================================
def summarize_selected_residues(
    residue_list: List[int],
    atom_records: List[Tuple[int, str, str, int]],
    coords: np.ndarray,
) -> np.ndarray:
    """
    Report and calculate centroid coordinates for selected residues.
    """
    centroids: List[np.ndarray] = []
    messages: List[str] = []

    for resid in residue_list:
        sel = [
            (r, rn, a, aid, c)
            for (r, rn, a, aid), c in zip(atom_records, coords)
            if r == resid
        ]

        if not sel:
            messages.append(f"⚠ Residue {resid} not found.")
            continue

        loc_coords = np.array([c for *_, c in sel], dtype=float)
        centroids.append(loc_coords.mean(axis=0))
        messages.append(f"• Residue {resid:4d}")

    print("\n=== Selected Anchor Residues ===")
    print("\n".join(messages))
    print("================================\n")

    return np.array(centroids, dtype=float)


# ======================================================================
# ORIENTATION ENGINE
# ======================================================================
def auto_orient_from_anchor_residues(
    enzyme_coords: np.ndarray,
    anchor_centroids: np.ndarray,
    surface_coords: np.ndarray,
    target_z: float,
) -> np.ndarray:
    """
    Orient the enzyme so that its anchor residues face the surface.
    """
    anchors = anchor_centroids.copy()
    centroid = anchors.mean(axis=0)

    A = anchors - centroid
    _, _, vh = np.linalg.svd(A)
    normal = vh[2] / np.linalg.norm(vh[2])

    target_normal = np.array([0.0, 0.0, -1.0])
    v = np.cross(normal, target_normal)
    s = np.linalg.norm(v)
    c = float(np.dot(normal, target_normal))

    if s < 1e-6:
        R = np.eye(3)
    else:
        vx = np.array(
            [[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], dtype=float
        )
        R = np.eye(3) + vx + vx @ vx * ((1 - c) / (s**2))

    center = enzyme_coords.mean(axis=0)
    rot = (R @ (enzyme_coords - center).T).T + center
    anchors_rot = (R @ (anchors - center).T).T + center

    if anchors_rot[:, 2].min() > rot[:, 2].mean():
        Rflip = np.diag([1.0, 1.0, -1.0])
        rot = (Rflip @ (rot - center).T).T + center

    z_anchor = anchors_rot[:, 2].min()
    z_surface = surface_coords[:, 2].max()
    dz = (z_surface + target_z) - z_anchor
    rot[:, 2] += dz

    xy_shift = surface_coords[:, :2].mean(axis=0) - rot[:, :2].mean(axis=0)
    rot[:, 0] += xy_shift[0]
    rot[:, 1] += xy_shift[1]

    return rot


# ======================================================================
# SAVE SYSTEM
# ======================================================================
def save_full_system(
    output_gro: str,
    surf_atoms: List[Tuple[int, str, str, int]],
    surf_coords: np.ndarray,
    enz_atoms: List[Tuple[int, str, str, int]],
    enz_coords: np.ndarray,
) -> None:
    """
    Write the final enzyme+surface system to a GRO file.
    """
    total = len(surf_coords) + len(enz_coords)

    with open(output_gro, "w") as fh:
        fh.write("MartiniSurf oriented system\n")
        fh.write(f"{total:5d}\n")

        next_resid = 1
        next_atomID = 1
        enz_map: dict[int, int] = {}

        for (r, rn, a, aid), (x, y, z) in zip(enz_atoms, enz_coords):
            if r not in enz_map:
                enz_map[r] = next_resid
                next_resid += 1

            fh.write(
                f"{enz_map[r]:5d}{rn:<5}{a:>5}{next_atomID:5d}"
                f"{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}\n"
            )
            next_atomID += 1

        surf_map: dict[int, int] = {}

        for (r, rn, a, aid), (x, y, z) in zip(surf_atoms, surf_coords):
            if r not in surf_map:
                surf_map[r] = next_resid
                next_resid += 1

            fh.write(
                f"{surf_map[r]:5d}{rn:<5}{a:>5}{next_atomID:5d}"
                f"{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}\n"
            )
            next_atomID += 1

        x_box = np.ptp(surf_coords[:, 0]) / 10.0
        y_box = np.ptp(surf_coords[:, 1]) / 10.0
        z_box = np.ptp(enz_coords[:, 2]) / 10.0 + 10.0

        fh.write(f"{x_box:12.5f}{y_box:12.5f}{z_box:12.5f}\n")

    print(f"✔ Saved oriented system → {output_gro}")


# ======================================================================
# MAIN PROGRAM
# ======================================================================
def main(argv: Sequence[str] | None = None) -> None:
    """
    Command-line interface for the MartiniSurf Orientation Assistant.
    """
    print(
        "\n==============================================================\n"
        "               MartiniSurf – Orientation Assistant\n"
        "==============================================================\n"
    )

    parser = argparse.ArgumentParser(
        description="Automatic CG-enzyme orientation on a flat surface."
    )

    parser.add_argument("--surface", required=True)
    parser.add_argument("--enzyme", required=True)
    parser.add_argument("--out", default="Enzyme_Surface.gro")

    parser.add_argument("--resA", nargs="+", type=int)
    parser.add_argument("--resB", nargs="+", type=int)
    parser.add_argument("--anchor", nargs="+", type=int)

    parser.add_argument("--dist", type=float, default=10.0)
    parser.add_argument("--display", choices=["on", "off"], default="off")

    args = parser.parse_args(argv)

    # Anchor selection
    if args.anchor:
        anchors = args.anchor
        print(f"• Using anchor list: {anchors}")
    else:
        anchors = []
        if args.resA:
            anchors.extend(args.resA)
        if args.resB:
            anchors.extend(args.resB)

        if not anchors:
            raise ValueError("You must supply --anchor OR (--resA/--resB).")

        print(f"• Using combined anchors: {anchors}")

    # Auto PDB → GRO
    if args.enzyme.lower().endswith(".pdb"):
        new_gro = args.enzyme.replace(".pdb", ".gro")
        print(f"⚠ Converting PDB → GRO: {new_gro}")
        convert_pdb_to_gro(args.enzyme, new_gro)
        args.enzyme = new_gro

    surf_coords, surf_atoms = load_gro_coords(args.surface)
    enz_coords, enz_atoms = load_gro_coords(args.enzyme)

    print(f"• Loaded surface: {len(surf_coords)} atoms")
    print(f"• Loaded enzyme : {len(enz_coords)} atoms")

    centroids = summarize_selected_residues(anchors, enz_atoms, enz_coords)

    oriented = auto_orient_from_anchor_residues(
        enz_coords, centroids, surf_coords, args.dist
    )

    save_full_system(args.out, surf_atoms, surf_coords, enz_atoms, oriented)


# ======================================================================
if __name__ == "__main__":
    main()
