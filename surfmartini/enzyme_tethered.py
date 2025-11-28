#!/usr/bin/env python3
"""
MartiniSurf Orientation Assistant

This module provides tools to:
- Load GROMACS .gro coordinates
- Convert PDB → GRO
- Compute anchor-residue centroids
- Automatically orient CG enzymes on flat surfaces
- Assemble the final enzyme + surface system into a .gro file

All functions are documented for scientific and HPC workflows.
"""

import argparse
from typing import Iterable, List, Tuple

import numpy as np


# ======================================================================
# GRO LOADING
# ======================================================================
def load_gro_coords(gro_file: str) -> Tuple[np.ndarray, List[Tuple]]:
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
    Coordinates from .gro are in nm; internally we store positions in Å.
    """
    coords, atoms = [], []

    with open(gro_file) as fh:
        lines = fh.readlines()[2:-1]  # skip title + atom count

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
        except Exception:
            continue

    return np.array(coords), atoms


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
        Output path for the generated .gro file.

    Returns
    -------
    str
        Path to the written GRO file.
    """
    atoms = []
    coords = []

    with open(pdb_file) as fh:
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

    coords = np.array(coords)

    with open(gro_file, "w") as fh:
        fh.write("Converted from PDB by MartiniSurf\n")
        fh.write(f"{len(coords):5d}\n")

        resid_map = {}
        new_resid = 0

        for (resid, resname, atom, atomid), (x, y, z) in zip(atoms, coords):
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
# RESIDUE SUMMARY / ANCHOR CENTROIDS
# ======================================================================
def summarize_selected_residues(
    residue_list: List[int],
    atom_records: List[Tuple],
    coords: np.ndarray,
) -> np.ndarray:
    """
    Report and calculate centroid coordinates for selected residues.

    Parameters
    ----------
    residue_list : list of int
        Residue IDs to highlight and summarize.
    atom_records : list
        Atom metadata from load_gro_coords().
    coords : np.ndarray
        Cartesian coordinates for each atom (N, 3) in Å.

    Returns
    -------
    np.ndarray
        Array of centroid coordinates (M, 3) in Å.
    """
    summary_lines = []
    centroids = []

    for resid in residue_list:
        sel_atoms = [
            (r, resname, atom, aid, c)
            for (r, resname, atom, aid), c in zip(atom_records, coords)
            if r == resid
        ]

        if not sel_atoms:
            summary_lines.append(f"⚠ Residue {resid} not found.")
            continue

        loc_coords = np.array([c for *_, c in sel_atoms])
        centroids.append(loc_coords.mean(axis=0))

        summary_lines.append(f"• Residue {resid:4d}")

    print("\n=== Selected Anchor Residues ===")
    print("\n".join(summary_lines))
    print("================================\n")

    return np.array(centroids)


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
    Orient the enzyme so that anchor residues face the surface.

    Parameters
    ----------
    enzyme_coords : np.ndarray
        Coordinates of the enzyme atoms (Å).
    anchor_centroids : np.ndarray
        Centroids of selected anchor residues (Å).
    surface_coords : np.ndarray
        Coordinates of surface atoms (Å).
    target_z : float
        Vertical offset above the surface, in Å.

    Returns
    -------
    np.ndarray
        Rotated and translated enzyme coordinates (Å).
    """
    anchors = anchor_centroids.copy()
    centroid = anchors.mean(axis=0)

    # PCA-like normal estimation
    A = anchors - centroid
    _, _, vh = np.linalg.svd(A)
    normal = vh[2] / np.linalg.norm(vh[2])

    target_normal = np.array([0, 0, -1])
    v = np.cross(normal, target_normal)
    s = np.linalg.norm(v)
    c = np.dot(normal, target_normal)

    # Rodrigues rotation
    if s < 1e-6:
        R = np.eye(3)
    else:
        vx = np.array(
            [
                [0, -v[2], v[1]],
                [v[2], 0, -v[0]],
                [-v[1], v[0], 0],
            ]
        )
        R = np.eye(3) + vx + vx.dot(vx) * ((1 - c) / s**2)

    # Rotate around enzyme center
    center = enzyme_coords.mean(axis=0)
    rot = (R @ (enzyme_coords - center).T).T + center
    anchors_rot = (R @ (anchors - center).T).T + center

    # Flip check
    if anchors_rot[:, 2].min() > rot[:, 2].mean():
        Rflip = np.diag([1, 1, -1])
        rot = (Rflip @ (rot - center).T).T + center
        anchors_rot = (Rflip @ (anchors_rot - center).T).T + center

    # Align Z position
    z_anchor = anchors_rot[:, 2].min()
    z_surface = surface_coords[:, 2].max()
    dz = (z_surface + target_z) - z_anchor
    rot[:, 2] += dz

    # Center XY alignment
    xy_shift = surface_coords[:, :2].mean(axis=0) - rot[:, :2].mean(axis=0)
    rot[:, 0] += xy_shift[0]
    rot[:, 1] += xy_shift[1]

    return rot


# ======================================================================
# SAVE ORIENTED GRO
# ======================================================================
def save_full_system(
    output_gro: str,
    surf_atoms: List[Tuple],
    surf_coords: np.ndarray,
    enz_atoms: List[Tuple],
    enz_coords: np.ndarray,
) -> None:
    """
    Write the final enzyme + surface system to a GRO file.

    Parameters
    ----------
    output_gro : str
        Path to output .gro file.
    surf_atoms : list
        Surface atom metadata.
    surf_coords : np.ndarray
        Surface atom coordinates (Å).
    enz_atoms : list
        Enzyme atom metadata.
    enz_coords : np.ndarray
        Enzyme atom coordinates (Å).
    """
    total = len(surf_coords) + len(enz_coords)

    with open(output_gro, "w") as fh:
        fh.write("MartiniSurf oriented system\n")
        fh.write(f"{total:5d}\n")

        next_resid = 1
        next_atomID = 1
        enz_map = {}

        # Enzyme
        for (r, rn, a, aid), (x, y, z) in zip(enz_atoms, enz_coords):
            if r not in enz_map:
                enz_map[r] = next_resid
                next_resid += 1

            fh.write(
                f"{enz_map[r]:5d}{rn:<5}{a:>5}{next_atomID:5d}"
                f"{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}\n"
            )
            next_atomID += 1

        # Surface
        surf_map = {}
        for (r, rn, a, aid), (x, y, z) in zip(surf_atoms, surf_coords):
            if r not in surf_map:
                surf_map[r] = next_resid
                next_resid += 1

            fh.write(
                f"{surf_map[r]:5d}{rn:<5}{a:>5}{next_atomID:5d}"
                f"{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}\n"
            )
            next_atomID += 1

        # Box dimensions (NumPy 2.0+ compatible)
        x_box = np.ptp(surf_coords[:, 0]) / 10.0
        y_box = np.ptp(surf_coords[:, 1]) / 10.0
        z_box = np.ptp(enz_coords[:, 2]) / 10.0 + 10.0

        fh.write(f"{x_box:12.5f}{y_box:12.5f}{z_box:12.5f}\n")

    print(f"✔ Saved oriented system → {output_gro}")


# ======================================================================
# MAIN PROGRAM
# ======================================================================
def main(argv: Iterable[str] | None = None) -> None:
    """
    Command-line interface for the MartiniSurf Orientation Assistant.

    Parameters
    ----------
    argv : iterable of str or None
        Input arguments. Default: sys.argv.
    """
    banner = r"""
==============================================================
               MartiniSurf – Orientation Assistant
==============================================================
"""
    print(banner)

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

    # ===== Determine anchor residues =====
    if args.anchor:
        anchors = args.anchor
        print(f"• Using anchor list: {anchors}")

    elif args.resA or args.resB:
        anchors = []
        if args.resA:
            anchors.extend(args.resA)
        if args.resB:
            anchors.extend(args.resB)
        print(f"• Using combined anchors: {anchors}")

        if not anchors:
            raise ValueError("resA/resB provided but no residues found.")
    else:
        raise ValueError("Missing --anchor OR ( --resA / --resB ).")

    # ===== Auto-convert PDB → GRO =====
    if args.enzyme.lower().endswith(".pdb"):
        gro_out = args.enzyme.replace(".pdb", ".gro")
        print(f"⚠ Converting PDB → GRO: {gro_out}")
        convert_pdb_to_gro(args.enzyme, gro_out)
        args.enzyme = gro_out

    # ===== Load structures =====
    surf_coords, surf_atoms = load_gro_coords(args.surface)
    enz_coords, enz_atoms = load_gro_coords(args.enzyme)

    print(f"• Loaded surface: {len(surf_coords)} atoms")
    print(f"• Loaded enzyme : {len(enz_coords)} atoms")

    # ===== Anchor centroids =====
    centroids = summarize_selected_residues(anchors, enz_atoms, enz_coords)

    # ===== Orientation =====
    oriented = auto_orient_from_anchor_residues(
        enz_coords, centroids, surf_coords, args.dist
    )

    # ===== Save final system =====
    save_full_system(args.out, surf_atoms, surf_coords, enz_atoms, oriented)


# ======================================================================
if __name__ == "__main__":
    main()
