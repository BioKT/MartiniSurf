#!/usr/bin/env python3
"""
MartiniSurf Orientation Assistant

Stable Classical Mode
+ Multi-Linker Mode (Physically Correct)
+ Optional Random Surface Linkers (Above Surface Only)
"""

import argparse
import numpy as np
from typing import Tuple, List


def convert_pdb_to_gro(pdb_file: str, gro_file: str) -> None:
    """
    Lightweight PDB->GRO converter kept for backward compatibility with tests
    and legacy standalone usage.
    """
    atoms = []

    with open(pdb_file, "r") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                resid = int(line[22:26])
                resname = line[17:20].strip() or "MOL"
                atomname = line[12:16].strip() or "X"
                x = float(line[30:38]) / 10.0
                y = float(line[38:46]) / 10.0
                z = float(line[46:54]) / 10.0
            except ValueError:
                continue
            atoms.append((resid, resname, atomname, x, y, z))

    with open(gro_file, "w") as fh:
        fh.write("Converted from PDB\n")
        fh.write(f"{len(atoms):5d}\n")
        for i, (resid, resname, atomname, x, y, z) in enumerate(atoms, start=1):
            fh.write(f"{resid:5d}{resname:<5}{atomname:>5}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        fh.write("   5.00000   5.00000   5.00000\n")


# ================================================================
# GRO LOADER
# ================================================================
def load_gro_coords(gro_file: str) -> Tuple[np.ndarray, List[Tuple[int, str, str, int]]]:

    coords = []
    atoms = []

    with open(gro_file, "r") as fh:
        lines = fh.readlines()[2:-1]

    for line in lines:
        try:
            resid    = int(line[0:5])
            resname  = line[5:10].strip()
            atomname = line[10:15].strip()
            atomid   = int(line[15:20])

            x = float(line[20:28]) * 10.0
            y = float(line[28:36]) * 10.0
            z = float(line[36:44]) * 10.0

            atoms.append((resid, resname, atomname, atomid))
            coords.append([x, y, z])
        except:
            continue

    return np.array(coords, float), atoms


# ================================================================
# RESIDUE CENTROIDS
# ================================================================
def summarize_selected_residues(residue_list, atom_records, coords):

    centroids = []

    for resid in residue_list:
        sel = [c for (r,_,_,_), c in zip(atom_records, coords) if r == resid]

        if not sel:
            # Legacy behavior used by tests: skip missing residues.
            continue

        centroids.append(np.mean(sel, axis=0))

    return np.array(centroids)


def _available_resids(atom_records):
    return sorted({int(r) for (r, _, _, _) in atom_records})


def _centroids_or_error(residue_list, atom_records, coords, label):
    centroids = summarize_selected_residues(residue_list, atom_records, coords)
    if centroids.size == 0:
        available = _available_resids(atom_records)
        preview = ", ".join(str(x) for x in available[:12])
        if len(available) > 12:
            preview += ", ..."
        raise ValueError(
            f"No atoms found for {label}. Requested residues: {residue_list}. "
            f"Available residue ids (sample): [{preview}]"
        )
    return centroids


# ================================================================
# CLASSICAL ORIENTATION
# ================================================================
def auto_orient_from_anchor_residues(system_coords,
                                     anchor_centroids,
                                     surface_coords,
                                     target_z):

    if anchor_centroids.size == 0:
        raise ValueError("Anchor centroid selection is empty. Check --anchor/--linker-group residue ids.")

    anchors  = anchor_centroids.copy()
    centroid = anchors.mean(0)

    A = anchors - centroid
    _, _, vh = np.linalg.svd(A)
    normal = vh[2] / np.linalg.norm(vh[2])

    target_normal = np.array([0, 0, -1])

    v = np.cross(normal, target_normal)
    s = np.linalg.norm(v)
    c = float(np.dot(normal, target_normal))

    if s < 1e-6:
        R = np.eye(3)
    else:
        vx = np.array([[0,-v[2],v[1]],
                       [v[2],0,-v[0]],
                       [-v[1],v[0],0]])
        R = np.eye(3) + vx + vx@vx*((1-c)/(s*s))

    center = system_coords.mean(0)

    rot = (R @ (system_coords - center).T).T + center
    anchors_rot = (R @ (anchors - center).T).T + center

    if anchors_rot[:,2].mean() > center[2]:
        Rflip = np.diag([1,1,-1])
        rot = (Rflip @ (rot - center).T).T + center
        # Keep anchor coordinates consistent with the flipped system before
        # computing z-translation.
        anchors_rot = (Rflip @ (anchors_rot - center).T).T + center

    z_anchor  = anchors_rot[:,2].min()
    z_surface = surface_coords[:,2].max()

    rot[:,2] += (z_surface + target_z) - z_anchor

    # Safety guard: anchor-based placement can still leave parts of large/asymmetric
    # proteins below the surface if selected anchors are not the lowest points.
    # Enforce a minimal global clearance to avoid surface penetration artifacts.
    min_clearance = 1.0  # Angstrom
    min_system_z = rot[:, 2].min()
    required_min_z = z_surface + min_clearance
    if min_system_z < required_min_z:
        rot[:, 2] += (required_min_z - min_system_z)

    xy_shift = surface_coords[:,:2].mean(0) - rot[:,:2].mean(0)
    rot[:,0] += xy_shift[0]
    rot[:,1] += xy_shift[1]

    return rot


# ================================================================
# SAVE SYSTEM
# ================================================================
def save_full_system(output_gro,
                     surf_atoms,
                     surf_coords,
                     enz_atoms,
                     enz_coords,
                     box_line):

    total = len(surf_coords) + len(enz_coords)

    with open(output_gro, "w") as fh:

        fh.write("MartiniSurf oriented system\n")
        fh.write(f"{total:5d}\n")

        next_res = 1
        next_atom = 1
        res_map = {}

        for (r,rn,a,_),(x,y,z) in zip(enz_atoms, enz_coords):
            if r not in res_map:
                res_map[r] = next_res
                next_res += 1

            fh.write(f"{res_map[r]:5d}{rn:<5}{a:>5}{next_atom:5d}"
                     f"{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}\n")
            next_atom += 1

        for (r,rn,a,_),(x,y,z) in zip(surf_atoms, surf_coords):
            fh.write(f"{next_res:5d}{rn:<5}{a:>5}{next_atom:5d}"
                     f"{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}\n")
            next_atom += 1

        fh.write(box_line + "\n")


# ================================================================
# MAIN
# ================================================================
def main(argv=None):

    parser = argparse.ArgumentParser()

    parser.add_argument("--surface", required=True)
    parser.add_argument("--system", required=True)
    parser.add_argument("--out", default="system_Surface.gro")

    parser.add_argument("--anchor", nargs="+", action="append")
    parser.add_argument("--dist", type=float, default=10.0)

    parser.add_argument("--linker-gro")
    parser.add_argument("--linker-group", nargs="+", action="append")
    parser.add_argument("--linker-prot-dist", type=float, default=3.0)
    parser.add_argument("--linker-surf-dist", type=float, default=3.0)
    parser.add_argument("--invert-linker", action="store_true")

    parser.add_argument("--surface-linkers", type=int, default=0)
    parser.add_argument("--surface-min-dist", type=float, default=3.0)

    args = parser.parse_args(argv)

    surf_coords, surf_atoms = load_gro_coords(args.surface)
    sys_coords, sys_atoms   = load_gro_coords(args.system)
    box_line = open(args.surface).readlines()[-1].strip()

    # ============================================================
    # CLASSICAL MODE
    # ============================================================
    if not args.linker_gro:

        if not args.anchor:
            raise ValueError("No anchors provided.")

        all_res = sorted({int(r) for g in args.anchor for r in g[1:]})

        centroids = _centroids_or_error(
            all_res,
            sys_atoms,
            sys_coords,
            "anchor groups"
        )

        oriented = auto_orient_from_anchor_residues(
            sys_coords,
            centroids,
            surf_coords,
            args.dist
        )

        final_atoms = sys_atoms


    # ============================================================
    # MULTI-LINKER MODE
    # ============================================================
    else:

        if not args.linker_group:
            raise ValueError("Provide --linker-group.")

        linker_coords, linker_atoms = load_gro_coords(args.linker_gro)
        if args.invert_linker:
            linker_coords = linker_coords[::-1].copy()
            linker_atoms = linker_atoms[::-1]

        all_res = sorted({int(r) for g in args.linker_group for r in g[1:]})

        centroids = _centroids_or_error(
            all_res,
            sys_atoms,
            sys_coords,
            "linker groups"
        )

        oriented_protein = auto_orient_from_anchor_residues(
            sys_coords,
            centroids,
            surf_coords,
            target_z=0.0
        )

        merged_coords = oriented_protein
        merged_atoms  = sys_atoms
        linker_tips = []

        # ONE linker per group
        for group in args.linker_group:

            residues = [int(r) for r in group[1:]]
            if not residues:
                raise ValueError(
                    f"Linker group {group[0]} has no residues. Use: --linker-group <group_id> <resid ...>"
                )

            group_centroid = _centroids_or_error(
                residues,
                sys_atoms,
                oriented_protein,
                f"linker group {group[0]}"
            ).mean(axis=0)

            axis = linker_coords[-1] - linker_coords[0]
            axis_norm = np.linalg.norm(axis)
            if axis_norm < 1e-12:
                raise ValueError(
                    "Linker has zero-length axis (first and last bead overlap). "
                    "Provide a linker with at least two distinct beads."
                )
            axis /= axis_norm

            target = np.array([0,0,-1])

            v = np.cross(axis, target)
            s = np.linalg.norm(v)
            c = np.dot(axis, target)

            if s < 1e-8:
                R = np.eye(3)
            else:
                vx = np.array([[0,-v[2],v[1]],
                               [v[2],0,-v[0]],
                               [-v[1],v[0],0]])
                R = np.eye(3) + vx + vx@vx*((1-c)/(s*s))

            rotated_linker = (R @ (linker_coords - linker_coords[0]).T).T
            rotated_linker += group_centroid + np.array([0,0,-args.linker_prot_dist])

            linker_tips.append(rotated_linker[-1])

            merged_coords = np.vstack([merged_coords, rotated_linker])
            merged_atoms  = merged_atoms + linker_atoms

        # Translate to surface
        z_surface = surf_coords[:,2].max()
        min_tip_z = min([tip[2] for tip in linker_tips])
        dz = (z_surface + args.linker_surf_dist) - min_tip_z
        merged_coords[:,2] += dz

        xy_shift = surf_coords[:,:2].mean(0) - merged_coords[:,:2].mean(0)
        merged_coords[:,0] += xy_shift[0]
        merged_coords[:,1] += xy_shift[1]

        # ---------------------------------------------------------
        # RANDOM SURFACE LINKERS (FIXED ABOVE SURFACE)
        # ---------------------------------------------------------
        if args.surface_linkers > 0:

            xmin, xmax = surf_coords[:,0].min(), surf_coords[:,0].max()
            ymin, ymax = surf_coords[:,1].min(), surf_coords[:,1].max()

            z_surface = surf_coords[:,2].max()

            for _ in range(args.surface_linkers):

                axis = linker_coords[-1] - linker_coords[0]
                axis_norm = np.linalg.norm(axis)
                if axis_norm < 1e-12:
                    raise ValueError(
                        "Linker has zero-length axis (first and last bead overlap). "
                        "Cannot place random surface linkers."
                    )
                axis /= axis_norm

                target = np.array([0,0,1])  # UPWARD

                v = np.cross(axis, target)
                s = np.linalg.norm(v)
                c = np.dot(axis, target)

                if s < 1e-8:
                    R = np.eye(3)
                else:
                    vx = np.array([[0,-v[2],v[1]],
                                   [v[2],0,-v[0]],
                                   [-v[1],v[0],0]])
                    R = np.eye(3) + vx + vx@vx*((1-c)/(s*s))

                vertical_linker = (R @ (linker_coords - linker_coords[0]).T).T

                rand_x = np.random.uniform(xmin, xmax)
                rand_y = np.random.uniform(ymin, ymax)

                # 🔥 ALWAYS ABOVE SURFACE
                vertical_linker += np.array([
                    rand_x,
                    rand_y,
                    z_surface + args.surface_min_dist
                ])

                merged_coords = np.vstack([merged_coords, vertical_linker])
                merged_atoms  = merged_atoms + linker_atoms

        oriented = merged_coords
        final_atoms = merged_atoms


    save_full_system(
        args.out,
        surf_atoms,
        surf_coords,
        final_atoms,
        oriented,
        box_line
    )

    print(f"✔ Saved oriented system → {args.out}")


if __name__ == "__main__":
    main()
