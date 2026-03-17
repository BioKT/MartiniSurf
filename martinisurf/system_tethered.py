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


def _parse_group_residues(raw_groups, option_name):
    parsed = []
    for group in raw_groups or []:
        if len(group) < 2:
            raise ValueError(
                f"{option_name} entries must include a group id and at least one residue."
            )
        try:
            residues = [int(r) for r in group[1:]]
        except ValueError as exc:
            raise ValueError(
                f"{option_name} residues must be integers. Received: {group}"
            ) from exc
        parsed.append((str(group[0]), residues))
    return parsed


def _candidate_pairs_or_all(centroids: np.ndarray) -> list[np.ndarray]:
    n = len(centroids)
    if n <= 2:
        return [np.array(centroids, float, copy=True)]
    out = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            out.append(np.array([centroids[i], centroids[j]], float))
    return out


def _select_best_two_group_z_subsets(group_a: np.ndarray, group_b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Pick one subset per anchor group that is globally most consistent in Z.
    This is used only in pre-config mode to ignore outlier residues during
    orientation while keeping all original residues for pull groups/topology.
    """
    cand_a = _candidate_pairs_or_all(group_a)
    cand_b = _candidate_pairs_or_all(group_b)
    best = None
    best_score = None
    for a in cand_a:
        for b in cand_b:
            all_z = np.concatenate([a[:, 2], b[:, 2]])
            score = (
                float(np.std(all_z)),
                float(abs(float(np.mean(a[:, 2])) - float(np.mean(b[:, 2])))),
                float(np.std(a[:, 2])) if len(a) > 1 else 0.0,
                float(np.std(b[:, 2])) if len(b) > 1 else 0.0,
            )
            if best_score is None or score < best_score:
                best_score = score
                best = (a, b)
    if best is None:
        return group_a, group_b
    return best


def _anchor_landmarks_for_groups(group_defs, atom_records, coords, mode, prefilter_z_consistency=False):
    residue_centroids = []
    for group_id, residues in group_defs:
        residue_centroids.append(
            _centroids_or_error(
                residues,
                atom_records,
                coords,
                f"anchor group {group_id}",
            )
        )

    if prefilter_z_consistency and mode == "residue" and len(residue_centroids) == 2:
        a, b = _select_best_two_group_z_subsets(residue_centroids[0], residue_centroids[1])
        residue_centroids = [a, b]

    if mode == "group":
        return np.array([centroids.mean(axis=0) for centroids in residue_centroids], float)

    return np.vstack(residue_centroids)


def _rotation_matrix_from_vectors(source, target):
    src = np.array(source, float)
    dst = np.array(target, float)
    src /= np.linalg.norm(src)
    dst /= np.linalg.norm(dst)

    cross = np.cross(src, dst)
    s = np.linalg.norm(cross)
    c = float(np.dot(src, dst))

    if s < 1e-8:
        if c > 0:
            return np.eye(3)
        ortho = np.array([1.0, 0.0, 0.0])
        if abs(src[0]) > 0.9:
            ortho = np.array([0.0, 1.0, 0.0])
        axis = np.cross(src, ortho)
        axis /= np.linalg.norm(axis)
        x, y, z = axis
        K = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
        return np.eye(3) + 2.0 * (K @ K)

    vx = np.array([[0, -cross[2], cross[1]],
                   [cross[2], 0, -cross[0]],
                   [-cross[1], cross[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))


def _rotation_matrix_about_axis(axis, angle):
    axis = np.array(axis, float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    c = np.cos(angle)
    s = np.sin(angle)
    C = 1.0 - c
    return np.array([
        [c + x*x*C, x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s, c + y*y*C, y*z*C - x*s],
        [z*x*C - y*s, z*y*C + x*s, c + z*z*C],
    ])


def _apply_rotation(coords, rotation, center):
    return (rotation @ (coords - center).T).T + center


def _finalize_anchor_pose(
    system_coords,
    anchor_coords,
    surface_coords,
    target_z,
    reference_coords,
    min_reference_dist=1.0,
):
    z_surface = surface_coords[:, 2].max()
    target_anchor_z = z_surface + target_z

    rot = np.array(system_coords, float, copy=True)
    anchors_rot = np.array(anchor_coords, float, copy=True)
    ref_rot = np.array(reference_coords, float, copy=True)

    dz = target_anchor_z - anchors_rot[:, 2].mean()
    rot[:, 2] += dz
    anchors_rot[:, 2] += dz
    ref_rot[:, 2] += dz

    required_min_z = z_surface + float(min_reference_dist)
    min_ref_z = ref_rot[:, 2].min()
    if min_ref_z < required_min_z:
        shift = required_min_z - min_ref_z
        rot[:, 2] += shift
        anchors_rot[:, 2] += shift
        ref_rot[:, 2] += shift

    xy_shift = surface_coords[:, :2].mean(0) - ref_rot[:, :2].mean(0)
    rot[:, 0] += xy_shift[0]
    rot[:, 1] += xy_shift[1]
    anchors_rot[:, 0] += xy_shift[0]
    anchors_rot[:, 1] += xy_shift[1]
    ref_rot[:, 0] += xy_shift[0]
    ref_rot[:, 1] += xy_shift[1]

    return rot, anchors_rot, ref_rot


def _lowest_z_subset(coords, fraction):
    if len(coords) == 0:
        return np.array(coords, float, copy=True)
    frac = float(fraction)
    count = int(np.ceil(len(coords) * frac))
    count = max(1, min(len(coords), count))
    if count >= len(coords):
        return np.array(coords, float, copy=True)
    idx = np.argpartition(coords[:, 2], count - 1)[:count]
    return coords[idx]


def _score_two_anchor_pose(reference_coords, balance_low_z=False, low_z_fraction=0.2):
    if not balance_low_z or len(reference_coords) < 3:
        return (float(reference_coords[:, 2].mean()),)

    low_z = _lowest_z_subset(reference_coords, low_z_fraction)
    # Prefer flatter low-Z subsets (more parallel to the surface),
    # then keep them as close as possible to the target surface distance.
    return (float(np.std(low_z[:, 2])), float(low_z[:, 2].mean()))


# ================================================================
# CLASSICAL ORIENTATION
# ================================================================
def auto_orient_from_anchor_residues(system_coords,
                                     anchor_centroids,
                                     surface_coords,
                                     target_z,
                                     reference_coords=None,
                                     min_reference_dist=1.0,
                                     balance_low_z=False,
                                     balance_low_z_fraction=0.2,
                                     orient_single_anchor_up=False):

    if anchor_centroids.size == 0:
        raise ValueError("Anchor centroid selection is empty. Check --anchor/--linker-group residue ids.")
    if balance_low_z and not (0.0 < float(balance_low_z_fraction) <= 1.0):
        raise ValueError("--balance-low-z-fraction must be in the interval (0, 1].")

    anchors = np.array(anchor_centroids, float, copy=True)
    ref = np.array(reference_coords if reference_coords is not None else system_coords, float, copy=True)

    if len(anchors) == 1:
        if orient_single_anchor_up:
            center = anchors[0]
            ref_vec = ref.mean(0) - center
            if np.linalg.norm(ref_vec) < 1e-8:
                ref_vec = np.array([0.0, 0.0, 1.0])
            target_vec = np.array([0.0, 0.0, 1.0])
            R_up = _rotation_matrix_from_vectors(ref_vec, target_vec)
            system_coords = _apply_rotation(system_coords, R_up, center)
            anchors = _apply_rotation(anchors, R_up, center)
            ref = _apply_rotation(ref, R_up, center)
        return _finalize_anchor_pose(
            system_coords,
            anchors,
            surface_coords,
            target_z,
            ref,
            min_reference_dist=min_reference_dist,
        )[0]

    if len(anchors) == 2:
        center = anchors.mean(0)
        axis = anchors[1] - anchors[0]
        axis_norm = np.linalg.norm(axis)
        if axis_norm < 1e-8:
            raise ValueError("Anchor groups collapse onto the same point; cannot orient the system.")

        target_axis = axis.copy()
        target_axis[2] = 0.0
        if np.linalg.norm(target_axis) < 1e-8:
            target_axis = np.array([1.0, 0.0, 0.0])

        R_align = _rotation_matrix_from_vectors(axis, target_axis)
        base_system = _apply_rotation(system_coords, R_align, center)
        base_anchors = _apply_rotation(anchors, R_align, center)
        base_ref = _apply_rotation(ref, R_align, center)

        roll_axis = base_anchors[1] - base_anchors[0]
        roll_axis /= np.linalg.norm(roll_axis)
        roll_center = base_anchors.mean(0)

        best_pose = None
        best_anchors = None
        best_ref = None
        best_score = None
        for angle in np.linspace(0.0, 2.0 * np.pi, 181, endpoint=False):
            R_roll = _rotation_matrix_about_axis(roll_axis, angle)
            trial_system = _apply_rotation(base_system, R_roll, roll_center)
            trial_anchors = _apply_rotation(base_anchors, R_roll, roll_center)
            trial_ref = _apply_rotation(base_ref, R_roll, roll_center)
            final_system, final_anchors, final_ref = _finalize_anchor_pose(
                trial_system,
                trial_anchors,
                surface_coords,
                target_z,
                trial_ref,
                min_reference_dist=min_reference_dist,
            )
            score = _score_two_anchor_pose(
                final_ref,
                balance_low_z=balance_low_z,
                low_z_fraction=balance_low_z_fraction,
            )
            if best_score is None or score < best_score:
                best_score = score
                best_pose = final_system
                best_anchors = final_anchors
                best_ref = final_ref

        # Auto-fix mirror solutions in balanced two-anchor mode:
        # if the selected orientation leaves the protein body below its anchors,
        # mirror along Z and re-apply anchor/surface constraints.
        if (
            balance_low_z
            and best_pose is not None
            and best_anchors is not None
            and best_ref is not None
            and float(best_ref[:, 2].mean()) < float(best_anchors[:, 2].mean())
        ):
            flip_center = best_anchors.mean(0)
            Rflip = np.diag([1.0, 1.0, -1.0])
            flip_system = _apply_rotation(best_pose, Rflip, flip_center)
            flip_anchors = _apply_rotation(best_anchors, Rflip, flip_center)
            flip_ref = _apply_rotation(best_ref, Rflip, flip_center)
            fixed_system, fixed_anchors, fixed_ref = _finalize_anchor_pose(
                flip_system,
                flip_anchors,
                surface_coords,
                target_z,
                flip_ref,
                min_reference_dist=min_reference_dist,
            )
            flip_score = _score_two_anchor_pose(
                fixed_ref,
                balance_low_z=balance_low_z,
                low_z_fraction=balance_low_z_fraction,
            )
            if flip_score < best_score:
                return fixed_system

        return best_pose

    centroid = anchors.mean(0)
    A = anchors - centroid
    _, _, vh = np.linalg.svd(A)
    normal = vh[2] / np.linalg.norm(vh[2])
    target_normal = np.array([0.0, 0.0, -1.0])
    center = ref.mean(0)

    R = _rotation_matrix_from_vectors(normal, target_normal)
    rot = _apply_rotation(system_coords, R, center)
    anchors_rot = _apply_rotation(anchors, R, center)
    ref_rot = _apply_rotation(ref, R, center)

    if ref_rot[:, 2].mean() < anchors_rot[:, 2].mean():
        Rflip = np.diag([1.0, 1.0, -1.0])
        rot = _apply_rotation(rot, Rflip, center)
        anchors_rot = _apply_rotation(anchors_rot, Rflip, center)
        ref_rot = _apply_rotation(ref_rot, Rflip, center)

    return _finalize_anchor_pose(
        rot,
        anchors_rot,
        surface_coords,
        target_z,
        ref_rot,
        min_reference_dist=min_reference_dist,
    )[0]


# ================================================================
# SAVE SYSTEM
# ================================================================
def save_full_system(output_gro,
                     surf_atoms,
                     surf_coords,
                     enz_atoms,
                     enz_coords,
                     box_line,
                     dna_mode: bool = False):
    surf_out = np.array(surf_coords, dtype=float, copy=True)
    enz_out = np.array(enz_coords, dtype=float, copy=True)
    total = len(surf_out) + len(enz_out)

    # Keep original XY when possible, but expand/shift when needed so all coordinates
    # fit inside the declared box and avoid PBC-splitting artifacts downstream.
    # GRO coordinates are written in nm; internal coords here are in Angstrom.
    box_tokens = box_line.split()
    if len(box_tokens) >= 3:
        try:
            box_x = float(box_tokens[0])
            box_y = float(box_tokens[1])
            box_z = float(box_tokens[2])
        except ValueError:
            box_x = box_y = box_z = 0.0
    else:
        box_x = box_y = box_z = 0.0

    xy_parts = []
    if len(surf_out) > 0:
        xy_parts.append(surf_out[:, :2])
    if len(enz_out) > 0:
        xy_parts.append(enz_out[:, :2])
    if xy_parts:
        xy = np.vstack(xy_parts)
        min_x = float(np.min(xy[:, 0]))
        min_y = float(np.min(xy[:, 1]))
        shift_x = -min_x if min_x < 0 else 0.0
        shift_y = -min_y if min_y < 0 else 0.0
        if shift_x or shift_y:
            if len(surf_out) > 0:
                surf_out[:, 0] += shift_x
                surf_out[:, 1] += shift_y
            if len(enz_out) > 0:
                enz_out[:, 0] += shift_x
                enz_out[:, 1] += shift_y
        box_x = max(box_x, float(np.max(xy[:, 0] + shift_x)) / 10.0)
        box_y = max(box_y, float(np.max(xy[:, 1] + shift_y)) / 10.0)

    if len(enz_out) > 0:
        if dna_mode:
            dna_resnames = {"DA", "DT", "DG", "DC"}
            dna_idx = [i for i, (_, rn, _, _) in enumerate(enz_atoms) if str(rn).strip() in dna_resnames]
            if dna_idx:
                z_req_nm = float(np.max(enz_out[dna_idx, 2])) / 10.0 + 3.0
            else:
                z_req_nm = float(np.max(enz_out[:, 2])) / 10.0 + 3.0
        else:
            z_req_nm = float(np.max(enz_out[:, 2])) / 10.0 + 3.0
    else:
        z_req_nm = 3.0
    if len(surf_out) > 0:
        z_req_nm = max(z_req_nm, float(np.max(surf_out[:, 2])) / 10.0 + 0.5)
    if dna_mode:
        box_z = z_req_nm
    else:
        box_z = max(box_z, z_req_nm)
    final_box_line = f"{box_x:10.5f}{box_y:10.5f}{box_z:10.5f}"

    with open(output_gro, "w") as fh:

        fh.write("MartiniSurf oriented system\n")
        fh.write(f"{total:5d}\n")

        next_res = 1
        next_atom = 1
        res_map = {}

        for (r,rn,a,_),(x,y,z) in zip(enz_atoms, enz_out):
            if r not in res_map:
                res_map[r] = next_res
                next_res += 1

            fh.write(f"{res_map[r]:5d}{rn:<5}{a:>5}{next_atom:5d}"
                     f"{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}\n")
            next_atom += 1

        for (r,rn,a,_),(x,y,z) in zip(surf_atoms, surf_out):
            fh.write(f"{next_res:5d}{rn:<5}{a:>5}{next_atom:5d}"
                     f"{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}\n")
            next_atom += 1

        fh.write(final_box_line + "\n")


# ================================================================
# MAIN
# ================================================================
def main(argv=None):

    parser = argparse.ArgumentParser()

    parser.add_argument("--surface", required=True)
    parser.add_argument("--system", required=True)
    parser.add_argument("--out", default="system_Surface.gro")

    parser.add_argument("--anchor", nargs="+", action="append")
    parser.add_argument(
        "--anchor-landmark-mode",
        choices=["residue", "group"],
        default="residue",
    )
    parser.add_argument(
        "--prefilter-anchor-z-consistency",
        action="store_true",
        help="Pre-config mode helper: pick the most Z-consistent anchor residues for orientation only.",
    )
    parser.add_argument("--dist", type=float, default=10.0)
    parser.add_argument("--reference-exclude-resname", action="append")

    parser.add_argument("--linker-gro")
    parser.add_argument("--linker-group", nargs="+", action="append")
    parser.add_argument("--linker-prot-dist", type=float, default=3.0)
    parser.add_argument("--linker-surf-dist", type=float, default=3.0)
    parser.add_argument("--invert-linker", action="store_true")

    parser.add_argument("--surface-linkers", type=int, default=0)
    parser.add_argument("--surface-min-dist", type=float, default=3.0)
    parser.add_argument("--dna-mode", action="store_true")
    parser.add_argument("--min-reference-z-dist", type=float, default=1.0)
    parser.add_argument(
        "--balance-low-z",
        action="store_true",
        help="In two-anchor mode, pick the roll angle that flattens the lowest-Z reference region.",
    )
    parser.add_argument(
        "--balance-low-z-fraction",
        type=float,
        default=0.2,
        help="Fraction (0,1] of lowest-Z reference beads used to score flatness in --balance-low-z mode.",
    )

    args = parser.parse_args(argv)
    if not (0.0 < args.balance_low_z_fraction <= 1.0):
        raise ValueError("--balance-low-z-fraction must be in the interval (0, 1].")

    surf_coords, surf_atoms = load_gro_coords(args.surface)
    sys_coords, sys_atoms   = load_gro_coords(args.system)
    box_line = open(args.surface).readlines()[-1].strip()
    reference_coords = sys_coords
    if args.reference_exclude_resname:
        excluded = {name.strip() for name in args.reference_exclude_resname if name.strip()}
        mask = np.array([resname not in excluded for (_, resname, _, _) in sys_atoms], dtype=bool)
        if mask.any():
            reference_coords = sys_coords[mask]

    # ============================================================
    # CLASSICAL MODE
    # ============================================================
    if not args.linker_gro:

        if not args.anchor:
            raise ValueError("No anchors provided.")

        anchor_groups = _parse_group_residues(args.anchor, "--anchor")
        centroids = _anchor_landmarks_for_groups(
            anchor_groups,
            sys_atoms,
            sys_coords,
            args.anchor_landmark_mode,
            prefilter_z_consistency=args.prefilter_anchor_z_consistency,
        )

        oriented = auto_orient_from_anchor_residues(
            sys_coords,
            centroids,
            surf_coords,
            args.dist,
            reference_coords=reference_coords,
            min_reference_dist=args.min_reference_z_dist,
            balance_low_z=args.balance_low_z,
            balance_low_z_fraction=args.balance_low_z_fraction,
            orient_single_anchor_up=(len(centroids) == 1),
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
            target_z=0.0,
            reference_coords=reference_coords,
            min_reference_dist=args.min_reference_z_dist,
            balance_low_z=args.balance_low_z,
            balance_low_z_fraction=args.balance_low_z_fraction,
            orient_single_anchor_up=True,
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
        box_line,
        dna_mode=args.dna_mode,
    )

    print(f"✔ Saved oriented system → {args.out}")


if __name__ == "__main__":
    main()
