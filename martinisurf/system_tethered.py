#!/usr/bin/env python3
"""
MartiniSurf Orientation Assistant
"""

import argparse
from typing import List, Tuple, Sequence
import numpy as np

# ======================================================================
# GRO LOADING
# ======================================================================
def load_gro_coords(gro_file: str) -> Tuple[np.ndarray, List[Tuple[int, str, str, int]]]:
    coords = []
    atoms = []
    with open(gro_file, "r") as fh:
        lines = fh.readlines()[2:-1]

    for line in lines:
        try:
            resid = int(line[0:5])
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            atomid = int(line[15:20])
            x = float(line[20:28]) * 10.0
            y = float(line[28:36]) * 10.0
            z = float(line[36:44]) * 10.0
            atoms.append((resid, resname, atomname, atomid))
            coords.append([x, y, z])
        except ValueError:
            continue

    return np.array(coords,float), atoms

# ======================================================================
# PDB → GRO
# ======================================================================
def convert_pdb_to_gro(pdb_file: str, gro_file: str) -> str:
    atoms = []
    coords = []
    with open(pdb_file, "r") as fh:
        for line in fh:
            if line.startswith(("ATOM","HETATM")):
                atomname = line[12:16].strip()
                resname = line[17:20].strip()
                resid = int(line[22:26])
                x = float(line[30:38])/10
                y = float(line[38:46])/10
                z = float(line[46:54])/10
                atoms.append((resid,resname,atomname,len(atoms)+1))
                coords.append([x,y,z])

    coords = np.array(coords)

    with open(gro_file,"w") as fh:
        fh.write("Converted from PDB by MartiniSurf\n")
        fh.write(f"{len(coords):5d}\n")
        resid_map={}
        new=0

        for (res,resn,a,aid),(x,y,z) in zip(atoms,coords):
            if res not in resid_map:
                new+=1
                resid_map[res]=new
            fh.write(f"{resid_map[res]:5d}{resn:<5}{a:>5}{aid:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")

        fh.write("   10.00000   10.00000   10.00000\n")

    return gro_file

# ======================================================================
# RESID CENTROIDS
# ======================================================================
def summarize_selected_residues(residue_list, atom_records, coords):
    centroids=[]
    msgs=[]
    for resid in residue_list:
        sel=[(r,rn,a,aid,c) for (r,rn,a,aid),c in zip(atom_records,coords) if r==resid]
        if not sel:
            msgs.append(f"⚠ Residue {resid} not found.")
            continue
        loc=np.array([c for *_,c in sel])
        centroids.append(loc.mean(0))
        msgs.append(f"• Residue {resid:4d}")

    print("\n=== Selected Anchor Residues ===")
    print("\n".join(msgs))
    print("================================\n")
    return np.array(centroids)

# ======================================================================
# ORIENTATION ENGINE
# ======================================================================
def auto_orient_from_anchor_residues(system_coords, anchor_centroids, surface_coords, target_z):
    anchors=anchor_centroids.copy()
    centroid=anchors.mean(0)

    A=anchors-centroid
    _,_,vh=np.linalg.svd(A)
    normal=vh[2]/np.linalg.norm(vh[2])

    target_normal=np.array([0,0,-1])
    v=np.cross(normal,target_normal)
    s=np.linalg.norm(v)
    c=float(np.dot(normal,target_normal))

    if s<1e-6:
        R=np.eye(3)
    else:
        vx=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
        R=np.eye(3)+vx+vx@vx*((1-c)/(s*s))

    center=system_coords.mean(0)
    rot=(R@(system_coords-center).T).T+center
    anchors_rot=(R@(anchors-center).T).T+center

    if anchors_rot[:,2].mean()>center[2]:
        Rflip=np.diag([1,1,-1])
        rot=(Rflip@(rot-center).T).T+center
        anchors_rot=(Rflip@(anchors_rot-center).T).T+center

    z_anchor=anchors_rot[:,2].min()
    z_surface=surface_coords[:,2].max()
    rot[:,2]+=(z_surface+target_z)-z_anchor

    xy_shift=surface_coords[:,:2].mean(0)-rot[:,:2].mean(0)
    rot[:,0]+=xy_shift[0]
    rot[:,1]+=xy_shift[1]
    return rot

### READ SURFACE BOX DIMENSIONS
def read_box_from_gro(gro_file):
    with open(gro_file) as fh:
        lines = fh.readlines()
    return lines[-1].strip()

# ======================================================================
# SAVE SYSTEM
# ======================================================================
def save_full_system(output_gro, surf_atoms, surf_coords,
                     enz_atoms, enz_coords, box_line):

    total = len(surf_coords) + len(enz_coords)

    with open(output_gro, "w") as fh:
        fh.write("MartiniSurf oriented system\n")
        fh.write(f"{total:5d}\n")

        next_resid = 1
        next_atom = 1
        enz_map = {}

        for (r,rn,a,aid),(x,y,z) in zip(enz_atoms,enz_coords):
            if r not in enz_map:
                enz_map[r] = next_resid
                next_resid += 1
            fh.write(
                f"{enz_map[r]:5d}{rn:<5}{a:>5}{next_atom:5d}"
                f"{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}\n"
            )
            next_atom += 1

        surf_map = {}
        for (r,rn,a,aid),(x,y,z) in zip(surf_atoms,surf_coords):
            if r not in surf_map:
                surf_map[r] = next_resid
                next_resid += 1
            fh.write(
                f"{surf_map[r]:5d}{rn:<5}{a:>5}{next_atom:5d}"
                f"{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}\n"
            )
            next_atom += 1

        fh.write(box_line + "\n")

    print(f"✔ Saved oriented system → {output_gro}")

# ======================================================================
# MAIN
# ======================================================================
def main(argv=None):

    print("\n==============================================================")
    print("               MartiniSurf – Orientation Assistant")
    print("==============================================================\n")

    parser=argparse.ArgumentParser(description="Automatic CG-system orientation")

    parser.add_argument("--surface",required=True)
    parser.add_argument("--system",required=True)
    parser.add_argument("--out",default="system_Surface.gro")

    # FIXED: REMOVE type=int → now list of strings ALWAYS
    parser.add_argument("--anchor", nargs="+", action="append",
        help="Use: --anchor GROUP RES RES RES")

    parser.add_argument("--dist",type=float,default=10.0)

    args=parser.parse_args(argv)

    # ======================================================
    # PARSE ANCHOR GROUPS
    # ======================================================
    anchor_groups=[]

    if args.anchor:
        for group_def in args.anchor:

            if not isinstance(group_def,(list,tuple)):
                raise ValueError(f"❌ Invalid anchor input: {group_def}")

            if len(group_def)<2:
                raise ValueError(f"❌ Anchor requires GROUP RESID... Got {group_def}")

            group_id = group_def[0]       # ignored for orientation
            residue_ids = group_def[1:]   # real residues

            try:
                residue_ids=[int(x) for x in residue_ids]
            except:
                raise ValueError(f"❌ Invalid residues: {group_def}")

            anchor_groups.append(residue_ids)

    if not anchor_groups:
        raise ValueError("❌ No anchors provided")

    print("\n=== Anchor groups===")
    for i,g in enumerate(anchor_groups,1):
        print(f"Anchor_{i}: {g}")

    all_res = sorted({r for g in anchor_groups for r in g})
    print("Flattened:",all_res,"\n")

    # ======================================================
    # LOAD SYSTEM
    # ======================================================
    if args.system.endswith(".pdb"):
        out = args.system.replace(".pdb",".gro")
        convert_pdb_to_gro(args.system,out)
        args.system=out

    surf_coords,surf_atoms = load_gro_coords(args.surface)
    enz_coords, enz_atoms = load_gro_coords(args.system)

    box_line = read_box_from_gro(args.surface)

    centroids = summarize_selected_residues(all_res, enz_atoms, enz_coords)

    oriented = auto_orient_from_anchor_residues(
        enz_coords, centroids, surf_coords, args.dist
    )

    save_full_system(
    args.out,
    surf_atoms, surf_coords,
    enz_atoms, oriented,
    box_line
    )


# ======================================================================
if __name__=="__main__":
    main()
