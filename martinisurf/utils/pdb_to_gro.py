#!/usr/bin/env python3
"""
Convert PDB → GRO using MDTraj (no GROMACS required).

Usage:
    from pdb_to_gro import pdb_to_gro
    pdb_to_gro("input.pdb", "output.gro")
"""

import mdtraj as md


def pdb_to_gro(pdb_path, gro_path):
    """
    Convert a PDB file into a GRO file.

    Parameters
    ----------
    pdb_path : str
        Input PDB filename
    gro_path : str
        Output GRO filename
    """

    print(f"🔄 Converting PDB → GRO (MDTraj): {pdb_path} → {gro_path}")

    # Load PDB
    traj = md.load(pdb_path)
    top = traj.topology

    # MDTraj writes GRO natively → perfect for Martini workflows
    traj.save_gro(gro_path)

    print("✔ GRO written:", gro_path)
    return gro_path
