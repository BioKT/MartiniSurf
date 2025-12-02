#!/usr/bin/env python3
"""
GoMartini System Builder with Anchoring Restraints

This module prepares a complete GōMartini simulation directory:
- Detects the enzyme-on-surface system
- Creates the 0_topology / 1_mdp / 2_system structure
- Copies Martini/GōMartini .itp files
- Builds system.top and system_anchor.top
- Generates Active_anchor.itp with position restraints
- Generates index.ndx (Residues_A / Residues_B)

Final step of the MartiniSurf pipeline.
"""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path
from typing import List, Sequence

import MDAnalysis as mda

import martinisurf


# ======================================================================
# UTILITIES
# ======================================================================
def ensure_dir(path: str | Path) -> None:
    """
    Create a directory if it does not already exist.

    Parameters
    ----------
    path : str or Path
        Directory path to create.
    """
    Path(path).mkdir(parents=True, exist_ok=True)


def write_list(values: List[int], fh, chunk: int = 15) -> None:
    """
    Write integers grouped in fixed-size chunks.

    Parameters
    ----------
    values : list[int]
        Atom IDs to write.
    fh : file-like
        Output file handle.
    chunk : int
        Number of entries per line.
    """
    for i in range(0, len(values), chunk):
        fh.write(" ".join(str(v) for v in values[i: i + chunk]) + "\n")


# ======================================================================
# MAIN WORKFLOW
# ======================================================================
def main(argv: Sequence[str] | None = None) -> None:
    """
    Build a GōMartini simulation system with anchoring restraints.

    Parameters
    ----------
    argv : sequence[str] or None
        Command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Build GoMartini system with anchoring restraints."
    )

    parser.add_argument("--moltype", required=True)
    #parser.add_argument("--resA", nargs="+", type=int)
    #parser.add_argument("--resB", nargs="+", type=int)
    #parser.add_argument("--anchor", nargs="+", type=int)
    parser.add_argument("--outdir", default="Simulation")
    # Multi-anchor input system 
    parser.add_argument(
    "--anchor",
    nargs="+",
    action="append",
    metavar=("GROUP", "RESID"),
    help=(
        "Define an anchor group as: GROUP RESID [RESID ...]. "
        "Repeat this flag to define multiple groups."
    ),
)

    args = parser.parse_args(argv) if argv is not None else parser.parse_args()

    # ===============================================================
    # Collect unlimited anchor groups
    # ===============================================================
 
    if args.anchor is None:
        raise ValueError("❌ No anchors provided. Use --anchor GROUP RES ... (repeatable)")

    # Convert to dict: { group_id: [res1, res2...] }
    anchor_groups = {}

    for group in args.anchor:
        group_id = int(group[0])
        residues = list(map(int, group[1:]))

        if group_id in anchor_groups:
            raise ValueError(f"❌ Duplicate group ID {group_id} in --anchor")

        anchor_groups[group_id] = residues

    # Sort numerically by group ID
    anchor_groups = dict(sorted(anchor_groups.items()))

    print("\n=== Anchor Groups Provided ===")
    for gid, reslist in anchor_groups.items():
        print(f"  Anchor_{gid}: {reslist}")
    print()


    # ==================================================================
    # Detect package paths
    # ==================================================================
    package_dir = Path(martinisurf.__file__).parent
    activeitp_pkg = package_dir / "system_itp"
    mdp_pkg = package_dir / "mdp_templates"

    print(f"• system_itp path:     {activeitp_pkg}")
    print(f"• MDP templates path: {mdp_pkg}")

    # ==================================================================
    # Locate immobilized_system.gro
    # ==================================================================
    cwd = Path.cwd()

    candidates = [
        cwd / "immobilized_system.gro",
        cwd / "2_system" / "immobilized_system.gro",
        Path(args.outdir).resolve() / "2_system" / "immobilized_system.gro",
    ]

    input_gro: Path | None = None
    for cand in candidates:
        if cand.exists():
            input_gro = cand
            break

    if input_gro is None:
        raise FileNotFoundError(
            "Could not find immobilized_system.gro in any expected location."
        )

    output_root = input_gro.parent.parent.resolve()

    print(f"\n• Using input system: {input_gro}")
    print(f"• Output directory:   {output_root}\n")

    u = mda.Universe(str(input_gro))

    # ===============================================================
    # Build atom lists for each anchor group
    # ===============================================================
    anchor_atoms = {}  # {group_id: [atom_ids]}

    for group_id, residues in anchor_groups.items():

        atom_list = [
            int(a.index + 1)
            for a in u.atoms
            if (a.resnum in residues and a.name != "CA")
        ]

        anchor_atoms[group_id] = atom_list
        print(f"  → Anchor_{group_id}: {len(atom_list)} atoms")

    # ==================================================================
    # Create folder structure
    # ==================================================================
    topo_dir = output_root / "0_topology"
    mdp_dir = output_root / "1_mdp"
    sys_dir = output_root / "2_system"

    ensure_dir(topo_dir)
    ensure_dir(mdp_dir)
    ensure_dir(sys_dir)

    shutil.copy(input_gro, sys_dir / "system.gro")

    # ==================================================================
    # Copy ActiveITP files
    # ==================================================================
    dst_active = topo_dir / "system_itp"
    ensure_dir(dst_active)

    print("• Copying system_itp files...")
    for fname in os.listdir(activeitp_pkg):
        src = activeitp_pkg / fname
        dst = dst_active / fname

        if src.is_file() and not dst.exists():
            shutil.copy(src, dst)
            print(f"  ✔ Copied {fname}")

    #Ensure Active.itp exists
    #source_active = dst_active / "martini_v3.0.0_Active.itp"

    #if active_alias.exists():
    #print("  ✔ Found existing Active.itp")
    #else:
    #shutil.copy(source_active, active_alias)
    #print("  ✔ Created Active.itp from martini_v3.0.0_Active.itp")

    # ==================================================================
    # Build system.top
    # ==================================================================
    print("\n• Generating system.top ...")

    u2 = mda.Universe(str(input_gro))

    n_active = 1
    n_surf = len({a.resid for a in u2.atoms if a.resname == "SRF"})
    n_water = len({a.resid for a in u2.atoms if a.resname in ["W", "WAT"]})
    n_ions = len({a.resid for a in u2.atoms if a.resname.upper() in ["NA", "CL"]})
    mol = args.moltype

    topfile = topo_dir / "system.top"
    with open(topfile, "w") as ftop:
        ftop.write("#define GO_VIRT\n\n")
        ftop.write('#include "system_itp/martini_v3.0.0_Active.itp"\n')
        ftop.write('#include "system_itp/martini_v3.0.0_solvents_v1.itp"\n')
        ftop.write('#include "system_itp/martini_v3.0.0_ions_v1.itp"\n')
        ftop.write(f'#include "system_itp/{mol}.itp"\n')
        ftop.write('#include "system_itp/surface.itp"\n\n')
        ftop.write("[ system ]\nGoMartini Surface Simulation\n\n")
        ftop.write("[ molecules ]\n")
        ftop.write(f"{mol}   {n_active}\n")
        ftop.write(f"SRF      {n_surf}\n")
        ftop.write(f"W        {n_water}\n")
        ftop.write(f"Na       {n_ions}\n")

    # ==================================================================
    # Build Active_anchor.itp
    # ==================================================================
    print("• Generating Active_anchor.itp ...")

    active_src = dst_active / f"{mol}.itp"
    active_res = dst_active / f"{mol}_anchor.itp"

    new_content: List[str] = []
    inside_posres = False

    with open(active_src) as fin:
        for line in fin:
            if "[ position_restraints" in line:
                inside_posres = True
                continue

            if inside_posres:
                if line.strip().startswith("["):
                    inside_posres = False
                    new_content.append(line)
                continue

            new_content.append(line)

    new_content.append("\n[ position_restraints ]\n")
    new_content.append("#ifdef POSRES\n")

    new_content.append("\n[ position_restraints ]\n")
    new_content.append("#ifdef POSRES\n")

    for group_id, atom_list in anchor_atoms.items():
        for atomid in atom_list:
            new_content.append(f"{atomid} 1 1000 1000 0\n")

    new_content.append("#endif\n")

    with open(active_res, "w") as fout:
        fout.write("".join(new_content))

    # ==================================================================
    # Build system_anchor.top
    # ==================================================================
    print("• Generating system_anchor.top ...")

    topfile_res = topo_dir / "system_anchor.top"
    with open(topfile_res, "w") as ftop:
        ftop.write("#define GO_VIRT\n\n")
        ftop.write('#include "system_itp/martini_v3.0.0_Active.itp"\n')
        ftop.write('#include "system_itp/martini_v3.0.0_solvents_v1.itp"\n')
        ftop.write('#include "system_itp/martini_v3.0.0_ions_v1.itp"\n')
        ftop.write(f'#include "system_itp/{mol}_anchor.itp"\n')
        ftop.write('#include "system_itp/surface.itp"\n\n')

        ftop.write("[ system ]\nGoMartini Surface Simulation (restr)\n\n")
        ftop.write("[ molecules ]\n")
        ftop.write(f"{mol}   {n_active}\n")
        ftop.write(f"SRF      {n_surf}\n")
        ftop.write(f"W        {n_water}\n")
        ftop.write(f"Na       {n_ions}\n")

    # ==================================================================
    # Build index.ndx
    # ==================================================================
    print("• Generating index.ndx ...")

    indexfile = topo_dir / "index.ndx"
    with open(indexfile, "w") as ndx:

        # First: all anchors merged
        ndx.write("[ Anchor_All ]\n")
        write_list([atom for lst in anchor_atoms.values() for atom in lst], ndx)

        # Then: individual groups
        for group_id, atom_list in anchor_atoms.items():
            ndx.write(f"\n[ Anchor_{group_id} ]\n")
            write_list(atom_list, ndx)

    # ==================================================================
    # Copy MDP files
    # ==================================================================
    print("• Copying MDP templates ...")

    for fname in [
    "gromacs_workflow",
    "minimization.mdp",
    "nvt.mdp",
    "npt.mdp",
    "deposition.mdp",
    "production.mdp",
    ]:

        src = mdp_pkg / fname
        dst = mdp_dir / fname
        if src.exists():
            shutil.copy(src, dst)
        else:
            print(f"  ⚠ Missing MDP template: {src}")

    # ==================================================================
    print("\n========================================")
    print("✔ GoMartini SYSTEM SUCCESSFULLY BUILT")
    print(f"  Output directory: {output_root}")
    print("========================================\n")


# ======================================================================
if __name__ == "__main__":
    main()
