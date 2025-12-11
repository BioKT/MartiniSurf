#!/usr/bin/env python3
"""
GoMartini System Builder with Anchoring Restraints

Automatically detects if the system is PROTEIN (Martini3) or DNA (Martini2.1-dna)
and builds the correct topology.

This module prepares a complete GōMartini simulation directory:
- Detects the molecule type based on residue names
- Creates the 0_topology / 1_mdp / 2_system structure
- Copies proper Martini/GōMartini .itp files
- Builds system.top and system_anchor.top
- Generates Active_anchor.itp with position restraints
- Builds index.ndx

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
    Path(path).mkdir(parents=True, exist_ok=True)


def write_list(values: List[int], fh, chunk: int = 15) -> None:
    for i in range(0, len(values), chunk):
        fh.write(" ".join(str(v) for v in values[i: i + chunk]) + "\n")


# ======================================================================
# MAIN WORKFLOW
# ======================================================================
def main(argv: Sequence[str] | None = None) -> None:

    parser = argparse.ArgumentParser(description="Build GoMartini system with anchoring restraints.")
    parser.add_argument("--moltype", required=True)
    parser.add_argument("--outdir", default="Simulation")

    parser.add_argument(
        "--anchor",
        nargs="+",
        action="append",
        metavar=("GROUP", "RESID"),
        help="Define anchor groups as: GROUP RESID [RESID ...]. Repeatable.",
    )

    args = parser.parse_args(argv) if argv else parser.parse_args()

    # ===============================================================
    # Collect anchor groups
    # ===============================================================
    if args.anchor is None:
        raise ValueError("❌ No anchors provided. Use --anchor GROUP RES ... (repeatable)")

    anchor_groups = {}
    for group in args.anchor:
        gid = int(group[0])
        residues = list(map(int, group[1:]))
        anchor_groups[gid] = residues

    anchor_groups = dict(sorted(anchor_groups.items()))

    print("\n=== Anchor Groups Provided ===")
    for gid, res in anchor_groups.items():
        print(f"  Anchor_{gid}: {res}")
    print()

    # ===============================================================
    # Find immobilized_system.gro
    # ===============================================================
    cwd = Path.cwd()
    candidates = [
        cwd / "immobilized_system.gro",
        cwd / "2_system" / "immobilized_system.gro",
        Path(args.outdir).resolve() / "2_system" / "immobilized_system.gro",
    ]

    input_gro = next((c for c in candidates if c.exists()), None)
    if input_gro is None:
        raise FileNotFoundError("❌ immobilized_system.gro not found.")

    output_root = input_gro.parent.parent.resolve()

    print(f"• Using input system: {input_gro}")
    print(f"• Output directory:   {output_root}\n")

    u = mda.Universe(str(input_gro))

    # ===============================================================
    # AUTO-DETECT DNA vs PROTEIN
    # ===============================================================
    dna_resnames = {"DA", "DT", "DC", "DG", "DAN", "THY", "GUA", "CYT"}

    is_dna = any(res in dna_resnames for res in u.residues.resnames)

    print("🔍 Molecule type detected:", "DNA" if is_dna else "PROTEIN/MARTINI3")

    # ===============================================================
    # Build atom lists for each anchor group
    # ===============================================================
    anchor_atoms = {}

    for gid, residues in anchor_groups.items():
        atom_list = [
            int(a.index + 1)
            for a in u.atoms
            if a.resnum in residues
        ]
        anchor_atoms[gid] = atom_list
        print(f"  → Anchor_{gid}: {len(atom_list)} atoms")

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

    # ===============================================================
    # Copy Martini .itp files
    # ===============================================================
    pkg_dir = Path(martinisurf.__file__).parent
    src_itp_dir = pkg_dir / "system_itp"

    dst_itp_dir = topo_dir / "system_itp"
    ensure_dir(dst_itp_dir)

    print("\n• Copying system_itp files...")

    if is_dna:
        required = [
            "martini_v2.1-dna.itp",
            "martini_v2.1P-dna.itp",
            "martini_v2.0_ions.itp",
        ]
    else:
        required = [
            "martini_v3.0.0.itp",
            "martini_v3.0.0_Active.itp",
            "martini_v3.0.0_solvents_v1.itp",
            "martini_v3.0.0_ions_v1.itp",
        ]

    for fname in required:
        src = src_itp_dir / fname
        if src.exists():
            shutil.copy(src, dst_itp_dir / fname)
            print(f"  ✔ Copied {fname}")
        else:
            print(f"  ⚠ Missing expected ITP: {fname}")

    # ===============================================================
    # Locate molecule ITP
    # ===============================================================
    print("\n• Detecting molecule ITP file...")

    # Protein: simple
    mol = args.moltype

    if not is_dna:
        possible_itps = list(dst_itp_dir.glob(f"{mol}.itp"))
    else:
        possible_itps = (
            list(dst_itp_dir.glob("*DNA*.itp"))
            + list(dst_itp_dir.glob("*Nucleic*.itp"))
        )

    if not possible_itps:
        raise FileNotFoundError(f"❌ No ITP found for molecule {mol} in {dst_itp_dir}")

    mol_itp = possible_itps[0]
    print(f"✔ Using molecule ITP: {mol_itp.name}")

    # ===============================================================
    # Build system.top
    # ===============================================================
    print("\n• Generating system.top ...")

    u2 = mda.Universe(str(input_gro))

    n_mol = 1
    n_surf = len({a.resid for a in u2.atoms if a.resname == "SRF"})
    n_water = len({a.resid for a in u2.atoms if a.resname in ["W", "WAT"]})
    n_ions = len({a.resid for a in u2.atoms if a.resname.upper() in ["NA", "CL"]})

    with open(topo_dir / "system.top", "w") as f:

        if is_dna:
            f.write('#include "system_itp/martini_v2.1-dna.itp"\n')
            if (dst_itp_dir / "martini_v2.1P-dna.itp").exists():
                f.write('#include "system_itp/martini_v2.1P-dna.itp"\n')
            f.write('#include "system_itp/martini_v2.0_ions.itp"\n')
        else:
            f.write('#include "system_itp/martini_v3.0.0.itp"\n')
            f.write('#include "system_itp/martini_v3.0.0_Active.itp"\n')
            f.write('#include "system_itp/martini_v3.0.0_solvents_v1.itp"\n')
            f.write('#include "system_itp/martini_v3.0.0_ions_v1.itp"\n')

        f.write(f'#include "system_itp/{mol_itp.name}"\n')
        f.write('#include "system_itp/surface.itp"\n\n')

        f.write("[ system ]\nGoMartini Surface Simulation\n\n")
        f.write("[ molecules ]\n")
        f.write(f"{args.moltype}   {n_mol}\n")
        f.write(f"SRF      {n_surf}\n")
        f.write(f"W        {n_water}\n")
        f.write(f"Na       {n_ions}\n")

    # ===============================================================
    # Build Active_anchor.itp
    # ===============================================================
    print("\n• Generating Active_anchor.itp ...")

    active_anchor = dst_itp_dir / f"{mol_itp.stem}_anchor.itp"

    new_lines = []
    inside_posres = False

    with open(mol_itp) as fin:
        for line in fin:
            if "[ position_restraints" in line:
                inside_posres = True
                continue
            if inside_posres:
                if line.strip().startswith("["):
                    inside_posres = False
                    new_lines.append(line)
                continue
            new_lines.append(line)

    new_lines.append("\n[ position_restraints ]\n#ifdef POSRES\n")

    for gid, atoms in anchor_atoms.items():
        for atom in atoms:
            new_lines.append(f"{atom} 1 1000 1000 0\n")

    new_lines.append("#endif\n")

    with open(active_anchor, "w") as fout:
        fout.write("".join(new_lines))

    # ===============================================================
    # Build system_anchor.top
    # ===============================================================
    print("• Generating system_anchor.top ...")

    with open(topo_dir / "system_anchor.top", "w") as f:
        if is_dna:
            f.write('#include "system_itp/martini_v2.1-dna.itp"\n')
            if (dst_itp_dir / "martini_v2.1P-dna.itp").exists():
                f.write('#include "system_itp/martini_v2.1P-dna.itp"\n')
            f.write('#include "system_itp/martini_v2.0_ions.itp"\n')
        else:
            f.write('#include "system_itp/martini_v3.0.0.itp"\n')
            f.write('#include "system_itp/martini_v3.0.0_Active.itp"\n')
            f.write('#include "system_itp/martini_v3.0.0_solvents_v1.itp"\n')
            f.write('#include "system_itp/martini_v3.0.0_ions_v1.itp"\n')

        f.write(f'#include "system_itp/{active_anchor.name}"\n')
        f.write('#include "system_itp/surface.itp"\n\n')

        f.write("[ system ]\nGoMartini Surface Simulation (restr)\n\n")
        f.write("[ molecules ]\n")
        f.write(f"{args.moltype}   {n_mol}\n")
        f.write(f"SRF      {n_surf}\n")
        f.write(f"W        {n_water}\n")
        f.write(f"Na       {n_ions}\n")

    # ===============================================================
    # Build index.ndx
    # ===============================================================
    print("• Generating index.ndx ...")

    with open(topo_dir / "index.ndx", "w") as ndx:
        for gid, atoms in anchor_atoms.items():
            ndx.write(f"\n[ Anchor_{gid} ]\n")
            write_list(atoms, ndx)

    # ===============================================================
    # Copy MDP templates
    # ===============================================================
    print("• Copying MDP templates ...")

    mdp_pkg = pkg_dir / "mdp_templates"
    for fname in [
        "gromacs_workflow",
        "minimization.mdp",
        "nvt.mdp",
        "npt.mdp",
        "deposition.mdp",
        "production.mdp",
    ]:
        src = mdp_pkg / fname
        if src.exists():
            shutil.copy(src, mdp_dir / fname)
        else:
            print(f"⚠ Missing MDP template: {fname}")

    # ===============================================================
    print("\n========================================")
    print("✔ GoMartini SYSTEM SUCCESSFULLY BUILT")
    print(f"  Output directory: {output_root}")
    print("========================================\n")


if __name__ == "__main__":
    main()
