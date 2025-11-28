#!/usr/bin/env python3
"""
GoMartini System Builder with Anchoring Restraints

This module prepares a complete GōMartini simulation directory:
- Detects input enzyme-on-surface system
- Creates standard 0_topology / 1_mdp / 2_system layout
- Copies the correct Martini GoITP and surface topologies
- Generates system.top and system_res.top
- Creates anchoring position restraints via Active_res.itp
- Builds index.ndx groups for A/B anchor residues

This is the final step of the MartiniSurf pipeline.
"""

import os
import shutil
import argparse
import MDAnalysis as mda
import surfmartini
from pathlib import Path
from typing import Iterable, List


# ======================================================================
# UTILITIES
# ======================================================================
def ensure_dir(path: str) -> None:
    """
    Create directory if it does not already exist.

    Parameters
    ----------
    path : str
        Directory path to ensure exists.
    """
    os.makedirs(path, exist_ok=True)


def write_list(values: List[int], fh, chunk: int = 15) -> None:
    """
    Write integer values grouped in fixed-size chunks.

    Parameters
    ----------
    values : list of int
        Atom IDs or indices to write.
    fh : file-like
        Output file handle.
    chunk : int
        Number of values per line.
    """
    for i in range(0, len(values), chunk):
        fh.write(" ".join(str(x) for x in values[i : i + chunk]) + "\n")


# ======================================================================
# MAIN WORKFLOW
# ======================================================================
def main(argv: Iterable[str] | None = None) -> None:
    """
    Build GoMartini simulation system with anchoring restraints.

    Parameters
    ----------
    argv : iterable of str, optional
        Command-line arguments. Default: sys.argv.
    """
    parser = argparse.ArgumentParser(
        description="Build GoMartini system with anchoring restraints."
    )

    parser.add_argument("--resA", nargs="+", type=int)
    parser.add_argument("--resB", nargs="+", type=int)
    parser.add_argument("--anchor", nargs="+", type=int)
    parser.add_argument("--outdir", default="Simulation")

    args = parser.parse_args(argv) if argv is not None else parser.parse_args()

    # ==================================================================
    # Determine anchor residues
    # ==================================================================
    if args.resA or args.resB:
        print("\n• Using explicit resA/resB selections")
        resA = args.resA or []
        resB = args.resB or []
    elif args.anchor:
        print("\n• Using unified --anchor list for both A and B")
        resA = resB = list(args.anchor)
    else:
        raise ValueError("You must provide --resA/--resB or --anchor.")

    print(f"  → Residues A = {resA}")
    print(f"  → Residues B = {resB}\n")

    # ==================================================================
    # Detect package paths
    # ==================================================================
    package_dir = Path(surfmartini.__file__).parent
    activeitp_pkg = package_dir / "ActiveITP"
    mdp_pkg = package_dir / "mdp_templates"

    print(f"• ActiveITP path:    {activeitp_pkg}")
    print(f"• MDP templates path:{mdp_pkg}")

    # ==================================================================
    # Locate the Enzyme_Surface.gro file
    # ==================================================================
    cwd = Path(os.getcwd())

    candidates = [
        cwd / "Enzyme_Surface.gro",
        cwd / "2_system" / "Enzyme_Surface.gro",
        Path(args.outdir).resolve() / "2_system" / "Enzyme_Surface.gro",
    ]

    for cand in candidates:
        if cand.exists():
            input_gro = cand
            break
    else:
        raise FileNotFoundError(
            "Could not find Enzyme_Surface.gro in any expected location."
        )

    output_root = input_gro.parent.parent.resolve()

    print(f"\n• Using input system: {input_gro}")
    print(f"• Output directory:   {output_root}\n")

    u = mda.Universe(str(input_gro))

    # ==================================================================
    # Map residues to atom IDs (excluding backbone CA)
    # ==================================================================
    residues_a_atoms = [
        a.index + 1 for a in u.atoms if (a.resnum in resA and a.name != "CA")
    ]
    residues_b_atoms = [
        a.index + 1 for a in u.atoms if (a.resnum in resB and a.name != "CA")
    ]

    print(f"  → Group A atoms: {len(residues_a_atoms)}")
    print(f"  → Group B atoms: {len(residues_b_atoms)}\n")

    # ==================================================================
    # Create standard folder structure
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
    dst_active = topo_dir / "ActiveITP"
    ensure_dir(dst_active)

    print("• Copying ActiveITP files...")
    for fname in os.listdir(activeitp_pkg):
        src = activeitp_pkg / fname
        dst = dst_active / fname

        if src.is_file() and not dst.exists():
            shutil.copy(src, dst)
            print(f"  ✔ Copied {fname}")

    # Ensure Active.itp exists (tests + real workflows)
    active_alias = os.path.join(dst_active, "Active.itp")
    source_active = os.path.join(dst_active, "martini_v3.0.0_Active.itp")

    if os.path.exists(active_alias):
        print("  ✔ Found existing Active.itp")
    else:
        shutil.copy(source_active, active_alias)
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

    topfile = topo_dir / "system.top"
    with open(topfile, "w") as ftop:
        ftop.write("#define GO_VIRT\n\n")
        ftop.write('#include "ActiveITP/martini_v3.0.0_Active.itp"\n')
        ftop.write('#include "ActiveITP/martini_v3.0.0_solvents_v1.itp"\n')
        ftop.write('#include "ActiveITP/martini_v3.0.0_ions_v1.itp"\n')
        ftop.write('#include "ActiveITP/Active.itp"\n')
        ftop.write('#include "ActiveITP/surface.itp"\n\n')

        ftop.write("[ system ]\nGoMartini Surface Simulation\n\n")
        ftop.write("[ molecules ]\n")
        ftop.write(f"Active   {n_active}\n")
        ftop.write(f"SRF      {n_surf}\n")
        ftop.write(f"W        {n_water}\n")
        ftop.write(f"Na       {n_ions}\n")

    # ==================================================================
    # Build Active_res.itp
    # ==================================================================
    print("• Generating Active_res.itp ...")

    active_src = dst_active / "Active.itp"
    active_res = dst_active / "Active_res.itp"

    inside_posres = False
    new_content = []

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

    for atom in residues_a_atoms + residues_b_atoms:
        new_content.append(f"{atom} 1 1000 1000 0\n")

    new_content.append("#endif\n")

    with open(active_res, "w") as fout:
        fout.write("".join(new_content))

    # ==================================================================
    # Build system_res.top
    # ==================================================================
    print("• Generating system_res.top ...")

    topfile_res = topo_dir / "system_res.top"
    with open(topfile_res, "w") as ftop:
        ftop.write("#define GO_VIRT\n\n")
        ftop.write('#include "ActiveITP/martini_v3.0.0_Active.itp"\n')
        ftop.write('#include "ActiveITP/martini_v3.0.0_solvents_v1.itp"\n')
        ftop.write('#include "ActiveITP/martini_v3.0.0_ions_v1.itp"\n')
        ftop.write('#include "ActiveITP/Active_res.itp"\n')
        ftop.write('#include "ActiveITP/surface.itp"\n\n')

        ftop.write("[ system ]\nGoMartini Surface Simulation (restr)\n\n")
        ftop.write("[ molecules ]\n")
        ftop.write(f"Active   {n_active}\n")
        ftop.write(f"SRF      {n_surf}\n")
        ftop.write(f"W        {n_water}\n")
        ftop.write(f"Na       {n_ions}\n")

    # ==================================================================
    # Build index.ndx
    # ==================================================================
    print("• Generating index.ndx ...")

    indexfile = topo_dir / "index.ndx"
    with open(indexfile, "w") as ndx:
        ndx.write("[ Residues_A ]\n")
        write_list(residues_a_atoms, ndx)
        ndx.write("\n[ Residues_B ]\n")
        write_list(residues_b_atoms, ndx)

    # ==================================================================
    # Copy MDP files
    # ==================================================================
    print("• Copying MDP templates ...")

    for fname in ["nvt.mdp", "npt.mdp", "deposition.mdp", "production.mdp"]:
        src = mdp_pkg / fname
        dst = mdp_dir / fname
        if src.exists():
            shutil.copy(src, dst)
        else:
            print(f"  ⚠ Template not found: {src}")

    print("\n========================================")
    print("✔ GoMartini SYSTEM SUCCESSFULLY BUILT")
    print(f"  Output directory: {output_root}")
    print("========================================\n")


# ======================================================================
if __name__ == "__main__":
    main()
