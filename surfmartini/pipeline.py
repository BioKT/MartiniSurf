#!/usr/bin/env python3
"""
MartiniSurf Full Pipeline

This script automates the complete workflow:

1. Martinize2  → Coarse-grained enzyme model
2. Surface     → Build hexagonal Martini surface
3. Orientation → Position enzyme above surface using anchor residues
4. GoMartini   → Build full simulation system (topology, restraints, mdp)

This is the highest-level entry point of MartiniSurf.
"""

import argparse
import shutil
import subprocess
from pathlib import Path


# ======================================================================
# Utility: run shell command
# ======================================================================
def run(cmd: list[str]) -> None:
    """
    Run a shell command and ensure it finishes successfully.

    Parameters
    ----------
    cmd : list of str
        Command + arguments.

    Raises
    ------
    RuntimeError
        If the command returns a non-zero exit code.
    """
    print(f"\n▶ Running command:\n  {' '.join(cmd)}\n")
    res = subprocess.run(cmd)

    if res.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")

    print("✔ Done\n")


# ======================================================================
# MAIN PIPELINE
# ======================================================================
def main(argv=None) -> None:
    """
    Execute the complete MartiniSurf workflow.

    Parameters
    ----------
    argv : iterable or None
        Optional command-line arguments (for testing).
    """
    parser = argparse.ArgumentParser(
        description=(
            "Full MartiniSurf pipeline: "
            "martinize2 → surface → orient → GoMartini system"
        )
    )

    parser.add_argument("--pdb", required=True, help="Input all-atom PDB")

    # Martinize2
    parser.add_argument("--go", nargs="?", default="auto")
    parser.add_argument("--eps", type=float, default=9.414)
    parser.add_argument("--dssp", default="mkdssp")
    parser.add_argument("--ff", default="martini3001")
    parser.add_argument("--moltype", default="Active")
    parser.add_argument("--merge", default=None)

    # Surface
    parser.add_argument("--surface-bead", default="P4")
    parser.add_argument("--dx", type=float, default=0.47)
    parser.add_argument("--lx", type=float, required=True)
    parser.add_argument("--ly", type=float, required=True)

    # Orientation
    parser.add_argument("--resA", nargs="+", type=int)
    parser.add_argument("--resB", nargs="+", type=int)
    parser.add_argument("--anchor", nargs="+", type=int)
    parser.add_argument("--dist", type=float, default=10.0)

    parser.add_argument("--outdir", default="Simulation_Files")

    args = parser.parse_args(argv)

    # ==================================================================
    # Prepare simulation directory structure
    # ==================================================================
    simdir = Path(args.outdir).resolve()

    print(f"\n📁 Creating simulation directory: {simdir}\n")

    if simdir.exists():
        shutil.rmtree(simdir)

    (simdir / "0_topology" / "ActiveITP").mkdir(parents=True)
    (simdir / "1_mdp").mkdir()
    (simdir / "2_system").mkdir()

    active_itp_dir = simdir / "0_topology" / "ActiveITP"
    system_dir = simdir / "2_system"

    # ==================================================================
    # 1. Martinize2 (always run in CWD)
    # ==================================================================
    tmpdir = simdir / "_martinize_tmp"
    tmpdir.mkdir()

    enzyme_cg_out = (tmpdir / "enzyme_cg.pdb").resolve()
    topfile_out = (tmpdir / "enzyme_cg.top").resolve()

    pdb_abs = Path(args.pdb).resolve()

    martinize_cmd = [
        "martinize2",
        "-f",
        str(pdb_abs),
        "-x",
        str(enzyme_cg_out),
        "-o",
        str(topfile_out),
        "-dssp",
        args.dssp,
        "-ff",
        args.ff,
        "-name",
        args.moltype,
    ]

    if args.go == "auto":
        martinize_cmd += ["-go"]
    elif args.go != "none":
        martinize_cmd += ["-go", args.go]

    martinize_cmd += ["-go-eps", str(args.eps)]

    if args.merge:
        martinize_cmd += ["-merge", args.merge]

    run(martinize_cmd)

    # Move all .itp and .ssd files generated in CWD
    for file in Path(".").glob("*"):
        if file.suffix in [".itp", ".ssd"]:
            shutil.move(str(file), active_itp_dir / file.name)

    shutil.copy(enzyme_cg_out, system_dir / "enzyme_cg.pdb")

    # ==================================================================
    # 2. Surface builder
    # ==================================================================
    print("\n🌐 Building surface...\n")
    import surfmartini.surface_builder as sb

    sb.main(
        [
            "--lx",
            str(args.lx),
            "--ly",
            str(args.ly),
            "--dx",
            str(args.dx),
            "--bead",
            args.surface_bead,
            "--output",
            str(system_dir / "surface"),
        ]
    )

    surface_itp_src = system_dir / "surface.itp"
    surface_itp_dst = active_itp_dir / "surface.itp"

    if surface_itp_src.exists():
        shutil.move(surface_itp_src, surface_itp_dst)
        print(f"✔ Moved surface.itp → {surface_itp_dst}")
    else:
        print("⚠ surface.itp not found in 2_system")

    # ==================================================================
    # 3. Orientation
    # ==================================================================
    print("\n🧬 Orienting enzyme on surface...\n")

    import surfmartini.enzyme_tethered as orient_mod
    from surfmartini.enzyme_tethered import convert_pdb_to_gro

    enzyme_cg_pdb = system_dir / "enzyme_cg.pdb"
    enzyme_cg_gro = system_dir / "enzyme_cg.gro"

    print("• Converting enzyme_cg.pdb → enzyme_cg.gro")
    convert_pdb_to_gro(str(enzyme_cg_pdb), str(enzyme_cg_gro))

    orient_args = [
        "--surface",
        str(system_dir / "surface.gro"),
        "--enzyme",
        str(enzyme_cg_gro),
        "--out",
        str(system_dir / "Enzyme_Surface.gro"),
        "--dist",
        str(args.dist),
    ]

    if args.resA:
        orient_args += ["--resA"] + list(map(str, args.resA))
    if args.resB:
        orient_args += ["--resB"] + list(map(str, args.resB))
    if args.anchor:
        orient_args += ["--anchor"] + list(map(str, args.anchor))

    orient_mod.main(orient_args)

    # ==================================================================
    # 4. GoMartini system builder
    # ==================================================================
    print("\n📦 Building GoMartini system...\n")
    import surfmartini.gomartini_system as gsm

    gsm_args = ["--outdir", str(simdir)]

    if args.resA:
        gsm_args += ["--resA"] + list(map(str, args.resA))
    if args.resB:
        gsm_args += ["--resB"] + list(map(str, args.resB))
    if args.anchor:
        gsm_args += ["--anchor"] + list(map(str, args.anchor))

    gsm.main(gsm_args)

    # ==================================================================
    # 5. Cleanup
    # ==================================================================
    shutil.rmtree(tmpdir)

    print("\n============================================")
    print("Finished Successfully!")
    print(f"  Final directory: {simdir}/")
    print("============================================\n")


# ======================================================================
if __name__ == "__main__":
    main()
