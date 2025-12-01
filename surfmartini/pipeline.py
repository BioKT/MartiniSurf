#!/usr/bin/env python3
"""
MartiniSurf Full Pipeline

This script automates the complete workflow:

1. Martinize2  → Coarse-grained enzyme model
2. Surface     → Build hexagonal Martini surface
3. Orientation → Position enzyme above surface using anchor residues
4. GoMartini   → Build full simulation system (topology, restraints, mdp)
"""

import argparse
import shutil
import subprocess
from pathlib import Path

from surfmartini.utils.pdb_rscb import load_clean_pdb


# ======================================================================
# 1) Parser exposed to main CLI (martinisurf -h)
# ======================================================================
def build_parser():
    parser = argparse.ArgumentParser(
        description="Full MartiniSurf pipeline: martinize2 → surface → orient → GoMartini system",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --------------------------
    # INPUTPDB
    # --------------------------
    parser.add_argument(
        "--pdb", required=True,
        help="Input structure: path/to/file.pdb or RCSB ID (e.g. 1A2B)"
    )

    parser.add_argument("--moltype", required=True,
        help="Name assigned to molecule in topology files")

    parser.add_argument("--merge", default=None,
        help="Merge chains: e.g. 'A,B' or 'all'")

    # --------------------------
    # FORCE FIELD
    # --------------------------
    parser.add_argument("--ff", default="martini3001",
        help="Martini force field")

    # --------------------------
    # POSITION RESTRAINTS
    # --------------------------
    parser.add_argument("--p", choices=["none", "all", "backbone"], default="backbone",
        help="Position restraints")

    parser.add_argument("--pf", type=float, default=1000,
        help="Position restraint force constant (kJ/mol/nm^2)")

    # --------------------------
    # DSSP
    # --------------------------
    parser.add_argument("--dssp", nargs="?", const="mkdssp", default="mkdssp",
        help="DSSP executable or mdtraj fallback")

    # --------------------------
    # ELASTIC NETWORK
    # --------------------------
    parser.add_argument("--elastic", action="store_true",
        help="Enable elastic network")

    parser.add_argument("--ef", type=float, default=700,
        help="Elastic force constant (kJ/mol/nm^2)")

    # --------------------------
    # GoMARTINI
    # --------------------------
    parser.add_argument("--go", nargs="?", const="auto",
        help="Enable GoMartini; optional contact map file")

    parser.add_argument("--go-eps", type=float, default=9.414,
        help="GoMartini interaction energy")

    parser.add_argument("--go-low", type=float, default=0.3,
        help="Min Go contact distance (nm)")

    parser.add_argument("--go-up", type=float, default=1.1,
        help="Max Go contact distance (nm)")

    parser.add_argument("--go-write-file", nargs="?", const=True, default=False,
        help="Write Go contact map")

    # --------------------------
    # PROTEIN OPTIONS
    # --------------------------
    parser.add_argument("--cys", default="auto",
        help="Cysteine bridges (auto/none/custom)")

    parser.add_argument("--mutate", nargs="+", default=[],
        help="Residue mutations: e.g. A-THR6:ALA")

    # --------------------------
    # SURFACE
    # --------------------------
    parser.add_argument("--surface-bead", default="P4",
        help="Martini bead used for the surface")

    parser.add_argument("--dx", type=float, default=0.47,
        help="Surface lattice spacing (nm)")

    parser.add_argument("--lx", type=float, required=True,
        help="Surface X length (nm)")

    parser.add_argument("--ly", type=float, required=True,
        help="Surface Y length (nm)")

    # --------------------------
    # ORIENTATION
    # --------------------------
    parser.add_argument("--resA", nargs="+", type=int,
        help="Anchor residues (set A)")

    parser.add_argument("--resB", nargs="+", type=int,
        help="Anchor residues (set B)")

    parser.add_argument("--anchor", nargs="+", type=int,
        help="Explicit anchors (override resA/resB)")

    parser.add_argument("--dist", type=float, default=10.0,
        help="Vertical enzyme–surface distance (Å)")

    # --------------------------
    # OUTPUT
    # --------------------------
    parser.add_argument("--outdir", default="Simulation_Files",
        help="Output directory")

    # --------------------------
    # Raw passthrough
    # --------------------------
    parser.add_argument("--m2-args", nargs=argparse.REMAINDER,
        help="Extra raw arguments passed directly to Martinize2")

    return parser


# ======================================================================
# Utility: run command
# ======================================================================
def run(cmd):
    print("\n▶ Running command:\n ", " ".join(cmd), "\n")
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Command failed.")
    print("✔ Done\n")


# ======================================================================
# MAIN pipeline
# ======================================================================
def main(argv=None):

    # Use the unified parser
    parser = build_parser()
    args = parser.parse_args(argv)

    # Create directories
    simdir = Path(args.outdir).resolve()
    if simdir.exists():
        shutil.rmtree(simdir)

    (simdir / "0_topology" / "system_itp").mkdir(parents=True)
    (simdir / "1_mdp").mkdir()
    (simdir / "2_system").mkdir()

    active_itp_dir = simdir / "0_topology" / "system_itp"
    system_dir = simdir / "2_system"

    # =========================================================================
    # 1) CLEAN OR DOWNLOAD PDB
    # =========================================================================
    pdb_abs = load_clean_pdb(args.pdb, workdir=simdir)

    tmpdir = simdir / "_martinize_tmp"
    tmpdir.mkdir()

    enzyme_cg_out = tmpdir / "enzyme_cg.pdb"
    topfile_out   = tmpdir / "enzyme_cg.top"

    # =========================================================================
    # 2) MARTINIZE2 COMMAND
    # =========================================================================
    martinize_cmd = [
        "martinize2",
        "-f", str(pdb_abs),
        "-x", str(enzyme_cg_out),
        "-o", str(topfile_out),
        "-dssp", args.dssp,
        "-ff", args.ff,
        "-name", args.moltype,
    ]

    # Position restraints
    if args.p != "none":
        martinize_cmd += ["-p", args.p]

    martinize_cmd += ["-pf", str(args.pf)]

    # Elastic
    if args.elastic:
        martinize_cmd.append("-elastic")

    if args.elastic:
        martinize_cmd += ["-ef", str(args.ef)]

    # GoMartini
    if args.go is not None:
        martinize_cmd += ["-go"] if args.go == "auto" else ["-go", args.go]

    martinize_cmd += [
        "-go-eps", str(args.go_eps),
        "-go-low", str(args.go_low),
        "-go-up",  str(args.go_up)
    ]

    if args.go_write_file:
        martinize_cmd.append("-go-write-file")

    # MERGE CHAINS
    if args.merge:
        martinize_cmd += ["-merge", args.merge]

    # Cysteines
    martinize_cmd += ["-cys", args.cys]

    # Mutations
    for mut in args.mutate:
        martinize_cmd += ["-mutate", mut]

    # Extra Martinize2 arguments
    if args.m2_args:
        martinize_cmd += args.m2_args

    # Debug print
    #print("\n🧩 Martinize2 flags used:")
    #for a in martinize_cmd:
    #    print("   ", a)

    # Run Martinize2
    run(martinize_cmd)

    # Move generated itps
    for f in Path(".").glob("*.itp"):
        shutil.move(str(f), active_itp_dir / f.name)
    
    # Move generated ssd
    for f in Path(".").glob("*.ssd"):
        shutil.move(str(f), active_itp_dir / f.name)

    # ============================================
    # FIX: Rename mol_0.itp → mol.itp
    # ============================================
    mol = args.moltype
    src = active_itp_dir / f"{mol}_0.itp"
    dst = active_itp_dir / f"{mol}.itp"

    if src.exists():
        src.rename(dst)
        print(f"✔ Renamed {src.name} → {dst.name}")

    shutil.copy(enzyme_cg_out, system_dir / "enzyme_cg.pdb")

    # =========================================================================
    # 3) Build surface
    # =========================================================================
    import surfmartini.surface_builder as sb

    sb.main([
        "--lx", str(args.lx),
        "--ly", str(args.ly),
        "--dx", str(args.dx),
        "--bead", args.surface_bead,
        "--output", str(system_dir / "surface"),
    ])

    if (system_dir / "surface.itp").exists():
        shutil.move(system_dir / "surface.itp", active_itp_dir / "surface.itp")

    # =========================================================================
    # 4) Orientation
    # =========================================================================
    import surfmartini.enzyme_tethered as orient_mod
    from surfmartini.enzyme_tethered import convert_pdb_to_gro

    enzyme_gro = system_dir / "enzyme_cg.gro"
    convert_pdb_to_gro(str(system_dir / "enzyme_cg.pdb"), str(enzyme_gro))

    orient_args = [
        "--surface", str(system_dir / "surface.gro"),
        "--enzyme",  str(enzyme_gro),
        "--out",     str(system_dir / "immobilized_system.gro"),
        "--dist",    str(args.dist),
    ]

    if args.resA:
        orient_args += ["--resA"] + [str(x) for x in args.resA]

    if args.resB:
        orient_args += ["--resB"] + [str(x) for x in args.resB]

    if args.anchor:
        orient_args += ["--anchor"] + [str(x) for x in args.anchor]

    orient_mod.main(orient_args)

    # =========================================================================
    # 5) Final system
    # =========================================================================
    import surfmartini.gomartini_system as gsm

    final_args = ["--outdir", str(simdir)]

    if args.resA:
        final_args += ["--resA"] + [str(x) for x in args.resA]

    if args.resB:
        final_args += ["--resB"] + [str(x) for x in args.resB]

    if args.anchor:
        final_args += ["--anchor"] + [str(x) for x in args.anchor]

    final_args += ["--moltype", args.moltype]


    gsm.main(final_args)

    shutil.rmtree(tmpdir)

    print("\n=====================================")
    print(" Finished Successfully!")
    print(" Output:", simdir)
    print("=====================================\n")


# ======================================================================
if __name__ == "__main__":
    main()
