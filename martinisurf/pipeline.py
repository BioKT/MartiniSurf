#!/usr/bin/env python3
"""
MartiniSurf Full Pipeline (Protein + DNA Compatible)

This script automates the complete workflow:

1. martinize2 or martinize-dna.py → Coarse-grained model
2. Surface                        → Build hexagonal Martini surface
3. Orientation                    → Position molecule above surface using anchor residues
4. GoMartini                      → Build full simulation system (topology, restraints, mdp)
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from martinisurf.utils.pdb_generation import load_clean_pdb
from martinisurf.utils.pdb_to_gro import pdb_to_gro

# ======================================================================
# PARSER
# ======================================================================
def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "MartiniSurf: Full pipeline for Protein/DNA coarse-graining, "
            "surface construction, orientation, and Gō–Martini system building."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --------------------------------------------------------------
    # 🧬 DNA MODE
    # --------------------------------------------------------------
    parser.add_argument(
        "--dna",
        action="store_true",
        help="Enable DNA workflow. Uses martinize-dna.py (Python 2.7) instead of martinize2."
    )

    parser.add_argument(
        "--dnatype",
        default="ds-stiff",
        help=(
            "Type of DNA elastic topology: "
            "ss | ds | ds-soft | ds-stiff | ss-stiff."
        )
    )

    # --------------------------------------------------------------
    # 📥 INPUT STRUCTURES
    # --------------------------------------------------------------
    parser.add_argument(
        "--pdb",
        required=True,
        help=(
            "Input structure. Can be: local file, RCSB ID (4 letters), "
            "or UniProt ID. Automatically cleaned before use."
        )
    )

    parser.add_argument(
        "--moltype",
        required=False,
        help="Name of the molecule (used for generated filenames and topology)."
    )

    parser.add_argument(
        "--merge",
        default=None,
        help="Merge chain identifiers (martinize-compatible syntax: e.g., A,B)."
    )

    # --------------------------------------------------------------
    # ⚗️ MARTINI FORCE FIELD
    # --------------------------------------------------------------
    parser.add_argument(
        "--ff",
        default="martini3001",
        help="Martini force field for protein martinization (ignored for DNA)."
    )

    # --------------------------------------------------------------
    # 🔗 POSITION RESTRAINTS
    # --------------------------------------------------------------
    parser.add_argument(
        "--p",
        choices=["none", "all", "backbone"],
        default="backbone",
        help="Which atoms to restrain during martinization."
    )

    parser.add_argument(
        "--pf",
        type=float,
        default=1000,
        help="Force constant for position restraints (kJ/mol/nm²)."
    )

    # --------------------------------------------------------------
    # 🧬 SECONDARY STRUCTURE
    # --------------------------------------------------------------
    parser.add_argument(
        "--dssp",
        action="store_true",
        help="Enable DSSP for secondary structure assignment."
    )

    parser.add_argument(
        "--no-dssp",
        dest="dssp",
        action="store_false",
        help="Disable DSSP (useful for DNA or custom ss assignment)."
    )

    parser.add_argument(
        "--ss",
        required=False,
        help="Provide your own secondary-structure file/string (overrides DSSP)."
    )

    parser.set_defaults(dssp=True)

    # --------------------------------------------------------------
    # 🧱 ELASTIC NETWORK
    # --------------------------------------------------------------
    parser.add_argument(
        "--elastic",
        action="store_true",
        help="Enable elastic network restraints during martinization."
    )

    parser.add_argument(
        "--ef",
        type=float,
        default=700,
        help="Elastic network force constant (kJ/mol/nm²)."
    )

    # --------------------------------------------------------------
    # 🧲 Gō–Martini OPTIONS (Protein only)
    # --------------------------------------------------------------
    parser.add_argument(
        "--go",
        nargs="?",
        const="auto",
        help=(
            "Enable the Gō–Martini potential. "
            "Use --go or --go auto for automatic contact detection."
        )
    )

    parser.add_argument(
        "--go-eps",
        type=float,
        default=9.414,
        help="Gō–Martini epsilon parameter (interaction strength)."
    )

    parser.add_argument(
        "--go-low",
        type=float,
        default=0.3,
        help="Lower bound for Gō–Martini contact distances (nm)."
    )

    parser.add_argument(
        "--go-up",
        type=float,
        default=1.1,
        help="Upper bound for Gō–Martini contact distances (nm)."
    )

    parser.add_argument(
        "--go-write-file",
        nargs="?",
        const=True,
        default=False,
        help="Export Gō contact map into an external file."
    )


    parser.add_argument(
        "--go-backbone",
        default="BB",
        help="Backbone bead name used for Gō–Martini (default: BB)."
    )

    parser.add_argument(
        "--go-atomname",
        default="CA",
        help="Atom name used for native contacts in Gō–Martini (default: CA)."
    )

    # --------------------------------------------------------------
    # 🔧 PROTEIN MODIFICATIONS
    # --------------------------------------------------------------
    parser.add_argument(
        "--cys",
        default="auto",
        help="Cysteine bridge detection: auto or explicit pairs."
    )

    parser.add_argument(
        "--mutate",
        nargs="+",
        default=[],
        help="Apply amino acid mutations before CG conversion (Protein only)."
    )

    parser.add_argument(
        "--maxwarn",
        type=int,
        default=0,
        help="Maximum warnings allowed from martinize2."
    )

    # --------------------------------------------------------------
    # 🌐 SURFACE GENERATION
    # --------------------------------------------------------------
    parser.add_argument(
        "--surface-bead",
        default="C1",
        help="Martini bead type used for the surface."
    )

    parser.add_argument(
        "--charge",
        type=int,
        default=0,
        help="Net charge per surface bead (integer)."
    )

    parser.add_argument(
        "--dx",
        type=float,
        default=0.47,
        help="Surface bead spacing (nm). Controls surface density."
    )

    parser.add_argument("--lx", type=float, help="Surface size in x direction (nm).")
    parser.add_argument("--ly", type=float, help="Surface size in y direction (nm).")

    parser.add_argument(
        "--surface",
        help="Provide an external surface GRO file instead of generating one."
    )

    # --------------------------------------------------------------
    # 📐 ORIENTATION
    # --------------------------------------------------------------
    parser.add_argument(
        "--anchor",
        nargs="+",
        action="append",
        metavar=("GROUP", "RESID"),
        help=(
            "Define anchor groups for orientation. "
            "Example: --anchor 1 8 10 11  (Group 1 uses residues 8,10,11). "
            "Repeat to define multiple groups."
        )
    )

    parser.add_argument(
        "--dist",
        type=float,
        default=10.0,
        help="Distance (nm) between the system and the surface."
    )

    # --------------------------------------------------------------
    # 📁 OUTPUT CONTROL
    # --------------------------------------------------------------
    parser.add_argument(
        "--outdir",
        default="Simulation_Files",
        help="Directory where all output files will be written."
    )

    parser.add_argument(
        "--m2-args",
        nargs=argparse.REMAINDER,
        help="Additional raw arguments passed directly to martinize2."
    )

    return parser

# ======================================================================
# UTILITY RUNNER
# ======================================================================
def run(cmd):
    print("\n▶ Running command:\n ", " ".join(cmd), "\n")
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Command failed.")
    print("✔ Done\n")


# ======================================================================
# MAIN PIPELINE
# ======================================================================
def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    # ---------------------------------------------------------
    # SETUP
    # ---------------------------------------------------------
    simdir = Path(args.outdir).resolve()
    if simdir.exists():
        shutil.rmtree(simdir)

    (simdir / "0_topology" / "system_itp").mkdir(parents=True)
    (simdir / "1_mdp").mkdir()
    (simdir / "2_system").mkdir()

    active_itp_dir = simdir / "0_topology" / "system_itp"
    system_dir     = simdir / "2_system"


    tmpdir = simdir / "_martinize_tmp"
    tmpdir.mkdir()

    mol = args.moltype if args.moltype else "DNA"

    system_cg_out = tmpdir / f"{mol}_cg.pdb"
    topfile_out   = tmpdir / f"{mol}_cg.top"

    # ---------------------------------------------------------
    # 1) CLEAN OR DOWNLOAD PDB
    # ---------------------------------------------------------
    pdb_abs = load_clean_pdb(args.pdb, workdir=simdir)

    # ---------------------------------------------------------
    # DNA branch → convert PDB → GRO (required by martinize-dna.py)
    # ---------------------------------------------------------
    if args.dna:
        dna_input_gro = tmpdir / "dna_input.gro"
        pdb_to_gro(str(pdb_abs), str(dna_input_gro))
        pdb_abs = dna_input_gro   # important: from now on use GRO file
        print(f"🧬 DNA mode: martinize-dna.py will use {pdb_abs}")

    # =====================================================================
    # 1) BUILD MARTINIZE COMMAND (Protein vs DNA)
    # =====================================================================
    if args.dna:
        print("🧬 DNA mode enabled → using martinize-dna.py")

        dna_script = Path(__file__).resolve().parent / "utils" / "martinize-dna.py"
        if not dna_script.exists():
            raise FileNotFoundError(f"martinize-dna.py not found at: {dna_script}")

        from martinisurf.utils.use_python2 import find_python2
        python2 = find_python2()

        martinize_cmd = [
            python2,
            str(dna_script),
            "-f", str(pdb_abs),
            "-x", str(system_cg_out),
            "-o", str(topfile_out),
            "-dnatype", args.dnatype,
            "-p", args.p.capitalize(),
            "-pf", str(args.pf),
        ]

        if args.elastic:
            martinize_cmd += ["-elastic", "-ef", str(args.ef)]

    else:
        print("🧬 Protein mode → using martinize2")

        martinize_cmd = [
            "martinize2",
            "-f", str(pdb_abs),
            "-x", str(system_cg_out),
            "-o", str(topfile_out),
            "-ff", args.ff,
            "-name", args.moltype,
            "-maxwarn", str(args.maxwarn),
        ]

        if args.p != "none":
            martinize_cmd += ["-p", args.p]
        martinize_cmd += ["-pf", str(args.pf)]

        if args.elastic:
            martinize_cmd += ["-elastic", "-ef", str(args.ef)]

    # =====================================================================
    # 2) DSSP / SS HANDLING (allowed for both Protein + DNA)
    # =====================================================================
    if args.ss:
        martinize_cmd += ["-ss", args.ss]

    elif args.dssp:
        dssp_bin = Path(__file__).resolve().parent / "dssp" / "mkdssp"
        if not dssp_bin.exists():
            raise FileNotFoundError(f"DSSP requested but binary not found: {dssp_bin}")
        print(f"✔ DSSP enabled → using: {dssp_bin}")
        martinize_cmd += ["-dssp", str(dssp_bin)]


    # =====================================================================
    # 3) MERGE, CYS (Protein + DNA) | MUTATE (Protein only)
    # =====================================================================
    if args.merge:
        martinize_cmd += ["-merge", args.merge]

    martinize_cmd += ["-cys", args.cys]

    if args.mutate:
        if args.dna:
            print("⚠ WARNING: --mutate ignored for DNA (not supported).")
        else:
            for mut in args.mutate:
                martinize_cmd += ["-mutate", mut]

    if not args.dna and args.m2_args:
        martinize_cmd += args.m2_args


    # =====================================================================
    # 4) GoMartini (Protein only)
    # =====================================================================
    if args.dna:
        print("🧬 DNA mode: GoMartini skipped.")
    else:
        if args.go is not None:
            martinize_cmd += ["-go"] if args.go == "auto" else ["-go", args.go]

            martinize_cmd += [
                "-go-eps", str(args.go_eps),
                "-go-low", str(args.go_low),
                "-go-up",  str(args.go_up),
                "-go-backbone", args.go_backbone,
                "-go-atomname", args.go_atomname,
            ]

            if args.go_write_file:
                martinize_cmd.append("-go-write-file")


    # =====================================================================
    # 5) RUN MARTINIZATION
    # =====================================================================
    run(martinize_cmd)

    for f in Path(".").glob("*.itp"):
        shutil.move(str(f), active_itp_dir / f.name)

    for f in Path(".").glob("*.ssd"):
        shutil.move(str(f), active_itp_dir / f.name)

    # Rename mol_0.itp → mol.itp (protein only)
    if not args.dna:
        src = active_itp_dir / f"{args.moltype}_0.itp"
        dst = active_itp_dir / f"{args.moltype}.itp"
        if src.exists():
            src.rename(dst)
            print(f"✔ Renamed {src.name} → {dst.name}")

    shutil.copy(system_cg_out, system_dir / f"{mol}_cg.pdb")


    # =====================================================================
    # 6) SURFACE GENERATION
    # =====================================================================
    surface_gro = system_dir / "surface.gro"

    if args.surface:
        shutil.copy(args.surface, surface_gro)
        possible_itp = Path(args.surface).with_suffix(".itp")
        if possible_itp.exists():
            shutil.copy(possible_itp, active_itp_dir / "surface.itp")
    else:
        print("\n🌐 Generating surface…\n")
        import martinisurf.surface_builder as sb
        sb.main([
            "--lx", str(args.lx),
            "--ly", str(args.ly),
            "--dx", str(args.dx),
            "--bead", args.surface_bead,
            "--charge", str(args.charge),
            "--output", str(system_dir / "surface"),
        ])

        surface_itp = system_dir / "surface.itp"
        if surface_itp.exists():
            shutil.move(surface_itp, active_itp_dir / "surface.itp")


    # =====================================================================
    # 7) ORIENTATION
    # =====================================================================
    import martinisurf.system_tethered as orient_mod
    from martinisurf.system_tethered import convert_pdb_to_gro

    system_gro = system_dir / f"{mol}_cg.gro"
    convert_pdb_to_gro(str(system_dir / f"{mol}_cg.pdb"), str(system_gro))

    orient_args = [
        "--surface", str(surface_gro),
        "--system",  str(system_gro),
        "--out",     str(system_dir / "immobilized_system.gro"),
        "--dist",    str(args.dist),
    ]

    if args.anchor:
        for group in args.anchor:
            orient_args += ["--anchor"] + [str(x) for x in group]

    orient_mod.main(orient_args)


    # =====================================================================
    # 8) FINAL GoMartini SYSTEM ASSEMBLY
    # =====================================================================
    import martinisurf.gromacs_inputs as gsm

    final_args = ["--outdir", str(simdir)]

    if not args.dna:
        final_args += ["--moltype", args.moltype]

    if args.anchor:
        for group in args.anchor:
            final_args += ["--anchor"] + [str(x) for x in group]

    gsm.main(final_args)

    shutil.rmtree(tmpdir)

    print("\n=====================================")
    print(" Finished Successfully!")
    print(" Output:", simdir)
    print("=====================================\n")


if __name__ == "__main__":
    main()
