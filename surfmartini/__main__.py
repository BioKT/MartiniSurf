#!/usr/bin/env python3
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        prog="martinisurf",
        description="SurfMartini — Toolkit for building Martini/GōMartini" 
        "surface–enzyme systems",
    )

    subparsers = parser.add_subparsers(dest="tool", required=True)

    # --------------------------
    # surface
    # --------------------------
    p_surface = subparsers.add_parser("surface")
    p_surface.add_argument("--bead", default="P4")
    p_surface.add_argument("--dx", type=float, default=0.47)
    p_surface.add_argument("--lx", type=float, required=True)
    p_surface.add_argument("--ly", type=float, required=True)
    p_surface.add_argument("--lz", type=float, default=5.0)
    p_surface.add_argument("--resname", default="SRF")
    p_surface.add_argument("--output", default="surface")
    p_surface.add_argument("--charge", type=float, default=0.0)

    # --------------------------
    # orient
    # --------------------------
    p_orient = subparsers.add_parser("orient")
    p_orient.add_argument("--surface", required=True)
    p_orient.add_argument("--enzyme", required=True)
    p_orient.add_argument("--out", default="Enzyme_Surface.gro")
    p_orient.add_argument("--resA", nargs="+", type=int)
    p_orient.add_argument("--resB", nargs="+", type=int)
    p_orient.add_argument("--anchor", nargs="+", type=int)
    p_orient.add_argument("--dist", type=float, default=10.0)
    p_orient.add_argument("--display", choices=["on", "off"], default="off")

    # --------------------------
    # system
    # --------------------------
    p_system = subparsers.add_parser("system")
    p_system.add_argument("--resA", nargs="+", type=int)
    p_system.add_argument("--resB", nargs="+", type=int)
    p_system.add_argument("--anchor", nargs="+", type=int)
    p_system.add_argument("--outdir", type=str, default="GoMartini_System")

    # --------------------------
    # build (pipeline)
    # --------------------------
    p_build = subparsers.add_parser("build")
    p_build.add_argument("--pdb", required=True)
    p_build.add_argument("--go", nargs="?", default="auto")
    p_build.add_argument("--eps", type=float, default=9.414)
    p_build.add_argument("--dssp", default="mkdssp")
    p_build.add_argument("--ff", default="martini3001")
    p_build.add_argument("--moltype", default="Active")
    p_build.add_argument("--merge", default=None)
    p_build.add_argument("--surface-bead", default="P4")
    p_build.add_argument("--dx", type=float, default=0.47)
    p_build.add_argument("--lx", type=float, required=True)
    p_build.add_argument("--ly", type=float, required=True)
    p_build.add_argument("--resA", nargs="+", type=int)
    p_build.add_argument("--resB", nargs="+", type=int)
    p_build.add_argument("--anchor", nargs="+", type=int)
    p_build.add_argument("--dist", type=float, default=10.0)
    p_build.add_argument("--outdir", default="SurfMartini_System")

    # --------------------------
    # Parse only the subcommand
    # --------------------------
    args = sys.argv[1:]  # everything except script name

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    tool = args[0]  # subcommand
    subcmd_args = args[1:]  # all remaining args go to the module

    # --------------------------
    # route to modules
    # --------------------------
    if tool == "surface":
        import surfmartini.surface_builder as module

    elif tool == "orient":
        import surfmartini.enzyme_tethered as module

    elif tool == "system":
        import surfmartini.gomartini_system as module

    elif tool == "build":
        import surfmartini.pipeline as module

    # --------------------------
    # FIX sys.argv for submodule
    # --------------------------
    sys.argv = ["martinisurf " + tool] + subcmd_args

    # --------------------------
    # Run module
    # --------------------------
    module.main()

if __name__ == "__main__":
    main()
