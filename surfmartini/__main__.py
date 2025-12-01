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
    # Force field
    p_build.add_argument("--ff", default="martini3001")
    # Merge chains
    p_build.add_argument("--merge", default=None)
    # Position restraints
    p_build.add_argument("--p", choices=["none", "all", "backbone"], default="none")
    p_build.add_argument("--pf", type=float, default=1000)
    # DSSP
    p_build.add_argument("--dssp", nargs="?", const="mkdssp", default="mkdssp")
    # Elastic network
    p_build.add_argument("--elastic", action="store_true")
    p_build.add_argument("--ef", type=float, default=700)
    # GoMartini options
    p_build.add_argument("--go", nargs="?", const="auto")
    p_build.add_argument("--go-eps", type=float, default=9.414)
    p_build.add_argument("--go-low", type=float, default=0.3)
    p_build.add_argument("--go-up", type=float, default=1.1)
    p_build.add_argument("--go-write-file", nargs="?", const=True, default=False)
    # Protein description
    p_build.add_argument("--cys", default="auto")
    p_build.add_argument("--mutate", nargs="+", default=[])
    # Passthrough for advanced Martinize2 args
    p_build.add_argument("--surface-bead", default="P4")
    p_build.add_argument("--dx", type=float, default=0.47)
    p_build.add_argument("--lx", type=float, required=True)
    p_build.add_argument("--ly", type=float, required=True)
    p_build.add_argument("--resA", nargs="+", type=int)
    p_build.add_argument("--resB", nargs="+", type=int)
    p_build.add_argument("--anchor", nargs="+", type=int)
    p_build.add_argument("--dist", type=float, default=10.0)
    p_build.add_argument("--outdir", default="SurfMartini_System")
    p_build.add_argument(
    "--m2-args",
    nargs=argparse.REMAINDER,
    help="Extra arguments passed directly to Martinize2"
    )


    # --------------------------
    # Parse tool with BUILD as default
    # --------------------------
    args = sys.argv[1:]

    # No arguments → show help
    if len(args) == 0:
        parser.print_help()
        return

    # Help flag at top level → show help
    if args[0] in ("-h", "--help"):
        parser.print_help()
        return

    # If first argument starts with "--", this means implicit build mode
    if args[0].startswith("--"):
        tool = "build"
        subcmd_args = args
    else:
        tool = args[0]
        subcmd_args = args[1:]


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
