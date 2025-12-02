#!/usr/bin/env python3
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        prog="martinisurf",
        description=(
            "MartiniSurf — Toolkit for building Martini/GōMartini "
            "surface–enzyme systems"
        ),
        add_help=False,
    )

    # ---------------------------------------------------------------
    # TOP-LEVEL HELP
    # ---------------------------------------------------------------
    parser.add_argument(
        "-h", "--help",
        action="store_true",
        help="Show detailed help (build flags)"
    )

    subparsers = parser.add_subparsers(dest="tool")

    # ================================================================
    # surface
    # ================================================================
    p_surface = subparsers.add_parser("surface")
    p_surface.add_argument("--bead", default="P4")
    p_surface.add_argument("--dx", type=float, default=0.47)
    p_surface.add_argument("--lx", type=float, required=True)
    p_surface.add_argument("--ly", type=float, required=True)
    p_surface.add_argument("--lz", type=float, default=5.0)
    p_surface.add_argument("--resname", default="SRF")
    p_surface.add_argument("--output", default="surface")
    p_surface.add_argument("--charge", type=float, default=0.0)

    # ================================================================
    # orient
    # ================================================================
    p_orient = subparsers.add_parser("orient")
    p_orient.add_argument("--surface", required=True)
    p_orient.add_argument("--enzyme", required=True)
    p_orient.add_argument("--out", default="Enzyme_Surface.gro")
    p_orient.add_argument("--dist", type=float, default=10.0)
    p_orient.add_argument("--display", choices=["on", "off"], default="off")

    # ---- MULTI-ANCHOR INPUT ----
    p_orient.add_argument(
        "--anchor",
        nargs="+",
        metavar=("GROUP", "RESID"),
        action="append",
        help=(
            "Define an anchor group as: GROUP RESID [RESID ...]. "
            "Repeat this flag to define multiple groups."
        ),
    )

    # ================================================================
    # system
    # ================================================================
    p_system = subparsers.add_parser("system")
    p_system.add_argument("--outdir", type=str, default="GoMartini_System")

    # ---- MULTI-ANCHOR INPUT ----
    p_system.add_argument(
        "--anchor",
        nargs="+",
        metavar=("GROUP", "RESID"),
        action="append",
        help=(
            "Define an anchor group as: GROUP RESID [RESID ...]. "
            "Repeat this flag to define multiple groups."
        ),
    )

    # ================================================================
    # build (pipeline)
    # ================================================================
    p_build = subparsers.add_parser("build")
    p_build.add_argument("--pdb", required=True)
    p_build.add_argument("--m2-args", nargs=argparse.REMAINDER)

    # ================================================================
    # PARSE TOP-LEVEL ARGS
    # ================================================================
    args = sys.argv[1:]

    # CASE 1 → top-level help shows the full BUILD parser
    if "-h" in args or "--help" in args:
        from martinisurf.pipeline import build_parser
        build_parser().print_help()
        return

    # CASE 2 → no arguments
    if len(args) == 0:
        parser.print_help()
        return

    # CASE 3 → subcommand + help
    if args[0] in ("surface", "orient", "system", "build") and (
        "-h" in args or "--help" in args
    ):
        parser.parse_args(args)
        return

    # CASE 4 → implicit "build"
    if args[0].startswith("--"):
        tool = "build"
        subcmd_args = args
    else:
        tool = args[0]
        subcmd_args = args[1:]

    # ================================================================
    # ROUTE MODULE
    # ================================================================
    if tool == "surface":
        import martinisurf.surface_builder as module
    elif tool == "orient":
        import martinisurf.enzyme_tethered as module
    elif tool == "system":
        import martinisurf.gomartini_system as module
    elif tool == "build":
        import martinisurf.pipeline as module
    else:
        print(f"❌ Unknown command: {tool}")
        return

    # FIX argv for module
    sys.argv = ["martinisurf " + tool] + subcmd_args

    module.main()


if __name__ == "__main__":
    main()
