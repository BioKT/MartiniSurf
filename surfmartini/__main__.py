#!/usr/bin/env python3
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        prog="martinisurf",
        description="SurfMartini — Toolkit for building Martini/GōMartini surface–enzyme systems",
        add_help=False
    )

    # Custom help: top-level -h shows BUILD flags
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
    p_orient.add_argument("--resA", nargs="+", type=int)
    p_orient.add_argument("--resB", nargs="+", type=int)
    p_orient.add_argument("--anchor", nargs="+", type=int)
    p_orient.add_argument("--dist", type=float, default=10.0)
    p_orient.add_argument("--display", choices=["on", "off"], default="off")

    # ================================================================
    # system
    # ================================================================
    p_system = subparsers.add_parser("system")
    p_system.add_argument("--resA", nargs="+", type=int)
    p_system.add_argument("--resB", nargs="+", type=int)
    p_system.add_argument("--anchor", nargs="+", type=int)
    p_system.add_argument("--outdir", type=str, default="GoMartini_System")

    # ================================================================
    # build (pipeline) — minimal parser, full help inside pipeline
    # ================================================================
    p_build = subparsers.add_parser("build")

    # The build parser here ONLY defines basic structure.
    # Full help and full argument list is loaded from pipeline.build_parser().
    p_build.add_argument("--pdb", required=True)
    p_build.add_argument("--m2-args", nargs=argparse.REMAINDER)

    # ================================================================
    # Parse top-level CLI
    # ================================================================
    args = sys.argv[1:]

    # CASE 1 → top level help → show FULL BUILD HELP
    if "-h" in args or "--help" in args:
        from surfmartini.pipeline import build_parser
        build_parser().print_help()
        return

    # CASE 2 → no arguments → general help
    if len(args) == 0:
        parser.print_help()
        return

    # CASE 3 → Explicit help for subcommands
    if args[0] in ("surface", "orient", "system", "build") and (
        "-h" in args or "--help" in args
    ):
        # Normal argparse help for subcommands
        parser.parse_args(args)  # will print help automatically
        return

    # CASE 4 → If the first arg begins with "--", this is an implicit "build"
    if args[0].startswith("--"):
        tool = "build"
        subcmd_args = args
    else:
        tool = args[0]
        subcmd_args = args[1:]

    # ================================================================
    # Route execution to correct module
    # ================================================================
    if tool == "surface":
        import surfmartini.surface_builder as module

    elif tool == "orient":
        import surfmartini.enzyme_tethered as module

    elif tool == "system":
        import surfmartini.gomartini_system as module

    elif tool == "build":
        import surfmartini.pipeline as module

    else:
        print(f"❌ Unknown command: {tool}")
        return

    # ================================================================
    # Fix sys.argv for submodule execution
    # ================================================================
    sys.argv = ["martinisurf " + tool] + subcmd_args

    # Call module
    module.main()


if __name__ == "__main__":
    main()
