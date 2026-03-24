#!/usr/bin/env python3
import argparse
import sys


def _print_compact_package_help():
    print(
        "MartiniSurf CLI\n"
        "Compact package help.\n"
    )
    print("Usage:")
    print("  martinisurf [build flags]")
    print("  martinisurf build [build flags]")
    print("  martinisurf surface [surface flags]")
    print("  martinisurf orient [orientation flags]")
    print("  martinisurf system [topology/index flags]")
    print("")
    print("Main command groups:")
    print("  build   Main pipeline (protein, DNA, linker, complex-config, solvation/ionization).")
    print("  surface Surface generator.")
    print("  orient  Orientation module.")
    print("  system  Topology/index generation module.")
    print("")
    print("Detailed options by module:")
    print("  martinisurf build -h")
    print("  martinisurf surface -h")
    print("  martinisurf orient -h")
    print("  martinisurf system -h")
    print("")
    print("Optional verbose package help:")
    print("  martinisurf --full-help")


def _print_full_package_help():
    from martinisurf.pipeline import build_parser
    import martinisurf.surface_builder as surface_builder
    import martinisurf.system_tethered as system_tethered
    import martinisurf.gromacs_inputs as gromacs_inputs

    print("=== BUILD (main pipeline) ===")
    build_parser().print_help()
    print("")
    print("=== SURFACE MODULE ===")
    try:
        surface_builder.main(["-h"])
    except SystemExit:
        pass
    print("")
    print("=== ORIENT MODULE ===")
    try:
        system_tethered.main(["-h"])
    except SystemExit:
        pass
    print("")
    print("=== SYSTEM MODULE ===")
    try:
        gromacs_inputs.main(["-h"])
    except SystemExit:
        pass


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
        help="Show compact package help"
    )
    parser.add_argument(
        "--full-help",
        action="store_true",
        help="Show full verbose help for all modules"
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

    p_orient.add_argument(
        "--anchor",
        nargs="+",
        metavar=("GROUP", "RESID"),
        action="append",
        help="Define an anchor group as: GROUP RESID [RESID ...]. Repeat flag to add more groups."
    )

    # ================================================================
    # system
    # ================================================================
    p_system = subparsers.add_parser("system")
    p_system.add_argument("--outdir", type=str, default="GoMartini_System")

    p_system.add_argument(
        "--anchor",
        nargs="+",
        metavar=("GROUP", "RESID"),
        action="append",
        help="Define anchor groups for final assembly."
    )

    # ================================================================
    # build
    # ================================================================
    p_build = subparsers.add_parser("build")
    p_build.add_argument("--pdb", required=True)
    p_build.add_argument(
        "build_args",
        nargs=argparse.REMAINDER,
        help="Additional flags forwarded to the build pipeline"
    )

    # ================================================================
    # PARSE TOP-LEVEL ARGS
    # ================================================================
    args = sys.argv[1:]

    # CASE 1 → top-level help
    if "--full-help" in args:
        _print_full_package_help()
        return

    # CASE 2 → explicit subcommand help
    if args and args[0] in ("surface", "orient", "system", "build") and (
        "-h" in args or "--help" in args
    ):
        # Route to the module help (not top-level compact help).
        pass
    elif "-h" in args or "--help" in args:
        _print_compact_package_help()
        return

    # CASE 3 → no args
    if len(args) == 0:
        _print_compact_package_help()
        return

    # CASE 4 → implicit build mode
    if args[0].startswith("--"):
        tool = "build"
        subcmd_args = args[:]      # CRITICAL: forward ALL flags intact
    else:
        tool = args[0]
        subcmd_args = args[1:]

    # ================================================================
    # ROUTE MODULE
    # ================================================================
    if tool == "surface":
        import martinisurf.surface_builder as module

    elif tool == "orient":
        import martinisurf.system_tethered as module

    elif tool == "system":
        import martinisurf.gromacs_inputs as module

    elif tool == "build":
        import martinisurf.pipeline as module

    else:
        print(f"❌ Unknown command: {tool}")
        return

    # ================================================================
    # FIX argv FOR SUBMODULE
    # ================================================================
    sys.argv = ["martinisurf " + tool] + subcmd_args

    module.main()


if __name__ == "__main__":
    main()
