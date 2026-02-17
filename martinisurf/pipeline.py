#!/usr/bin/env python3
"""
MartiniSurf Full Pipeline (Protein + DNA Compatible)

Stable architecture:
• Protein → martinize2
• DNA     → martinize-dna.py
• Classical anchor mode
• Multi-linker mode
• Optional random surface linkers
"""

import argparse
import importlib.util
import os
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional

from martinisurf.utils.pdb_generation import load_clean_pdb
from martinisurf.utils.pdb_to_gro import pdb_to_gro


def _read_gro_first_resname(gro_path: str) -> str | None:
    with open(gro_path, "r") as fh:
        lines = fh.readlines()
    for line in lines[2:-1]:
        if len(line) >= 10:
            return line[5:10].strip() or None
    return None


def _read_gro_first_atomname(gro_path: str) -> str | None:
    with open(gro_path, "r") as fh:
        lines = fh.readlines()
    for line in lines[2:-1]:
        if len(line) >= 15:
            return line[10:15].strip() or None
    return None


def _read_gro_last_atomname(gro_path: str) -> str | None:
    with open(gro_path, "r") as fh:
        lines = fh.readlines()[2:-1]
    for line in reversed(lines):
        if len(line) >= 15:
            name = line[10:15].strip()
            if name:
                return name
    return None


def _read_gro_atomname_for_resid(gro_path: str, resid: int) -> str | None:
    with open(gro_path, "r") as fh:
        lines = fh.readlines()[2:-1]
    for line in lines:
        if len(line) < 15:
            continue
        try:
            gro_resid = int(line[0:5])
        except ValueError:
            continue
        if gro_resid == resid:
            name = line[10:15].strip()
            if name:
                return name
    return None


def _read_gro_atom_count(gro_path: str) -> int | None:
    with open(gro_path, "r") as fh:
        lines = fh.readlines()
    if len(lines) < 2:
        return None
    try:
        return int(lines[1].strip())
    except ValueError:
        return None


def _bead_size_class(bead_name: Optional[str]) -> str:
    if not bead_name:
        return "R"
    name = bead_name.strip().upper()
    if name.startswith("T"):
        return "T"
    if name.startswith("S"):
        return "S"
    return "R"


def _sigma_nm(is_dna: bool, ff_name: str, class_a: str, class_b: str) -> float:
    if is_dna or "martini2" in ff_name.lower():
        return 0.47

    pair = tuple(sorted((class_a, class_b)))
    table = {
        ("R", "R"): 0.47,
        ("R", "S"): 0.43,
        ("R", "T"): 0.395,
        ("S", "S"): 0.41,
        ("S", "T"): 0.365,
        ("T", "T"): 0.34,
    }
    return table.get(pair, 0.47)


def _normalize_merge_groups(merge_values: Optional[list[str]]) -> list[str]:
    if not merge_values:
        return []
    groups: list[str] = []
    for raw in merge_values:
        if raw is None:
            continue
        item = str(raw).strip()
        if not item:
            continue
        groups.append(item)
    return groups


def _validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if not args.surface and (args.lx is None or args.ly is None):
        parser.error("When --surface is not provided, both --lx and --ly are required.")

    if args.linker and not args.linker_group:
        parser.error("Linker mode requires at least one --linker-group.")

    if not args.linker and not args.anchor:
        parser.error("Provide --anchor in classical mode, or use --linker mode.")

    if args.anchor and args.linker:
        print("⚠ Both --anchor and --linker were provided. Linker mode will be used.")


def _print_config_summary(args: argparse.Namespace) -> None:
    mode = "DNA" if args.dna else "Protein"
    orient_mode = "linker" if args.linker else "anchor"
    print("\n=== MartiniSurf Configuration ===")
    print(f"Mode:            {mode}")
    print(f"PDB/Input:        {args.pdb}")
    print(f"Orientation:      {orient_mode}")
    print(f"Output:           {args.outdir}")
    if args.go and not args.dna:
        print("Go model:         enabled")
    if not args.dna:
        print(f"Max warnings:     {args.maxwarn}")
    if args.linker:
        group_count = len(args.linker_group) if args.linker_group else 0
        print(f"Linker groups:    {group_count}")
        print(f"Invert linker:    {args.invert_linker}")
    elif args.anchor:
        print(f"Anchor groups:    {len(args.anchor)}")
    if args.merge:
        print(f"Merge groups:     {', '.join(args.merge)}")
    print("=================================\n")


# ======================================================================
# PARSER
# ======================================================================

def build_parser():
    parser = argparse.ArgumentParser(
        prog="martinisurf",
        description=(
            "Build complete MartiniSurf systems for Protein/DNA on surfaces.\n"
            "Use ONE orientation mode: classical anchors (--anchor) or linker mode (--linker)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Examples:\n"
            "  Anchor mode:\n"
            "    martinisurf --pdb 1RJW --moltype Protein --lx 20 --ly 20 "
            "--anchor 1 8 10 --anchor 2 1025 1027\n"
            "  Linker mode:\n"
            "    martinisurf --dna --pdb 4C64.pdb --surface surface.gro "
            "--linker linker.gro --linker-group 1 1\n"
        ),
    )

    input_group = parser.add_argument_group("Input And Molecule")
    input_group.add_argument("--pdb", required=True, help="Local PDB path, RCSB ID (4 chars), or UniProt ID (6 chars).")
    input_group.add_argument("--moltype", help="Molecule name for protein topology output.")
    input_group.add_argument(
        "--go",
        action="store_true",
        help="Enable Go model in martinize2 protein mode.",
    )
    input_group.add_argument("--ff", default="martini3001", help="Force field name for martinize2 (protein mode).")
    input_group.add_argument("--dna", action="store_true", help="Enable DNA mode (uses martinize-dna.py).")
    input_group.add_argument("--dnatype", default="ds-stiff", help="DNA type for martinize-dna.py.")
    input_group.add_argument(
        "--merge",
        action="append",
        metavar="CHAINS",
        help=(
            "Merge chains during martinization. Example: --merge A,B,C,D "
            "(can be repeated). Use --merge all to merge every chain."
        ),
    )

    martinize_group = parser.add_argument_group("Martinization Controls")
    martinize_group.add_argument("--p", choices=["none", "all", "backbone"], default="backbone", help="Position restraints selection.")
    martinize_group.add_argument("--pf", type=float, default=1000, help="Position restraints force constant.")
    martinize_group.add_argument("--maxwarn", type=int, default=0, help="Allowed martinize2 warnings before abort.")
    martinize_group.add_argument("--dssp", action="store_true", help="Use DSSP during protein martinization.")
    martinize_group.add_argument("--no-dssp", dest="dssp", action="store_false", help="Disable DSSP during protein martinization.")
    martinize_group.add_argument("--elastic", action="store_true", help="Enable elastic network.")
    martinize_group.add_argument("--ef", type=float, default=700, help="Elastic network force constant.")
    parser.set_defaults(dssp=True)

    surface_group = parser.add_argument_group("Surface (Required if --surface is omitted)")
    surface_group.add_argument("--surface", help="Existing surface .gro file. If omitted, a surface is generated.")
    surface_group.add_argument(
        "--surface-mode",
        choices=["2-1", "4-1"],
        default="2-1",
        help="Surface lattice mode for generated surfaces.",
    )
    surface_group.add_argument("--lx", type=float, help="Surface size in X (nm) for generated surface.")
    surface_group.add_argument("--ly", type=float, help="Surface size in Y (nm) for generated surface.")
    surface_group.add_argument("--dx", type=float, default=0.47, help="Surface bead spacing (nm).")
    surface_group.add_argument("--surface-bead", default="C1", help="Surface bead type for generated surface.")
    surface_group.add_argument("--charge", type=int, default=0, help="Surface bead charge for generated surface.")

    anchor_group = parser.add_argument_group("Orientation: Classical Anchor Mode")
    anchor_group.add_argument(
        "--anchor",
        nargs="+",
        action="append",
        metavar=("GROUP", "RESID"),
        help="Anchor group: GROUP RESID [RESID ...]. Repeat to define multiple groups.",
    )
    anchor_group.add_argument("--dist", type=float, default=10.0, help="Anchor-to-surface target distance (A).")

    linker_group = parser.add_argument_group("Orientation: Linker Mode")
    linker_group.add_argument("--linker", help="Linker .gro file.")
    linker_group.add_argument(
        "--linker-group",
        nargs="+",
        action="append",
        metavar=("GROUP", "RESID"),
        help="Residue group(s) to attach linker(s): GROUP RESID [RESID ...].",
    )
    linker_group.add_argument("--linker-prot-dist", type=float, help="Linker-to-protein/DNA distance (A). Auto if omitted.")
    linker_group.add_argument("--linker-surf-dist", type=float, help="Linker-to-surface distance (A). Auto if omitted.")
    linker_group.add_argument("--invert-linker", action="store_true", help="Reverse linker bead order before attachment.")
    linker_group.add_argument("--surface-linkers", type=int, default=0, help="Add random extra linkers on the surface.")

    output_group = parser.add_argument_group("Output")
    output_group.add_argument("--outdir", default="Simulation_Files", help="Output directory for generated system files.")

    return parser


# ======================================================================
# RUNNER
# ======================================================================

def run(cmd, cwd=None):
    print("\n▶ Running:\n ", " ".join(cmd), "\n")
    res = subprocess.run(cmd, cwd=cwd)
    if res.returncode != 0:
        raise RuntimeError("Command failed.")
    print("✔ Done\n")


def _parse_version_triplet(text: str) -> tuple[int, int, int] | None:
    match = re.search(r"(\d+)\.(\d+)\.(\d+)", text or "")
    if not match:
        return None
    return int(match.group(1)), int(match.group(2)), int(match.group(3))


def _version_leq(found: tuple[int, int, int], limit: tuple[int, int, int]) -> bool:
    return found <= limit


def _is_dssp_binary_compatible(binary_path: str) -> bool:
    try:
        res = subprocess.run(
            [binary_path, "--version"],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return False

    version_text = f"{res.stdout}\n{res.stderr}"
    parsed = _parse_version_triplet(version_text)
    if parsed is None:
        # If version cannot be parsed, avoid hard-failing and let martinize2 retry logic handle it.
        print(f"⚠ Could not parse DSSP version for {binary_path}.")
        return True

    compatible = _version_leq(parsed, (3, 1, 4))
    if not compatible:
        print(
            f"⚠ DSSP at {binary_path} is version {parsed[0]}.{parsed[1]}.{parsed[2]}, "
            "but martinize2 is only compatible with DSSP <= 3.1.4."
        )
    return compatible


def _select_dssp_flags() -> list[str]:
    # Martinize2 recommendation: prefer mdtraj for secondary structure assignment.
    if importlib.util.find_spec("mdtraj") is not None:
        return ["-dssp"]

    # Fallback to binary DSSP only when compatible.
    dssp_env = os.environ.get("DSSP", "").strip()
    if dssp_env and _is_dssp_binary_compatible(dssp_env):
        return ["-dssp", dssp_env]

    system_mkdssp = shutil.which("mkdssp")
    if system_mkdssp and _is_dssp_binary_compatible(system_mkdssp):
        return ["-dssp", system_mkdssp]

    dssp_bin = Path(__file__).resolve().parent / "dssp" / "mkdssp"
    if dssp_bin.exists() and _is_dssp_binary_compatible(str(dssp_bin)):
        return ["-dssp", str(dssp_bin)]

    print("⚠ DSSP requested but no compatible setup found. Continuing without DSSP.")
    print("  Install `mdtraj` (recommended) or provide DSSP <= 3.1.4 via $DSSP.")
    return []


def _backup_existing_output_dir(simdir: Path) -> Path | None:
    if not simdir.exists():
        return None

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = simdir.with_name(f"{simdir.name}_backup_{timestamp}")
    suffix = 1
    while backup.exists():
        backup = simdir.with_name(f"{simdir.name}_backup_{timestamp}_{suffix}")
        suffix += 1

    shutil.move(str(simdir), str(backup))
    print(f"ℹ Existing output moved to backup: {backup}")
    return backup


# ======================================================================
# MAIN
# ======================================================================

def main(argv=None):

    parser = build_parser()
    args = parser.parse_args(argv)
    _validate_args(parser, args)
    _print_config_summary(args)
    merge_groups = _normalize_merge_groups(args.merge)

    simdir = Path(args.outdir).resolve()
    _backup_existing_output_dir(simdir)

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

    # ===============================================================
    # 1) CLEAN INPUT
    # ===============================================================
    pdb_abs = load_clean_pdb(args.pdb, workdir=simdir)

    # ===============================================================
    # 2) MARTINIZATION
    # ===============================================================

    if args.dna:
        print("🧬 DNA mode → using martinize-dna.py")

        dna_script = Path(__file__).resolve().parent / "utils" / "martinize-dna.py"
        from martinisurf.utils.use_python2 import find_python2
        python2 = find_python2()

        dna_input_gro = tmpdir / "dna_input.gro"
        pdb_to_gro(str(pdb_abs), str(dna_input_gro))

        martinize_cmd = [
            python2,
            str(dna_script),
            "-f", str(dna_input_gro),
            "-x", str(system_cg_out),
            "-o", str(topfile_out),
            "-dnatype", args.dnatype,
            "-p", args.p.capitalize(),
            "-pf", str(args.pf),
        ]
        for group in merge_groups:
            martinize_cmd += ["-merge", group]

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
            "-name", mol,
            "-maxwarn", str(args.maxwarn),
        ]
        for group in merge_groups:
            martinize_cmd += ["-merge", group]

        if args.p != "none":
            martinize_cmd += ["-p", args.p]

        martinize_cmd += ["-pf", str(args.pf)]

        if args.elastic:
            martinize_cmd += ["-elastic", "-ef", str(args.ef)]

        if args.go:
            martinize_cmd += ["-go"]

        if args.dssp:
            martinize_cmd += _select_dssp_flags()

    # martinize2 in some Colab/runtime setups fails when DSSP binary is present
    # but not functional. In that case retry once without DSSP.
    try:
        run(martinize_cmd, cwd=tmpdir)
    except RuntimeError:
        if (not args.dna) and args.dssp and "-dssp" in martinize_cmd:
            print("⚠ martinize2 failed with DSSP. Retrying without DSSP...")
            retry_cmd = martinize_cmd[:]
            dssp_idx = retry_cmd.index("-dssp")
            # Remove -dssp and optional binary path if present.
            del retry_cmd[dssp_idx]
            if dssp_idx < len(retry_cmd) and not retry_cmd[dssp_idx].startswith("-"):
                del retry_cmd[dssp_idx]
            run(retry_cmd, cwd=tmpdir)
        else:
            raise

    # Move ITP files
    for f in tmpdir.glob("*.itp"):
        shutil.move(str(f), active_itp_dir / f.name)

    shutil.copy(system_cg_out, system_dir / f"{mol}_cg.pdb")

    # ===============================================================
    # 3) SURFACE
    # ===============================================================

    surface_gro = system_dir / "surface.gro"

    if args.surface:

        # Copy GRO
        shutil.copy(args.surface, surface_gro)

        # 🔹 Copy surface.itp if exists
        possible_itp = Path(args.surface).with_suffix(".itp")
        if possible_itp.exists():
            shutil.copy(possible_itp, active_itp_dir / "surface.itp")
            print("✔ Copied surface.itp from provided surface file")
        else:
            print("⚠ No surface.itp found next to provided surface.gro")

    else:
        import martinisurf.surface_builder as sb

        sb.main([
            "--mode", args.surface_mode,
            "--lx", str(args.lx),
            "--ly", str(args.ly),
            "--dx", str(args.dx),
            "--bead", args.surface_bead,
            "--charge", str(args.charge),
            "--output", str(system_dir / "surface"),
        ])

        # 🔹 Move generated surface.itp into topology
        generated_itp = system_dir / "surface.itp"
        if generated_itp.exists():
            shutil.move(generated_itp, active_itp_dir / "surface.itp")
            print("✔ Moved generated surface.itp into topology")
        else:
            print("⚠ surface.itp not generated by surface_builder")

    # ===============================================================
    # 4) ORIENTATION
    # ===============================================================

    import martinisurf.system_tethered as orient_mod

    system_gro = system_dir / f"{mol}_cg.gro"
    pdb_to_gro(str(system_dir / f"{mol}_cg.pdb"), str(system_gro))

    orient_args = [
        "--surface", str(surface_gro),
        "--system", str(system_gro),
        "--out", str(system_dir / "immobilized_system.gro"),
    ]

    if args.linker:
        first_group_resid = None
        if args.linker_group and len(args.linker_group[0]) > 1:
            try:
                first_group_resid = int(args.linker_group[0][1])
            except ValueError:
                first_group_resid = None

        linker_first = _read_gro_first_atomname(args.linker)
        linker_last = _read_gro_last_atomname(args.linker)
        # Protein side is the first linker bead by default; swap when inverted.
        linker_head = linker_last if args.invert_linker else linker_first
        linker_tail = linker_first if args.invert_linker else linker_last
        target_atom = _read_gro_atomname_for_resid(str(system_gro), first_group_resid) if first_group_resid else None
        surface_atom = _read_gro_first_atomname(str(surface_gro))

        prot_sigma_nm = _sigma_nm(
            is_dna=args.dna,
            ff_name=args.ff,
            class_a=_bead_size_class(linker_head),
            class_b=_bead_size_class(target_atom),
        )
        surf_sigma_nm = _sigma_nm(
            is_dna=args.dna,
            ff_name=args.ff,
            class_a=_bead_size_class(linker_tail),
            class_b=_bead_size_class(surface_atom or args.surface_bead),
        )

        linker_prot_dist_ang = args.linker_prot_dist if args.linker_prot_dist is not None else prot_sigma_nm * 10.0
        linker_surf_dist_ang = args.linker_surf_dist if args.linker_surf_dist is not None else surf_sigma_nm * 10.0

        print(
            f"ℹ Linker distances (A): prot={linker_prot_dist_ang:.3f} "
            f"surf={linker_surf_dist_ang:.3f}"
        )

        orient_args += [
            "--linker-gro", str(args.linker),
            "--linker-prot-dist", str(linker_prot_dist_ang),
            "--linker-surf-dist", str(linker_surf_dist_ang),
        ]
        if args.invert_linker:
            orient_args += ["--invert-linker"]

        if args.surface_linkers > 0:
            orient_args += ["--surface-linkers", str(args.surface_linkers)]

        if args.linker_group:
            for group in args.linker_group:
                orient_args += ["--linker-group"] + [str(x) for x in group]

    elif args.anchor:
        orient_args += ["--dist", str(args.dist)]
        for group in args.anchor:
            orient_args += ["--anchor"] + [str(x) for x in group]

    orient_mod.main(orient_args)

    # ===============================================================
    # 5) FINAL GMX SYSTEM (CRITICAL FIX)
    # ===============================================================

    import martinisurf.gromacs_inputs as gsm

    final_args = ["--outdir", str(simdir)]

    if not args.dna:
        final_args += ["--moltype", mol]

    if args.linker:
        final_args += ["--use-linker"]
        linker_resname = _read_gro_first_resname(args.linker)
        if linker_resname:
            final_args += ["--linker-resname", linker_resname]
        linker_size = _read_gro_atom_count(args.linker)
        if linker_size and linker_size > 0:
            final_args += ["--linker-size", str(linker_size)]
        final_args += ["--linker-pull-init-prot", str(prot_sigma_nm)]
        final_args += ["--linker-pull-init-surf", str(surf_sigma_nm)]

        linker_itp = Path(args.linker).with_suffix(".itp")
        if not linker_itp.exists():
            raise FileNotFoundError(
                f"Linker mode requires the matching .itp next to linker GRO: {linker_itp}"
            )
        shutil.copy(linker_itp, active_itp_dir / linker_itp.name)
        final_args += ["--linker-itp-name", linker_itp.name]
        print(f"✔ Copied {linker_itp.name} into topology")

    # If classical anchor mode
    if args.anchor:
        for group in args.anchor:
            final_args += ["--anchor"] + [str(x) for x in group]

    # If linker mode → reuse linker groups as anchors
    elif args.linker_group:
        for group in args.linker_group:
            final_args += ["--anchor"] + [str(x) for x in group]

    gsm.main(final_args)

    shutil.rmtree(tmpdir)

    print("\n=====================================")
    print("✔ MartiniSurf Finished Successfully")
    print("=====================================\n")


if __name__ == "__main__":
    main()
