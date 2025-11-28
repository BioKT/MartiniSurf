#!/usr/bin/env python3
"""
Generate a hexagonal surface and a minimal GROMACS topology (.itp).

This module creates a 2D hexagonal lattice of beads to represent
a flat surface for Martini or GōMartini simulations. It outputs:

- A .gro coordinate file
- A minimal .itp topology file
- A clean automatic copy of surface.itp into ActiveITP/
"""

import argparse
import math
import os
import shutil
from pathlib import Path
from typing import Iterable, List, Tuple


# ======================================================================
# Main entry point
# ======================================================================
def main(argv: Iterable[str] | None = None) -> None:
    """
    Entry point for the hexagonal-surface generator.

    Parameters
    ----------
    argv : iterable of str or None
        Command-line arguments. If None, argparse parses sys.argv.

    Notes
    -----
    This function also handles correct placement of the generated
    `surface.itp` file depending on whether the script is run in:
        - Standalone mode
        - MartiniSurf pipeline mode (inside Simulation/<...>)
    """
    parser = argparse.ArgumentParser(
        description="Generate a hexagonal surface and minimal surface.itp"
    )

    parser.add_argument("--bead", type=str, default="P4")
    parser.add_argument("--dx", type=float, default=0.47)
    parser.add_argument("--lx", type=float, required=True)
    parser.add_argument("--ly", type=float, required=True)
    parser.add_argument("--lz", type=float, default=5.0)
    parser.add_argument("--resname", type=str, default="SRF")
    parser.add_argument(
        "--output",
        type=str,
        default="surface",
        help="Basename or full path for the output files",
    )
    parser.add_argument("--charge", type=float, default=0.0)

    args = parser.parse_args(argv)

    # ---------------------------------------------------------
    # Resolve output directory and file base name
    # ---------------------------------------------------------
    outdir = os.path.dirname(args.output) or "."
    basename = os.path.basename(args.output)
    os.makedirs(outdir, exist_ok=True)

    # ---------------------------------------------------------
    # Generate hexagonal lattice parameters
    # ---------------------------------------------------------
    scale = args.dx / 0.142  # scaling from 0.142 nm reference
    a = 0.246 * scale  # primitive lattice constant

    # Unit-cell atom positions
    atoms_unit: List[Tuple[float, float, float]] = [
        (0.0, 0.0, 0.0),
        (a, 0.0, 0.0),
        (2 * a, 0.0, 0.0),
        (0.5 * a, (math.sqrt(3) / 2) * a, 0.0),
        (1.5 * a, (math.sqrt(3) / 2) * a, 0.0),
        (2.5 * a, (math.sqrt(3) / 2) * a, 0.0),
    ]

    Lx_cell = 3.0 * a
    Ly_cell = math.sqrt(3) * a

    nx = max(1, round(args.lx / Lx_cell))
    ny = max(1, round(args.ly / Ly_cell))

    Lx = nx * Lx_cell
    Ly = ny * Ly_cell

    print(f"• Building hexagonal surface: {Lx:.3f} × {Ly:.3f} × {args.lz:.3f} nm")

    # ---------------------------------------------------------
    # Build all atom coordinates
    # ---------------------------------------------------------
    atoms: List[Tuple[float, float, float]] = []
    for i in range(nx):
        for j in range(ny):
            dx = i * Lx_cell
            dy = j * Ly_cell
            for x, y, z in atoms_unit:
                atoms.append((x + dx, y + dy, z))

    # ---------------------------------------------------------
    # Write .gro file
    # ---------------------------------------------------------
    gro_path = Path(outdir, f"{basename}.gro")

    with gro_path.open("w") as fgro:
        fgro.write(f"Hexagonal surface {args.lx}x{args.ly} nm (d={args.dx:.2f} nm)\n")
        fgro.write(f"{len(atoms):5d}\n")
        for i, (x, y, z) in enumerate(atoms, start=1):
            fgro.write(
                f"{i:5d}{args.resname:<4}{args.bead:>4}{i:6d}"
                f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
            )
        fgro.write(f"{Lx:10.5f}{Ly:10.5f}{args.lz:10.5f}\n")

    print(f"✔ Wrote {gro_path}  ({len(atoms)} beads)")

    # ---------------------------------------------------------
    # Write minimal ITP
    # ---------------------------------------------------------
    itp_path = Path(outdir, f"{basename}.itp")

    with itp_path.open("w") as fitp:
        fitp.write(";;;;;; Minimal surface topology\n\n")
        fitp.write("[ moleculetype ]\n")
        fitp.write("; molname   nrexcl\n")
        fitp.write(f"  {args.resname}        1\n\n")
        fitp.write("[ atoms ]\n")
        fitp.write("; id  type     resnr  residu  atom  cgnr  charge\n")
        fitp.write(
            f"  1   {args.bead:<6}   1   {args.resname:<4}   C     1     {args.charge:.3f}\n"
        )

    print(f"✔ Wrote {itp_path}")

    # ---------------------------------------------------------
    # Copy surface.itp to ActiveITP/
    # ---------------------------------------------------------
    final_dst = resolve_activeitp_destination(outdir)
    final_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(itp_path, final_dst)

    print(f"✔ Copied surface.itp → {final_dst}")


# ======================================================================
# Helper: detect correct output destination for surface.itp
# ======================================================================
def resolve_activeitp_destination(outdir: str) -> Path:
    """
    Determine the correct destination path for surface.itp depending on
    standalone or pipeline execution.

    Parameters
    ----------
    outdir : str
        Directory where output files were generated.

    Returns
    -------
    Path
        Final destination for ActiveITP/surface.itp.
    """
    outdir_path = Path(outdir)

    # Default: standalone mode
    default_dst = outdir_path / "ActiveITP" / "surface.itp"

    # Pipeline detection
    if "Simulation" in str(outdir_path):
        parts = list(outdir_path.parts)
        if "Simulation" in parts:
            idx = parts.index("Simulation")
            simulation_root = Path(*parts[: idx + 1])

            topology_dir = simulation_root / "0_topology" / "ActiveITP"
            if topology_dir.exists():
                return topology_dir / "surface.itp"

    return default_dst


# ======================================================================
if __name__ == "__main__":
    main()
